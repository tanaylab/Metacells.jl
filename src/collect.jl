"""
Collect metacells.
"""
module Collect

export collect_metacells!

using ..Contracts
using ..Downsample

using DataAxesFormats
using Statistics
using StatsBase
using TanayLabUtilities

"""
    collect_metacells(
        daf::DafWriter;
        UMIs_regularization::AbstractFloat = $(DEFAULT.UMIs_regularization),
    )::Nothing

Given an assignment of cells to metacell, compute an estimation of the fraction of UMIs of each gene for each metacell.

The naive way to do this would be to just take the total UMIs of the gene out of the total UMIs of the metacell.
However, this method has a weakness; a single strong cell with a "very different" fraction of the gene will dominate
(that is, the method is sensitive to outliers).

Instead, we take the geometric mean of the fractions of the gene in the cells. This however raises a few issues we need
to deal with:

  - Genes with zero UMIs are a problem; to combat this, we add a `UMIs_regularization` factor when computing the
    fractions, take the geomean, and subtract the regularization at the end (so all-zero genes will still get a zero
    overall fraction).

  - We want to give more weight to cells with more UMIs, so we use a scaled geomean. The weight we give to each cell is
    the log of the total number of UMIs in it.

  - The geomean fractions of all the genes in a metacell do not sum to one, so we scale them. This has the unfortunate
    side effect that the end result does not obey a nice relation to the linear fraction of the gene in each of the
    cells; in particular, the final result might be higher than all of these per-cell fractions (ouch).

This raises an important point about the whole "fraction of UMIs of a gene in a cell" concept, which is that it is
highly dependent on the set of genes you pick to compute the fraction out of (that is, the denominator). This places
restrictions on how you should use these fractions:

  - You can't just compare these fractions between two arbitrary data sets, as the denominators aren't the same.
    You must identify the set of common genes, and renormalize the fractions to sum to one in this subset, in both data
    sets. Only then can you meaningfully compare the results.

  - Even if you compare two data sets with the same set of genes (or even two subsets of the same data set, such as two
    "cell types"), if a specific "gene program" has a very high total expression in only one of them, then all other
    genes will artificially appear to be lower.

    This is why we recommend excluding ribosomal genes from the data sets; they can take up anything between almost none
    to over two thirds of the total UMIs, and this varies between "cell types". This means that if they are included in
    the denominator, the fractions of otherwise "identical" genes will appear to differ by a factor of up to 3X

    Luckily, most gene programs we deal with have a total expression of a few percent at most, so this effect is
    negligible. However if you identify a gene program with a total expression higher than, say, 10%, you should
    consider its effect on the rest of the genes.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        gene_cell_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_total_UMIs_vector(GuaranteedOutput),
        gene_metacell_total_UMIs_matrix(GuaranteedOutput),
        gene_metacell_fraction_matrix(GuaranteedOutput),
    ],
) function collect_metacells!(
    daf::DafWriter;
    UMIs_regularization::AbstractFloat = 1 / 16,
    min_downsamples::Integer = 750, # TODOX function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = 0.05, # TODOX function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = 0.5, # TODOX function_default(downsamples, :max_downsamples_quantile),
    rng::AbstractRNG = default_rng(),
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    metacell_name_per_cell = get_vector(daf, "cell", "metacell").array
    name_per_metacell = axis_vector(daf, "metacell")

    UMIs_per_gene_per_cell = get_matrix(daf, "gene", "cell", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    total_UMIs_per_gene_per_metacell = zeros(UInt32, n_genes, n_metacells)
    total_UMIs_per_metacell = zeros(UInt32, n_metacells)
    fraction_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)

    parallel_loop_with_rng(n_metacells; rng) do metacell_index, rng
        metacell_name = name_per_metacell[metacell_index]
        index_per_cell_in_metacell = findall(metacell_name_per_cell .== metacell_name)
        n_cells_in_metacell = length(index_per_cell_in_metacell)

        @views UMIs_per_gene_per_cell_in_metacell = UMIs_per_gene_per_cell[:, index_per_cell_in_metacell]

        total_UMIs_in_metacell_per_gene = vec(sum(UMIs_per_gene_per_cell_in_metacell; dims = 2))
        @assert_vector(total_UMIs_in_metacell_per_gene, n_genes)
        total_UMIs_per_gene_per_metacell[:, metacell_index] .= total_UMIs_in_metacell_per_gene

        @views total_UMIs_in_metacell_per_cell_in_metacell = total_UMIs_per_cell[index_per_cell_in_metacell]

        total_UMIs_in_metacell = sum(total_UMIs_in_metacell_per_cell_in_metacell)
        @assert sum(total_UMIs_in_metacell_per_gene) == total_UMIs_in_metacell

        total_UMIs_per_metacell[metacell_index] = total_UMIs_in_metacell

        downsample_UMIs_of_cells = downsamples(
            total_UMIs_in_metacell_per_cell_in_metacell;
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        downsampled_UMIs_per_gene_per_cell_in_metacell =
            downsample(UMIs_per_gene_per_cell_in_metacell, downsample_UMIs_of_cells; dims = 2, rng)

        total_downsampled_UMIs_per_cell_in_metacell = vec(sum(downsampled_UMIs_per_gene_per_cell_in_metacell; dims = 1))
        @assert_vector(total_downsampled_UMIs_per_cell_in_metacell, n_cells_in_metacell)
        @assert all(total_downsampled_UMIs_per_cell_in_metacell .<= downsample_UMIs_of_cells)
        @assert all(total_downsampled_UMIs_per_cell_in_metacell .<= total_UMIs_in_metacell_per_cell_in_metacell)

        weight_per_cell_in_metacell = log2.(total_UMIs_in_metacell_per_cell_in_metacell)
        fraction_per_gene_per_cell_in_metacell =
            downsampled_UMIs_per_gene_per_cell_in_metacell ./ total_downsampled_UMIs_per_cell_in_metacell

        regularization_per_cell_in_metacell = UMIs_regularization ./ total_downsampled_UMIs_per_cell_in_metacell
        fraction_per_gene_per_cell_in_metacell .+= regularization_per_cell_in_metacell

        fraction_per_gene =
            weighted_geomean(fraction_per_gene_per_cell_in_metacell, weight_per_cell_in_metacell; dims = 2)
        @assert_vector(fraction_per_gene, n_genes)

        regularization_in_metacell = weighted_geomean(regularization_per_cell_in_metacell, weight_per_cell_in_metacell)
        fraction_per_gene .-= regularization_in_metacell

        total_downsampled_UMIs_of_genes = vec(sum(downsampled_UMIs_per_gene_per_cell_in_metacell; dims = 2))
        @assert all(isapprox.(fraction_per_gene[total_downsampled_UMIs_of_genes .== 0], 0))
        fraction_per_gene[total_downsampled_UMIs_of_genes .== 0] .= 0
        @assert all(fraction_per_gene .>= 0)
        fraction_per_gene ./= sum(fraction_per_gene)

        fraction_per_gene_per_metacell[:, metacell_index] .= fraction_per_gene
    end

    set_matrix!(daf, "gene", "metacell", "total_UMIs", total_UMIs_per_gene_per_metacell)
    set_vector!(daf, "metacell", "total_UMIs", total_UMIs_per_metacell)
    set_matrix!(daf, "gene", "metacell", fraction_per_gene_per_metacell)

    return nothing
end

function weighted_geomean(values::AbstractVector{<:Real}, weights::AbstractVector{<:Real})::AbstractFloat
    return 2 .^ (mean(log2.(values) .* weights) ./ sum(weights))
end

function weighted_geomean(
    values::AbstractMatrix{<:Real},
    weights::AbstractVector{<:Real};
    dims::Integer,
)::AbstractVector{<:AbstractFloat}
    if dims == 1
        return vec(2 .^ (mean(log2.(values) .* weights; dims) ./ sum(weights)))
    else
        return vec(2 .^ (mean(log2.(values) .* transpose(weights); dims) ./ sum(weights)))
    end
end

end
