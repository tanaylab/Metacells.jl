"""
Do simple metacells analysis.
"""
module AnalyzeMetacells

export compute_cells_types_by_metacells!
export compute_metacells_covered_UMIs!
export compute_metacells_euclidean_distances!
export compute_metacells_genes_covered_fractions!
export compute_metacells_genes_geomean_fractions!
export compute_metacells_genes_linear_fractions!
export compute_metacells_genes_log_covered_fractions!
export compute_metacells_genes_log_geomean_fractions!
export compute_metacells_genes_log_linear_fractions!
export compute_metacells_genes_UMIs!
export compute_metacells_max_skeleton_fold_distances!
export compute_metacells_n_cells!
export compute_metacells_total_UMIs!
export compute_metacells_types_by_cells!

using Base.Threads
using DataAxesFormats
using Distances
using Distributions
using TanayLabUtilities
using Random
using StatsBase

using ..Defaults
using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.cell_type_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_metacell_euclidean_skeleton_distance
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.metacell_type_vector

"""
    function compute_metacells_n_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of metacells per metacell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [cell_metacell_vector(RequiredInput), metacell_n_cells_vector(GuaranteedOutput)],
) function compute_metacells_n_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells_per_metacell = daf["/ cell : index @ metacell ! %> Count"].array
    set_vector!(daf, "metacell", "n_cells", n_cells_per_metacell; overwrite)
    return nothing
end

"""
    function compute_metacells_genes_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

An estimated linear fraction of the UMIs of each gene in each metacell. This is just the total UMIs of the gene
in the metacell divided by the total UMIs of the metacell, which is the "best" estimate assuming multinomial sampling
noise. However, this is sensitive to a few cells with very high expression levels ("bursty" genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_gene_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_linear_fractions(daf; qualifier = "linear", genes_mask = "!is_excluded", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_log_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of each gene in each metacell. This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "linear", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

An estimated linear fraction of the UMIs of each covered gene in each metacell. By considering only the covered genes
this avoid the impact of highly-expressed lateral genes (e.g., cell cycle). Otherwise is similar to
[`compute_metacells_genes_linear_fractions!`](@ref).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_gene_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_linear_fractions(daf; qualifier = "covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_log_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each metacell. This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_gene_log_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "covered", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of UMIs used to estimate the fraction of each gene in each metacell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        metacell_gene_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")
    n_genes = axis_length(daf, "gene")

    UMIs_per_gene_per_cell = get_matrix(daf, "gene", "cell", "UMIs").array

    UMIs_per_gene_per_metacell = Matrix{UInt32}(undef, n_genes, n_metacells)
    metacell_per_cell = get_vector(daf, "cell", "metacell").array

    @threads :greedy for metacell_index in 1:n_metacells
        metacell_name = name_per_metacell[metacell_index]
        metacell_cell_indices = findall(metacell_per_cell .== metacell_name)
        @assert length(metacell_cell_indices) > 0

        @views UMIs_per_gene_per_metacell_cell = UMIs_per_gene_per_cell[:, metacell_cell_indices]
        UMIs_per_gene_per_metacell[:, metacell_index] .= vec(sum(UMIs_per_gene_per_metacell_cell; dims = 2))
    end

    set_matrix!(daf, "gene", "metacell", "UMIs", bestify(UMIs_per_gene_per_metacell); overwrite)  # NOJET
    return nothing
end

"""
    function compute_metacells_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of UMIs used to estimate the fraction of all the genes in each metacell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_total_UMIs_vector(GuaranteedOutput),
    ],
) function compute_metacells_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_UMIs(daf; qualifier = "total", genes_mask = "!is_excluded", overwrite)
    return nothing
end

"""
    function compute_metacells_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of covered genes per metacell.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_metacells_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_UMIs(daf; qualifier = "covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    compute_metacells_genes_geomean_fractions!(
        daf::DafWriter;
        UMIs_regularization::AbstractFloat = $(DEFAULT.UMIs_regularization),
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given an assignment of cells to metacell, compute an geomean estimation of the fraction of UMIs of each gene for each
metacell.

The linear way to do this would be to just take the total UMIs of the gene out of the total UMIs of the metacell.
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
        gene_is_excluded_vector(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_gene_geomean_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_geomean_fractions!(
    daf::DafWriter;
    UMIs_regularization::AbstractFloat = 1 / 16,
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    metacell_name_per_cell = get_vector(daf, "cell", "metacell").array
    name_per_metacell = axis_vector(daf, "metacell")

    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    UMIs_per_included_gene_per_cell = daf["/ gene &! is_excluded / cell : UMIs"].array

    included_gene_indices = daf["/ gene &! is_excluded : index"].array
    n_included_genes = length(included_gene_indices)

    geomean_fraction_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)

    parallel_loop_with_rng(1:n_metacells; rng) do metacell_index, rng
        metacell_name = name_per_metacell[metacell_index]
        index_per_metacell_cell = findall(metacell_name_per_cell .== metacell_name)
        n_metacell_cells = length(index_per_metacell_cell)

        @views UMIs_per_included_gene_per_metacell_cell = UMIs_per_included_gene_per_cell[:, index_per_metacell_cell]

        total_UMIs_in_metacell_per_included_gene = vec(sum(UMIs_per_included_gene_per_metacell_cell; dims = 2))
        @assert_vector(total_UMIs_in_metacell_per_included_gene, n_included_genes)

        @views total_UMIs_in_metacell_per_metacell_cell = total_UMIs_per_cell[index_per_metacell_cell]

        total_UMIs_in_metacell = sum(total_UMIs_in_metacell_per_metacell_cell)
        @assert sum(total_UMIs_in_metacell_per_included_gene) == total_UMIs_in_metacell

        downsample_UMIs_of_cells = downsamples(
            total_UMIs_in_metacell_per_metacell_cell;
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        downsampled_UMIs_per_included_gene_per_metacell_cell =
            downsample(UMIs_per_included_gene_per_metacell_cell, downsample_UMIs_of_cells; dims = 1, rng)

        total_downsampled_UMIs_per_metacell_cell =
            vec(sum(downsampled_UMIs_per_included_gene_per_metacell_cell; dims = 1))
        @assert_vector(total_downsampled_UMIs_per_metacell_cell, n_metacell_cells)

        weight_per_metacell_cell = log2.(total_UMIs_in_metacell_per_metacell_cell)
        fraction_per_included_gene_per_metacell_cell =
            downsampled_UMIs_per_included_gene_per_metacell_cell ./ transpose(total_downsampled_UMIs_per_metacell_cell)

        regularization_per_metacell_cell = UMIs_regularization ./ total_downsampled_UMIs_per_metacell_cell
        fraction_per_included_gene_per_metacell_cell .+= transpose(regularization_per_metacell_cell)

        geomean_fraction_per_included_gene =
            weighted_geomean(fraction_per_included_gene_per_metacell_cell, weight_per_metacell_cell; dims = 2)
        @assert_vector(geomean_fraction_per_included_gene, n_included_genes)

        regularization_in_metacell = weighted_geomean(regularization_per_metacell_cell, weight_per_metacell_cell)
        geomean_fraction_per_included_gene .-= regularization_in_metacell

        total_downsampled_UMIs_of_included_genes =
            vec(sum(downsampled_UMIs_per_included_gene_per_metacell_cell; dims = 2))
        @assert all(
            isapprox.(
                geomean_fraction_per_included_gene[total_downsampled_UMIs_of_included_genes .== 0],
                0;
                atol = 1e-6,
            ),
        )
        geomean_fraction_per_included_gene[total_downsampled_UMIs_of_included_genes .== 0] .= 0
        @assert all(geomean_fraction_per_included_gene .>= 0)
        geomean_fraction_per_included_gene ./= sum(geomean_fraction_per_included_gene)

        geomean_fraction_per_gene_per_metacell[included_gene_indices, metacell_index] .=
            geomean_fraction_per_included_gene
        return nothing
    end

    set_matrix!(daf, "gene", "metacell", "geomean_fraction", bestify(geomean_fraction_per_gene_per_metacell); overwrite)

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

"""
    function compute_metacells_genes_log_geomean_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the estimated geomean fraction of the UMIs of each gene in each metacell. This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_geomean_fraction_matrix(RequiredInput),
        metacell_gene_log_geomean_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_geomean_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "geomean", overwrite)
    return nothing
end

function do_compute_metacells_genes_linear_fractions(
    daf::DafWriter;
    qualifier::AbstractString,
    genes_mask::AbstractString,
    overwrite::Bool,
)::Nothing
    UMIs_per_metacell_per_mask_gene = daf["/ metacell / gene &$(genes_mask) : UMIs"].array
    total_mask_UMIs_per_metacell = sum(UMIs_per_metacell_per_mask_gene; dims = 2)

    n_metacells = axis_length(daf, "metacell")
    n_genes = axis_length(daf, "gene")
    mask_gene_indices = daf["/ gene &$(genes_mask) : index"].array
    linear_fraction_per_metacell_per_gene = zeros(Float32, n_metacells, n_genes)
    linear_fraction_per_metacell_per_gene[:, mask_gene_indices] .=
        UMIs_per_metacell_per_mask_gene ./ total_mask_UMIs_per_metacell

    set_matrix!(  # NOJET
        daf,
        "metacell",
        "gene",
        "$(qualifier)_fraction",
        bestify(linear_fraction_per_metacell_per_gene);
        overwrite,
    )  # NOJET

    return nothing
end

function do_compute_metacells_genes_log_fractions(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat,
    qualifier::AbstractString,
    overwrite::Bool,
)::Nothing
    if qualifier == ""
        name = "fraction"
    else
        name = "$(qualifier)_fraction"
    end
    fraction_per_gene_per_metacell = daf["/ gene / metacell : $(name)"].array

    empty_dense_matrix!(
        daf,
        "gene",
        "metacell",
        "log_$(name)",
        Float32;
        overwrite,
    ) do log_fraction_per_gene_per_metacell
        log_fraction_per_gene_per_metacell .= log2.(fraction_per_gene_per_metacell .+ gene_fraction_regularization)
        return nothing
    end

    return nothing
end

function do_compute_metacells_UMIs(
    daf::DafWriter;
    qualifier::AbstractString,
    genes_mask::AbstractString,
    scale_per_gene::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
    overwrite::Bool,
)::Nothing
    UMIs_per_gene_per_metacell = daf["/ gene &$(genes_mask) / metacell : UMIs"].array

    if scale_per_gene !== nothing
        UMIs_per_gene_per_metacell = round.(UMIs_per_gene_per_metacell .* scale_per_gene)
    end

    UMIs_per_metacell = UInt32.(vec(sum(UMIs_per_gene_per_metacell; dims = 1)))
    set_vector!(daf, "metacell", "$(qualifier)_UMIs", UMIs_per_metacell; overwrite)

    return nothing
end

"""
    function compute_cells_types_by_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The type of each cell, based on the type of the metacell it belongs to.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_type_vector(RequiredInput),
        cell_type_vector(GuaranteedOutput),
    ],
) function compute_cells_types_by_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_cell = daf["/ cell : metacell ?? '' => type"].array
    set_vector!(daf, "cell", "type", type_per_cell; overwrite)
    return nothing
end

"""
    function compute_metacells_types_by_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The type of each metacell, based on is cell types. This assumes that each cell was assigned type, either because we
imported the data from a cell-based atlas, or from previous computed metacells.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_type_vector(RequiredInput),
        metacell_type_vector(GuaranteedOutput),
    ],
) function compute_metacells_types_by_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_metacell = daf["/ cell : type @ metacell ! %> Mode"].array
    set_vector!(daf, "metacell", "type", type_per_metacell; overwrite)
    return nothing
end

"""
    function compute_metacells_euclidean_distances!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The Euclidean distance between the log of the covered fraction of the skeleton genes between the metacells.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        metacell_gene_log_covered_fraction_matrix(RequiredInput),
        metacell_metacell_euclidean_skeleton_distance(GuaranteedOutput),
    ],
) function compute_metacells_euclidean_distances!(daf::DafWriter; overwrite::Bool = false)::Nothing
    log_covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : log_covered_fraction"].array
    distances_between_metacells = pairwise(Euclidean(), log_covered_fraction_per_skeleton_per_metacell)  # NOJET
    set_matrix!(daf, "metacell", "metacell", "euclidean_skeleton_distance", distances_between_metacells; overwrite)
    return nothing
end

"""
    function compute_metacells_max_skeleton_fold_distances!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The maximal significant fold factor between the covered fraction of skeleton genes between the metacells. This uses
heuristics to require the fold factor be based on a sufficient number of UMIs to be robust.

The fold factor is log (base 2) of the gene expression using the `gene_fraction_regularization`. For computing this fold
factor, we ignore "insignificant" genes whose total UMIs in the compared metacells isn't at least
`min_significant_gene_UMIs`. We also we reduce the distance using the `fold_confidence` based on the number of UMIs used
to estimate the expression in the metacells.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_covered_UMIs_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_metacell_max_skeleton_fold_distance(GuaranteedOutput),
    ],
) function compute_metacells_max_skeleton_fold_distances!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_significant_gene_UMIs >= 0
    @assert 0 < fold_confidence < 1

    covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : covered_fraction"].array
    covered_UMIs_per_metacell = get_vector(daf, "metacell", "covered_UMIs").array

    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_covered_fractions_per_skeleton_per_metacells =  # NOJET
        confidence_stdevs .* sqrt.(covered_fraction_per_skeleton_per_metacell .* transpose(covered_UMIs_per_metacell)) ./
        transpose(covered_UMIs_per_metacell)

    low_log_covered_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            max.(
                covered_fraction_per_skeleton_per_metacell .- confidence_covered_fractions_per_skeleton_per_metacells,
                0.0,
            ) .+ gene_fraction_regularization,
        )

    high_log_covered_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            covered_fraction_per_skeleton_per_metacell .+ confidence_covered_fractions_per_skeleton_per_metacells .+
            gene_fraction_regularization,
        )

    UMIs_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : UMIs"].array

    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)
    distances_between_metacells[1, 1] = 0.0

    @threads :greedy for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views base_metacell_UMIs_per_skeleton = vec(UMIs_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_low_log_covered_fraction_per_skeleton =
            vec(low_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_high_log_covered_fraction_per_skeleton =
            vec(high_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])

        n_other_metacells = base_metacell_index - 1
        other_metacells_indices = 1:(base_metacell_index - 1)
        @views UMIs_per_skeleton_per_other_metacells = UMIs_per_skeleton_per_metacell[:, other_metacells_indices]
        @views low_log_covered_fraction_per_skeleton_per_other_metacells =
            low_log_covered_fraction_per_skeleton_per_metacell[:, other_metacells_indices]
        @views high_log_covered_fraction_per_skeleton_per_other_metacells =
            high_log_covered_fraction_per_skeleton_per_metacell[:, other_metacells_indices]

        significant_fold_per_skeleton_per_other_metacell = confident_gene_distance.(
            min_significant_gene_UMIs,
            base_metacell_UMIs_per_skeleton,
            base_metacell_low_log_covered_fraction_per_skeleton,
            base_metacell_high_log_covered_fraction_per_skeleton,
            UMIs_per_skeleton_per_other_metacells,
            low_log_covered_fraction_per_skeleton_per_other_metacells,
            high_log_covered_fraction_per_skeleton_per_other_metacells,
        )
        @assert_matrix(significant_fold_per_skeleton_per_other_metacell, n_skeletons, n_other_metacells, Columns)

        distances_between_base_and_other_metacells =
            vec(maximum(significant_fold_per_skeleton_per_other_metacell; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, n_other_metacells)

        distances_between_metacells[other_metacells_indices, base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, other_metacells_indices] .=
            distances_between_base_and_other_metacells
    end

    set_matrix!(daf, "metacell", "metacell", "max_skeleton_fold_distance", distances_between_metacells; overwrite)
    return nothing
end

@inline function confident_gene_distance(  # UNTESTED
    min_significant_gene_UMIs::Integer,
    base_metacell_total_UMIs_of_gene::Integer,
    base_metacell_low_log_covered_fraction_of_gene::AbstractFloat,
    base_metacell_high_log_covered_fraction_of_gene::AbstractFloat,
    other_metacell_total_UMIs_of_gene::Integer,
    other_metacell_low_log_covered_fraction_of_gene::AbstractFloat,
    other_metacell_high_log_covered_fraction_of_gene::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = base_metacell_total_UMIs_of_gene + other_metacell_total_UMIs_of_gene
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = base_metacell_high_log_covered_fraction_of_gene < other_metacell_high_log_covered_fraction_of_gene

    highest_low_log_covered_fraction_of_gene =
        is_base_low * base_metacell_high_log_covered_fraction_of_gene +
        !is_base_low * other_metacell_high_log_covered_fraction_of_gene

    lowest_high_log_covered_fraction_of_gene =
        is_base_low * other_metacell_low_log_covered_fraction_of_gene +
        !is_base_low * base_metacell_low_log_covered_fraction_of_gene

    confident_gap = lowest_high_log_covered_fraction_of_gene - highest_low_log_covered_fraction_of_gene

    return (is_significant * max(confident_gap, 0.0))
end

end  # module
