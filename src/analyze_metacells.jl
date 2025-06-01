"""
Do simple metacells analysis.
"""
module AnalyzeMetacells

export compute_cells_linear_covered_metacell_cross_entropy!
export compute_cells_linear_covered_metacell_kl_divergence!
export compute_cells_linear_metacell_cross_entropy!
export compute_cells_linear_metacell_kl_divergence!
export compute_cells_types_by_metacells!
export compute_metacells_covered_UMIs!
export compute_metacells_genes_geomean_fractions!
export compute_metacells_genes_linear_covered_fractions!
export compute_metacells_genes_linear_fractions!
export compute_metacells_genes_log_geomean_fractions!
export compute_metacells_genes_log_linear_covered_fractions!
export compute_metacells_genes_log_linear_fractions!
export compute_metacells_genes_log_scaled_linear_covered_fractions!
export compute_metacells_genes_log_scaled_linear_fractions!
export compute_metacells_genes_scaled_linear_covered_fractions!
export compute_metacells_genes_scaled_linear_fractions!
export compute_metacells_genes_UMIs!
export compute_metacells_mean_cells_linear_covered_cross_entropy!
export compute_metacells_mean_cells_linear_covered_kl_divergence!
export compute_metacells_mean_cells_linear_cross_entropy!
export compute_metacells_mean_cells_linear_kl_divergence!
export compute_metacells_n_cells!
export compute_metacells_scaled_covered_UMIs!
export compute_metacells_scaled_total_UMIs!
export compute_metacells_total_UMIs!
export compute_metacells_types_by_cells!

using Base.Threads
using DataAxesFormats
using TanayLabUtilities
using Random
using Statistics

using ..Defaults
using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_linear_covered_metacell_cross_entropy_vector
import Metacells.Contracts.cell_linear_covered_metacell_kl_divergence_vector
import Metacells.Contracts.cell_linear_metacell_cross_entropy_vector
import Metacells.Contracts.cell_linear_metacell_kl_divergence_vector
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.cell_type_vector
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_scaled_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_scaled_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_mean_cells_linear_covered_cross_entropy_vector
import Metacells.Contracts.metacell_mean_cells_linear_covered_kl_divergence_vector
import Metacells.Contracts.metacell_mean_cells_linear_cross_entropy_vector
import Metacells.Contracts.metacell_mean_cells_linear_kl_divergence_vector
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.metacell_scaled_covered_UMIs_vector
import Metacells.Contracts.metacell_scaled_total_UMIs_vector
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
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "linear", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_scaled_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each gene in each metacell, scaled by divergence. We apply this scaling to
reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_scaled_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_scaled_fractions(daf; qualifier = "linear", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_log_scaled_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each gene in each metacell, scaled by divergence. We apply this scaling to
reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_scaled_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_scaled_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "scaled_linear", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_linear_covered_fractions!(
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
        metacell_gene_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_linear_fractions(daf; qualifier = "linear_covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_log_linear_covered_fractions!(
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
        metacell_gene_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "linear_covered", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_scaled_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each covered gene in each metacell, scaled by divergence. We apply this
scaling to reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        metacell_gene_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_scaled_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_scaled_fractions(daf; qualifier = "linear_covered", overwrite)
    return nothing
end

"""
    function compute_metacells_genes_log_scaled_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each covered gene in each metacell, scaled by divergence. We apply this
scaling to reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_scaled_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_log_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_log_scaled_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    do_compute_metacells_genes_log_fractions(
        daf;
        gene_fraction_regularization,
        qualifier = "scaled_linear_covered",
        overwrite,
    )
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
    function compute_metacells_scaled_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of UMIs used to estimate the fraction of all the genes in each metacell, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        gene_divergence_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_scaled_total_UMIs_vector(GuaranteedOutput),
    ],
) function compute_metacells_scaled_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_gene = daf["/ gene &! is_excluded : divergence"].array
    do_compute_metacells_UMIs(
        daf;
        qualifier = "scaled_total",
        genes_mask = "!is_excluded",
        scale_per_gene = 1.0 .- divergence_per_gene,
        overwrite,
    )
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
    function compute_metacells_scaled_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of covered genes per metacell, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_scaled_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_metacells_scaled_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_gene = daf["/ gene & is_covered : divergence"].array
    do_compute_metacells_UMIs(
        daf;
        qualifier = "scaled_covered",
        genes_mask = "is_covered",
        scale_per_gene = 1.0 .- divergence_per_gene,
        overwrite,
    )
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

function do_compute_metacells_genes_scaled_fractions(
    daf::DafWriter;
    qualifier::AbstractString,
    overwrite::Bool,
)::Nothing
    name = "$(qualifier)_fraction"

    fraction_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", name).array
    divergence_per_gene = get_vector(daf, "gene", "divergence").array

    scaled_fraction_per_gene_per_metacell = fraction_per_gene_per_metacell .* (1.0 .- divergence_per_gene)

    set_matrix!(
        daf,
        "gene",
        "metacell",
        "scaled_$(name)",
        bestify(scaled_fraction_per_gene_per_metacell; eltype = Float32);
        overwrite,
    )
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
    compute_cells_linear_metacell_cross_entropy!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
        overwrite::Bool = false,
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_linear_metacell_cross_entropy_vector(GuaranteedOutput),
    ],
) function compute_cells_linear_metacell_cross_entropy!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_cells_linear_metacell_cross_entropy(
        daf;
        gene_fraction_regularization,
        genes_mask = "!is_excluded",
        qualifier = "",
        overwrite,
    )
    return nothing
end

"""
    compute_cells_linear_covered_metacell_cross_entropy!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION
        overwrite::Bool = false,
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        cell_linear_covered_metacell_cross_entropy_vector(GuaranteedOutput),
    ],
) function compute_cells_linear_covered_metacell_cross_entropy!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_cells_linear_metacell_cross_entropy(
        daf;
        gene_fraction_regularization,
        genes_mask = "is_covered",
        qualifier = "covered",
        overwrite,
    )
    return nothing
end

function do_compute_cells_linear_metacell_cross_entropy(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat,
    genes_mask::AbstractString,
    qualifier::AbstractString,
    overwrite::Bool,
)::Nothing
    n_cells = axis_length(daf, "cell")
    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")

    linear_metacell_cross_entropy_per_cell = zeros(Float32, n_cells)

    @threads :greedy for metacell_index in 1:n_metacells
        metacell_name = name_per_metacell[metacell_index]

        UMIs_per_masked_gene_per_metacell_cell =
            daf["/ gene &$(genes_mask) / cell & metacell = $(metacell_name) : UMIs"].array
        n_masked_genes, n_metacell_cells = size(UMIs_per_masked_gene_per_metacell_cell)

        total_UMIs_per_metacell_cell = vec(sum(UMIs_per_masked_gene_per_metacell_cell; dims = 1))
        @assert_vector(total_UMIs_per_metacell_cell, n_metacell_cells)

        total_UMIs_in_metacell = sum(total_UMIs_per_metacell_cell)

        total_UMIs_per_masked_gene = vec(sum(UMIs_per_masked_gene_per_metacell_cell; dims = 2))
        @assert_vector(total_UMIs_per_masked_gene, n_masked_genes)

        fraction_per_masked_gene_per_metacell_cell =
            UMIs_per_masked_gene_per_metacell_cell ./ transpose(total_UMIs_per_metacell_cell)
        @assert_matrix(fraction_per_masked_gene_per_metacell_cell, n_masked_genes, n_metacell_cells)

        UMIs_per_masked_gene_per_other_cell = total_UMIs_per_masked_gene .- UMIs_per_masked_gene_per_metacell_cell  # NOJET
        @assert_matrix(UMIs_per_masked_gene_per_other_cell, n_masked_genes, n_metacell_cells)

        total_UMIs_per_other_cell = total_UMIs_in_metacell .- total_UMIs_per_metacell_cell  # NOJET # TODOX

        fraction_per_masked_gene_per_other_cell =
            UMIs_per_masked_gene_per_other_cell ./ transpose(total_UMIs_per_other_cell)
        @assert_matrix(fraction_per_masked_gene_per_other_cell, n_masked_genes, n_metacell_cells)

        cross_entropy_per_metacell_cell = vec(
            .- sum(fraction_per_masked_gene_per_metacell_cell .* log2.(fraction_per_masked_gene_per_other_cell .+ gene_fraction_regularization); dims = 1),
        )
        @assert_vector(cross_entropy_per_metacell_cell, n_metacell_cells)

        indices_of_metacell_cells = daf["/ cell & metacell = $(metacell_name) : index"].array
        linear_metacell_cross_entropy_per_cell[indices_of_metacell_cells] .= cross_entropy_per_metacell_cell
    end

    if qualifier == ""
        name = "linear_metacell_cross_entropy"
        @debug "Mean cells linear cross_entropy: $(mean(linear_metacell_cross_entropy_per_cell[linear_metacell_cross_entropy_per_cell .!= 0]))"
    else
        name = "linear_$(qualifier)_metacell_cross_entropy"
        @debug "Mean cells linear $(qualifier) cross_entropy: $(mean(linear_metacell_cross_entropy_per_cell[linear_metacell_cross_entropy_per_cell .!= 0]))"
    end
    set_vector!(daf, "cell", name, linear_metacell_cross_entropy_per_cell; overwrite)
    return nothing
end

"""
    compute_metacells_mean_cells_linear_cross_entropy!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_linear_metacell_cross_entropy_vector(RequiredInput),
        metacell_mean_cells_linear_cross_entropy_vector(GuaranteedOutput),
    ],
) function compute_metacells_mean_cells_linear_cross_entropy!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_cells_linear_cross_entropy_per_metacell =
        daf["/ cell : linear_metacell_cross_entropy @ metacell ! %> Mean"].array
    set_vector!(
        daf,
        "metacell",
        "mean_cells_linear_cross_entropy",
        mean_cells_linear_cross_entropy_per_metacell;
        overwrite,
    )
    @debug "Mean metacells linear cross_entropy: $(mean(mean_cells_linear_cross_entropy_per_metacell))"
    return nothing
end

"""
    compute_metacells_mean_cells_linear_covered_cross_entropy!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_linear_covered_metacell_cross_entropy_vector(RequiredInput),
        metacell_mean_cells_linear_covered_cross_entropy_vector(GuaranteedOutput),
    ],
) function compute_metacells_mean_cells_linear_covered_cross_entropy!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_cells_linear_covered_cross_entropy_per_metacell =
        daf["/ cell : linear_covered_metacell_cross_entropy @ metacell ! %> Mean"].array
    set_vector!(
        daf,
        "metacell",
        "mean_cells_linear_covered_cross_entropy",
        mean_cells_linear_covered_cross_entropy_per_metacell;
        overwrite,
    )
    @debug "Mean metacells linear covered cross_entropy: $(mean(mean_cells_linear_covered_cross_entropy_per_metacell))"
    return nothing
end

"""
    compute_cells_linear_metacell_kl_divergence!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
        overwrite::Bool = false,
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_linear_metacell_kl_divergence_vector(GuaranteedOutput),
    ],
) function compute_cells_linear_metacell_kl_divergence!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_cells_linear_metacell_kl_divergence(
        daf;
        gene_fraction_regularization,
        genes_mask = "!is_excluded",
        qualifier = "",
        overwrite,
    )
    return nothing
end

"""
    compute_cells_linear_covered_metacell_kl_divergence!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION
        overwrite::Bool = false,
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        cell_linear_covered_metacell_kl_divergence_vector(GuaranteedOutput),
    ],
) function compute_cells_linear_covered_metacell_kl_divergence!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_cells_linear_metacell_kl_divergence(
        daf;
        gene_fraction_regularization,
        genes_mask = "is_covered",
        qualifier = "covered",
        overwrite,
    )
    return nothing
end

function do_compute_cells_linear_metacell_kl_divergence(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat,
    genes_mask::AbstractString,
    qualifier::AbstractString,
    overwrite::Bool,
)::Nothing
    n_cells = axis_length(daf, "cell")
    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")

    linear_metacell_kl_divergence_per_cell = zeros(Float32, n_cells)

    @threads :greedy for metacell_index in 1:n_metacells
        metacell_name = name_per_metacell[metacell_index]

        UMIs_per_masked_gene_per_metacell_cell =
            daf["/ gene &$(genes_mask) / cell & metacell = $(metacell_name) : UMIs"].array
        n_masked_genes, n_metacell_cells = size(UMIs_per_masked_gene_per_metacell_cell)

        total_UMIs_per_metacell_cell = vec(sum(UMIs_per_masked_gene_per_metacell_cell; dims = 1))
        @assert_vector(total_UMIs_per_metacell_cell, n_metacell_cells)

        total_UMIs_in_metacell = sum(total_UMIs_per_metacell_cell)

        total_UMIs_per_masked_gene = vec(sum(UMIs_per_masked_gene_per_metacell_cell; dims = 2))
        @assert_vector(total_UMIs_per_masked_gene, n_masked_genes)

        fraction_per_masked_gene_per_metacell_cell =
            UMIs_per_masked_gene_per_metacell_cell ./ transpose(total_UMIs_per_metacell_cell)
        @assert_matrix(fraction_per_masked_gene_per_metacell_cell, n_masked_genes, n_metacell_cells)

        UMIs_per_masked_gene_per_other_cell = total_UMIs_per_masked_gene .- UMIs_per_masked_gene_per_metacell_cell
        @assert_matrix(UMIs_per_masked_gene_per_other_cell, n_masked_genes, n_metacell_cells)

        total_UMIs_per_other_cell = total_UMIs_in_metacell .- total_UMIs_per_metacell_cell  # NOJET # TODOX

        fraction_per_masked_gene_per_other_cell =
            UMIs_per_masked_gene_per_other_cell ./ transpose(total_UMIs_per_other_cell)
        @assert_matrix(fraction_per_masked_gene_per_other_cell, n_masked_genes, n_metacell_cells)

        kl_divergence_per_metacell_cell = vec(
            sum(
                fraction_per_masked_gene_per_metacell_cell .* (
                    log2.(fraction_per_masked_gene_per_metacell_cell .+ gene_fraction_regularization) .-
                    log2.(fraction_per_masked_gene_per_other_cell .+ gene_fraction_regularization)
                );
                dims = 1,
            ),
        )
        @assert_vector(kl_divergence_per_metacell_cell, n_metacell_cells)

        indices_of_metacell_cells = daf["/ cell & metacell = $(metacell_name) : index"].array
        linear_metacell_kl_divergence_per_cell[indices_of_metacell_cells] .= kl_divergence_per_metacell_cell
    end

    if qualifier == ""
        name = "linear_metacell_kl_divergence"
        @debug "Mean cells linear kl_divergence: $(mean(linear_metacell_kl_divergence_per_cell[linear_metacell_kl_divergence_per_cell .!= 0]))"
    else
        name = "linear_$(qualifier)_metacell_kl_divergence"
        @debug "Mean cells linear $(qualifier) kl_divergence: $(mean(linear_metacell_kl_divergence_per_cell[linear_metacell_kl_divergence_per_cell .!= 0]))"
    end
    set_vector!(daf, "cell", name, linear_metacell_kl_divergence_per_cell; overwrite)
    return nothing
end

"""
    compute_metacells_mean_cells_linear_kl_divergence!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_linear_metacell_kl_divergence_vector(RequiredInput),
        metacell_mean_cells_linear_kl_divergence_vector(GuaranteedOutput),
    ],
) function compute_metacells_mean_cells_linear_kl_divergence!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_cells_linear_kl_divergence_per_metacell =
        daf["/ cell : linear_metacell_kl_divergence @ metacell ! %> Mean"].array
    set_vector!(
        daf,
        "metacell",
        "mean_cells_linear_kl_divergence",
        mean_cells_linear_kl_divergence_per_metacell;
        overwrite,
    )
    @debug "Mean metacells linear kl_divergence: $(mean(mean_cells_linear_kl_divergence_per_metacell))"
    return nothing
end

"""
    compute_metacells_mean_cells_linear_covered_kl_divergence!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_linear_covered_metacell_kl_divergence_vector(RequiredInput),
        metacell_mean_cells_linear_covered_kl_divergence_vector(GuaranteedOutput),
    ],
) function compute_metacells_mean_cells_linear_covered_kl_divergence!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_cells_linear_covered_kl_divergence_per_metacell =
        daf["/ cell : linear_covered_metacell_kl_divergence @ metacell ! %> Mean"].array
    set_vector!(
        daf,
        "metacell",
        "mean_cells_linear_covered_kl_divergence",
        mean_cells_linear_covered_kl_divergence_per_metacell;
        overwrite,
    )
    @debug "Mean metacells linear covered kl_divergence: $(mean(mean_cells_linear_covered_kl_divergence_per_metacell))"
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

end  # module
