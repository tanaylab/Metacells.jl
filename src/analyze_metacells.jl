"""
Do simple metacells analysis.
"""
module AnalyzeMetacells

export compute_cells_types_by_metacells!
export compute_metacells_cells_correlations!
export compute_metacells_cells_projected_correlations!
export compute_metacells_euclidean_distances!
export compute_metacells_genes_UMIs!
export compute_metacells_genes_correlations!
export compute_metacells_genes_geomean_fractions!
export compute_metacells_genes_linear_fractions!
export compute_metacells_genes_log_geomean_fractions!
export compute_metacells_genes_log_linear_fractions!
export compute_metacells_max_skeleton_fold_distances!
export compute_metacells_n_cells!
export compute_metacells_regulators_correlations!
export compute_metacells_total_UMIs!
export compute_metacells_types_by_cells!
export compute_projected_metacells_genes_correlations!

using Base.Threads
using DataAxesFormats
using Distances
using Distributions
using MultipleTesting
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
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.gene_is_regulator_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_gene_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_geomean_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_metacell_euclidean_skeleton_distance
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
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

An estimated linear fraction of the UMIs of each non-excluded gene in each metacell. This is just the total UMIs of the
gene in the metacell divided by the total UMIs of the metacell, which is the "best" estimate assuming multinomial
sampling noise. However, this is sensitive to a few cells with very high expression levels ("bursty" genes).

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
        index_per_metacell_cell = findall(metacell_per_cell .== metacell_name)
        @assert length(index_per_metacell_cell) > 0

        @views UMIs_per_gene_per_metacell_cell = UMIs_per_gene_per_cell[:, index_per_metacell_cell]
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

The total number of UMIs used to estimate the fraction of all the non-excluded genes in each metacell.

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
    compute_metacells_genes_geomean_fractions!(
        daf::DafWriter;
        UMIs_regularization::AbstractFloat = $(DEFAULT.UMIs_regularization),
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given an assignment of cells to metacell, compute an geomean estimation of the fraction of UMIs of each non-excluded
gene for each metacell.

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

    index_per_included_gene = daf["/ gene &! is_excluded : index"].array
    n_included_genes = length(index_per_included_gene)

    geomean_fraction_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)

    progress_counter = Atomic{Int}(0)
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

        geomean_fraction_per_gene_per_metacell[index_per_included_gene, metacell_index] .=
            geomean_fraction_per_included_gene

        counter = atomic_add!(progress_counter, 1)
        print("\r$(progress_counter[]) ($(percent(counter + 1, n_metacells))) ...")
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

The log base 2 of the estimated geomean fraction of the UMIs of each non-excluded gene in each metacell. This adds the
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
    index_per_mask_gene = daf["/ gene &$(genes_mask) : index"].array
    linear_fraction_per_metacell_per_gene = zeros(Float32, n_metacells, n_genes)
    linear_fraction_per_metacell_per_gene[:, index_per_mask_gene] .=
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

The Euclidean distance between the log of the linear fraction of the skeleton genes between the metacells.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        metacell_metacell_euclidean_skeleton_distance(GuaranteedOutput),
    ],
) function compute_metacells_euclidean_distances!(daf::DafWriter; overwrite::Bool = false)::Nothing
    log_linear_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : log_linear_fraction"].array
    distances_between_metacells = pairwise(Euclidean(), log_linear_fraction_per_skeleton_per_metacell)  # NOJET
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

The maximal significant fold factor between the fraction of skeleton genes between the metacells. This uses heuristics
to require the fold factor be based on a sufficient number of UMIs to be robust.

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
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
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

    linear_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : linear_fraction"].array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    confidence_stds = quantile(Normal(), fold_confidence)

    confidence_linear_fractions_per_skeleton_per_metacells =  # NOJET
        confidence_stds .* sqrt.(linear_fraction_per_skeleton_per_metacell .* transpose(total_UMIs_per_metacell)) ./
        transpose(total_UMIs_per_metacell)

    low_log_linear_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            max.(
                linear_fraction_per_skeleton_per_metacell .- confidence_linear_fractions_per_skeleton_per_metacells,
                0.0,
            ) .+ gene_fraction_regularization,
        )

    high_log_linear_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            linear_fraction_per_skeleton_per_metacell .+ confidence_linear_fractions_per_skeleton_per_metacells .+
            gene_fraction_regularization,
        )

    UMIs_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : UMIs"].array

    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)
    distances_between_metacells[1, 1] = 0.0

    @threads :greedy for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views base_metacell_UMIs_per_skeleton = vec(UMIs_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_low_log_linear_fraction_per_skeleton =
            vec(low_log_linear_fraction_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_high_log_linear_fraction_per_skeleton =
            vec(high_log_linear_fraction_per_skeleton_per_metacell[:, base_metacell_index])

        n_other_metacells = base_metacell_index - 1
        index_per_other_metacell = 1:(base_metacell_index - 1)
        @views UMIs_per_skeleton_per_other_metacells = UMIs_per_skeleton_per_metacell[:, index_per_other_metacell]
        @views low_log_linear_fraction_per_skeleton_per_other_metacells =
            low_log_linear_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]
        @views high_log_linear_fraction_per_skeleton_per_other_metacells =
            high_log_linear_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]

        significant_fold_per_skeleton_per_other_metacell = confident_gene_distance.(
            min_significant_gene_UMIs,
            base_metacell_UMIs_per_skeleton,
            base_metacell_low_log_linear_fraction_per_skeleton,
            base_metacell_high_log_linear_fraction_per_skeleton,
            UMIs_per_skeleton_per_other_metacells,
            low_log_linear_fraction_per_skeleton_per_other_metacells,
            high_log_linear_fraction_per_skeleton_per_other_metacells,
        )
        @assert_matrix(significant_fold_per_skeleton_per_other_metacell, n_skeletons, n_other_metacells, Columns)

        distances_between_base_and_other_metacells =
            vec(maximum(significant_fold_per_skeleton_per_other_metacell; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, n_other_metacells)

        distances_between_metacells[index_per_other_metacell, base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, index_per_other_metacell] .=
            distances_between_base_and_other_metacells
    end

    set_matrix!(daf, "metacell", "metacell", "max_skeleton_fold_distance", distances_between_metacells; overwrite)
    return nothing
end

@inline function confident_gene_distance(  # UNTESTED
    min_significant_gene_UMIs::Integer,
    base_metacell_total_UMIs_of_gene::Integer,
    base_metacell_low_log_linear_fraction_of_gene::AbstractFloat,
    base_metacell_high_log_linear_fraction_of_gene::AbstractFloat,
    other_metacell_total_UMIs_of_gene::Integer,
    other_metacell_low_log_linear_fraction_of_gene::AbstractFloat,
    other_metacell_high_log_linear_fraction_of_gene::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = base_metacell_total_UMIs_of_gene + other_metacell_total_UMIs_of_gene
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = base_metacell_high_log_linear_fraction_of_gene < other_metacell_high_log_linear_fraction_of_gene

    highest_low_log_linear_fraction_of_gene =
        is_base_low * base_metacell_high_log_linear_fraction_of_gene +
        !is_base_low * other_metacell_high_log_linear_fraction_of_gene

    lowest_high_log_linear_fraction_of_gene =
        is_base_low * other_metacell_low_log_linear_fraction_of_gene +
        !is_base_low * base_metacell_low_log_linear_fraction_of_gene

    confident_gap = lowest_high_log_linear_fraction_of_gene - highest_low_log_linear_fraction_of_gene

    return (is_significant * max(confident_gap, 0.0))
end

"""
    function compute_metacells_genes_correlations!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = (DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the correlation between cells and metacells (marker) gene expression levels. TODOX.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        gene_is_excluded_vector(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        gene_correlation_between_cells_and_metacells_vector(GuaranteedOutput),
    ],
) function compute_metacells_genes_correlations!(  # UNTESTED
    daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    return adapter(
        daf;
        input_axes = ["cell" => "/ cell & metacell", "metacell" => "=", "gene" => "/ gene &! is_excluded"],
        input_data = [
            ("cell", "gene", "UMIs") => "=",
            ("cell", "total_UMIs") => "=",
            ("cell", "metacell") => "=",
            ("metacell", "gene", "UMIs") => "=",
            ("metacell", "total_UMIs") => "=",
        ],
        output_axes = ["gene" => "="],
        output_data = [("gene", "correlation_between_cells_and_metacells") => "="],
        empty = Dict(("gene", "correlation_between_cells_and_metacells") => 0.0),
        overwrite,
    ) do adapted
        n_genes = axis_length(adapted, "gene")

        UMIs_per_cell_per_gene = get_matrix(adapted, "cell", "gene", "UMIs").array
        UMIs_per_metacell_per_gene = get_matrix(adapted, "metacell", "gene", "UMIs").array
        total_UMIs_per_cell = get_vector(adapted, "cell", "total_UMIs").array
        total_UMIs_per_metacell = get_vector(adapted, "metacell", "total_UMIs").array
        metacell_index_per_cell = adapted["/ cell : metacell => index"].array
        total_punctuated_metacell_UMIs_per_cell =
            total_UMIs_per_metacell[metacell_index_per_cell] .- total_UMIs_per_cell

        progress_counter = Atomic{Int}(0)
        correlation_between_cells_and_metacells_per_gene = Vector{Float32}(undef, n_genes)
        @threads :greedy for gene_index in 1:n_genes
            @views UMIs_per_cell = UMIs_per_cell_per_gene[:, gene_index]
            @views metacell_UMIs_per_metacell = UMIs_per_metacell_per_gene[:, gene_index]
            cell_log_fraction_per_cell = log2.(UMIs_per_cell ./ total_UMIs_per_cell .+ gene_cell_fraction_regularization)
            punctuated_metacell_log_fraction_per_cell = log2.(
                (
                    (metacell_UMIs_per_metacell[metacell_index_per_cell] .- UMIs_per_cell) ./
                    total_punctuated_metacell_UMIs_per_cell
                ) .+ gene_cell_fraction_regularization,
            )
            correlation_between_cells_and_metacells_per_gene[gene_index] =
                cor(cell_log_fraction_per_cell, punctuated_metacell_log_fraction_per_cell)
            if isnan(correlation_between_cells_and_metacells_per_gene[gene_index])
                correlation_between_cells_and_metacells_per_gene[gene_index] = 0.0
            end
            counter = atomic_add!(progress_counter, 1)
            if counter % 100 == 0
                print("\r$(progress_counter[]) ($(percent(counter + 1, n_genes))) ...")
            end
        end

        set_vector!(
            adapted,
            "gene",
            "correlation_between_cells_and_metacells",
            correlation_between_cells_and_metacells_per_gene;
            overwrite,
        )
        is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
        is_pertinent_marker_per_gene = .!get_vector(daf, "gene", "is_lateral").array .& is_marker_per_gene
        @debug "Mean correlation of pertinent marker genes (between cells and their punctuated metacells): $(mean(correlation_between_cells_and_metacells_per_gene[is_pertinent_marker_per_gene]))"
        return nothing
    end
end

"""
    function compute_metacells_regulators_correlations!(
        daf::DafWriter;
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute covariance of regulator genes with other genes in each metacell. In each metacell, we downsample the
(non-excluded)) gene UMIs in each cell using `min_downsamples`, `min_downsamples` and `max_downsamples_quantile`, then
correlate the regulator genes with the rest of the genes. For each regulator gene, we look for the gene with the highest
absolute correlation, only if the total UMIs of both the anchor gene and correlated gene is at least
`min_significant_gene_UMIs`.

TODOX

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
    ],
    data = [
        gene_is_lateral_vector(RequiredInput),
        gene_is_excluded_vector(RequiredInput),
        gene_is_regulator_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_gene_is_environment_marker_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_gene_most_p_value_matrix(GuaranteedOutput),
        metacell_gene_most_q_value_matrix(GuaranteedOutput),
        metacell_gene_low_p_value_matrix(GuaranteedOutput),
        metacell_gene_low_q_value_matrix(GuaranteedOutput),
        metacell_gene_most_significant_correlated_gene_matrix(GuaranteedOutput),
        metacell_gene_most_significant_correlation_matrix(GuaranteedOutput),
    ],
) function compute_metacells_regulators_correlations!(  # UNTESTED
    daf::DafWriter;
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = 0.25,
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    min_significant_gene_UMIs::Integer = 20,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")

    index_per_included = daf["/ gene &! is_excluded : index"].array
    is_regulator_per_included = daf["/ gene &! is_excluded : is_regulator"].array
    is_lateral_per_included = daf["/ gene &! is_excluded : is_lateral"].array

    block_index_per_metacell = daf["/ metacell : block => index"].array
    is_environment_marker_per_included_gene_per_block =
        daf["/ gene &! is_excluded / block : is_environment_marker"].array .& .!is_lateral_per_included

    most_p_value_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)
    most_q_value_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)
    most_significant_correlation_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)
    most_significant_correlated_gene_per_gene_per_metacell = fill("", n_genes, n_metacells)
    low_p_value_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)
    low_q_value_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_metacells; rng) do metacell_index, rng
        metacell_name = name_per_metacell[metacell_index]
        total_UMIs_per_metacell_cell = daf["/ cell & metacell = $(metacell_name) : total_UMIs"].array

        UMIs_per_included_per_metacell_cell =
            daf["/ gene &! is_excluded / cell & metacell = $(metacell_name) : UMIs"].array

        downsamples_UMIs = downsamples(
            total_UMIs_per_metacell_cell;
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        is_significant_per_metacell_cell = total_UMIs_per_metacell_cell .>= downsamples_UMIs
        n_significant_metacell_cells = sum(is_significant_per_metacell_cell)

        @views UMIs_per_included_per_significant_metacell_cell =
            UMIs_per_included_per_metacell_cell[:, is_significant_per_metacell_cell]
        downsampled_UMIs_per_included_per_significant_metacell_cell =
            downsample(UMIs_per_included_per_significant_metacell_cell, downsamples_UMIs; dims = 2, rng)

        total_downsampled_UMIs_per_significant_metacell_cell =
            vec(sum(downsampled_UMIs_per_included_per_significant_metacell_cell; dims = 1))
        @assert_vector(total_downsampled_UMIs_per_significant_metacell_cell, n_significant_metacell_cells)
        @assert all(total_downsampled_UMIs_per_significant_metacell_cell .== downsamples_UMIs)

        block_index = block_index_per_metacell[metacell_index]
        @views is_environment_marker_per_included = is_environment_marker_per_included_gene_per_block[:, block_index]
        n_environment_markers = sum(is_environment_marker_per_included)

        @views downsampled_UMIs_per_environment_marker_per_significant_metacell_cell =
            downsampled_UMIs_per_included_per_significant_metacell_cell[is_environment_marker_per_included, :]

        total_downsampled_UMIs_per_environment_marker =
            vec(sum(downsampled_UMIs_per_environment_marker_per_significant_metacell_cell; dims = 2))
        @assert_vector(total_downsampled_UMIs_per_environment_marker, n_environment_markers)

        is_significant_per_environment_marker =
            total_downsampled_UMIs_per_environment_marker .>= min_significant_gene_UMIs

        @views is_regulator_per_environment_marker = is_regulator_per_included[is_environment_marker_per_included]
        @views index_per_environment_marker = index_per_included[is_environment_marker_per_included]
        @views index_per_significant_environment_marker =
            index_per_environment_marker[is_significant_per_environment_marker]
        n_significant_environment_markers = length(index_per_significant_environment_marker)

        is_significant_regulator_per_environment_marker =
            is_regulator_per_environment_marker .& is_significant_per_environment_marker

        index_per_significant_regulator_marker =
            index_per_environment_marker[findall(is_significant_regulator_per_environment_marker)]
        n_significant_regulator_markers = length(index_per_significant_regulator_marker)

        if n_significant_environment_markers == 1
            return nothing
        end

        @views downsampled_UMIs_per_significant_environment_marker_per_significant_metacell_cell =
            downsampled_UMIs_per_environment_marker_per_significant_metacell_cell[
                is_significant_per_environment_marker,
                :,
            ]
        @views downsampled_UMIs_per_significant_regulator_marker_per_significant_metacell_cell =
            downsampled_UMIs_per_environment_marker_per_significant_metacell_cell[
                is_significant_regulator_per_environment_marker,
                :,
            ]

        correlation_between_significant_and_regulator_genes, p_value_between_significant_and_regulator_genes =
            spearman_with_pvalue(
                flip(downsampled_UMIs_per_significant_environment_marker_per_significant_metacell_cell),
                flip(downsampled_UMIs_per_significant_regulator_marker_per_significant_metacell_cell),
            )
        @assert_matrix(
            correlation_between_significant_and_regulator_genes,
            n_significant_environment_markers,
            n_significant_regulator_markers,
            Columns
        )
        @assert_matrix(
            p_value_between_significant_and_regulator_genes,
            n_significant_environment_markers,
            n_significant_regulator_markers,
            Columns
        )

        q_value_between_significant_and_regulator_genes =
            Matrix{Float32}(undef, n_significant_environment_markers, n_significant_regulator_markers)
        is_other_per_significant_environment_marker = fill(true, n_significant_environment_markers)
        for regulator_gene_position in 1:n_significant_regulator_markers
            regulator_gene_index = index_per_significant_regulator_marker[regulator_gene_position]
            significant_gene_position = findfirst(index_per_significant_environment_marker .== regulator_gene_index)
            is_other_per_significant_environment_marker[significant_gene_position] = false

            if isapprox(
                p_value_between_significant_and_regulator_genes[significant_gene_position, regulator_gene_position],
                0,
            )
                p_value_between_significant_and_regulator_genes[significant_gene_position, regulator_gene_position] = 1
            end
            @assert isapprox(
                correlation_between_significant_and_regulator_genes[significant_gene_position, regulator_gene_position],
                1,
            )

            @views q_value_between_significant_and_regulator_gene =
                q_value_between_significant_and_regulator_genes[:, regulator_gene_position]

            q_value_between_significant_and_regulator_gene[significant_gene_position] = 1

            p_value_between_other_significant_and_regulator_genes = p_value_between_significant_and_regulator_genes[
                is_other_per_significant_environment_marker,
                regulator_gene_position,
            ]
            q_value_between_significant_and_regulator_gene[is_other_per_significant_environment_marker] .=
                q_values_of_p_values(p_value_between_other_significant_and_regulator_genes)

            is_other_per_significant_environment_marker[significant_gene_position] = true
        end

        most_significant_position_per_regulator = argmin(q_value_between_significant_and_regulator_genes; dims = 1)
        for regulator_most_significant_position in most_significant_position_per_regulator
            significant_gene_position, regulator_gene_position = regulator_most_significant_position.I

            significant_gene_index = index_per_significant_environment_marker[significant_gene_position]
            regulator_gene_index = index_per_significant_regulator_marker[regulator_gene_position]
            if significant_gene_index == regulator_gene_index
                @assert isapprox(
                    q_value_between_significant_and_regulator_genes[regulator_most_significant_position],
                    1,
                )
                continue
            end

            most_p_value_per_gene_per_metacell[regulator_gene_index, metacell_index] =
                p_value_between_significant_and_regulator_genes[regulator_most_significant_position]
            most_q_value_per_gene_per_metacell[regulator_gene_index, metacell_index] =
                q_value_between_significant_and_regulator_genes[regulator_most_significant_position]
            most_significant_correlation_per_gene_per_metacell[regulator_gene_index, metacell_index] =
                correlation_between_significant_and_regulator_genes[regulator_most_significant_position]
            most_significant_correlated_gene_per_gene_per_metacell[regulator_gene_index, metacell_index] =
                name_per_gene[significant_gene_index]
        end

        for regulator_gene_position in 1:n_significant_regulator_markers
            regulator_gene_index = index_per_significant_regulator_marker[regulator_gene_position]
            @views p_value_with_regulator = p_value_between_significant_and_regulator_genes[:, regulator_gene_position]
            @views q_value_with_regulator = q_value_between_significant_and_regulator_genes[:, regulator_gene_position]
            if n_significant_regulator_markers > 10
                low_p_value = quantile(p_value_with_regulator, 10 / n_significant_environment_markers) # TODOX
                low_q_value = quantile(q_value_with_regulator, 10 / n_significant_environment_markers) # TODOX
            else
                low_p_value = maximum(p_value_with_regulator) # TODOX
                low_q_value = maximum(q_value_with_regulator) # TODOX
            end
            low_p_value_per_gene_per_metacell[regulator_gene_index, metacell_index] = low_p_value
            low_q_value_per_gene_per_metacell[regulator_gene_index, metacell_index] = low_q_value
        end

        # @assert false  # TODOX
        counter = atomic_add!(progress_counter, 1)
        print("\r$(progress_counter[]) ($(percent(counter + 1, n_metacells))) ...")
        return nothing
    end

    set_matrix!(daf, "gene", "metacell", "most_p_value", sparsify(most_p_value_per_gene_per_metacell); overwrite)
    set_matrix!(daf, "gene", "metacell", "most_q_value", sparsify(most_q_value_per_gene_per_metacell); overwrite)
    set_matrix!(daf, "gene", "metacell", "low_p_value", sparsify(low_p_value_per_gene_per_metacell); overwrite)
    set_matrix!(daf, "gene", "metacell", "low_q_value", sparsify(low_q_value_per_gene_per_metacell); overwrite)

    set_matrix!(
        daf,
        "gene",
        "metacell",
        "most_significant_correlation",
        most_significant_correlation_per_gene_per_metacell;
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "metacell",
        "gene.most_significant_correlated",
        most_significant_correlated_gene_per_gene_per_metacell;
        overwrite,
    )

    return nothing
end

function spearman_with_pvalue(
    left_matrix::AbstractMatrix{<:Real},
    right_matrix::AbstractMatrix{<:Real},
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}}
    n_left_columns = size(left_matrix, 2)
    n_right_columns = size(right_matrix, 2)

    correlations = Matrix{Float64}(undef, n_left_columns, n_right_columns)
    p_values = Matrix{Float64}(undef, n_left_columns, n_right_columns)

    for left_column in 1:n_left_columns
        for right_column in 1:n_right_columns
            correlations[left_column, right_column], p_values[left_column, right_column] =
                spearman_with_pvalue(left_matrix[:, left_column], right_matrix[:, right_column])
        end
    end

    return correlations, p_values
end

function spearman_with_pvalue(
    left_values::AbstractVector{<:Real},
    right_values::AbstractVector{<:Real},
)::Tuple{AbstractFloat, AbstractFloat}
    @assert length(left_values) == length(right_values)
    n_values = length(left_values)
    if n_values <= 2
        return (1.0, 1.0)
    end

    # Compute tied ranks
    left_ranks = tiedrank(left_values)
    right_ranks = tiedrank(right_values)

    # Compute Spearman correlation
    correlation = cor(left_ranks, right_ranks)

    # Compute p-value using t-distribution approximation
    # t = correlation * sqrt((n_values-2)/(1-correlation^2))
    if abs(correlation) == 1.0
        p_value = 0.0
    elseif isnan(correlation)
        correlation = 1.0
        p_value = 1.0
    else
        t_statistic = correlation * sqrt((n_values - 2) / (1 - correlation^2))
        t_distribution = TDist(n_values - 2)
        p_value = 2 * (1 - cdf(t_distribution, abs(t_statistic)))
    end

    return (correlation, p_value)
end

function q_values_of_p_values(p_values::AbstractVector{<:AbstractFloat})::AbstractVector{<:AbstractFloat}
    if !(p_values isa Vector)
        p_values = Vector(p_values)
    end
    return adjust(p_values, BenjaminiHochberg())
end

"""
    function compute_metacells_cells_correlations!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the correlation between cells and metacells (marker) gene expression levels. TODOX.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        gene_is_marker_vector(RequiredInput),
        # gene_is_regulator_vector(RequiredInput),  TODOX
        gene_is_lateral_vector(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        # cell_genes_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # cell_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        cell_pertinent_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # cell_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # cell_pertinent_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # metacell_mean_cells_genes_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # metacell_mean_cells_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        metacell_mean_cells_pertinent_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # metacell_mean_cells_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        # metacell_mean_cells_pertinent_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
    ],
) function compute_metacells_cells_correlations!(  # UNTESTED
    daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(daf, "cell")

    correlation_per_cell = [
        zeros(Float32, n_cells),
        zeros(Float32, n_cells),
        zeros(Float32, n_cells),
        zeros(Float32, n_cells),
        zeros(Float32, n_cells),
    ]

    total_UMIs_per_cell = daf["/ cell : total_UMIs"].array
    total_UMIs_per_metacell = daf["/ metacell : total_UMIs"].array
    metacell_index_per_cell = daf["/ cell : metacell ?? 0 => index"].array
    UMIs_per_cell_per_gene = daf["/ cell / gene : UMIs"].array
    UMIs_per_gene_per_metacell = daf["/ gene / metacell : UMIs"].array

    is_included_per_gene = .!get_vector(daf, "gene", "is_excluded").array
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array .& is_included_per_gene
    is_pertinent_marker_per_gene = .!get_vector(daf, "gene", "is_lateral").array .& is_marker_per_gene
    # is_regulator_per_gene = get_vector(daf, "gene", "is_regulator").array .& is_included_per_gene
    # is_pertinent_regulator_per_gene = .!get_vector(daf, "gene", "is_lateral").array .& is_regulator_per_gene

    progress_counter = Atomic{Int}(0)
    @threads for cell_index in 1:n_cells
        metacell_index = metacell_index_per_cell[cell_index]
        if metacell_index > 0
            total_cell_UMIs = total_UMIs_per_cell[cell_index]
            total_metacell_UMIs = total_UMIs_per_metacell[metacell_index]
            total_punctuated_metacell_UMIs = total_metacell_UMIs - total_cell_UMIs

            for (correlation_index, genes_mask) in enumerate((
                # is_included_per_gene,
                # is_marker_per_gene,
                is_pertinent_marker_per_gene,
                # is_regulator_per_gene,
                # is_pertinent_regulator_per_gene,
            ))
                cell_UMIs_per_masked_gene = UMIs_per_cell_per_gene[cell_index, genes_mask]
                metacell_UMIs_per_masked_gene = UMIs_per_gene_per_metacell[genes_mask, metacell_index]

                cell_log_fraction_per_masked_gene =
                    log2.(cell_UMIs_per_masked_gene ./ total_cell_UMIs .+ gene_cell_fraction_regularization)
                punctuated_metacell_log_fraction_per_masked_gene = log2.(
                    ((metacell_UMIs_per_masked_gene .- cell_UMIs_per_masked_gene) ./ total_punctuated_metacell_UMIs) .+ gene_cell_fraction_regularization,
                )

                correlation = cor(cell_log_fraction_per_masked_gene, punctuated_metacell_log_fraction_per_masked_gene)
                if isnan(correlation)
                    correlation = 0.0
                end
                correlation_per_cell[correlation_index][cell_index] = correlation
            end
        end
        counter = atomic_add!(progress_counter, 1)
        if counter % 100 == 0
            print("\r$(progress_counter[]) ($(percent(counter + 1, n_cells))) ...")
        end
    end
    println()

    is_in_metacell_per_cell = metacell_index_per_cell .> 0
    for (correlation_index, (cells_property, title)) in enumerate((
        #"genes_correlation_with_punctuated_metacells",
        #"markers_correlation_with_punctuated_metacells",
        ("pertinent_markers_correlation_with_punctuated_metacells", "pertinent markers"),
        #"regulators_correlation_with_punctuated_metacells",
        #"pertinent_regulators_correlation_with_punctuated_metacells",
    ))
        @debug "Mean correlation of cells with their punctuated metacells (of $(title)): $(mean(correlation_per_cell[correlation_index][is_in_metacell_per_cell]))"
        set_vector!(daf, "cell", cells_property, correlation_per_cell[correlation_index]; overwrite)

        correlation_per_metacell = daf["/ cell : $(cells_property) @ metacell ! %> Mean"].array
        @debug "Mean correlation of punctuated metacells and their cells (of $(title)): $(mean(correlation_per_metacell))"
        set_vector!(daf, "metacell", "mean_cells_$(cells_property)", correlation_per_metacell; overwrite)
    end

    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
    ],
) Contract(;
    axes = [gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_projected_metacell_vector(RequiredInput),
        cell_metacell_vector(OptionalInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        cell_pertinent_markers_correlation_with_projected_metacells_vector(GuaranteedOutput),
    ],
) function compute_metacells_cells_projected_correlations!(  # UNTESTED
    modules_daf::DafReader,
    cells_daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(cells_daf, "cell")

    correlation_per_cell = Vector{Float32}(undef, n_cells)

    total_UMIs_per_cell = cells_daf["/ cell : total_UMIs"].array
    metacell_per_cell = cells_daf["/ cell : metacell.projected"].array
    UMIs_per_cell_per_gene = cells_daf["/ cell / gene : UMIs"].array

    metacell_index_per_cell = axis_indices(modules_daf, "metacell", metacell_per_cell)
    total_UMIs_per_metacell = modules_daf["/ metacell : total_UMIs"].array
    UMIs_per_gene_per_metacell = modules_daf["/ gene / metacell : UMIs"].array

    pertinent_genes = modules_daf["/ gene & is_marker &! is_lateral"]
    indices_of_pertinent_genes = axis_indices(cells_daf, "gene", pertinent_genes)

    progress_counter = Atomic{Int}(0)
    @threads for cell_index in 1:n_cells
        total_cell_UMIs = total_UMIs_per_cell[cell_index]
        cell_UMIs_per_pertinent_marker = UMIs_per_cell_per_gene[cell_index, indices_of_pertinent_genes]
        cell_log_fraction_per_pertinent_marker =
            log2.(cell_UMIs_per_pertinent_marker ./ total_cell_UMIs .+ gene_cell_fraction_regularization)

        metacell_index = metacell_index_per_cell[cell_index]
        total_metacell_UMIs = total_UMIs_per_metacell[metacell_index]
        metacell_UMIs_per_pertinent_marker = UMIs_per_gene_per_metacell[indices_of_pertinent_genes, metacell_index]

        metacell_log_fraction_per_pertinent_marker =
            log2.(metacell_UMIs_per_pertinent_marker ./ total_metacell_UMIs .+ gene_cell_fraction_regularization)

        correlation = cor(cell_log_fraction_per_pertinent_marker, metacell_log_fraction_per_pertinent_marker)
        if isnan(correlation)
            correlation = 0.0
        end
        correlation_per_cell[cell_index] = correlation

        counter = atomic_add!(progress_counter, 1)
        if counter % 100 == 0
            print("\r$(progress_counter[]) ($(percent(counter + 1, n_cells))) ...")
        end
    end

    set_vector!(cells_daf, "cell", "pertinent_markers_correlation_with_projected_metacells", correlation_per_cell; overwrite)
    @debug "Mean correlation of cells and their projected metacells (of pertinent markers): $(mean(correlation_per_cell))"

    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
    ],
) Contract(
    axes = [gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_projected_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        gene_correlation_between_cells_and_projected_metacells_vector(GuaranteedOutput),
    ],
) function compute_projected_metacells_genes_correlations!(  # UNTESTED
    modules_daf::DafReader,
    cells_daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    n_genes_in_cells = axis_length(cells_daf, "gene")
    cells_genes = Set(axis_vector(cells_daf, "gene"))

    included_genes = modules_daf["/ gene &! is_excluded"]
    included_genes = filter(included_genes) do gene
        return gene in cells_genes
    end

    indices_in_modules_per_included_gene = axis_indices(modules_daf, "gene", included_genes)
    indices_in_cells_per_included_gene = axis_indices(cells_daf, "gene", included_genes)
    n_included_genes = length(indices_in_cells_per_included_gene)

    UMIs_per_cell_per_gene = get_matrix(cells_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(cells_daf, "cell", "total_UMIs").array

    UMIs_per_metacell_per_gene = get_matrix(modules_daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_metacell = get_vector(modules_daf, "metacell", "total_UMIs").array

    metacell_per_cell = cells_daf["/ cell : metacell.projected"].array
    metacell_index_per_cell = axis_indices(modules_daf, "metacell", metacell_per_cell)

    progress_counter = Atomic{Int}(0)
    correlation_between_cells_and_projected_metacells_per_gene = zeros(Float32, n_genes_in_cells)

    @threads :greedy for included_gene_position in 1:n_included_genes
        gene_index_in_cells = indices_in_cells_per_included_gene[included_gene_position]
        @views UMIs_per_cell = UMIs_per_cell_per_gene[:, gene_index_in_cells]
        cell_log_fraction_per_cell = log2.(UMIs_per_cell ./ total_UMIs_per_cell .+ gene_cell_fraction_regularization)

        gene_index_in_modules = indices_in_cells_per_included_gene[included_gene_position]
        @views metacell_UMIs_per_metacell = UMIs_per_metacell_per_gene[:, gene_index_in_modules]
        metacell_log_fraction_per_cell = log2.(
            metacell_UMIs_per_metacell[metacell_index_per_cell] ./ total_UMIs_per_metacell[metacell_index_per_cell] .+
            gene_cell_fraction_regularization,
        )

        correlation_between_cells_and_projected_metacells_per_gene[gene_index_in_cells] =
            cor(cell_log_fraction_per_cell, metacell_log_fraction_per_cell)

        if isnan(correlation_between_cells_and_projected_metacells_per_gene[gene_index_in_cells])
            correlation_between_cells_and_projected_metacells_per_gene[gene_index_in_cells] = 0.0
        end
        counter = atomic_add!(progress_counter, 1)
        if counter % 100 == 0
            print("\r$(progress_counter[]) ($(percent(counter + 1, n_included_genes))) ...")
        end
    end

    set_vector!(
        cells_daf,
        "gene",
        "correlation_between_cells_and_projected_metacells",
        bestify(correlation_between_cells_and_projected_metacells_per_gene);
        overwrite,
    )

    pertinent_markers = Set(modules_daf["/ gene & is_marker &! is_lateral"])
    filter!(included_genes) do gene
        return gene in pertinent_markers
    end
    indices_in_cells_per_pertinent_marker = axis_indices(cells_daf, "gene", included_genes)
    @debug "Mean correlation of pertinent marker genes (between cells and their projected metacells): $(mean(correlation_between_cells_and_projected_metacells_per_gene[indices_in_cells_per_pertinent_marker]))"
    return nothing
end

end  # module
