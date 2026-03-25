"""
Do simple metacells analysis.
"""
module AnalyzeMetacells

export compute_matrix_of_UMIs_per_gene_per_metacell!
export compute_matrix_of_correlation_per_gene_per_gene_of_subset_of_metacells!
export compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!
export compute_matrix_of_linear_fraction_per_gene_per_metacell!
export compute_matrix_of_log_linear_fraction_per_gene_per_metacell!
export compute_matrix_of_max_skeleton_fold_distance_between_metacells!
export compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!
export compute_vector_of_n_cells_per_metacell!
export compute_vector_of_total_UMIs_per_metacell!
export compute_vector_of_type_per_cell_by_metacells!
export compute_vector_of_type_per_metacell_by_cells!

using Base.Threads
using DataAxesFormats
using Distances
using Distributions
using LoopVectorization
using MultipleTesting
using TanayLabUtilities
using Random
using StatsBase

using ..Defaults
using ..Contracts

import Base.Threads.maxthreadid
import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_correlation_between_cells_and_punctuated_metacells_per_gene
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell
import Metacells.Contracts.vector_of_type_per_cell
import Metacells.Contracts.vector_of_type_per_metacell

"""
    function compute_vector_of_n_cells_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_cells_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [vector_of_metacell_per_cell(RequiredInput), vector_of_n_cells_per_metacell(CreatedOutput)],
) function compute_vector_of_n_cells_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells_per_metacell = daf["@ cell / metacell =@ >> Count"].array
    set_vector!(daf, "metacell", "n_cells", n_cells_per_metacell; overwrite)
    @debug "Mean cells in metacell: $(mean(n_cells_per_metacell))" _group = :mcs_results
    return nothing
end

"""
    function compute_matrix_of_UMIs_per_gene_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_UMIs_per_gene_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_UMIs_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    UMIs_per_gene_per_metacell = daf["@ gene @ cell :: UMIs |/ metacell =@ >| Sum"].array
    set_matrix!(daf, "gene", "metacell", "UMIs", bestify(UMIs_per_gene_per_metacell); overwrite)  # NOJET
    return nothing
end

"""
    function compute_vector_of_total_UMIs_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_total_UMIs_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_total_UMIs_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_metacell = daf["@ metacell @ gene [ ! is_excluded ] :: UMIs >| Sum"].array
    set_vector!(daf, "metacell", "total_UMIs", total_UMIs_per_metacell; overwrite)
    @debug "Mean UMIs in metacell: $(mean(total_UMIs_per_metacell))" _group = :mcs_results
    return nothing
end

"""
    function compute_matrix_of_linear_fraction_per_gene_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_linear_fraction_per_gene_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_linear_fraction_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    indices_of_included_genes = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    UMIs_per_metacell_per_included_gene = daf["@ metacell @ gene [ ! is_excluded ] :: UMIs"].array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    linear_fraction_per_metacell_per_gene = zeros(Float32, n_metacells, n_genes)
    @. linear_fraction_per_metacell_per_gene[:, indices_of_included_genes] =
        UMIs_per_metacell_per_included_gene / total_UMIs_per_metacell

    set_matrix!(  # NOJET
        daf,
        "metacell",
        "gene",
        "linear_fraction",
        bestify(linear_fraction_per_metacell_per_gene);
        overwrite,
    )  # NOJET

    return nothing
end

"""
    function compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_log_linear_fraction_per_gene_per_metacell`](@ref). This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    fraction_per_metacell_per_gene =
        mutable_array(densify(get_matrix(daf, "metacell", "gene", "linear_fraction").array))
    empty_dense_matrix!(
        daf,
        "metacell",
        "gene",
        "log_linear_fraction",
        Float32;
        overwrite,
    ) do log_fraction_per_metacell_per_gene
        @assert LoopVectorization.check_args(log_fraction_per_metacell_per_gene) "check_args failed in compute_matrix_of_log_linear_fraction_per_gene_per_metacell!\nfor log_fraction_per_metacell_per_gene: $(brief(log_fraction_per_metacell_per_gene))"
        @assert LoopVectorization.check_args(fraction_per_metacell_per_gene) "check_args failed in compute_matrix_of_log_linear_fraction_per_gene_per_metacell!\nfor fraction_per_metacell_per_gene: $(brief(fraction_per_metacell_per_gene))"
        n_metacells, n_genes = size(log_fraction_per_metacell_per_gene)
        parallel_loop_wo_rng(  # NOJET
            1:n_genes;
            name = "log_linear_fraction_per_gene_per_metacell",
            progress = DebugProgress(n_genes; group = :mcs_loops, desc = "log_linear_fraction_per_gene_per_metacell"),
            progress_chunk = 100,
        ) do gene_index
            @turbo for metacell_index in 1:n_metacells
                log_fraction_per_metacell_per_gene[metacell_index, gene_index] =
                    log2(fraction_per_metacell_per_gene[metacell_index, gene_index] + gene_fraction_regularization)
            end
            return nothing
        end
        return nothing
    end
    return nothing
end

"""
    function compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_euclidean_skeleton_fold_distance_between_metacells`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(CreatedOutput),
    ],
) function compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    log_linear_fraction_per_metacell_per_skeleton =
        daf["@ metacell @ gene [ is_skeleton ] :: log_linear_fraction"].array
    log_linear_fraction_per_skeleton_per_metacell = flipped(log_linear_fraction_per_metacell_per_skeleton)
    n_metacells = axis_length(daf, "metacell")
    distances_between_metacells = parallel_pairwise(  # NOJET
        Euclidean(),
        log_linear_fraction_per_skeleton_per_metacell;
        dims = 2,
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "euclidean_skeleton_fold_distance_between_metacells",
        ),
    )
    set_matrix!(daf, "metacell", "metacell", "euclidean_skeleton_fold_distance", distances_between_metacells; overwrite)
    return nothing
end

"""
    function compute_matrix_of_max_skeleton_fold_distance_between_metacells!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_max_skeleton_fold_distance_between_metacells`](@ref). For computing this fold factor, we add
the `gene_fraction_regularization` to deal with zero fractions. We ignore "insignificant" genes whose total UMIs in the
compared metacells isn't at least `min_significant_gene_UMIs`. We also we reduce the distance using the
`fold_confidence` based on the number of UMIs used to estimate the expression in the metacells. That is, when we have a
high fold factor, we are pretty certain the metacells truly have a significant difference in the expression of the
skeleton gene.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_max_skeleton_fold_distance_between_metacells(CreatedOutput),
    ],
) function compute_matrix_of_max_skeleton_fold_distance_between_metacells!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_significant_gene_UMIs >= 0
    @assert 0 < fold_confidence < 1

    linear_fraction_per_metacell_per_skeleton = daf["@ metacell @ gene [ is_skeleton ] :: linear_fraction"].array
    linear_fraction_per_skeleton_per_metacell = Matrix(flip(linear_fraction_per_metacell_per_skeleton))
    total_UMIs_per_metacell = mutable_array(densify(get_vector(daf, "metacell", "total_UMIs").array))

    confidence_stds = quantile(Normal(), fold_confidence)

    confidence_linear_fractions_per_skeleton_per_metacells =  # NOJET
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell))
    @assert LoopVectorization.check_args(confidence_linear_fractions_per_skeleton_per_metacells) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (confidence)\nfor confidence_linear_fractions_per_skeleton_per_metacells: $(brief(confidence_linear_fractions_per_skeleton_per_metacells))"
    @assert LoopVectorization.check_args(linear_fraction_per_skeleton_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (confidence)\nfor linear_fraction_per_skeleton_per_metacell: $(brief(linear_fraction_per_skeleton_per_metacell))"
    @assert LoopVectorization.check_args(total_UMIs_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (confidence)\nfor total_UMIs_per_metacell: $(brief(total_UMIs_per_metacell))"
    n_skeletons, n_metacells = size(linear_fraction_per_skeleton_per_metacell)
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "confidence_linear_fractions_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "confidence_linear_fractions_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index] =
                confidence_stds * sqrt(
                    linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] *
                    total_UMIs_per_metacell[metacell_index],
                ) / total_UMIs_per_metacell[metacell_index]
        end
        return nothing
    end

    low_log_linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell)) # NOJET
    @assert LoopVectorization.check_args(low_log_linear_fraction_per_skeleton_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (low)\nfor low_log_linear_fraction_per_skeleton_per_metacell: $(brief(low_log_linear_fraction_per_skeleton_per_metacell))"
    @assert LoopVectorization.check_args(linear_fraction_per_skeleton_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (low)\nfor linear_fraction_per_skeleton_per_metacell: $(brief(linear_fraction_per_skeleton_per_metacell))"
    @assert LoopVectorization.check_args(confidence_linear_fractions_per_skeleton_per_metacells) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (low)\nfor confidence_linear_fractions_per_skeleton_per_metacells: $(brief(confidence_linear_fractions_per_skeleton_per_metacells))"
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "low_log_linear_fraction_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "low_log_linear_fraction_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            low_log_linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] = log2(
                max(
                    linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] -
                    confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index],
                    0.0,
                ) + gene_fraction_regularization,
            )
        end
        return nothing
    end

    high_log_linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell)) # NOJET
    @assert LoopVectorization.check_args(high_log_linear_fraction_per_skeleton_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (high)\nfor high_log_linear_fraction_per_skeleton_per_metacell: $(brief(high_log_linear_fraction_per_skeleton_per_metacell))"
    @assert LoopVectorization.check_args(linear_fraction_per_skeleton_per_metacell) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (high)\nfor linear_fraction_per_skeleton_per_metacell: $(brief(linear_fraction_per_skeleton_per_metacell))"
    @assert LoopVectorization.check_args(confidence_linear_fractions_per_skeleton_per_metacells) "check_args failed in compute_matrix_of_max_skeleton_fold_distance_between_metacells! (high)\nfor confidence_linear_fractions_per_skeleton_per_metacells: $(brief(confidence_linear_fractions_per_skeleton_per_metacells))"
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "high_log_linear_fraction_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "high_log_linear_fraction_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            high_log_linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] = log2(
                linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] +
                confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index] +
                gene_fraction_regularization,
            )
        end
        return nothing
    end

    UMIs_per_metacell_per_skeleton = daf["@ metacell @ gene [ is_skeleton ] :: UMIs"].array
    UMIs_per_skeleton_per_metacell = flip(UMIs_per_metacell_per_skeleton)

    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)
    distances_between_metacells[1, 1] = 0.0

    parallel_loop_wo_rng(
        reverse(2:(n_metacells));
        progress = DebugProgress(
            n_metacells - 1;
            group = :mcs_loops,
            desc = "max_skeleton_fold_distance_between_metacells",
        ),
    ) do base_metacell_index
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
        return nothing
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
    function compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = (DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_correlation_between_cells_and_punctuated_metacells_per_gene`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_correlation_between_cells_and_punctuated_metacells_per_gene(CreatedOutput),
    ],
) function compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    gene_fraction_regularization = Float32(gene_fraction_regularization)
    n_genes = axis_length(daf, "gene")

    index_per_included_gene = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    n_included_genes = length(index_per_included_gene)

    indices_of_grouped_cells = get_query(daf, "@ cell [ metacell ] : index").array
    n_grouped_cells = length(indices_of_grouped_cells)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array

    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    total_UMIs_per_grouped_cell = total_UMIs_per_cell[indices_of_grouped_cells]

    metacell_index_per_grouped_cell = daf["@ cell : metacell ?? : index"].array

    total_punctuated_metacell_UMIs_per_grouped_cell =
        total_UMIs_per_metacell[metacell_index_per_grouped_cell] .- total_UMIs_per_grouped_cell

    correlation_between_cells_and_punctuated_metacells_per_gene = Vector{Float32}(undef, n_genes)

    gene_cell_log_fraction_per_grouped_cell_per_thread =
        [Vector{Float32}(undef, n_grouped_cells) for _ in 1:maxthreadid()]
    gene_punctuated_metacell_log_fraction_per_grouped_cell_per_thread =
        [Vector{Float32}(undef, n_grouped_cells) for _ in 1:maxthreadid()]

    todox_all_zero = Atomic{Int32}(0)
    parallel_loop_wo_rng(
        1:n_included_genes;
        policy = :static,
        progress = DebugProgress(
            n_included_genes;
            group = :mcs_loops,
            desc = "correlation_between_cells_and_punctuated_metacells_per_gene",
        ),
        progress_chunk = 100,
    ) do included_gene_position
        gene_index = index_per_included_gene[included_gene_position]

        gene_cell_log_fraction_per_grouped_cell = gene_cell_log_fraction_per_grouped_cell_per_thread[threadid()]
        gene_punctuated_metacell_log_fraction_per_grouped_cell =
            gene_punctuated_metacell_log_fraction_per_grouped_cell_per_thread[threadid()]

        all_zero_grouped_cell_UMIs = true
        for grouped_cell_position in 1:n_grouped_cells
            if UMIs_per_cell_per_gene[indices_of_grouped_cells[grouped_cell_position], gene_index] > 0
                all_zero_grouped_cell_UMIs = false
                break
            end
        end
        if all_zero_grouped_cell_UMIs
            atomic_add!(todox_all_zero, Int32(1))
            correlation_between_cells_and_punctuated_metacells_per_gene[gene_index] = 0.0f0
            return nothing
        end

        all_zero_grouped_punctuated_metacell_UMIs = true
        for grouped_cell_position in 1:n_grouped_cells
            cell_index = indices_of_grouped_cells[grouped_cell_position]
            metacell_index = metacell_index_per_grouped_cell[grouped_cell_position]
            if UMIs_per_metacell_per_gene[metacell_index, gene_index] > UMIs_per_cell_per_gene[cell_index, gene_index]
                all_zero_grouped_punctuated_metacell_UMIs = false
                break
            end
        end
        if all_zero_grouped_punctuated_metacell_UMIs
            atomic_add!(todox_all_zero, Int32(1))
            correlation_between_cells_and_punctuated_metacells_per_gene[gene_index] = 0.0f0
            return nothing
        end

        for (grouped_cell_position, cell_index) in enumerate(indices_of_grouped_cells)
            metacell_index = metacell_index_per_grouped_cell[grouped_cell_position]
            gene_cell_UMIs = UMIs_per_cell_per_gene[cell_index, gene_index]
            gene_metacell_UMIs = UMIs_per_metacell_per_gene[metacell_index, gene_index]
            gene_cell_log_fraction_per_grouped_cell[grouped_cell_position] = gene_cell_UMIs
            gene_punctuated_metacell_log_fraction_per_grouped_cell[grouped_cell_position] =
                gene_metacell_UMIs - gene_cell_UMIs
        end
        @assert LoopVectorization.check_args(gene_cell_log_fraction_per_grouped_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!\nfor gene_cell_log_fraction_per_grouped_cell: $(brief(gene_cell_log_fraction_per_grouped_cell))"
        @assert LoopVectorization.check_args(total_UMIs_per_grouped_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!\nfor total_UMIs_per_grouped_cell: $(brief(total_UMIs_per_grouped_cell))"
        @assert LoopVectorization.check_args(gene_punctuated_metacell_log_fraction_per_grouped_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!\nfor gene_punctuated_metacell_log_fraction_per_grouped_cell: $(brief(gene_punctuated_metacell_log_fraction_per_grouped_cell))"
        @assert LoopVectorization.check_args(total_punctuated_metacell_UMIs_per_grouped_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!\nfor total_punctuated_metacell_UMIs_per_grouped_cell: $(brief(total_punctuated_metacell_UMIs_per_grouped_cell))"
        @turbo for grouped_cell_position in 1:n_grouped_cells
            gene_cell_log_fraction_per_grouped_cell[grouped_cell_position] = log2(
                gene_cell_log_fraction_per_grouped_cell[grouped_cell_position] /
                total_UMIs_per_grouped_cell[grouped_cell_position] + gene_fraction_regularization,
            )
            gene_punctuated_metacell_log_fraction_per_grouped_cell[grouped_cell_position] = log2(
                gene_punctuated_metacell_log_fraction_per_grouped_cell[grouped_cell_position] /
                total_punctuated_metacell_UMIs_per_grouped_cell[grouped_cell_position] +
                gene_fraction_regularization,
            )
        end

        correlation_between_cells_and_punctuated_metacells_per_gene[gene_index] = zero_cor_between_vectors(
            gene_cell_log_fraction_per_grouped_cell,
            gene_punctuated_metacell_log_fraction_per_grouped_cell,
        )

        return nothing
    end

    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(n_included_genes) ($(percent(todox_all_zero[], n_included_genes)))" _group =
        :todox
    set_vector!(
        daf,
        "gene",
        "correlation_between_cells_and_punctuated_metacells",
        correlation_between_cells_and_punctuated_metacells_per_gene;
        overwrite,
    )

    is_pertinent_marker_per_gene =
        get_vector(daf, "gene", "is_marker").array .& .!get_vector(daf, "gene", "is_lateral").array
    @debug (
        "Mean correlation of pertinent marker genes between cells and their punctuated metacells: " *
        "$(mean(correlation_between_cells_and_punctuated_metacells_per_gene[is_pertinent_marker_per_gene]))"
    ) _group = :mcs_results

    return nothing
end

"""
    function compute_vector_of_type_per_cell_by_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_cell`](@ref) by using the type of the metacell each cell belongs to. Cells that do
not belong to any metacells (including excluded cells) are given an empty type. This method is used when metacells are
annotated with types (which is typically a supervised process).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_metacell(RequiredInput),
        vector_of_type_per_cell(CreatedOutput),
    ],
) function compute_vector_of_type_per_cell_by_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_cell = daf["@ cell : metacell ?? '' : type"].array
    set_vector!(daf, "cell", "type", type_per_cell; overwrite)
    return nothing
end

"""
    function compute_vector_of_type_per_metacell_by_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_metacell`](@ref) by using the type of the cells grouped into each metacell. The
most frequent ("mode") of the type is used for each metacell. This method is used when cells are annotated with types
(typically when importing single-cell data with type annotations).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_cell(RequiredInput),
        vector_of_type_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_type_per_metacell_by_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_metacell = daf["@ cell : type / metacell =@ >> Mode"].array
    set_vector!(daf, "metacell", "type", type_per_metacell; overwrite)
    return nothing
end

"""
    compute_matrix_of_correlation_per_gene_per_gene_of_subset_of_metacells!(
        daf::DafWriter;
        metacells_subset::Union{QueryString, BitVector, AbstractVector{Bool}},
        matrix_name::AbstractString,
        overwrite::Bool = false,
    )::Nothing

Correlate the marker genes in an arbitrary `metacells_subset`. This subset can be either an explicit mask or a query
that masks the `metacell` axis. The result is stored in a per-gene-per-gene `matrix_name`.
"""
@logged :mc_ops @computation Contract(;
    is_relaxed = true,
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_marker_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
    ],
) function compute_matrix_of_correlation_per_gene_per_gene_of_subset_of_metacells!(
    daf::DafWriter;
    metacells_subset::Union{QueryString, BitVector, AbstractVector{Bool}},
    matrix_name::AbstractString,
    overwrite::Bool = false,
)::Nothing
    if metacells_subset isa QueryString
        @assert is_axis_query(metacells_subset)
        @assert query_axis_name(metacells_subset) == "metacell"
        metacells_subset = get_query(daf, metacells_subset)
    end
    indices_of_subset_metacells = findall(metacells_subset)
    n_subset_metacells = length(indices_of_subset_metacells)

    linear_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "linear_fraction").array
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)
    n_marker_genes = length(indices_of_markers)

    linear_fraction_per_subset_metacell_per_marker_gene = Matrix{Float32}(undef, n_subset_metacells, n_marker_genes)
    parallel_loop_wo_rng(
        1:n_marker_genes;
        progress = DebugProgress(
            n_marker_genes,
            group = :mcs_loops,
            desc = "linear_fraction_per_subset_metacell_per_marker_gene",
        ),
    ) do marker_gene_position
        marker_gene_index = indices_of_markers[marker_gene_position]
        @views linear_fraction_per_subset_metacell =
            linear_fraction_per_subset_metacell_per_marker_gene[:, marker_gene_position]
        return linear_fraction_per_subset_metacell .=
            linear_fraction_per_metacell_per_gene[indices_of_subset_metacells, marker_gene_index]
    end

    correlation_between_marker_genes =
        zero_cor_between_matrix_columns(linear_fraction_per_subset_metacell_per_marker_gene)
    correlation_between_genes = embed_dense_matrix_in_sparse_matrix(
        correlation_between_marker_genes;
        rows_indices = indices_of_markers,
        n_rows = n_marker_genes,
        columns_indices = indices_of_markers,
        n_columns = n_marker_genes,
    )
    set_matrix!(daf, "gene", "gene", matrix_name, correlation_between_genes; overwrite)
    return nothing
end

end  # module
