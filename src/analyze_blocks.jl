"""
Do simple blocks analysis.
"""
module AnalyzeBlocks

export compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!
export compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!
export compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!
export compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!
export compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!
export compute_matrix_of_is_in_neighborhood_per_block_per_block!
export compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!
export compute_matrix_of_is_neighborhood_marker_per_gene_per_block!
export compute_matrix_of_linear_fraction_per_gene_per_block!
export compute_matrix_of_log_linear_fraction_per_gene_per_block!
export compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!
export compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!
export compute_matrix_of_UMIs_per_gene_per_block!
export compute_vector_of_block_closest_by_pertinent_markers_per_cell!
export compute_vector_of_n_cells_per_block!
export compute_vector_of_n_metacells_per_block!
export compute_vector_of_n_neighborhood_blocks_per_block!
export compute_vector_of_n_neighborhood_cells_per_block!
export compute_vector_of_n_neighborhood_metacells_per_block!
export compute_vector_of_total_neighborhood_UMIs_per_block!
export compute_vector_of_total_UMIs_per_block!
export compute_vector_of_type_per_block_by_cells!
export compute_vector_of_type_per_block_by_metacells!
export compute_vector_of_type_per_cell_by_blocks!
export compute_vector_of_type_per_metacell_by_blocks!

using Base.Threads
using DataAxesFormats
using Distances
using LinearAlgebra
using StatsBase
using TanayLabUtilities

using ..AnalyzeGenes
using ..Defaults
using ..Contracts

import ..AnalyzeGenes.fill_vector_of_is_correlated_with_skeleton_per_gene!

# Needed because of JET:
import Metacells.Contracts.base_block_axis
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_block
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block
import Metacells.Contracts.matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_distinct_per_gene_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_marker_per_gene_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks
import Metacells.Contracts.matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_block_closest_by_pertinent_markers_per_cell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_block
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_n_metacells_per_block
import Metacells.Contracts.vector_of_n_neighborhood_blocks_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_n_neighborhood_metacells_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell
import Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
import Metacells.Contracts.vector_of_type_per_block
import Metacells.Contracts.vector_of_type_per_metacell

"""
    function compute_vector_of_n_metacells_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_metacells_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [vector_of_block_per_metacell(RequiredInput), vector_of_n_metacells_per_block(CreatedOutput)],
) function compute_vector_of_n_metacells_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells_per_block = daf["@ metacell / block =@ >> Count || 0"].array
    set_vector!(daf, "block", "n_metacells", n_metacells_per_block; overwrite)
    @debug "Mean metacells in block: $(mean(n_metacells_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    function compute_vector_of_n_cells_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_cells_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_n_cells_per_metacell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_cells_per_block(CreatedOutput),
    ],
) function compute_vector_of_n_cells_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells_per_block = daf["@ metacell : n_cells / block =@ >> Sum || 0"].array
    set_vector!(daf, "block", "n_cells", n_cells_per_block; overwrite)
    @debug "Mean cells in block: $(mean(n_cells_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    function compute_matrix_of_UMIs_per_gene_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_UMIs_per_gene_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_UMIs_per_gene_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    UMIs_per_gene_per_block = daf["@ gene @ metacell :: UMIs |/ block =@ >| Sum || 0"].array
    set_matrix!(daf, "gene", "block", "UMIs", bestify(UMIs_per_gene_per_block); overwrite)
    return nothing
end

"""
    function compute_vector_of_total_UMIs_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of non-excluded genes per block.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_block(RequiredInput),
        vector_of_total_UMIs_per_block(CreatedOutput),
    ],
) function compute_vector_of_total_UMIs_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_block = daf["@ block @ gene [ ! is_excluded ] :: UMIs >| Sum"].array
    set_vector!(daf, "block", "total_UMIs", total_UMIs_per_block; overwrite)
    @debug "Mean total UMIs in block: $(mean(total_UMIs_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    function compute_matrix_of_linear_fraction_per_gene_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_linear_fraction_per_gene_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_block(RequiredInput),
        vector_of_total_UMIs_per_block(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_linear_fraction_per_gene_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    indices_of_included_genes = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    UMIs_per_block_per_included_gene = daf["@ block @ gene [ ! is_excluded ] :: UMIs"].array
    total_included_UMIs_per_block = get_vector(daf, "block", "total_UMIs").array

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    linear_fraction_per_block_per_gene = zeros(Float32, n_blocks, n_genes)
    @views linear_fraction_per_block = linear_fraction_per_block_per_gene[:, indices_of_included_genes]
    @. linear_fraction_per_block = UMIs_per_block_per_included_gene / total_included_UMIs_per_block

    set_matrix!(daf, "block", "gene", "linear_fraction", bestify(linear_fraction_per_block_per_gene); overwrite)  # NOJET

    return nothing
end

"""
    function compute_matrix_of_log_linear_fraction_per_gene_per_block!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_log_linear_fraction_per_gene_per_block`](@ref). This adds the `gene_fraction_regularization`
to deal with zero fractions.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_log_linear_fraction_per_gene_per_block!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    overwrite::Bool = false,
)::Nothing
    fraction_per_gene_per_block = get_matrix(daf, "gene", "block", "linear_fraction").array
    empty_dense_matrix!(
        daf,
        "gene",
        "block",
        "log_linear_fraction",
        Float32;
        overwrite,
    ) do log_fraction_per_gene_per_block
        @. log_fraction_per_gene_per_block = log2(fraction_per_gene_per_block + gene_fraction_regularization)
        return nothing
    end
    return nothing
end

"""
    function compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The maximal significant skeleton genes fractions fold factor between metacells of the blocks.
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks(CreatedOutput),
    ],
) function compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")

    distances_between_metacells = get_matrix(daf, "metacell", "metacell", "euclidean_skeleton_fold_distance").array

    block_index_per_metacell = get_query(daf, "@ metacell : block : index").array
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)

    mean_distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    parallel_loop_wo_rng(
        reverse(1:n_blocks);
        progress = DebugProgress(
            n_blocks;
            group = :mcs_loops,
            desc = "mean_euclidean_skeleton_fold_distance_between_blocks",
        ),
    ) do base_block_index
        indices_of_base_metacells = metacell_indices_per_block[base_block_index]
        @views distance_per_metacell_per_base_metacell = distances_between_metacells[:, indices_of_base_metacells]

        mean_distance_from_base_per_metacell = vec(mean(distance_per_metacell_per_base_metacell; dims = 2))  # NOLINT
        @assert_vector(mean_distance_from_base_per_metacell, n_metacells)

        for other_block_index in 1:base_block_index
            indices_of_other_metacells = metacell_indices_per_block[other_block_index]
            n_other_block_metacells = length(indices_of_other_metacells)
            @assert n_other_block_metacells > 0

            mean_distance_between_base_and_other_block = 0.0
            for other_metacell_index in 1:n_other_block_metacells
                mean_distance_between_base_and_other_block += mean_distance_from_base_per_metacell[other_metacell_index]
            end
            mean_distance_between_base_and_other_block /= n_other_block_metacells

            mean_distances_between_blocks[base_block_index, other_block_index] =
                mean_distance_between_base_and_other_block
            mean_distances_between_blocks[other_block_index, base_block_index] =
                mean_distance_between_base_and_other_block
        end

        return nothing
    end

    set_matrix!(
        daf,
        "block",
        "block",
        "mean_euclidean_skeleton_fold_distance",
        mean_distances_between_blocks;
        overwrite,
    )
    return nothing
end

"""
    compute_vector_of_block_closest_by_pertinent_markers_per_cell!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = false
    )::Nothing

Compute and set [`vector_of_block_closest_by_pertinent_markers_per_cell`](@ref).
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), cell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_block_closest_by_pertinent_markers_per_cell(CreatedOutput),
    ],
) function compute_vector_of_block_closest_by_pertinent_markers_per_cell!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_cells = axis_length(daf, "cell")

    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    UMIs_per_cell_per_pertinent_marker = daf["@ cell @ gene [ is_marker & ! is_lateral ] :: UMIs"].array
    UMIs_per_pertinent_marker_per_cell = flipped(UMIs_per_cell_per_pertinent_marker)
    log_linear_fraction_per_pertinent_marker_per_cell = log2.(
        densify(UMIs_per_pertinent_marker_per_cell) ./ transpose(total_UMIs_per_cell) .+ gene_fraction_regularization,
    )

    linear_fraction_per_block_per_pertinent_marker =
        daf["@ block @ gene [ is_marker & ! is_lateral ] :: linear_fraction"].array
    linear_fraction_per_pertinent_marker_per_block = flipped(linear_fraction_per_block_per_pertinent_marker)
    log_linear_fraction_per_pertinent_marker_per_block =
        log2.(linear_fraction_per_pertinent_marker_per_block .+ gene_fraction_regularization)

    distances_between_blocks_and_cells = parallel_pairwise(
        Euclidean(),
        log_linear_fraction_per_pertinent_marker_per_block,
        log_linear_fraction_per_pertinent_marker_per_cell,
        dims = 2,
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "distances_between_blocks_and_cells"),
    )

    name_per_block = axis_vector(daf, "block")
    closest_block_name_per_cell =
        name_per_block[first.(Tuple.(vec(argmin(distances_between_blocks_and_cells; dims = 1))))]
    @assert_vector(closest_block_name_per_cell, n_cells)

    set_vector!(daf, "cell", "block.closest_by_pertinent_markers", closest_block_name_per_cell; overwrite)

    return nothing
end

"""
    compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!(daf::DafWriter; overwrite::Bool = false)::Nothing

Compute and set [`matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block`](@ref).
$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_block_closest_by_pertinent_markers_per_cell(RequiredInput),
        matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block(CreatedOutput),
    ],
) function compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")

    closest_block_index_per_grouped_cell = daf["@ cell [ metacell ] : block.closest_by_pertinent_markers : index"].array
    original_block_index_per_grouped_cell = get_query(daf, "@ cell : metacell ?? : block : index").array
    n_grouped_cells = length(closest_block_index_per_grouped_cell)

    confusion_per_other_block_per_original_block = zeros(UInt32, n_blocks, n_blocks)
    for (closest_block_index, original_block_index) in
        zip(closest_block_index_per_grouped_cell, original_block_index_per_grouped_cell)
        confusion_per_other_block_per_original_block[closest_block_index, original_block_index] += 1
    end

    set_matrix!(
        daf,
        "block",
        "block",
        "confusion_by_closest_by_pertinent_markers",
        confusion_per_other_block_per_original_block;
        overwrite,
    )

    n_stable_cells =
        sum(confusion_per_other_block_per_original_block[diagind(confusion_per_other_block_per_original_block)])

    @debug "Stable cells: $(n_stable_cells) ($(percent(n_stable_cells, n_grouped_cells)))" _group = :mcs_results

    return nothing
end

"""
    compute_matrix_of_is_in_neighborhood_per_block_per_block!(
        daf::DafWriter;
        min_neighbour_confusion_fractions::AbstractFloat = $(DEFAULT.min_neighbour_confusion_fractions),
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        min_total_UMIs_in_neighborhood::Integer = $(DEFAULT.min_total_UMIs_in_neighborhood),
        max_blocks_in_neighborhood::Integer = $(DEFAULT.max_blocks_in_neighborhood),
        max_metacells_in_neighborhood::Integer = $(DEFAULT.max_metacells_in_neighborhood)
        max_total_UMIs_in_neighborhood::Integer = $(DEFAULT.max_total_UMIs_in_neighborhood),
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_is_in_neighborhood_per_block_per_block`](@ref).

  - We normalize the confusion matrix between blocks so that each column sums to 1 (that is, we look at the fraction of
    the cells of each block that are confused with each other block).
  - We then make it symmetric by adding it to its transpose.

For each block, we sort all blocks by the normalized symmetric confusion score. To break ties (especially for zero
confusion), we use the mean distance between the blocks. We then pick the first few blocks to be in the neighborhood.

  - All blocks that have a symmetric normalized confusion fraction of at least `min_neighbour_confusion_fractions` are
    included, but we stop if we exceed `max_blocks_in_neighborhood`, `max_metacells_in_neighborhood`, or
    `max_total_UMIs_in_neighborhood`.
  - If we didn't reach `min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` or
    `min_total_UMIs_in_neighborhood`, we keep adding blocks even if their symmetric normalized confusion fraction is
    low (even zero).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        vector_of_n_cells_per_block(RequiredInput),
        matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block(RequiredInput),
        vector_of_total_UMIs_per_block(RequiredInput),
        vector_of_n_metacells_per_block(RequiredInput),
        matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(CreatedOutput),
    ],
) function compute_matrix_of_is_in_neighborhood_per_block_per_block!(
    daf::DafWriter;
    min_neighbour_confusion_fractions::AbstractFloat = 0.01,
    min_blocks_in_neighborhood::Integer = 5,
    min_metacells_in_neighborhood::Integer = 20,
    min_total_UMIs_in_neighborhood::Integer = 2_000_000,
    max_blocks_in_neighborhood::Integer = 100,
    max_metacells_in_neighborhood::Integer = 400,
    max_total_UMIs_in_neighborhood::Integer = 4_000_000_0,
    overwrite::Bool = false,
)::Nothing
    @assert 0 < min_blocks_in_neighborhood < max_blocks_in_neighborhood
    @assert 0 < min_metacells_in_neighborhood < max_metacells_in_neighborhood
    @assert 0 < min_total_UMIs_in_neighborhood < max_total_UMIs_in_neighborhood

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    n_cells_per_block = get_vector(daf, "block", "n_cells").array

    mean_euclidean_skeleton_fold_distance_per_block_per_block =
        get_matrix(daf, "block", "block", "mean_euclidean_skeleton_fold_distance")
    n_metacells_per_block = get_vector(daf, "block", "n_metacells")
    total_UMIs_per_block = get_vector(daf, "block", "total_UMIs")
    confusion_per_other_block_per_original_block =
        get_matrix(daf, "block", "block", "confusion_by_closest_by_pertinent_markers").array

    confusion_fraction_per_block_per_block = Float32.(confusion_per_other_block_per_original_block ./ n_cells_per_block)
    confusion_fraction_per_block_per_block .+= transpose(confusion_fraction_per_block_per_block)

    is_in_neighborhood_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)
    priority_per_other_block_per_thread = [Vector{Tuple{Float32, Float32}}(undef, n_blocks) for _ in 1:nthreads()]
    ordered_block_indices_per_thread = [Vector{Int}(undef, n_blocks) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "is_in_neighborhood_per_block_per_block"),
        policy = :static,
    ) do base_block_index
        priority_per_other_block = priority_per_other_block_per_thread[threadid()]
        ordered_block_indices = ordered_block_indices_per_thread[threadid()]
        @views confusion_fraction_per_other_block = confusion_fraction_per_block_per_block[:, base_block_index]
        @views fold_distance_per_other_block =
            mean_euclidean_skeleton_fold_distance_per_block_per_block[:, base_block_index]

        for i in 1:n_blocks
            priority_per_other_block[i] =
                (-confusion_fraction_per_other_block[i], Float32(fold_distance_per_other_block[i]))
        end
        priority_per_other_block[base_block_index] = (-Inf32, -Inf32)
        sortperm!(ordered_block_indices, priority_per_other_block)

        neighborhood_n_confused_blocks = 0
        neighborhood_n_disjoint_blocks = 0
        neighborhood_n_overflow_blocks = 0

        neighborhood_n_blocks = 0
        neighborhood_n_metacells = 0
        neighborhood_total_UMIs = 0

        while true
            next_block_index = ordered_block_indices[neighborhood_n_blocks + 1]
            include_next_block = false

            if neighborhood_n_blocks == 0
                include_next_block = true
                neighborhood_n_confused_blocks += 1
            elseif confusion_fraction_per_other_block[next_block_index] >= min_neighbour_confusion_fractions
                if neighborhood_n_blocks <= max_blocks_in_neighborhood &&
                   neighborhood_n_metacells <= max_metacells_in_neighborhood &&
                   neighborhood_total_UMIs <= max_total_UMIs_in_neighborhood
                    include_next_block = true
                    neighborhood_n_confused_blocks += 1
                end
            else
                if neighborhood_n_blocks < min_blocks_in_neighborhood ||
                   neighborhood_n_metacells < min_metacells_in_neighborhood ||
                   neighborhood_total_UMIs < min_total_UMIs_in_neighborhood
                    include_next_block = true
                    neighborhood_n_disjoint_blocks += 1
                end
            end

            if !include_next_block
                break
            end

            neighborhood_n_blocks += 1
            is_in_neighborhood_per_other_block_per_base_block[next_block_index, base_block_index] = true
            neighborhood_n_metacells += n_metacells_per_block[next_block_index]
            neighborhood_total_UMIs += total_UMIs_per_block[next_block_index]
        end

        while neighborhood_n_blocks + neighborhood_n_overflow_blocks < n_blocks
            next_block_index = ordered_block_indices[neighborhood_n_blocks + neighborhood_n_overflow_blocks + 1]
            if confusion_fraction_per_other_block[next_block_index] >= min_neighbour_confusion_fractions
                neighborhood_n_overflow_blocks += 1
            else
                break
            end
        end

        @debug (
            "Neighborhood: $(name_per_block[base_block_index])" *
            " Blocks: $(neighborhood_n_blocks)" *
            " Confused: $(neighborhood_n_confused_blocks)" *
            " Disjoint: $(neighborhood_n_disjoint_blocks)" *
            " Overflow: $(neighborhood_n_overflow_blocks)" *
            " Metacells: $(neighborhood_n_metacells)" *
            " Covered M-UMIs: $(neighborhood_total_UMIs / 1e6)"
        ) _group = :mcs_details
        @assert neighborhood_n_disjoint_blocks == 0 || neighborhood_n_overflow_blocks == 0
        return nothing
    end

    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_neighborhood",
        is_in_neighborhood_per_other_block_per_base_block;
        overwrite,
    )
    return nothing
end

"""
    function compute_vector_of_n_neighborhood_blocks_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_neighborhood_blocks_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_n_neighborhood_blocks_per_block(CreatedOutput),
    ],
) function compute_vector_of_n_neighborhood_blocks_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_neighborhood_blocks_per_block = daf["@ block @ block :: is_in_neighborhood >- Sum"].array
    set_vector!(daf, "block", "n_neighborhood_blocks", n_neighborhood_blocks_per_block; overwrite)
    @debug "Mean blocks in neighborhood: $(mean(n_neighborhood_blocks_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    function compute_vector_of_n_neighborhood_metacells_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_neighborhood_metacells_per_block`](@ref)

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        vector_of_n_metacells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_n_neighborhood_metacells_per_block(CreatedOutput),
    ],
) function compute_vector_of_n_neighborhood_metacells_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    compute_vector_of_neighborhood_something_per_block(
        daf;
        overwrite,
        vector_property = "n_metacells",
        matrix_property = "n_neighborhood_metacells",
        result_name = "metacells",
    )
    return nothing
end

"""
    function compute_vector_of_n_neighborhood_cells_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of cells in the metacells of the blocks of the neighborhood centered at a block.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        vector_of_n_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(CreatedOutput),
    ],
) function compute_vector_of_n_neighborhood_cells_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    compute_vector_of_neighborhood_something_per_block(
        daf;
        overwrite,
        vector_property = "n_cells",
        matrix_property = "n_neighborhood_cells",
        result_name = "cells",
    )
    return nothing
end

"""
    function compute_vector_of_total_neighborhood_UMIs_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_total_neighborhood_UMIs_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        vector_of_total_UMIs_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_total_neighborhood_UMIs_per_block(CreatedOutput),
    ],
) function compute_vector_of_total_neighborhood_UMIs_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    compute_vector_of_neighborhood_something_per_block(
        daf;
        overwrite,
        vector_property = "total_UMIs",
        matrix_property = "total_neighborhood_UMIs",
        result_name = "total UMIs",
    )
    return nothing
end

function compute_vector_of_neighborhood_something_per_block(
    daf::DafWriter;
    overwrite::Bool = false,
    vector_property::AbstractString,
    matrix_property::AbstractString,
    result_name::AbstractString,
)::Nothing
    n_blocks = axis_length(daf, "block")

    value_per_block = get_vector(daf, "block", vector_property).array
    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array

    neighborhood_value_per_block = Vector{UInt32}(undef, n_blocks)
    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "$(matrix_property)_per_block"),
    ) do block_index
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        neighborhood_value_per_block[block_index] = dot(is_in_neighborhood_per_other_block, value_per_block)
        return nothing
    end

    set_vector!(daf, "block", matrix_property, neighborhood_value_per_block; overwrite)
    @debug "Mean $(result_name) in neighborhood: $(mean(neighborhood_value_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    compute_matrix_of_is_neighborhood_marker_per_gene_per_block!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = 2 ^ -13.5,
        min_marker_gene_range_fold::Real = 1.0,
        min_marker_quantile::Real = 0.1,
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_is_neighborhood_marker_per_gene_per_block`](@ref).

We only consider genes which are markers in the overall population. We then call
[`compute_vector_of_is_marker_per_gene!`](@ref) but only looking at the metacells in each neighborhood.
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_is_neighborhood_marker_per_gene_per_block!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 2 ^ -13.5,
    min_marker_gene_range_fold::Real = 1.0,
    min_marker_quantile::Real = 0.1,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_neighborhood_marker_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "is_neighborhood_marker_per_gene_per_block"),
    ) do block_index
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = [
                "metacell" => "@ metacell [ block :: is_in_neighborhood @| $(block_name) ]",
                "gene" => "@ gene [ is_marker ]",
            ],
            input_data = [
                ("metacell", "gene", "linear_fraction") => "=",
                ("metacell", "gene", "log_linear_fraction") => "=",
                ("gene", "full_index") => "index",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            compute_vector_of_is_marker_per_gene!(
                adapted;
                min_marker_gene_max_fraction,
                min_marker_gene_range_fold,
                min_marker_quantile,
            )
            full_index_per_gene = get_vector(adapted, "gene", "full_index").array
            is_neighborhood_marker_per_gene = get_vector(adapted, "gene", "is_marker").array
            is_neighborhood_marker_per_gene_per_block[full_index_per_gene, block_index] =
                is_neighborhood_marker_per_gene
            return nothing
        end
        return nothing
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_neighborhood_marker",
        bestify(is_neighborhood_marker_per_gene_per_block);
        overwrite,
    )
    @debug "Mean markers in neighborhood: $(sum(is_neighborhood_marker_per_gene_per_block) / n_blocks)" _group =
        :mcs_results
    return nothing
end

"""
    compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!(
        daf::DafWriter;
        min_distinct_gene_max_fraction::AbstractFloat = 2 ^ -14.5,
        min_distinct_gene_mean_fold::Real = 2,
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_is_in_neighborhood_per_block_per_block`](@ref). A gene is distinct in the neighborhood
if:

  - It has a maximal expression level of at least `min_distinct_gene_max_fraction`.
  - When comparing its fold factor compared to its median expression level in the population, then the mean fold factor
    in the neighborhood is at least `min_distinct_gene_mean_fold`. This is not the absolute fold, that is, we only look
    for genes which are stronger than the population as a whole, not weaker.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_is_neighborhood_distinct_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!(
    daf::DafWriter;
    min_distinct_gene_max_fraction::AbstractFloat = 2 ^ -14.5,
    min_distinct_gene_mean_fold::Real = 2,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")

    linear_fraction_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", "linear_fraction").array
    log_linear_fraction_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", "log_linear_fraction").array
    median_log_linear_fraction_per_gene = daf["@ metacell @ gene :: log_linear_fraction >- Median"].array
    @assert_vector(median_log_linear_fraction_per_gene, n_genes)

    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array
    block_index_per_metacell = daf["@ metacell : block : index"].array

    n_metacells = length(block_index_per_metacell)
    is_neighborhood_distinct_per_gene_per_block = zeros(Bool, n_genes, n_blocks)
    is_in_neighborhood_per_metacell_per_thread = [BitVector(undef, n_metacells) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "is_neighborhood_distinct_per_gene_per_block"),
        policy = :static,
    ) do block_index
        is_in_neighborhood_per_metacell = is_in_neighborhood_per_metacell_per_thread[threadid()]
        @views is_neighborhood_distinct_per_gene = is_neighborhood_distinct_per_gene_per_block[:, block_index]
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]

        is_in_neighborhood_per_metacell .= getindex.(Ref(is_in_neighborhood_per_other_block), block_index_per_metacell)

        linear_fraction_per_gene_per_neighborhood_metacell =
            linear_fraction_per_gene_per_metacell[:, is_in_neighborhood_per_metacell]
        log_linear_fraction_per_gene_per_neighborhood_metacell =
            log_linear_fraction_per_gene_per_metacell[:, is_in_neighborhood_per_metacell]
        max_linear_fraction_per_gene = vec(maximum(linear_fraction_per_gene_per_neighborhood_metacell; dims = 2))

        fold_per_gene_per_neighborhood_metacell =
            log_linear_fraction_per_gene_per_neighborhood_metacell .- median_log_linear_fraction_per_gene
        mean_fold_per_gene = vec(mean(fold_per_gene_per_neighborhood_metacell; dims = 2))  # NOLINT

        is_neighborhood_distinct_per_gene .=
            (mean_fold_per_gene .>= min_distinct_gene_mean_fold) .&
            (max_linear_fraction_per_gene .>= min_distinct_gene_max_fraction)

        return nothing
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "is_neighborhood_distinct",
        bestify(is_neighborhood_distinct_per_gene_per_block);
        overwrite,
    )

    @debug "Mean distincts in neighborhood: $(sum(is_neighborhood_distinct_per_gene_per_block) / n_blocks)" _group =
        :mcs_results

    return nothing
end

"""
    compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!(
        daf::DafWriter;
        min_gene_correlation::AbstractFloat = function_default(
            compute_vector_of_is_correlated_with_skeleton_per_gene!,
            :min_gene_correlation,
        ),
        min_gene_correlation_quantile::AbstractFloat = function_default(
            compute_vector_of_is_correlated_with_skeleton_per_gene!,
            :min_gene_correlation_quantile,
        ),
        genes_correlation_window::Integer = function_default(
            compute_vector_of_is_correlated_with_skeleton_per_gene!,
            :genes_correlation_window,
        ),
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block`](@ref). This just invokes
[`compute_vector_of_is_correlated_with_skeleton_per_gene!`](@ref) for the metacells in each neighborhood.
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_n_neighborhood_metacells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!(
    daf::DafWriter;
    min_gene_correlation::AbstractFloat = function_default(
        compute_vector_of_is_correlated_with_skeleton_per_gene!,
        :min_gene_correlation,
    ),
    min_gene_correlation_quantile::AbstractFloat = function_default(
        compute_vector_of_is_correlated_with_skeleton_per_gene!,
        :min_gene_correlation_quantile,
    ),
    genes_correlation_window::Integer = function_default(
        compute_vector_of_is_correlated_with_skeleton_per_gene!,
        :genes_correlation_window,
    ),
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    n_metacells = axis_length(daf, "metacell")

    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)
    n_markers = length(indices_of_markers)

    is_skeletons_per_gene = get_vector(daf, "gene", "is_skeleton").array
    n_skeletons = sum(is_skeletons_per_gene)

    log_fraction_per_metacell_per_skeleton = daf["@ metacell @ gene [ is_skeleton ] :: log_linear_fraction"].array
    log_fraction_per_metacell_per_marker = daf["@ metacell @ gene [ is_marker ] :: log_linear_fraction"].array

    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array
    block_index_per_metacell = daf["@ metacell : block : index"].array

    is_correlated_with_skeleton_in_neighborhood_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    n_neighborhood_metacells_per_block = get_vector(daf, "block", "n_neighborhood_metacells").array
    max_n_neighborhood_metacells = maximum(n_neighborhood_metacells_per_block)

    is_in_neighborhood_per_metacell_per_thread = [BitVector(undef, n_metacells) for _ in 1:nthreads()]
    tmp_per_marker_per_thread = Matrix{Float32}(undef, n_markers, nthreads())
    correlation_per_skeleton_per_marker_per_thread =
        [Matrix{Float32}(undef, n_skeletons, n_markers) for _ in 1:nthreads()]
    scratch_per_max_metacell_per_skeleton_per_thread =
        [Matrix{Float32}(undef, max_n_neighborhood_metacells, n_skeletons) for _ in 1:nthreads()]
    scratch_per_max_metacell_per_marker_per_thread =
        [Matrix{Float32}(undef, max_n_neighborhood_metacells, n_markers) for _ in 1:nthreads()]
    scratch_per_skeleton_per_thread = [Vector{Float32}(undef, n_skeletons) for _ in 1:nthreads()]
    scratch_per_marker_per_thread = [Vector{Float32}(undef, n_markers) for _ in 1:nthreads()]
    scratch_markers_window_per_thread = [Vector{Float32}(undef, genes_correlation_window) for _ in 1:nthreads()]

    log_fraction_per_max_neighborhood_metacell_per_skeleton_per_thread =
        [Matrix{Float32}(undef, max_n_neighborhood_metacells, n_skeletons) for _ in 1:nthreads()]
    log_fraction_per_max_neighborhood_metacell_per_marker_per_thread =
        [Matrix{Float32}(undef, max_n_neighborhood_metacells, n_markers) for _ in 1:nthreads()]
    indices_of_neighborhood_metacells_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_metacells) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_blocks;
        policy = :static,
        progress = DebugProgress(
            n_blocks;
            group = :mcs_loops,
            desc = "is_correlated_with_skeleton_in_neighborhood_per_gene_per_block",
        ),
    ) do block_index
        @views tmp_per_marker = tmp_per_marker_per_thread[:, threadid()]
        correlation_per_skeleton_per_marker = correlation_per_skeleton_per_marker_per_thread[threadid()]
        scratch_per_max_metacell_per_skeleton = scratch_per_max_metacell_per_skeleton_per_thread[threadid()]
        scratch_per_max_metacell_per_marker = scratch_per_max_metacell_per_marker_per_thread[threadid()]
        scratch_per_skeleton = scratch_per_skeleton_per_thread[threadid()]
        scratch_per_marker = scratch_per_marker_per_thread[threadid()]
        scratch_marker_window = scratch_markers_window_per_thread[threadid()]

        log_fraction_per_max_neighborhood_metacell_per_skeleton =
            log_fraction_per_max_neighborhood_metacell_per_skeleton_per_thread[threadid()]
        log_fraction_per_max_neighborhood_metacell_per_marker =
            log_fraction_per_max_neighborhood_metacell_per_marker_per_thread[threadid()]

        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_metacell = is_in_neighborhood_per_metacell_per_thread[threadid()]
        is_in_neighborhood_per_metacell .= getindex.(Ref(is_in_neighborhood_per_other_block), block_index_per_metacell)

        indices_of_max_neighborhood_metacells = indices_of_neighborhood_metacells_per_thread[threadid()]
        n_neighborhood_metacells = 0
        @foreach_true_index is_in_neighborhood_per_metacell metacell_index begin
            n_neighborhood_metacells += 1
            indices_of_max_neighborhood_metacells[n_neighborhood_metacells] = metacell_index
        end
        @assert n_neighborhood_metacells == n_neighborhood_metacells_per_block[block_index]
        @views indices_of_neighborhood_metacells = indices_of_max_neighborhood_metacells[1:n_neighborhood_metacells]

        @views scratch_per_neighborhood_metacell_per_skeleton =
            scratch_per_max_metacell_per_skeleton[1:n_neighborhood_metacells, :]
        @views scratch_per_neighborhood_metacell_per_marker =
            scratch_per_max_metacell_per_marker[1:n_neighborhood_metacells, :]

        @views log_fraction_per_neighborhood_metacell_per_skeleton =
            log_fraction_per_max_neighborhood_metacell_per_skeleton[1:n_neighborhood_metacells, :]
        @views log_fraction_per_neighborhood_metacell_per_marker =
            log_fraction_per_max_neighborhood_metacell_per_marker[1:n_neighborhood_metacells, :]

        for skeleton_index in 1:n_skeletons
            @views destination = log_fraction_per_neighborhood_metacell_per_skeleton[:, skeleton_index]
            @views source = log_fraction_per_metacell_per_skeleton[:, skeleton_index]
            destination .= getindex.(Ref(source), indices_of_neighborhood_metacells)
        end
        for marker_index in 1:n_markers
            @views destination = log_fraction_per_neighborhood_metacell_per_marker[:, marker_index]
            @views source = log_fraction_per_metacell_per_marker[:, marker_index]
            destination .= getindex.(Ref(source), indices_of_neighborhood_metacells)
        end

        @views is_correlated_with_skeleton_in_neighborhood_per_gene =
            is_correlated_with_skeleton_in_neighborhood_per_gene_per_block[:, block_index]

        fill_vector_of_is_correlated_with_skeleton_per_gene!(;
            min_gene_correlation,
            min_gene_correlation_quantile,
            genes_correlation_window,
            indices_of_markers,
            log_fraction_per_metacell_per_skeleton = log_fraction_per_neighborhood_metacell_per_skeleton,
            log_fraction_per_metacell_per_marker = log_fraction_per_neighborhood_metacell_per_marker,
            is_correlated_with_skeleton_per_gene = is_correlated_with_skeleton_in_neighborhood_per_gene,
            correlation_per_skeleton_per_marker,
            max_correlation_per_marker = tmp_per_marker,
            scratch_per_metacell_per_skeleton = scratch_per_neighborhood_metacell_per_skeleton,
            scratch_per_metacell_per_marker = scratch_per_neighborhood_metacell_per_marker,
            scratch_per_skeleton,
            scratch_per_marker,
            scratch_marker_window,
        )

        return nothing
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "is_correlated_with_skeleton_in_neighborhood",
        bestify(is_correlated_with_skeleton_in_neighborhood_per_gene_per_block);
        overwrite,
    )

    @debug (
        "Mean markers correlated with skeleton in neighborhood: " *
        "$(sum(is_correlated_with_skeleton_in_neighborhood_per_gene_per_block) / n_blocks)"
    ) _group = :mcs_results

    return nothing
end

"""
    compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool,
    )::Nothing

Compute and set [`matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block(CreatedOutput),
        vector_of_is_lateral_per_gene(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(RequiredInput),
    ],
) function compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(daf, "cell")
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")

    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_neighborhood_marker_per_gene_per_block = get_matrix(daf, "gene", "block", "is_neighborhood_marker").array

    indices_of_included_genes = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    n_included_genes = length(indices_of_included_genes)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array

    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    metacell_index_per_cell = get_query(daf, "@ cell : metacell ?? 0 : index").array
    block_index_per_cell = daf["@ cell : metacell ?? 0 : block : index"].array

    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array

    is_neighborhood_pertinent_marker_per_gene = BitVector(undef, n_genes)

    correlation_per_gene_per_block = zeros(Float32, n_genes, n_blocks)
    mean_correlation_per_block = Vector{Float32}(undef, n_blocks)

    progress = DebugProgress(
        n_blocks * n_included_genes;
        group = :mcs_loops,
        desc = "correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block",
    )

    n_neighborhood_cells_per_block = get_vector(daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    cell_log_fraction_per_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    punctuated_metacell_log_fraction_per_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    is_in_neighborhood_per_cell = BitVector(undef, n_cells)

    metacell_index_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_UMIs_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_metacell_UMIs_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_punctuated_metacell_UMIs_per_max_neighborhood_cell = total_metacell_UMIs_per_max_neighborhood_cell
    indices_of_max_neighborhood_cells = Vector{Int}(undef, max_n_neighborhood_cells)

    todox_all_zero = Atomic{Int32}(0)
    for block_index in 1:n_blocks
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_cell .=
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1))

        n_neighborhood_cells = 0
        @foreach_true_index is_in_neighborhood_per_cell cell_index begin
            n_neighborhood_cells += 1
            indices_of_max_neighborhood_cells[n_neighborhood_cells] = cell_index
        end
        @assert n_neighborhood_cells == n_neighborhood_cells_per_block[block_index]
        @views indices_of_neighborhood_cells = indices_of_max_neighborhood_cells[1:n_neighborhood_cells]

        @views metacell_index_per_neighborhood_cell = metacell_index_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_UMIs_per_neighborhood_cell = total_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_metacell_UMIs_per_neighborhood_cell =
            total_metacell_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_punctuated_metacell_UMIs_per_neighborhood_cell =
            total_punctuated_metacell_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]

        metacell_index_per_neighborhood_cell .= getindex.(Ref(metacell_index_per_cell), indices_of_neighborhood_cells)
        total_UMIs_per_neighborhood_cell .= getindex.(Ref(total_UMIs_per_cell), indices_of_neighborhood_cells)
        total_metacell_UMIs_per_neighborhood_cell .=
            getindex.(Ref(total_UMIs_per_metacell), metacell_index_per_neighborhood_cell)
        @. total_punctuated_metacell_UMIs_per_neighborhood_cell =
            total_metacell_UMIs_per_neighborhood_cell .- total_UMIs_per_neighborhood_cell

        parallel_loop_wo_rng(1:n_included_genes; progress, policy = :static, progress_chunk = 100) do included_gene_position
            gene_index = indices_of_included_genes[included_gene_position]

            cell_log_fraction_per_max_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell_per_thread[threadid()]
            punctuated_metacell_log_fraction_per_max_neighborhood_cell =
                punctuated_metacell_log_fraction_per_max_neighborhood_cell_per_thread[threadid()]

            @views cell_log_fraction_per_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell[1:n_neighborhood_cells]
            @views punctuated_metacell_log_fraction_per_neighborhood_cell =
                punctuated_metacell_log_fraction_per_max_neighborhood_cell[1:n_neighborhood_cells]

            all_zero_neighborhood_cell_UMIs = true
            for neighborhood_cell_position in 1:n_neighborhood_cells
                if UMIs_per_cell_per_gene[indices_of_neighborhood_cells[neighborhood_cell_position], gene_index] > 0
                    all_zero_neighborhood_cell_UMIs = false
                    break
                end
            end
            if all_zero_neighborhood_cell_UMIs
                atomic_add!(todox_all_zero, Int32(1))
                correlation_per_gene_per_block[gene_index, block_index] = 0.0f0
                return nothing
            end

            all_zero_neighborhood_punctuated_metacell_UMIs = true
            for neighborhood_cell_position in 1:n_neighborhood_cells
                cell_index = indices_of_neighborhood_cells[neighborhood_cell_position]
                metacell_index = metacell_index_per_neighborhood_cell[neighborhood_cell_position]
                if UMIs_per_metacell_per_gene[metacell_index, gene_index] > UMIs_per_cell_per_gene[cell_index, gene_index]
                    all_zero_neighborhood_punctuated_metacell_UMIs = false
                    break
                end
            end
            if all_zero_neighborhood_punctuated_metacell_UMIs
                atomic_add!(todox_all_zero, Int32(1))
                correlation_per_gene_per_block[gene_index, block_index] = 0.0f0
                return nothing
            end

            for (neighborhood_cell_position, cell_index) in enumerate(indices_of_neighborhood_cells)
                metacell_index = metacell_index_per_neighborhood_cell[neighborhood_cell_position]
                UMIs_in_cell = Float32(UMIs_per_cell_per_gene[cell_index, gene_index])
                UMIs_in_metacell = Float32(UMIs_per_metacell_per_gene[metacell_index, gene_index])
                total_UMIs_in_cell = Float32(total_UMIs_per_neighborhood_cell[neighborhood_cell_position])
                total_UMIs_in_punctuated_metacell =
                    Float32(total_punctuated_metacell_UMIs_per_neighborhood_cell[neighborhood_cell_position])
                cell_log_fraction_per_neighborhood_cell[neighborhood_cell_position] =
                    log2(UMIs_in_cell / total_UMIs_in_cell + gene_fraction_regularization)
                punctuated_metacell_log_fraction_per_neighborhood_cell[neighborhood_cell_position] = log2(
                    (UMIs_in_metacell - UMIs_in_cell) / total_UMIs_in_punctuated_metacell +
                    gene_fraction_regularization,
                )
            end

            cell_log_fraction_per_neighborhood_cell .-= mean(cell_log_fraction_per_neighborhood_cell)  # NOLINT
            punctuated_metacell_log_fraction_per_neighborhood_cell .-=
                mean(punctuated_metacell_log_fraction_per_neighborhood_cell)  # NOLINT

            cell_log_fraction_per_neighborhood_cell ./= std(cell_log_fraction_per_neighborhood_cell)  # NOLINT
            punctuated_metacell_log_fraction_per_neighborhood_cell ./=
                std(punctuated_metacell_log_fraction_per_neighborhood_cell)  # NOLINT

            correlation =
                dot(cell_log_fraction_per_neighborhood_cell, punctuated_metacell_log_fraction_per_neighborhood_cell)
            correlation /= n_neighborhood_cells - 1
            if isnan(correlation)
                correlation = 0.0
            end

            correlation_per_gene_per_block[gene_index, block_index] = correlation

            return nothing
        end

        @views is_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene_per_block[:, block_index]
        @. is_neighborhood_pertinent_marker_per_gene = is_neighborhood_marker_per_gene & !is_lateral_per_gene
        mean_correlation_per_block[block_index] =
            mean(correlation_per_gene_per_block[is_neighborhood_pertinent_marker_per_gene, block_index])  # NOLINT
    end
    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(n_blocks * n_included_genes) ($(percent(todox_all_zero[], n_blocks * n_included_genes)))" _group = :todox

    set_matrix!(
        daf,
        "gene",
        "block",
        "correlation_between_neighborhood_cells_and_punctuated_metacells",
        bestify(correlation_per_gene_per_block);
        overwrite,
    )

    @debug (
        "Mean correlation of neighborhood pertinent marker genes between neighborhood cells and their punctuated metacells: " *
        "$(mean(mean_correlation_per_block))"  # NOLINT
    ) _group = :mcs_results
    return nothing
end

"""
    compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!(;
        other_daf::DafWriter,
        base_daf::DafReader,
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set
`[matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block`](@ref). This
will also copy [`base_block_axis`](@ref) from the `base_daf` into the `other_daf` if needed.

# Other

$(CONTRACT1)

# Base

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(
    name = "other_daf",
    axes = [
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        base_block_axis(GuaranteedOutput),
    ],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            CreatedOutput,
        ),
    ],
) Contract(
    name = "base_daf",
    axes = [
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
    ],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
    ],
) function compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!(;
    other_daf::DafWriter,
    base_daf::DafReader,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    bin::Maybe{Integer} = nothing,
    overwrite::Bool = false,
)::Nothing
    n_base_blocks = axis_length(base_daf, "block")
    n_cells = axis_length(base_daf, "cell")
    n_genes = axis_length(base_daf, "gene")

    @assert axis_vector(base_daf, "gene") == axis_vector(other_daf, "gene")
    @assert axis_vector(base_daf, "cell") == axis_vector(other_daf, "cell")
    @assert get_vector(base_daf, "gene", "is_excluded").array == get_vector(other_daf, "gene", "is_excluded").array
    @assert get_vector(base_daf, "gene", "is_lateral").array == get_vector(other_daf, "gene", "is_lateral").array

    if has_axis(other_daf, "base_block")
        @assert axis_vector(other_daf, "base_block") == axis_vector(base_daf, "block")
    else
        copy_axis!(source = base_daf, destination = other_daf, axis = "block", rename = "base_block"; overwrite)
    end

    indices_of_included_genes = get_query(base_daf, "@ gene [ ! is_excluded ] : index").array
    n_included_genes = length(indices_of_included_genes)

    is_in_neighborhood_per_other_block_per_base_block =
        get_matrix(base_daf, "block", "block", "is_in_neighborhood").array

    base_block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array

    is_lateral_per_gene = get_vector(base_daf, "gene", "is_lateral").array
    is_base_neighborhood_marker_per_gene_per_base_block =
        get_matrix(base_daf, "gene", "block", "is_neighborhood_marker").array

    other_metacell_index_per_cell = get_query(other_daf, "@ cell : metacell ?? 0 : index").array

    total_UMIs_per_cell = get_vector(other_daf, "cell", "total_UMIs").array
    total_UMIs_per_other_metacell = get_vector(other_daf, "metacell", "total_UMIs").array
    UMIs_per_cell_per_gene = get_matrix(other_daf, "cell", "gene", "UMIs").array
    UMIs_per_other_metacell_per_gene = get_matrix(other_daf, "metacell", "gene", "UMIs").array

    correlation_per_gene_per_base_block = zeros(Float32, n_genes, n_base_blocks)

    is_base_neighborhood_pertinent_marker_per_gene = BitVector(undef, n_genes)

    is_relevant_per_gene = nothing
    is_in_bin_per_gene = nothing

    mean_correlation_per_base_block = nothing
    in_bin_mean_correlation_per_base_block = nothing
    out_bin_mean_correlation_per_base_block = nothing

    correlation_per_max_relevant_gene = Vector{Float32}(undef, n_genes)

    if bin === nothing
        mean_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
    else
        if base_daf isa DataAxesFormats.Contracts.ContractDaf
            daf = base_daf.daf
        else
            daf = base_daf
        end
        is_relevant_per_gene = BitVector(undef, n_genes)
        is_in_bin_per_gene = get_vector(daf, "gene", "bin").array .== bin
        @assert any(is_in_bin_per_gene)
        @assert !all(is_in_bin_per_gene)
        in_bin_mean_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
        out_bin_mean_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
    end

    n_neighborhood_cells_per_block = get_vector(base_daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    is_in_neighborhood_per_cell = BitVector(undef, n_cells)

    other_metacell_index_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_UMIs_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_other_metacell_UMIs_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    total_punctuated_metacell_UMIs_per_max_neighborhood_cell = Vector{UInt32}(undef, max_n_neighborhood_cells)
    indices_of_max_neighborhood_cells = Vector{Int}(undef, max_n_neighborhood_cells)

    cell_log_fraction_per_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    punctuated_metacell_log_fraction_per_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]

    progress = DebugProgress(
        n_base_blocks * n_included_genes;
        group = :mcs_loops,
        desc = "correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block",
    )

    todox_all_zero = Atomic{Int32}(0)
    for base_block_index in 1:n_base_blocks
        @views is_in_neighborhood_per_other_block =
            is_in_neighborhood_per_other_block_per_base_block[:, base_block_index]
        is_in_neighborhood_per_cell .=
            (other_metacell_index_per_cell .> 0) .& (base_block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(base_block_index_per_cell, 1))

        n_neighborhood_cells = 0
        @foreach_true_index is_in_neighborhood_per_cell cell_index begin
            n_neighborhood_cells += 1
            indices_of_max_neighborhood_cells[n_neighborhood_cells] = cell_index
        end
        @assert n_neighborhood_cells <= n_neighborhood_cells_per_block[base_block_index]
        @views indices_of_neighborhood_cells = indices_of_max_neighborhood_cells[1:n_neighborhood_cells]

        @views other_metacell_index_per_neighborhood_cell =
            other_metacell_index_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_UMIs_per_neighborhood_cell = total_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_other_metacell_UMIs_per_neighborhood_cell =
            total_other_metacell_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]
        @views total_punctuated_metacell_UMIs_per_neighborhood_cell =
            total_punctuated_metacell_UMIs_per_max_neighborhood_cell[1:n_neighborhood_cells]

        other_metacell_index_per_neighborhood_cell .=
            getindex.(Ref(other_metacell_index_per_cell), indices_of_neighborhood_cells)
        total_UMIs_per_neighborhood_cell .= getindex.(Ref(total_UMIs_per_cell), indices_of_neighborhood_cells)
        total_other_metacell_UMIs_per_neighborhood_cell .=
            getindex.(Ref(total_UMIs_per_other_metacell), other_metacell_index_per_neighborhood_cell)
        @. total_punctuated_metacell_UMIs_per_neighborhood_cell =
            total_other_metacell_UMIs_per_neighborhood_cell .- total_UMIs_per_neighborhood_cell

        parallel_loop_wo_rng(1:n_included_genes; progress, policy = :static, progress_chunk = 100) do included_gene_position
            gene_index = indices_of_included_genes[included_gene_position]

            cell_log_fraction_per_max_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell_per_thread[threadid()]
            punctuated_metacell_log_fraction_per_max_neighborhood_cell =
                punctuated_metacell_log_fraction_per_max_neighborhood_cell_per_thread[threadid()]

            @views cell_log_fraction_per_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell[1:n_neighborhood_cells]
            @views punctuated_metacell_log_fraction_per_neighborhood_cell =
                punctuated_metacell_log_fraction_per_max_neighborhood_cell[1:n_neighborhood_cells]

            all_zero_neighborhood_cell_UMIs = true
            for neighborhood_cell_position in 1:n_neighborhood_cells
                if UMIs_per_cell_per_gene[indices_of_neighborhood_cells[neighborhood_cell_position], gene_index] > 0
                    all_zero_neighborhood_cell_UMIs = false
                    break
                end
            end
            if all_zero_neighborhood_cell_UMIs
                atomic_add!(todox_all_zero, Int32(1))
                correlation_per_gene_per_base_block[gene_index, base_block_index] = 0.0f0
                return nothing
            end

            all_zero_neighborhood_punctuated_metacell_UMIs = true
            for neighborhood_cell_position in 1:n_neighborhood_cells
                cell_index = indices_of_neighborhood_cells[neighborhood_cell_position]
                other_metacell_index = other_metacell_index_per_neighborhood_cell[neighborhood_cell_position]
                if UMIs_per_other_metacell_per_gene[other_metacell_index, gene_index] > UMIs_per_cell_per_gene[cell_index, gene_index]
                    all_zero_neighborhood_punctuated_metacell_UMIs = false
                    break
                end
            end
            if all_zero_neighborhood_punctuated_metacell_UMIs
                atomic_add!(todox_all_zero, Int32(1))
                correlation_per_gene_per_base_block[gene_index, base_block_index] = 0.0f0
                return nothing
            end

            for (neighborhood_cell_position, cell_index) in enumerate(indices_of_neighborhood_cells)
                other_metacell_index = other_metacell_index_per_neighborhood_cell[neighborhood_cell_position]
                UMIs_in_cell = Float32(UMIs_per_cell_per_gene[cell_index, gene_index])
                UMIs_in_other_metacell = Float32(UMIs_per_other_metacell_per_gene[other_metacell_index, gene_index])
                total_UMIs_in_cell = Float32(total_UMIs_per_neighborhood_cell[neighborhood_cell_position])
                total_UMIs_in_punctuated_metacell =
                    Float32(total_punctuated_metacell_UMIs_per_neighborhood_cell[neighborhood_cell_position])
                cell_log_fraction_per_neighborhood_cell[neighborhood_cell_position] =
                    log2(UMIs_in_cell / total_UMIs_in_cell + gene_fraction_regularization)
                punctuated_metacell_log_fraction_per_neighborhood_cell[neighborhood_cell_position] = log2(
                    (UMIs_in_other_metacell - UMIs_in_cell) / total_UMIs_in_punctuated_metacell +
                    gene_fraction_regularization,
                )
            end

            cell_log_fraction_per_neighborhood_cell .-= mean(cell_log_fraction_per_neighborhood_cell)  # NOLINT
            punctuated_metacell_log_fraction_per_neighborhood_cell .-=
                mean(punctuated_metacell_log_fraction_per_neighborhood_cell)  # NOLINT

            cell_log_fraction_per_neighborhood_cell ./= std(cell_log_fraction_per_neighborhood_cell)  # NOLINT
            punctuated_metacell_log_fraction_per_neighborhood_cell ./=
                std(punctuated_metacell_log_fraction_per_neighborhood_cell)  # NOLINT

            correlation =
                dot(cell_log_fraction_per_neighborhood_cell, punctuated_metacell_log_fraction_per_neighborhood_cell)
            correlation /= n_neighborhood_cells - 1
            if isnan(correlation)
                correlation = 0.0
            end

            correlation_per_gene_per_base_block[gene_index, base_block_index] = correlation

            return nothing
        end

        @views is_base_neighborhood_marker_per_gene =
            is_base_neighborhood_marker_per_gene_per_base_block[:, base_block_index]
        @. is_base_neighborhood_pertinent_marker_per_gene = is_base_neighborhood_marker_per_gene & !is_lateral_per_gene

        if mean_correlation_per_base_block !== nothing
            n_relevant_genes = sum(is_base_neighborhood_pertinent_marker_per_gene)
            @views correlation_per_relevant_gene = correlation_per_max_relevant_gene[1:n_relevant_genes]
            correlation_per_relevant_gene .=
                correlation_per_gene_per_base_block[is_base_neighborhood_pertinent_marker_per_gene, base_block_index]
            mean_correlation_per_base_block[base_block_index] = mean(correlation_per_relevant_gene)  # NOLINT

        else
            @. is_relevant_per_gene = is_base_neighborhood_pertinent_marker_per_gene & is_in_bin_per_gene  # NOJET
            n_relevant_genes = sum(is_relevant_per_gene)  # NOJET
            @views correlation_per_relevant_gene = correlation_per_max_relevant_gene[1:n_relevant_genes]
            correlation_per_relevant_gene .= correlation_per_gene_per_base_block[is_relevant_per_gene, base_block_index]
            in_bin_mean_correlation_per_base_block[base_block_index] = mean(correlation_per_relevant_gene)  # NOJET # NOLINT

            @. is_relevant_per_gene = is_base_neighborhood_pertinent_marker_per_gene & !is_in_bin_per_gene
            n_relevant_genes = sum(is_relevant_per_gene)
            @views correlation_per_relevant_gene = correlation_per_max_relevant_gene[1:n_relevant_genes]
            correlation_per_relevant_gene .= correlation_per_gene_per_base_block[is_relevant_per_gene, base_block_index]
            out_bin_mean_correlation_per_base_block[base_block_index] = mean(correlation_per_relevant_gene)  # NOJET # NOLINT
        end
    end
    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(n_base_blocks * n_included_genes) ($(percent(todox_all_zero[], n_base_blocks * n_included_genes)))" _group = :todox

    set_matrix!(
        other_daf,
        "gene",
        "base_block",
        "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
        bestify(correlation_per_gene_per_base_block);
        overwrite,
    )

    if mean_correlation_per_base_block !== nothing
        @debug (
            "Mean correlation of base neighborhood pertinent marker genes between base neighborhood cells and their punctuated metacells: " *
            "$(mean(mean_correlation_per_base_block))"  # NOLINT
        ) _group = :mcs_results
    end
    if in_bin_mean_correlation_per_base_block !== nothing
        @debug (
            "Mean correlation of base neighborhood bin pertinent marker genes between base neighborhood cells and their punctuated metacells: " *
            "$(mean(in_bin_mean_correlation_per_base_block))"  # NOLINT
        ) _group = :mcs_results
    end
    if out_bin_mean_correlation_per_base_block !== nothing
        @debug (
            "Mean correlation of base neighborhood !bin pertinent marker genes between base neighborhood cells and their punctuated metacells: " *
            "$(mean(out_bin_mean_correlation_per_base_block))"  # NOLINT
        ) _group = :mcs_results
    end

    return nothing
end

"""
    compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block`](@ref). This takes a long time. Still better than
doing full two-way cross validation, though, and only needed to be done once (for some base metacells).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(RequiredInput),
        matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(CreatedOutput),
    ],
) function compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= gene_fraction_regularization <= 1

    n_cells = axis_length(daf, "cell")
    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    n_blocks = axis_length(daf, "block")

    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    block_index_per_cell = daf["@ cell : metacell ?? 0 : block : index"].array

    is_neighborhood_marker_per_gene_per_block = get_matrix(daf, "gene", "block", "is_neighborhood_marker").array
    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array

    n_neighborhood_markers_per_block = vec(sum(is_neighborhood_marker_per_gene_per_block; dims = 1))
    @assert_vector(n_neighborhood_markers_per_block, n_blocks)
    max_n_neighborhood_markers = maximum(n_neighborhood_markers_per_block)

    correlation_between_max_neighborhood_markers =
        Matrix{Float32}(undef, max_n_neighborhood_markers, max_n_neighborhood_markers)

    n_neighborhood_cells_per_block = get_vector(daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    most_correlated_gene_per_gene_per_block = Matrix{AbstractString}(undef, n_genes, n_blocks)
    fill!(most_correlated_gene_per_gene_per_block, "")

    total_pertinent_neighborhood_marker_genes = sum(is_neighborhood_marker_per_gene_per_block[.!is_lateral_per_gene, :])

    is_in_neighborhood_per_cell = BitVector(undef, n_cells)
    is_pertinent_neighborhood_marker_per_gene = BitVector(undef, n_genes)
    indices_of_max_neighborhood_cells = Vector{Int}(undef, max_n_neighborhood_cells)
    indices_of_max_neighborhood_markers = Vector{Int}(undef, max_n_neighborhood_markers)

    log_fraction_per_max_neighborhood_cell_per_max_neighborhood_marker =
        Matrix{Float32}(undef, max_n_neighborhood_cells, max_n_neighborhood_markers)

    progress = DebugProgress(
        3 * total_pertinent_neighborhood_marker_genes - 1;
        group = :mcs_loops,
        desc = "most_correlated_pertinent_neighborhood_markers_per_gene_per_block",
    )

    for block_index in 1:n_blocks
        @views most_correlated_gene_per_gene = most_correlated_gene_per_gene_per_block[:, block_index]

        @views is_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene_per_block[:, block_index]
        @. is_pertinent_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene & !is_lateral_per_gene

        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_cell .=
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1))

        n_pertinent_neighborhood_markers = 0
        @foreach_true_index is_pertinent_neighborhood_marker_per_gene gene_index begin
            n_pertinent_neighborhood_markers += 1
            indices_of_max_neighborhood_markers[n_pertinent_neighborhood_markers] = gene_index
        end
        @views indices_of_pertinent_neighborhood_markers =
            indices_of_max_neighborhood_markers[1:n_pertinent_neighborhood_markers]

        n_neighborhood_cells = 0
        @foreach_true_index is_in_neighborhood_per_cell cell_index begin
            n_neighborhood_cells += 1
            indices_of_max_neighborhood_cells[n_neighborhood_cells] = cell_index
        end
        @assert n_neighborhood_cells == n_neighborhood_cells_per_block[block_index]
        @views indices_of_neighborhood_cells = indices_of_max_neighborhood_cells[1:n_neighborhood_cells]

        @views log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker =
            log_fraction_per_max_neighborhood_cell_per_max_neighborhood_marker[
                1:n_neighborhood_cells,
                1:n_pertinent_neighborhood_markers,
            ]

        @views correlation_between_pertinent_neighborhood_markers = correlation_between_max_neighborhood_markers[
            1:n_pertinent_neighborhood_markers,
            1:n_pertinent_neighborhood_markers,
        ]
        correlation_between_pertinent_neighborhood_markers[1, 1] = 0

        parallel_loop_wo_rng(
            1:n_pertinent_neighborhood_markers;
            name = "compute_log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker",
            progress,
            progress_chunk = 100,
        ) do pertinent_neighborhood_marker_position
            gene_index = indices_of_pertinent_neighborhood_markers[pertinent_neighborhood_marker_position]
            for (neighborhood_cell_position, cell_index) in enumerate(indices_of_neighborhood_cells)
                log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker[
                    neighborhood_cell_position,
                    pertinent_neighborhood_marker_position,
                ] = log2(
                    UMIs_per_cell_per_gene[cell_index, gene_index] ./ total_UMIs_per_cell[cell_index] +
                    gene_fraction_regularization,
                )
            end
        end

        parallel_loop_wo_rng(
            1:(n_pertinent_neighborhood_markers - 1);
            name = "compute_correlation_between_pertinent_neighborhood_markers",
            progress,
            progress_chunk = 100,
        ) do base_position_index
            base_position = n_pertinent_neighborhood_markers + 1 - base_position_index
            @assert 2 <= base_position <= n_pertinent_neighborhood_markers
            correlation_between_pertinent_neighborhood_markers[base_position, base_position] = 0
            @views base_column = log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker[:, base_position]
            for other_position in 1:(base_position - 1)
                @views other_column =
                    log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker[:, other_position]
                correlation = zero_cor_between_vectors(base_column, other_column)
                correlation_between_pertinent_neighborhood_markers[base_position, other_position] = correlation
                correlation_between_pertinent_neighborhood_markers[other_position, base_position] = correlation
            end
            return nothing
        end

        parallel_loop_wo_rng(
            1:n_pertinent_neighborhood_markers;
            name = "compute_most_correlated_pertinent_neighborhood_markers",
            progress,
            progress_chunk = 100,
        ) do pertinent_neighborhood_marker_position
            gene_index = indices_of_pertinent_neighborhood_markers[pertinent_neighborhood_marker_position]
            @views correlation_with_pertinent_neighborhood_marker =
                correlation_between_pertinent_neighborhood_markers[:, pertinent_neighborhood_marker_position]
            most_correlated_pertinent_neighborhood_marker_position =
                argmax(correlation_with_pertinent_neighborhood_marker)
            most_correlated_gene_index =
                indices_of_pertinent_neighborhood_markers[most_correlated_pertinent_neighborhood_marker_position]
            most_correlated_gene_per_gene[gene_index] = name_per_gene[most_correlated_gene_index]
            return nothing
        end
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "most_correlated_gene_in_neighborhood",
        most_correlated_gene_per_gene_per_block;
        overwrite,
    )

    return nothing
end

"""
    compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!(;
        other_daf::DafWriter,
        base_daf::DafReader;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set
[`matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block`](@ref). This
will also copy [`base_block_axis`](@ref) from the `base_daf` into the `other_daf` if needed.

The `cell` and `gene` axes must be identical in both `base_daf` and `other_daf`.

That is, we evaluate the quality of the `other_daf` metacells by correlating the genes in each cell in each neighborhood
of `base_daf` with the expression of their most correlated gene in their `other_daf` metacells. This works as a poor
man's cross-validation; not as good as real cross-validation but much less expensive to compute.

# Other

$(CONTRACT1)

# Base

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(
    name = "other_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        gene_axis(RequiredInput),
        base_block_axis(GuaranteedOutput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block(
            CreatedOutput,
        ),
    ],
) Contract(
    name = "base_daf",
    is_relaxed = true,
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(RequiredInput),
    ],
) function compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!(;
    other_daf::DafWriter,
    base_daf::DafReader,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= gene_fraction_regularization <= 1

    @assert axis_vector(base_daf, "gene") == axis_vector(other_daf, "gene")
    @assert axis_vector(base_daf, "cell") == axis_vector(other_daf, "cell")

    n_cells = axis_length(base_daf, "cell")
    n_genes = axis_length(base_daf, "gene")

    linear_fraction_per_other_metacell_per_gene = get_matrix(other_daf, "metacell", "gene", "linear_fraction").array
    other_metacell_index_per_cell = other_daf["@ cell : metacell ?? 0 : index"].array

    n_base_blocks = axis_length(base_daf, "block")

    if has_axis(other_daf, "base_block")
        @assert axis_vector(other_daf, "base_block") == axis_vector(base_daf, "block")
    else
        copy_axis!(source = base_daf, destination = other_daf, axis = "block", rename = "base_block"; overwrite)
    end

    UMIs_per_cell_per_gene = get_matrix(base_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(base_daf, "cell", "total_UMIs").array

    base_block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array

    is_in_base_neighborhood_per_other_block_per_base_block =
        get_matrix(base_daf, "block", "block", "is_in_neighborhood").array

    most_correlated_gene_in_base_neighborhood_per_gene_per_base_block =
        get_matrix(base_daf, "gene", "block", "most_correlated_gene_in_neighborhood").array
    n_most_correlated_gene_per_base_block =
        vec(sum(most_correlated_gene_in_base_neighborhood_per_gene_per_base_block .!= ""; dims = 1))
    @assert_vector(n_most_correlated_gene_per_base_block, n_base_blocks)
    n_max_correlated_genes = maximum(n_most_correlated_gene_per_base_block)

    correlation_with_most_per_gene_per_base_block = zeros(Float32, n_genes, n_base_blocks)

    mean_correlation_with_most_per_base_block = Vector{Float32}(undef, n_base_blocks)

    total_correlated_genes = sum(most_correlated_gene_in_base_neighborhood_per_gene_per_base_block .!= "")

    is_correlated_per_gene = BitVector(undef, n_genes)
    is_in_base_neighborhood_per_cell = BitVector(undef, n_cells)

    n_neighborhood_cells_per_block = get_vector(base_daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)
    indices_of_max_correlated_genes = Vector{Int}(undef, n_max_correlated_genes)
    indices_of_max_base_neighborhood_cells = Vector{Int}(undef, max_n_neighborhood_cells)

    cell_log_fraction_per_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    other_metacell_log_fraction_per_other_max_neighborhood_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    correlation_per_max_correlated_gene_per_thread =
        [Vector{Float32}(undef, n_max_correlated_genes) for _ in 1:nthreads()]

    progress = DebugProgress(
        total_correlated_genes,
        group = :mcs_loops,
        desc = "correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block",
    )

    todox_all_zero = Atomic{Int32}(0)
    for base_block_index in 1:n_base_blocks
        @views most_correlated_gene_in_base_neighborhood_per_gene =
            most_correlated_gene_in_base_neighborhood_per_gene_per_base_block[:, base_block_index]
        @. is_correlated_per_gene = most_correlated_gene_in_base_neighborhood_per_gene != ""

        @views is_in_base_neighborhood_per_other_block =
            is_in_base_neighborhood_per_other_block_per_base_block[:, base_block_index]
        is_in_base_neighborhood_per_cell .=
            (other_metacell_index_per_cell .> 0) .& (base_block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_base_neighborhood_per_other_block), max.(base_block_index_per_cell, 1))

        n_correlated_genes = 0
        @foreach_true_index is_correlated_per_gene gene_index begin
            n_correlated_genes += 1
            indices_of_max_correlated_genes[n_correlated_genes] = gene_index
        end
        @views indices_of_correlated_genes = indices_of_max_correlated_genes[1:n_correlated_genes]

        n_base_neighborhood_cells = 0
        @foreach_true_index is_in_base_neighborhood_per_cell cell_index begin
            n_base_neighborhood_cells += 1
            indices_of_max_base_neighborhood_cells[n_base_neighborhood_cells] = cell_index
        end
        @views indices_of_base_neighborhood_cells = indices_of_max_base_neighborhood_cells[1:n_base_neighborhood_cells]

        parallel_loop_wo_rng(1:n_correlated_genes; progress, policy = :static, progress_chunk = 100) do correlated_gene_position
            gene_index = indices_of_correlated_genes[correlated_gene_position]

            cell_log_fraction_per_max_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell_per_thread[threadid()]
            other_metacell_log_fraction_per_other_max_neighborhood_cell =
                other_metacell_log_fraction_per_other_max_neighborhood_cell_per_thread[threadid()]

            @views cell_log_fraction_per_base_neighborhood_cell =
                cell_log_fraction_per_max_neighborhood_cell[1:n_base_neighborhood_cells]
            @views other_metacell_log_fraction_per_other_base_neighborhood_cell =
                other_metacell_log_fraction_per_other_max_neighborhood_cell[1:n_base_neighborhood_cells]

            all_zero_base_neighborhood_cell_UMIs = true
            for base_neighborhood_cell_position in 1:n_base_neighborhood_cells
                if UMIs_per_cell_per_gene[indices_of_base_neighborhood_cells[base_neighborhood_cell_position], gene_index] > 0
                    all_zero_base_neighborhood_cell_UMIs = false
                    break
                end
            end
            if all_zero_base_neighborhood_cell_UMIs
                atomic_add!(todox_all_zero, Int32(1))
                correlation_with_most_per_gene_per_base_block[gene_index, base_block_index] = 0.0f0
                return nothing
            end

            all_zero_base_neighborhood_other_metacell_fractions = true
            for base_neighborhood_cell_position in 1:n_base_neighborhood_cells
                cell_index = indices_of_base_neighborhood_cells[base_neighborhood_cell_position]
                other_metacell_index = other_metacell_index_per_cell[cell_index]
                if linear_fraction_per_other_metacell_per_gene[other_metacell_index, gene_index] > 0
                    all_zero_base_neighborhood_other_metacell_fractions = false
                    break
                end
            end
            if all_zero_base_neighborhood_other_metacell_fractions
                atomic_add!(todox_all_zero, Int32(1))
                correlation_with_most_per_gene_per_base_block[gene_index, base_block_index] = 0.0f0
                return nothing
            end

            for (base_neighborhood_cell_position, cell_index) in enumerate(indices_of_base_neighborhood_cells)
                other_metacell_index = other_metacell_index_per_cell[cell_index]
                cell_log_fraction_per_base_neighborhood_cell[base_neighborhood_cell_position] = log2(
                    UMIs_per_cell_per_gene[cell_index, gene_index] / total_UMIs_per_cell[cell_index] +
                    gene_fraction_regularization,
                )
                other_metacell_log_fraction_per_other_base_neighborhood_cell[base_neighborhood_cell_position] =
                    log2(linear_fraction_per_other_metacell_per_gene[other_metacell_index, gene_index] + gene_fraction_regularization)
            end

            return correlation_with_most_per_gene_per_base_block[gene_index, base_block_index] =
                zero_cor_between_vectors(
                    cell_log_fraction_per_base_neighborhood_cell,
                    other_metacell_log_fraction_per_other_base_neighborhood_cell,
                )
        end

        correlation_per_max_correlated_gene = correlation_per_max_correlated_gene_per_thread[threadid()]
        @views correlation_per_relevant_gene = correlation_per_max_correlated_gene[1:n_correlated_genes]
        correlation_per_relevant_gene .=
            correlation_with_most_per_gene_per_base_block[is_correlated_per_gene, base_block_index]
        mean_correlation_with_most = mean(correlation_per_relevant_gene)  # NOLINT
        @assert mean_correlation_with_most > 0
        mean_correlation_with_most_per_base_block[base_block_index] = mean_correlation_with_most
    end
    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(total_correlated_genes) ($(percent(todox_all_zero[], total_correlated_genes)))" _group = :todox

    set_matrix!(
        other_daf,
        "gene",
        "base_block",
        "correlation_with_most_between_base_neighborhood_cells_and_metacells",
        bestify(correlation_with_most_per_gene_per_base_block);
        overwrite,
    )

    @debug (
        "Mean correlation of base neighborhood pertinent marker genes with friends between base neighborhood cells and their metacells:" *
        " $(mean(mean_correlation_with_most_per_base_block))"  # NOLINT
    ) _group = :mcs_results

    return nothing
end

"""
    compute_vector_of_type_per_block_by_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_block`](@ref) by using the type of the metacells grouped into each block. The most
frequent ("mode") of the type is used for each block. This method is used when metacells are annotated with types (which
is typically a supervised process).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_metacell(RequiredInput),
        vector_of_type_per_block(CreatedOutput),
    ],
) function compute_vector_of_type_per_block_by_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_block = daf["@ metacell : type / block =@ >> Mode"].array
    set_vector!(daf, "block", "type", type_per_block; overwrite)
    return nothing
end

"""
    compute_vector_of_type_per_block_by_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_block`](@ref) by using the type of the metacells grouped into each block. The most
frequent ("mode") of the type is used for each block. This method is used when cells are annotated with types (typically
when importing single-cell data with type annotations).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_metacell(RequiredInput),
        vector_of_type_per_block(CreatedOutput),
    ],
) function compute_vector_of_type_per_block_by_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_block = daf["@ cell : metacell ?? : type / block =@ >> Mode"].array
    set_vector!(daf, "block", "type", type_per_block; overwrite)
    return nothing
end

"""
    compute_vector_of_type_per_cell_by_blocks!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_cell`](@ref) by using the type of the blocks the metacells of the cells belong to.
This method is used when blocks are annotated with types (which is typically a supervised process).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_block(RequiredInput),
        vector_of_type_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_type_per_cell_by_blocks!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_cell = daf["@ cell : metacell ?? '' : block : type"].array
    set_vector!(daf, "cell", "type", type_per_cell; overwrite)
    return nothing
end

"""
    compute_vector_of_type_per_metacell_by_blocks!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_metacell`](@ref) by using the type of the blocks the metacells belong to. This
method is used when blocks are annotated with types (which is typically a supervised process).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_block(RequiredInput),
        vector_of_type_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_type_per_metacell_by_blocks!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_metacell = daf["@ metacell : block : type"].array
    set_vector!(daf, "metacell", "type", type_per_metacell; overwrite)
    return nothing
end

end  # module
