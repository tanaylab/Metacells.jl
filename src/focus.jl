"""
Reduce a large set of metacells to a smaller, non-redundant focus set, for analyzing large data sets.
"""
module Focus

export compute_focus_metacells!
export copy_focus_metacells_data!

using DataAxesFormats
using Statistics
using TanayLabUtilities

using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..ComputeBlocks
using ..Contracts
using ..Defaults

import ..AnalyzeMetacells.nearest_neighbors_per_point
import ..ComputeBlocks.compute_confidence_intervals
import ..ComputeBlocks.confident_gene_distance

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell

"""
    compute_focus_metacells!(;
        full_daf::DafReader,
        focus_daf::DafWriter,
        max_redundant_fold::Real = $(DEFAULT.max_redundant_fold),
        n_neighbors::Integer = $(DEFAULT.n_neighbors),
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        max_fold_gene_fraction_regularization::AbstractFloat = $(DEFAULT.max_fold_gene_fraction_regularization),
        max_fold_min_significant_gene_UMIs::Integer = $(DEFAULT.max_fold_min_significant_gene_UMIs),
        max_fold_confidence::AbstractFloat = $(DEFAULT.max_fold_confidence),
        max_fold_min_confidence_UMIs::Real = $(DEFAULT.max_fold_min_confidence_UMIs),
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        marker_min_significant_gene_UMIs::Integer = $(DEFAULT.marker_min_significant_gene_UMIs),
        marker_gene_fraction_regularization::AbstractFloat = $(DEFAULT.marker_gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Choose the focus metacells of the `full_daf`. Create the focus [`metacell_axis`](@ref) of the `focus_daf` with the
focus metacells, and its [`cell_axis`](@ref) with their cells. Create its [`gene_axis`](@ref) with all the genes,
unless it already has one, e.g. from a base repository of gene masks. Set its [`vector_of_is_skeleton_per_gene`](@ref)
to the skeleton genes that are still markers of the focus metacells.

This works in passes, until a pass drops no metacell. Each pass:

  - Keeps only the skeleton genes that are still markers of the remaining metacells, using the
    `min_marker_gene_max_fraction`, `min_marker_gene_range_fold`, `min_marker_quantile`,
    `marker_min_significant_gene_UMIs` and `marker_gene_fraction_regularization` (see
    [`compute_vector_of_is_marker_per_gene!`](@ref)).
  - Finds the `n_neighbors` nearest neighbors of each remaining metacell, using the Euclidean distance between the log
    of the fractions of the skeleton genes, adding the `gene_fraction_regularization` (see
    [`compute_matrix_of_knn_euclidean_skeleton_fold_distance_per_metacell_per_metacell!`](@ref)). Unlike other fold
    factors, this does not require each gene to be based on at least `max_fold_min_significant_gene_UMIs`. Testing all
    the pairs would be too expensive for a large number of metacells. Instead, the high `gene_fraction_regularization`
    reduces the effect of genes with few UMIs. This only chooses the pairs to test. The test itself does use the floor.
  - Computes the maximal confident fold of the skeleton genes between each metacell and each of its neighbors, using
    the `max_fold_gene_fraction_regularization`, `max_fold_min_significant_gene_UMIs`, `max_fold_confidence` and
    `max_fold_min_confidence_UMIs` (see [`compute_metacells_blocks!`](@ref)).
  - Considers the pairs whose maximal fold is at most `max_redundant_fold`, in ascending order of the fold. If both
    metacells of a pair are still there, drops the one with fewer cells. On a tie, it drops the one with fewer UMIs,
    and on a tie, the one with the higher index.

This only creates the axes and the skeleton genes mask of the `focus_daf`. To populate the rest of the data, use
[`copy_focus_metacells_data!`](@ref).

# Full

$(CONTRACT1)

# Focus

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(
    name = "full_daf",
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_n_cells_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
    ],
) Contract(
    name = "focus_daf",
    axes = [metacell_axis(CreatedOutput), cell_axis(CreatedOutput), gene_axis(GuaranteedOutput)],
    data = [vector_of_is_skeleton_per_gene(CreatedOutput)],
) function compute_focus_metacells!(;  # UNTESTED
    full_daf::DafReader,
    focus_daf::DafWriter,
    max_redundant_fold::Real = 0.5,
    n_neighbors::Integer = function_default(
        compute_matrix_of_knn_euclidean_skeleton_fold_distance_per_metacell_per_metacell!,
        :n_neighbors,
    ),
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_matrix_of_knn_euclidean_skeleton_fold_distance_per_metacell_per_metacell!,
        :gene_fraction_regularization,
    ),
    max_fold_gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_blocks!,
        :max_fold_gene_fraction_regularization,
    ),
    max_fold_min_significant_gene_UMIs::Integer = function_default(
        compute_metacells_blocks!,
        :max_fold_min_significant_gene_UMIs,
    ),
    max_fold_confidence::AbstractFloat = function_default(compute_metacells_blocks!, :max_fold_confidence),
    max_fold_min_confidence_UMIs::Real = function_default(compute_metacells_blocks!, :max_fold_min_confidence_UMIs),
    min_marker_gene_max_fraction::AbstractFloat = function_default(
        compute_vector_of_is_marker_per_gene!,
        :min_marker_gene_max_fraction,
    ),
    min_marker_gene_range_fold::Real = function_default(
        compute_vector_of_is_marker_per_gene!,
        :min_marker_gene_range_fold,
    ),
    min_marker_quantile::Real = function_default(compute_vector_of_is_marker_per_gene!, :min_marker_quantile),
    marker_min_significant_gene_UMIs::Integer = function_default(
        compute_vector_of_is_marker_per_gene!,
        :min_significant_gene_UMIs,
    ),
    marker_gene_fraction_regularization::AbstractFloat = function_default(
        compute_vector_of_is_marker_per_gene!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert max_redundant_fold >= 0
    @assert 0 < n_neighbors
    @assert 0 <= min_marker_quantile < 1
    @assert marker_min_significant_gene_UMIs >= 0
    @assert marker_gene_fraction_regularization > 0

    is_skeleton_per_gene = Vector{Bool}(get_vector(full_daf, "gene", "is_skeleton").array)
    indices_of_skeletons = findall(is_skeleton_per_gene)
    n_skeletons = length(indices_of_skeletons)

    linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(full_daf["@ gene [ is_skeleton ] @ metacell :: linear_fraction"].array)
    log_linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(full_daf["@ gene [ is_skeleton ] @ metacell :: log_linear_fraction"].array)
    UMIs_per_skeleton_per_metacell = Matrix{UInt32}(full_daf["@ gene [ is_skeleton ] @ metacell :: UMIs"].array)
    total_UMIs_per_metacell = plain_array(get_vector(full_daf, "metacell", "total_UMIs"))
    n_cells_per_metacell = plain_array(get_vector(full_daf, "metacell", "n_cells"))
    n_metacells = length(n_cells_per_metacell)

    confidence_intervals = compute_confidence_intervals(;
        linear_fraction_per_skeleton_per_metacell,
        total_UMIs_per_metacell,
        gene_fraction_regularization = max_fold_gene_fraction_regularization,
        fold_confidence = max_fold_confidence,
        min_confidence_UMIs = max_fold_min_confidence_UMIs,
    )

    is_focus_per_metacell = trues(n_metacells)
    is_focus_per_skeleton = trues(n_skeletons)
    pass = 0
    while true
        pass += 1
        indices_of_focus_metacells = findall(is_focus_per_metacell)

        # Keep only the skeleton genes that are still markers of the remaining metacells.
        for skeleton_position in findall(is_focus_per_skeleton)
            @views linear_fraction_per_focus_metacell =
                linear_fraction_per_skeleton_per_metacell[skeleton_position, indices_of_focus_metacells]
            @views log_linear_fraction_per_focus_metacell =
                log_linear_fraction_per_skeleton_per_metacell[skeleton_position, indices_of_focus_metacells]
            if min_marker_quantile == 0
                log_minimal_fraction = minimum(log_linear_fraction_per_focus_metacell)
            else
                log_minimal_fraction = quantile(log_linear_fraction_per_focus_metacell, min_marker_quantile)
            end
            minimal_fraction = 2^log_minimal_fraction - marker_gene_fraction_regularization

            # An insignificant fold counts as zero.
            range_fold = 0.0f0
            for metacell_index in indices_of_focus_metacells
                if UMIs_per_skeleton_per_metacell[skeleton_position, metacell_index] +
                   minimal_fraction * total_UMIs_per_metacell[metacell_index] >= marker_min_significant_gene_UMIs
                    range_fold = max(
                        range_fold,
                        log_linear_fraction_per_skeleton_per_metacell[skeleton_position, metacell_index] -
                        log_minimal_fraction,
                    )
                end
            end

            is_focus_per_skeleton[skeleton_position] =
                maximum(linear_fraction_per_focus_metacell) >= min_marker_gene_max_fraction &&
                range_fold >= min_marker_gene_range_fold
        end
        focus_skeleton_positions = findall(is_focus_per_skeleton)
        @assert !isempty(focus_skeleton_positions)

        # The maximal confident fold between each remaining metacell and each of its nearest neighbors.
        log_linear_fraction_per_focus_skeleton_per_focus_metacell = log2.(
            linear_fraction_per_skeleton_per_metacell[focus_skeleton_positions, indices_of_focus_metacells] .+
            Float32(gene_fraction_regularization),
        )
        neighbor_position_per_edge, _ =
            nearest_neighbors_per_point(log_linear_fraction_per_focus_skeleton_per_focus_metacell, n_neighbors)
        max_fold_per_edge = Vector{Float32}(undef, length(neighbor_position_per_edge))
        parallel_loop_wo_rng(
            eachindex(indices_of_focus_metacells);
            name = "focus_max_fold_per_edge",
            progress = DebugProgress(
                length(indices_of_focus_metacells);
                group = :mcs_loops,
                desc = "focus_max_fold_per_edge",
            ),
            progress_chunk = 100,
        ) do focus_position
            metacell_index = indices_of_focus_metacells[focus_position]
            for edge_index in ((focus_position - 1) * n_neighbors + 1):(focus_position * n_neighbors)
                neighbor_index = indices_of_focus_metacells[neighbor_position_per_edge[edge_index]]
                max_fold = 0.0f0
                for skeleton_position in focus_skeleton_positions
                    max_fold = max(
                        max_fold,
                        confident_gene_distance(
                            max_fold_min_significant_gene_UMIs,
                            UMIs_per_skeleton_per_metacell[skeleton_position, metacell_index],
                            confidence_intervals.low_log_linear_fraction_per_skeleton_per_metacell[
                                skeleton_position,
                                metacell_index,
                            ],
                            confidence_intervals.high_log_linear_fraction_per_skeleton_per_metacell[
                                skeleton_position,
                                metacell_index,
                            ],
                            UMIs_per_skeleton_per_metacell[skeleton_position, neighbor_index],
                            confidence_intervals.low_log_linear_fraction_per_skeleton_per_metacell[
                                skeleton_position,
                                neighbor_index,
                            ],
                            confidence_intervals.high_log_linear_fraction_per_skeleton_per_metacell[
                                skeleton_position,
                                neighbor_index,
                            ],
                        ),
                    )
                end
                max_fold_per_edge[edge_index] = max_fold
            end
            return nothing
        end

        # Drop the smaller metacell of each redundant pair, in ascending order of the fold.
        redundant_pairs = Tuple{Float32, Int, Int}[]
        for (edge_index, max_fold) in enumerate(max_fold_per_edge)
            if max_fold <= max_redundant_fold
                metacell_index = indices_of_focus_metacells[div(edge_index - 1, n_neighbors) + 1]
                neighbor_index = indices_of_focus_metacells[neighbor_position_per_edge[edge_index]]
                push!(
                    redundant_pairs,
                    (max_fold, min(metacell_index, neighbor_index), max(metacell_index, neighbor_index)),
                )
            end
        end
        sort!(unique!(redundant_pairs))

        n_dropped = 0
        for (_, first_index, second_index) in redundant_pairs
            if is_focus_per_metacell[first_index] && is_focus_per_metacell[second_index]
                is_focus_per_metacell[dropped_metacell_index(
                    first_index,
                    second_index,
                    n_cells_per_metacell,
                    total_UMIs_per_metacell,
                )] = false
                n_dropped += 1
            end
        end

        @debug (
            "Focus pass: $(pass) metacells: $(length(indices_of_focus_metacells)) dropped: $(n_dropped)" *
            " skeletons: $(length(focus_skeleton_positions))"
        ) _group = :mcs_results

        if n_dropped == 0
            break
        end
    end

    is_skeleton_per_gene[indices_of_skeletons[.!is_focus_per_skeleton]] .= false

    name_per_metacell = axis_vector(full_daf, "metacell")
    metacell_index_per_cell = full_daf["@ cell : metacell ?? 0 : index"].array
    is_focus_per_cell =
        [metacell_index > 0 && is_focus_per_metacell[metacell_index] for metacell_index in metacell_index_per_cell]

    add_axis!(focus_daf, "metacell", name_per_metacell[is_focus_per_metacell]; overwrite)
    add_axis!(focus_daf, "cell", axis_vector(full_daf, "cell")[is_focus_per_cell]; overwrite)
    if has_axis(focus_daf, "gene")
        @assert axis_vector(focus_daf, "gene") == axis_vector(full_daf, "gene")
    else
        add_axis!(focus_daf, "gene", axis_vector(full_daf, "gene"))
    end
    set_vector!(focus_daf, "gene", "is_skeleton", is_skeleton_per_gene; overwrite)
    return nothing
end

"""
    copy_focus_metacells_data!(;
        full_daf::DafReader,
        focus_daf::DafWriter,
        excluded_contract::Contract,
    )::Nothing

Copy the data of the `full_daf` into the `focus_daf`, restricted to the focus metacells and cells that
[`compute_focus_metacells!`](@ref) chose. Data the `focus_daf` already has is kept. This includes the focus skeleton
genes mask, and the gene masks if the `focus_daf` rests on their repository.

The outputs of the `excluded_contract` are not copied. These are typically the outputs of the computations that are
redone for the focus metacells, because they depend on which metacells there are (e.g., the marker genes). Both layouts
of an excluded matrix are excluded.

# Full

$(CONTRACT1)

# Focus

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(
    name = "full_daf",
    is_relaxed = true,
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput), gene_axis(RequiredInput)],
) Contract(
    name = "focus_daf",
    is_relaxed = true,
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput), gene_axis(RequiredInput)],
) function copy_focus_metacells_data!(;  # UNTESTED
    full_daf::DafReader,
    focus_daf::DafWriter,
    excluded_contract::Contract,
)::Nothing
    view_data = Pair[VIEW_ALL_DATA...]
    for (key, (expectation, _, _)) in something(excluded_contract.data, [])
        if expectation in (CreatedOutput, GuaranteedOutput, OptionalOutput)
            push!(view_data, key => nothing)
            if key isa Tuple && length(key) == 3
                rows_axis, columns_axis, name = key
                push!(view_data, (columns_axis, rows_axis, name) => nothing)
            end
        end
    end
    copy_all!(; source = viewer(full_daf; data = view_data), destination = focus_daf, insist = false)
    return nothing
end

# Which of two redundant metacells to drop: the one with fewer cells, then the one with fewer UMIs, then the one with
# the higher index.
function dropped_metacell_index(  # UNTESTED
    first_index::Integer,
    second_index::Integer,
    n_cells_per_metacell::AbstractVector{<:Integer},
    total_UMIs_per_metacell::AbstractVector{<:Real},
)::Integer
    if n_cells_per_metacell[first_index] != n_cells_per_metacell[second_index]
        return n_cells_per_metacell[first_index] < n_cells_per_metacell[second_index] ? first_index : second_index
    elseif total_UMIs_per_metacell[first_index] != total_UMIs_per_metacell[second_index]
        return total_UMIs_per_metacell[first_index] < total_UMIs_per_metacell[second_index] ? first_index : second_index
    else
        return max(first_index, second_index)
    end
end

end  # module
