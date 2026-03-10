"""
Group "very close" metacells into disjoint blocks.
"""
module ComputeBlocks

export compute_metacells_blocks!

using Base.Threads
using Clustering
using DataAxesFormats
using Distributions
using LinearAlgebra
using MultivariateStats
using Random
using StatsBase
using TanayLabUtilities

using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_block_per_metacell

@kwdef struct ConfidenceIntervals
    low_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_metacells_blocks!(
        daf::DafWriter;
        num_blocks::Integer = $(DEFAULT.num_blocks),
        max_block_span::Real = $(DEFAULT.max_block_span),
        prefix::AbstractString = $(DEFAULT.prefix),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set the [`block_axis`](@ref) and [`vector_of_block_per_metacell`](@ref). We group the metacells into blocks
where the metacells in each one are "similar" in all the skeleton genes. That is, each block is an approximation of a
single cell state, assuming the skeleton genes adequately predict the rest of the genes.

If `num_blocks` is not specified, we decide on the number of blocks by using hierarchical clustering based on the
maximal skeleton fold distance between the metacells in each cluster of `max_block_span` (that is, `hclust` with
`:complete` linkage using the [`matrix_of_max_skeleton_fold_distance_between_metacells`](@ref).

The actual blocks are computed using hierarchical clustering with `:ward` linkage based on the euclidean skeleton
distance between the metacells (that is, using the [`matrix_of_euclidean_skeleton_fold_distance_between_metacells`](@ref)).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(CreatedOutput)],
    data = [
        matrix_of_max_skeleton_fold_distance_between_metacells(RequiredInput),
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(RequiredInput),
        vector_of_block_per_metacell(CreatedOutput),
    ],
) function compute_metacells_blocks!(
    daf::DafWriter;
    num_blocks::Maybe{Integer} = nothing,
    max_block_span::Real = function_default(compute_vector_of_is_marker_per_gene!, :min_marker_gene_range_fold),
    prefix::AbstractString = "B",
    overwrite::Bool = false,
)::Nothing
    if num_blocks === nothing
        distances_between_metacells = get_matrix(daf, "metacell", "metacell", "max_skeleton_fold_distance").array
        clusters = hclust(distances_between_metacells; linkage = :complete)
        block_index_per_metacell = cutree(clusters; h = max_block_span)
        num_blocks = maximum(block_index_per_metacell)
    end

    distances_between_metacells = get_matrix(daf, "metacell", "metacell", "euclidean_skeleton_fold_distance").array
    clusters = hclust(distances_between_metacells; linkage = :ward)  # NOJET
    block_index_per_metacell = cutree(clusters; k = num_blocks)

    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    n_blocks = length(metacell_indices_per_block)
    name_per_block = group_names(axis_vector(daf, "metacell"), metacell_indices_per_block; prefix)

    add_axis!(daf, "block", name_per_block)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)

    @debug "Blocks: $(n_blocks)" _group = :mcs_results

    return nothing
end

end
