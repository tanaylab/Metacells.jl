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
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_skeleton_per_gene

@kwdef struct ConfidenceIntervals
    low_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_metacells_blocks!(
        daf::DafWriter;
        n_blocks::Integer = $(DEFAULT.n_blocks),
        max_block_span::Real = $(DEFAULT.max_block_span),
        prefix::AbstractString = $(DEFAULT.prefix),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set the [`block_axis`](@ref) and [`vector_of_block_per_metacell`](@ref). We group the metacells into blocks
where the metacells in each one are "similar" in all the skeleton genes. That is, each block is an approximation of a
single cell state, assuming the skeleton genes adequately predict the rest of the genes.

If `n_blocks` is not specified, we decide on the number of blocks by using hierarchical clustering based on the
maximal skeleton fold distance between the metacells in each cluster of `max_block_span` (that is, `hclust` with
`:complete` linkage using the [`matrix_of_max_skeleton_fold_distance_between_metacells`](@ref).

The actual blocks are computed using K-means on the log of the skeleton genes fractions of the metacells, using the same
regularization as [`matrix_of_euclidean_skeleton_fold_distance_between_metacells`](@ref). If `n_blocks` is not
specified, K-means starts from the centroids of the hierarchical clustering blocks. Otherwise, it starts from `:kmpp`
seeds chosen using the `rng`. Therefore, specifying `n_blocks` gives different blocks than the automatic method, even
if the number of blocks is the same. Blocks which end up empty are dropped.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), block_axis(CreatedOutput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_max_skeleton_fold_distance_between_metacells(RequiredInput),
        vector_of_block_per_metacell(CreatedOutput),
    ],
) function compute_metacells_blocks!(  # UNTESTED
    daf::DafWriter;
    n_blocks::Maybe{Integer} = nothing,
    max_block_span::Real = function_default(compute_vector_of_is_marker_per_gene!, :min_marker_gene_range_fold),
    prefix::AbstractString = "B",
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    linear_fraction_per_metacell_per_skeleton = plain_array(daf["@ metacell @ gene [ is_skeleton ] :: linear_fraction"])
    log_linear_fraction_per_skeleton_per_metacell =
        log2.(flipped(linear_fraction_per_metacell_per_skeleton) .+ GENE_FRACTION_REGULARIZATION_FOR_CELLS)

    if n_blocks === nothing
        distances_between_metacells = plain_array(get_matrix(daf, "metacell", "metacell", "max_skeleton_fold_distance"))
        clusters = hclust(distances_between_metacells; linkage = :complete)
        seed_block_index_per_metacell = cutree(clusters; h = max_block_span)
        mean_log_linear_fraction_per_skeleton_per_block = reduce(
            hcat,
            [
                vec(mean(view(log_linear_fraction_per_skeleton_per_metacell, :, metacell_indices); dims = 2)) for  # NOLINT
                metacell_indices in collect_group_members(seed_block_index_per_metacell)
            ],
        )
    else
        seed_metacell_index_per_block = initseeds(:kmpp, log_linear_fraction_per_skeleton_per_metacell, n_blocks; rng)
        mean_log_linear_fraction_per_skeleton_per_block =
            log_linear_fraction_per_skeleton_per_metacell[:, seed_metacell_index_per_block]
    end

    block_index_per_metacell = assignments(
        kmeans!(log_linear_fraction_per_skeleton_per_metacell, mean_log_linear_fraction_per_skeleton_per_block; rng),
    )
    n_blocks = compact_groups!(block_index_per_metacell)
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    name_per_block = group_names(axis_vector(daf, "metacell"), metacell_indices_per_block; prefix)

    add_axis!(daf, "block", name_per_block; overwrite)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)

    @debug "Blocks: $(n_blocks)" _group = :mcs_results

    return nothing
end

end
