"""
Group "very close" metacells into disjoint blocks, and blocks into overlapping vicinities, for constructing a local
linear approximation to the manifold.
"""
module ComputeBlocks

export compute_blocks!
export compute_blocks_is_in_neighborhood!
export compute_blocks_is_in_environment!
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
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_block_max_skeleton_fold_distance
import Metacells.Contracts.block_block_mean_euclidean_skeleton_distance_matrix
import Metacells.Contracts.block_total_UMIs_vector
import Metacells.Contracts.block_n_metacells_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_metacell_euclidean_skeleton_distance
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance

@kwdef struct ConfidenceIntervals
    low_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_blocks!(
        daf::DafWriter;
        max_block_span::Real = $(DEFAULT.max_block_span),
        num_blocks::Integer = $(DEFAULT.num_blocks),
        prefix::AbstractString = $(DEFAULT.prefix),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group the metacells into blocks where the metacells in each one are "similar" in all the skeleton genes. That is, each
block is an approximation of a single cell state, assuming the skeleton genes adequately predict the rest of the genes.

TODOX

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        metacell_metacell_max_skeleton_fold_distance(OptionalInput),
        metacell_metacell_euclidean_skeleton_distance(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
    ],
) function compute_blocks!(
    daf::DafWriter;
    max_block_span::Real = function_default(identify_marker_genes!, :min_marker_gene_range_fold),
    num_blocks::Maybe{Integer} = nothing,
    prefix::AbstractString = "B",
    overwrite::Bool = false,
)::Nothing
    if num_blocks === nothing
        distances_between_metacells = get_matrix(daf, "metacell", "metacell", "max_skeleton_fold_distance").array
        clusters = hclust(distances_between_metacells; linkage = :complete)
        block_index_per_metacell = cutree(clusters; h = max_block_span)
        num_blocks = maximum(block_index_per_metacell)
        @debug "num_blocks: $(num_blocks)"
    end

    distances_between_metacells = get_matrix(daf, "metacell", "metacell", "euclidean_skeleton_distance").array
    clusters = hclust(distances_between_metacells; linkage = :ward)  # NOJET
    block_index_per_metacell = cutree(clusters; k = num_blocks)

    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    n_blocks = length(metacell_indices_per_block)
    name_per_block = group_names(daf, "metacell", metacell_indices_per_block; prefix)

    @info "Blocks: $(n_blocks) Mean metacells: $(axis_length(daf, "metacell") / n_blocks)"
    add_axis!(daf, "block", name_per_block)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)
    return nothing
end

"""
    function compute_blocks_is_in_neighborhood!(
        daf::DafWriter;
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        min_total_UMIs_in_neighborhood::Integer = $(DEFAULT.min_total_UMIs_in_neighborhood),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

Group the blocks into overlapping small tight neighborhoods, and larger environments, such that each block's environment
can be reasonably approximated using a local linear model, evaluated using cross-validation on the RMSE in the
neighborhood at its core.

Using the block-block distances, we define a tight neighborhood of each block containing at least
`min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` and `min_total_UMIs_in_neighborhood` (note that the
latter only counts UMIs of skeleton genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        block_total_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_mean_euclidean_skeleton_distance_matrix(OptionalInput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_neighborhood!(
    daf::DafWriter;
    min_blocks_in_neighborhood::Integer = 5,
    min_metacells_in_neighborhood::Integer = 20,
    min_total_UMIs_in_neighborhood::Integer = 2_000_000,
    overwrite::Bool = false,
)::Nothing
    @assert min_blocks_in_neighborhood > 0
    @assert min_metacells_in_neighborhood > 0
    @assert min_total_UMIs_in_neighborhood > 0

    compute_blocks_is_in_vicinity!(
        daf::DafWriter;
        name = "neighborhood",
        min_blocks_in_vicinity = min_blocks_in_neighborhood,
        min_metacells_in_vicinity = min_metacells_in_neighborhood,
        min_total_UMIs_in_vicinity = min_total_UMIs_in_neighborhood,
        overwrite,
    )
    return nothing
end

"""
    function compute_blocks_is_in_environment!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

Group the blocks into overlapping small tight neighborhoods, and larger environments, such that each block's environment
can be reasonably approximated using a local linear model, evaluated using cross-validation on the RMSE in the
neighborhood at its core.

Using the block-block distances, we define a tight neighborhood of each block containing at least
`min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` and `min_total_UMIs_in_neighborhood` (note that the
latter only counts UMIs of skeleton genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        block_total_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_mean_euclidean_skeleton_distance_matrix(OptionalInput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_environment!(
    daf::DafWriter;
    min_blocks_in_environment::Integer = function_default(
        compute_blocks_is_in_neighborhood!,
        :min_blocks_in_neighborhood,
    ) * 2,
    min_metacells_in_environment::Integer = function_default(
        compute_blocks_is_in_neighborhood!,
        :min_metacells_in_neighborhood,
    ) * 2,
    min_total_UMIs_in_environment::Integer = function_default(
        compute_blocks_is_in_neighborhood!,
        :min_total_UMIs_in_neighborhood,
    ) * 2,
    overwrite::Bool = false,
)::Nothing
    @assert min_blocks_in_environment > 0
    @assert min_metacells_in_environment > 0
    @assert min_total_UMIs_in_environment > 0

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    compute_blocks_is_in_vicinity!(
        daf::DafWriter;
        name = "environment",
        min_blocks_in_vicinity = min_blocks_in_environment,
        min_metacells_in_vicinity = min_metacells_in_environment,
        min_total_UMIs_in_vicinity = min_total_UMIs_in_environment,
        overwrite,
    )
    return nothing
end

function compute_blocks_is_in_vicinity!(
    daf::DafWriter;
    min_blocks_in_vicinity::Integer,
    min_metacells_in_vicinity::Integer,
    min_total_UMIs_in_vicinity::Integer,
    name::AbstractString,
    overwrite::Bool = false,
)::Nothing
    @assert min_blocks_in_vicinity > 0
    @assert min_total_UMIs_in_vicinity >= 0

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    mean_euclidean_skeleton_distance_per_block_per_block =
        get_matrix(daf, "block", "block", "mean_euclidean_skeleton_distance")
    n_metacells_per_block = get_vector(daf, "block", "n_metacells")
    total_UMIs_per_block = get_vector(daf, "block", "total_UMIs")
    is_in_vicinity_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)

    min_blocks_in_vicinity = min(min_blocks_in_vicinity, n_blocks)

    parallel_loop_wo_rng(1:n_blocks) do base_block_index
        @views distance_from_base_per_block = mean_euclidean_skeleton_distance_per_block_per_block[:, base_block_index]
        ordered_block_indices = sortperm(distance_from_base_per_block)
        if ordered_block_indices[1] != base_block_index
            distance_from_base_per_block = copy_array(distance_from_base_per_block)
            distance_from_base_per_block[base_block_index] = 0
            ordered_block_indices = sortperm(distance_from_base_per_block)
            @assert ordered_block_indices[1] == base_block_index
        end

        vicinity_n_blocks = 0
        vicinity_n_metacells = 0
        vicinity_total_UMIs = 0

        while vicinity_n_blocks < min_blocks_in_vicinity ||
                  vicinity_n_metacells < min_metacells_in_vicinity ||
                  vicinity_total_UMIs < min_total_UMIs_in_vicinity
            vicinity_n_blocks += 1
            next_block_index = ordered_block_indices[vicinity_n_blocks]
            distance_from_base = distance_from_base_per_block[next_block_index]
            is_in_vicinity_per_other_block_per_base_block[next_block_index, base_block_index] = true
            vicinity_n_metacells += n_metacells_per_block[next_block_index]
            vicinity_total_UMIs += total_UMIs_per_block[next_block_index]
        end

        @debug (
            "$(uppercasefirst(name)): $(name_per_block[base_block_index])" *
            " Blocks: $(vicinity_n_blocks)" *
            " Metacells: $(vicinity_n_metacells)" *
            " Covered M-UMIs: $(vicinity_total_UMIs / 1e6)"
        )
        return nothing
    end

    set_matrix!(  # NOJET
        daf,
        "block",
        "block",
        "is_in_$(name)",
        bestify(is_in_vicinity_per_other_block_per_base_block);
        overwrite,
    )

    return nothing
end

end
