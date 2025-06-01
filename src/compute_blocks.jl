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
using Distances
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
import Metacells.Contracts.block_component_base_covered_fraction_matrix
import Metacells.Contracts.block_component_gene_covered_coefficient_tensor
import Metacells.Contracts.block_component_is_found_matrix
import Metacells.Contracts.block_component_n_genes_matrix
import Metacells.Contracts.block_components_neighborhood_RMSE_vector
import Metacells.Contracts.block_components_neighborhood_XRMSE_vector
import Metacells.Contracts.block_covered_UMIs_vector
import Metacells.Contracts.block_gene_base_covered_fraction_matrix
import Metacells.Contracts.block_gene_component_matrix
import Metacells.Contracts.block_gene_environment_marker_rank_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_is_environment_skeleton_matrix
import Metacells.Contracts.block_n_found_components_vector
import Metacells.Contracts.block_n_metacells_vector
import Metacells.Contracts.component_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance

@kwdef struct ConfidenceIntervals
    low_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_blocks!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        max_block_span::Real = $(DEFAULT.max_block_span),
        num_blocks::Integer = $(DEFAULT.num_blocks),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group the metacells into blocks where the metacells in each one are "similar" in all the skeleton genes. That is, each
block is an approximation of a single cell state, assuming the skeleton genes adequately predict the rest of the genes.

First, we compute metacell-metacell distances based on the maximal fold factor between `is_skeleton` genes. The fold
factor is log (base 2) of the gene expression using the `gene_fraction_regularization`. For computing this fold factor,
we ignore "insignificant" genes whose total UMIs in the compared metacells isn't at least `min_significant_gene_UMIs`.
We also we reduce the distance using the `fold_confidence` based on the number of UMIs used to estimate the expression
in the metacells. Two metacells can only belong to the same block if the final fold factor is at most `max_block_span`
in all the (significant, skeleton) genes.

We then compute block-block distances which are the maximal distance between the metacells of the blocks.

!!! note

    This uses the virtual [`metacell_gene_covered_fraction_matrix`](@ref). You will need an `adapter` to map these to
    concrete fractions (geomean, linear, scaled, ...).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        metacell_gene_covered_fraction_matrix(OptionalInput),
        metacell_gene_log_covered_fraction_matrix(OptionalInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_covered_UMIs_vector(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        metacell_metacell_max_skeleton_fold_distance(OptionalOutput),
        metacell_metacell_euclidean_skeleton_distance(OptionalOutput),
        block_block_max_skeleton_fold_distance(OptionalOutput),
        block_block_mean_euclidean_skeleton_distance(OptionalOutput),
    ],
) function compute_blocks!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    max_block_span::Real = function_default(identify_marker_genes!, :min_marker_gene_range_fold),
    num_blocks::Integer = 431,  # TODOX
    overwrite::Bool = false,
    metacells_distance_method = :euclidean, # TODOX
    block_distance_aggregation::Function = mean,  # TODOX
    block_hclust_linkage::Symbol = :ward,  # TODOX
)::Nothing
    @assert 0 < fold_confidence < 1
    @assert max_block_span > 0

    n_metacells = axis_length(daf, "metacell")

    UMIs_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : UMIs"].array
    covered_UMIs_per_metacell = get_vector(daf, "metacell", "covered_UMIs").array

    if metacells_distance_method == :euclidean
        log_covered_fraction_per_skeleton_per_metacell =
            daf["/ gene & is_skeleton / metacell : log_covered_fraction"].array
        distances_between_metacells = pairwise(Euclidean(), log_covered_fraction_per_skeleton_per_metacell)
        set_matrix!(daf, "metacell", "metacell", "euclidean_skeleton_distance", distances_between_metacells)
    else
        @assert metacells_distance_method == :max_fold
        covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : covered_fraction"].array
        confidence_intervals = compute_confidence_intervals(;
            gene_fraction_regularization,
            covered_fraction_per_skeleton_per_metacell,
            covered_UMIs_per_metacell,
            fold_confidence,
        )

        distances_between_metacells = compute_maximal_distances_between_metacells(;
            confidence_intervals,
            UMIs_per_skeleton_per_metacell,
            min_significant_gene_UMIs,
        )
        set_matrix!(daf, "metacell", "metacell", "max_skeleton_fold_distance", distances_between_metacells)
    end

    clusters = hclust(distances_between_metacells; linkage = block_hclust_linkage)  # NOJET
    if block_hclust_linkage == :complete
        block_index_per_metacell = cutree(clusters; h = max_block_span)
    else
        block_index_per_metacell = cutree(clusters; k = num_blocks)
    end
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    n_blocks = length(metacell_indices_per_block)
    name_per_block = group_names(daf, "metacell", metacell_indices_per_block; prefix = "B")

    @info "Blocks: $(n_blocks) Mean metacells: $(n_metacells / n_blocks)"
    add_axis!(daf, "block", name_per_block)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)

    distances_between_blocks = compute_distances_between_blocks(;
        distances_between_metacells,
        metacell_indices_per_block,
        block_distance_aggregation,
        max_block_span,
    )
    if metacells_distance_method == :euclidean
        set_matrix!(daf, "block", "block", "mean_euclidean_skeleton_distance", distances_between_blocks)
    else
        set_matrix!(daf, "block", "block", "max_skeleton_fold_distance", distances_between_blocks)
    end
    return nothing
end

function compute_confidence_intervals(;  # UNTESTED
    gene_fraction_regularization::Real,
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    covered_UMIs_per_metacell::AbstractVector{<:Integer},
    fold_confidence::AbstractFloat,
)::ConfidenceIntervals
    n_skeletons, n_metacells = size(covered_fraction_per_skeleton_per_metacell)
    @assert_matrix(covered_fraction_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert_vector(covered_UMIs_per_metacell, n_metacells)

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

    return ConfidenceIntervals(;
        low_log_covered_fraction_per_skeleton_per_metacell,
        high_log_covered_fraction_per_skeleton_per_metacell,
    )
end

function compute_euclidean_distances_between_metacells(;  # UNTESTED
    confidence_intervals::ConfidenceIntervals,
    UMIs_per_skeleton_per_metacell::AbstractMatrix{<:Integer},
    min_significant_gene_UMIs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)
    @assert_matrix(UMIs_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
end

function compute_maximal_distances_between_metacells(;  # UNTESTED
    confidence_intervals::ConfidenceIntervals,
    UMIs_per_skeleton_per_metacell::AbstractMatrix{<:Integer},
    min_significant_gene_UMIs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)
    @assert_matrix(UMIs_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert_matrix(
        confidence_intervals.high_log_covered_fraction_per_skeleton_per_metacell,
        n_skeletons,
        n_metacells,
        Columns
    )
    @assert_matrix(
        confidence_intervals.low_log_covered_fraction_per_skeleton_per_metacell,
        n_skeletons,
        n_metacells,
        Columns
    )

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)
    distances_between_metacells[1, 1] = 0.0

    @threads :greedy for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views base_metacell_UMIs_per_skeleton = vec(UMIs_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_low_log_covered_fraction_per_skeleton =
            vec(confidence_intervals.low_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_high_log_covered_fraction_per_skeleton =
            vec(confidence_intervals.high_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])

        n_other_metacells = base_metacell_index - 1
        other_metacells_indices = 1:(base_metacell_index - 1)
        @views UMIs_per_skeleton_per_other_metacells = UMIs_per_skeleton_per_metacell[:, other_metacells_indices]
        @views low_log_covered_fraction_per_skeleton_per_other_metacells =
            confidence_intervals.low_log_covered_fraction_per_skeleton_per_metacell[:, other_metacells_indices]
        @views high_log_covered_fraction_per_skeleton_per_other_metacells =
            confidence_intervals.high_log_covered_fraction_per_skeleton_per_metacell[:, other_metacells_indices]

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

    return distances_between_metacells
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

function compute_distances_between_blocks(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    block_distance_aggregation::Function,
    max_block_span::Real,
)::AbstractMatrix{<:AbstractFloat}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacell_indices_per_block)
    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    @threads :greedy for base_block_index in reverse(1:n_blocks)
        base_metacells_indices = metacell_indices_per_block[base_block_index]

        @views distance_per_metacell_per_base_metacell = distances_between_metacells[:, base_metacells_indices]

        aggregated_distance_from_base_per_metacell =
            vec(block_distance_aggregation(distance_per_metacell_per_base_metacell; dims = 2))
        @assert length(aggregated_distance_from_base_per_metacell) == n_metacells

        for other_block_index in 1:base_block_index
            other_metacells_indices = metacell_indices_per_block[other_block_index]

            aggregated_distance_from_base_per_other_metacell =
                aggregated_distance_from_base_per_metacell[other_metacells_indices]

            aggregated_distance_between_base_and_other_block =
                block_distance_aggregation(aggregated_distance_from_base_per_other_metacell)

            distances_between_blocks[base_block_index, other_block_index] =
                aggregated_distance_between_base_and_other_block
            distances_between_blocks[other_block_index, base_block_index] =
                aggregated_distance_between_base_and_other_block
        end

        # TODOX @assert distances_between_blocks[base_block_index, base_block_index] <= max_block_span
        if base_block_index > 1
            @views distances_between_base_and_other_blocks =
                distances_between_blocks[1:(base_block_index - 1), base_block_index]
            # TODOX @assert all(distances_between_base_and_other_blocks .> max_block_span)
        end
    end

    return distances_between_blocks
end

"""
    function compute_blocks_is_in_neighborhood!(
        daf::DafWriter;
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        min_covered_UMIs_in_neighborhood::Integer = $(DEFAULT.min_covered_UMIs_in_neighborhood),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

Group the blocks into overlapping small tight neighborhoods, and larger environments, such that each block's environment
can be reasonably approximated using a local linear model, evaluated using cross-validation on the RMSE in the
neighborhood at its core.

Using the block-block distances, we define a tight neighborhood of each block containing at least
`min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` and `min_covered_UMIs_in_neighborhood` (note that the
latter only counts UMIs of skeleton genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        block_covered_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_max_skeleton_fold_distance(OptionalInput),
        block_block_mean_euclidean_skeleton_distance(OptionalInput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_neighborhood!(
    daf::DafWriter;
    min_blocks_in_neighborhood::Integer = 5,
    min_metacells_in_neighborhood::Integer = 20,
    min_covered_UMIs_in_neighborhood::Integer = 2_000_000,
    overwrite::Bool = false,
)::Nothing
    compute_blocks_is_in_vicinity!(
        daf::DafWriter;
        name = "neighborhood",
        min_blocks_in_vicinity = min_blocks_in_neighborhood,
        min_metacells_in_vicinity = min_metacells_in_neighborhood,
        min_covered_UMIs_in_vicinity = min_covered_UMIs_in_neighborhood,
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
`min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` and `min_covered_UMIs_in_neighborhood` (note that the
latter only counts UMIs of skeleton genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        block_covered_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_max_skeleton_fold_distance(OptionalInput),
        block_block_mean_euclidean_skeleton_distance(OptionalInput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_environment!(
    daf::DafWriter;
    min_blocks_in_environment::Integer = function_default(compute_blocks_is_in_neighborhood!, :min_blocks_in_neighborhood) * 2,
    min_metacells_in_environment::Integer = function_default(compute_blocks_is_in_neighborhood!, :min_metacells_in_neighborhood) * 2,
    min_covered_UMIs_in_environment::Integer = function_default(compute_blocks_is_in_neighborhood!, :min_covered_UMIs_in_neighborhood) * 2,
    overwrite::Bool = false,
)::Nothing
    compute_blocks_is_in_vicinity!(
        daf::DafWriter;
        name = "environment",
        min_blocks_in_vicinity = min_blocks_in_environment,
        min_metacells_in_vicinity = min_metacells_in_environment,
        min_covered_UMIs_in_vicinity = min_covered_UMIs_in_environment,
        overwrite,
    )
    return nothing
end

function compute_blocks_is_in_vicinity!(
    daf::DafWriter;
    min_blocks_in_vicinity::Integer,
    min_metacells_in_vicinity::Integer,
    min_covered_UMIs_in_vicinity::Integer,
    name::AbstractString,
    overwrite::Bool = false,
)::Nothing
    @assert min_blocks_in_vicinity > 0
    @assert min_covered_UMIs_in_vicinity >= 0

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    if has_matrix(daf, "block", "block", "max_skeleton_fold_distance")
        max_skeleton_fold_distance_per_block_per_block = get_matrix(daf, "block", "block", "max_skeleton_fold_distance")
    else
        max_skeleton_fold_distance_per_block_per_block = get_matrix(daf, "block", "block", "mean_euclidean_skeleton_distance")
    end
    n_metacells_per_block = get_vector(daf, "block", "n_metacells")
    covered_UMIs_per_block = get_vector(daf, "block", "covered_UMIs")
    is_in_vicinity_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)

    min_blocks_in_vicinity = min(min_blocks_in_vicinity, n_blocks)

    @threads :greedy for base_block_index in 1:n_blocks
        @views distance_from_base_per_block = max_skeleton_fold_distance_per_block_per_block[:, base_block_index]

        ordered_block_indices = sortperm(distance_from_base_per_block)
        @assert ordered_block_indices[1] == base_block_index

        vicinity_n_blocks = 0
        vicinity_n_metacells = 0
        vicinity_covered_UMIs = 0

        while vicinity_n_blocks < min_blocks_in_vicinity ||
                  vicinity_n_metacells < min_metacells_in_vicinity ||
                  vicinity_covered_UMIs < min_covered_UMIs_in_vicinity
            vicinity_n_blocks += 1
            next_block_index = ordered_block_indices[vicinity_n_blocks]
            is_in_vicinity_per_other_block_per_base_block[base_block_index, next_block_index] = true
            vicinity_n_metacells += n_metacells_per_block[next_block_index]
            vicinity_covered_UMIs += covered_UMIs_per_block[next_block_index]
        end

        @debug (
            "$(uppercasefirst(name)): $(name_per_block[base_block_index])" *
            " Blocks: $(vicinity_n_blocks)" *
            " Metacells: $(vicinity_n_metacells)" *
            " Covered M-UMIs: $(vicinity_covered_UMIs / 1e6)"
        )
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

end  # component
