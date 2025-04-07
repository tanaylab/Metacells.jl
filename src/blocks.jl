"""
Approximate the manifold of actual cell states (captured by metacells) using linear programs in each local region.
"""
module Blocks

export compute_covered_fractions!
export compute_blocks!
export sharpen_metacells!

using Base.Threads
using Clustering
using DataAxesFormats
using Distributions
using MultivariateStats
using Random
using SparseArrays
using StatsBase
using TanayLabUtilities

using ..Contracts
using ..Defaults
using ..IdentifyGenes

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_gene_mean_scaled_log_covered_fraction_matrix
import Metacells.Contracts.block_linear_RMSE_vector
import Metacells.Contracts.block_linear_XRMSE_vector
import Metacells.Contracts.block_n_principal_components_vector
import Metacells.Contracts.block_principal_component_gene_covered_coefficient_tensor
import Metacells.Contracts.block_principal_component_gene_skeleton_coefficient_tensor
import Metacells.Contracts.block_principal_component_is_used_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_new_metacell_vector
import Metacells.Contracts.cell_old_metacell_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_fraction_matrix
import Metacells.Contracts.metacell_gene_scaled_log_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_total_UMIs_matrix
import Metacells.Contracts.new_metacell_axis
import Metacells.Contracts.new_metacell_block_vector
import Metacells.Contracts.old_metacell_axis
import Metacells.Contracts.old_metacell_block_vector
import Metacells.Contracts.principal_component_axis

"""
    function compute_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute normalized fractions for each metacell based only on the covered genes. By changing the denominator to only this
smaller set of genes, we remove effects of possibly high-expression "irrelevant" gene programs (e.g., cell cycle). In
general whenever we talk about gene "fractions" we have to keep in mind what is the total the fractions are out of, and
ensure it is fit for purpose (and, of course, ensure the same divergence factors are used for all the data). Here we are
going to be using these fractions for creating a distance measure in the manifold, which will be used for both deciding
what is a "local" region and for estimating errors of reconsctructing the manifold, taking into consideration only the
covered genes.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        gene_divergence_vector(OptionalInput),
        metacell_gene_total_UMIs_matrix(RequiredInput),
        metacell_gene_fraction_matrix(RequiredInput),
        metacell_gene_covered_fraction_matrix(GuaranteedOutput),
        metacell_gene_scaled_log_covered_fraction_matrix(GuaranteedOutput),
        metacell_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0

    adapter(  # NOJET
        daf;
        input_axes = ["metacell" => "=", "gene" => "/ gene & is_covered"],
        input_data = [
            ("gene", "divergence") => "=",
            ("gene", "metacell", "total_UMIs") => "=",
            ("gene", "metacell", "fraction") => "=",
        ],
        output_axes = ["gene" => "=", "metacell" => "="],
        output_data = [
            ("metacell", "gene", "covered_fraction") => "=",
            ("metacell", "gene", "scaled_log_covered_fraction") => "=",
            ("metacell", "covered_UMIs") => "=",
        ],
        empty = Dict([
            ("metacell", "gene", "covered_fraction") => 0,
            ("metacell", "gene", "scaled_log_covered_fraction") => 0,
        ]),
        overwrite,
    ) do adapted
        n_metacells = axis_length(adapted, "metacell")

        covered_UMIs_per_metacell = adapted["/ gene / metacell : total_UMIs %> Sum"].array
        @assert_vector(covered_UMIs_per_metacell, n_metacells)
        set_vector!(adapted, "metacell", "covered_UMIs", covered_UMIs_per_metacell)

        fraction_per_covered_per_metacell = get_matrix(adapted, "gene", "metacell", "fraction").array
        covered_fraction_per_metacell = vec(sum(fraction_per_covered_per_metacell; dims = 1))
        @assert_vector(covered_fraction_per_metacell, n_metacells)

        covered_fraction_per_gene_per_metacell =
            fraction_per_covered_per_metacell .* transpose(1 ./ covered_fraction_per_metacell)

        scaled_log_covered_fraction_per_gene_per_metacell =
            log2.(covered_fraction_per_gene_per_metacell .+ gene_fraction_regularization)

        if has_vector(adapted, "gene", "divergence")
            divergence_per_covered = get_vector(adapted, "gene", "divergence").array
            scaled_log_covered_fraction_per_gene_per_metacell .*= (1 .- divergence_per_covered)
        end

        set_matrix!(
            adapted,
            "gene",
            "metacell",
            "covered_fraction",
            covered_fraction_per_gene_per_metacell;
            eltype = Float32,
        )

        return set_matrix!(
            adapted,
            "gene",
            "metacell",
            "scaled_log_covered_fraction",
            scaled_log_covered_fraction_per_gene_per_metacell;
            eltype = Float32,
        )
    end

    return nothing
end

@kwdef struct BlocksData
    n_blocks::Integer
    metacell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}}
    distances_between_blocks::AbstractMatrix{<:AbstractFloat}
    covered_UMIs_per_metacell::AbstractVector{<:Integer}
end

@kwdef struct ConfidenceIntervals
    low_scaled_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_scaled_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

@kwdef mutable struct LocalModel
    n_principal_components::Integer
    #
    mean_scaled_log_covered_fraction_per_skeleton::AbstractVector{<:AbstractFloat}
    coefficient_per_skeleton_per_principal_component::AbstractMatrix{<:AbstractFloat}
    #
    mean_scaled_log_covered_fraction_per_covered::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    coefficient_per_principal_component_per_covered::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    #
    XRMSE::Maybe{AbstractFloat} = nothing
    RMSE::Maybe{AbstractFloat} = nothing
end

@kwdef struct CrossValidationIndices
    train_indices_per_part::AbstractVector{<:AbstractVector{<:Integer}}
    test_indices_per_part::AbstractVector{<:AbstractVector{<:Integer}}
end

"""
    function compute_blocks!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        max_block_span::Real = $(DEFAULT.max_block_span),
        max_principal_components = $(DEFAULT.max_principal_components),
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        min_covered_UMIs_in_neighborhood::Integer = $(DEFAULT.min_covered_UMIs_in_neighborhood),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group the metacells into blocks where the metacells in each one are "similar" in all the skeleton genes. That is, each
block is an approximation of a single cell state, assuming the skeleton genes adequately predict the rest of the genes.
Then, group these blocks into overlapping small tight neighborhoods, and larger environments, such that each block's
environment can be reasonably approximated using a local linear model, evaluated using cross-validation on the RMSE in
the neighborhood at its core.

First, we compute metacell-metacell distances based on the maximal fold factor between `is_skeleton` genes. The fold
factor is log (base 2) of the gene expression using the `gene_fraction_regularization`. For computing this fold factor,
we ignore "insignificant" genes whose total UMIs in the compared metacells isn't at least `min_significant_gene_UMIs`.
We also we reduce the distance using the `fold_confidence` based on the number of UMIs used to estimate the expression
in the metacells, Two metacells can only belong to the same block if the final fold factor is at most `max_block_span`
in all the (significant, skeleton) genes.

We then compute block-block distances which are the mean distance between the metacells of the blocks. Using these
distances, we define a tight neighborhood of each block containing at least `min_blocks_in_neighborhood`,
`min_metacells_in_neighborhood` and `min_covered_UMIs_in_neighborhood` (note that the latter only counts UMIs of
skeleton genes).

Having computed the neighborhoods, we expand them (using the same block-block distances) as long as computing a local
linear model for the environment gives a better approximation of the metacells in the neighborhood. This local model
starts with computing `max_principal_components`. These are assumed to over-fit the solution, so we use cross-validation
(in `cross_validation_parts`) to pick a subset of these principal components, which hopefully approximates the true
dimensionality of the local data.

!!! note

    There is no attempt to reconcile the local model of neighboring blocks. In general the linear model computed here
    isn't very useful for interpreting the data as principal components for noisy high-dimensional data are pretty
    opaque.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(GuaranteedOutput),
        principal_component_axis(GuaranteedOutput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        gene_divergence_vector(RequiredInput),
        metacell_gene_total_UMIs_matrix(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_gene_scaled_log_covered_fraction_matrix(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
        block_n_principal_components_vector(GuaranteedOutput),
        block_principal_component_is_used_matrix(GuaranteedOutput),
        block_gene_mean_scaled_log_covered_fraction_matrix(GuaranteedOutput),
        block_principal_component_gene_skeleton_coefficient_tensor(GuaranteedOutput),
        block_principal_component_gene_covered_coefficient_tensor(GuaranteedOutput),
        block_linear_RMSE_vector(GuaranteedOutput),
        block_linear_XRMSE_vector(GuaranteedOutput),
    ],
) function compute_blocks!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    max_block_span::Real = function_default(identify_marker_genes!, :min_marker_gene_range_fold),
    max_principal_components::Integer = 40,
    min_blocks_in_neighborhood::Integer = 4,
    min_metacells_in_neighborhood::Integer = 20,
    min_covered_UMIs_in_neighborhood::Integer = 2_000_000,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
    block_distance_aggregation::Function = mean,  # TODOX
    leaky_pca_cross_validation::Bool = false,  # TODOX
)::Nothing
    @assert gene_fraction_regularization > 0
    @assert 0 < fold_confidence < 1
    @assert max_block_span > 0
    @assert max_principal_components > 0
    @assert min_blocks_in_neighborhood > 0
    @assert min_covered_UMIs_in_neighborhood >= 0
    @assert cross_validation_parts > 0
    @assert min_metacells_in_neighborhood >= cross_validation_parts

    blocks = compute_coherent_blocks!(
        daf;
        gene_fraction_regularization,
        min_significant_gene_UMIs,
        fold_confidence,
        max_block_span,
        overwrite,
        block_distance_aggregation,
    )

    is_in_neighborhood_per_other_block_per_base_block = compute_neighborhoods!(
        daf;
        min_blocks_in_neighborhood,
        min_metacells_in_neighborhood,
        min_covered_UMIs_in_neighborhood,
        blocks,
        overwrite,
    )

    add_axis!(daf, "principal_component", ["PC$(index)" for index in 1:max_principal_components]; overwrite)

    compute_environments!(
        daf;
        max_principal_components,
        cross_validation_parts,
        is_in_neighborhood_per_other_block_per_base_block,
        blocks,
        rng,
        overwrite,
        leaky_pca_cross_validation,
    )

    return nothing
end

function compute_coherent_blocks!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat,
    min_significant_gene_UMIs::Integer,
    fold_confidence::AbstractFloat,
    max_block_span::Real,
    overwrite::Bool,
    block_distance_aggregation::Function,
)::BlocksData
    n_metacells = axis_length(daf, "metacell")

    covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : covered_fraction"].array
    total_UMIs_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : total_UMIs"].array
    covered_UMIs_per_metacell = daf["/ gene & is_covered / metacell : total_UMIs %> Sum"].array

    divergene_per_skeleton = daf["/ gene & is_skeleton : divergence"].array

    confidence_intervals = compute_confidence_intervals(;
        covered_fraction_per_skeleton_per_metacell,
        covered_UMIs_per_metacell,
        divergence_per_gene = divergene_per_skeleton,
        gene_fraction_regularization,
        fold_confidence,
    )

    distances_between_metacells = compute_distances_between_metacells(;
        confidence_intervals,
        total_UMIs_per_gene_per_metacell = total_UMIs_per_skeleton_per_metacell,
        min_significant_gene_UMIs,
    )

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    block_index_per_metacell = cutree(clusters; h = max_block_span)
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    n_blocks = length(metacell_indices_per_block)
    name_per_block = group_names(daf, "metacell", metacell_indices_per_block; prefix = "B")

    @info (
        "Blocks: $(n_blocks) Mean metacells: $(n_metacells / n_blocks)" *
        " Mean covered M-UMIs: $(sum(covered_UMIs_per_metacell) / (1e6 * n_blocks))"
    )
    add_axis!(daf, "block", name_per_block)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)

    distances_between_blocks = compute_distances_between_blocks(;
        distances_between_metacells,
        metacell_indices_per_block,
        block_distance_aggregation,
    )

    return BlocksData(; n_blocks, metacell_indices_per_block, distances_between_blocks, covered_UMIs_per_metacell)
end

function compute_confidence_intervals(;  # UNTESTED
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    covered_UMIs_per_metacell::AbstractVector{<:Integer},
    divergence_per_gene::AbstractVector{<:AbstractFloat},
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
)::ConfidenceIntervals
    n_genes, n_metacells = size(covered_fraction_per_skeleton_per_metacell)
    @assert_matrix(covered_fraction_per_skeleton_per_metacell, n_genes, n_metacells, Columns)
    @assert_vector(covered_UMIs_per_metacell, n_metacells)
    @assert_vector(divergence_per_gene, n_genes)

    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_fractions_per_gene_per_metacells =  # NOJET
        confidence_stdevs .*
        sqrt.(transpose(covered_UMIs_per_metacell) .* covered_fraction_per_skeleton_per_metacell) ./
        transpose(covered_UMIs_per_metacell)

    low_scaled_log_covered_fraction_per_skeleton_per_metacell =  # NOJET
        log2.(
            max.(covered_fraction_per_skeleton_per_metacell .- confidence_fractions_per_gene_per_metacells, 0.0) .+
            gene_fraction_regularization
        ) .* (1.0 .- divergence_per_gene)

    high_scaled_log_covered_fraction_per_skeleton_per_metacell =
        log2.(
            covered_fraction_per_skeleton_per_metacell .+ confidence_fractions_per_gene_per_metacells .+
            gene_fraction_regularization
        ) .* (1.0 .- divergence_per_gene)

    return ConfidenceIntervals(;
        low_scaled_log_covered_fraction_per_skeleton_per_metacell,
        high_scaled_log_covered_fraction_per_skeleton_per_metacell,
    )
end

function compute_distances_between_metacells(;  # UNTESTED
    confidence_intervals::ConfidenceIntervals,
    total_UMIs_per_gene_per_metacell::AbstractMatrix{<:Integer},
    min_significant_gene_UMIs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_genes, n_metacells = size(confidence_intervals.low_scaled_log_covered_fraction_per_skeleton_per_metacell)
    @assert_matrix(
        confidence_intervals.high_scaled_log_covered_fraction_per_skeleton_per_metacell,
        n_genes,
        n_metacells,
        Columns
    )
    @assert_matrix(
        confidence_intervals.low_scaled_log_covered_fraction_per_skeleton_per_metacell,
        n_genes,
        n_metacells,
        Columns
    )
    @assert_matrix(total_UMIs_per_gene_per_metacell, n_genes, n_metacells, Columns)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)

    distances_between_metacells[1, 1] = 0.0

    @threads :greedy for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views base_metacell_total_UMIs_per_gene = vec(total_UMIs_per_gene_per_metacell[:, base_metacell_index])
        @views base_metacell_low_scaled_log_covered_fraction_per_gene =
            vec(confidence_intervals.low_scaled_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_high_scaled_log_covered_fraction_per_gene =
            vec(confidence_intervals.high_scaled_log_covered_fraction_per_skeleton_per_metacell[:, base_metacell_index])

        n_other_metacells = base_metacell_index - 1
        index_per_other_metacell = 1:(base_metacell_index - 1)
        @views total_UMIs_per_gene_per_other_metacells = total_UMIs_per_gene_per_metacell[:, index_per_other_metacell]
        @views low_scaled_log_covered_fraction_per_gene_per_other_metacells =
            confidence_intervals.low_scaled_log_covered_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]
        @views high_scaled_log_covered_fraction_per_gene_per_other_metacells =
            confidence_intervals.high_scaled_log_covered_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]

        significant_fold_per_gene_per_other_metacell =
            confident_gene_distance.(
                min_significant_gene_UMIs,
                base_metacell_total_UMIs_per_gene,
                base_metacell_low_scaled_log_covered_fraction_per_gene,
                base_metacell_high_scaled_log_covered_fraction_per_gene,
                total_UMIs_per_gene_per_other_metacells,
                low_scaled_log_covered_fraction_per_gene_per_other_metacells,
                high_scaled_log_covered_fraction_per_gene_per_other_metacells,
            )
        @assert_matrix(significant_fold_per_gene_per_other_metacell, n_genes, n_other_metacells, Columns)

        distances_between_base_and_other_metacells =
            vec(maximum(significant_fold_per_gene_per_other_metacell; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, n_other_metacells)

        distances_between_metacells[index_per_other_metacell, base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, index_per_other_metacell] .=
            distances_between_base_and_other_metacells
    end

    return distances_between_metacells
end

@inline function confident_gene_distance(  # UNTESTED
    min_significant_gene_UMIs::Integer,
    base_metacell_total_UMIs_of_gene::Integer,
    base_metacell_low_scaled_log_covered_fraction_of_gene::AbstractFloat,
    base_metacell_high_scaled_log_covered_fraction_of_gene::AbstractFloat,
    other_metacell_total_UMIs_of_gene::Integer,
    other_metacell_low_scaled_log_covered_fraction_of_gene::AbstractFloat,
    other_metacell_high_scaled_log_covered_fraction_of_gene::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = base_metacell_total_UMIs_of_gene + other_metacell_total_UMIs_of_gene
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low =
        base_metacell_high_scaled_log_covered_fraction_of_gene < other_metacell_high_scaled_log_covered_fraction_of_gene

    highest_low_scaled_log_covered_fraction_of_gene =
        is_base_low * base_metacell_high_scaled_log_covered_fraction_of_gene +
        !is_base_low * other_metacell_high_scaled_log_covered_fraction_of_gene

    lowest_high_scaled_log_covered_fraction_of_gene =
        is_base_low * other_metacell_low_scaled_log_covered_fraction_of_gene +
        !is_base_low * base_metacell_low_scaled_log_covered_fraction_of_gene

    confident_gap = lowest_high_scaled_log_covered_fraction_of_gene - highest_low_scaled_log_covered_fraction_of_gene

    return (is_significant * max(confident_gap, 0.0))
end

function compute_distances_between_blocks(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    block_distance_aggregation::Function,
)::AbstractMatrix{<:AbstractFloat}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacell_indices_per_block)
    aggregated_distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    aggregated_distances_between_blocks[1, 1] = -1.0
    @threads :greedy for base_block_index in reverse(2:n_blocks)
        aggregated_distances_between_blocks[base_block_index, base_block_index] = -1.0
        index_per_base_metacells = metacell_indices_per_block[base_block_index]

        distance_per_metacell_per_base_metacell = distances_between_metacells[:, index_per_base_metacells]

        aggregated_distance_from_base_per_metacell =
            vec(block_distance_aggregation(distance_per_metacell_per_base_metacell; dims = 2))
        @assert length(aggregated_distance_from_base_per_metacell) == n_metacells

        for other_block_index in 1:(base_block_index - 1)
            index_per_other_metacell = metacell_indices_per_block[other_block_index]

            @views aggregated_distance_from_base_per_other_metacell =
                aggregated_distance_from_base_per_metacell[index_per_other_metacell]

            aggregated_distance_between_base_and_other_block =
                block_distance_aggregation(aggregated_distance_from_base_per_other_metacell)

            aggregated_distances_between_blocks[base_block_index, other_block_index] =
                aggregated_distance_between_base_and_other_block
            aggregated_distances_between_blocks[other_block_index, base_block_index] =
                aggregated_distance_between_base_and_other_block
        end
    end

    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)
    @threads :greedy for base_block_index in 1:n_blocks
        @views aggregated_distance_between_others_and_base_block =
            aggregated_distances_between_blocks[:, base_block_index]
        rank_of_others_for_base_block = invperm(sortperm(aggregated_distance_between_others_and_base_block))
        distances_between_blocks[:, base_block_index] .= rank_of_others_for_base_block
    end

    distances_between_blocks .*= transpose(distances_between_blocks)
    distances_between_blocks .-= 1
    distances_between_blocks ./= maximum(distances_between_blocks)
    @assert minimum(distances_between_blocks) == 0

    return distances_between_blocks
end

function compute_neighborhoods!(
    daf::DafWriter;
    min_blocks_in_neighborhood::Integer,
    min_metacells_in_neighborhood::Integer,
    min_covered_UMIs_in_neighborhood::Integer,
    blocks::BlocksData,
    overwrite::Bool,
)::AbstractMatrix{<:Bool}
    name_per_block = axis_vector(daf, "block")
    @assert_vector(blocks.metacell_indices_per_block, blocks.n_blocks)
    @assert_matrix(blocks.distances_between_blocks, blocks.n_blocks, blocks.n_blocks, Columns)
    is_in_neighborhood_per_other_block_per_base_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)

    min_blocks_in_neighborhood = min(min_blocks_in_neighborhood, blocks.n_blocks)

    @threads :greedy for base_block_index in 1:(blocks.n_blocks)
        @views distance_from_base_per_block = blocks.distances_between_blocks[:, base_block_index]

        ordered_block_indices = sortperm(distance_from_base_per_block)
        @assert ordered_block_indices[1] == base_block_index

        neighborhood_n_blocks = 0
        neighborhood_n_metacells = 0
        neighborhood_covered_UMIs = 0

        while neighborhood_n_blocks < min_blocks_in_neighborhood ||
                  neighborhood_n_metacells < min_metacells_in_neighborhood ||
                  neighborhood_covered_UMIs < min_covered_UMIs_in_neighborhood
            neighborhood_n_blocks += 1
            next_block_index = ordered_block_indices[neighborhood_n_blocks]
            is_in_neighborhood_per_other_block_per_base_block[base_block_index, next_block_index] = true
            is_in_neighborhood_per_other_block_per_base_block[next_block_index, base_block_index] = true
            neighborhood_n_metacells += length(blocks.metacell_indices_per_block[next_block_index])
            neighborhood_covered_UMIs +=
                sum(blocks.covered_UMIs_per_metacell[blocks.metacell_indices_per_block[next_block_index]])
        end

        @debug (
            "Neighborhood: $(name_per_block[base_block_index])" *
            " Initial Blocks: $(neighborhood_n_blocks)" *
            " Metacells: $(neighborhood_n_metacells)" *
            " Covered M-UMIs: $(neighborhood_covered_UMIs / 1e6)"
        )
    end

    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_neighborhood",
        SparseMatrixCSC(is_in_neighborhood_per_other_block_per_base_block);
        overwrite,
    )

    total_blocks_in_neighborhoods = 0
    total_metacells_in_neighborhoods = 0
    covered_UMIs_in_neighborhoods = 0

    for base_block_index in 1:(blocks.n_blocks)
        @views is_in_neighborhood_of_base_per_block =
            is_in_neighborhood_per_other_block_per_base_block[:, base_block_index]
        block_indices_of_neighborhood = findall(is_in_neighborhood_of_base_per_block)
        neighborhood_n_blocks = length(block_indices_of_neighborhood)
        neighborhood_n_metacells = sum(length.(blocks.metacell_indices_per_block[block_indices_of_neighborhood]))
        neighborhood_covered_UMIs = sum([
            sum(blocks.covered_UMIs_per_metacell[blocks.metacell_indices_per_block[block_index]]) for
            block_index in block_indices_of_neighborhood
        ])
        @debug "Neighborhood: $(name_per_block[base_block_index]) Total Blocks: $(neighborhood_n_blocks) Metacells: $(neighborhood_n_metacells) Covered M-UMIs: $(neighborhood_covered_UMIs / 1e6)"
        @debug (
            "Neighborhood: $(name_per_block[base_block_index])" *
            " Final Blocks: $(neighborhood_n_blocks)" *
            " Metacells: $(neighborhood_n_metacells)" *
            " Covered M-UMIs: $(neighborhood_covered_UMIs / 1e6)"
        )

        total_blocks_in_neighborhoods += neighborhood_n_blocks
        total_metacells_in_neighborhoods += neighborhood_n_metacells
        covered_UMIs_in_neighborhoods += neighborhood_covered_UMIs
    end

    @info (
        "Neighborhoods mean Blocks: $(total_blocks_in_neighborhoods[] / blocks.n_blocks)" *
        " Metacells: $(total_metacells_in_neighborhoods[] / blocks.n_blocks)" *
        " Covered M-UMIs: $(covered_UMIs_in_neighborhoods[] / (1e6 * blocks.n_blocks))"
    )

    return is_in_neighborhood_per_other_block_per_base_block
end

@kwdef struct ReusableMatrices
    skeleton_coefficient_per_gene_per_principal_component::Matrix{Float32}
    covered_coefficient_per_principal_component_per_gene::Matrix{Float32}
end

function reset_reusable_matrices(reusable_matrices::ReusableMatrices)::Nothing
    reusable_matrices.skeleton_coefficient_per_gene_per_principal_component .= 0
    reusable_matrices.covered_coefficient_per_principal_component_per_gene .= 0
    return nothing
end

function compute_environments!(
    daf::DafWriter;
    max_principal_components::Integer,
    cross_validation_parts::Integer,
    is_in_neighborhood_per_other_block_per_base_block::AbstractMatrix{<:Bool},
    blocks::BlocksData,
    rng::AbstractRNG,
    overwrite::Bool,
    leaky_pca_cross_validation::Bool,
)::Nothing
    n_genes = axis_length(daf, "gene")

    index_per_skeleton = daf["/ gene & is_skeleton : index"].array
    index_per_covered = daf["/ gene & is_covered : index"].array

    scaled_log_covered_fraction_per_skeleton_per_metacell =
        daf["/ gene & is_skeleton / metacell : scaled_log_covered_fraction"].array
    scaled_log_covered_fraction_per_covered_per_metacell =
        daf["/ gene & is_covered / metacell : scaled_log_covered_fraction"].array

    is_in_environment_per_other_block_per_base_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    linear_RMSE_per_block = Vector{Float32}(undef, blocks.n_blocks)
    linear_XRMSE_per_block = Vector{Float32}(undef, blocks.n_blocks)
    n_principal_components_per_block = Vector{UInt32}(undef, blocks.n_blocks)
    is_used_per_principal_component_per_block = zeros(Bool, max_principal_components, blocks.n_blocks)

    environment_n_principal_components_per_block = zeros(Int32, blocks.n_blocks)
    mean_scaled_log_covered_fraction_per_gene_per_block = zeros(Float32, n_genes, blocks.n_blocks)

    name_per_block = axis_vector(daf, "block")

    total_blocks_in_environments = Atomic{Int}(0)
    total_metacells_in_environments = Atomic{Int}(0)
    covered_UMIs_in_environments = Atomic{UInt64}(0)

    reusable_storage = ReusableStorage(reset_reusable_matrices) do
        return ReusableMatrices(;
            skeleton_coefficient_per_gene_per_principal_component = zeros(Float32, n_genes, max_principal_components),
            covered_coefficient_per_principal_component_per_gene = zeros(Float32, max_principal_components, n_genes),
        )
    end

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(blocks.n_blocks; rng) do base_block_index, rng
        @views is_in_neighborhood_of_base_per_block =
            is_in_neighborhood_per_other_block_per_base_block[:, base_block_index]
        @views distance_from_base_per_block = blocks.distances_between_blocks[:, base_block_index]

        ordered_block_indices = sortperm(distance_from_base_per_block)
        @assert ordered_block_indices[1] == base_block_index

        block_indices_of_neighborhood = findall(is_in_neighborhood_of_base_per_block)
        neighborhood_n_blocks = length(block_indices_of_neighborhood)

        additional_ordered_block_indices = filter(ordered_block_indices) do block_index
            return !(block_index in block_indices_of_neighborhood)
        end
        @assert length(additional_ordered_block_indices) + neighborhood_n_blocks == blocks.n_blocks

        metacell_indices_of_neighborhood = vcat(blocks.metacell_indices_per_block[block_indices_of_neighborhood]...)

        neighborhood_cross_validation_indices = pick_cross_validation_indices(;
            full_indices = metacell_indices_of_neighborhood,
            cross_validation_parts,
            rng,
        )

        environment_n_blocks, _, environment_model = minimize_cost(;
            minimal_value = neighborhood_n_blocks,
            maximal_value = blocks.n_blocks,
        ) do environment_m_blocks
            if environment_m_blocks == neighborhood_n_blocks
                additional_metacell_indices_of_environment = Int[]
            else
                @views additional_block_indices_of_environment =
                    additional_ordered_block_indices[1:(environment_m_blocks - neighborhood_n_blocks)]
                additional_metacell_indices_of_environment =
                    vcat(blocks.metacell_indices_per_block[additional_block_indices_of_environment]...)
            end

            @assert isempty(
                intersect(Set(additional_metacell_indices_of_environment), Set(metacell_indices_of_neighborhood)),
            ) # TODOX

            environment_model = compute_environment_model(
                max_principal_components,
                metacell_indices_of_neighborhood,
                neighborhood_cross_validation_indices,
                additional_metacell_indices_of_environment,
                scaled_log_covered_fraction_per_skeleton_per_metacell,
                scaled_log_covered_fraction_per_covered_per_metacell,
                leaky_pca_cross_validation,
            )

            #@debug "TODOX BLKs: $(environment_m_blocks) PCS: $(environment_model.n_principal_components) RMSE: $(environment_model.RMSE) XRMSE: $(environment_model.XRMSE)"
            return (environment_model.XRMSE, environment_model)
        end

        @views block_indices_of_environment = ordered_block_indices[1:environment_n_blocks]
        is_in_environment_per_other_block_per_base_block[block_indices_of_environment, base_block_index] .= true
        environment_n_principal_components_per_block[base_block_index] = environment_model.n_principal_components

        linear_RMSE_per_block[base_block_index] = environment_model.RMSE
        linear_XRMSE_per_block[base_block_index] = environment_model.XRMSE
        n_principal_components_per_block[base_block_index] = environment_model.n_principal_components
        is_used_per_principal_component_per_block[1:(environment_model.n_principal_components), base_block_index] .=
            true

        with_reusable(reusable_storage) do reusable_matrices
            reusable_matrices.skeleton_coefficient_per_gene_per_principal_component[
                index_per_skeleton,
                1:(environment_model.n_principal_components),
            ] .= environment_model.coefficient_per_skeleton_per_principal_component

            reusable_matrices.covered_coefficient_per_principal_component_per_gene[
                1:(environment_model.n_principal_components),
                index_per_covered,
            ] .= SparseMatrixCSC(environment_model.coefficient_per_principal_component_per_covered)

            mean_scaled_log_covered_fraction_per_gene_per_block[index_per_covered, base_block_index] =
                environment_model.mean_scaled_log_covered_fraction_per_covered

            environment_n_metacells = sum(length.(blocks.metacell_indices_per_block[block_indices_of_environment]))
            environment_covered_UMIs = sum([
                sum(blocks.covered_UMIs_per_metacell[blocks.metacell_indices_per_block[block_index]]) for
                block_index in block_indices_of_environment
            ])

            sparse_skeleton_coefficient_per_gene_per_principal_component =
                SparseMatrixCSC(reusable_matrices.skeleton_coefficient_per_gene_per_principal_component)
            sparse_covered_coefficient_per_principal_component_per_gene =
                SparseMatrixCSC(reusable_matrices.covered_coefficient_per_principal_component_per_gene)

            set_matrix!(
                daf,
                "gene",
                "principal_component",
                "$(name_per_block[base_block_index])_skeleton_coefficient",
                sparse_skeleton_coefficient_per_gene_per_principal_component;
                overwrite,
            )

            set_matrix!(
                daf,
                "principal_component",
                "gene",
                "$(name_per_block[base_block_index])_covered_coefficient",
                sparse_covered_coefficient_per_principal_component_per_gene;
                overwrite,
            )

            atomic_add!(total_blocks_in_environments, environment_n_blocks)
            atomic_add!(total_metacells_in_environments, environment_n_metacells)
            atomic_add!(covered_UMIs_in_environments, UInt64(environment_covered_UMIs))

            counter = atomic_add!(progress_counter, 1)
            @debug (
                "- Environment: $(name_per_block[base_block_index]) ($(percent(counter + 1, blocks.n_blocks)))" *
                " Blocks: $(environment_n_blocks)" *
                " Metacells: $(environment_n_metacells)" *
                " Covered M-UMIs: $(environment_covered_UMIs / 1e6)" *
                " Principal components: $(environment_model.n_principal_components)" *
                " RMSE: $(environment_model.RMSE)" *
                " XRMSE: $(environment_model.XRMSE)"
            )

            return nothing
        end
    end

    @info (
        "Environments mean Blocks: $(total_blocks_in_environments[] / blocks.n_blocks)" *
        " Metacells: $(total_metacells_in_environments[] / blocks.n_blocks)" *
        " Covered M-UMIs: $(covered_UMIs_in_environments[] / (1e6 * blocks.n_blocks))" *
        " Principal components: $(mean(n_principal_components_per_block))" *
        " RMSE: $(mean(linear_RMSE_per_block))" *
        " XRMSE: $(mean(linear_XRMSE_per_block))"
    )

    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_environment",
        SparseMatrixCSC(is_in_environment_per_other_block_per_base_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "mean_scaled_log_covered_fraction",
        SparseMatrixCSC(mean_scaled_log_covered_fraction_per_gene_per_block);
        overwrite,
    )
    set_vector!(daf, "block", "n_principal_components", n_principal_components_per_block; overwrite)
    set_vector!(daf, "block", "linear_RMSE", linear_RMSE_per_block; overwrite)
    set_vector!(daf, "block", "linear_XRMSE", linear_XRMSE_per_block; overwrite)
    set_matrix!(
        daf,
        "principal_component",
        "block",
        "is_used",
        SparseMatrixCSC(is_used_per_principal_component_per_block);
        overwrite,
    )

    return nothing
end

function compute_environment_model(
    max_principal_components::Integer,
    metacell_indices_of_neighborhood::AbstractVector{<:Integer},
    neighborhood_cross_validation_indices::CrossValidationIndices,
    additional_metacell_indices_of_environment::AbstractVector{<:Integer},
    scaled_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    scaled_log_covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    leaky_pca_cross_validation::Bool = false,
)::LocalModel
    metacell_indices_of_environment = vcat(metacell_indices_of_neighborhood, additional_metacell_indices_of_environment)

    full_pca_model = compute_local_model(;
        scaled_log_covered_fraction_per_skeleton_per_metacell,
        scaled_log_covered_fraction_per_covered_per_metacell,
        metacell_indices = metacell_indices_of_environment,
        max_principal_components,
        base_model = nothing,
        only_pca = true,
    )

    trained_model_per_part = [
        compute_local_model(;
            scaled_log_covered_fraction_per_skeleton_per_metacell,
            scaled_log_covered_fraction_per_covered_per_metacell,
            metacell_indices = vcat(train_indices, additional_metacell_indices_of_environment),
            max_principal_components,
            base_model = leaky_pca_cross_validation ? full_pca_model : nothing,
            only_pca = false,
        ) for train_indices in neighborhood_cross_validation_indices.train_indices_per_part
    ]

    max_principal_components =
        minimum([trained_model.n_principal_components for trained_model in trained_model_per_part])

    environment_n_principal_components, XRMSE, _ = minimize_cost(;
        minimal_value = 1,
        maximal_value = max_principal_components,
    ) do environment_m_principal_components
        RMSE_per_part = [
            RMSE_of_model(;
                local_model = limited_model(trained_model, environment_m_principal_components),
                scaled_log_covered_fraction_per_skeleton_per_metacell,
                scaled_log_covered_fraction_per_covered_per_metacell,
                metacell_indices = test_indices,
            ) for (trained_model, test_indices) in
            zip(trained_model_per_part, neighborhood_cross_validation_indices.test_indices_per_part)
        ]
        #@debug "TODOX PCs: $(environment_m_principal_components) XRMSE: $(mean(RMSE_per_part))"
        return (mean(RMSE_per_part), nothing)
    end

    full_model = compute_local_model(;
        scaled_log_covered_fraction_per_skeleton_per_metacell,
        scaled_log_covered_fraction_per_covered_per_metacell,
        metacell_indices = metacell_indices_of_environment,
        max_principal_components = environment_n_principal_components,
        base_model = limited_model(full_pca_model, environment_n_principal_components),
        only_pca = false,
    )

    full_model.XRMSE = XRMSE
    full_model.RMSE = RMSE_of_model(;
        local_model = full_model,
        scaled_log_covered_fraction_per_skeleton_per_metacell,
        scaled_log_covered_fraction_per_covered_per_metacell,
        metacell_indices = metacell_indices_of_environment,
    )

    return full_model
end

function limited_model(local_model::LocalModel, n_principal_components::Integer)::LocalModel
    @assert 1 <= n_principal_components <= local_model.n_principal_components

    if n_principal_components == local_model.n_principal_components
        return local_model
    end

    mean_scaled_log_covered_fraction_per_skeleton = local_model.mean_scaled_log_covered_fraction_per_skeleton

    coefficient_per_skeleton_per_principal_component =
        local_model.coefficient_per_skeleton_per_principal_component[:, 1:n_principal_components]

    mean_scaled_log_covered_fraction_per_covered = local_model.mean_scaled_log_covered_fraction_per_covered

    if local_model.coefficient_per_principal_component_per_covered === nothing
        coefficient_per_principal_component_per_covered = nothing
    else
        coefficient_per_principal_component_per_covered =
            local_model.coefficient_per_principal_component_per_covered[1:n_principal_components, :]
    end

    return LocalModel(;
        n_principal_components,
        mean_scaled_log_covered_fraction_per_skeleton,
        coefficient_per_skeleton_per_principal_component,
        mean_scaled_log_covered_fraction_per_covered,
        coefficient_per_principal_component_per_covered,
    )
end

function compute_local_model(;
    scaled_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    scaled_log_covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    metacell_indices::AbstractVector{<:Integer},
    max_principal_components::Integer,
    base_model::Maybe{LocalModel},
    only_pca::Bool,
)::LocalModel
    n_metacells = size(scaled_log_covered_fraction_per_skeleton_per_metacell, 2)
    n_skeleton = size(scaled_log_covered_fraction_per_skeleton_per_metacell, 1)
    n_covered = size(scaled_log_covered_fraction_per_covered_per_metacell, 1)
    n_model_metacells = length(metacell_indices)

    @assert_matrix(scaled_log_covered_fraction_per_skeleton_per_metacell, n_skeleton, n_metacells, Columns)
    @assert_matrix(scaled_log_covered_fraction_per_covered_per_metacell, n_covered, n_metacells, Columns)

    @views scaled_log_covered_fraction_per_skeleton_per_model_metacell =
        scaled_log_covered_fraction_per_skeleton_per_metacell[:, metacell_indices]
    @views scaled_log_covered_fraction_per_covered_per_model_metacell =
        scaled_log_covered_fraction_per_covered_per_metacell[:, metacell_indices]

    if base_model !== nothing
        n_principal_components = base_model.n_principal_components
        mean_scaled_log_covered_fraction_per_skeleton = base_model.mean_scaled_log_covered_fraction_per_skeleton
        offset_scaled_log_covered_fraction_per_skeleton_per_model_metacell, _ = centralize(
            scaled_log_covered_fraction_per_skeleton_per_model_metacell,
            mean_scaled_log_covered_fraction_per_skeleton,
        )
        coefficient_per_skeleton_per_principal_component = base_model.coefficient_per_skeleton_per_principal_component

    else
        offset_scaled_log_covered_fraction_per_skeleton_per_model_metacell,
        mean_scaled_log_covered_fraction_per_skeleton =
            centralize(scaled_log_covered_fraction_per_skeleton_per_model_metacell)

        pca = fit(
            PCA,
            offset_scaled_log_covered_fraction_per_skeleton_per_model_metacell;
            mean = 0,
            maxoutdim = max_principal_components,
        )
        coefficient_per_skeleton_per_principal_component = loadings(pca)
        n_principal_components = outdim(pca)
        @assert_matrix(coefficient_per_skeleton_per_principal_component, n_skeleton, n_principal_components, Columns)
    end

    if only_pca
        mean_scaled_log_covered_fraction_per_covered = nothing
        coefficient_per_principal_component_per_covered = nothing

    else
        principal_component_per_model_metacell =
            transpose(coefficient_per_skeleton_per_principal_component) *
            offset_scaled_log_covered_fraction_per_skeleton_per_model_metacell
        @assert_matrix(principal_component_per_model_metacell, n_principal_components, n_model_metacells, Columns)

        offset_principal_component_per_model_metacell, mean_per_principal_component =
            centralize(principal_component_per_model_metacell)
        if base_model === nothing
            max_mean_principal_component = maximum(abs.(mean_per_principal_component))
            max_offset_principal_component = maximum(abs.(offset_principal_component_per_model_metacell))
            max_mean_principal_component_ratio = max_mean_principal_component / max_offset_principal_component
            @assert max_mean_principal_component_ratio < 1e-4
        end

        offset_scaled_log_covered_fraction_per_covered_per_model_metacell,
        mean_scaled_log_covered_fraction_per_covered =
            centralize(scaled_log_covered_fraction_per_covered_per_model_metacell)

        coefficient_per_principal_component_per_covered =
            transpose(offset_principal_component_per_model_metacell) \
            transpose(offset_scaled_log_covered_fraction_per_covered_per_model_metacell)
        @assert_matrix(coefficient_per_principal_component_per_covered, n_principal_components, n_covered, Columns)
    end

    return LocalModel(;
        n_principal_components,
        mean_scaled_log_covered_fraction_per_skeleton,
        coefficient_per_skeleton_per_principal_component,
        mean_scaled_log_covered_fraction_per_covered,
        coefficient_per_principal_component_per_covered,
    )
end

function RMSE_of_model(;
    local_model::LocalModel,
    scaled_log_covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    scaled_log_covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    metacell_indices::AbstractVector{<:Integer},
)::AbstractFloat
    n_metacells = size(scaled_log_covered_fraction_per_skeleton_per_metacell, 2)
    n_skeletons = size(scaled_log_covered_fraction_per_skeleton_per_metacell, 1)
    n_covered = size(scaled_log_covered_fraction_per_covered_per_metacell, 1)

    @assert_matrix(scaled_log_covered_fraction_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert_matrix(scaled_log_covered_fraction_per_covered_per_metacell, n_covered, n_metacells, Columns)

    principal_component_per_model_metacell = compute_profiles_principal_components(;
        local_model,
        scaled_log_covered_fraction_per_skeleton_per_profile = scaled_log_covered_fraction_per_skeleton_per_metacell,
        profile_indices = metacell_indices,
    )

    @views scaled_log_covered_fraction_per_covered_per_model_metacell =
        scaled_log_covered_fraction_per_covered_per_metacell[:, metacell_indices]

    squared_error_per_covered_per_model_metacell =
        (
            (
                (
                    transpose(local_model.coefficient_per_principal_component_per_covered) *
                    principal_component_per_model_metacell
                ) .+ local_model.mean_scaled_log_covered_fraction_per_covered
            ) .- scaled_log_covered_fraction_per_covered_per_model_metacell
        ) .^ 2

    return sqrt(mean(squared_error_per_covered_per_model_metacell))
end

function centralize(
    value_per_gene_per_metacell::AbstractMatrix{<:AbstractFloat},
    mean_per_gene::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractVector{<:AbstractFloat}}
    n_genes, n_metacells = size(value_per_gene_per_metacell)
    @assert_matrix(value_per_gene_per_metacell, n_genes, n_metacells)
    if mean_per_gene === nothing
        mean_per_gene = vec(mean(value_per_gene_per_metacell; dims = 2))  # NOLINT
    end
    @assert_vector(mean_per_gene, n_genes)

    return (value_per_gene_per_metacell .- mean_per_gene, mean_per_gene)
end

function pick_cross_validation_indices(;
    full_indices::AbstractVector{<:Integer},
    cross_validation_parts::Integer,
    rng::AbstractRNG,
)::CrossValidationIndices
    n_full = length(full_indices)
    @assert n_full >= cross_validation_parts

    train_indices_per_part = Vector{Int}[]
    test_indices_per_part = Vector{Int}[]

    shuffled_full_indices = shuffle(rng, full_indices)
    parts_size = n_full / cross_validation_parts

    for part_index in 1:cross_validation_parts
        first_test_position = Int(round((part_index - 1) * parts_size)) + 1
        last_test_position = Int(round(part_index * parts_size))
        test_positions = first_test_position:last_test_position

        train_indices = vcat(
            shuffled_full_indices[1:(first_test_position - 1)],
            shuffled_full_indices[(last_test_position + 1):end],
        )
        test_indices = shuffled_full_indices[test_positions]
        @assert length(test_indices) + length(train_indices) == n_full

        push!(train_indices_per_part, train_indices)
        push!(test_indices_per_part, test_indices)
    end

    return CrossValidationIndices(; train_indices_per_part, test_indices_per_part)
end

function minimize_cost(  # UNTESTED
    compute_cost_result::Function;
    minimal_value::Integer,
    maximal_value::Integer,
)::Tuple{Integer, AbstractFloat, Any}
    min_significant_cost = 1e-3

    @assert min_significant_cost > 0
    @assert 0 < minimal_value <= maximal_value

    if minimal_value == maximal_value
        return (minimal_value, compute_cost_result(minimal_value)...)
    end

    sample_cost, sample_result = compute_cost_result(minimal_value)
    sampled_values = [minimal_value]
    minimal_cost = sample_cost
    sampled_costs = [sample_cost]
    sampled_results = Any[sample_result]

    offset = 1
    while sampled_values[end] != maximal_value
        sample_value = min(minimal_value + offset, maximal_value)
        offset *= 2
        sample_cost, sample_result = compute_cost_result(sample_value)
        push!(sampled_values, sample_value)
        push!(sampled_costs, sample_cost)
        minimal_cost = min(minimal_cost, sample_cost)
        push!(sampled_results, sample_result)

        if sample_value == maximal_value ||
           sample_cost == 0.0 ||
           (length(sampled_values) > 4 && sampled_costs[end] - minimal_cost > min_significant_cost * 2)
            break
        end
    end

    while true
        best_index = nothing
        best_cost = nothing
        for (index, cost) in enumerate(sampled_costs)
            if cost <= minimal_cost + min_significant_cost
                best_index = index
                best_cost = cost
                break
            end
        end
        @assert best_index !== nothing
        @assert best_cost !== nothing

        best_value = sampled_values[best_index]
        sample_value = nothing
        if best_index == 1
            next_value = sampled_values[2]
            if next_value > best_value + 1
                sample_value = div(best_value + next_value, 2)
            end
        elseif best_index == length(sampled_costs)
            prev_value = sampled_values[end - 1]
            if prev_value < best_value - 1
                sample_value = div(best_value + prev_value, 2)
            end
        else
            next_value = sampled_values[best_index + 1]
            prev_value = sampled_values[best_index - 1]
            next_gap = next_value - best_value
            prev_gap = best_value - prev_value
            if prev_gap > 1 && prev_gap >= next_gap
                sample_value = div(best_value + prev_value, 2)
            elseif next_gap > 1
                sample_value = div(best_value + next_value, 2)
            end
        end

        if sample_value === nothing
            return (best_value, sampled_costs[best_index], sampled_results[best_index])
        end

        sample_cost, sample_result = compute_cost_result(sample_value)
        sample_index = searchsortedfirst(sampled_values, sample_value)
        insert!(sampled_values, sample_index, sample_value)
        insert!(sampled_costs, sample_index, sample_cost)
        minimal_cost = min(minimal_cost, sample_cost)
        insert!(sampled_results, sample_index, sample_result)
    end
end

"""
TODOX
"""
@logged @computation Contract(
    axes = [
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
        old_metacell_axis(RequiredInput),
        new_metacell_axis(GuaranteedOutput),
        block_axis(RequiredInput),
        principal_component_axis(RequiredInput),
    ],
    data = [
        gene_divergence_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_old_metacell_vector(RequiredInput),
        cell_new_metacell_vector(GuaranteedOutput),
        old_metacell_block_vector(RequiredInput),
        new_metacell_block_vector(GuaranteedOutput),
        block_n_principal_components_vector(RequiredInput),
        block_principal_component_is_used_matrix(RequiredInput),
        block_gene_mean_scaled_log_covered_fraction_matrix(RequiredInput),
        block_principal_component_gene_skeleton_coefficient_tensor(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
    ],
) function sharpen_metacells!(
    daf::DafWriter;
    gene_UMIs_regularization::AbstractFloat = 1 / 16,
    min_migrated_cells_fraction::AbstractFloat = 0.01,
    min_metacell_size::Integer = 12,
    kmeans_rounds::Integer = 10,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_UMIs_regularization > 0
    @assert min_migrated_cells_fraction < 1
    @assert min_metacell_size >= 0
    @assert kmeans_rounds > 0

    n_cells = axis_length(daf, "cell")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")
    n_principal_components_per_block = get_vector(daf, "block", "n_principal_components").array

    mean_metacell_cells_per_block = [
        daf["/ cell & old_metacell ?? => block => is_in_neighborhood ;= $(base_block_name) : index %> Count"] /
        daf["/ old_metacell & block => is_in_neighborhood ;= $(base_block_name) : index %> Count"] for
        base_block_name in name_per_block
    ]

    block_index_per_cell = daf["/ cell : old_metacell ?? 0 => block => index"].array
    block_indices_per_neighborhood = [
        parent(daf["/ block & is_in_neighborhood ;= $(name_per_block[base_block_index]) : index"].array) for
        base_block_index in 1:n_blocks
    ]

    scaled_log_covered_fraction_per_skeleton_per_cell =
        compute_scaled_log_covered_fraction_per_skeleton_per_cell(daf; gene_UMIs_regularization)

    block_index_per_cell_per_round = AbstractVector{<:Integer}[block_index_per_cell]
    n_migrated = n_cells
    while n_migrated >= n_cells * min_migrated_cells_fraction
        preferred_block_index_per_cell_per_block = compute_preferred_block_index_per_cell_per_block(
            daf;
            kmeans_rounds,
            scaled_log_covered_fraction_per_skeleton_per_cell,
            name_per_block,
            n_principal_components_per_block,
            block_index_per_cell,
            mean_metacell_cells_per_block,
            block_indices_per_neighborhood,
            rng,
        )

        block_index_per_cell, n_migrated = compute_preferred_block_index_of_cells(;
            block_index_per_cell,
            preferred_block_index_per_cell_per_block,
            block_index_per_cell_per_round,
        )
        push!(block_index_per_cell_per_round, block_index_per_cell)
    end

    new_metacell_index_per_cell, block_index_per_new_metacell = compute_new_metacells(
        daf;
        min_metacell_size,
        kmeans_rounds,
        scaled_log_covered_fraction_per_skeleton_per_cell,
        name_per_block,
        n_principal_components_per_block,
        block_index_per_cell,
        mean_metacell_cells_per_block,
        rng,
    )

    cell_indices_per_new_metacell = collect_group_members(new_metacell_index_per_cell)
    name_per_new_metacell = group_names(daf, "cell", cell_indices_per_new_metacell; prefix = "M")
    add_axis!(daf, "new_metacell", name_per_new_metacell)

    new_metacell_name_per_cell = [
        new_metacell_index == 0 ? "" : name_per_new_metacell[new_metacell_index] for
        new_metacell_index in new_metacell_index_per_cell
    ]
    set_vector!(daf, "cell", "new_metacell", new_metacell_name_per_cell; overwrite)

    block_name_per_new_metacell = name_per_block[block_index_per_new_metacell]
    set_vector!(daf, "new_metacell", "block", block_name_per_new_metacell; overwrite)

    return nothing
end

function compute_scaled_log_covered_fraction_per_skeleton_per_cell(
    daf::DafReader;
    gene_UMIs_regularization::AbstractFloat,
)::AbstractMatrix{<:AbstractFloat}
    n_cells = axis_length(daf, "cell")
    n_skeletons = sum(daf["/ gene : is_skeleton"].array)

    covered_UMIs_per_cell = densify(daf["/ gene & is_covered / cell : UMIs %> Sum"].array)
    @assert_vector(covered_UMIs_per_cell, n_cells)

    UMIs_per_skeleton_per_cell = densify(daf["/ gene & is_skeleton / cell : UMIs"].array)
    @assert_matrix(UMIs_per_skeleton_per_cell, n_skeletons, n_cells, Columns)

    divergence_per_skeleton = densify(daf["/ gene & is_skeleton : divergence"].array)
    @assert_vector(divergence_per_skeleton, n_skeletons)

    scaled_log_covered_fraction_per_skeleton_per_cell =
        log2.(UMIs_per_skeleton_per_cell ./ transpose(covered_UMIs_per_cell) .+ gene_UMIs_regularization) .*
        (1.0 .- divergence_per_skeleton)

    return scaled_log_covered_fraction_per_skeleton_per_cell
end

function compute_preferred_block_index_per_cell_per_block(
    daf::DafReader;
    kmeans_rounds::Integer,
    scaled_log_covered_fraction_per_skeleton_per_cell::AbstractMatrix{<:AbstractFloat},
    name_per_block::AbstractVector{<:AbstractString},
    n_principal_components_per_block::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    block_indices_per_neighborhood::AbstractVector{<:AbstractVector{<:Integer}},
    rng::AbstractRNG,
)::Vector{Maybe{SparseVector{<:Integer}}}
    n_cells = length(block_index_per_cell)
    n_blocks = length(name_per_block)

    cell_indices_per_block = [findall(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    preferred_block_index_per_cell_per_block = Vector{Maybe{SparseVector{<:Integer}}}(undef, n_blocks)
    preferred_block_index_per_cell_per_block .= nothing

    #progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(n_blocks; rng) do base_block_index, rng
        #base_block_name = name_per_block[base_block_index]

        neighborhood_cell_indices = vcat(cell_indices_per_block[block_indices_per_neighborhood[base_block_index]]...)
        n_neighborhood_cells = length(neighborhood_cell_indices)
        if n_neighborhood_cells == 0
            return nothing
        end

        local_model = load_local_model(daf; block_index = base_block_index, only_pca = true)
        @assert local_model.n_principal_components == n_principal_components_per_block[base_block_index]

        principal_component_per_neighborhood_cell = compute_profiles_principal_components(;
            local_model,
            scaled_log_covered_fraction_per_skeleton_per_profile = scaled_log_covered_fraction_per_skeleton_per_cell,
            profile_indices = neighborhood_cell_indices,
        )

        n_neighborhood_clusters =
            max(Int(round(n_neighborhood_cells / mean_metacell_cells_per_block[base_block_index])), 1)

        kmeans_result =
            kmeans_in_rounds(principal_component_per_neighborhood_cell, n_neighborhood_clusters; kmeans_rounds, rng)
        cluster_index_per_neighborhood_cell = assignments(kmeans_result)
        n_clusters = length(counts(kmeans_result))

        block_index_per_neighborhood_cell = vcat(
            [
                fill(block_index, length(cell_indices_per_block[block_index])) for
                block_index in block_indices_per_neighborhood[base_block_index]
            ]...,
        )
        @assert_vector(block_index_per_neighborhood_cell, n_neighborhood_cells)

        preferred_block_index_per_cluster = [
            mode(block_index_per_neighborhood_cell[cluster_index_per_neighborhood_cell .== cluster_index]) for
            cluster_index in 1:n_clusters
        ]
        preferred_block_index_per_neighborhood_cell =
            preferred_block_index_per_cluster[cluster_index_per_neighborhood_cell]

        reorder = sortperm(neighborhood_cell_indices)

        return preferred_block_index_per_cell_per_block[base_block_index] = SparseVector(
            n_cells,
            neighborhood_cell_indices[reorder],
            preferred_block_index_per_neighborhood_cell[reorder],
        )

        # counter = atomic_add!(progress_counter, 1)
        # n_stable_cells = sum(preferred_block_index_per_neighborhood_cell .== base_block_index)
        # @debug (
        #     "- Neighborhood: $(base_block_name) ($(percent(counter + 1, n_blocks)))" *
        #     " Cells: $(n_neighborhood_cells)" *
        #     " Stable: $(percent(n_stable_cells, n_neighborhood_cells))"
        # )
    end

    return preferred_block_index_per_cell_per_block
end

function compute_preferred_block_index_of_cells(;
    block_index_per_cell::AbstractVector{<:Integer},
    preferred_block_index_per_cell_per_block::Vector{Maybe{SparseVector{<:Integer}}},
    block_index_per_cell_per_round::AbstractVector{<:AbstractVector{<:Integer}},
)::Tuple{Vector{<:Integer}, Integer}
    n_cells = length(block_index_per_cell)

    n_quiescent = Atomic{Int}(0)
    n_stationary = Atomic{Int}(0)
    n_orbited = Atomic{Int}(0)
    n_migrated = Atomic{Int}(0)

    preferred_blocks_of_cells = zeros(UInt32, n_cells)
    @threads :greedy for cell_index in 1:n_cells
        original_block_index_of_cell = block_index_per_cell[cell_index]
        if original_block_index_of_cell == 0
            continue
        end

        preferred_block_index_per_block = preferred_block_index_per_cell_per_block[original_block_index_of_cell]
        if preferred_block_index_per_block === nothing
            preferred_blocks_of_cells[cell_index] = original_block_index_of_cell
            atomic_add!(n_quiescent, 1)
            continue
        end

        preferred_block_index_of_cell = preferred_block_index_per_block[cell_index]
        @assert preferred_block_index_of_cell > 0

        preferred_block_index_per_other_block = preferred_block_index_per_cell_per_block[preferred_block_index_of_cell]
        if preferred_block_index_per_other_block === nothing
            preferred_blocks_of_cells[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
            continue
        end

        back_preferred_block_index_of_cell = preferred_block_index_per_other_block[cell_index]
        @assert back_preferred_block_index_of_cell > 0

        if preferred_block_index_of_cell == original_block_index_of_cell ||
           back_preferred_block_index_of_cell != preferred_block_index_of_cell
            preferred_blocks_of_cells[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
        else
            preferred_blocks_of_cells[cell_index] = preferred_block_index_of_cell
            orbited = false
            for old_block_index_per_cell in block_index_per_cell_per_round
                if preferred_block_index_of_cell == old_block_index_per_cell[cell_index]
                    orbited = true
                    break
                end
            end
            if orbited
                atomic_add!(n_orbited, 1)
            else
                atomic_add!(n_migrated, 1)
            end
        end
    end

    @debug (
        "Cells: $(n_cells)" *
        " Stationary: $(n_stationary[]) ($(percent(n_stationary[], n_cells)))" *
        " Orbited: $(n_orbited[]) ($(percent(n_orbited[], n_cells)))" *
        " Migrated: $(n_migrated[]) ($(percent(n_migrated[], n_cells)))"
    )

    return (preferred_blocks_of_cells, n_migrated[])
end

@kwdef struct LocalClusters
    block_cell_indices::AbstractVector{<:Integer}
    cluster_index_per_block_cell::AbstractVector{<:Integer}
    is_too_small_per_cluster::Union{AbstractVector{Bool}, BitVector}
end

function compute_new_metacells(
    daf::DafReader;
    min_metacell_size::Integer,
    kmeans_rounds::Integer,
    scaled_log_covered_fraction_per_skeleton_per_cell::AbstractMatrix{<:AbstractFloat},
    name_per_block::AbstractVector{<:AbstractString},
    n_principal_components_per_block::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    rng::AbstractRNG,
)::Tuple{AbstractVector{<:Integer}, AbstractVector{<:Integer}}
    n_cells = size(scaled_log_covered_fraction_per_skeleton_per_cell, 2)

    local_clusters_per_block = compute_local_clusters(
        daf;
        min_metacell_size,
        kmeans_rounds,
        scaled_log_covered_fraction_per_skeleton_per_cell,
        name_per_block,
        n_principal_components_per_block,
        block_index_per_cell,
        mean_metacell_cells_per_block,
        rng,
    )

    n_old_metacells = axis_length(daf, "old_metacell")
    return combine_local_clusters(; local_clusters_per_block, n_old_metacells, n_cells)
end

function compute_local_clusters(
    daf::DafReader;
    min_metacell_size::Integer,
    kmeans_rounds::Integer,
    scaled_log_covered_fraction_per_skeleton_per_cell::AbstractMatrix{<:AbstractFloat},
    name_per_block::AbstractVector{<:AbstractString},
    n_principal_components_per_block::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    rng::AbstractRNG,
)::AbstractVector{Maybe{LocalClusters}}
    n_blocks = length(name_per_block)

    local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)

    #   progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(n_blocks; rng) do block_index, rng
        #block_name = name_per_block[block_index]

        block_cell_indices = findall(block_index_per_cell .== block_index)
        n_block_cells = length(block_cell_indices)
        if n_block_cells == 0
            local_clusters_per_block[block_index] = nothing
            #counter = atomic_add!(progress_counter, 1)
            #@debug "- Block: $(block_name) ($(percent(counter + 1, n_blocks))) Empty"
            return nothing
        end

        local_model = load_local_model(daf; block_index, only_pca = true)
        @assert local_model.n_principal_components == n_principal_components_per_block[block_index]

        principal_component_per_block_cell = compute_profiles_principal_components(;
            local_model,
            scaled_log_covered_fraction_per_skeleton_per_profile = scaled_log_covered_fraction_per_skeleton_per_cell,
            profile_indices = block_cell_indices,
        )

        n_block_clusters = max(Int(round(n_block_cells / mean_metacell_cells_per_block[block_index])), 1)

        kmeans_result = kmeans_with_sizes(
            principal_component_per_block_cell,
            n_block_clusters;
            min_cluster_size = min_metacell_size,
            max_cluster_size = mean_metacell_cells_per_block[block_index],
            kmeans_rounds,
            rng,
        )

        cluster_sizes = counts(kmeans_result)
        return local_clusters_per_block[block_index] = LocalClusters(;
            block_cell_indices,
            cluster_index_per_block_cell = assignments(kmeans_result),
            is_too_small_per_cluster = cluster_sizes .< min_metacell_size,
        )

        #       n_clusters = length(cluster_sizes)
        #       n_large_clusters = sum(cluster_sizes .> 2 * mean_metacell_cells_per_block[block_index])
        #       n_small_clusters = sum(cluster_sizes .< min_metacell_size)
        #       counter = atomic_add!(progress_counter, 1)
        #       @debug (
        #           "- Block: $(block_name) ($(percent(counter + 1, n_blocks)))" *
        #           " Clusters: $(n_clusters)" *
        #           " Min: $(minimum(cluster_sizes))" *
        #           " Mean: $(mean(cluster_sizes))" *
        #           " Max: $(maximum(cluster_sizes))" *
        #           " Too small: $(n_small_clusters)" *
        #           " Too large: $(n_large_clusters)"
        #       )
    end

    return local_clusters_per_block
end

function combine_local_clusters(;
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    n_old_metacells::Integer,
    n_cells::Integer,
)::Tuple{AbstractVector{<:Integer}, <:AbstractVector{<:Integer}}
    n_total_new_metacells = 0
    n_outlier_new_metacells = 0
    n_new_outlier_cells = 0
    for local_clusters in local_clusters_per_block
        if local_clusters !== nothing
            n_outlier_new_metacells += sum(local_clusters.is_too_small_per_cluster)
            n_total_new_metacells += length(local_clusters.is_too_small_per_cluster)
            for too_small_cluster_index in findall(local_clusters.is_too_small_per_cluster)
                n_new_outlier_cells += sum(local_clusters.cluster_index_per_block_cell .== too_small_cluster_index)
            end
        end
    end

    new_metacell_index_per_cell = zeros(UInt32, n_cells)
    block_index_per_new_metacell = Vector{UInt32}(undef, n_total_new_metacells - n_outlier_new_metacells)
    n_new_metacells = 0

    for (block_index, local_clusters) in enumerate(local_clusters_per_block)
        if local_clusters !== nothing
            n_clusters = length(local_clusters.is_too_small_per_cluster)
            for cluster_index in 1:n_clusters
                if !local_clusters.is_too_small_per_cluster[cluster_index]
                    n_new_metacells += 1
                    @views cell_indices_of_new_metacell =
                        local_clusters.block_cell_indices[local_clusters.cluster_index_per_block_cell .== cluster_index]
                    new_metacell_index_per_cell[cell_indices_of_new_metacell] .= n_new_metacells
                    block_index_per_new_metacell[n_new_metacells] = block_index
                end
            end
        end
    end

    @assert n_new_metacells == n_total_new_metacells - n_outlier_new_metacells

    @debug (
        "Metacells Old: $(n_old_metacells)" *
        " New: $(n_old_metacells)" *
        " Outliers: $(n_outlier_new_metacells)" *
        " Cells: $(n_new_outlier_cells)" *
        " ($(percent(n_outlier_new_metacells, n_cells)))"
    )

    return (new_metacell_index_per_cell, block_index_per_new_metacell)
end

function load_local_model(daf::DafReader; block_index::Integer, only_pca::Bool)::LocalModel
    block_name = axis_vector(daf, "block")[block_index]
    n_principal_components = daf["/ block = $(block_name) : n_principal_components"]

    mean_scaled_log_covered_fraction_per_skeleton =
        daf["/ block = $(block_name) / gene & is_skeleton : mean_scaled_log_covered_fraction"].array

    coefficient_per_skeleton_per_principal_component =
        daf["/ gene & is_skeleton / principal_component & is_used ; block = $(block_name) : $(block_name)_skeleton_coefficient"].array

    if only_pca
        mean_scaled_log_covered_fraction_per_covered = nothing
        coefficient_per_principal_component_per_covered = nothing
    else
        mean_scaled_log_covered_fraction_per_covered =
            daf["/ block = $(block_name) / gene & is_covered : mean_scaled_log_covered_fraction"].array

        coefficient_per_principal_component_per_covered =
            daf["/ principal_component & is_used ; block = $(block_name) / gene & is_covered : $(block_name)_covered_coefficient"].array
    end

    return LocalModel(;
        n_principal_components,
        mean_scaled_log_covered_fraction_per_skeleton,
        coefficient_per_skeleton_per_principal_component,
        mean_scaled_log_covered_fraction_per_covered,
        coefficient_per_principal_component_per_covered,
    )
end

function compute_profiles_principal_components(;
    local_model::LocalModel,
    scaled_log_covered_fraction_per_skeleton_per_profile::AbstractMatrix{<:AbstractFloat},
    profile_indices::AbstractVector{<:Integer},
)::AbstractMatrix{<:AbstractFloat}
    n_profiles = length(profile_indices)

    @views scaled_log_covered_fraction_per_skeleton_per_profile =
        scaled_log_covered_fraction_per_skeleton_per_profile[:, profile_indices]

    offset_scaled_log_covered_fraction_per_skeleton_per_profile =
        scaled_log_covered_fraction_per_skeleton_per_profile .-
        local_model.mean_scaled_log_covered_fraction_per_skeleton

    principal_component_per_profile =
        transpose(local_model.coefficient_per_skeleton_per_principal_component) *
        offset_scaled_log_covered_fraction_per_skeleton_per_profile

    @assert_matrix(principal_component_per_profile, local_model.n_principal_components, n_profiles)
    return principal_component_per_profile
end

function kmeans_with_sizes(
    values_of_points::AbstractMatrix{<:AbstractFloat},
    initial_k::Integer;
    min_cluster_size::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    rng::AbstractRNG,
)::KmeansResult
    kmeans_result = nothing
    max_k = size(values_of_points, 2)
    k = min(initial_k, max_k)
    centers = nothing
    cluster_sizes = nothing

    while true
        kmeans_result = kmeans_in_rounds(values_of_points, k; centers, kmeans_rounds, rng)
        cluster_sizes = counts(kmeans_result)

        @assert length(cluster_sizes) == k
        largest_cluster_size = maximum(cluster_sizes)
        if largest_cluster_size <= max_cluster_size
            break
        end

        n_large_clusters = sum(cluster_sizes .> max_cluster_size)
        n_split_clusters = max(1, div(n_large_clusters + 1, 2))
        indices_of_split_clusters = partialsortperm(cluster_sizes, 1:n_split_clusters; rev = true)

        clusters_of_points = assignments(kmeans_result)
        centers = kmeans_result.centers

        if k == max_k
            break
        end

        new_centers = Vector{AbstractVector{<:AbstractFloat}}()
        for split_cluster_index in indices_of_split_clusters
            indices_of_points_in_split_cluster = findall(clusters_of_points .== split_cluster_index)
            Random.seed!(rng, 123456)
            split_result = kmeans(values_of_points[:, indices_of_points_in_split_cluster], 2; rng)
            centers[:, split_cluster_index] .= split_result.centers[:, 1]
            push!(new_centers, split_result.centers[:, 2])
            k += 1
        end
        centers = hcat(centers, new_centers...)
    end

    while k > 1
        @assert length(cluster_sizes) == k
        smallest_cluster_index = argmin(cluster_sizes)
        smallest_cluster_size = cluster_sizes[smallest_cluster_index]

        if smallest_cluster_size >= min_cluster_size
            return kmeans_result
        end

        centers = kmeans_result.centers[:, 1:k .!= smallest_cluster_index]

        k -= 1
        merged_kmeans_result = kmeans_in_rounds(values_of_points, k; centers, kmeans_rounds, rng)
        cluster_sizes = counts(merged_kmeans_result)
        largest_cluster_size = maximum(cluster_sizes)
        if largest_cluster_size > max_cluster_size
            k += 1
            break
        end

        kmeans_result = merged_kmeans_result
    end

    @assert length(counts(kmeans_result)) == k
    return kmeans_result
end

function kmeans_in_rounds(
    values_of_points::AbstractMatrix{<:AbstractFloat},
    k::Integer;
    centers::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing,
    kmeans_rounds::Integer,
    rng::AbstractRNG,
)::KmeansResult
    best_kmeans_result = nothing

    for _ in 1:kmeans_rounds
        if centers === nothing
            kmeans_result = kmeans(values_of_points, k; rng)
        else
            kmeans_result = kmeans!(values_of_points, copy_array(centers); rng)
        end

        if best_kmeans_result === nothing || kmeans_result.totalcost < best_kmeans_result.totalcost
            best_kmeans_result = kmeans_result
        end
    end

    @assert best_kmeans_result !== nothing
    return best_kmeans_result
end

end  # module
