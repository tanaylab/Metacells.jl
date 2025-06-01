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
import Metacells.Contracts.block_gene_base_covered_fraction_matrix
import Metacells.Contracts.block_covered_UMIs_vector
import Metacells.Contracts.block_n_metacells_vector
import Metacells.Contracts.block_n_used_principal_components_vector
import Metacells.Contracts.block_pca_RMSE_vector
import Metacells.Contracts.block_pca_XRMSE_vector
import Metacells.Contracts.block_principal_component_gene_covered_coefficient_tensor
import Metacells.Contracts.block_principal_component_gene_skeleton_coefficient_tensor
import Metacells.Contracts.block_principal_component_is_used_matrix
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.principal_component_axis

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
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_covered_UMIs_vector(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        metacell_metacell_max_skeleton_fold_distance(GuaranteedOutput),
        block_block_max_skeleton_fold_distance(GuaranteedOutput),
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
    overwrite::Bool = false,
    block_distance_aggregation::Function = maximum,  # TODOX
)::Nothing
    @assert 0 < fold_confidence < 1
    @assert max_block_span > 0

    n_metacells = axis_length(daf, "metacell")

    covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : covered_fraction"].array
    UMIs_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : UMIs"].array
    covered_UMIs_per_metacell = get_vector(daf, "metacell", "covered_UMIs").array

    confidence_intervals = compute_confidence_intervals(;
        gene_fraction_regularization,
        covered_fraction_per_skeleton_per_metacell,
        covered_UMIs_per_metacell,
        fold_confidence,
    )

    distances_between_metacells = compute_distances_between_metacells(;
        confidence_intervals,
        UMIs_per_skeleton_per_metacell,
        min_significant_gene_UMIs,
    )
    set_matrix!(daf, "metacell", "metacell", "max_skeleton_fold_distance", distances_between_metacells)

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    block_index_per_metacell = cutree(clusters; h = max_block_span)
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
    set_matrix!(daf, "block", "block", "max_skeleton_fold_distance", distances_between_blocks)
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

function compute_distances_between_metacells(;  # UNTESTED
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

        @assert distances_between_blocks[base_block_index, base_block_index] <= max_block_span
        if base_block_index > 1
            @views distances_between_base_and_other_blocks =
                distances_between_blocks[1:(base_block_index - 1), base_block_index]
            @assert all(distances_between_base_and_other_blocks .> max_block_span)
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

Group the blocks into overlapping small tight neighborhoods, and larger environments, such that each block's environment
can be reasonably approximated using a local linear model, evaluated using cross-validation on the RMSE in the
neighborhood at its core.

Using the block-block distances, we define a tight neighborhood of each block containing at least
`min_blocks_in_neighborhood`, `min_metacells_in_neighborhood` and `min_covered_UMIs_in_neighborhood` (note that the
latter only counts UMIs of skeleton genes).

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
    axes = [block_axis(RequiredInput)],
    data = [
        block_covered_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_max_skeleton_fold_distance(RequiredInput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_neighborhood!(
    daf::DafWriter;
    min_blocks_in_neighborhood::Integer = 4,
    min_metacells_in_neighborhood::Integer = 20,
    min_covered_UMIs_in_neighborhood::Integer = 2_000_000,
    overwrite::Bool = false,
)::Nothing
    @assert min_blocks_in_neighborhood > 0
    @assert min_covered_UMIs_in_neighborhood >= 0

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    max_skeleton_fold_distance_per_block_per_block = get_matrix(daf, "block", "block", "max_skeleton_fold_distance")
    n_metacells_per_block = get_vector(daf, "block", "n_metacells")
    covered_UMIs_per_block = get_vector(daf, "block", "covered_UMIs")
    is_in_neighborhood_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)

    min_blocks_in_neighborhood = min(min_blocks_in_neighborhood, n_blocks)

    @threads :greedy for base_block_index in 1:n_blocks
        @views distance_from_base_per_block = max_skeleton_fold_distance_per_block_per_block[:, base_block_index]

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
            neighborhood_n_metacells += n_metacells_per_block[next_block_index]
            neighborhood_covered_UMIs += covered_UMIs_per_block[next_block_index]
        end

        @debug (
            "Neighborhood: $(name_per_block[base_block_index])" *
            " Blocks: $(neighborhood_n_blocks)" *
            " Metacells: $(neighborhood_n_metacells)" *
            " Covered M-UMIs: $(neighborhood_covered_UMIs / 1e6)"
        )
    end

    set_matrix!(  # NOJET
        daf,
        "block",
        "block",
        "is_in_neighborhood",
        bestify(is_in_neighborhood_per_other_block_per_base_block);
        overwrite,
    )

    return nothing
end

@kwdef struct EmbeddingModel
    n_full_principal_components::Integer
    base_covered_fraction_per_skeleton::AbstractVector{<:AbstractFloat}
    coefficient_per_skeleton_per_full_principal_component::AbstractMatrix{<:AbstractFloat}
end

@kwdef struct PredictiveModel
    n_used_principal_components::Integer
    base_covered_fraction_per_covered::AbstractVector{<:AbstractFloat} = nothing
    coefficient_per_used_principal_component_per_covered::AbstractMatrix{<:AbstractFloat} = nothing
end

@kwdef struct LocalModel
    embedding_model::EmbeddingModel
    predictive_model::PredictiveModel
    XRMSE::AbstractFloat
    RMSE::AbstractFloat
end

@kwdef struct ReusableMatrices
    skeleton_coefficient_per_gene_per_principal_component::Matrix{Float32}
    covered_coefficient_per_principal_component_per_gene::Matrix{Float32}
end

function TanayLabUtilities.reset_reusable_storage!(reusable_matrices::ReusableMatrices)::Nothing
    reusable_matrices.skeleton_coefficient_per_gene_per_principal_component .= 0
    reusable_matrices.covered_coefficient_per_principal_component_per_gene .= 0
    return nothing
end

"""
    function compute_blocks_is_in_environment!(
        daf::DafWriter;
        max_principal_components = $(DEFAULT.max_principal_components),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Expand the block neighborhoods to larger linear environment. We expand the neighborhood as long as computing a local
linear model for the environment gives a better approximation of the metacells in the neighborhood. This local model
with computing `max_principal_components`. These are assumed to over-fit the solution, so we use cross-validation (in
`cross_validation_parts`) to pick a subset of these principal components to actually use, which hopefully approximates
the true dimensionality of the local data.

!!! note

    This uses the virtual [`metacell_gene_covered_fraction_matrix`](@ref). You will need an `adapter` to map these to
    concrete fractions (geomean, linear, scaled, ...).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        principal_component_axis(GuaranteedOutput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_block_max_skeleton_fold_distance(RequiredInput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
        block_principal_component_is_used_matrix(GuaranteedOutput),
        block_n_used_principal_components_vector(GuaranteedOutput),
        block_gene_base_covered_fraction_matrix(GuaranteedOutput),
        block_principal_component_gene_skeleton_coefficient_tensor(GuaranteedOutput),
        block_principal_component_gene_covered_coefficient_tensor(GuaranteedOutput),
        block_pca_RMSE_vector(GuaranteedOutput),
        block_pca_XRMSE_vector(GuaranteedOutput),
    ],
) function compute_blocks_is_in_environment!(
    daf::DafWriter;
    max_principal_components::Integer = 40,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert max_principal_components > 0
    @assert cross_validation_parts > 1

    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    block_index_per_metacell = daf["/ metacell : block => index"].array
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)

    max_skeleton_fold_distance_per_block_per_block = get_matrix(daf, "block", "block", "max_skeleton_fold_distance")
    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood")

    covered_indices = daf["/ gene & is_covered : index"].array
    skeleton_indices = daf["/ gene & is_skeleton : index"].array

    covered_fraction_per_covered_per_metacell = daf["/ gene & is_covered / metacell : covered_fraction"].array
    covered_fraction_per_skeleton_per_metacell = daf["/ gene & is_skeleton / metacell : covered_fraction"].array

    if overwrite
        delete_axis!(daf, "principal_component")
    end
    add_axis!(daf, "principal_component", ["PC$(index)" for index in 1:max_principal_components]; overwrite)

    pca_RMSE_per_block = Vector{Float32}(undef, n_blocks)
    pca_XRMSE_per_block = Vector{Float32}(undef, n_blocks)
    n_used_principal_components_per_block = Vector{UInt32}(undef, n_blocks)
    is_in_environment_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)
    is_used_per_principal_component_per_block = zeros(Bool, max_principal_components, n_blocks)
    base_covered_fraction_per_gene_per_block = zeros(Float32, n_genes, n_blocks)

    reusable_storage = ReusableStorage() do
        return ReusableMatrices(;
            skeleton_coefficient_per_gene_per_principal_component = zeros(Float32, n_genes, max_principal_components),
            covered_coefficient_per_principal_component_per_gene = zeros(Float32, max_principal_components, n_genes),
        )
    end

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng) do base_block_index, rng
        block_name = name_per_block[base_block_index]

        @views is_in_neighborhood_of_base_per_block =
            is_in_neighborhood_per_other_block_per_base_block[:, base_block_index]

        @views distance_from_base_per_block = max_skeleton_fold_distance_per_block_per_block[:, base_block_index]
        ordered_block_indices = sortperm(distance_from_base_per_block)
        @assert ordered_block_indices[1] == base_block_index

        block_indices_of_neighborhood = findall(is_in_neighborhood_of_base_per_block)
        neighborhood_n_blocks = length(block_indices_of_neighborhood)

        metacell_indices_of_neighborhood = vcat(metacell_indices_per_block[block_indices_of_neighborhood]...)
        neighborhood_cross_validation_indices = pick_cross_validation_indices(;
            full_indices = metacell_indices_of_neighborhood,
            cross_validation_parts,
            rng,
        )

        additional_ordered_block_indices = filter(ordered_block_indices) do block_index
            return !(block_index in block_indices_of_neighborhood)
        end
        @assert length(additional_ordered_block_indices) + neighborhood_n_blocks == n_blocks

        environment_n_blocks, _, environment_model = minimize_cost(;
            minimal_value = neighborhood_n_blocks,
            maximal_value = n_blocks,
            min_significant_cost = GENE_FRACTION_REGULARIZATION / 10,
        ) do environment_m_blocks
            if environment_m_blocks == neighborhood_n_blocks
                additional_block_indices_of_environment = UInt32[]
                additional_metacell_indices_of_environment = UInt32[]
            else
                @views additional_block_indices_of_environment =
                    additional_ordered_block_indices[1:(environment_m_blocks - neighborhood_n_blocks)]
                additional_metacell_indices_of_environment =
                    vcat(metacell_indices_per_block[additional_block_indices_of_environment]...)
            end
            metacell_indices_of_environment =
                vcat(metacell_indices_of_neighborhood, additional_metacell_indices_of_environment)

            # @debug "BASE: $(block_name) BLKS: $(environment_m_blocks)..."

            environment_model = compute_environment_model(
                max_principal_components,
                metacell_indices_of_environment,
                metacell_indices_of_neighborhood,
                neighborhood_cross_validation_indices,
                additional_metacell_indices_of_environment,
                covered_fraction_per_skeleton_per_metacell,
                covered_fraction_per_covered_per_metacell,
            )

            # @debug (
            #     "BASE: $(block_name)" *
            #     " BLKs: $(environment_m_blocks)" *
            #     " PCs: $(environment_model.predictive_model.n_used_principal_components)" *
            #     " RMSE: $(environment_model.RMSE)" *
            #     " XRMSE: $(environment_model.XRMSE)"
            # )
            return (environment_model.XRMSE, environment_model)
        end

        @views block_indices_of_environment = ordered_block_indices[1:environment_n_blocks]

        pca_RMSE_per_block[base_block_index] = environment_model.RMSE
        pca_XRMSE_per_block[base_block_index] = environment_model.XRMSE
        n_used_principal_components_per_block[base_block_index] =
            environment_model.predictive_model.n_used_principal_components
        is_in_environment_per_other_block_per_base_block[block_indices_of_environment, base_block_index] .= true
        is_used_per_principal_component_per_block[
            1:(environment_model.predictive_model.n_used_principal_components),
            base_block_index,
        ] .= true
        base_covered_fraction_per_gene_per_block[covered_indices, base_block_index] =
            environment_model.predictive_model.base_covered_fraction_per_covered

        with_reusable(reusable_storage) do reusable_matrices
            reusable_matrices.skeleton_coefficient_per_gene_per_principal_component[
                skeleton_indices,
                1:(environment_model.embedding_model.n_full_principal_components),
            ] .= environment_model.embedding_model.coefficient_per_skeleton_per_full_principal_component

            reusable_matrices.covered_coefficient_per_principal_component_per_gene[
                1:(environment_model.predictive_model.n_used_principal_components),
                covered_indices,
            ] .= environment_model.predictive_model.coefficient_per_used_principal_component_per_covered

            set_matrix!(
                daf,
                "gene",
                "principal_component",
                "$(block_name)_skeleton_coefficient",
                bestify(reusable_matrices.skeleton_coefficient_per_gene_per_principal_component);
                overwrite,
            )

            set_matrix!(
                daf,
                "principal_component",
                "gene",
                "$(block_name)_covered_coefficient",
                bestify(reusable_matrices.covered_coefficient_per_principal_component_per_gene);
                overwrite,
            )

            counter = atomic_add!(progress_counter, 1)
            @debug (
                "- Environment: $(block_name) ($(percent(counter + 1, n_blocks)))" *
                " Blocks: $(environment_n_blocks) = $(neighborhood_n_blocks) + $(environment_n_blocks - neighborhood_n_blocks)" *
                " Principal components: $(environment_model.predictive_model.n_used_principal_components)" *
                " RMSE: $(environment_model.RMSE)" *
                " XRMSE: $(environment_model.XRMSE)"
            )

            return nothing
        end
    end

    @debug "Mean RMSE: $(mean(pca_RMSE_per_block)) XRMSE: $(mean(pca_XRMSE_per_block))"
    set_vector!(daf, "block", "pca_RMSE", pca_RMSE_per_block; overwrite)
    set_vector!(daf, "block", "pca_XRMSE", pca_XRMSE_per_block; overwrite)
    set_vector!(daf, "block", "n_used_principal_components", n_used_principal_components_per_block; overwrite)
    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_environment",
        bestify(is_in_environment_per_other_block_per_base_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "principal_component",
        "block",
        "is_used",
        bestify(is_used_per_principal_component_per_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "base_covered_fraction",
        bestify(base_covered_fraction_per_gene_per_block);
        overwrite,
    )

    return nothing
end

function minimize_cost(  # UNTESTED
    compute_cost_result::Function;
    minimal_value::Integer,
    maximal_value::Integer,
    min_significant_cost::Real,
)::Tuple{Integer, AbstractFloat, Any}
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

function compute_environment_model(
    max_principal_components::Integer,
    metacell_indices_of_environment::AbstractVector{<:Integer},
    metacell_indices_of_neighborhood::AbstractVector{<:Integer},
    neighborhood_cross_validation_indices::CrossValidationIndices,
    additional_metacell_indices_of_environment::AbstractVector{<:Integer},
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
)::LocalModel
    expanded_train_indices_per_part = [
        vcat(train_indices, additional_metacell_indices_of_environment) for
        train_indices in neighborhood_cross_validation_indices.train_indices_per_part
    ]

    trained_embedding_model_per_part = [
        compute_embedding_model(;
            covered_fraction_per_skeleton_per_metacell,
            metacell_indices = expanded_train_indices_per_part[part_index],
            max_principal_components,
        ) for part_index in 1:neighborhood_cross_validation_indices.n_parts
    ]

    max_principal_components =
        minimum([trained_model.n_full_principal_components for trained_model in trained_embedding_model_per_part])

    environment_n_used_principal_components, XRMSE, _ = minimize_cost(;
        minimal_value = 1,
        maximal_value = max_principal_components,
        min_significant_cost = GENE_FRACTION_REGULARIZATION / 10,
    ) do environment_m_used_principal_components
        XRMSE = mean([
            compute_RMSE_of_model(;
                embedding_model = trained_embedding_model_per_part[part_index],
                predictive_model = compute_predictive_model(;
                    embedding_model = trained_embedding_model_per_part[part_index],
                    covered_fraction_per_skeleton_per_metacell,
                    covered_fraction_per_covered_per_metacell,
                    metacell_indices = expanded_train_indices_per_part[part_index],
                    n_used_principal_components = environment_m_used_principal_components,
                ),
                covered_fraction_per_skeleton_per_metacell,
                covered_fraction_per_covered_per_metacell,
                metacell_indices = neighborhood_cross_validation_indices.test_indices_per_part[part_index],
            ) for part_index in 1:neighborhood_cross_validation_indices.n_parts
        ])
        # @debug "- PCs: $(environment_m_used_principal_components) XRMSE: $(XRMSE)"
        return (XRMSE, nothing)
    end

    # @debug "PCs: $(environment_n_used_principal_components) XRMSE: $(XRMSE)"

    embedding_model = compute_embedding_model(;
        covered_fraction_per_skeleton_per_metacell,
        metacell_indices = metacell_indices_of_environment,
        max_principal_components = environment_n_used_principal_components,
    )
    predictive_model = compute_predictive_model(;
        embedding_model,
        covered_fraction_per_skeleton_per_metacell,
        covered_fraction_per_covered_per_metacell,
        metacell_indices = metacell_indices_of_environment,
        n_used_principal_components = embedding_model.n_full_principal_components,
    )
    RMSE = compute_RMSE_of_model(;
        embedding_model,
        predictive_model,
        covered_fraction_per_skeleton_per_metacell,
        covered_fraction_per_covered_per_metacell,
        metacell_indices = metacell_indices_of_neighborhood,
    )

    return LocalModel(; embedding_model, predictive_model, XRMSE, RMSE)
end

function compute_embedding_model(;
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    metacell_indices::AbstractVector{<:Integer},
    max_principal_components::Integer,
)::EmbeddingModel
    n_skeletons, n_metacells = size(covered_fraction_per_skeleton_per_metacell)
    n_model_metacells = length(metacell_indices)
    n_skeletons = size(covered_fraction_per_skeleton_per_metacell, 1)
    @assert_matrix(covered_fraction_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert n_model_metacells <= n_metacells

    @views covered_fraction_per_skeleton_per_model_metacell =
        covered_fraction_per_skeleton_per_metacell[:, metacell_indices]

    offset_covered_fraction_per_skeleton_per_model_metacell, mean_covered_fraction_per_skeleton =
        centralize(covered_fraction_per_skeleton_per_model_metacell)

    pca = fit(
        PCA,
        offset_covered_fraction_per_skeleton_per_model_metacell;
        mean = 0,
        maxoutdim = max_principal_components,
    )

    coefficient_per_skeleton_per_full_principal_component = projection(pca)
    n_full_principal_components = outdim(pca)
    @assert_matrix(
        coefficient_per_skeleton_per_full_principal_component,
        n_skeletons,
        n_full_principal_components,
        Columns
    )

    return EmbeddingModel(;
        n_full_principal_components,
        base_covered_fraction_per_skeleton = mean_covered_fraction_per_skeleton,
        coefficient_per_skeleton_per_full_principal_component,
    )
end

function compute_predictive_model(;
    embedding_model::EmbeddingModel,
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    metacell_indices::AbstractVector{<:Integer},
    n_used_principal_components::Integer,
)::PredictiveModel
    n_skeletons = size(covered_fraction_per_skeleton_per_metacell, 1)
    n_covered = size(covered_fraction_per_covered_per_metacell, 1)
    n_metacells = size(covered_fraction_per_skeleton_per_metacell, 2)
    @assert_matrix(covered_fraction_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert_matrix(covered_fraction_per_covered_per_metacell, n_covered, n_metacells, Columns)
    @assert n_used_principal_components <= embedding_model.n_full_principal_components

    used_principal_component_per_model_metacell = compute_embedding(;
        embedding_model,
        covered_fraction_per_skeleton_per_profile = covered_fraction_per_skeleton_per_metacell,
        profile_indices = metacell_indices,
        n_used_principal_components,
    )

    @views covered_fraction_per_covered_per_model_metacell =
        covered_fraction_per_covered_per_metacell[:, metacell_indices]

    offset_covered_fraction_per_covered_per_model_metacell, mean_covered_fraction_per_covered_per_model_metacell =
        centralize(covered_fraction_per_covered_per_model_metacell)

    coefficient_per_used_principal_component_per_covered =
        transpose(used_principal_component_per_model_metacell) \
        transpose(offset_covered_fraction_per_covered_per_model_metacell)
    @assert_matrix(
        coefficient_per_used_principal_component_per_covered,
        n_used_principal_components,
        n_covered,
        Columns
    )

    return PredictiveModel(;
        n_used_principal_components,
        base_covered_fraction_per_covered = mean_covered_fraction_per_covered_per_model_metacell,
        coefficient_per_used_principal_component_per_covered,
    )
end

function compute_RMSE_of_model(;
    embedding_model::EmbeddingModel,
    predictive_model::PredictiveModel,
    covered_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    metacell_indices::AbstractVector{<:Integer},
)::AbstractFloat
    n_skeletons = size(covered_fraction_per_skeleton_per_metacell, 1)
    n_covered = size(covered_fraction_per_covered_per_metacell, 1)
    n_metacells = size(covered_fraction_per_skeleton_per_metacell, 2)
    @assert_matrix(covered_fraction_per_skeleton_per_metacell, n_skeletons, n_metacells, Columns)
    @assert_matrix(covered_fraction_per_covered_per_metacell, n_covered, n_metacells, Columns)
    @assert predictive_model.n_used_principal_components <= embedding_model.n_full_principal_components

    predicted_covered_fraction_per_covered_per_model_metacell = compute_prediction(;
        embedding_model,
        predictive_model,
        covered_fraction_per_skeleton_per_profile = covered_fraction_per_skeleton_per_metacell,
        covered_fraction_per_covered_per_profile = covered_fraction_per_covered_per_metacell,
        profile_indices = metacell_indices,
    )

    @views covered_fraction_per_covered_per_model_metacell =
        covered_fraction_per_covered_per_metacell[:, metacell_indices]

    return sqrt(
        mean(
            (
                predicted_covered_fraction_per_covered_per_model_metacell .-
                covered_fraction_per_covered_per_model_metacell
            ) .^ 2,
        ),
    )
end

function centralize(
    value_per_gene_per_metacell::AbstractMatrix{<:AbstractFloat},
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractVector{<:AbstractFloat}}
    n_genes, n_metacells = size(value_per_gene_per_metacell)
    @assert_matrix(value_per_gene_per_metacell, n_genes, n_metacells)

    mean_per_gene = vec(mean(value_per_gene_per_metacell; dims = 2))  # NOLINT
    @assert_vector(mean_per_gene, n_genes)

    return (value_per_gene_per_metacell .- mean_per_gene, mean_per_gene)
end

function compute_embedding(;
    embedding_model::EmbeddingModel,
    covered_fraction_per_skeleton_per_profile::AbstractMatrix{<:AbstractFloat},
    profile_indices::AbstractVector{<:Integer},
    n_used_principal_components::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_skeletons = size(covered_fraction_per_skeleton_per_profile, 1)
    n_profiles = size(covered_fraction_per_skeleton_per_profile, 2)
    n_model_profiles = length(profile_indices)
    @assert_matrix(covered_fraction_per_skeleton_per_profile, n_skeletons, n_profiles, Columns)
    @assert n_used_principal_components <= embedding_model.n_full_principal_components
    @assert n_model_profiles <= n_profiles

    @views covered_fraction_per_skeleton_per_model_profile =
        covered_fraction_per_skeleton_per_profile[:, profile_indices]

    offset_covered_fraction_per_skeleton_per_model_profile =
        covered_fraction_per_skeleton_per_model_profile .- embedding_model.base_covered_fraction_per_skeleton

    @views coefficient_per_skeleton_per_used_principal_component =
        embedding_model.coefficient_per_skeleton_per_full_principal_component[:, 1:n_used_principal_components]

    used_principal_component_per_model_profile =
        transpose(coefficient_per_skeleton_per_used_principal_component) *
        offset_covered_fraction_per_skeleton_per_model_profile

    @assert_matrix(used_principal_component_per_model_profile, n_used_principal_components, n_model_profiles)
    return used_principal_component_per_model_profile
end

function compute_prediction(;
    embedding_model::EmbeddingModel,
    predictive_model::PredictiveModel,
    covered_fraction_per_skeleton_per_profile::AbstractMatrix{<:AbstractFloat},
    covered_fraction_per_covered_per_profile::AbstractMatrix{<:AbstractFloat},
    profile_indices::AbstractVector{<:Integer},
)::AbstractMatrix{<:AbstractFloat}
    n_skeletons = size(covered_fraction_per_skeleton_per_profile, 1)
    n_covered = size(covered_fraction_per_covered_per_profile, 1)
    n_profiles = size(covered_fraction_per_skeleton_per_profile, 2)
    n_model_profiles = length(profile_indices)
    @assert_matrix(covered_fraction_per_skeleton_per_profile, n_skeletons, n_profiles, Columns)
    @assert_matrix(covered_fraction_per_covered_per_profile, n_covered, n_profiles, Columns)
    @assert predictive_model.n_used_principal_components <= embedding_model.n_full_principal_components
    @assert n_model_profiles <= n_profiles

    used_principal_component_per_model_profile = compute_embedding(;
        embedding_model,
        covered_fraction_per_skeleton_per_profile = covered_fraction_per_skeleton_per_profile,
        profile_indices,
        n_used_principal_components = predictive_model.n_used_principal_components,
    )

    @views covered_fraction_per_covered_per_model_profile = covered_fraction_per_covered_per_profile[:, profile_indices]

    predicted_covered_fraction_per_covered_per_model_profile =
        (
            transpose(predictive_model.coefficient_per_used_principal_component_per_covered) *
            used_principal_component_per_model_profile
        ) .+ predictive_model.base_covered_fraction_per_covered

    @assert_matrix(predicted_covered_fraction_per_covered_per_model_profile, n_covered, n_model_profiles)
    return predicted_covered_fraction_per_covered_per_model_profile
end

end  # module
