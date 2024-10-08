"""
Approximate the manifold of actual cell states (captured by metacells) using linear programs in each local region.
"""
module Programs

export compute_blocks!

using ..Contracts
using ..Defaults
using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using DataAxesFormats
using DataAxesFormats.GenericLogging
using DataAxesFormats.GenericTypes
using Distributions
using MultivariateStats
using NonNegLeastSquares
using Printf
using Random
using SparseArrays
using Statistics

import Random.default_rng

Indices = AbstractVector{<:Integer}
Mask = Union{AbstractVector{Bool}, BitVector}
Subset = Union{Indices, Mask}

@kwdef struct Blocks
    blocks_of_metacells::Indices
    metacells_of_blocks::AbstractVector{<:Indices}
    n_blocks::Integer
    distances_between_blocks::AbstractMatrix{<:AbstractFloat}
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
        min_blocks_in_approximate_environment::Integer = $(DEFAULT.min_blocks_in_approximate_environment),
        min_metacells_in_approximate_environment::Integer = $(DEFAULT.min_metacells_in_approximate_environment),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group the metacells into blocks where the metacells in each one are "similar" in all the transcription factors. That is,
each block is an approximation of a single cell state, assuming the transcription factors control the expression of the
rest of the genes. Then, group these blocks into overlapping small tight neighborhoods, and larger environments, such
that each block's environment can be reasonably approximated using a linear model (which is not computed), evaluated
using cross-validation on the RMSE in the neighborhood at its core.

First, we compute metacell-metacell distances based on the maximal fold factor between `is_transcription_factor` genes,
which aren't `is_forbidden_factor`. The fold factor is log (base 2) of the gene expression using the
`gene_fraction_regularization`. For computing this fold factor, we ignore genes whose total UMIs in the compared
metacells isn't at least `min_significant_gene_UMIs`. We also we reduce the distance using the `fold_confidence` based
on the number of UMIs used to estimate the expression in the metacells, and the `divergence` of the gene. Two metacells
can only belong to the same block if the final fold factor is at most `max_block_span` in all the genes.

We then compute block-block distance which is the mean distance between the metacells of the blocks. Using this
distances, we define a tight neighborhood of each block containing at least `min_blocks_in_neighborhood` and
`min_metacells_in_neighborhood`. We also compute an approximate environment of at least
`min_blocks_in_approximate_environment` and `min_metacells_in_approximate_environment`. Using this approximate
environment, we compute `max_principal_components` for the transcription factor genes, and figure out how many of them
are useful to predict the value of the rest of the genes using `cross_validation_parts` minimizing the RMSE
in the neighborhood.

Having determined the dimensionality of the environment of the block, we search for the optimal environment such that
the RMSE of the cross validation of a linear model with this humber of principal components (computed on the transcription
factor genes) for predicting the value of the rest of the genes.

$(CONTRACT)
"""
@logged @computation Contract(;
    is_relaxed = true,
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        metacell_total_UMIs_vector(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_transcription_factor_vector(RequiredInput),
        gene_is_forbidden_factor_vector(OptionalInput),
        gene_metacell_fraction_matrix(RequiredInput),
        gene_metacell_total_UMIs_matrix(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
    ],
) function compute_blocks!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    max_block_span::Real = max(function_default(identify_marker_genes!, :min_marker_gene_range_fold) - 1, 1),
    max_principal_components::Integer = 40,
    min_blocks_in_neighborhood::Integer = 4,
    min_metacells_in_neighborhood::Integer = 50,
    min_blocks_in_approximate_environment::Integer = 20,
    min_metacells_in_approximate_environment::Integer = 500,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_significant_gene_UMIs > 0
    @assert 0 <= fold_confidence <= 1
    @assert max_block_span > 0
    @assert max_principal_components > 0
    @assert min_blocks_in_neighborhood > 0
    @assert min_metacells_in_neighborhood > 0
    @assert min_blocks_in_approximate_environment > 0
    @assert min_metacells_in_approximate_environment > 0
    @assert cross_validation_parts > 1

    n_genes = axis_length(daf, "gene")
    @assert n_genes > 0
    @debug "Genes: $(n_genes)"

    names_of_genes = axis_array(daf, "gene")

    n_metacells = axis_length(daf, "metacell")
    @assert n_metacells > 0
    @debug "Metacells: $(n_metacells)"

    divergence_of_genes = get_vector(daf, "gene", "divergence").array

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    factor_genes_mask =
        get_vector(daf, "gene", "is_transcription_factor") .&
        .!get_vector(daf, "gene", "is_forbidden_factor"; default = false)
    factor_genes_indices = findall(factor_genes_mask)
    n_factor_genes = length(factor_genes_indices)
    @debug "Factors $(length(factor_genes_indices)): [ $(join(sort(names_of_genes[factor_genes_indices]), ", ")) ]"
    @assert 0 < n_factor_genes < n_genes

    measured_genes_mask = .!factor_genes_mask
    n_measured_genes = n_genes - n_factor_genes
    @assert n_measured_genes > 0

    total_UMIs_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    total_UMIs_of_metacells = get_vector(daf, "metacell", "total_UMIs").array

    log_fractions_of_genes_in_metacells = log2.(fractions_of_genes_in_metacells .+ gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    scale_of_genes = 1.0 .- divergence_of_genes

    total_UMIs_of_factor_genes_in_metacells = total_UMIs_of_genes_in_metacells[factor_genes_indices, :]
    log_decreased_fractions_of_factor_genes_in_metacells, log_increased_fractions_of_factor_genes_in_metacells =
        compute_confidence(;
            fractions_of_genes_in_metacells = fractions_of_genes_in_metacells[factor_genes_indices, :],
            total_UMIs_of_genes_in_metacells = total_UMIs_of_factor_genes_in_metacells,
            total_UMIs_of_metacells = total_UMIs_of_metacells,
            scale_of_genes = scale_of_genes[factor_genes_indices],
            gene_fraction_regularization = gene_fraction_regularization,
            fold_confidence = fold_confidence,
        )

    blocks = compute_blocks_by_confidence(;
        log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_factor_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_factor_genes_in_metacells,
        total_UMIs_of_genes_in_metacells = total_UMIs_of_factor_genes_in_metacells,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        max_block_span = max_block_span,
    )
    @debug "Blocks: $(blocks.n_blocks)"
    block_names = group_names(daf, "metacell", blocks.metacells_of_blocks; prefix = "B")

    add_axis!(daf, "block", block_names)
    set_vector!(daf, "metacell", "block", block_names[blocks.blocks_of_metacells]; overwrite = overwrite)

    scaled_log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells .* scale_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_metacells, n_genes, n_metacells)

    scaled_log_fractions_of_factor_genes_in_metacells =
        scaled_log_fractions_of_genes_in_metacells[factor_genes_indices, :]
    @assert_matrix(scaled_log_fractions_of_factor_genes_in_metacells, n_factor_genes, n_metacells)

    scaled_log_fractions_of_measured_genes_in_metacells =
        scaled_log_fractions_of_genes_in_metacells[measured_genes_mask, :]
    @assert_matrix(scaled_log_fractions_of_measured_genes_in_metacells, n_measured_genes, n_metacells)

    scaled_log_fractions_of_genes_in_blocks =
        transpose(
            daf["/ metacell / gene : fraction % Log base 2 eps $(gene_fraction_regularization) @ block %> Mean"].array,
        ) .* scale_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_blocks, n_genes, blocks.n_blocks)

    scaled_log_fractions_of_factor_genes_in_blocks = scaled_log_fractions_of_genes_in_blocks[factor_genes_indices, :]
    @assert_matrix(scaled_log_fractions_of_factor_genes_in_blocks, n_factor_genes, blocks.n_blocks)

    scaled_log_fractions_of_measured_genes_in_blocks = scaled_log_fractions_of_genes_in_blocks[measured_genes_mask, :]
    @assert_matrix(scaled_log_fractions_of_measured_genes_in_blocks, n_measured_genes, blocks.n_blocks)

    is_in_neighborhood_of_block_of_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    is_in_environment_of_block_of_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    expansion_dimensions_of_blocks = Vector{UInt16}(undef, blocks.n_blocks)
    expansion_rmse_of_blocks = Vector{Float32}(undef, blocks.n_blocks)

    @threads for block_index in 1:(blocks.n_blocks)
        @views is_in_neighborhood_of_block = is_in_neighborhood_of_block_of_block[:, block_index]
        @views is_in_environment_of_block = is_in_environment_of_block_of_block[:, block_index]
        expansion_dimensions_of_blocks[block_index], expansion_rmse_of_blocks[block_index] =
            compute_pcs_vicinity_of_block!(;
                scaled_log_fractions_of_factor_genes_in_metacells = scaled_log_fractions_of_factor_genes_in_metacells,
                scaled_log_fractions_of_measured_genes_in_metacells = scaled_log_fractions_of_measured_genes_in_metacells,
                distances_between_blocks = blocks.distances_between_blocks,
                min_blocks_in_neighborhood = min_blocks_in_neighborhood,
                min_metacells_in_neighborhood = min_metacells_in_neighborhood,
                min_blocks_in_approximate_environment = min_blocks_in_approximate_environment,
                min_metacells_in_approximate_environment = min_metacells_in_approximate_environment,
                max_principal_components = max_principal_components,
                n_parts = cross_validation_parts,
                rng = rng,
                blocks = blocks,
                block_index = block_index,
                is_in_neighborhood_of_block = is_in_neighborhood_of_block,
                is_in_environment_of_block = is_in_environment_of_block,
            )
    end

    set_vector!(daf, "block", "expansion_dimensions", expansion_dimensions_of_blocks)
    set_vector!(daf, "block", "expansion_rmse", expansion_rmse_of_blocks)

    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_neighborhood",
        SparseMatrixCSC(is_in_neighborhood_of_block_of_block);
        overwrite = overwrite,
    )
    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_environment",
        SparseMatrixCSC(is_in_environment_of_block_of_block);
        overwrite = overwrite,
    )

    return nothing
end

function compute_confidence(;  # untested
    fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    total_UMIs_of_metacells::AbstractVector{<:Unsigned},
    scale_of_genes::AbstractVector{<:AbstractFloat},
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}}
    n_genes, n_metacells = size(fractions_of_genes_in_metacells)
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)
    @assert_matrix(total_UMIs_of_genes_in_metacells, n_genes, n_metacells, Columns)
    @assert_vector(total_UMIs_of_metacells, n_metacells)
    @assert_vector(scale_of_genes, n_genes)

    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_fractions_of_genes_in_metacells =  # NOJET
        confidence_stdevs .* sqrt.(transpose(total_UMIs_of_metacells) .* fractions_of_genes_in_metacells) ./
        transpose(total_UMIs_of_metacells)

    log_decreased_fractions_of_genes_in_metacells =  # NOJET
        log2.(
            max.(fractions_of_genes_in_metacells .- confidence_fractions_of_genes_in_metacells, 0.0) .+
            gene_fraction_regularization
        ) .* scale_of_genes

    log_increased_fractions_of_genes_in_metacells =
        log2.(
            fractions_of_genes_in_metacells .+ confidence_fractions_of_genes_in_metacells .+
            gene_fraction_regularization
        ) .* scale_of_genes

    @assert_matrix(log_decreased_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns,)
    @assert_matrix(log_increased_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns,)

    return (log_decreased_fractions_of_genes_in_metacells, log_increased_fractions_of_genes_in_metacells)
end

function compute_blocks_by_confidence(;  # untested
    log_decreased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    log_increased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    min_significant_gene_UMIs::Integer,
    max_block_span::Real,
)::Blocks
    distances_between_metacells = compute_distances_between_metacells(;
        log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_genes_in_metacells,
        total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
    )

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    blocks_of_metacells = Vector{UInt32}(cutree(clusters; h = max_block_span))
    metacells_of_blocks = collect_group_members(blocks_of_metacells)
    n_blocks = length(metacells_of_blocks)

    distances_between_blocks = compute_distances_between_blocks(;
        distances_between_metacells = distances_between_metacells,
        metacells_of_blocks = metacells_of_blocks,
    )

    return Blocks(;
        blocks_of_metacells = blocks_of_metacells,
        metacells_of_blocks = metacells_of_blocks,
        n_blocks = n_blocks,
        distances_between_blocks = distances_between_blocks,
    )
end

function compute_distances_between_metacells(;  # untested
    log_decreased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    log_increased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    min_significant_gene_UMIs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_genes, n_metacells = size(log_decreased_fractions_of_genes_in_metacells)
    @assert_matrix(log_decreased_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)
    @assert_matrix(log_increased_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)

    distances_between_metacells[1, 1] = 0.0
    @threads for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views total_UMIs_of_genes_in_base_metacell = vec(total_UMIs_of_genes_in_metacells[:, base_metacell_index])
        @views log_decreased_fractions_of_genes_in_base_metacell =
            vec(log_decreased_fractions_of_genes_in_metacells[:, base_metacell_index])
        @views log_increased_fractions_of_genes_in_base_metacell =
            vec(log_increased_fractions_of_genes_in_metacells[:, base_metacell_index])

        @views total_UMIs_of_genes_in_other_metacells = total_UMIs_of_genes_in_metacells[:, 1:(base_metacell_index - 1)]
        @views log_decreased_fractions_of_genes_in_other_metacells =
            log_decreased_fractions_of_genes_in_metacells[:, 1:(base_metacell_index - 1)]
        @views log_increased_fractions_of_genes_in_other_metacells =
            log_increased_fractions_of_genes_in_metacells[:, 1:(base_metacell_index - 1)]

        significant_folds_of_genes_in_other_metacells =
            gene_distance.(
                min_significant_gene_UMIs,
                total_UMIs_of_genes_in_base_metacell,
                log_decreased_fractions_of_genes_in_base_metacell,
                log_increased_fractions_of_genes_in_base_metacell,
                total_UMIs_of_genes_in_other_metacells,
                log_decreased_fractions_of_genes_in_other_metacells,
                log_increased_fractions_of_genes_in_other_metacells,
            )
        @assert_matrix(significant_folds_of_genes_in_other_metacells, n_genes, base_metacell_index - 1, Columns,)

        distances_between_base_and_other_metacells =
            vec(maximum(significant_folds_of_genes_in_other_metacells; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, base_metacell_index - 1)

        distances_between_metacells[1:(base_metacell_index - 1), base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, 1:(base_metacell_index - 1)] .=
            distances_between_base_and_other_metacells
    end

    return distances_between_metacells
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Integer,
    total_UMIs_of_gene_in_base_metacell::Integer,
    log_decreased_fractions_of_gene_in_base_metacell::AbstractFloat,
    log_increased_fractions_of_gene_in_base_metacell::AbstractFloat,
    total_UMIs_of_genes_in_other_metacell::Integer,
    log_decreased_fractions_of_gene_in_other_metacell::AbstractFloat,
    log_increased_fractions_of_gene_in_other_metacell::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = total_UMIs_of_gene_in_base_metacell + total_UMIs_of_genes_in_other_metacell
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = log_increased_fractions_of_gene_in_base_metacell < log_increased_fractions_of_gene_in_other_metacell

    log_increased_low_fractions_of_gene =
        is_base_low * log_increased_fractions_of_gene_in_base_metacell +
        !is_base_low * log_increased_fractions_of_gene_in_other_metacell

    log_decreased_high_fractions_of_gene =
        is_base_low * log_decreased_fractions_of_gene_in_other_metacell +
        !is_base_low * log_decreased_fractions_of_gene_in_base_metacell

    return (is_significant * max.((log_decreased_high_fractions_of_gene - log_increased_low_fractions_of_gene), 0.0))
end

function compute_distances_between_blocks(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacells_of_blocks::AbstractVector{<:Indices},
)::AbstractMatrix{<:AbstractFloat}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacells_of_blocks)
    mean_distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    mean_distances_between_blocks[1, 1] = -1.0
    @threads for base_block_index in reverse(2:n_blocks)
        mean_distances_between_blocks[base_block_index, base_block_index] = -1.0
        metacells_of_base_block = metacells_of_blocks[base_block_index]

        distance_of_metacells_from_base_block_metacells = distances_between_metacells[:, metacells_of_base_block]

        distance_of_metacells_from_base_block = vec(mean(distance_of_metacells_from_base_block_metacells; dims = 2))
        @assert length(distance_of_metacells_from_base_block) == n_metacells

        for other_block_index in 1:(base_block_index - 1)
            metacells_of_other_block = metacells_of_blocks[other_block_index]

            @views distance_of_other_block_metacells_from_base_block =
                distance_of_metacells_from_base_block[metacells_of_other_block]

            distance_between_other_and_base_block = mean(distance_of_other_block_metacells_from_base_block)

            mean_distances_between_blocks[base_block_index, other_block_index] = distance_between_other_and_base_block
            mean_distances_between_blocks[other_block_index, base_block_index] = distance_between_other_and_base_block
        end
    end

    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)
    @threads for base_block_index in 1:n_blocks
        @views mean_distance_between_others_and_base_block = mean_distances_between_blocks[:, base_block_index]
        rank_of_others_for_base_block = invperm(sortperm(mean_distance_between_others_and_base_block))
        distances_between_blocks[:, base_block_index] .= rank_of_others_for_base_block
    end

    distances_between_blocks .*= transpose(distances_between_blocks)
    distances_between_blocks .-= 1
    distances_between_blocks ./= maximum(distances_between_blocks)
    @assert minimum(distances_between_blocks) == 0

    return distances_between_blocks
end

function compute_pcs_vicinity_of_block!(;  # untested
    scaled_log_fractions_of_factor_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    scaled_log_fractions_of_measured_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    distances_between_blocks::AbstractMatrix{<:AbstractFloat},
    min_blocks_in_neighborhood::Integer,
    min_metacells_in_neighborhood::Integer,
    min_blocks_in_approximate_environment::Integer,
    min_metacells_in_approximate_environment::Integer,
    max_principal_components::Integer,
    n_parts::Integer,
    rng::AbstractRNG,
    blocks::Blocks,
    block_index::Integer,
    is_in_neighborhood_of_block::AbstractVector{Bool},
    is_in_environment_of_block::AbstractVector{Bool},
)::Tuple{Integer, AbstractFloat}
    distances_between_others_and_block = distances_between_blocks[:, block_index]
    @assert_vector(distances_between_others_and_block, blocks.n_blocks)

    ordered_block_indices = sortperm(distances_between_others_and_block)
    @assert ordered_block_indices[1] == block_index

    @debug "Block: $(block_index) metacells: $(sum(region_metacells_mask(blocks, ordered_block_indices[1:1])))"

    n_blocks_in_neighborhood = compute_vicinity_of_block(;
        blocks = blocks,
        min_blocks_in_vicinity = min_blocks_in_neighborhood,
        min_metacells_in_vicinity = min_metacells_in_neighborhood,
        ordered_block_indices = ordered_block_indices,
    )
    neighborhood_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood])
    @debug "- Approximate neighborhood blocks: $(n_blocks_in_neighborhood) metacells: $(sum(neighborhood_metacells_mask))"

    n_blocks_in_environment = compute_vicinity_of_block(;
        blocks = blocks,
        min_blocks_in_vicinity = min_blocks_in_approximate_environment,
        min_metacells_in_vicinity = min_metacells_in_approximate_environment,
        ordered_block_indices = ordered_block_indices,
    )
    environment_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_environment])
    @debug "- Approximate environment blocks: $(n_blocks_in_environment) metacells: $(sum(environment_metacells_mask))"

    pca, n_principal_components = compute_environment_principal_components(;
        scaled_log_fractions_of_factor_genes_in_metacells = scaled_log_fractions_of_factor_genes_in_metacells,
        scaled_log_fractions_of_measured_genes_in_metacells = scaled_log_fractions_of_measured_genes_in_metacells,
        neighborhood_metacells_mask = neighborhood_metacells_mask,
        environment_metacells_mask = environment_metacells_mask,
        max_principal_components = max_principal_components,
        n_parts = n_parts,
        rng = rng,
    )

    n_blocks_in_environment, rmse = quad_minimize_cost(;
        minimal_value = n_blocks_in_neighborhood,
        maximal_value = blocks.n_blocks,
        max_uncertainty = 1,
    ) do m_blocks_in_environment
        environment_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:m_blocks_in_environment])
        scaled_log_fractions_of_factor_genes_in_environment_metacells =
            scaled_log_fractions_of_factor_genes_in_metacells[:, environment_metacells_mask]
        values_of_principal_components_of_environment_metacells =
            predict(pca, scaled_log_fractions_of_factor_genes_in_environment_metacells)

        rmse = compute_cross_validation_rmse_of_pca(;
            values_of_included_profiles_of_measured_genes = transposer(
                scaled_log_fractions_of_measured_genes_in_metacells[:, environment_metacells_mask],
            ),
            values_of_principal_components_of_included_profiles = values_of_principal_components_of_environment_metacells,
            core_profiles_mask = neighborhood_metacells_mask[environment_metacells_mask],
            n_principal_components = n_principal_components,
            n_parts = n_parts,
            rng = rng,
        )
        # @debug "  Environment blocks: $(m_blocks_in_environment) RMSE: $(rmse)"
        return rmse
    end

    environment_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_environment])
    @debug "- Refined environment blocks: $(n_blocks_in_environment) metacells: $(sum(environment_metacells_mask)) RMSE: $(rmse)"

    is_in_neighborhood_of_block[ordered_block_indices[1:n_blocks_in_neighborhood]] .= true
    is_in_environment_of_block[ordered_block_indices[1:n_blocks_in_environment]] .= true

    return (n_principal_components, rmse)
end

function compute_vicinity_of_block(;  # untested
    blocks::Blocks,
    min_blocks_in_vicinity::Real,
    min_metacells_in_vicinity::Real,
    ordered_block_indices::Indices,
)::Integer
    n_metacells_in_vicinity = 0
    n_blocks_in_vicinity = 0
    vicinity_metacells_mask = nothing

    while n_blocks_in_vicinity < min_blocks_in_vicinity || n_metacells_in_vicinity < min_metacells_in_vicinity
        n_blocks_in_vicinity += 1
        vicinity_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_vicinity])
        @assert vicinity_metacells_mask !== nothing
        n_metacells_in_vicinity = sum(vicinity_metacells_mask)
    end

    return n_blocks_in_vicinity
end

function compute_environment_principal_components(;  # untested
    scaled_log_fractions_of_factor_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    scaled_log_fractions_of_measured_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    neighborhood_metacells_mask::Union{AbstractVector{Bool}, BitVector},
    environment_metacells_mask::Union{AbstractVector{Bool}, BitVector},
    max_principal_components::Integer,
    n_parts::Integer,
    rng::AbstractRNG,
)::Tuple{PCA, Integer}
    n_metacells = length(neighborhood_metacells_mask)
    n_factor_genes = size(scaled_log_fractions_of_factor_genes_in_metacells, 1)
    n_measured_genes = size(scaled_log_fractions_of_measured_genes_in_metacells, 1)

    @assert_matrix(scaled_log_fractions_of_factor_genes_in_metacells, n_factor_genes, n_metacells, Columns)
    @assert_matrix(scaled_log_fractions_of_measured_genes_in_metacells, n_measured_genes, n_metacells, Columns)

    scaled_log_fractions_of_factor_genes_in_environment_metacells =
        scaled_log_fractions_of_factor_genes_in_metacells[:, environment_metacells_mask]
    n_environment_metacells = size(scaled_log_fractions_of_factor_genes_in_environment_metacells, 2)

    pca = fit(PCA, scaled_log_fractions_of_factor_genes_in_environment_metacells; maxoutdim = max_principal_components)
    @assert indim(pca) == n_factor_genes
    @assert outdim(pca) <= max_principal_components
    m_principal_components = outdim(pca)

    @debug "- Max principal components: $(m_principal_components)"

    values_of_principal_components_of_environment_metacells =
        predict(pca, scaled_log_fractions_of_factor_genes_in_environment_metacells)
    @assert_matrix(
        values_of_principal_components_of_environment_metacells,
        m_principal_components,
        n_environment_metacells
    )

    scaled_log_fractions_of_environment_metacells_of_measured_genes =
        transposer(scaled_log_fractions_of_measured_genes_in_metacells[:, environment_metacells_mask])

    n_principal_components, rmse = quad_minimize_cost(;
        minimal_value = 1,
        maximal_value = m_principal_components,
        max_uncertainty = 1,
    ) do m_principal_components
        rmse = compute_cross_validation_rmse_of_pca(;
            values_of_included_profiles_of_measured_genes = scaled_log_fractions_of_environment_metacells_of_measured_genes,
            values_of_principal_components_of_included_profiles = values_of_principal_components_of_environment_metacells,
            core_profiles_mask = neighborhood_metacells_mask[environment_metacells_mask],
            n_principal_components = m_principal_components,
            n_parts = n_parts,
            rng = rng,
        )
        # @debug "  Used principal components: $(m_principal_components) RMSE: $(rmse)"
        return rmse
    end

    @debug "- Best principal components: $(n_principal_components) RMSE: $(rmse)"
    return pca, n_principal_components
end

function compute_cross_validation_rmse_of_pca(;  # untested
    values_of_included_profiles_of_measured_genes::AbstractMatrix{<:AbstractFloat},
    values_of_principal_components_of_included_profiles::AbstractMatrix{<:AbstractFloat},
    core_profiles_mask::Union{AbstractVector{Bool}, BitVector},
    n_principal_components::Integer,
    n_parts::Integer,
    rng::AbstractRNG,
)::AbstractFloat
    n_included_profiles, n_measured_genes = size(values_of_included_profiles_of_measured_genes)
    shuffled_core_profiles_indices = shuffle!(rng, findall(core_profiles_mask))
    n_core_profiles = length(shuffled_core_profiles_indices)
    parts_size = n_core_profiles / n_parts

    rmse_of_parts = Vector{Float32}(undef, n_parts)

    @threads for part_index in 1:n_parts
        first_left_out_position = Int(round((part_index - 1) * parts_size)) + 1
        last_left_out_position = Int(round(part_index * parts_size))
        left_out_positions = first_left_out_position:last_left_out_position

        test_profiles_indices = shuffled_core_profiles_indices[left_out_positions]
        n_test_profiles = length(test_profiles_indices)
        train_profiles_mask = fill(true, n_included_profiles)
        train_profiles_mask[test_profiles_indices] .= false

        values_of_principal_components_of_train_profiles =
            values_of_principal_components_of_included_profiles[1:n_principal_components, train_profiles_mask]

        values_of_train_profiles_of_measured_genes =
            values_of_included_profiles_of_measured_genes[train_profiles_mask, :]

        mean_values_of_principal_components_of_train_profiles =
            vec(mean(values_of_principal_components_of_train_profiles; dims = 2))
        @assert_vector(mean_values_of_principal_components_of_train_profiles, n_principal_components)

        mean_values_in_train_profiles_of_measured_genes =
            vec(mean(values_of_included_profiles_of_measured_genes; dims = 1))
        @assert_vector(mean_values_in_train_profiles_of_measured_genes, n_measured_genes)

        values_of_principal_components_of_train_profiles .-= mean_values_of_principal_components_of_train_profiles
        values_of_train_profiles_of_measured_genes .-= transpose(mean_values_in_train_profiles_of_measured_genes)

        coefficients_of_principal_components_of_measured_genes =
            transpose(values_of_principal_components_of_train_profiles) \ values_of_train_profiles_of_measured_genes
        @assert_matrix(coefficients_of_principal_components_of_measured_genes, n_principal_components, n_measured_genes)

        values_of_principal_components_of_test_profiles =
            values_of_principal_components_of_included_profiles[1:n_principal_components, test_profiles_indices]

        values_of_measured_genes_in_test_profiles =
            transposer(values_of_included_profiles_of_measured_genes[test_profiles_indices, :])

        predicted_values_of_measured_genes_in_test_profiles =
            transpose(coefficients_of_principal_components_of_measured_genes) *
            values_of_principal_components_of_test_profiles
        @assert_matrix(predicted_values_of_measured_genes_in_test_profiles, n_measured_genes, n_test_profiles, Columns)

        predicted_values_of_measured_genes_in_test_profiles .+= mean_values_in_train_profiles_of_measured_genes
        predicted_values_of_measured_genes_in_test_profiles .-= values_of_measured_genes_in_test_profiles
        predicted_values_of_measured_genes_in_test_profiles .*= predicted_values_of_measured_genes_in_test_profiles

        rmse_of_measured_genes = sqrt.(vec(mean(predicted_values_of_measured_genes_in_test_profiles; dims = 2)))
        @assert_vector(rmse_of_measured_genes, n_measured_genes)
        rmse_of_parts[part_index] = mean(rmse_of_measured_genes)
    end

    return mean(rmse_of_parts)
end

function region_metacells_mask(blocks::Blocks, block_indices::Indices)::Mask  # untested
    metacells_mask = zeros(Bool, length(blocks.blocks_of_metacells))
    for block_index in block_indices
        metacells_mask[blocks.metacells_of_blocks[block_index]] .= true
    end
    return metacells_mask
end

function quad_minimize_cost(  # untested
    compute_cost::Function;
    minimal_value::Integer,
    maximal_value::Integer,
    max_uncertainty::Integer,
)::Tuple{Integer, AbstractFloat}
    @assert 0 < minimal_value < maximal_value

    max_not_u_shaped_samples = 10
    sample_radius_fraction = 0.1
    min_estimate_samples = 10

    @assert 0 < sample_radius_fraction < 1

    sample_values = [minimal_value]
    sample_costs = [compute_cost(minimal_value)]

    is_u_shaped = false

    estimated_value = nothing
    estimated_cost = nothing

    offset = 1
    while !is_u_shaped
        sample_value = minimal_value + offset
        if sample_value <= maximal_value
            offset *= 2
        else
            sample_value = rand(minimal_value:maximal_value)
        end
        sample_cost = compute_cost(sample_value)

        push!(sample_values, sample_value)
        push!(sample_costs, sample_cost)

        if length(sample_values) >= 5
            estimated_value, estimated_cost, is_u_shaped = estimate_peak(sample_values, sample_costs)
            if is_u_shaped
                break
            end
            if length(sample_values) > max_not_u_shaped_samples
                error("minimized function is not U-shaped")
            end
        end
    end

    estimated_value = min(estimated_value, maximal_value)
    estimated_value = max(estimated_value, minimal_value)

    sample_radius = max(max_uncertainty, Int(round((maximal_value - minimal_value + 1) * sample_radius_fraction)))

    index = 0
    while true
        previous_estimated_value = estimated_value
        weights = 1.0 ./ (1.0 .+ abs.(sample_values .- previous_estimated_value))
        estimated_value, estimated_cost, is_u_shaped = estimate_peak(sample_values, sample_costs, weights)

        if !is_u_shaped
            index += 1
            if index > max_not_u_shaped_samples
                error("minimized function is not U-shaped")
            end
            sample_value = rand(minimal_value:maximal_value)
            sample_cost = compute_cost(sample_value)
            push!(sample_values, sample_value)
            push!(sample_costs, sample_cost)
            continue
        end

        estimated_value = min(estimated_value, maximal_value)
        estimated_value = max(estimated_value, minimal_value)

        if abs.(estimated_value - previous_estimated_value) <= max_uncertainty
            n_high_samples = sum(sample_values .>= estimated_value)
            n_low_samples = sum(sample_values .<= estimated_value)
            if n_high_samples >= min_estimate_samples / 2 && n_low_samples >= min_estimate_samples / 2
                return (estimated_value, estimated_cost)
            end
        end

        min_sample_value = max(estimated_value - sample_radius, minimal_value)
        max_sample_value = min(estimated_value + sample_radius, maximal_value)
        sample_value = rand(min_sample_value:max_sample_value)
        sample_cost = compute_cost(sample_value)

        push!(sample_values, sample_value)
        push!(sample_costs, sample_cost)
    end
end

function estimate_peak(  # untested
    xs::AbstractVector{<:Real},
    ys::AbstractVector{<:Real},
    ws::Maybe{AbstractVector{<:Real}} = nothing,
)::Tuple{Integer, AbstractFloat, Bool}
    if ws === nothing
        sum_x = sum(xs)
        sum_y = sum(ys)
        sum_x2 = sum(xs .* xs)
        sum_xy = sum(xs .* ys)  # NOJET
        sum_x3 = sum(xs .* xs .* xs)  # NOJET
        sum_x2y = sum(xs .* xs .* ys)
        sum_x4 = sum(xs .* xs .* xs .* xs)
        n = length(xs)
    else
        sum_x = sum(xs .* ws)
        sum_y = sum(ys .* ws)
        sum_x2 = sum(xs .* xs .* ws)
        sum_xy = sum(xs .* ys .* ws)
        sum_x3 = sum(xs .* xs .* xs .* ws)
        sum_x2y = sum(xs .* xs .* ys .* ws)
        sum_x4 = sum(xs .* xs .* xs .* xs .* ws)
        n = sum(ws)
    end

    sxx = sum_x2 - sum_x^2 / n
    sxy = sum_xy - sum_x * sum_y / n
    sxx2 = sum_x3 - sum_x2 * sum_x / n
    sx2y = sum_x2y - sum_x2 * sum_y / n
    sx2x2 = sum_x4 - sum_x2^2 / n

    d = sxx * sx2x2 - sxx2^2
    b = (sxy * sx2x2 - sx2y * sxx2) / d
    a = (sx2y * sxx - sxy * sxx2) / d
    c = sum_y / n - b * sum_x / n - a * sum_x2 / n

    estimate = -b / (2 * a)
    cost = a * estimate * estimate + b * estimate + c

    return (Int(round(estimate)), cost, a > 0)
end

end  # module
