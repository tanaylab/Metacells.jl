"""
Approximate the manifold of actual cell states (captured by metacells) using linear programs in each local region.
"""
module Programs

export compute_blocks!
export compute_block_programs

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
using NMF
using NonNegLeastSquares
using Printf
using Random
using SparseArrays
using Statistics

import Random.default_rng

@kwdef struct Blocks
    blocks_of_metacells::AbstractVector{<:Integer}
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}}
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
            gene_fraction_regularization,
            fold_confidence,
        )

    blocks = compute_blocks_by_confidence(;
        log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_factor_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_factor_genes_in_metacells,
        total_UMIs_of_genes_in_metacells = total_UMIs_of_factor_genes_in_metacells,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        max_block_span,
    )
    @debug "Blocks: $(blocks.n_blocks)"
    block_names = group_names(daf, "metacell", blocks.metacells_of_blocks; prefix = "B")

    add_axis!(daf, "block", block_names)
    set_vector!(daf, "metacell", "block", block_names[blocks.blocks_of_metacells]; overwrite)

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
                scaled_log_fractions_of_factor_genes_in_metacells,
                scaled_log_fractions_of_measured_genes_in_metacells,
                distances_between_blocks = blocks.distances_between_blocks,
                min_blocks_in_neighborhood,
                min_metacells_in_neighborhood,
                min_blocks_in_approximate_environment,
                min_metacells_in_approximate_environment,
                max_principal_components,
                cross_validation_parts,
                rng,
                blocks,
                block_index,
                is_in_neighborhood_of_block,
                is_in_environment_of_block,
                factor_genes_indices,
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
        overwrite,
    )
    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_environment",
        SparseMatrixCSC(is_in_environment_of_block_of_block);
        overwrite,
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
        log_decreased_fractions_of_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells,
        total_UMIs_of_genes_in_metacells,
        min_significant_gene_UMIs,
    )

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    blocks_of_metacells = Vector{UInt32}(cutree(clusters; h = max_block_span))
    metacells_of_blocks = collect_group_members(blocks_of_metacells)
    n_blocks = length(metacells_of_blocks)

    distances_between_blocks = compute_distances_between_blocks(; distances_between_metacells, metacells_of_blocks)

    return Blocks(; blocks_of_metacells, metacells_of_blocks, n_blocks, distances_between_blocks)
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
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}},
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
    cross_validation_parts::Integer,
    rng::AbstractRNG,
    blocks::Blocks,
    block_index::Integer,
    is_in_neighborhood_of_block::AbstractVector{Bool},
    is_in_environment_of_block::AbstractVector{Bool},
    factor_genes_indices::Any,
)::Tuple{Integer, AbstractFloat}
    n_metacells = size(scaled_log_fractions_of_factor_genes_in_metacells, 2)
    n_factor_genes = length(factor_genes_indices)
    n_measured_genes = size(scaled_log_fractions_of_measured_genes_in_metacells, 1)

    @assert_matrix(scaled_log_fractions_of_factor_genes_in_metacells, n_factor_genes, n_metacells, Columns)
    @assert_matrix(scaled_log_fractions_of_measured_genes_in_metacells, n_measured_genes, n_metacells, Columns)

    distances_between_others_and_block = distances_between_blocks[:, block_index]
    @assert_vector(distances_between_others_and_block, blocks.n_blocks)

    ordered_block_indices = sortperm(distances_between_others_and_block)
    @assert ordered_block_indices[1] == block_index

    @debug "Block: $(block_index) metacells: $(sum(collect_metacells_region_mask(blocks, ordered_block_indices[1:1])))"

    n_blocks_in_neighborhood = compute_vicinity_of_block(;
        blocks,
        min_blocks_in_vicinity = min_blocks_in_neighborhood,
        min_metacells_in_vicinity = min_metacells_in_neighborhood,
        ordered_block_indices,
    )
    neighborhood_metacells_mask =
        collect_metacells_region_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood])
    @debug "- Neighborhood blocks: $(n_blocks_in_neighborhood) metacells: $(sum(neighborhood_metacells_mask))"

    n_blocks_in_approximate_environment = compute_vicinity_of_block(;
        blocks,
        min_blocks_in_vicinity = min_blocks_in_approximate_environment,
        min_metacells_in_vicinity = min_metacells_in_approximate_environment,
        ordered_block_indices,
    )
    approximate_environment_metacells_mask =
        collect_metacells_region_mask(blocks, ordered_block_indices[1:n_blocks_in_approximate_environment])
    scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells =
        scaled_log_fractions_of_factor_genes_in_metacells[:, approximate_environment_metacells_mask]
    n_approximate_environment_metacells =
        size(scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells, 2)
    @debug "- Approximate environment blocks: $(n_blocks_in_approximate_environment) metacells: $(n_approximate_environment_metacells)"

    pca = fit(
        PCA,
        scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells;
        maxoutdim = max_principal_components,
    )
    @assert indim(pca) == n_factor_genes
    @assert outdim(pca) <= max_principal_components
    max_n_principal_components = outdim(pca)
    @debug "- Max principal components: $(max_n_principal_components)"

    coefficients_of_factor_genes_in_principal_components = loadings(pca)  # NOJET
    @assert_matrix(
        coefficients_of_factor_genes_in_principal_components,
        n_factor_genes,
        max_n_principal_components,
        Columns
    )
    mean_scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells =
        vec(mean(scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells; dims = 2))
    @assert_vector(mean_scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells, n_factor_genes)
    scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells .-=
        mean_scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells

    values_of_principal_components_of_approximate_environment_metacells =
        transpose(coefficients_of_factor_genes_in_principal_components) *
        scaled_log_fractions_of_factor_genes_in_approximate_environment_metacells
    @assert_matrix(
        values_of_principal_components_of_approximate_environment_metacells,
        max_n_principal_components,
        n_approximate_environment_metacells
    )

    scaled_log_fractions_of_measured_genes_of_approximate_environment_metacells =
        scaled_log_fractions_of_measured_genes_in_metacells[:, approximate_environment_metacells_mask]

    neighborhood_metacells_mask_in_approximate_environment =
        neighborhood_metacells_mask[approximate_environment_metacells_mask]
    neighborhood_metacells_indices_in_approximate_environment =
        findall(neighborhood_metacells_mask_in_approximate_environment)

    n_principal_components, rmse, _ = fast_minimize_cost(;
        minimal_value = 1,
        maximal_value = max_n_principal_components,
        linear_init = false,
    ) do m_principal_components
        rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = values_of_principal_components_of_approximate_environment_metacells[
                1:m_principal_components,
                :,
            ],
            measured_in_profiles = scaled_log_fractions_of_measured_genes_of_approximate_environment_metacells,
            core_profiles_indices = neighborhood_metacells_indices_in_approximate_environment,
            cross_validation_parts,
            rng,
        )

        @debug "  Used principal components: $(m_principal_components) RMSE: $(rmse)"
        return (rmse, nothing)
    end
    @debug "- Predictive principal components: $(n_principal_components) RMSE: $(rmse)"

    coefficients_of_factor_genes_in_principal_components =
        coefficients_of_factor_genes_in_principal_components[:, 1:n_principal_components]

    n_blocks_in_environment, rmse, _ = fast_minimize_cost(;
        minimal_value = n_blocks_in_neighborhood,
        maximal_value = blocks.n_blocks,
        linear_init = false,
    ) do m_blocks_in_environment
        environment_metacells_mask =
            collect_metacells_region_mask(blocks, ordered_block_indices[1:m_blocks_in_environment])
        scaled_log_fractions_of_factor_genes_in_environment_metacells =
            scaled_log_fractions_of_factor_genes_in_metacells[:, environment_metacells_mask]
        n_environment_metacells = size(scaled_log_fractions_of_factor_genes_in_environment_metacells, 2)

        values_of_principal_components_of_environment_metacells =
            transpose(coefficients_of_factor_genes_in_principal_components) *
            scaled_log_fractions_of_factor_genes_in_environment_metacells
        @assert_matrix(
            values_of_principal_components_of_environment_metacells,
            n_principal_components,
            n_environment_metacells
        )

        scaled_log_fractions_of_measured_genes_in_environment_metacells =
            scaled_log_fractions_of_measured_genes_in_metacells[:, environment_metacells_mask]

        neighborhood_metacells_mask_in_environment = neighborhood_metacells_mask[environment_metacells_mask]
        neighborhood_metacells_indices_in_environment = findall(neighborhood_metacells_mask_in_environment)

        rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = values_of_principal_components_of_environment_metacells,
            measured_in_profiles = scaled_log_fractions_of_measured_genes_in_environment_metacells,
            core_profiles_indices = neighborhood_metacells_indices_in_environment,
            cross_validation_parts,
            rng,
        )
        #@debug "  Environment blocks: $(m_blocks_in_environment) RMSE: $(rmse)"
        return (rmse, nothing)
    end

    environment_metacells_mask = collect_metacells_region_mask(blocks, ordered_block_indices[1:n_blocks_in_environment])
    @debug "- Refined environment blocks: $(n_blocks_in_environment) metacells: $(sum(environment_metacells_mask)) RMSE: $(rmse)"

    is_in_neighborhood_of_block[ordered_block_indices[1:n_blocks_in_neighborhood]] .= true
    is_in_environment_of_block[ordered_block_indices[1:n_blocks_in_environment]] .= true

    return (n_principal_components, rmse)
end

function compute_vicinity_of_block(;  # untested
    blocks::Blocks,
    min_blocks_in_vicinity::Real,
    min_metacells_in_vicinity::Real,
    ordered_block_indices::AbstractVector{<:Integer},
)::Integer
    n_metacells_in_vicinity = 0
    n_blocks_in_vicinity = 0
    vicinity_metacells_mask = nothing

    while n_blocks_in_vicinity < min_blocks_in_vicinity || n_metacells_in_vicinity < min_metacells_in_vicinity
        n_blocks_in_vicinity += 1
        vicinity_metacells_mask = collect_metacells_region_mask(blocks, ordered_block_indices[1:n_blocks_in_vicinity])
        @assert vicinity_metacells_mask !== nothing
        n_metacells_in_vicinity = sum(vicinity_metacells_mask)
    end

    return n_blocks_in_vicinity
end

function collect_metacells_region_mask(
    blocks::Blocks,
    block_indices::AbstractVector{<:Integer},
)::Union{AbstractVector{Bool}, BitVector}  # untested
    region_metacells_mask = zeros(Bool, length(blocks.blocks_of_metacells))
    for block_index in block_indices
        region_metacells_mask[blocks.metacells_of_blocks[block_index]] .= true
    end
    return region_metacells_mask
end

function fast_minimize_cost(  # untested
    compute_cost_result::Function;
    minimal_value::Integer,
    maximal_value::Integer,
    linear_init::Bool,
)::Tuple{Integer, AbstractFloat, Any}
    min_significant_rmse = 1e-3

    @assert min_significant_rmse > 0
    @assert 0 < minimal_value < maximal_value

    if linear_init
        linear_init_steps = 5
        @assert linear_init_steps > 1

        step_size = (maximal_value - minimal_value) / linear_init_steps
        sampled_values = [minimal_value]
        for step_index in 1:linear_init_steps
            sample_value = Int(round(minimal_value + step_size * step_index))
            if sample_value !== sampled_values[end]
                @assert sample_value > sampled_values[end]
                push!(sampled_values, sample_value)
            end
        end
        @assert sampled_values[end] == maximal_value

        sampled_computed = compute_cost_result.(sampled_values)  # NOJET
        sampled_costs = [sample_cost for (sample_cost, _) in sampled_computed]
        sampled_results = [sample_result for (_, sample_result) in sampled_computed]
    else
        sample_cost, sample_result = compute_cost_result(minimal_value)
        sampled_values = [minimal_value]
        sampled_costs = [sample_cost]
        sampled_results = [sample_result]

        offset = 1
        while sampled_values[end] != maximal_value
            sample_value = min(minimal_value + offset, maximal_value)
            offset *= 2
            sample_cost, sample_result = compute_cost_result(sample_value)
            push!(sampled_values, sample_value)
            push!(sampled_costs, sample_cost)
            push!(sampled_results, sample_result)

            if sample_value == maximal_value ||
               (length(sampled_values) > 4 && sampled_costs[end] - minimum(sampled_costs) > min_significant_rmse * 2)
                break
            end
        end
    end

    while true
        best_index = 1
        best_cost = sampled_costs[1]
        n_sampled = length(sampled_costs)
        for index in 2:n_sampled
            cost = sampled_costs[index]
            if cost + min_significant_rmse < best_cost
                best_index = index
                best_cost = cost
            end
        end

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
        insert!(sampled_results, sample_index, sample_result)
    end
end

@kwdef mutable struct BlockContext
    block_index::Integer
    n_factor_genes::Integer
    n_measured_genes::Integer
    n_environment_metacells::Integer
    n_neighborhood_metacells::Integer
    scaled_log_fractions_of_factor_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    central_scaled_log_fractions_of_factor_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    mean_scaled_log_fractions_of_factor_genes_in_environment_metacells::AbstractVector{<:AbstractFloat}
    scaled_log_fractions_of_measured_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    central_scaled_log_fractions_of_measured_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    mean_scaled_log_fractions_of_measured_genes_in_environment_metacells::AbstractVector{<:AbstractFloat}
    neighborhood_metacells_mask_in_environment::Union{AbstractVector{Bool}, BitVector}
    neighborhood_metacells_indices_in_environment::AbstractVector{<:Integer}
    environment_metacells_indices::AbstractVector{<:Integer}
    neighborhood_metacells_indices::AbstractVector{<:Integer}
    factor_genes_indices::AbstractVector{<:Integer}
    measured_genes_indices::AbstractVector{<:Integer}
    names_of_genes::AbstractVector{<:AbstractString}
    cross_validation_parts::Integer
    rng::AbstractRNG
end

@kwdef mutable struct BlockModel
    context::BlockContext
    method::AbstractString
    parameters::AbstractString
    filter::Maybe{AbstractFloat}
    modules::Bool
    names_of_programs::Maybe{AbstractVector{<:AbstractString}} = nothing
    n_programs::Maybe{Integer} = nothing
    n_positive_programs::Maybe{Integer} = nothing
    n_used_factors::Maybe{Integer} = nothing
    n_used_coefficients::Maybe{Integer} = nothing
    coefficients_of_factor_genes_in_programs::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    values_of_programs_in_environment_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    central_values_of_programs_in_environment_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    mean_values_of_programs_in_environment_metacells::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    coefficients_of_programs_of_measured_genes::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    predicted_measured_genes_in_neighborhood_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    g_rmse_of_measured_genes::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    g_rmse::Maybe{AbstractFloat} = nothing
    x_rmse::Maybe{AbstractFloat} = nothing
end

@logged function compute_block_programs(  # untested
    daf::DafReader;
    block_index::Integer,
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    max_principal_components::Integer = 40,
    min_marker_gene_range_fold::Real = max(
        function_default(identify_marker_genes!, :min_marker_gene_range_fold) - 1,
        1,
    ),
    min_marker_gene_max_fraction::AbstractFloat = function_default(
        identify_marker_genes!,
        :min_marker_gene_max_fraction,
    ),
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert max_principal_components > 0
    @assert cross_validation_parts > 1
    @assert 0 < block_index <= axis_length(daf, "block")

    mkdir("csvs/block.$(block_index)")

    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    divergence_of_genes = get_vector(daf, "gene", "divergence").array
    names_of_genes = axis_array(daf, "gene")

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    names_of_blocks = axis_array(daf, "block")
    block_name = names_of_blocks[block_index]

    neighborhood_metacells_indices = daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"]
    environment_metacells_indices = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"]
    n_neighborhood_metacells = length(neighborhood_metacells_indices)
    n_environment_metacells = length(environment_metacells_indices)

    neighborhood_metacells_mask = zeros(Bool, n_metacells)
    environment_metacells_mask = zeros(Bool, n_metacells)

    neighborhood_metacells_mask[neighborhood_metacells_indices] .= true
    environment_metacells_mask[environment_metacells_indices] .= true
    @debug "Metacells in neighborhood: $(n_neighborhood_metacells)"
    @debug "Metacells in environment: $(n_environment_metacells)"

    environment_marker_genes_mask = marker_genes_of_environment(;
        daf,
        gene_fraction_regularization,
        min_marker_gene_range_fold,
        min_marker_gene_max_fraction,
        environment_metacells_mask,
    )
    @debug "Marker genes in environment: $(sum(environment_marker_genes_mask))"

    factor_genes_mask =
        get_vector(daf, "gene", "is_transcription_factor") .&
        .!get_vector(daf, "gene", "is_forbidden_factor"; default = false) .&
        environment_marker_genes_mask
    factor_genes_indices = findall(factor_genes_mask)
    n_factor_genes = length(factor_genes_indices)
    @debug "Factors $(length(factor_genes_indices)): [ $(join(sort(names_of_genes[factor_genes_indices]), ", ")) ]"
    @assert 0 < n_factor_genes < n_genes

    measured_genes_mask = .!factor_genes_mask .& environment_marker_genes_mask
    n_measured_genes = sum(measured_genes_mask)
    measured_genes_indices = findall(measured_genes_mask)
    n_measured_genes = length(measured_genes_indices)
    @debug "Measured genes in environment: $(n_measured_genes)"
    @assert n_measured_genes > 0

    log_fractions_of_genes_in_metacells = log2.(fractions_of_genes_in_metacells .+ gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    scale_of_genes = 1.0 .- divergence_of_genes
    scaled_log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells .* scale_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_metacells, n_genes, n_metacells)

    open("csvs/block.$(block_index)/scale_of_genes.csv", "w") do file
        println(file, "gene,scale")
        for (name_of_gene, scale_of_gene) in zip(names_of_genes, scale_of_genes)
            println(file, "$(name_of_gene),$(scale_of_gene)")
        end
    end

    open("csvs/block.$(block_index)/measured_in_neighborhood.csv", "w") do file
        println(file, "metacell,measured_gene,value")
        for metacell_index in neighborhood_metacells_indices
            for gene_index in measured_genes_indices
                println(
                    file,
                    "$(metacell_index),$(names_of_genes[gene_index]),$(scaled_log_fractions_of_genes_in_metacells[gene_index, metacell_index])",
                )
            end
        end
    end

    open("csvs/block.$(block_index)/factors_in_neighborhood.csv", "w") do file
        println(file, "metacell,factor_gene,value")
        for metacell_index in neighborhood_metacells_indices
            for gene_index in factor_genes_indices
                println(
                    file,
                    "$(metacell_index),$(names_of_genes[gene_index]),$(scaled_log_fractions_of_genes_in_metacells[gene_index, metacell_index])",
                )
            end
        end
    end

    scaled_log_fractions_of_factor_genes_in_environment_metacells =
        scaled_log_fractions_of_genes_in_metacells[factor_genes_indices, environment_metacells_mask]
    @assert_matrix(
        scaled_log_fractions_of_factor_genes_in_environment_metacells,
        n_factor_genes,
        n_environment_metacells
    )

    central_scaled_log_fractions_of_factor_genes_in_environment_metacells,
    mean_scaled_log_fractions_of_factor_genes_in_environment_metacells =
        centralize(scaled_log_fractions_of_factor_genes_in_environment_metacells)

    scaled_log_fractions_of_measured_genes_in_environment_metacells =
        scaled_log_fractions_of_genes_in_metacells[measured_genes_mask, environment_metacells_mask]
    @assert_matrix(
        scaled_log_fractions_of_measured_genes_in_environment_metacells,
        n_measured_genes,
        n_environment_metacells
    )

    central_scaled_log_fractions_of_measured_genes_in_environment_metacells,
    mean_scaled_log_fractions_of_measured_genes_in_environment_metacells =
        centralize(scaled_log_fractions_of_measured_genes_in_environment_metacells)

    neighborhood_metacells_mask_in_environment = neighborhood_metacells_mask[environment_metacells_mask]
    neighborhood_metacells_indices_in_environment = findall(neighborhood_metacells_mask_in_environment)

    context = BlockContext(;
        block_index,
        n_factor_genes,
        n_measured_genes,
        n_environment_metacells,
        n_neighborhood_metacells,
        scaled_log_fractions_of_factor_genes_in_environment_metacells,
        central_scaled_log_fractions_of_factor_genes_in_environment_metacells,
        mean_scaled_log_fractions_of_factor_genes_in_environment_metacells,
        scaled_log_fractions_of_measured_genes_in_environment_metacells,
        central_scaled_log_fractions_of_measured_genes_in_environment_metacells,
        mean_scaled_log_fractions_of_measured_genes_in_environment_metacells,
        neighborhood_metacells_mask_in_environment,
        neighborhood_metacells_indices_in_environment,
        environment_metacells_indices,
        neighborhood_metacells_indices,
        factor_genes_indices,
        measured_genes_indices,
        names_of_genes,
        cross_validation_parts,
        rng,
    )
    store_block_context(context)

    # PCA

    pca = fit(PCA, scaled_log_fractions_of_factor_genes_in_environment_metacells; maxoutdim = max_principal_components)
    @assert indim(pca) == n_factor_genes
    n_pca_principal_components = outdim(pca)
    @assert n_pca_principal_components <= max_principal_components

    pca_coefficients_of_factor_genes_in_principal_components = loadings(pca)
    @assert_matrix(pca_coefficients_of_factor_genes_in_principal_components, n_factor_genes, n_pca_principal_components)

    pca_model = BlockModel(;
        context,
        method = "pca",
        parameters = "all",
        filter = nothing,
        modules = false,
        coefficients_of_factor_genes_in_programs = pca_coefficients_of_factor_genes_in_principal_components,
        n_programs = n_pca_principal_components,
    )
    evaluate_block_model!(pca_model)

    store_block_model(pca_model)
    modules_model(pca_model)
    filter_model(pca_model)

    _, _, subset_pca_model = fast_minimize_cost(;
        minimal_value = 1,
        maximal_value = n_pca_principal_components,
        linear_init = true,
    ) do m_used_principal_components
        subset_pca_model = BlockModel(;
            context,
            method = "pca",
            parameters = "top=$(m_used_principal_components)",
            filter = nothing,
            modules = false,
            coefficients_of_factor_genes_in_programs = pca_coefficients_of_factor_genes_in_principal_components[
                :,
                1:m_used_principal_components,
            ],
            n_programs = m_used_principal_components,
        )
        evaluate_block_model!(subset_pca_model)
        return (subset_pca_model.x_rmse, subset_pca_model)
    end

    store_block_model(subset_pca_model)
    modules_model(subset_pca_model)
    filter_model(subset_pca_model)

    # POLARITY

    polarity_coefficients_of_factor_genes_in_principal_components = hcat(
        max.(pca_coefficients_of_factor_genes_in_principal_components, 0),
        .-min.(pca_coefficients_of_factor_genes_in_principal_components, 0),
    )
    @assert_matrix(
        polarity_coefficients_of_factor_genes_in_principal_components,
        n_factor_genes,
        2 * n_pca_principal_components
    )

    polarity_model = BlockModel(;
        context,
        method = "polarity",
        parameters = "all",
        filter = nothing,
        modules = false,
        coefficients_of_factor_genes_in_programs = polarity_coefficients_of_factor_genes_in_principal_components,
        n_programs = 2 * n_pca_principal_components,
        n_positive_programs = n_pca_principal_components,
    )
    evaluate_block_model!(polarity_model)

    store_block_model(polarity_model)
    modules_model(polarity_model)
    filter_model(polarity_model)

    _, _, subset_polarity_model = fast_minimize_cost(;
        minimal_value = 1,
        maximal_value = n_pca_principal_components,
        linear_init = true,
    ) do n_positive_principal_components
        _, _, subset_polarity_model = fast_minimize_cost(;
            minimal_value = 1,
            maximal_value = n_pca_principal_components,
            linear_init = true,
        ) do n_negative_principal_components
            subset_polarity_model = BlockModel(;
                context,
                method = "polarity",
                parameters = "pos=$(n_positive_principal_components).neg=$(n_negative_principal_components)",
                filter = nothing,
                modules = false,
                coefficients_of_factor_genes_in_programs = polarity_coefficients_of_factor_genes_in_principal_components[
                    :,
                    vcat(
                        collect(1:n_positive_principal_components),
                        collect(
                            (n_pca_principal_components + 1):(n_pca_principal_components + n_negative_principal_components),
                        ),
                    ),
                ],
                n_programs = n_positive_principal_components + n_negative_principal_components,
            )
            evaluate_block_model!(subset_polarity_model)
            return (subset_polarity_model.x_rmse, subset_polarity_model)
        end
        return (subset_polarity_model.x_rmse, subset_polarity_model)
    end

    store_block_model(subset_polarity_model)
    modules_model(subset_polarity_model)
    filter_model(subset_polarity_model)

    # NMF

    rebased_scaled_log_fractions_of_factor_genes_in_environment_metacells =
        scaled_log_fractions_of_factor_genes_in_environment_metacells .- log2(gene_fraction_regularization)

    max_nmf_programs = min(max_principal_components, context.n_factor_genes)
    result = nnmf(rebased_scaled_log_fractions_of_factor_genes_in_environment_metacells, max_nmf_programs)
    nmf_coefficients_of_factor_genes_in_programs = result.W
    @assert_matrix(nmf_coefficients_of_factor_genes_in_programs, n_factor_genes, max_nmf_programs, Columns)

    nmf_model = BlockModel(;
        context,
        method = "nmf",
        parameters = "all",
        filter = nothing,
        modules = false,
        coefficients_of_factor_genes_in_programs = nmf_coefficients_of_factor_genes_in_programs,
        n_programs = max_nmf_programs,
    )
    evaluate_block_model!(nmf_model)

    store_block_model(nmf_model)
    modules_model(nmf_model)
    filter_model(nmf_model)

    _, _, nmf_model =
        fast_minimize_cost(; minimal_value = 1, maximal_value = max_nmf_programs, linear_init = false) do n_nmf_programs
            result = nnmf(rebased_scaled_log_fractions_of_factor_genes_in_environment_metacells, n_nmf_programs)  # NOJET
            nmf_coefficients_of_factor_genes_in_programs = result.W
            @assert_matrix(nmf_coefficients_of_factor_genes_in_programs, n_factor_genes, n_nmf_programs, Columns)

            nmf_model = BlockModel(;
                context,
                method = "nmf",
                parameters = "few=$(n_nmf_programs)",
                filter = nothing,
                modules = false,
                coefficients_of_factor_genes_in_programs = nmf_coefficients_of_factor_genes_in_programs,
                n_programs = n_nmf_programs,
            )
            evaluate_block_model!(nmf_model)
            return (nmf_model.x_rmse, nmf_model)
        end

    store_block_model(nmf_model)
    modules_model(nmf_model)
    filter_model(nmf_model)

    return nothing
end

function filter_model(base_model::BlockModel)::Nothing
    @assert base_model.filter === nothing
    @assert !base_model.modules

    BAD_RMSE = 10.0

    _, _, filtered_model =
        fast_minimize_cost(; minimal_value = 1, maximal_value = 1000, linear_init = true) do threshold_k
            threshold = threshold_k / 1000
            filtered_coefficients_of_factor_genes_in_programs =
                copy_array(base_model.coefficients_of_factor_genes_in_programs)
            filtered_coefficients_of_factor_genes_in_programs[abs.(
                filtered_coefficients_of_factor_genes_in_programs
            ) .< threshold] .= 0

            maximal_coefficients_in_programs =
                vec(maximum(abs.(filtered_coefficients_of_factor_genes_in_programs); dims = 1))
            @assert_vector(maximal_coefficients_in_programs, base_model.n_programs)
            nonzero_programs_mask = maximal_coefficients_in_programs .> 0
            filtered_coefficients_of_factor_genes_in_programs =
                filtered_coefficients_of_factor_genes_in_programs[:, nonzero_programs_mask]
            n_programs = size(filtered_coefficients_of_factor_genes_in_programs, 2)
            if n_programs == 0
                return (BAD_RMSE, nothing)
            end
            if base_model.n_positive_programs === nothing
                n_positive_programs = nothing
            else
                n_positive_programs = sum(nonzero_programs_mask[1:(base_model.n_positive_programs)])  # NOJET
            end

            filtered_model = BlockModel(;
                context = base_model.context,
                method = base_model.method,
                parameters = base_model.parameters,
                filter = threshold,
                modules = false,
                coefficients_of_factor_genes_in_programs = filtered_coefficients_of_factor_genes_in_programs,
                n_programs,
                n_positive_programs,
            )

            evaluate_block_model!(filtered_model)
            return (filtered_model.x_rmse, filtered_model)
        end

    store_block_model(filtered_model)
    return modules_model(filtered_model)
end

function modules_model(base_model::BlockModel)::Nothing
    @assert !base_model.modules

    modules_coefficients_of_factor_genes_in_programs = copy_array(base_model.coefficients_of_factor_genes_in_programs)
    zero_programs_mask = fill(true, base_model.n_programs)  # NOJET
    for factor_gene_position in 1:(base_model.context.n_factor_genes)
        program_index = argmax(vec(abs.(modules_coefficients_of_factor_genes_in_programs[factor_gene_position, :])))
        zero_programs_mask[program_index] = false
        modules_coefficients_of_factor_genes_in_programs[factor_gene_position, zero_programs_mask] .= 0
        zero_programs_mask[program_index] = true
    end

    maximal_coefficients_in_programs = vec(maximum(abs.(modules_coefficients_of_factor_genes_in_programs); dims = 1))
    @assert_vector(maximal_coefficients_in_programs, base_model.n_programs)
    nonzero_programs_mask = maximal_coefficients_in_programs .> 0
    modules_coefficients_of_factor_genes_in_programs =
        modules_coefficients_of_factor_genes_in_programs[:, nonzero_programs_mask]
    n_programs = size(modules_coefficients_of_factor_genes_in_programs, 2)

    if n_programs > 0
        if base_model.n_positive_programs === nothing
            n_positive_programs = nothing
        else
            n_positive_programs = sum(nonzero_programs_mask[1:(base_model.n_positive_programs)])  # NOJET
        end
        modules_model = BlockModel(;
            context = base_model.context,
            method = base_model.method,
            parameters = base_model.parameters,
            filter = base_model.filter,
            modules = true,
            coefficients_of_factor_genes_in_programs = modules_coefficients_of_factor_genes_in_programs,
            n_programs,
            n_positive_programs,
        )
        evaluate_block_model!(modules_model)
        store_block_model(modules_model)
    end

    return nothing
end

function evaluate_block_model!(model::BlockModel)::Nothing
    @assert_matrix(model.coefficients_of_factor_genes_in_programs, model.context.n_factor_genes, model.n_programs)

    if model.n_used_factors === nothing
        factors_coefficients_mask = vec(maximum(abs.(model.coefficients_of_factor_genes_in_programs); dims = 2) .> 0)
        @assert_vector(factors_coefficients_mask, model.context.n_factor_genes)
        model.n_used_factors = sum(factors_coefficients_mask)
    end

    if model.n_used_coefficients === nothing
        model.n_used_coefficients = sum(model.coefficients_of_factor_genes_in_programs .!= 0)
    end

    if model.values_of_programs_in_environment_metacells === nothing
        model.values_of_programs_in_environment_metacells =
            transpose(model.coefficients_of_factor_genes_in_programs) *
            model.context.scaled_log_fractions_of_factor_genes_in_environment_metacells
    end
    @assert_matrix(
        model.values_of_programs_in_environment_metacells,
        model.n_programs,
        model.context.n_environment_metacells
    )

    if model.central_values_of_programs_in_environment_metacells === nothing ||  # NOJET
       model.mean_values_of_programs_in_environment_metacells
        model.central_values_of_programs_in_environment_metacells,  # NOJET
        model.mean_values_of_programs_in_environment_metacells =
            centralize(model.values_of_programs_in_environment_metacells)
    end
    @assert_matrix(
        model.central_values_of_programs_in_environment_metacells,
        model.n_programs,
        model.context.n_environment_metacells
    )

    if model.coefficients_of_programs_of_measured_genes === nothing
        model.coefficients_of_programs_of_measured_genes =
            transpose(model.central_values_of_programs_in_environment_metacells) \
            transpose(model.context.central_scaled_log_fractions_of_measured_genes_in_environment_metacells)
    end
    @assert_matrix(model.coefficients_of_programs_of_measured_genes, model.n_programs, model.context.n_measured_genes)

    if model.predicted_measured_genes_in_neighborhood_metacells === nothing
        model.predicted_measured_genes_in_neighborhood_metacells =
            transpose(model.coefficients_of_programs_of_measured_genes) *
            model.central_values_of_programs_in_environment_metacells[
                :,
                model.context.neighborhood_metacells_indices_in_environment,
            ]
        model.predicted_measured_genes_in_neighborhood_metacells .+=  # NOJET
            model.context.mean_scaled_log_fractions_of_measured_genes_in_environment_metacells
    end
    @assert_matrix(
        model.predicted_measured_genes_in_neighborhood_metacells,
        model.context.n_measured_genes,
        model.context.n_neighborhood_metacells,
    )

    if model.g_rmse_of_measured_genes === nothing
        error_of_measured_genes_in_neighborhood_metacells =
            model.predicted_measured_genes_in_neighborhood_metacells .-
            model.context.scaled_log_fractions_of_measured_genes_in_environment_metacells[
                :,
                model.context.neighborhood_metacells_indices_in_environment,
            ]
        error_of_measured_genes_in_neighborhood_metacells .*= error_of_measured_genes_in_neighborhood_metacells
        model.g_rmse_of_measured_genes = sqrt.(vec(mean(error_of_measured_genes_in_neighborhood_metacells; dims = 2)))
    end
    @assert_vector(model.g_rmse_of_measured_genes, model.context.n_measured_genes)

    if model.g_rmse === nothing
        model.g_rmse = mean(model.g_rmse_of_measured_genes)
    end

    if model.x_rmse === nothing
        model.x_rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = model.values_of_programs_in_environment_metacells,
            measured_in_profiles = model.context.scaled_log_fractions_of_measured_genes_in_environment_metacells,
            core_profiles_indices = model.context.neighborhood_metacells_indices_in_environment,
            cross_validation_parts = model.context.cross_validation_parts,
            rng = model.context.rng,
        )
    end

    return nothing
end

function store_block_context(context::BlockContext)::Nothing
    csv_path_prefix = "csvs/block.$(context.block_index)/"

    open(csv_path_prefix * "factors_in_environment.csv", "w") do file
        println(file, "metacell,factor_gene,value")
        for (metacell_position, metacell_index) in enumerate(context.environment_metacells_indices)
            for (gene_position, gene_index) in enumerate(context.factor_genes_indices)
                gene_name = context.names_of_genes[gene_index]
                value = context.scaled_log_fractions_of_factor_genes_in_environment_metacells[
                    gene_position,
                    metacell_position,
                ]
                println(file, "$(metacell_index),$(gene_name),$(value)")
            end
        end
    end

    open(csv_path_prefix * "measured_in_environment.csv", "w") do file
        println(file, "metacell,measured_gene,value")
        for (metacell_position, metacell_index) in enumerate(context.environment_metacells_indices)
            for (gene_position, gene_index) in enumerate(context.measured_genes_indices)
                gene_name = context.names_of_genes[gene_index]
                value = context.scaled_log_fractions_of_measured_genes_in_environment_metacells[
                    gene_position,
                    metacell_position,
                ]
                println(file, "$(metacell_index),$(gene_name),$(value)")
            end
        end
    end

    return nothing
end

function store_block_model(model::BlockModel)::Nothing
    if !isfile("csvs/methods.csv")
        open("csvs/methods.csv", "a") do file
            return println(
                file,
                "block" *
                ",method" *
                ",params" *
                ",filter" *
                ",unique" *
                ",programs" *
                ",used_factors" *
                ",used_coeffs" *
                ",g_rmse" *
                ",x_rmse",
            )
        end
    end

    open("csvs/methods.csv", "a") do file
        return println(
            file,
            "$(model.context.block_index)" *
            ",$(model.method)" *
            ",$(model.parameters)" *
            ",$(model.filter === nothing ? "1.0" : model.filter)" *
            ",$(model.modules ? "modules" : "programs")" *
            ",$(model.n_programs)" *
            ",$(model.n_used_factors)" *
            ",$(model.n_used_coefficients)" *
            ",$(model.g_rmse)" *
            ",$(model.x_rmse)",
        )
    end

    prefix = "$(model.method).$(model.parameters)"
    if model.filter !== nothing
        prefix *= ".filter=$(model.filter)"
    end
    if model.modules
        prefix *= ".modules"
    end

    @debug "$(prefix) X-RMSE: $(model.x_rmse) programs: $(model.n_programs) factors: $(model.n_used_factors) coeffs: $(model.n_used_coefficients)"

    csv_path_prefix = "csvs/block.$(model.context.block_index)/$(prefix)"
    open(csv_path_prefix * ".factor_programs.csv", "w") do file
        println(file, "block,program,factor_gene,coefficient")
        for program_index in 1:(model.n_programs)  # NOJET
            if model.names_of_programs === nothing
                program_name = "$(program_index)"
            else
                program_name = model.names_of_programs[program_index]
            end
            for gene_position in 1:(model.context.n_factor_genes)
                gene_name = model.context.names_of_genes[model.context.factor_genes_indices[gene_position]]
                value = model.coefficients_of_programs_of_measured_genes[program_index, gene_position]
                println(file, "$(model.context.block_index),$(program_name),$(gene_name),$(value)")
            end
        end
    end

    open(csv_path_prefix * ".measured_programs.csv", "w") do file
        println(file, "block,program,measured_gene,coefficient")
        for program_index in 1:(model.n_programs)  # NOJET
            if model.names_of_programs === nothing
                program_name = "$(program_index)"
            else
                program_name = model.names_of_programs[program_index]
            end
            for gene_position in 1:(model.context.n_measured_genes)
                gene_name = model.context.names_of_genes[model.context.measured_genes_indices[gene_position]]
                value = model.coefficients_of_programs_of_measured_genes[program_index, gene_position]
                println(file, "$(model.context.block_index),$(program_name),$(gene_name),$(value)")
            end
        end
    end

    open(csv_path_prefix * ".g_rmse_of_measured.csv", "w") do file
        println(file, "block,gene,g_rmse")
        for gene_position in 1:(model.context.n_measured_genes)
            gene_name = model.context.names_of_genes[model.context.measured_genes_indices[gene_position]]
            value = model.g_rmse_of_measured_genes[gene_position]
            println(file, "$(model.context.block_index),$(gene_name),$(value)")
        end
    end

    open(csv_path_prefix * ".programs_in_environment.csv", "w") do file
        println(file, "block,metacell,program,value")
        for (metacell_position, metacell_index) in enumerate(model.context.environment_metacells_indices)
            for program_index in 1:(model.n_programs)  # NOJET
                if model.names_of_programs === nothing
                    program_name = "$(program_index)"
                else
                    program_name = model.names_of_programs[program_index]
                end
                value = model.values_of_programs_in_environment_metacells[program_index, metacell_position]
                println(file, "$(model.context.block_index),$(metacell_index),$(program_name),$(value)")
            end
        end
    end

    open(csv_path_prefix * ".programs_in_neighborhood.csv", "w") do file
        println(file, "block,metacell,program,value")
        for metacell_position in model.context.neighborhood_metacells_indices_in_environment
            for program_index in 1:(model.n_programs)  # NOJET
                if model.names_of_programs === nothing
                    program_name = "$(program_index)"
                else
                    program_name = model.names_of_programs[program_index]
                end
                metacell_index = model.context.environment_metacells_indices[metacell_position]
                value = model.values_of_programs_in_environment_metacells[program_index, metacell_position]
                println(file, "$(model.context.block_index),$(metacell_index),$(program_name),$(value)")
            end
        end
    end

    open(csv_path_prefix * ".predicted_in_neighborhood.csv", "w") do file
        println(file, "block,metacell,measured_gene,value")
        for (metacell_position, metacell_index) in enumerate(model.context.neighborhood_metacells_indices)
            for (gene_position, gene_index) in enumerate(model.context.measured_genes_indices)
                metacell_index = model.context.environment_metacells_indices[metacell_position]
                value = model.predicted_measured_genes_in_neighborhood_metacells[gene_position, metacell_position]
                gene_name = model.context.names_of_genes[gene_index]
                println(file, "$(model.context.block_index),$(metacell_index),$(gene_name),$(value)")
            end
        end
    end

    return nothing
end

function centralize(
    values_of_profiles::AbstractMatrix{<:AbstractFloat},
    means_of_values::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractVector{<:AbstractFloat}}
    n_values = size(values_of_profiles, 1)
    if means_of_values === nothing
        means_of_values = vec(mean(values_of_profiles; dims = 2))
    end
    @assert_vector(means_of_values, n_values)
    return (values_of_profiles .- means_of_values, means_of_values)
end

function solve_least_squares(
    central_predictors_in_profiles::AbstractMatrix{<:AbstractFloat},
    central_measured_in_profiles::AbstractMatrix{<:AbstractFloat},
)::AbstractMatrix{<:AbstractFloat}
    n_predictors, n_profiles = size(central_predictors_in_profiles)
    n_measured = size(central_measured_in_profiles, 1)
    @assert_matrix(central_predictors_in_profiles, n_predictors, n_profiles, Columns)
    @assert_matrix(central_measured_in_profiles, n_measured, n_profiles, Columns)
    coefficients_of_predictors_of_measured =
        transpose(central_predictors_in_profiles) \ transpose(central_measured_in_profiles)
    @assert_matrix(coefficients_of_predictors_of_measured, n_predictors, n_measured, Columns)
    return coefficients_of_predictors_of_measured
end

function marker_genes_of_environment(;  # untested
    daf::DafReader,
    gene_fraction_regularization::AbstractFloat,
    min_marker_gene_range_fold::Real,
    min_marker_gene_max_fraction::AbstractFloat,
    environment_metacells_mask::Union{AbstractVector{Bool}, BitVector},
)::Union{AbstractVector{Bool}, BitVector}
    chain = chain_writer([daf, MemoryDaf(; name = "environment")]; name = "mask_chain")
    set_vector!(chain, "metacell", "is_in_environment", environment_metacells_mask; overwrite = true)
    adapter(  # NOJET
        chain;
        input_axes = ["metacell" => "/metacell & is_in_environment", "gene" => "="],
        input_data = [("gene", "divergence") => "=", ("metacell", "gene", "fraction") => "="],
        output_axes = ["gene" => "="],
        output_data = [("gene", "is_marker") => "="],
        overwrite = true,
    ) do adapted
        return identify_marker_genes!(
            adapted;
            gene_fraction_regularization = gene_fraction_regularization,
            min_marker_gene_range_fold = min_marker_gene_range_fold,
            min_marker_gene_max_fraction = min_marker_gene_max_fraction,
            overwrite = true,
        )
    end
    return get_vector(chain, "gene", "is_marker").array
end

function compute_cross_validation_of_prediction(
    solve::Function;
    predictors_in_profiles::AbstractMatrix{<:AbstractFloat},
    measured_in_profiles::AbstractMatrix{<:AbstractFloat},
    core_profiles_indices::AbstractVector{<:Integer},
    cross_validation_parts::Integer,
    rng::AbstractRNG,
)::AbstractFloat
    n_predictors, n_profiles = size(predictors_in_profiles)
    n_measured = size(measured_in_profiles, 1)
    n_core_profiles = length(core_profiles_indices)

    @assert_matrix(predictors_in_profiles, n_predictors, n_profiles, Columns)
    @assert_matrix(measured_in_profiles, n_measured, n_profiles, Columns)
    @assert_vector(core_profiles_indices, n_core_profiles)

    parts_size = n_core_profiles / cross_validation_parts
    rmse_of_parts = Vector{Float32}(undef, cross_validation_parts)
    shuffled_core_profiles_indices = shuffle(rng, core_profiles_indices)
    train_profiles_mask = Vector{Bool}(undef, n_profiles)

    for part_index in 1:cross_validation_parts
        first_left_out_position = Int(round((part_index - 1) * parts_size)) + 1
        last_left_out_position = Int(round(part_index * parts_size))
        left_out_positions = first_left_out_position:last_left_out_position

        test_profiles_indices = shuffled_core_profiles_indices[left_out_positions]
        n_test_profiles = length(test_profiles_indices)

        train_profiles_mask .= true
        train_profiles_mask[test_profiles_indices] .= false
        n_train_profiles = n_profiles - n_test_profiles

        central_predictors_in_train_profiles, means_of_predictors_in_train_profiles =
            centralize(predictors_in_profiles[:, train_profiles_mask])
        @assert_matrix(central_predictors_in_train_profiles, n_predictors, n_train_profiles, Columns)
        @assert_vector(means_of_predictors_in_train_profiles, n_predictors)

        central_measured_in_train_profiles, means_of_measured_in_train_profiles =
            centralize(measured_in_profiles[:, train_profiles_mask])
        @assert_matrix(central_measured_in_train_profiles, n_measured, n_train_profiles, Columns)
        @assert_vector(means_of_measured_in_train_profiles, n_measured)

        coefficients_of_predictors_of_measured =
            solve(central_predictors_in_train_profiles, central_measured_in_train_profiles)
        @assert_matrix(coefficients_of_predictors_of_measured, n_predictors, n_measured, Columns)

        central_predictors_in_test_profiles, _ =
            centralize(predictors_in_profiles[:, test_profiles_indices], means_of_predictors_in_train_profiles)
        @assert_matrix(central_predictors_in_test_profiles, n_predictors, n_test_profiles, Columns)

        predicted_measured_of_test_profiles =
            transpose(coefficients_of_predictors_of_measured) * central_predictors_in_test_profiles
        @assert_matrix(predicted_measured_of_test_profiles, n_measured, n_test_profiles, Columns)
        predicted_measured_of_test_profiles .+= means_of_measured_in_train_profiles

        predicted_measured_of_test_profiles .-= measured_in_profiles[:, test_profiles_indices]
        predicted_measured_of_test_profiles .*= predicted_measured_of_test_profiles
        rmse_of_measured = sqrt.(vec(mean(predicted_measured_of_test_profiles; dims = 2)))
        @assert_vector(rmse_of_measured, n_measured)

        rmse_of_parts[part_index] = mean(rmse_of_measured)
    end

    return mean(rmse_of_parts)
end

end  # module
