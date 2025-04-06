"""
Approximate the manifold of actual cell states (captured by metacells) using linear programs in each local region.
"""
module Programs

export compute_blocks!
export compute_global_programs!
export compute_local_programs!

using ..Contracts
using ..Defaults
using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using DataAxesFormats
using Distances
using Distributions
using MultivariateStats
using NMF
using NonNegLeastSquares
using Printf
using Random
using SparseArrays
using Statistics
using TanayLabUtilities

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
    is_relaxed = true, # TODOX
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        metacell_total_UMIs_vector(RequiredInput),
        gene_divergence_vector(RequiredInput),
        #gene_is_transcription_factor_vector(RequiredInput),
        #gene_is_forbidden_factor_vector(OptionalInput),
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

    names_of_genes = axis_vector(daf, "gene")

    n_metacells = axis_length(daf, "metacell")
    @assert n_metacells > 0
    @debug "Metacells: $(n_metacells)"

    divergence_of_genes = get_vector(daf, "gene", "divergence").array

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    global_factor_genes_mask =
        get_vector(daf, "gene", "is_transcription_factor").array .&
        .!get_vector(daf, "gene", "is_forbidden_factor"; default = false).array
    global_factor_genes_indices = findall(global_factor_genes_mask)
    n_global_factor_genes = length(global_factor_genes_indices)
    @debug "Factors $(length(global_factor_genes_indices)): [ $(join(sort(names_of_genes[global_factor_genes_indices]), ", ")) ]"
    @assert 0 < n_global_factor_genes < n_genes

    global_measured_genes_mask = .!global_factor_genes_mask
    global_measured_genes_indices = findall(global_measured_genes_mask)
    n_global_measured_genes = length(global_measured_genes_indices)
    @assert n_global_measured_genes > 0

    total_UMIs_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    total_UMIs_of_metacells = get_vector(daf, "metacell", "total_UMIs").array

    log_fractions_of_genes_in_metacells = log2.(fractions_of_genes_in_metacells .+ gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    scales_of_genes = 1.0 .- divergence_of_genes

    total_UMIs_of_global_factor_genes_in_metacells = total_UMIs_of_genes_in_metacells[global_factor_genes_indices, :]
    log_decreased_fractions_of_global_factor_genes_in_metacells,
    log_increased_fractions_of_global_factor_genes_in_metacells = compute_confidence(;
        fractions_of_genes_in_metacells = fractions_of_genes_in_metacells[global_factor_genes_indices, :],
        total_UMIs_of_metacells,
        scales_of_genes = scales_of_genes[global_factor_genes_indices],
        gene_fraction_regularization,
        fold_confidence,
    )

    blocks = compute_blocks_by_confidence(;
        log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_global_factor_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_global_factor_genes_in_metacells,
        total_UMIs_of_genes_in_metacells = total_UMIs_of_global_factor_genes_in_metacells,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        max_block_span,
    )
    @debug "Blocks: $(blocks.n_blocks)"
    block_names = group_names(daf, "metacell", blocks.metacells_of_blocks; prefix = "B")

    add_axis!(daf, "block", block_names)
    set_vector!(daf, "metacell", "block", block_names[blocks.blocks_of_metacells]; overwrite)

    add_axis!(daf, "factor", names_of_genes[global_factor_genes_mask])
    set_vector!(daf, "factor", "gene", names_of_genes[global_factor_genes_mask])

    scaled_log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells .* scales_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_metacells, n_genes, n_metacells)

    scaled_log_fractions_of_global_factor_genes_in_metacells =
        scaled_log_fractions_of_genes_in_metacells[global_factor_genes_indices, :]
    @assert_matrix(scaled_log_fractions_of_global_factor_genes_in_metacells, n_global_factor_genes, n_metacells)

    scaled_log_fractions_of_global_measured_genes_in_metacells =
        scaled_log_fractions_of_genes_in_metacells[global_measured_genes_mask, :]
    @assert_matrix(scaled_log_fractions_of_global_measured_genes_in_metacells, n_global_measured_genes, n_metacells)

    scaled_log_fractions_of_genes_in_blocks =
        transpose(
            daf["/ metacell / gene : fraction % Log base 2 eps $(gene_fraction_regularization) @ block %> Mean"].array,
        ) .* scales_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_blocks, n_genes, blocks.n_blocks)

    scaled_log_fractions_of_global_factor_genes_in_blocks =
        scaled_log_fractions_of_genes_in_blocks[global_factor_genes_indices, :]
    @assert_matrix(scaled_log_fractions_of_global_factor_genes_in_blocks, n_global_factor_genes, blocks.n_blocks)

    scaled_log_fractions_of_global_measured_genes_in_blocks =
        scaled_log_fractions_of_genes_in_blocks[global_measured_genes_mask, :]
    @assert_matrix(scaled_log_fractions_of_global_measured_genes_in_blocks, n_global_measured_genes, blocks.n_blocks)

    is_in_neighborhood_of_block_of_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    is_in_environment_of_block_of_block = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    expansion_pca_components_of_blocks = Vector{UInt16}(undef, blocks.n_blocks)
    expansion_x_rmse_of_blocks = Vector{Float32}(undef, blocks.n_blocks)
    is_marker_of_genes_of_blocks = Matrix{Bool}(undef, n_genes, blocks.n_blocks)

    #TODOX @threads
    for block_index in 1:(blocks.n_blocks)
        @views is_in_neighborhood_of_block = is_in_neighborhood_of_block_of_block[:, block_index]
        @views is_in_environment_of_block = is_in_environment_of_block_of_block[:, block_index]
        @views is_marker_of_genes_of_block = is_marker_of_genes_of_blocks[:, block_index]
        expansion_pca_components_of_blocks[block_index], expansion_x_rmse_of_blocks[block_index] =
            compute_expanded_vicinity_of_block!(;
                daf,
                gene_fraction_regularization,
                scaled_log_fractions_of_global_factor_genes_in_metacells,
                scaled_log_fractions_of_global_measured_genes_in_metacells,
                distances_between_blocks = blocks.distances_between_blocks,
                min_blocks_in_neighborhood,
                min_metacells_in_neighborhood,
                min_blocks_in_approximate_environment,
                min_metacells_in_approximate_environment,
                max_principal_components,
                min_marker_gene_range_fold,
                min_marker_gene_max_fraction,
                cross_validation_parts,
                rng,
                blocks,
                block_index,
                is_in_neighborhood_of_block,
                is_in_environment_of_block,
                is_marker_of_genes_of_block,
                global_factor_genes_indices,
            )
    end

    set_vector!(daf, "block", "expansion_pca_components", expansion_pca_components_of_blocks)
    set_vector!(daf, "block", "expansion_pca_x_rmse", expansion_x_rmse_of_blocks)

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
    set_matrix!(daf, "gene", "block", "is_marker", SparseMatrixCSC(is_marker_of_genes_of_blocks))

    return nothing
end

function compute_confidence(;  # untested
    fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_metacells::AbstractVector{<:Unsigned},
    scales_of_genes::AbstractVector{<:AbstractFloat},
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}}
    n_genes, n_metacells = size(fractions_of_genes_in_metacells)
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)
    @assert_vector(total_UMIs_of_metacells, n_metacells)
    @assert_vector(scales_of_genes, n_genes)

    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_fractions_of_genes_in_metacells =  # NOJET
        confidence_stdevs .* sqrt.(transpose(total_UMIs_of_metacells) .* fractions_of_genes_in_metacells) ./
        transpose(total_UMIs_of_metacells)

    log_decreased_fractions_of_genes_in_metacells =  # NOJET
        log2.(
            max.(fractions_of_genes_in_metacells .- confidence_fractions_of_genes_in_metacells, 0.0) .+
            gene_fraction_regularization
        ) .* scales_of_genes

    log_increased_fractions_of_genes_in_metacells =
        log2.(
            fractions_of_genes_in_metacells .+ confidence_fractions_of_genes_in_metacells .+
            gene_fraction_regularization
        ) .* scales_of_genes

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
    #TODOX @threads
    for base_metacell_index in reverse(2:(n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views total_UMIs_of_genes_in_base_metacell = vec(total_UMIs_of_genes_in_metacells[:, base_metacell_index])
        @views log_decreased_fractions_of_genes_in_base_metacell =  # NOJET
            vec(log_decreased_fractions_of_genes_in_metacells[:, base_metacell_index])
        @views log_increased_fractions_of_genes_in_base_metacell =
            vec(log_increased_fractions_of_genes_in_metacells[:, base_metacell_index])  # NOJET

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
    #TODOX @threads
    for base_block_index in reverse(2:n_blocks)
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
    #TODOX @threads
    for base_block_index in 1:n_blocks
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

function compute_expanded_vicinity_of_block!(;  # untested
    daf::DafReader,
    gene_fraction_regularization::AbstractFloat,
    scaled_log_fractions_of_global_factor_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    scaled_log_fractions_of_global_measured_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    distances_between_blocks::AbstractMatrix{<:AbstractFloat},
    min_blocks_in_neighborhood::Integer,
    min_metacells_in_neighborhood::Integer,
    min_blocks_in_approximate_environment::Integer,
    min_metacells_in_approximate_environment::Integer,
    max_principal_components::Integer,
    min_marker_gene_range_fold::Real,
    min_marker_gene_max_fraction::AbstractFloat,
    cross_validation_parts::Integer,
    rng::AbstractRNG,
    blocks::Blocks,
    block_index::Integer,
    is_in_neighborhood_of_block::AbstractVector{Bool},
    is_in_environment_of_block::AbstractVector{Bool},
    is_marker_of_genes_of_block::AbstractVector{Bool},
    global_factor_genes_indices::Any,
)::Tuple{Integer, AbstractFloat}
    n_metacells = size(scaled_log_fractions_of_global_factor_genes_in_metacells, 2)
    n_global_factor_genes = length(global_factor_genes_indices)
    n_global_measured_genes = size(scaled_log_fractions_of_global_measured_genes_in_metacells, 1)

    @assert_matrix(
        scaled_log_fractions_of_global_factor_genes_in_metacells,
        n_global_factor_genes,
        n_metacells,
        Columns
    )
    @assert_matrix(
        scaled_log_fractions_of_global_measured_genes_in_metacells,
        n_global_measured_genes,
        n_metacells,
        Columns
    )

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
    scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells =
        scaled_log_fractions_of_global_factor_genes_in_metacells[:, approximate_environment_metacells_mask]
    n_approximate_environment_metacells =
        size(scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells, 2)
    @debug "- Approximate environment blocks: $(n_blocks_in_approximate_environment) metacells: $(n_approximate_environment_metacells)"

    pca = fit(
        PCA,
        scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells;
        maxoutdim = max_principal_components,
    )
    @assert indim(pca) == n_global_factor_genes
    @assert outdim(pca) <= max_principal_components
    max_n_principal_components = outdim(pca)
    @debug "- Max principal components: $(max_n_principal_components)"

    coefficients_of_global_factor_genes_in_principal_components = loadings(pca)  # NOJET
    @assert_matrix(
        coefficients_of_global_factor_genes_in_principal_components,
        n_global_factor_genes,
        max_n_principal_components,
        Columns
    )
    mean_scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells =
        vec(mean(scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells; dims = 2))
    @assert_vector(
        mean_scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells,
        n_global_factor_genes
    )
    scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells .-=
        mean_scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells

    values_of_principal_components_of_approximate_environment_metacells =
        transpose(coefficients_of_global_factor_genes_in_principal_components) *
        scaled_log_fractions_of_global_factor_genes_in_approximate_environment_metacells
    @assert_matrix(
        values_of_principal_components_of_approximate_environment_metacells,
        max_n_principal_components,
        n_approximate_environment_metacells
    )

    scaled_log_fractions_of_global_measured_genes_of_approximate_environment_metacells =
        scaled_log_fractions_of_global_measured_genes_in_metacells[:, approximate_environment_metacells_mask]

    neighborhood_metacells_mask_in_approximate_environment =
        neighborhood_metacells_mask[approximate_environment_metacells_mask]
    neighborhood_metacells_indices_in_approximate_environment =
        findall(neighborhood_metacells_mask_in_approximate_environment)

    n_principal_components, x_rmse, _ = fast_minimize_cost(;
        minimal_value = 1,
        maximal_value = max_n_principal_components,
        linear_init = true,
        prefer_smaller = true,
    ) do m_principal_components
        x_rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = values_of_principal_components_of_approximate_environment_metacells[
                1:m_principal_components,
                :,
            ],
            actual_in_profiles = scaled_log_fractions_of_global_measured_genes_of_approximate_environment_metacells,
            core_profiles_indices = neighborhood_metacells_indices_in_approximate_environment,
            cross_validation_parts,
            rng,
        )

        @debug "TODOX Predictive principal components: $(m_principal_components) X-RMSE: $(x_rmse)"
        return (x_rmse, nothing)
    end
    @debug "- Predictive principal components: $(n_principal_components) X-RMSE: $(x_rmse)"

    coefficients_of_global_factor_genes_in_principal_components =
        coefficients_of_global_factor_genes_in_principal_components[:, 1:n_principal_components]

    n_blocks_in_environment, x_rmse, _ = fast_minimize_cost(;
        minimal_value = n_blocks_in_neighborhood,
        maximal_value = blocks.n_blocks,
        linear_init = true,
        prefer_smaller = true,
    ) do m_blocks_in_environment
        environment_metacells_mask =
            collect_metacells_region_mask(blocks, ordered_block_indices[1:m_blocks_in_environment])
        scaled_log_fractions_of_global_factor_genes_in_environment_metacells =
            scaled_log_fractions_of_global_factor_genes_in_metacells[:, environment_metacells_mask]
        n_environment_metacells = size(scaled_log_fractions_of_global_factor_genes_in_environment_metacells, 2)

        values_of_principal_components_of_environment_metacells =
            transpose(coefficients_of_global_factor_genes_in_principal_components) *
            scaled_log_fractions_of_global_factor_genes_in_environment_metacells
        @assert_matrix(
            values_of_principal_components_of_environment_metacells,
            n_principal_components,
            n_environment_metacells
        )

        scaled_log_fractions_of_global_measured_genes_in_environment_metacells =
            scaled_log_fractions_of_global_measured_genes_in_metacells[:, environment_metacells_mask]

        neighborhood_metacells_mask_in_environment = neighborhood_metacells_mask[environment_metacells_mask]
        neighborhood_metacells_indices_in_environment = findall(neighborhood_metacells_mask_in_environment)

        x_rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = values_of_principal_components_of_environment_metacells,
            actual_in_profiles = scaled_log_fractions_of_global_measured_genes_in_environment_metacells,
            core_profiles_indices = neighborhood_metacells_indices_in_environment,
            cross_validation_parts,
            rng,
        )
        @debug "TODOX Refined environment blocks: $(m_blocks_in_environment) metacells: $(sum(environment_metacells_mask)) X-RMSE: $(x_rmse)"
        return (x_rmse, nothing)  # TODOX: G-RMSE
    end

    environment_metacells_mask = collect_metacells_region_mask(blocks, ordered_block_indices[1:n_blocks_in_environment])
    @debug "- Refined environment blocks: $(n_blocks_in_environment) metacells: $(sum(environment_metacells_mask)) X-RMSE: $(x_rmse)"

    is_in_neighborhood_of_block[ordered_block_indices[1:n_blocks_in_neighborhood]] .= true
    is_in_environment_of_block[ordered_block_indices[1:n_blocks_in_environment]] .= true

    environment_marker_genes_mask = marker_genes_of_environment(;
        daf,
        gene_fraction_regularization,
        min_marker_gene_range_fold,
        min_marker_gene_max_fraction,
        environment_metacells_mask,
    )
    is_marker_of_genes_of_block[:] .= environment_marker_genes_mask

    return (n_principal_components, x_rmse)
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
    prefer_smaller::Bool,
)::Tuple{Integer, AbstractFloat, Any}
    min_significant_rmse = 1e-3

    @assert min_significant_rmse > 0
    @assert 0 < minimal_value < maximal_value

    sample_cost, sample_result = compute_cost_result(minimal_value)
    sampled_values = [minimal_value]
    minimal_cost = sample_cost
    sampled_costs = [sample_cost]
    sampled_results = Any[sample_result]

    if linear_init
        linear_init_steps = 5
        @assert linear_init_steps > 1

        step_size = (maximal_value - minimal_value) / linear_init_steps
        for step_index in 1:linear_init_steps
            sample_value = Int(round(minimal_value + step_size * step_index))
            if sample_value !== sampled_values[end]
                @assert sample_value > sampled_values[end]
                push!(sampled_values, sample_value)

                sample_cost, sample_result = compute_cost_result(sample_value)
                push!(sampled_costs, sample_cost)
                minimal_cost = min(minimal_cost, sample_cost)
                push!(sampled_results, sample_result)

                if sample_cost == 0
                    break
                end
            end
        end
    else
        @assert prefer_smaller
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
               (length(sampled_values) > 4 && sampled_costs[end] - minimal_cost > min_significant_rmse * 2)
                break
            end
        end
    end

    while true
        best_index = nothing
        best_cost = nothing
        if prefer_smaller
            for (index, cost) in enumerate(sampled_costs)
                if cost <= minimal_cost + min_significant_rmse
                    best_index = index
                    best_cost = cost
                    break
                end
            end
        else
            for (index, cost) in reverse!(collect(enumerate(sampled_costs)))
                if cost <= minimal_cost + min_significant_rmse
                    best_index = index
                    best_cost = cost
                    break
                end
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

@kwdef mutable struct GlobalContext
    n_genes::Integer
    n_metacells::Integer
    n_blocks::Integer
    n_global_factor_genes::Integer
    n_global_measured_genes::Integer
    global_factor_genes_mask::Union{AbstractVector{Bool}, BitVector}
    global_factor_genes_indices::AbstractVector{<:Integer}
    global_measured_genes_mask::Union{AbstractVector{Bool}, BitVector}
    global_measured_genes_indices::AbstractVector{<:Integer}
    is_marker_of_genes_of_blocks::Matrix{Bool}
    scaled_log_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
    names_of_genes::AbstractVector{<:AbstractString}
    names_of_blocks::AbstractVector{<:AbstractString}
    names_of_global_factor_genes::AbstractVector{<:AbstractString}
    gene_fraction_regularization::AbstractFloat
    cross_validation_parts::Integer
    rng::AbstractRNG
end

@kwdef mutable struct LocalContext
    global_context::GlobalContext
    block_index::Integer
    block_name::AbstractString
    n_local_factor_genes::Integer
    n_local_measured_genes::Integer
    n_environment_metacells::Integer
    n_neighborhood_metacells::Integer
    scaled_log_fractions_of_local_factor_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    scaled_log_fractions_of_local_measured_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    central_scaled_log_fractions_of_local_measured_genes_in_environment_metacells::AbstractMatrix{<:AbstractFloat}
    mean_scaled_log_fractions_of_local_measured_genes_in_environment_metacells::AbstractVector{<:AbstractFloat}
    neighborhood_metacells_mask_in_environment::Union{AbstractVector{Bool}, BitVector}
    neighborhood_metacells_indices_in_environment::AbstractVector{<:Integer}
    environment_metacells_indices::AbstractVector{<:Integer}
    neighborhood_metacells_indices::AbstractVector{<:Integer}
    local_factor_genes_indices::AbstractVector{<:Integer}
    local_factor_genes_indices_in_global_factor_genes::AbstractVector{<:Integer}
    local_measured_genes_indices::AbstractVector{<:Integer}
    local_measured_genes_indices_in_global_measured_genes::AbstractVector{<:Integer}
end

@kwdef mutable struct LocalModel
    local_context::LocalContext
    method::AbstractString
    parameters::AbstractString
    filter::Maybe{AbstractFloat}
    modules::Bool
    names_of_programs::Maybe{AbstractVector{<:AbstractString}} = nothing
    n_programs::Maybe{Integer} = nothing
    n_positive_programs::Maybe{Integer} = nothing
    n_used_factors::Maybe{Integer} = nothing
    n_used_coefficients::Maybe{Integer} = nothing
    coefficients_of_local_factor_genes_in_programs::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    values_of_programs_in_environment_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    central_values_of_programs_in_environment_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    mean_values_of_programs_in_environment_metacells::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    coefficients_of_programs_of_local_measured_genes::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    predicted_local_measured_genes_in_neighborhood_metacells::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing
    g_rmse_of_local_measured_genes::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    g_rmse::Maybe{AbstractFloat} = nothing
    x_rmse::Maybe{AbstractFloat} = nothing
end

function evaluate_local_model!(model::LocalModel)::Nothing
    @assert_matrix(
        model.coefficients_of_local_factor_genes_in_programs,
        model.local_context.n_local_factor_genes,
        model.n_programs
    )

    if model.n_used_factors === nothing
        used_factors_mask = vec(maximum(abs.(model.coefficients_of_local_factor_genes_in_programs); dims = 2) .> 0)
        model.n_used_factors = sum(used_factors_mask)
    end

    if model.n_used_coefficients === nothing
        model.n_used_coefficients = sum(model.coefficients_of_local_factor_genes_in_programs .!= 0)
    end

    if model.values_of_programs_in_environment_metacells === nothing
        model.values_of_programs_in_environment_metacells =
            transpose(model.coefficients_of_local_factor_genes_in_programs) *
            model.local_context.scaled_log_fractions_of_local_factor_genes_in_environment_metacells
    end
    @assert_matrix(
        model.values_of_programs_in_environment_metacells,
        model.n_programs,
        model.local_context.n_environment_metacells
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
        model.local_context.n_environment_metacells
    )

    if model.coefficients_of_programs_of_local_measured_genes === nothing
        model.coefficients_of_programs_of_local_measured_genes =
            transpose(model.central_values_of_programs_in_environment_metacells) \
            transpose(model.local_context.central_scaled_log_fractions_of_local_measured_genes_in_environment_metacells)
    end
    @assert_matrix(
        model.coefficients_of_programs_of_local_measured_genes,
        model.n_programs,
        model.local_context.n_local_measured_genes
    )

    if model.predicted_local_measured_genes_in_neighborhood_metacells === nothing
        model.predicted_local_measured_genes_in_neighborhood_metacells =
            transpose(model.coefficients_of_programs_of_local_measured_genes) *
            model.central_values_of_programs_in_environment_metacells[
                :,
                model.local_context.neighborhood_metacells_indices_in_environment,
            ]
        model.predicted_local_measured_genes_in_neighborhood_metacells .+=  # NOJET
            model.local_context.mean_scaled_log_fractions_of_local_measured_genes_in_environment_metacells
    end
    @assert_matrix(
        model.predicted_local_measured_genes_in_neighborhood_metacells,
        model.local_context.n_local_measured_genes,
        model.local_context.n_neighborhood_metacells,
    )

    if model.g_rmse_of_local_measured_genes === nothing
        error_of_local_measured_genes_in_neighborhood_metacells =
            model.predicted_local_measured_genes_in_neighborhood_metacells .-
            model.local_context.scaled_log_fractions_of_local_measured_genes_in_environment_metacells[
                :,
                model.local_context.neighborhood_metacells_indices_in_environment,
            ]
        error_of_local_measured_genes_in_neighborhood_metacells .*=
            error_of_local_measured_genes_in_neighborhood_metacells
        model.g_rmse_of_local_measured_genes =
            sqrt.(vec(mean(error_of_local_measured_genes_in_neighborhood_metacells; dims = 2)))
    end
    @assert_vector(model.g_rmse_of_local_measured_genes, model.local_context.n_local_measured_genes)

    if model.g_rmse === nothing
        model.g_rmse = mean(model.g_rmse_of_local_measured_genes)
    end

    if model.x_rmse === nothing
        model.x_rmse = compute_cross_validation_of_prediction(
            solve_least_squares;
            predictors_in_profiles = model.values_of_programs_in_environment_metacells,
            actual_in_profiles = model.local_context.scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
            core_profiles_indices = model.local_context.neighborhood_metacells_indices_in_environment,
            cross_validation_parts = model.local_context.global_context.cross_validation_parts,
            rng = model.local_context.global_context.rng,
        )
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
    central_actual_in_profiles::AbstractMatrix{<:AbstractFloat},
)::AbstractMatrix{<:AbstractFloat}
    n_predictors, n_profiles = size(central_predictors_in_profiles)
    n_actual = size(central_actual_in_profiles, 1)
    @assert_matrix(central_predictors_in_profiles, n_predictors, n_profiles, Columns)
    @assert_matrix(central_actual_in_profiles, n_actual, n_profiles, Columns)
    coefficients_of_predictors_of_actual =
        transpose(central_predictors_in_profiles) \ transpose(central_actual_in_profiles)
    @assert_matrix(coefficients_of_predictors_of_actual, n_predictors, n_actual, Columns)
    return coefficients_of_predictors_of_actual
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
    actual_in_profiles::AbstractMatrix{<:AbstractFloat},
    core_profiles_indices::AbstractVector{<:Integer},
    cross_validation_parts::Integer,
    rng::AbstractRNG,
)::AbstractFloat
    n_predictors, n_profiles = size(predictors_in_profiles)
    n_actual = size(actual_in_profiles, 1)
    n_core_profiles = length(core_profiles_indices)

    @assert_matrix(predictors_in_profiles, n_predictors, n_profiles, Columns)
    @assert_matrix(actual_in_profiles, n_actual, n_profiles, Columns)
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

        central_actual_in_train_profiles, means_of_actual_in_train_profiles =
            centralize(actual_in_profiles[:, train_profiles_mask])
        @assert_matrix(central_actual_in_train_profiles, n_actual, n_train_profiles, Columns)
        @assert_vector(means_of_actual_in_train_profiles, n_actual)

        coefficients_of_predictors_of_actual =
            solve(central_predictors_in_train_profiles, central_actual_in_train_profiles)
        @assert_matrix(coefficients_of_predictors_of_actual, n_predictors, n_actual, Columns)

        central_predictors_in_test_profiles, _ =
            centralize(predictors_in_profiles[:, test_profiles_indices], means_of_predictors_in_train_profiles)
        @assert_matrix(central_predictors_in_test_profiles, n_predictors, n_test_profiles, Columns)

        predicted_actual_of_test_profiles =
            transpose(coefficients_of_predictors_of_actual) * central_predictors_in_test_profiles
        @assert_matrix(predicted_actual_of_test_profiles, n_actual, n_test_profiles, Columns)
        predicted_actual_of_test_profiles .+= means_of_actual_in_train_profiles

        predicted_actual_of_test_profiles .-= actual_in_profiles[:, test_profiles_indices]
        predicted_actual_of_test_profiles .*= predicted_actual_of_test_profiles
        rmse_of_actual = sqrt.(vec(mean(predicted_actual_of_test_profiles; dims = 2)))
        @assert_vector(rmse_of_actual, n_actual)

        rmse_of_parts[part_index] = mean(rmse_of_actual)
    end

    return mean(rmse_of_parts)
end

function solve_nmf(;
    values_in_profiles::AbstractMatrix{<:AbstractFloat},
    n_programs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    n_values = size(values_in_profiles, 1)
    result = nnmf(values_in_profiles, n_programs)  # NOJET
    coefficients_of_values_in_programs = result.W
    @assert_matrix(coefficients_of_values_in_programs, n_values, n_programs, Columns)

    scale_of_programs = vec(sum(coefficients_of_values_in_programs; dims = 1))
    @assert_vector(scale_of_programs, n_programs)
    scale_of_programs[scale_of_programs .== 0.0] .= 1.0

    coefficients_of_values_in_programs ./= transpose(scale_of_programs)
    return coefficients_of_values_in_programs
end

@logged function compute_local_programs!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    max_principal_components::Integer = 40,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert max_principal_components > 0
    @assert cross_validation_parts > 1

    global_context = read_global_context(; daf, gene_fraction_regularization, cross_validation_parts, rng)

    coefficients_of_global_factors_in_nmf_programs = Vector{Vector{Float32}}()
    nmf_programs_of_blocks = Vector{UInt16}(undef, global_context.n_blocks)
    nmf_g_rmse_of_blocks = Vector{Float32}(undef, global_context.n_blocks)
    nmf_x_rmse_of_blocks = Vector{Float32}(undef, global_context.n_blocks)

    #TODOX @threads
    for block_index in 1:(global_context.n_blocks)
        local_context =
            read_local_context(; daf, global_context, block_index, force_use_all_global_factor_genes = false)

        _, _, (nmf_model, density) = fast_minimize_cost(;
            minimal_value = 1,
            maximal_value = min(max_principal_components, local_context.n_local_factor_genes),
            linear_init = true,
            prefer_smaller = true,
        ) do m_nmf_programs
            nmf_coefficients_of_local_factor_genes_in_programs = solve_nmf(;
                values_in_profiles = local_context.rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells,
                n_programs = m_nmf_programs,
            )
            @assert_matrix(
                nmf_coefficients_of_local_factor_genes_in_programs,
                local_context.n_local_factor_genes,
                m_nmf_programs,
                Columns
            )

            nmf_model = LocalModel(;
                local_context,
                method = "local_nmf",
                parameters = "few=$(m_nmf_programs)",
                filter = nothing,
                modules = false,
                coefficients_of_local_factor_genes_in_programs = nmf_coefficients_of_local_factor_genes_in_programs,
                n_programs = m_nmf_programs,
            )
            evaluate_local_model!(nmf_model)

            density = depict_percent(
                sum(nmf_model.coefficients_of_local_factor_genes_in_programs .!= 0),
                length(nmf_model.coefficients_of_local_factor_genes_in_programs),
            )

            @debug "TODOX NMF programs: $(nmf_model.n_programs) Density: $(density) G-RMSE: $(nmf_model.g_rmse) X-RMSE: $(nmf_model.x_rmse)"
            return (nmf_model.x_rmse, (nmf_model, density))
        end

        @debug "- NMF programs: $(nmf_model.n_programs) G-RMSE: $(nmf_model.g_rmse) X-RMSE: $(nmf_model.x_rmse)"

        unfiltered_g_rmse = nmf_model.g_rmse
        unfiltered_x_rmse = nmf_model.x_rmse
        unfiltered_density = density

        threshold_k, x_rmse, (filtered_nmf_model, density) = fast_minimize_cost(;
            minimal_value = 1,
            maximal_value = 1001,
            linear_init = true,
            prefer_smaller = false,
        ) do threshold_k
            threshold = (threshold_k - 1) / 1000
            filtered_mask = abs.(nmf_model.coefficients_of_local_factor_genes_in_programs) .< threshold
            if !any(filtered_mask)
                @debug "TODOX Filtered NMF programs: $(nmf_model.n_programs) Threshold: $(threshold) Density: $(unfiltered_density) G-RMSE: $(unfiltered_g_rmse) X-RMSE: $(unfiltered_x_rmse)"
                return (unfiltered_x_rmse, (nmf_model, unfiltered_density))
            else
                filtered_coefficients_of_local_factors_in_local_nmf_programs =
                    copy_array(nmf_model.coefficients_of_local_factor_genes_in_programs)
                filtered_coefficients_of_local_factors_in_local_nmf_programs[filtered_mask] .= 0

                filtered_nmf_model = LocalModel(;
                    local_context,
                    method = "local_nmf",
                    parameters = "few=$(nmf_model.n_programs)",
                    filter = threshold,
                    modules = false,
                    coefficients_of_local_factor_genes_in_programs = filtered_coefficients_of_local_factors_in_local_nmf_programs,
                    n_programs = nmf_model.n_programs,
                )
                evaluate_local_model!(filtered_nmf_model)

                density = depict_percent(
                    sum(filtered_nmf_model.coefficients_of_local_factor_genes_in_programs .!= 0),
                    length(filtered_nmf_model.coefficients_of_local_factor_genes_in_programs),
                )

                @debug "TODOX Filtered NMF programs: $(nmf_model.n_programs) Threshold: $(threshold) Density: $(density) G-RMSE: $(filtered_nmf_model.g_rmse) X-RMSE: $(filtered_nmf_model.x_rmse)"
                return (filtered_nmf_model.x_rmse, (filtered_nmf_model, density))
            end
        end

        threshold = (threshold_k - 1) / 1000
        @debug "TODOX Filtered NMF programs: $(filtered_nmf_model.n_programs) Threshold: $(threshold) Density: $(density) G-RMSE: $(filtered_nmf_model.g_rmse) X-RMSE: $(filtered_nmf_model.x_rmse)"

        nmf_coefficients_of_global_factor_genes_in_programs =
            zeros(Float32, global_context.n_global_factor_genes, filtered_nmf_model.n_programs)
        nmf_coefficients_of_global_factor_genes_in_programs[
            local_context.local_factor_genes_indices_in_global_factor_genes,
            :,
        ] .= filtered_nmf_model.coefficients_of_local_factor_genes_in_programs

        for program_index in 1:(filtered_nmf_model.n_programs)
            global_factor_coefficients = nmf_coefficients_of_global_factor_genes_in_programs[:, program_index]
            @assert_vector(global_factor_coefficients, global_context.n_global_factor_genes)
            push!(coefficients_of_global_factors_in_nmf_programs, global_factor_coefficients)
        end

        nmf_programs_of_blocks[block_index] = filtered_nmf_model.n_programs
        nmf_g_rmse_of_blocks[block_index] = filtered_nmf_model.g_rmse
        nmf_x_rmse_of_blocks[block_index] = filtered_nmf_model.x_rmse
    end

    n_nmf_programs = length(coefficients_of_global_factors_in_nmf_programs)

    names_of_nmf_programs = Vector{AbstractString}(undef, n_nmf_programs)
    blocks_of_nmf_programs = Vector{AbstractString}(undef, n_nmf_programs)

    next_nmf_program_position = 1
    for block_index in 1:(global_context.n_blocks)
        block_name = global_context.names_of_blocks[block_index]
        for nmf_program_index in 1:nmf_programs_of_blocks[block_index]
            names_of_nmf_programs[next_nmf_program_position] = "B$(block_index)_N$(nmf_program_index)"
            blocks_of_nmf_programs[next_nmf_program_position] = block_name
            next_nmf_program_position += 1
        end
    end
    @assert next_nmf_program_position == n_nmf_programs + 1

    if overwrite && has_axis(daf, "local_nmf_program")
        delete_axis!(daf, "local_nmf_program")
    end
    add_axis!(daf, "local_nmf_program", names_of_nmf_programs)

    set_vector!(daf, "local_nmf_program", "block", blocks_of_nmf_programs; overwrite = overwrite)

    coefficients_of_global_factors_in_nmf_programs = hcat(coefficients_of_global_factors_in_nmf_programs...)
    @assert_matrix(coefficients_of_global_factors_in_nmf_programs, global_context.n_global_factor_genes, n_nmf_programs)
    set_matrix!(
        daf,
        "factor",
        "local_nmf_program",
        "coefficient",
        coefficients_of_global_factors_in_nmf_programs;
        overwrite = overwrite,
    )

    set_vector!(daf, "block", "local_nmf_programs", nmf_programs_of_blocks; overwrite = overwrite)
    set_vector!(daf, "block", "local_nmf_x_rmse", nmf_x_rmse_of_blocks; overwrite = overwrite)
    set_vector!(daf, "block", "local_nmf_g_rmse", nmf_g_rmse_of_blocks; overwrite = overwrite)

    return nothing
end

function read_global_context(;
    daf::DafReader,
    gene_fraction_regularization::AbstractFloat,
    cross_validation_parts::Integer,
    rng::AbstractRNG,
)::GlobalContext
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    n_global_factor_genes = axis_length(daf, "factor")

    names_of_global_factor_genes = axis_vector(daf, "factor")
    global_factor_genes_indices = axis_indices(daf, "gene", names_of_global_factor_genes)

    global_factor_genes_mask = zeros(Bool, n_genes)
    global_factor_genes_mask[global_factor_genes_indices] .= true

    global_measured_genes_mask = .!global_factor_genes_mask
    global_measured_genes_indices = findall(global_measured_genes_mask)
    n_global_measured_genes = length(global_measured_genes_indices)

    names_of_blocks = axis_vector(daf, "block")
    names_of_genes = axis_vector(daf, "gene")

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    log_fractions_of_genes_in_metacells = log2.(fractions_of_genes_in_metacells .+ gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    divergence_of_genes = get_vector(daf, "gene", "divergence").array
    scales_of_genes = 1.0 .- divergence_of_genes

    scaled_log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells .* scales_of_genes
    @assert_matrix(scaled_log_fractions_of_genes_in_metacells, n_genes, n_metacells)

    names_of_genes = axis_vector(daf, "gene")
    names_of_global_factor_genes = axis_vector(daf, "factor")
    @assert all(names_of_genes[global_factor_genes_indices] .== names_of_global_factor_genes)

    names_of_blocks = axis_vector(daf, "block")

    is_marker_of_genes_of_blocks = get_matrix(daf, "gene", "block", "is_marker")

    return GlobalContext(;
        n_genes,
        n_metacells,
        n_blocks,
        n_global_factor_genes,
        n_global_measured_genes,
        global_factor_genes_mask,
        global_factor_genes_indices,
        global_measured_genes_mask,
        global_measured_genes_indices,
        is_marker_of_genes_of_blocks,
        scaled_log_fractions_of_genes_in_metacells,
        names_of_genes,
        names_of_global_factor_genes,
        names_of_blocks,
        gene_fraction_regularization,
        cross_validation_parts,
        rng,
    )
end

function read_local_context(;
    daf::DafReader,
    global_context::GlobalContext,
    block_index::Integer,
    force_use_all_global_factor_genes::Bool,
)::LocalContext
    block_name = global_context.names_of_blocks[block_index]

    if !force_use_all_global_factor_genes
        @debug "Block: $(block_name)"
    end

    neighborhood_metacells_indices = daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"]
    environment_metacells_indices = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"]
    n_neighborhood_metacells = length(neighborhood_metacells_indices)
    n_environment_metacells = length(environment_metacells_indices)

    neighborhood_metacells_mask = zeros(Bool, global_context.n_metacells)
    environment_metacells_mask = zeros(Bool, global_context.n_metacells)

    neighborhood_metacells_mask[neighborhood_metacells_indices] .= true
    environment_metacells_mask[environment_metacells_indices] .= true

    if !force_use_all_global_factor_genes
        @debug "- Metacells in neighborhood: $(n_neighborhood_metacells)"
        @debug "- Metacells in environment: $(n_environment_metacells)"
    end

    environment_marker_genes_mask = global_context.is_marker_of_genes_of_blocks[:, block_index]
    n_local_marker_genes = sum(environment_marker_genes_mask)
    if !force_use_all_global_factor_genes
        @debug "- Marker genes in environment: $(n_local_marker_genes)"
    end

    if force_use_all_global_factor_genes
        local_factor_genes_mask = global_context.global_factor_genes_mask
        local_factor_genes_indices = global_context.global_factor_genes_indices
        n_local_factor_genes = global_context.n_global_factor_genes
        local_factor_genes_mask_in_global_factor_genes_mask = fill(true, n_local_factor_genes)
        local_factor_genes_indices_in_global_factor_genes = collect(1:n_local_factor_genes)
    else
        local_factor_genes_mask = global_context.global_factor_genes_mask .& environment_marker_genes_mask
        local_factor_genes_indices = findall(local_factor_genes_mask)
        n_local_factor_genes = length(local_factor_genes_indices)
        local_factor_genes_mask_in_global_factor_genes_mask =
            local_factor_genes_mask[global_context.global_factor_genes_mask]
        @assert sum(local_factor_genes_mask_in_global_factor_genes_mask) == n_local_factor_genes
        local_factor_genes_indices_in_global_factor_genes = findall(local_factor_genes_mask_in_global_factor_genes_mask)
        @debug "- Local factors genes: $(n_local_factor_genes) out of: $(global_context.n_global_factor_genes)"
    end

    local_measured_genes_mask = global_context.global_measured_genes_mask .& environment_marker_genes_mask
    local_measured_genes_indices = findall(local_measured_genes_mask)
    n_local_measured_genes = length(local_measured_genes_indices)

    local_measured_genes_mask_in_global_measured_genes_mask =
        local_measured_genes_mask[global_context.global_measured_genes_mask]
    @assert sum(local_measured_genes_mask_in_global_measured_genes_mask) == n_local_measured_genes
    local_measured_genes_indices_in_global_measured_genes =
        findall(local_measured_genes_mask_in_global_measured_genes_mask)

    if !force_use_all_global_factor_genes
        @debug "- Local measured genes in environment: $(n_local_measured_genes) out of: $(global_context.n_global_measured_genes)"
    end

    scaled_log_fractions_of_local_factor_genes_in_environment_metacells =
        global_context.scaled_log_fractions_of_genes_in_metacells[
            local_factor_genes_indices,
            environment_metacells_mask,
        ]
    @assert_matrix(
        scaled_log_fractions_of_local_factor_genes_in_environment_metacells,
        n_local_factor_genes,
        n_environment_metacells
    )

    rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells, _ =
        centralize(scaled_log_fractions_of_local_factor_genes_in_environment_metacells)
    rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells .=
        min.(rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells, 2.0)
    rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells .=
        max.(rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells, -2.0)
    rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells .+= 2.0

    #   rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells =
    #       scaled_log_fractions_of_local_factor_genes_in_environment_metacells .-
    #       log2(global_context.gene_fraction_regularization)

    scaled_log_fractions_of_local_measured_genes_in_environment_metacells =
        global_context.scaled_log_fractions_of_genes_in_metacells[local_measured_genes_mask, environment_metacells_mask]
    @assert_matrix(
        scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
        n_local_measured_genes,
        n_environment_metacells
    )

    central_scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
    mean_scaled_log_fractions_of_local_measured_genes_in_environment_metacells =
        centralize(scaled_log_fractions_of_local_measured_genes_in_environment_metacells)

    neighborhood_metacells_mask_in_environment = neighborhood_metacells_mask[environment_metacells_mask]
    neighborhood_metacells_indices_in_environment = findall(neighborhood_metacells_mask_in_environment)

    return LocalContext(;
        global_context,
        block_index,
        block_name,
        n_local_factor_genes,
        n_local_measured_genes,
        n_environment_metacells,
        n_neighborhood_metacells,
        scaled_log_fractions_of_local_factor_genes_in_environment_metacells,
        rebased_scaled_log_fractions_of_local_factor_genes_in_environment_metacells,
        scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
        central_scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
        mean_scaled_log_fractions_of_local_measured_genes_in_environment_metacells,
        neighborhood_metacells_mask_in_environment,
        neighborhood_metacells_indices_in_environment,
        environment_metacells_indices,
        neighborhood_metacells_indices,
        local_factor_genes_indices,
        local_factor_genes_indices_in_global_factor_genes,
        local_measured_genes_indices,
        local_measured_genes_indices_in_global_measured_genes,
    )
end

@logged function compute_global_programs!(
    daf::DafReader;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    global_context = read_global_context(; daf, gene_fraction_regularization, cross_validation_parts, rng)

    n_local_nmf_programs = axis_length(daf, "local_nmf_program")
    coefficients_of_global_factors_in_local_nmf_programs =
        get_matrix(daf, "factor", "local_nmf_program", "coefficient").array
    blocks_of_local_nmf_programs = get_vector(daf, "local_nmf_program", "block").array
    local_nmf_programs_in_block = get_vector(daf, "block", "local_nmf_programs").array
    max_local_nmf_programs_in_block = maximum(local_nmf_programs_in_block)

    g_rmse_of_local_nmf_programs = get_vector(daf, "block", "local_nmf_g_rmse").array
    x_rmse_of_local_nmf_programs = get_vector(daf, "block", "local_nmf_x_rmse").array
    @debug "Local NMF programs: $(n_local_nmf_programs) G-RMSE: $(mean(g_rmse_of_local_nmf_programs)) X-RMSE: $(mean(x_rmse_of_local_nmf_programs))"

    n_global_nmf_programs,
    x_rmse,
    (
        g_rmse,
        n_merged_nmf_programs,
        coefficients_of_global_factors_in_global_nmf_programs,
        global_nmf_program_indices_of_local_nmf_programs,
        density,
    ) = fast_minimize_cost(;
        minimal_value = max_local_nmf_programs_in_block,
        maximal_value = n_local_nmf_programs,
        linear_init = true,
        prefer_smaller = true,
    ) do m_global_nmf_programs
        result = kmeans(coefficients_of_global_factors_in_local_nmf_programs, m_global_nmf_programs; init = :kmpp)
        global_nmf_program_indices_of_local_nmf_programs = result.assignments
        @assert_vector(global_nmf_program_indices_of_local_nmf_programs, n_local_nmf_programs)

        coefficients_of_global_factors_in_global_nmf_programs = result.centers
        @assert_matrix(
            coefficients_of_global_factors_in_global_nmf_programs,
            global_context.n_global_factor_genes,
            m_global_nmf_programs,
            Columns
        )

        n_merged_nmf_programs = 0
        for block_index in 1:(global_context.n_blocks)
            block_name = global_context.names_of_blocks[block_index]
            global_nmf_program_indices_of_block_local_nmf_programs =
                global_nmf_program_indices_of_local_nmf_programs[blocks_of_local_nmf_programs .== block_name]

            n_merged_nmf_programs +=
                local_nmf_programs_in_block[block_index] -
                length(unique(global_nmf_program_indices_of_block_local_nmf_programs))
        end

        g_rmse_of_blocks, x_rmse_of_blocks = evaluate_global_nmf_programs(
            daf;
            global_context,
            blocks_of_local_nmf_programs,
            global_nmf_program_indices_of_local_nmf_programs,
            coefficients_of_global_factors_in_global_nmf_programs,
            gene_fraction_regularization,
            cross_validation_parts,
            rng,
        )

        g_rmse = mean(g_rmse_of_blocks)
        x_rmse = mean(x_rmse_of_blocks)

        density = depict_percent(
            sum(coefficients_of_global_factors_in_global_nmf_programs .!= 0),
            length(coefficients_of_global_factors_in_global_nmf_programs),
        )

        @debug "TODOX Global NMF programs: $(m_global_nmf_programs) Merged programs: $(n_merged_nmf_programs) Density: $(density) G-RMSE: $(g_rmse) X-RMSE: $(x_rmse)"

        return (
            x_rmse,
            (
                g_rmse,
                n_merged_nmf_programs,
                coefficients_of_global_factors_in_global_nmf_programs,
                global_nmf_program_indices_of_local_nmf_programs,
                density,
            ),
        )
    end

    @debug "Global NMF programs: $(n_global_nmf_programs) Merged programs: $(n_merged_nmf_programs) Density: $(density) G-RMSE: $(g_rmse) X-RMSE: $(x_rmse)"

    unfiltered_g_rmse = g_rmse
    unfiltered_x_rmse = x_rmse
    unfiltered_density = density

    threshold_k, x_rmse, (g_rmse_of_blocks, x_rmse_of_blocks, density) = fast_minimize_cost(;
        minimal_value = 1,
        maximal_value = 1001,
        linear_init = true,
        prefer_smaller = false,
    ) do threshold_k
        threshold = (threshold_k - 1) / 1000
        filtered_mask = abs.(coefficients_of_global_factors_in_global_nmf_programs) .< threshold
        if !any(filtered_mask)
            @debug "TODOX Global NMF programs: $(n_global_nmf_programs) Merged programs: $(n_merged_nmf_programs) Threshold: $(threshold) Density: $(unfiltered_density) G-RMSE: $(unfiltered_g_rmse) X-RMSE: $(unfiltered_x_rmse)"
            return (unfiltered_x_rmse, (nothing, nothing, nothing))
        else
            filtered_coefficients_of_global_factors_in_global_nmf_programs =
                copy_array(coefficients_of_global_factors_in_global_nmf_programs)
            filtered_coefficients_of_global_factors_in_global_nmf_programs[filtered_mask] .= 0
            g_rmse_of_blocks, x_rmse_of_blocks = evaluate_global_nmf_programs(
                daf;
                global_context,
                blocks_of_local_nmf_programs,
                global_nmf_program_indices_of_local_nmf_programs,
                coefficients_of_global_factors_in_global_nmf_programs = filtered_coefficients_of_global_factors_in_global_nmf_programs,
                gene_fraction_regularization,
                cross_validation_parts,
                rng,
            )
            density = depict_percent(
                sum(filtered_coefficients_of_global_factors_in_global_nmf_programs .!= 0),
                length(filtered_coefficients_of_global_factors_in_global_nmf_programs),
            )

            g_rmse = mean(g_rmse_of_blocks)
            x_rmse = mean(x_rmse_of_blocks)

            @debug "TODOX Global NMF programs: $(n_global_nmf_programs) Merged programs: $(n_merged_nmf_programs) Threshold: $(threshold) Density: $(density) G-RMSE: $(g_rmse) X-RMSE: $(x_rmse)"
            return (x_rmse, (g_rmse_of_blocks, x_rmse_of_blocks, density))
        end
    end

    g_rmse = mean(g_rmse_of_blocks)

    threshold = (threshold_k - 1) / 1000
    @debug "Global NMF programs: $(n_global_nmf_programs) Merged programs: $(n_merged_nmf_programs) Threshold: $(threshold) Density: $(density) G-RMSE: $(g_rmse) X-RMSE: $(x_rmse)"

    if overwrite && has_axis(daf, "global_nmf_program")
        delete_axis!(daf, "global_nmf_program")
    end

    global_nmf_programs_of_blocks = Vector{UInt16}(undef, global_context.n_blocks)
    for block_index in 1:(global_context.n_blocks)
        block_name = global_context.names_of_blocks[block_index]
        global_nmf_program_indices_of_block_local_nmf_programs =
            global_nmf_program_indices_of_local_nmf_programs[blocks_of_local_nmf_programs .== block_name]

        global_nmf_programs_of_blocks[block_index] =
            length(unique(global_nmf_program_indices_of_block_local_nmf_programs))
    end

    names_of_global_nmf_programs =
        ["G$(global_nmf_program_index)" for global_nmf_program_index in 1:n_global_nmf_programs]
    add_axis!(daf, "global_nmf_program", names_of_global_nmf_programs)
    filtered_coefficients_of_global_factors_in_global_nmf_programs =
        copy_array(coefficients_of_global_factors_in_global_nmf_programs)
    filtered_coefficients_of_global_factors_in_global_nmf_programs[abs.(
        filtered_coefficients_of_global_factors_in_global_nmf_programs
    ) .< threshold] .= 0
    set_matrix!(
        daf,
        "factor",
        "global_nmf_program",
        "coefficient",
        filtered_coefficients_of_global_factors_in_global_nmf_programs,
    )
    set_vector!(
        daf,
        "local_nmf_program",
        "global_nmf_program",
        names_of_global_nmf_programs[global_nmf_program_indices_of_local_nmf_programs],
    )
    set_vector!(daf, "block", "global_nmf_g_rmse", g_rmse_of_blocks)
    set_vector!(daf, "block", "global_nmf_x_rmse", x_rmse_of_blocks)
    set_vector!(daf, "block", "global_nmf_programs", global_nmf_programs_of_blocks)

    return nothing
end

function compute_modules_cluster_programs(;
    n_blocks::Integer,
    blocks_of_local_programs::AbstractVector{<:Integer},
    clusters_of_local_programs::AbstractVector{<:Integer},
    cluster_programs::AbstractMatrix{<:AbstractFloat},
)::AbstractMatrix{<:AbstractFloat}
    n_global_factor_genes, n_clusters = size(cluster_programs)
    overlap_of_clusters = zeros(Bool, n_clusters, n_clusters)
    for block_index in 1:n_blocks
        clustered_programs_of_block =
            sort!(unique(clusters_of_local_programs[blocks_of_local_programs .== block_index]))
        n_programs = length(clustered_programs_of_block)
        for high_position in 2:n_programs
            high_cluster = clustered_programs_of_block[high_position]
            for low_position in 1:high_position
                low_cluster = clustered_programs_of_block[low_position]
                overlap_of_clusters[high_cluster, low_cluster] = true
                overlap_of_clusters[low_cluster, high_cluster] = true
            end
        end
    end

    tmp_cluster_programs = copy_array(cluster_programs)
    modules_cluster_programs = copy_array(cluster_programs)
    for global_factor_gene_position in 1:n_global_factor_genes
        while true
            strongest_program = argmax(vec(abs.(tmp_cluster_programs[global_factor_gene_position, :])))
            strongest_coefficient = tmp_cluster_programs[global_factor_gene_position, strongest_program]
            if strongest_coefficient == 0
                break
            end
            tmp_cluster_programs[global_factor_gene_position, overlap_of_clusters[:, strongest_program]] .= 0
            tmp_cluster_programs[global_factor_gene_position, strongest_program] = 0
            modules_cluster_programs[global_factor_gene_position, overlap_of_clusters[:, strongest_program]] .= 0
        end
    end

    return modules_cluster_programs
end

function evaluate_global_nmf_programs(
    daf::DafReader;
    global_context::GlobalContext,
    blocks_of_local_nmf_programs::AbstractVector{<:AbstractString},
    global_nmf_program_indices_of_local_nmf_programs::AbstractVector{<:Integer},
    coefficients_of_global_factors_in_global_nmf_programs::AbstractMatrix{<:AbstractFloat},
    gene_fraction_regularization::AbstractFloat,
    cross_validation_parts::Integer,
    rng::AbstractRNG,
)::Tuple{Vector{Float32}, Vector{Float32}}
    g_rmse_of_blocks = Vector{Float32}(undef, global_context.n_blocks)
    x_rmse_of_blocks = Vector{Float32}(undef, global_context.n_blocks)

    # TODOX @threads
    for block_index in 1:(global_context.n_blocks)
        block_name = global_context.names_of_blocks[block_index]
        local_context = read_local_context(; daf, global_context, block_index, force_use_all_global_factor_genes = true)
        global_nmf_program_indices_of_block =
            sort!(unique(global_nmf_program_indices_of_local_nmf_programs[blocks_of_local_nmf_programs .== block_name]))
        n_block_global_nmf_programs = length(global_nmf_program_indices_of_block)
        coefficients_of_global_factor_genes_in_block_global_nmf_programs =
            coefficients_of_global_factors_in_global_nmf_programs[:, global_nmf_program_indices_of_block]
        @assert_matrix(
            coefficients_of_global_factor_genes_in_block_global_nmf_programs,
            global_context.n_global_factor_genes,
            n_block_global_nmf_programs,
            Columns,
        )
        model = LocalModel(;
            local_context,
            method = "global_nmf",
            parameters = "block=$(block_name)",
            filter = nothing,
            modules = false,
            coefficients_of_local_factor_genes_in_programs = coefficients_of_global_factor_genes_in_block_global_nmf_programs,
            n_programs = n_block_global_nmf_programs,
        )
        evaluate_local_model!(model)
        g_rmse_of_blocks[block_index] = model.g_rmse
        x_rmse_of_blocks[block_index] = model.x_rmse
    end

    return (g_rmse_of_blocks, x_rmse_of_blocks)
end

end  # module
