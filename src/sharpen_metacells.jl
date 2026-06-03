"""
Compute "better" metacells based on the local linear model approximating the manifold.
"""
module SharpenMetacells

export sharpen_metacells!

using Base.Threads
using Clustering
using DataAxesFormats
using ProgressMeter
using Random
using SparseArrays
using StatsBase
using TanayLabUtilities

using ..AnalyzeBlocks
using ..AnalyzeModules
using ..Contracts
using ..Defaults

import ..AnalyzeBlocks.compute_correlation_with_most_for_base_block!
import ..AnalyzeModules.maximal_cells_dispersion_of_modules!

import Base.Threads.maxthreadid
import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_is_strong_per_gene_per_block
import Metacells.Contracts.matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.vector_of_base_block_per_metacell
import Metacells.Contracts.vector_of_block_closest_by_pertinent_markers_per_cell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_base_outlier_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_block
import Metacells.Contracts.vector_of_n_metacells_per_block
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell

# TODOX: Get rid of unused members.
struct KmeansSizesBuffers{T <: AbstractFloat}
    kmeans_buffers::Tuple{KMeansBuffers{T}, KMeansBuffers{T}}
    max_best_assignments::Vector{Int}
    max_best_counts::Vector{Int}
    max_best_weight_per_cluster::Vector{Float64}
    max_best_centers::Matrix{T}
    max_best_is_too_small_per_cluster::BitVector
    max_best_is_too_large_per_cluster::BitVector
    max_current_assignments::Vector{Int}
    max_current_counts::Vector{Int}
    max_current_weight_per_cluster::Vector{Float64}
    max_current_centers::Matrix{T}
    max_current_is_too_small_per_cluster::BitVector
    max_current_is_too_large_per_cluster::BitVector
    max_candidate_assignments::Vector{Int}
    max_candidate_counts::Vector{Int}
    max_candidate_weight_per_cluster::Vector{Float64}
    max_candidate_centers::Matrix{T}
    max_candidate_is_too_small_per_cluster::BitVector
    max_candidate_is_too_large_per_cluster::BitVector
    max_split_point_indices::Vector{Int}
    max_split_points::Matrix{T}
    is_rejected_for_split_per_max_cluster::BitVector
    is_gene_in_module::BitVector
    total_UMIs_per_max_cluster_cell::Vector{Float32}
    normalized_factor_per_max_cluster_cell::Vector{Float32}
    max_best_cells_dispersion::Vector{Float32}
    max_current_cells_dispersion::Vector{Float32}
    max_candidate_cells_dispersion::Vector{Float32}
    cells_dispersion_per_max_module::Vector{Float32}
    mean_normalized_per_max_module::Vector{Float64}
end

function KmeansSizesBuffers{T}(;
    n_dims::Integer,
    max_k::Integer,
    n_points::Integer,
    n_genes::Integer,
    n_modules::Integer,
) where {T <: AbstractFloat}
    return KmeansSizesBuffers{T}(
        (KMeansBuffers{T}(; n_dims, max_k, n_points), KMeansBuffers{T}(; n_dims, max_k, n_points)),
        Vector{Int}(undef, n_points),         # max_best_assignments
        Vector{Int}(undef, max_k),            # max_best_counts
        Vector{Float64}(undef, max_k),        # max_best_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_best_centers
        BitVector(undef, max_k),              # max_best_is_too_small_per_cluster
        BitVector(undef, max_k),              # max_best_is_too_large_per_cluster
        Vector{Int}(undef, n_points),         # max_current_assignments
        Vector{Int}(undef, max_k),            # max_current_counts
        Vector{Float64}(undef, max_k),        # max_current_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_current_centers
        BitVector(undef, max_k),              # max_current_is_too_small_per_cluster
        BitVector(undef, max_k),              # max_current_is_too_large_per_cluster
        Vector{Int}(undef, n_points),         # max_candidate_assignments
        Vector{Int}(undef, max_k),            # max_candidate_counts
        Vector{Float64}(undef, max_k),        # max_candidate_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_candidate_centers
        BitVector(undef, max_k),              # max_candidate_is_too_small_per_cluster
        BitVector(undef, max_k),              # max_candidate_is_too_large_per_cluster
        Vector{Int}(undef, n_points),         # max_split_point_indices
        Matrix{T}(undef, n_dims, n_points),   # max_split_points
        BitVector(undef, max_k),              # is_rejected_for_split_per_max_cluster
        BitVector(undef, n_genes),            # is_gene_in_module
        Vector{Float32}(undef, n_points),     # total_UMIs_per_max_cluster_cell
        Vector{Float32}(undef, n_points),     # normalized_factor_per_max_cluster_cell
        Vector{Float32}(undef, max_k),        # max_best_cells_dispersion
        Vector{Float32}(undef, max_k),        # max_current_cells_dispersion
        Vector{Float32}(undef, max_k),        # max_candidate_cells_dispersion
        Vector{Float32}(undef, n_modules),    # cells_dispersion_per_max_module
        Vector{Float64}(undef, n_modules),    # mean_normalized_per_max_module
    )
end

@kwdef mutable struct KmeansSolution{T <: AbstractFloat}
    is_filled::Bool
    n_dims::Int
    n_points::Int
    k::Int
    n_too_tight::Int
    n_too_small::Int
    n_too_wide::Int
    n_too_large::Int
    centers::Matrix{T}
    assignments::Vector{Int}
    counts::Vector{Int}
    weight_per_cluster::Vector{Float64}
    cells_dispersion_per_cluster::Vector{Float32}
end

@kwdef struct LocalClusters
    n_clusters::Int
    block_cell_indices::AbstractVector{<:Integer}
    cluster_index_per_block_cell::AbstractVector{<:Integer}
    is_too_small_per_cluster::BitVector
end

@kwdef struct DispersionContext
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer}
    total_UMIs_per_cell::AbstractVector{<:Integer}
    block_cell_indices::AbstractVector{<:Integer}
    is_found_per_module::BitVector
    module_index_per_gene::AbstractVector{<:Integer}
    normalized_UMIs_quantile::AbstractFloat
    min_module_UMIs::Integer
    min_cells_dispersion::AbstractFloat
    max_cells_dispersion::AbstractFloat
end

# Phase 2 per-thread workspace for `compute_delta_correlation`. The global fields are shared across all walkable
# blocks (read-only). The per-walkable-block fields and the per-call scratch are private to this thread - the per-
# walkable-block fields get re-assigned for every work item before the call, and the scratch arrays are sized to the
# walkable blocks' maxima. One instance per thread is allocated up-front; the parallel loop body picks
# `workspace_per_thread[threadid()]` under `:static_greedy`.
@kwdef mutable struct DeltaCorrelationContext
    # Global baseline (shared across walkable blocks).
    grouped_per_base_block::AbstractVector{GroupedSeriesCorrelations{Float32}}
    baseline_mean_correlation_per_base_block::AbstractVector{Float32}
    group_index_per_block_per_base_block::AbstractMatrix{<:Integer}
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer}
    total_UMIs_per_cell::AbstractVector{<:Integer}
    gene_fraction_regularization::Float32
    # Per-walkable-block (assigned per work item before the `compute_delta_correlation` call).
    walkable_block_index::Int
    block_cell_indices::AbstractVector{<:Integer}
    affected_base_block_indices::AbstractVector{<:Integer}
    gene_index_per_friend_position::AbstractVector{<:Integer}
    friend_position_per_series_per_base_block::Dict{Int, Vector{Int}}
    block_cell_position_per_point_per_base_block::Dict{Int, Vector{Int}}
    # Per-thread scratch (sized to the walkable blocks' maxima).
    total_UMIs_per_max_cluster::Vector{Float64}
    UMIs_per_friend_position_per_max_cluster::Matrix{Float64}
    variable_per_block_cell_per_friend_position::Matrix{Float32}
    is_active_per_block_cell::Vector{Bool}
end

# Phase-1 result for a single block: a list of candidate K-solutions (owned, sized to each solution's `k`). Phase 1
# collects every perfect-K result encountered during its walk; if it found at least one perfect, all candidates are
# perfect and the baseline is the minimal-K perfect (the smallest-K candidate). If Phase 1 found no perfect, the list
# holds a single imperfect compromise (the best-by-penalty result seen during the walk). The candidates are sorted
# ascending by `k` at the end of Phase 1 - so `candidates[1]` is always the baseline.
struct SolutionCandidate
    k::Int
    counts::Vector{Int}
    weight_per_cluster::Vector{Float64}
    assignments::Vector{Int}
    is_perfect::Bool
end

struct SolutionCandidates
    n_points::Int
    candidates::Vector{SolutionCandidate}
end

"""
    function sharpen_metacells!(;
        sharp_daf::DafWriter,
        base_daf::DafReader,
        prefix::AbstractString = $(DEFAULT.prefix),
        min_cells_in_metacell::Integer = $(DEFAULT.min_cells_in_metacell),
        target_cell_total_UMIs_quantile::AbstractFloat = $(DEFAULT.target_cell_total_UMIs_quantile),
        min_metacell_total_UMIs::Integer = $(DEFAULT.min_metacell_total_UMIs),
        max_cells_dispersion_in_metacell::AbstractFloat = $(DEFAULT.max_cells_dispersion_in_metacell),
        min_cells_dispersion_in_metacell::AbstractFloat = $(DEFAULT.min_cells_dispersion_in_metacell),
        normalized_UMIs_quantile::AbstractFloat = $(DEFAULT.normalized_UMIs_quantile),
        min_module_UMIs::Integer = $(DEFAULT.min_module_UMIs),
        min_migration_likelihood::AbstractFloat = $(DEFAULT.min_migration_likelihood),
        kmeans_rounds::Integer = $(DEFAULT.kmeans_rounds),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given an `base_daf` metacells repository with a blocks structure and local gene modules that describe the cell state
manifold, compute a `sharp_daf` metacells repository, which hopefully more faithfully captures this manifold.

 1. We cluster using K-means all the cells in each neighborhood using the z-score of the expression of the modules of
    that neighborhood. The number of clusters is the number of base metacells in that neighborhood.
 2. We assign to each cluster the base block which is most frequent in the cluster cells.
    Cells of the base base block of the neighborhood, which belong to a cluster that is assigned to a different
    block, and which also belong to a cluster of that block in the neighborhood of that block, are migrated to that
    block, but only if the enrichment of the cells of that block in the cluster is at least `min_migration_likelihood`
    times what would be expected assuming random clustering based on the relative sizes of the base and other blocks.
 3. Having finalized the block to which each cell belongs to, we cluster all the cells in each block using K-means
    using the modules of the neighborhood of that block. We start with the expected number of metacells in that block
    (based on the mean number of cells per metacell in the base block) and adjust the number of clusters to try and
    enforce the sizes of the clusters - not more than twice that mean and no less than `min_cells_in_metacell`. A
    cluster is also considered too-large if its maximal `cells_dispersion` (across the block's modules) is above
    `max_cells_dispersion_in_metacell`, and too-small if its maximal `cells_dispersion` is below
    `min_cells_dispersion_in_metacell`. In edge cases we dissolve too-small clusters, so this can create new outlier
    cells. We increase the number of target metacells for every input metacell whose maximal `cells_dispersion` is
    above `max_cells_dispersion_in_metacell`.
 4. The final clusters are the sharp metacells. We name them using the `prefix`, the convention is to advance the letter
    for each sharpening round (`M` to `N` to `O` to ...).

Whenever we call K-means we repeat the call `kmeans_rounds` times and pick the best result.

# Sharpened Metacells

$(CONTRACT1)

# Original Metacells

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "sharp_daf",
    axes = [cell_axis(RequiredInput), metacell_axis(CreatedOutput)],
    data = [vector_of_metacell_per_cell(CreatedOutput), vector_of_base_block_per_metacell(CreatedOutput)],
) Contract(;
    name = "base_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_is_base_outlier_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_closest_by_pertinent_markers_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_cells_per_block(RequiredInput),
        vector_of_n_metacells_per_block(RequiredInput),
        matrix_of_cells_dispersion_per_metacell_per_module(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        vector_of_n_modules_per_block(RequiredInput),
        matrix_of_is_found_per_module_per_block(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block(RequiredInput),
    ],
) Contract(;
    name = "base_blocks_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(RequiredInput),
        matrix_of_is_strong_per_gene_per_block(RequiredInput),
    ],
) function sharpen_metacells!(;
    sharp_daf::DafWriter,
    base_daf::DafReader,
    base_blocks_daf::DafReader,
    prefix::AbstractString = "M",
    min_cells_in_metacell::Integer = 12,
    target_cell_total_UMIs_quantile::AbstractFloat = 0.25,
    min_metacell_total_UMIs::Integer = 40000,
    min_migration_likelihood::AbstractFloat = 1.5,
    max_cells_dispersion_in_metacell::AbstractFloat = 3.0,
    min_cells_dispersion_in_metacell::AbstractFloat = 1.0,
    normalized_UMIs_quantile::AbstractFloat = function_default(
        compute_matrix_of_cells_dispersion_per_metacell_per_module!,
        :normalized_UMIs_quantile,
    ),
    min_module_UMIs::Integer = function_default(
        compute_matrix_of_cells_dispersion_per_metacell_per_module!,
        :min_module_UMIs,
    ),
    kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert min_cells_in_metacell >= 0
    @assert 0 < target_cell_total_UMIs_quantile < 1
    @assert min_metacell_total_UMIs > 0
    @assert min_migration_likelihood > 0
    @assert 0 <= min_cells_dispersion_in_metacell <= max_cells_dispersion_in_metacell
    @assert 0 <= normalized_UMIs_quantile <= 1
    @assert min_module_UMIs >= 0
    @assert kmeans_rounds > 0

    n_cells = axis_length(base_daf, "cell")
    n_blocks = axis_length(base_daf, "block")
    n_modules = axis_length(base_daf, "module")
    n_base_metacells = axis_length(base_daf, "metacell")
    name_per_block = axis_vector(base_daf, "block")

    UMIs_per_cell_per_gene = get_matrix(base_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(base_daf, "cell", "total_UMIs").array

    mean_linear_fraction_in_neighborhood_cells_per_module_per_block =
        get_matrix(base_daf, "module", "block", "mean_linear_fraction_in_neighborhood_cells").array
    std_linear_fraction_in_neighborhood_cells_per_module_per_block =
        get_matrix(base_daf, "module", "block", "std_linear_fraction_in_neighborhood_cells").array

    is_in_neighborhood_per_other_block_per_base_block =
        get_matrix(base_daf, "block", "block", "is_in_neighborhood").array
    is_found_per_module_per_block = get_matrix(base_daf, "module", "block", "is_found").array
    module_index_per_gene_per_block = base_daf["@ gene @ block :: module ?? 0 : index"].array

    block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array
    is_base_outlier_per_cell = get_vector(base_daf, "cell", "is_base_outlier").array
    closest_block_index_per_cell = base_daf["@ cell : block.closest_by_pertinent_markers : index"].array
    # Grouped cells use their metacell's block; round outliers (lost their metacell in a previous sharpening round)
    # use the closest pertinent-markers block as a synthetic base; permanent (base) outliers stay at 0 and are skipped.
    block_index_per_cell = ifelse.(
        is_base_outlier_per_cell,
        zero(eltype(block_index_per_cell)),
        ifelse.(block_index_per_cell .> 0, block_index_per_cell, closest_block_index_per_cell),
    )
    n_cells_per_block = [Int(sum(block_index_per_cell .== block_index)) for block_index in 1:n_blocks]

    # Floor on metacell size (in cells) needed for each block to clear `min_metacell_total_UMIs`, assuming each cell
    # contributes at the `target_cell_total_UMIs_quantile` of the block's per-cell UMI distribution.
    target_min_cells_per_block = Vector{Int}(undef, n_blocks)
    parallel_loop_wo_rng(
        1:n_blocks;
        name = "target_min_cells_per_block",
        order = sortperm(n_cells_per_block; rev = true),  # NOJET
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "target_min_cells_per_block"),
    ) do block_index
        @views block_total_UMIs_per_cell = total_UMIs_per_cell[block_index_per_cell .== block_index]
        if isempty(block_total_UMIs_per_cell)
            target_min_cells_per_block[block_index] = min_cells_in_metacell
        else
            target_total_UMIs_per_cell = quantile(block_total_UMIs_per_cell, target_cell_total_UMIs_quantile)  # NOLINT
            n_cells_for_min_UMIs = Int(ceil(min_metacell_total_UMIs / target_total_UMIs_per_cell))
            target_min_cells_per_block[block_index] = max(n_cells_for_min_UMIs, min_cells_in_metacell)
        end
        return nothing
    end

    n_metacells_per_block = get_vector(base_daf, "block", "n_metacells").array
    block_index_per_metacell = base_daf["@ metacell : block : index"].array
    max_cells_dispersion_per_metacell = base_daf["@ module @ metacell :: cells_dispersion >- Max"].array
    @assert_vector(max_cells_dispersion_per_metacell, n_base_metacells)

    n_target_metacells_per_block = copy_array(n_metacells_per_block)
    n_added_metacells = 0
    for metacell_index in 1:n_base_metacells
        if max_cells_dispersion_per_metacell[metacell_index] > max_cells_dispersion_in_metacell
            block_index = block_index_per_metacell[metacell_index]
            if (n_target_metacells_per_block[block_index] + 1) * target_min_cells_per_block[block_index] * 2 <
               n_cells_per_block[block_index]
                n_added_metacells += 1
                n_target_metacells_per_block[block_index_per_metacell[metacell_index]] += 1
            end
        end
    end
    @debug "Adding metacells: $(n_added_metacells) to base: $(n_base_metacells) due to cells dispersion" _group =
        :mcs_results

    mean_metacell_cells_per_block = n_cells_per_block ./ max.(n_target_metacells_per_block, 1)

    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)
    n_neighborhood_cells_per_block = [
        sum(
            (block_index_per_cell .> 0) .& getindex.(
                Ref(is_in_neighborhood_per_other_block_per_base_block[:, block_index]),
                max.(block_index_per_cell, 1),
            ),
        ) for block_index in 1:n_blocks
    ]
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    n_genes = size(module_index_per_gene_per_block, 1)
    max_n_kmeans_points = max(max_n_neighborhood_cells, maximum(n_cells_per_block))
    max_n_kmeans_clusters =
        2 * maximum(
            max(
                Int(
                    round(
                        max(
                            n_neighborhood_cells_per_block[block_index],
                            n_cells_per_block[block_index],
                        ) / mean_metacell_cells_per_block[block_index],
                    ),
                ),
                1,
            ) for block_index in 1:n_blocks
        )
    kmeans_sizes_max_buffers_per_thread = [
        KmeansSizesBuffers{Float32}(;
            n_dims = max_n_block_modules,
            max_k = max_n_kmeans_clusters,
            n_points = max_n_kmeans_points,
            n_genes,
            n_modules,
        ) for _ in 1:maxthreadid()
    ]

    preferred_block_index_per_cell_per_block = compute_preferred_block_index_per_cell_per_block(;
        base_daf,
        kmeans_rounds,
        name_per_block,
        n_cells_per_block,
        n_neighborhood_cells_per_block,
        block_index_per_cell,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        mean_metacell_cells_per_block,
        is_in_neighborhood_per_other_block_per_base_block,
        is_found_per_module_per_block,
        module_index_per_gene_per_block,
        mean_linear_fraction_in_neighborhood_cells_per_module_per_block,
        std_linear_fraction_in_neighborhood_cells_per_module_per_block,
        min_migration_likelihood,
        kmeans_sizes_max_buffers_per_thread,
        rng,
    )

    block_index_per_cell =
        compute_preferred_block_index_of_cells(; block_index_per_cell, preferred_block_index_per_cell_per_block)

    # Post-migration blocks can have more cells than pre-migration; reallocate buffers if needed.
    post_migration_n_cells_per_block =
        [sum(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    post_max_n_kmeans_points = maximum(post_migration_n_cells_per_block)
    post_max_n_kmeans_clusters =
        2 * maximum(
            max(
                Int(
                    round(
                        post_migration_n_cells_per_block[block_index] / mean_metacell_cells_per_block[block_index],
                    ),
                ),
                1,
            ) for block_index in 1:n_blocks
        )
    if post_max_n_kmeans_points > max_n_kmeans_points || post_max_n_kmeans_clusters > max_n_kmeans_clusters
        max_n_kmeans_points = max(max_n_kmeans_points, post_max_n_kmeans_points)
        max_n_kmeans_clusters = max(max_n_kmeans_clusters, post_max_n_kmeans_clusters)
        kmeans_sizes_max_buffers_per_thread = [
            KmeansSizesBuffers{Float32}(;
                n_dims = max_n_block_modules,
                max_k = max_n_kmeans_clusters,
                n_points = max_n_kmeans_points,
                n_genes,
                n_modules,
            ) for _ in 1:maxthreadid()
        ]
    end

    @assert axis_vector(base_blocks_daf, "cell") == axis_vector(base_daf, "cell")
    @assert axis_vector(base_blocks_daf, "gene") == axis_vector(base_daf, "gene")

    base_block_index_per_cell = base_blocks_daf["@ cell : metacell ?? 0 : block : index"].array
    is_in_base_neighborhood_per_other_base_block_per_base_block =
        get_matrix(base_blocks_daf, "block", "block", "is_in_neighborhood").array
    most_correlated_gene_in_base_neighborhood_per_gene_per_base_block =
        get_matrix(base_blocks_daf, "gene", "block", "most_correlated_gene_in_neighborhood").array
    is_strong_per_gene_per_base_block = get_matrix(base_blocks_daf, "gene", "block", "is_strong").array
    gene_name_to_index = axis_dict(base_blocks_daf, "gene")
    n_base_blocks = axis_length(base_blocks_daf, "block")

    # Relevant gene set = union of strong genes and friend genes (the most-correlated-in-base-neighborhood partner)
    # across all base blocks. Per-block covariance scoring reads `log_fraction_per_relevant_gene_per_cell` at both
    # strong and friend positions for each (s, f) pair in the cell's base block.
    is_relevant_per_gene = falses(n_genes)
    is_strong_per_gene = BitVector(undef, n_genes)
    for base_block_index in 1:n_base_blocks
        @views most_correlated_gene_in_base_neighborhood_per_gene =
            most_correlated_gene_in_base_neighborhood_per_gene_per_base_block[:, base_block_index]
        is_strong_per_gene .= @view(is_strong_per_gene_per_base_block[:, base_block_index])
        @foreach_true_index is_strong_per_gene gene_index begin  # NOLINT
            name = most_correlated_gene_in_base_neighborhood_per_gene[gene_index]  # NOLINT
            if !isempty(name)
                is_relevant_per_gene[gene_index] = true                            # strong itself      # NOLINT
                is_relevant_per_gene[gene_name_to_index[name]] = true              # its friend
            end
        end
    end
    relevant_gene_position_per_gene = zeros(Int, n_genes)
    n_relevant_genes = 0
    @inbounds for gene_index in 1:n_genes
        if is_relevant_per_gene[gene_index]
            n_relevant_genes += 1
            relevant_gene_position_per_gene[gene_index] = n_relevant_genes
        end
    end
    gene_index_per_relevant_gene_position = Vector{Int}(undef, n_relevant_genes)
    @inbounds for gene_index in 1:n_genes
        relevant_position = relevant_gene_position_per_gene[gene_index]
        if relevant_position > 0
            gene_index_per_relevant_gene_position[relevant_position] = gene_index
        end
    end

    # Strong and friend genes as raw (global) gene indices - parallel vectors per base block. Consumed by
    # `build_grouped_for_base_block` (for reading `UMIs_per_cell_per_gene` columns) and by
    # `precompute_walkable_indirection` (to build each walkable block's friend-gene subspace).
    strong_gene_indices_per_base_block = Vector{Vector{Int}}(undef, n_base_blocks)
    friend_gene_indices_per_base_block = Vector{Vector{Int}}(undef, n_base_blocks)
    for base_block_index in 1:n_base_blocks
        @views most_correlated_gene_in_base_neighborhood_per_gene =
            most_correlated_gene_in_base_neighborhood_per_gene_per_base_block[:, base_block_index]
        is_strong_per_gene .= @view(is_strong_per_gene_per_base_block[:, base_block_index])
        strong_gene_indices = Int[]
        friend_gene_indices = Int[]
        @foreach_true_index is_strong_per_gene gene_index begin  # NOLINT
            name = most_correlated_gene_in_base_neighborhood_per_gene[gene_index]  # NOLINT
            if !isempty(name)
                friend_gene_index = gene_name_to_index[name]
                push!(strong_gene_indices, gene_index)                                          # NOLINT
                push!(friend_gene_indices, friend_gene_index)
            end
        end
        strong_gene_indices_per_base_block[base_block_index] = strong_gene_indices
        friend_gene_indices_per_base_block[base_block_index] = friend_gene_indices
    end

    local_clusters_per_block = compute_local_clusters(;
        base_daf,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        min_cells_in_metacell,
        min_metacell_total_UMIs,
        min_cells_dispersion_in_metacell,
        max_cells_dispersion_in_metacell,
        normalized_UMIs_quantile,
        min_module_UMIs,
        kmeans_rounds,
        is_found_per_module_per_block,
        module_index_per_gene_per_block,
        block_index_per_cell,
        mean_metacell_cells_per_block,
        mean_linear_fraction_in_neighborhood_cells_per_module_per_block,
        std_linear_fraction_in_neighborhood_cells_per_module_per_block,
        base_block_index_per_cell,
        is_in_base_neighborhood_per_other_base_block_per_base_block,
        strong_gene_indices_per_base_block,
        friend_gene_indices_per_base_block,
        relevant_gene_position_per_gene,
        gene_index_per_relevant_gene_position,
        gene_fraction_regularization,
        max_n_kmeans_clusters,
        kmeans_sizes_max_buffers_per_thread,
        rng,
    )

    cells_of_sharp_metacells, block_name_per_sharp_metacell = combine_local_clusters(;
        min_cells_in_metacell,
        local_clusters_per_block,
        name_per_block,
        n_base_metacells,
        n_cells,
    )

    name_per_sharp_metacell = group_names(axis_vector(base_daf, "cell"), cells_of_sharp_metacells; prefix)  # NOJET
    sharp_metacell_name_per_cell = fill("", n_cells)
    for (sharp_metacell_name, cells_of_sharp_metacell) in zip(name_per_sharp_metacell, cells_of_sharp_metacells)
        sharp_metacell_name_per_cell[cells_of_sharp_metacell] .= sharp_metacell_name
    end

    add_axis!(sharp_daf, "metacell", name_per_sharp_metacell; overwrite)
    set_vector!(sharp_daf, "cell", "metacell", sharp_metacell_name_per_cell; overwrite)
    set_vector!(sharp_daf, "metacell", "base_block", block_name_per_sharp_metacell; overwrite)

    return nothing
end

function compute_preferred_block_index_per_cell_per_block(;
    base_daf::DafReader,
    kmeans_rounds::Integer,
    name_per_block::AbstractVector{<:AbstractString},
    n_cells_per_block::AbstractVector{<:Integer},
    n_neighborhood_cells_per_block::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    is_in_neighborhood_per_other_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    min_migration_likelihood::AbstractFloat,
    is_found_per_module_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    mean_linear_fraction_in_neighborhood_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_linear_fraction_in_neighborhood_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    kmeans_sizes_max_buffers_per_thread::AbstractVector{<:KmeansSizesBuffers},
    rng::AbstractRNG,
)::Vector{Maybe{SparseVector{<:Integer}}}
    n_cells = length(block_index_per_cell)
    n_blocks = length(name_per_block)
    n_modules = axis_length(base_daf, "module")

    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)

    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    preferred_block_index_per_cell_per_block = Vector{Maybe{SparseVector{<:Integer}}}(undef, n_blocks)
    preferred_block_index_per_cell_per_block .= nothing

    n_genes = size(module_index_per_gene_per_block, 1)

    z_score_per_max_module_per_max_neighborhood_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_neighborhood_cells) for _ in 1:maxthreadid()]
    is_in_neighborhood_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]
    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:maxthreadid()]
    max_n_neighborhood_clusters = maximum(
        max(
            Int(round(n_neighborhood_cells_per_block[block_index] / mean_metacell_cells_per_block[block_index])),
            1,
        ) for block_index in 1:n_blocks
    )
    preferred_block_index_per_cluster_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_clusters) for _ in 1:maxthreadid()]
    block_count_scratch_per_thread = [zeros(Int, n_blocks) for _ in 1:maxthreadid()]
    preferred_block_index_per_max_neighborhood_cell_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_cells) for _ in 1:maxthreadid()]

    # Weight progress by per-block neighborhood-cell count; blocks vary in size, so block-count alone is misleading.
    progress = DebugProgress(
        sum(n_neighborhood_cells_per_block);
        group = :mcs_loops,
        desc = "preferred_block_index_per_cell_per_block",
    )

    parallel_loop_with_rng(  # NOJET
        1:n_blocks;
        rng,
        policy = :static_greedy,
        order = sortperm(n_neighborhood_cells_per_block; rev = true),
        name = "compute_preferred_block_index_per_cell_per_block",
        progress,
    ) do block_index, rng
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_cell = is_in_neighborhood_per_cell_per_thread[threadid()]
        is_in_neighborhood_per_cell .=
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1))
        indices_of_neighborhood_cells = findall(is_in_neighborhood_per_cell)
        n_neighborhood_cells = n_neighborhood_cells_per_block[block_index]
        @assert n_neighborhood_cells == length(indices_of_neighborhood_cells)
        @assert n_neighborhood_cells > 0

        is_found_per_module = is_found_per_module_per_thread[threadid()]
        is_found_per_module .= is_found_per_module_per_block[:, block_index]
        n_block_modules = n_modules_per_block[block_index]
        @assert n_block_modules == sum(is_found_per_module)

        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]

        z_score_per_max_module_per_max_neighborhood_cell =
            z_score_per_max_module_per_max_neighborhood_cell_per_thread[threadid()]
        @views z_score_per_found_module_per_neighborhood_cell =
            z_score_per_max_module_per_max_neighborhood_cell[1:n_block_modules, 1:n_neighborhood_cells]
        @views mean_linear_fraction_in_neighborhood_cells_per_module =
            mean_linear_fraction_in_neighborhood_cells_per_module_per_block[:, block_index]
        @views std_linear_fraction_in_neighborhood_cells_per_module =
            std_linear_fraction_in_neighborhood_cells_per_module_per_block[:, block_index]
        is_gene_in_module = is_gene_in_module_per_thread[threadid()]

        compute_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_neighborhood_cell,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            is_found_per_module,
            is_in_region_per_cell = is_in_neighborhood_per_cell,
            module_index_per_gene,
            is_gene_in_module,
            mean_linear_fraction_in_neighborhood_cells_per_module,
            std_linear_fraction_in_neighborhood_cells_per_module,
        )

        n_neighborhood_clusters = max(Int(round(n_neighborhood_cells / mean_metacell_cells_per_block[block_index])), 1)
        kmeans_result = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                z_score_per_found_module_per_neighborhood_cell,
                n_neighborhood_clusters;
                buffers = kmeans_sizes_max_buffers_per_thread[threadid()].kmeans_buffers,
                rounds = kmeans_rounds,
                rng,
            )
        end
        cluster_index_per_neighborhood_cell = assignments(kmeans_result)
        n_cells_per_cluster = counts(kmeans_result)

        preferred_block_index_per_cluster = preferred_block_index_per_cluster_per_thread[threadid()]
        count_per_block = block_count_scratch_per_thread[threadid()]

        preferred_block_index_per_max_neighborhood_cell =
            preferred_block_index_per_max_neighborhood_cell_per_thread[threadid()]

        preferred_block_index_per_neighborhood_cell = pick_preferred_block_index_per_neighborhood_cell(;
            block_index,
            min_migration_likelihood,
            n_cells_per_block,
            block_index_per_cell,
            indices_of_neighborhood_cells,
            cluster_index_per_neighborhood_cell,
            n_cells_per_cluster,
            preferred_block_index_per_cluster,
            count_per_block,
            preferred_block_index_per_max_neighborhood_cell,
        )

        preferred_block_index_per_cell_per_block[block_index] =
            SparseVector(n_cells, indices_of_neighborhood_cells, Vector(preferred_block_index_per_neighborhood_cell))

        # parallel_loop already ticks progress by 1 at iteration end; add the rest of this block's neighborhood-cell
        # weight so the bar reflects actual work done.
        if progress !== nothing && n_neighborhood_cells > 1
            next!(progress; step = n_neighborhood_cells - 1)  # NOJET
        end
        return nothing
    end

    return preferred_block_index_per_cell_per_block
end

function pick_preferred_block_index_per_neighborhood_cell(;
    block_index::Integer,
    min_migration_likelihood::AbstractFloat,
    n_cells_per_block::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    indices_of_neighborhood_cells::AbstractVector{<:Integer},
    cluster_index_per_neighborhood_cell::AbstractVector{<:Integer},
    n_cells_per_cluster::AbstractVector{<:Integer},
    preferred_block_index_per_cluster::AbstractVector{<:Integer},
    count_per_block::AbstractVector{<:Integer},
    preferred_block_index_per_max_neighborhood_cell::AbstractVector{<:Integer},
)::AbstractVector{<:Integer}
    n_block_cells = n_cells_per_block[block_index]
    n_clusters = length(n_cells_per_cluster)
    n_neighborhood_cells = length(indices_of_neighborhood_cells)

    @views preferred_block_index_per_cluster[1:n_clusters] .= 0

    for cluster_index in 1:n_clusters
        fill!(count_per_block, 0)
        n_block_cells_in_cluster = 0

        for i in 1:n_neighborhood_cells
            if cluster_index_per_neighborhood_cell[i] == cluster_index
                block_index = block_index_per_cell[indices_of_neighborhood_cells[i]]
                count_per_block[block_index] += 1
                if block_index == block_index
                    n_block_cells_in_cluster += 1
                end
            end
        end

        if n_block_cells_in_cluster > 0
            most_frequent_block_index = block_index
            most_frequent_count = n_block_cells_in_cluster
            for block_index in eachindex(count_per_block)
                if count_per_block[block_index] > most_frequent_count
                    most_frequent_count = count_per_block[block_index]
                    most_frequent_block_index = block_index
                end
            end

            if most_frequent_block_index == block_index
                preferred_block_index_per_cluster[cluster_index] = most_frequent_block_index
            else
                n_most_frequent_block_cells = n_cells_per_block[most_frequent_block_index]
                n_most_frequent_block_cells_in_cluster = most_frequent_count
                n_cluster_cells = n_cells_per_cluster[cluster_index]

                block_fraction_in_cluster = n_block_cells_in_cluster / n_cluster_cells
                most_frequent_block_fraction_in_cluster = n_most_frequent_block_cells_in_cluster / n_cluster_cells

                most_frequent_block_fraction_out_of_both =
                    most_frequent_block_fraction_in_cluster /
                    (block_fraction_in_cluster + most_frequent_block_fraction_in_cluster)

                neutral_block_fraction_out_of_both =
                    n_most_frequent_block_cells / (n_block_cells + n_most_frequent_block_cells)

                if most_frequent_block_fraction_out_of_both >=
                   min_migration_likelihood * neutral_block_fraction_out_of_both
                    preferred_block_index_per_cluster[cluster_index] = most_frequent_block_index
                end
            end
        end
    end

    @views preferred_block_index_per_neighborhood_cell =
        preferred_block_index_per_max_neighborhood_cell[1:n_neighborhood_cells]
    for neighborhood_cell_position in eachindex(preferred_block_index_per_neighborhood_cell)
        preferred_block_index_per_neighborhood_cell[neighborhood_cell_position] =
            preferred_block_index_per_cluster[cluster_index_per_neighborhood_cell[neighborhood_cell_position]]
    end
    return preferred_block_index_per_neighborhood_cell
end

function compute_preferred_block_index_of_cells(;
    block_index_per_cell::AbstractVector{<:Integer},
    preferred_block_index_per_cell_per_block::Vector{Maybe{SparseVector{<:Integer}}},
)::Vector{<:Integer}
    n_cells = length(block_index_per_cell)

    n_stationary = Atomic{Int}(0)
    n_restless = Atomic{Int}(0)
    n_migrated = Atomic{Int}(0)

    new_block_index_per_cell = zeros(UInt32, n_cells)
    parallel_loop_wo_rng(
        1:n_cells;
        name = "compute_preferred_block_index_of_cells",
        progress = DebugProgress(n_cells; group = :mcs_loops, desc = "preferred_block_index_of_cells"),
        progress_chunk = 100,
    ) do cell_index
        base_block_index_of_cell = block_index_per_cell[cell_index]
        if base_block_index_of_cell == 0
            return nothing
        end

        preferred_block_index_per_block = preferred_block_index_per_cell_per_block[base_block_index_of_cell]
        if preferred_block_index_per_block === nothing
            new_block_index_per_cell[cell_index] = base_block_index_of_cell
            atomic_add!(n_stationary, 1)
            return nothing
        end

        preferred_block_index_of_cell = preferred_block_index_per_block[cell_index]
        if preferred_block_index_of_cell == 0
            preferred_block_index_per_other_block = nothing
        else
            preferred_block_index_per_other_block =
                preferred_block_index_per_cell_per_block[preferred_block_index_of_cell]
        end

        if preferred_block_index_per_other_block === nothing
            new_block_index_per_cell[cell_index] = base_block_index_of_cell
            atomic_add!(n_stationary, 1)
            return nothing
        end

        back_preferred_block_index_of_cell = preferred_block_index_per_other_block[cell_index]
        if preferred_block_index_of_cell == base_block_index_of_cell
            new_block_index_per_cell[cell_index] = base_block_index_of_cell
            atomic_add!(n_stationary, 1)
        elseif back_preferred_block_index_of_cell != preferred_block_index_of_cell
            new_block_index_per_cell[cell_index] = base_block_index_of_cell
            atomic_add!(n_restless, 1)
        else
            new_block_index_per_cell[cell_index] = preferred_block_index_of_cell
            atomic_add!(n_migrated, 1)
        end
        return nothing
    end

    @debug (
        "Cells: $(n_cells)" *
        " Stationary: $(n_stationary[]) ($(percent(n_stationary[], n_cells)))" *
        " Restless: $(n_restless[]) ($(percent(n_restless[], n_cells)))" *
        " Migrated: $(n_migrated[]) ($(percent(n_migrated[], n_cells)))"
    ) _group = :mcs_results

    return new_block_index_per_cell
end

function compute_local_clusters(;
    base_daf::DafReader,
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    min_cells_in_metacell::Integer,
    min_metacell_total_UMIs::Integer,
    min_cells_dispersion_in_metacell::AbstractFloat,
    max_cells_dispersion_in_metacell::AbstractFloat,
    normalized_UMIs_quantile::AbstractFloat,
    min_module_UMIs::Integer,
    kmeans_rounds::Integer,
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    mean_linear_fraction_in_neighborhood_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_linear_fraction_in_neighborhood_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    base_block_index_per_cell::AbstractVector{<:Integer},
    is_in_base_neighborhood_per_other_base_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    strong_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    friend_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    relevant_gene_position_per_gene::AbstractVector{<:Integer},
    gene_index_per_relevant_gene_position::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat,
    max_n_kmeans_clusters::Integer,
    kmeans_sizes_max_buffers_per_thread::AbstractVector{<:KmeansSizesBuffers},
    rng::AbstractRNG,
    max_optimization_passes::Integer = 2,
)::AbstractVector{Maybe{LocalClusters}}
    n_cells = length(total_UMIs_per_cell)
    n_blocks = size(is_found_per_module_per_block, 2)
    n_modules = axis_length(base_daf, "module")
    n_genes = size(module_index_per_gene_per_block, 1)
    n_base_blocks = size(is_in_base_neighborhood_per_other_base_block_per_base_block, 1)

    local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)
    local_clusters_per_block .= nothing
    solution_candidates_per_block = Vector{Maybe{SolutionCandidates}}(undef, n_blocks)
    solution_candidates_per_block .= nothing

    n_cells_per_block = [sum(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    max_n_block_cells = maximum(n_cells_per_block)

    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)

    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]
    is_in_block_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:maxthreadid()]
    z_score_per_max_module_per_max_block_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_block_cells) for _ in 1:maxthreadid()]

    # Per-block setup: build the z-score matrix, dispersion context, and a few derived per-block values shared by
    # Phase 1 and Phase 2 worker bodies. Captures per-thread scratch from the enclosing scope. Returns
    # `(z_score_per_found_module_per_block_cell, dispersion_context, sizes_buffers, n_block_clusters,
    # weight_per_block_cell, block_cell_indices)` or `nothing` if the block has no cells.
    function setup_block_context(block_index::Integer)
        is_in_block_per_cell = is_in_block_per_cell_per_thread[threadid()]
        is_in_block_per_cell .= block_index_per_cell .== block_index
        block_cell_indices = findall(is_in_block_per_cell)
        n_block_cells = n_cells_per_block[block_index]
        @assert n_block_cells == length(block_cell_indices)
        if n_block_cells == 0
            return nothing
        end

        is_found_per_module = is_found_per_module_per_thread[threadid()]
        is_found_per_module .= is_found_per_module_per_block[:, block_index]
        n_block_modules = n_modules_per_block[block_index]
        @assert n_block_modules == sum(is_found_per_module)

        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        z_score_per_max_module_per_max_block_cell = z_score_per_max_module_per_max_block_cell_per_thread[threadid()]
        @views z_score_per_found_module_per_block_cell =
            z_score_per_max_module_per_max_block_cell[1:n_block_modules, 1:n_block_cells]
        @views mean_linear_fraction_in_neighborhood_cells_per_module =
            mean_linear_fraction_in_neighborhood_cells_per_module_per_block[:, block_index]
        @views std_linear_fraction_in_neighborhood_cells_per_module =
            std_linear_fraction_in_neighborhood_cells_per_module_per_block[:, block_index]

        compute_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_block_cell,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            is_found_per_module,
            is_in_region_per_cell = is_in_block_per_cell,
            module_index_per_gene,
            is_gene_in_module = is_gene_in_module_per_thread[threadid()],
            mean_linear_fraction_in_neighborhood_cells_per_module,
            std_linear_fraction_in_neighborhood_cells_per_module,
        )

        n_block_clusters = max(Int(round(n_block_cells / mean_metacell_cells_per_block[block_index])), 1)
        sizes_buffers = kmeans_sizes_max_buffers_per_thread[threadid()]
        @views weight_per_block_cell = total_UMIs_per_cell[block_cell_indices]

        dispersion_context = DispersionContext(;
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            block_cell_indices,
            is_found_per_module,
            module_index_per_gene,
            normalized_UMIs_quantile,
            min_module_UMIs,
            min_cells_dispersion = min_cells_dispersion_in_metacell,
            max_cells_dispersion = max_cells_dispersion_in_metacell,
        )

        return (
            z_score_per_found_module_per_block_cell,
            dispersion_context,
            sizes_buffers,
            n_block_clusters,
            weight_per_block_cell,
            block_cell_indices,
        )
    end

    # ===== Phase 1: per-block walk to a perfect K (or exhaustion) =====
    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        policy = :static_greedy,
        order = sortperm(n_cells_per_block; rev = true),  # NOJET
        name = "compute_local_clusters_phase_1",
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "local_clusters_phase_1"),
    ) do block_index, rng
        setup = setup_block_context(block_index)
        if setup === nothing
            return nothing
        end
        (z_score_per_found_module_per_block_cell, dispersion_context, sizes_buffers, n_block_clusters,
         weight_per_block_cell, block_cell_indices) = setup
        n_block_cells = length(block_cell_indices)

        if n_block_clusters == 1
            weight_total = sum(weight_per_block_cell)
            local_clusters_per_block[block_index] = LocalClusters(;
                n_clusters = 1,
                block_cell_indices,
                cluster_index_per_block_cell = fill(1, n_block_cells),
                is_too_small_per_cluster = BitVector([
                    n_block_cells < min_cells_in_metacell || weight_total < min_metacell_total_UMIs,
                ]),
            )
            return nothing
        end

        candidates = flame_timed("kmeans_with_sizes_phase_1") do
            return kmeans_with_sizes_phase_1(
                z_score_per_found_module_per_block_cell,
                weight_per_block_cell,
                n_block_clusters;
                max_k = 2 * n_block_clusters,
                min_cluster_size = min_cells_in_metacell,
                min_cluster_weight = min_metacell_total_UMIs,
                max_cluster_size = mean_metacell_cells_per_block[block_index] * 2,
                kmeans_rounds,
                sizes_buffers,
                dispersion_context,
                rng,
            )
        end
        solution_candidates_per_block[block_index] = candidates
        # Tentative LocalClusters from the baseline (candidates[1]); Phase 2 may overwrite for walkable blocks.
        local_clusters_per_block[block_index] = build_local_clusters_from_candidate(
            candidates.candidates[1],
            block_cell_indices,
            min_cells_in_metacell,
            min_metacell_total_UMIs,
        )
        return nothing
    end

    # `block_cell_indices_per_block` is constant across optimization passes (block membership is fixed by Phase 1;
    # only cluster assignments within a block change). Same for `n_points_per_base_block` (depends only on
    # `base_block_index_per_cell` and the neighborhood matrix).
    block_cell_indices_per_block = [
        local_clusters === nothing ? Int[] : Vector{Int}(local_clusters.block_cell_indices)
        for local_clusters in local_clusters_per_block
    ]

    n_cells_per_cell_base_block = zeros(Int, n_base_blocks)
    for cell_index in 1:n_cells
        base_block_of_cell = base_block_index_per_cell[cell_index]
        if base_block_of_cell > 0
            n_cells_per_cell_base_block[base_block_of_cell] += 1
        end
    end
    n_points_per_base_block =
        vec(is_in_base_neighborhood_per_other_base_block_per_base_block' * n_cells_per_cell_base_block)
    order_per_base_block = sortperm(n_points_per_base_block; rev = true)

    # Pre-allocate the grouped-correlations outputs; each pass overwrites them via `build_grouped_for_base_block`.
    grouped_per_base_block = Vector{GroupedSeriesCorrelations{Float32}}(undef, n_base_blocks)
    baseline_mean_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
    group_index_per_block_per_base_block = zeros(Int, n_blocks, n_base_blocks)

    # ===== Walkable blocks: those whose Phase 1 produced more than one candidate K =====
    walkable_block_indices = Int[]
    for block_index in 1:n_blocks
        solution_candidates = solution_candidates_per_block[block_index]
        if solution_candidates !== nothing && length(solution_candidates.candidates) > 1
            push!(walkable_block_indices, block_index)
        end
    end

    if isempty(walkable_block_indices)
        return local_clusters_per_block
    end

    (affected_base_block_indices_per_walkable_block,
     gene_index_per_friend_position_per_walkable_block,
     _,
     friend_position_per_series_per_base_block_per_walkable_block,
     block_cell_position_per_point_per_base_block_per_walkable_block) =
        precompute_walkable_indirection(
            walkable_block_indices,
            block_cell_indices_per_block,
            base_block_index_per_cell,
            is_in_base_neighborhood_per_other_base_block_per_base_block,
            friend_gene_indices_per_base_block,
            relevant_gene_position_per_gene,
            n_base_blocks,
            n_genes,
        )

    n_walkable_blocks = length(walkable_block_indices)

    # Flatten (walkable_position, candidate_index) into a single work list covering *all* candidates - whichever index
    # is the current baseline for a walkable block is short-circuited inside the Phase 2 body to `delta = 0`. Going
    # through all candidates lets the optimization-pass loop change which candidate is the baseline across passes
    # without rebuilding the list.
    work_item_start_per_walkable_block = Vector{Int}(undef, n_walkable_blocks + 1)
    work_item_start_per_walkable_block[1] = 1
    walkable_position_per_work_item = Int[]
    candidate_index_per_work_item = Int[]
    n_block_cells_per_work_item = Int[]
    max_n_friend_positions = 0
    max_n_walkable_block_cells = 0
    for walkable_position in 1:n_walkable_blocks
        block_index = walkable_block_indices[walkable_position]
        solution_candidates = solution_candidates_per_block[block_index]
        n_candidates = length(solution_candidates.candidates)
        n_block_cells = solution_candidates.n_points
        max_n_walkable_block_cells = max(max_n_walkable_block_cells, n_block_cells)
        max_n_friend_positions = max(
            max_n_friend_positions,
            length(gene_index_per_friend_position_per_walkable_block[walkable_position]),
        )
        for candidate_index in 1:n_candidates
            push!(walkable_position_per_work_item, walkable_position)
            push!(candidate_index_per_work_item, candidate_index)
            push!(n_block_cells_per_work_item, n_block_cells)
        end
        work_item_start_per_walkable_block[walkable_position + 1] =
            work_item_start_per_walkable_block[walkable_position] + n_candidates
    end
    n_work_items = length(walkable_position_per_work_item)
    order_per_work_item = sortperm(n_block_cells_per_work_item; rev = true)  # NOJET

    # One mutable `DeltaCorrelationContext` per thread - the per-walkable-block fields get re-assigned per work item
    # under top-level `:static_greedy` (stable `threadid()`), and the scratch fields stay sized to the global maxima.
    delta_context_per_thread = [
        DeltaCorrelationContext(;
            grouped_per_base_block,
            baseline_mean_correlation_per_base_block,
            group_index_per_block_per_base_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            gene_fraction_regularization = Float32(gene_fraction_regularization),
            # Per-walkable-block fields - re-assigned per work item.
            walkable_block_index = 0,
            block_cell_indices = Int[],
            affected_base_block_indices = Int[],
            gene_index_per_friend_position = Int[],
            friend_position_per_series_per_base_block = Dict{Int, Vector{Int}}(),
            block_cell_position_per_point_per_base_block = Dict{Int, Vector{Int}}(),
            # Per-thread scratch.
            total_UMIs_per_max_cluster = Vector{Float64}(undef, max_n_kmeans_clusters),
            UMIs_per_friend_position_per_max_cluster =
                Matrix{Float64}(undef, max_n_friend_positions, max_n_kmeans_clusters),
            variable_per_block_cell_per_friend_position =
                Matrix{Float32}(undef, max_n_walkable_block_cells, max_n_friend_positions),
            is_active_per_block_cell = Vector{Bool}(undef, max_n_walkable_block_cells),
        ) for _ in 1:maxthreadid()
    ]
    delta_per_work_item = Vector{Float64}(undef, n_work_items)

    # Per walkable block, the candidate index currently treated as the baseline. Initialized to `1` (the minimal-K
    # candidate, which is also what `local_clusters_per_block` was built from after Phase 1). Each optimization pass
    # updates this to whatever candidate it picks.
    chosen_candidate_index_per_walkable_block = fill(1, n_walkable_blocks)

    # Per-walkable-block scratches sized to that walkable block's own dimensions. Pass A (parallel argmax) populates
    # these for every changed walkable block; Pass B (parallel `replace_group!`) reads them per (changed_walkable,
    # affected base block) without re-populating - the two-pass split eliminates the populate redundancy that a
    # parallel-by-base-block commit would otherwise pay.
    variable_per_block_cell_per_friend_position_per_walkable_block = Vector{Matrix{Float32}}(undef, n_walkable_blocks)
    is_active_per_block_cell_per_walkable_block = Vector{Vector{Bool}}(undef, n_walkable_blocks)
    for walkable_position in 1:n_walkable_blocks
        block_index = walkable_block_indices[walkable_position]
        n_block_cells = solution_candidates_per_block[block_index].n_points
        n_friend_positions = length(gene_index_per_friend_position_per_walkable_block[walkable_position])
        variable_per_block_cell_per_friend_position_per_walkable_block[walkable_position] =
            Matrix{Float32}(undef, n_block_cells, n_friend_positions)
        is_active_per_block_cell_per_walkable_block[walkable_position] = Vector{Bool}(undef, n_block_cells)
    end

    # Per-base-block commit queues + locks for Pass A's append. Cleared per pass via `empty!` (no realloc), so this
    # whole structure is allocated once and reused across all optimization passes.
    commits_per_base_block = [Int[] for _ in 1:n_base_blocks]
    append_lock_per_base_block = [SpinLock() for _ in 1:n_base_blocks]

    # Run one optimization pass: evaluate every (walkable, candidate) work item's delta-correlation against the cached
    # groupeds (Phase 2), then in parallel per walkable argmax the deltas + populate the chosen candidate's persistent
    # scratch + queue commits (Pass A), then in parallel per base block drain queued commits via `replace_group!` and
    # refresh the cached baseline mean correlation (Pass B). Returns the number of walkable blocks whose chosen
    # candidate changed.
    function run_optimization_pass!(pass_index::Integer)::Int
        # Phase 2 flat parallel loop over all (walkable block, candidate) pairs. The current baseline candidate of
        # each walkable block short-circuits to `delta = 0` (it's the per-pass identity).
        progress_phase_2 = DebugProgress(
            sum(n_block_cells_per_work_item);
            group = :mcs_loops,
            desc = "local_clusters_phase_2 pass $pass_index",
        )
        parallel_loop_wo_rng(
            1:n_work_items;
            policy = :static_greedy,
            order = copy(order_per_work_item),
            name = "compute_local_clusters_phase_2.pass_$pass_index",
        ) do work_index
            walkable_position = walkable_position_per_work_item[work_index]
            candidate_index = candidate_index_per_work_item[work_index]
            if candidate_index == chosen_candidate_index_per_walkable_block[walkable_position]
                delta_per_work_item[work_index] = 0.0
            else
                block_index = walkable_block_indices[walkable_position]
                solution_candidates = solution_candidates_per_block[block_index]
                candidate = solution_candidates.candidates[candidate_index]
                delta_context = delta_context_per_thread[threadid()]
                delta_context.walkable_block_index = block_index
                delta_context.block_cell_indices = block_cell_indices_per_block[block_index]
                delta_context.affected_base_block_indices =
                    affected_base_block_indices_per_walkable_block[walkable_position]
                delta_context.gene_index_per_friend_position =
                    gene_index_per_friend_position_per_walkable_block[walkable_position]
                delta_context.friend_position_per_series_per_base_block =
                    friend_position_per_series_per_base_block_per_walkable_block[walkable_position]
                delta_context.block_cell_position_per_point_per_base_block =
                    block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position]
                delta_per_work_item[work_index] =
                    compute_delta_correlation(candidate, solution_candidates.n_points, delta_context)
            end
            if progress_phase_2 !== nothing
                next!(progress_phase_2; step = n_block_cells_per_work_item[work_index])  # NOJET
            end
            return nothing
        end

        # Clear last pass's per-base-block commit queues (preserves capacity, no realloc).
        for base_block_index in 1:n_base_blocks
            empty!(commits_per_base_block[base_block_index])
        end

        # Pass A: parallel over walkable blocks. Each task argmaxes deltas over its candidates; on change, populates
        # the walkable's persistent variable + mask scratch via `populate_candidate_scratches!`, rebuilds
        # `local_clusters_per_block[block_index]`, and lock-appends the walkable to each affected base block's commit
        # queue. The populate runs at most once per changed walkable across the whole pass.
        n_changes = Threads.Atomic{Int}(0)
        parallel_loop_wo_rng(
            1:n_walkable_blocks;
            policy = :static_greedy,
            name = "compute_local_clusters_phase_2_reduce.pass_$pass_index",
            progress = DebugProgress(
                n_walkable_blocks; group = :mcs_loops, desc = "phase_2_reduce pass $pass_index",
            ),
        ) do walkable_position
            block_index = walkable_block_indices[walkable_position]
            solution_candidates = solution_candidates_per_block[block_index]
            candidates = solution_candidates.candidates
            previous_candidate_index = chosen_candidate_index_per_walkable_block[walkable_position]
            best_candidate_index = previous_candidate_index
            best_delta = 0.0
            slice_start = work_item_start_per_walkable_block[walkable_position]
            slice_end = work_item_start_per_walkable_block[walkable_position + 1] - 1
            for work_index in slice_start:slice_end
                delta = delta_per_work_item[work_index]
                if delta > best_delta
                    best_delta = delta
                    best_candidate_index = candidate_index_per_work_item[work_index]
                end
            end
            if best_candidate_index == previous_candidate_index
                return nothing
            end
            chosen = candidates[best_candidate_index]
            @warn "TODOX - PHASE 2 PASS $pass_index CHOSEN K: $(chosen.k) (baseline K: $(candidates[previous_candidate_index].k); delta: $(round(best_delta; digits=4)); $(length(candidates)) candidates)"
            Threads.atomic_add!(n_changes, 1)
            chosen_candidate_index_per_walkable_block[walkable_position] = best_candidate_index
            local_clusters_per_block[block_index] = build_local_clusters_from_candidate(
                chosen,
                block_cell_indices_per_block[block_index],
                min_cells_in_metacell,
                min_metacell_total_UMIs,
            )

            # Populate this walkable's persistent variable + mask scratch under the chosen candidate.
            delta_context = delta_context_per_thread[threadid()]
            delta_context.walkable_block_index = block_index
            delta_context.block_cell_indices = block_cell_indices_per_block[block_index]
            delta_context.affected_base_block_indices = affected_base_block_indices_per_walkable_block[walkable_position]
            delta_context.gene_index_per_friend_position =
                gene_index_per_friend_position_per_walkable_block[walkable_position]
            delta_context.friend_position_per_series_per_base_block =
                friend_position_per_series_per_base_block_per_walkable_block[walkable_position]
            delta_context.block_cell_position_per_point_per_base_block =
                block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position]
            populate_candidate_scratches!(
                chosen,
                solution_candidates.n_points,
                delta_context,
                variable_per_block_cell_per_friend_position_per_walkable_block[walkable_position],
                is_active_per_block_cell_per_walkable_block[walkable_position],
            )

            # Lock-append the walkable to each affected base block's commit queue.
            for base_block_index in affected_base_block_indices_per_walkable_block[walkable_position]
                if group_index_per_block_per_base_block[block_index, base_block_index] == 0
                    continue
                end
                lock(append_lock_per_base_block[base_block_index]) do
                    return push!(commits_per_base_block[base_block_index], walkable_position)
                end
            end
            return nothing
        end

        # Pass B: parallel over base blocks. Each task drains its commit queue with `replace_group!` reusing the
        # per-walkable populated scratch (no re-populate), then refreshes the cached baseline mean correlation. Base
        # blocks with no commits skip the recompute - their mean is unchanged.
        parallel_loop_wo_rng(
            1:n_base_blocks;
            policy = :static_greedy,
            name = "compute_local_clusters_phase_2_commit.pass_$pass_index",
            progress = DebugProgress(
                n_base_blocks; group = :mcs_loops, desc = "phase_2_commit pass $pass_index",
            ),
        ) do base_block_index
            commits = commits_per_base_block[base_block_index]
            if isempty(commits)
                return nothing
            end
            grouped_for_base_block = grouped_per_base_block[base_block_index]
            for walkable_position in commits
                block_index = walkable_block_indices[walkable_position]
                group_index_in_base_block = group_index_per_block_per_base_block[block_index, base_block_index]
                block_cell_position_per_point =
                    block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position][base_block_index]
                friend_position_per_series =
                    friend_position_per_series_per_base_block_per_walkable_block[walkable_position][base_block_index]
                replace_group!(
                    grouped_for_base_block,
                    group_index_in_base_block,
                    variable_per_block_cell_per_friend_position_per_walkable_block[walkable_position],
                    is_active_per_block_cell_per_walkable_block[walkable_position],
                    block_cell_position_per_point,
                    friend_position_per_series,
                )
            end
            baseline_mean_correlation_per_base_block[base_block_index] =
                Float32(mean_correlation(grouped_for_base_block))
            return nothing
        end

        return n_changes[]
    end

    # One-time setup: build the per-base-block grouped correlations from the Phase-1 `local_clusters_per_block`. All
    # later updates happen in place via Pass B's `replace_group!`.
    (total_UMIs_per_baseline_metacell, UMIs_per_relevant_gene_per_baseline_metacell,
     baseline_metacell_index_per_cell, _) =
        compute_baseline_metacell_aggregates(
            local_clusters_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            gene_index_per_relevant_gene_position,
            n_cells,
        )

    parallel_loop_wo_rng(
        1:n_base_blocks;
        policy = :static_greedy,
        order = copy(order_per_base_block),
        name = "compute_local_clusters.build_grouped",
        progress = DebugProgress(n_base_blocks; group = :mcs_loops, desc = "build_grouped"),
    ) do base_block_index
        @views is_in_base_neighborhood_for_base_block =
            is_in_base_neighborhood_per_other_base_block_per_base_block[:, base_block_index]
        (grouped, group_index_per_block) = build_grouped_for_base_block(;
            n_blocks,
            n_cells,
            base_block_index_per_cell,
            block_index_per_cell,
            block_cell_indices_per_block,
            is_in_base_neighborhood_for_base_block,
            baseline_metacell_index_per_cell,
            total_UMIs_per_baseline_metacell,
            UMIs_per_relevant_gene_per_baseline_metacell,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            strong_gene_indices = strong_gene_indices_per_base_block[base_block_index],
            friend_gene_indices = friend_gene_indices_per_base_block[base_block_index],
            relevant_gene_position_per_gene,
            gene_fraction_regularization,
        )
        grouped_per_base_block[base_block_index] = grouped
        @views group_index_per_block_per_base_block[:, base_block_index] .= group_index_per_block
        baseline_mean_correlation_per_base_block[base_block_index] = Float32(mean_correlation(grouped))
        return nothing
    end

    # Run up to `max_optimization_passes` passes; skip remaining passes once one converges (no walkable block changed).
    for pass_index in 1:max_optimization_passes
        n_changes = run_optimization_pass!(pass_index)
        @warn "TODOX - PHASE 2 PASS $pass_index: $n_changes walkable blocks moved"
        if n_changes == 0
            break
        end
    end

    return local_clusters_per_block
end

function combine_local_clusters(;
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    name_per_block::AbstractVector{<:AbstractString},
    min_cells_in_metacell::Integer,
    n_base_metacells::Integer,
    n_cells::Integer,
)::Tuple{AbstractVector{<:AbstractVector{<:Integer}}, <:AbstractVector{<:AbstractString}}
    n_total_new_metacells = 0
    n_new_outlier_cells = 0
    for local_clusters in local_clusters_per_block
        if local_clusters !== nothing
            n_total_new_metacells += local_clusters.n_clusters
            @foreach_true_index local_clusters.is_too_small_per_cluster cluster_index begin  # NOLINT
                n_new_outlier_cells += sum(local_clusters.cluster_index_per_block_cell .== cluster_index)  # NOLINT
            end
        end
    end

    cells_of_new_metacells = Vector{AbstractVector{<:Integer}}()
    block_name_per_new_metacell = Vector{AbstractString}()

    for (block_index, local_clusters) in enumerate(local_clusters_per_block)
        if local_clusters !== nothing
            block_name = name_per_block[block_index]
            for cluster_index in 1:local_clusters.n_clusters
                if !local_clusters.is_too_small_per_cluster[cluster_index]
                    @views cell_indices_of_new_metacell =
                        local_clusters.block_cell_indices[local_clusters.cluster_index_per_block_cell .== cluster_index]
                    @assert length(cell_indices_of_new_metacell) >= min_cells_in_metacell
                    push!(cells_of_new_metacells, cell_indices_of_new_metacell)
                    push!(block_name_per_new_metacell, block_name)
                end
            end
        end
    end

    @debug (
        "Metacells Original: $(n_base_metacells)" *
        " Sharpened: $(length(block_name_per_new_metacell))" *
        " Outlier cells: $(n_new_outlier_cells)" *
        " ($(percent(n_new_outlier_cells, n_cells)))"
    ) _group = :mcs_results

    return (cells_of_new_metacells, block_name_per_new_metacell)
end

# Build one base block's baseline `GroupedSeriesCorrelations{Float32}`. Each group corresponds to one block whose
# cells contribute to the base block's neighborhood; the points within a group are the contributing cells laid out
# contiguously, preserving the per-block cell order. Per (point, base-block series) we store the cell's
# `log2(strong_UMIs / total_UMIs + reg)` as the fixed side and the punctuated baseline-metacell
# `log2((baseline_metacell_friend_UMIs - cell_friend_UMIs) / (baseline_metacell_total - cell_total) + reg)` as the
# variable side - with the mask `is_active_per_point` set false (and the variable value zeroed) when the cell's
# baseline metacell is missing or the punctuation denominator is non-positive. `baseline_metacell_index_per_cell` and
# the per-baseline-metacell aggregates (`total_UMIs_per_baseline_metacell`,
# `UMIs_per_relevant_gene_per_baseline_metacell`) encode "Phase 1's choice of K per block" - the baseline candidate's cluster
# assignments. The returned `group_index_per_block` tells the indirect-gather scoring API which group inside this base
# block corresponds to each contributing block.
function build_grouped_for_base_block(;
    n_blocks::Integer,
    n_cells::Integer,
    base_block_index_per_cell::AbstractVector{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    block_cell_indices_per_block::AbstractVector,
    is_in_base_neighborhood_for_base_block::AbstractVector{Bool},
    baseline_metacell_index_per_cell::AbstractVector{<:Integer},
    total_UMIs_per_baseline_metacell::AbstractVector{<:Integer},
    UMIs_per_relevant_gene_per_baseline_metacell::AbstractMatrix{<:Integer},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    strong_gene_indices::AbstractVector{<:Integer},
    friend_gene_indices::AbstractVector{<:Integer},
    relevant_gene_position_per_gene::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat,
)::Tuple{GroupedSeriesCorrelations{Float32}, Vector{Int}}
    n_series = length(strong_gene_indices)
    @assert length(friend_gene_indices) == n_series

    # Count how many cells each block contributes to this base block's neighborhood.
    n_points_per_block = zeros(Int, n_blocks)
    for cell_index in 1:n_cells
        base_block_of_cell = base_block_index_per_cell[cell_index]
        block_of_cell = block_index_per_cell[cell_index]
        if base_block_of_cell > 0 && block_of_cell > 0 && is_in_base_neighborhood_for_base_block[base_block_of_cell]
            n_points_per_block[block_of_cell] += 1
        end
    end

    # Assign group indices to contributing blocks in `block_index` order.
    group_index_per_block = zeros(Int, n_blocks)
    n_groups = 0
    for block_index in 1:n_blocks
        if n_points_per_block[block_index] > 0
            n_groups += 1
            group_index_per_block[block_index] = n_groups
        end
    end

    n_points = sum(n_points_per_block)
    n_points_per_group = Vector{Int}(undef, n_groups)
    for block_index in 1:n_blocks
        group = group_index_per_block[block_index]
        if group > 0
            n_points_per_group[group] = n_points_per_block[block_index]
        end
    end

    # Lay out cells contiguously by group, preserving the order in `block_cell_indices_per_block` within each block.
    cell_indices_in_order = Vector{Int}(undef, n_points)
    write_position = 1
    for block_index in 1:n_blocks
        if group_index_per_block[block_index] == 0
            continue
        end
        for cell_index in block_cell_indices_per_block[block_index]
            base_block_of_cell = base_block_index_per_cell[cell_index]
            if base_block_of_cell > 0 && is_in_base_neighborhood_for_base_block[base_block_of_cell]
                cell_indices_in_order[write_position] = cell_index
                write_position += 1
            end
        end
    end
    @assert write_position == n_points + 1

    fixed_per_point_per_series = Matrix{Float32}(undef, n_points, n_series)
    variable_per_point_per_series = Matrix{Float32}(undef, n_points, n_series)
    is_active_per_point = Vector{Bool}(undef, n_points)
    regularization = Float32(gene_fraction_regularization)
    for point_index in 1:n_points
        cell_index = cell_indices_in_order[point_index]
        cell_total = Int(total_UMIs_per_cell[cell_index])
        baseline_metacell_index = Int(baseline_metacell_index_per_cell[cell_index])
        baseline_metacell_total = baseline_metacell_index > 0 ? Int(total_UMIs_per_baseline_metacell[baseline_metacell_index]) : 0
        is_active = baseline_metacell_index > 0 && baseline_metacell_total > cell_total
        is_active_per_point[point_index] = is_active

        for series_index in 1:n_series
            strong_gene_index = strong_gene_indices[series_index]
            friend_gene_index = friend_gene_indices[series_index]
            friend_relevant_position = relevant_gene_position_per_gene[friend_gene_index]

            cell_strong_UMIs = Int(UMIs_per_cell_per_gene[cell_index, strong_gene_index])
            fixed_per_point_per_series[point_index, series_index] =
                log2(Float32(cell_strong_UMIs) / Float32(cell_total) + regularization)

            if is_active
                baseline_metacell_friend_UMIs = Int(UMIs_per_relevant_gene_per_baseline_metacell[baseline_metacell_index, friend_relevant_position])
                cell_friend_UMIs = Int(UMIs_per_cell_per_gene[cell_index, friend_gene_index])
                punctuated_friend = baseline_metacell_friend_UMIs - cell_friend_UMIs
                punctuated_total = baseline_metacell_total - cell_total
                variable_per_point_per_series[point_index, series_index] =
                    log2(Float32(punctuated_friend) / Float32(punctuated_total) + regularization)
            else
                variable_per_point_per_series[point_index, series_index] = 0.0f0
            end
        end
    end

    grouped = GroupedSeriesCorrelations(
        fixed_per_point_per_series,
        variable_per_point_per_series,
        n_points_per_group;
        is_active_per_point,
    )

    return (grouped, group_index_per_block)
end

# Given the per-block Phase-1 `LocalClusters`, compute the per-baseline-metacell aggregates
# (`total_UMIs_per_baseline_metacell` and `UMIs_per_relevant_gene_per_baseline_metacell`), the per-cell
# `baseline_metacell_index_per_cell` (the cell's baseline metacell index in the global array, or 0 if the cell's
# Phase 1 cluster is too-small), and `first_baseline_metacell_per_block` indicating where each block's baseline
# metacells begin in the global array. For a cell in some block whose Phase 1 cluster is `cluster_index`, the
# baseline metacell index is `first_baseline_metacell_per_block[block_index] + cluster_index - 1`. Blocks with
# `local_clusters === nothing` contribute no baseline metacells.
function compute_baseline_metacell_aggregates(
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    gene_index_per_relevant_gene_position::AbstractVector{<:Integer},
    n_cells::Integer,
)::Tuple{Vector{UInt32}, Matrix{UInt32}, Vector{UInt32}, Vector{Int}}
    n_blocks = length(local_clusters_per_block)
    n_relevant_genes = length(gene_index_per_relevant_gene_position)

    first_baseline_metacell_per_block = Vector{Int}(undef, n_blocks)
    n_baseline_metacells = 0
    for block_index in 1:n_blocks
        local_clusters = local_clusters_per_block[block_index]
        first_baseline_metacell_per_block[block_index] = n_baseline_metacells + 1
        if local_clusters !== nothing
            n_baseline_metacells += local_clusters.n_clusters
        end
    end

    total_UMIs_per_baseline_metacell = Vector{UInt32}(undef, n_baseline_metacells)
    UMIs_per_relevant_gene_per_baseline_metacell = Matrix{UInt32}(undef, n_baseline_metacells, n_relevant_genes)
    baseline_metacell_index_per_cell = zeros(UInt32, n_cells)
    cell_indices_per_baseline_metacell = [Int[] for _ in 1:n_baseline_metacells]

    for block_index in 1:n_blocks
        local_clusters = local_clusters_per_block[block_index]
        if local_clusters === nothing
            continue
        end
        block_first_baseline_metacell = first_baseline_metacell_per_block[block_index]
        @inbounds for block_cell_position in eachindex(local_clusters.cluster_index_per_block_cell)
            cluster_index = local_clusters.cluster_index_per_block_cell[block_cell_position]
            cell_index = local_clusters.block_cell_indices[block_cell_position]
            baseline_metacell_index = block_first_baseline_metacell + cluster_index - 1
            push!(cell_indices_per_baseline_metacell[baseline_metacell_index], cell_index)
            baseline_metacell_index_per_cell[cell_index] =
                local_clusters.is_too_small_per_cluster[cluster_index] ? UInt32(0) : UInt32(baseline_metacell_index)
        end
    end

    parallel_loop_wo_rng(
        1:n_baseline_metacells;
        policy = :greedy,
        name = "compute_baseline_metacell_aggregates",
        progress = DebugProgress(n_baseline_metacells; group = :mcs_loops, desc = "compute_baseline_metacell_aggregates"),
    ) do baseline_metacell_index
        cell_indices = cell_indices_per_baseline_metacell[baseline_metacell_index]
        baseline_metacell_total = zero(UInt32)
        @inbounds for cell_index in cell_indices
            baseline_metacell_total += UInt32(total_UMIs_per_cell[cell_index])
        end
        total_UMIs_per_baseline_metacell[baseline_metacell_index] = baseline_metacell_total
        @inbounds for relevant_gene_position in 1:n_relevant_genes
            gene_index = gene_index_per_relevant_gene_position[relevant_gene_position]
            gene_sum = zero(UInt32)
            for cell_index in cell_indices
                gene_sum += UInt32(UMIs_per_cell_per_gene[cell_index, gene_index])
            end
            UMIs_per_relevant_gene_per_baseline_metacell[baseline_metacell_index, relevant_gene_position] = gene_sum
        end
        return nothing
    end

    return (total_UMIs_per_baseline_metacell, UMIs_per_relevant_gene_per_baseline_metacell, baseline_metacell_index_per_cell, first_baseline_metacell_per_block)
end

# Per walkable block, precompute the indirection vectors the delta-correlation scoring will need: which base blocks'
# correlations the block can affect, the block's friend-gene subspace (column axis for the per-(block, K-candidate)
# log-fill cache), and per (block, affected base block) the (column position in the block's friend subspace, the
# block's `block_cell_position` per point) maps.
#   * `affected_base_block_indices_per_walkable_block[walkable_position]`: list of base blocks whose neighborhoods overlap any of the
#     block's cells.
#   * `gene_index_per_friend_position_per_walkable_block[walkable_position]`: the block's friend-gene subspace (global gene index per
#     column).
#   * `relevant_gene_position_per_friend_position_per_walkable_block[walkable_position]`: same column →
#     `UMIs_per_relevant_gene_per_baseline_metacell` column index for the baseline-metacell friend-gene lookup.
#   * `friend_position_per_series_per_base_block_per_walkable_block[walkable_position][base_block_index]`: per (block, base block), the
#     block's friend column for each of the base block's strong-friend series.
#   * `block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position][base_block_index]`: per (block, base block), the block's
#     `block_cell_position`s of the cells in the base block's neighborhood, in the block's canonical block-cell order.
#     Row indirection for the indirect-gather query.
function precompute_walkable_indirection(
    walkable_block_indices::AbstractVector{<:Integer},
    block_cell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    base_block_index_per_cell::AbstractVector{<:Integer},
    is_in_base_neighborhood_per_other_base_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    friend_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    relevant_gene_position_per_gene::AbstractVector{<:Integer},
    n_base_blocks::Integer,
    n_genes::Integer,
)::Tuple{
    Vector{Vector{Int}},
    Vector{Vector{Int}},
    Vector{Vector{Int}},
    Vector{Dict{Int, Vector{Int}}},
    Vector{Dict{Int, Vector{Int}}},
}
    n_walkable_blocks = length(walkable_block_indices)

    affected_base_block_indices_per_walkable_block = Vector{Vector{Int}}(undef, n_walkable_blocks)
    gene_index_per_friend_position_per_walkable_block = Vector{Vector{Int}}(undef, n_walkable_blocks)
    relevant_gene_position_per_friend_position_per_walkable_block = Vector{Vector{Int}}(undef, n_walkable_blocks)
    friend_position_per_series_per_base_block_per_walkable_block =
        Vector{Dict{Int, Vector{Int}}}(undef, n_walkable_blocks)
    block_cell_position_per_point_per_base_block_per_walkable_block =
        Vector{Dict{Int, Vector{Int}}}(undef, n_walkable_blocks)

    is_base_block_of_block_per_thread = [BitVector(undef, n_base_blocks) for _ in 1:maxthreadid()]
    is_affected_base_block_per_thread = [BitVector(undef, n_base_blocks) for _ in 1:maxthreadid()]
    is_friend_gene_for_walkable_block_per_thread = [BitVector(undef, n_genes) for _ in 1:maxthreadid()]
    friend_position_per_gene_per_thread = [zeros(Int, n_genes) for _ in 1:maxthreadid()]

    n_block_cells_per_walkable_block =
        [length(block_cell_indices_per_block[block_index]) for block_index in walkable_block_indices]

    parallel_loop_wo_rng(
        1:n_walkable_blocks;
        policy = :static_greedy,
        order = sortperm(n_block_cells_per_walkable_block; rev = true),  # NOJET
        name = "precompute_walkable_indirection",
        progress = DebugProgress(n_walkable_blocks; group = :mcs_loops, desc = "precompute_walkable_indirection"),
    ) do walkable_position
        block_index = walkable_block_indices[walkable_position]
        block_cell_indices_of_walkable_block = block_cell_indices_per_block[block_index]
        n_block_cells = length(block_cell_indices_of_walkable_block)

        # Affected base blocks for this walkable block.
        is_base_block_of_block = is_base_block_of_block_per_thread[threadid()]
        is_affected_base_block = is_affected_base_block_per_thread[threadid()]
        fill!(is_base_block_of_block, false)
        for cell_index in block_cell_indices_of_walkable_block
            base_block_index = base_block_index_per_cell[cell_index]
            if base_block_index > 0
                is_base_block_of_block[base_block_index] = true
            end
        end
        fill!(is_affected_base_block, false)
        @foreach_true_index is_base_block_of_block base_block_index begin  # NOLINT
            @views is_affected_base_block .|=
                is_in_base_neighborhood_per_other_base_block_per_base_block[base_block_index, :]  # NOLINT
        end
        affected_base_block_indices = findall(is_affected_base_block)
        affected_base_block_indices_per_walkable_block[walkable_position] = affected_base_block_indices

        # Friend-gene subspace for this walkable block: union of friend genes across its affected base blocks.
        is_friend_gene_for_walkable_block = is_friend_gene_for_walkable_block_per_thread[threadid()]
        friend_position_per_gene = friend_position_per_gene_per_thread[threadid()]
        fill!(is_friend_gene_for_walkable_block, false)
        for base_block_index in affected_base_block_indices
            for friend_gene_index in friend_gene_indices_per_base_block[base_block_index]
                is_friend_gene_for_walkable_block[friend_gene_index] = true
            end
        end
        fill!(friend_position_per_gene, 0)
        gene_index_per_friend_position = Int[]
        for gene_index in 1:n_genes
            if is_friend_gene_for_walkable_block[gene_index]
                push!(gene_index_per_friend_position, gene_index)
                friend_position_per_gene[gene_index] = length(gene_index_per_friend_position)
            end
        end
        gene_index_per_friend_position_per_walkable_block[walkable_position] = gene_index_per_friend_position
        relevant_gene_position_per_friend_position_per_walkable_block[walkable_position] =
            [relevant_gene_position_per_gene[gene_index] for gene_index in gene_index_per_friend_position]

        friend_position_per_series_per_base_block = Dict{Int, Vector{Int}}()
        block_cell_position_per_point_per_base_block = Dict{Int, Vector{Int}}()
        for base_block_index in affected_base_block_indices
            friend_position_per_series_per_base_block[base_block_index] = [
                friend_position_per_gene[friend_gene_index]
                for friend_gene_index in friend_gene_indices_per_base_block[base_block_index]
            ]
            @views is_in_base_neighborhood_for_base_block =
                is_in_base_neighborhood_per_other_base_block_per_base_block[:, base_block_index]
            block_cell_position_per_point = Int[]
            for block_cell_position in 1:n_block_cells
                cell_index = block_cell_indices_of_walkable_block[block_cell_position]
                base_block_of_cell = base_block_index_per_cell[cell_index]
                if base_block_of_cell > 0 && is_in_base_neighborhood_for_base_block[base_block_of_cell]
                    push!(block_cell_position_per_point, block_cell_position)
                end
            end
            block_cell_position_per_point_per_base_block[base_block_index] = block_cell_position_per_point
        end
        friend_position_per_series_per_base_block_per_walkable_block[walkable_position] =
            friend_position_per_series_per_base_block
        block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position] =
            block_cell_position_per_point_per_base_block
        return nothing
    end

    return (
        affected_base_block_indices_per_walkable_block,
        gene_index_per_friend_position_per_walkable_block,
        relevant_gene_position_per_friend_position_per_walkable_block,
        friend_position_per_series_per_base_block_per_walkable_block,
        block_cell_position_per_point_per_base_block_per_walkable_block,
    )
end

# TODOX: only used by the diagnostic NEXT / BETTER / PERFECT K logs in `kmeans_with_sizes`. Returns extrema for cells
# per cluster, total UMIs per cluster, and dispersion per cluster, taken across the solution's active K range. The
# dispersion extrema skip clusters whose cell count is below `min_cluster_size` (outliers) - those have undefined or
# degenerate dispersion and would otherwise drive the min to 0.
function todox_solution_extrema(
    solution::KmeansSolution,
    min_cluster_size::Real,
)::Tuple{Tuple{Int, Int}, Tuple{Float64, Float64}, Tuple{Float32, Float32}}
    cells_range = extrema(@view solution.counts[1:solution.k])
    UMIs_range = extrema(@view solution.weight_per_cluster[1:solution.k])
    dispersion_min = Float32(Inf)
    dispersion_max = -Float32(Inf)
    @inbounds for cluster_index in 1:solution.k
        if solution.counts[cluster_index] >= min_cluster_size
            dispersion = solution.cells_dispersion_per_cluster[cluster_index]
            dispersion_min = min(dispersion_min, dispersion)
            dispersion_max = max(dispersion_max, dispersion)
        end
    end
    return (cells_range, UMIs_range, (dispersion_min, dispersion_max))
end

# TODOX: dumps every non-outlier cluster's (count, weight, dispersion) so the offending cluster with dispersion = 0 is
# identifiable when the corresponding assertion in `kmeans_with_sizes` fires. For each cluster with `dispersion == 0`
# AND `count >= min_cluster_size`, also dumps per-found-module mean normalized UMIs (= what
# `maximal_cells_dispersion_of_modules!` computes and compares against `min_module_UMIs` to decide whether to skip the
# module). That distinguishes the three known causes of a zero non-outlier dispersion: (a) count < 2 (impossible if
# `min_cluster_size >= 2`); (b) `is_found_per_module` is all false (n_found = 0 in the dump); (c) every found module's
# mean is below `min_module_UMIs` (n_found > 0 and n_passing = 0 in the dump).
function todox_describe_solution_clusters(
    solution::KmeansSolution,
    min_cluster_size::Real,
    sizes_buffers::KmeansSizesBuffers,
    dispersion_context::DispersionContext,
    context::AbstractString,
)::Nothing
    @warn "TODOX - $(context) non-outlier clusters dump (min_cluster_size=$(min_cluster_size)):"
    for cluster_index in 1:solution.k
        count = solution.counts[cluster_index]
        if count < min_cluster_size
            continue
        end
        weight = solution.weight_per_cluster[cluster_index]
        dispersion = solution.cells_dispersion_per_cluster[cluster_index]
        @warn "  cluster $(cluster_index): count=$(count) weight=$(round(weight / 1000; digits=1))K dispersion=$(dispersion)"
        if dispersion == 0
            todox_describe_cluster_modules(solution, cluster_index, sizes_buffers, dispersion_context)
            @assert false  # TODOX
        end
    end
    return nothing
end

# TODOX: mirrors the per-found-module setup of `maximal_cells_dispersion_of_modules!` for a single cluster but logs the
# `mean_normalized_module_UMIs` it computes per found module (the value compared against `min_module_UMIs` to decide
# whether the module is `continue`d). Trashes `sizes_buffers` scratch (max_split_point_indices, total_UMIs_per_max_-
# cluster_cell, normalized_factor_per_max_cluster_cell, is_gene_in_module); only called immediately before an assertion
# fires so the trash is irrelevant.
function todox_describe_cluster_modules(
    solution::KmeansSolution,
    cluster_index::Integer,
    sizes_buffers::KmeansSizesBuffers,
    dispersion_context::DispersionContext,
)::Nothing
    n_block_cells = length(dispersion_context.block_cell_indices)
    n_cluster_cells = 0
    @inbounds for point_index in 1:n_block_cells
        if solution.assignments[point_index] == cluster_index
            n_cluster_cells += 1
            sizes_buffers.max_split_point_indices[n_cluster_cells] =
                dispersion_context.block_cell_indices[point_index]
        end
    end
    @views indices_of_cluster_cells = sizes_buffers.max_split_point_indices[1:n_cluster_cells]
    @views total_UMIs_per_cells = sizes_buffers.total_UMIs_per_max_cluster_cell[1:n_cluster_cells]
    @views normalized_factor_per_cells = sizes_buffers.normalized_factor_per_max_cluster_cell[1:n_cluster_cells]
    normalized_factor_per_cells .= getindex.(Ref(dispersion_context.total_UMIs_per_cell), indices_of_cluster_cells)
    total_UMIs_per_cells .= normalized_factor_per_cells
    normalized_total_UMIs = quantile!(total_UMIs_per_cells, dispersion_context.normalized_UMIs_quantile)
    @. normalized_factor_per_cells = normalized_total_UMIs / normalized_factor_per_cells

    n_modules = length(dispersion_context.is_found_per_module)
    n_genes = length(dispersion_context.module_index_per_gene)
    is_gene_in_module = sizes_buffers.is_gene_in_module
    n_found = 0
    n_passing = 0
    max_mean_for_found = 0.0
    for module_index in 1:n_modules
        if !dispersion_context.is_found_per_module[module_index]
            continue
        end
        n_found += 1
        @. is_gene_in_module = dispersion_context.module_index_per_gene == module_index
        sum_normalized = 0.0
        @inbounds for cell_position in 1:n_cluster_cells
            cell_index = indices_of_cluster_cells[cell_position]
            module_UMIs = 0
            for gene_index in 1:n_genes
                if is_gene_in_module[gene_index]
                    module_UMIs += Int(dispersion_context.UMIs_per_cell_per_gene[cell_index, gene_index])
                end
            end
            sum_normalized += module_UMIs * normalized_factor_per_cells[cell_position]
        end
        mean_normalized = sum_normalized / n_cluster_cells
        passing = mean_normalized >= dispersion_context.min_module_UMIs
        if passing
            n_passing += 1
        end
        max_mean_for_found = max(max_mean_for_found, mean_normalized)
        @warn "    module $(module_index): mean_normalized=$(round(mean_normalized; digits=3)) passing=$(passing)"
    end
    @warn "    summary: n_found=$(n_found) n_passing=$(n_passing) max_mean=$(round(max_mean_for_found; digits=3)) min_module_UMIs=$(dispersion_context.min_module_UMIs)"
    return nothing
end

function compute_z_score_per_found_module_per_region_cell!(;
    z_score_per_found_module_per_region_cell::AbstractMatrix{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    is_found_per_module::BitVector,
    is_in_region_per_cell::BitVector,
    module_index_per_gene::AbstractVector{<:Integer},
    is_gene_in_module::BitVector,
    mean_linear_fraction_in_neighborhood_cells_per_module::AbstractVector{<:AbstractFloat},
    std_linear_fraction_in_neighborhood_cells_per_module::AbstractVector{<:AbstractFloat},
)::Nothing
    z_score_per_found_module_per_region_cell .= 0
    @foreach_true_index_position is_found_per_module module_index found_module_position begin  # NOLINT
        @. is_gene_in_module = module_index_per_gene == module_index  # NOLINT
        @foreach_true_index_position is_in_region_per_cell cell_index region_cell_position begin  # NOLINT
            @foreach_true_index is_gene_in_module gene_index begin  # NOLINT
                z_score_per_found_module_per_region_cell[found_module_position, region_cell_position] +=  # NOLINT
                    UMIs_per_cell_per_gene[cell_index, gene_index]  # NOLINT
            end
            z_score_per_found_module_per_region_cell[found_module_position, region_cell_position] /=  # NOLINT
                total_UMIs_per_cell[cell_index]  # NOLINT
        end
        mean_linear_fraction_in_neighborhood_cells = mean_linear_fraction_in_neighborhood_cells_per_module[module_index]  # NOLINT
        std_linear_fraction_in_neighborhood_cells = std_linear_fraction_in_neighborhood_cells_per_module[module_index]  # NOLINT
        z_score_per_found_module_per_region_cell[found_module_position, :] .-=  # NOLINT
            mean_linear_fraction_in_neighborhood_cells
        z_score_per_found_module_per_region_cell[found_module_position, :] ./= std_linear_fraction_in_neighborhood_cells  # NOLINT
    end
    return nothing
end

# Copy all fields and the active-range data from source to target.
function copy_solution!(target::KmeansSolution{T}, source::KmeansSolution{T})::Nothing where {T <: AbstractFloat}
    @assert source.is_filled
    target.is_filled = true
    target.n_dims = source.n_dims
    target.n_points = source.n_points
    target.k = source.k
    target.n_too_tight = source.n_too_tight
    target.n_too_small = source.n_too_small
    target.n_too_wide = source.n_too_wide
    target.n_too_large = source.n_too_large
    @views target.centers[1:source.n_dims, 1:source.k] .= source.centers[1:source.n_dims, 1:source.k]
    @views target.assignments[1:source.n_points] .= source.assignments[1:source.n_points]
    @views target.counts[1:source.k] .= source.counts[1:source.k]
    @views target.weight_per_cluster[1:source.k] .= source.weight_per_cluster[1:source.k]
    @views target.cells_dispersion_per_cluster[1:source.k] .= source.cells_dispersion_per_cluster[1:source.k]
    return nothing
end

# Recompute weight_per_cluster from the active assignments.
function compute_weight_per_cluster!(solution::KmeansSolution, weight_per_point::AbstractVector{<:Real})::Nothing
    @views fill!(solution.weight_per_cluster[1:solution.k], 0.0)
    @inbounds for point_index in eachindex(weight_per_point)
        cluster_index = solution.assignments[point_index]
        solution.weight_per_cluster[cluster_index] += weight_per_point[point_index]
    end
    return nothing
end

# Recompute the maximal cells_dispersion (across the dispersion context's found modules) for a single cluster and store
# it in solution.cells_dispersion_per_cluster.
function compute_cluster_dispersion!(
    solution::KmeansSolution,
    cluster_index::Integer,
    sizes_buffers::KmeansSizesBuffers,
    dispersion_context::DispersionContext,
)::Nothing
    n_block_cells = length(dispersion_context.block_cell_indices)
    n_cluster_cells = 0
    @inbounds for point_index in 1:n_block_cells
        if solution.assignments[point_index] == cluster_index
            n_cluster_cells += 1
            sizes_buffers.max_split_point_indices[n_cluster_cells] = dispersion_context.block_cell_indices[point_index]
        end
    end

    @views indices_of_cluster_cells = sizes_buffers.max_split_point_indices[1:n_cluster_cells]
    solution.cells_dispersion_per_cluster[cluster_index] = maximal_cells_dispersion_of_modules!(;
        cells_dispersion_per_module = sizes_buffers.cells_dispersion_per_max_module,
        mean_normalized_per_module = sizes_buffers.mean_normalized_per_max_module,
        indices_of_cells = indices_of_cluster_cells,
        UMIs_per_cell_per_gene = dispersion_context.UMIs_per_cell_per_gene,
        total_UMIs_per_cell = dispersion_context.total_UMIs_per_cell,
        is_found_per_module = dispersion_context.is_found_per_module,
        module_index_per_gene = dispersion_context.module_index_per_gene,
        is_gene_in_module = sizes_buffers.is_gene_in_module,
        total_UMIs_per_max_cells = sizes_buffers.total_UMIs_per_max_cluster_cell,
        normalized_factor_per_max_cells = sizes_buffers.normalized_factor_per_max_cluster_cell,
        normalized_UMIs_quantile = dispersion_context.normalized_UMIs_quantile,
        min_module_UMIs = dispersion_context.min_module_UMIs,
    )
    return nothing
end

# Recompute the maximal cells_dispersion for every active cluster in the solution.
function compute_all_dispersions!(
    solution::KmeansSolution,
    sizes_buffers::KmeansSizesBuffers,
    dispersion_context::DispersionContext,
)::Nothing
    for cluster_index in 1:solution.k
        compute_cluster_dispersion!(solution, cluster_index, sizes_buffers, dispersion_context)
    end
    return nothing
end


# Fill `variable_per_block_cell_per_friend_position` and `is_active_per_block_cell` with `candidate`'s punctuated
# cluster-mean log-fractions over the walkable block's cells × friend genes; cells in too-small or zero-punctuated-
# total clusters are masked off and their variable values zeroed. Uses `delta_context`'s per-thread cluster-aggregate
# scratches. Allocation-free; the supplied buffers may be per-thread (for Phase 2 evaluation, overwritten per work
# item) or per-walkable-block (for commit, kept alive across the Phase 2 → commit hand-off).
function populate_candidate_scratches!(
    candidate::SolutionCandidate,
    n_block_cells::Integer,
    delta_context::DeltaCorrelationContext,
    variable_per_block_cell_per_friend_position::AbstractMatrix{Float32},
    is_active_per_block_cell::AbstractVector{Bool},
)::Nothing
    n_friend_positions = length(delta_context.gene_index_per_friend_position)
    k = candidate.k

    total_UMIs_per_max_cluster = delta_context.total_UMIs_per_max_cluster
    UMIs_per_friend_position_per_max_cluster = delta_context.UMIs_per_friend_position_per_max_cluster
    @views fill!(total_UMIs_per_max_cluster[1:k], 0.0)
    @views fill!(UMIs_per_friend_position_per_max_cluster[1:n_friend_positions, 1:k], 0.0)

    # Per K-candidate cluster: total UMIs of the block's cells in the cluster + per-friend-gene UMI sums (only over
    # the walkable block's friend-gene subspace, not the full relevant gene set).
    @inbounds for block_cell_position in 1:n_block_cells
        cluster_index = candidate.assignments[block_cell_position]
        cell_index = delta_context.block_cell_indices[block_cell_position]
        total_UMIs_per_max_cluster[cluster_index] += Float64(delta_context.total_UMIs_per_cell[cell_index])
        for friend_position in 1:n_friend_positions
            gene_index = delta_context.gene_index_per_friend_position[friend_position]
            UMIs_per_friend_position_per_max_cluster[friend_position, cluster_index] +=
                Float64(delta_context.UMIs_per_cell_per_gene[cell_index, gene_index])
        end
    end

    regularization = delta_context.gene_fraction_regularization
    @inbounds for block_cell_position in 1:n_block_cells
        cluster_index = candidate.assignments[block_cell_position]
        cell_index = delta_context.block_cell_indices[block_cell_position]
        cell_total = Float64(delta_context.total_UMIs_per_cell[cell_index])
        cluster_total_UMIs = total_UMIs_per_max_cluster[cluster_index]
        punctuated_total = cluster_total_UMIs - cell_total
        is_active = candidate.counts[cluster_index] >= 2 && punctuated_total > 0
        is_active_per_block_cell[block_cell_position] = is_active
        for friend_position in 1:n_friend_positions
            if is_active
                gene_index = delta_context.gene_index_per_friend_position[friend_position]
                cell_friend_UMIs = Float64(delta_context.UMIs_per_cell_per_gene[cell_index, gene_index])
                cluster_friend_UMIs = UMIs_per_friend_position_per_max_cluster[friend_position, cluster_index]
                punctuated_friend = cluster_friend_UMIs - cell_friend_UMIs
                variable_per_block_cell_per_friend_position[block_cell_position, friend_position] =
                    log2(Float32(punctuated_friend / punctuated_total) + regularization)
            else
                variable_per_block_cell_per_friend_position[block_cell_position, friend_position] = 0.0f0
            end
        end
    end
    return nothing
end

# Score one candidate K-solution for a walkable block: returns the delta of grouped mean correlation against the
# current baseline, summed across the block's `affected_base_block_indices`. Per affected base block, the mean
# correlation of the grouped correlations when the block's group is replaced (under `candidate`'s assignments) minus
# the cached baseline mean correlation for that base block. Uses `delta_context`'s per-thread variable + mask scratches
# for the gather, sized to the walkable blocks' global maxima.
function compute_delta_correlation(
    candidate::SolutionCandidate,
    n_block_cells::Integer,
    delta_context::DeltaCorrelationContext,
)::Float64
    variable_per_block_cell_per_friend_position = delta_context.variable_per_block_cell_per_friend_position
    is_active_per_block_cell = delta_context.is_active_per_block_cell
    populate_candidate_scratches!(
        candidate,
        n_block_cells,
        delta_context,
        variable_per_block_cell_per_friend_position,
        is_active_per_block_cell,
    )

    # Sum per-base-block delta: query each affected base block's grouped correlations with the block's candidate
    # values via the 6-arg indirect-gather form, subtract the cached baseline mean.
    delta = 0.0
    walkable_block_index = delta_context.walkable_block_index
    for base_block_index in delta_context.affected_base_block_indices
        group_index_in_base_block =
            delta_context.group_index_per_block_per_base_block[walkable_block_index, base_block_index]
        if group_index_in_base_block == 0
            continue
        end
        grouped_for_base_block = delta_context.grouped_per_base_block[base_block_index]
        block_cell_position_per_point = delta_context.block_cell_position_per_point_per_base_block[base_block_index]
        friend_position_per_series = delta_context.friend_position_per_series_per_base_block[base_block_index]

        new_mean_correlation = mean_correlation_if_group_replaced(
            grouped_for_base_block,
            group_index_in_base_block,
            variable_per_block_cell_per_friend_position,
            is_active_per_block_cell,
            block_cell_position_per_point,
            friend_position_per_series,
        )
        delta +=
            Float64(new_mean_correlation) -
            Float64(delta_context.baseline_mean_correlation_per_base_block[base_block_index])
    end

    return delta
end

# Recompute n_too_tight (clusters whose maximal cells_dispersion is below dispersion_context.min_cells_dispersion),
# n_too_small (clusters whose count is below min_cluster_size OR whose weight is below min_cluster_weight),
# n_too_wide (clusters whose maximal cells_dispersion is above dispersion_context.max_cells_dispersion), and
# n_too_large (clusters whose count is above max_cluster_size). Assumes solution.cells_dispersion_per_cluster is
# already populated for the active clusters.
function update_size_statistics!(
    solution::KmeansSolution,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,  # NOLINT
    dispersion_context::DispersionContext,
)::Nothing
    max_cluster_size = 2 * solution.n_points / solution.k  # TODOX
    n_too_tight = 0
    n_too_small = 0
    n_too_wide = 0
    n_too_large = 0
    @inbounds for cluster_index in 1:solution.k
        if solution.counts[cluster_index] < min_cluster_size ||
           solution.weight_per_cluster[cluster_index] < min_cluster_weight
            n_too_small += 1
        end
        if solution.counts[cluster_index] > max_cluster_size
            n_too_large += 1
        end
        max_dispersion = solution.cells_dispersion_per_cluster[cluster_index]
        if max_dispersion < dispersion_context.min_cells_dispersion
            n_too_tight += 1
        end
        if max_dispersion > dispersion_context.max_cells_dispersion
            n_too_wide += 1
        end
    end
    solution.n_too_tight = n_too_tight
    solution.n_too_small = n_too_small
    solution.n_too_wide = n_too_wide
    solution.n_too_large = n_too_large
    return nothing
end

function penalty(solution::KmeansSolution)::Tuple{Int, Int, Int}
    n_too_tight_or_small = solution.n_too_tight + solution.n_too_small
    n_too_wide_or_large = solution.n_too_wide + solution.n_too_large
    # Third slot is `K` for both perfect and imperfect: smaller K wins on the tiebreaker. Perfect candidates are
    # collected per block by Phase 1 and ranked downstream against the delta-correlation, not here.
    return (n_too_tight_or_small, n_too_wide_or_large, solution.k)
end

# Split a single cluster into two using a 2-cluster kmeans on its points. Replaces the split cluster's center with one
# sub-center, appends the other at a fresh last position (incrementing `candidate_solution.k`), reassigns the points,
# and updates the counts of both sub-clusters.
function split_one_cluster!(
    candidate_solution::KmeansSolution{T},
    old_cluster_index::Integer,
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Nothing where {T <: AbstractFloat}
    n_cluster_points = 0
    @inbounds for point_index in 1:candidate_solution.n_points
        if candidate_solution.assignments[point_index] == old_cluster_index
            n_cluster_points += 1
            sizes_buffers.max_split_point_indices[n_cluster_points] = point_index
        end
    end
    @assert n_cluster_points >= 2

    @views split_point_indices = sizes_buffers.max_split_point_indices[1:n_cluster_points]
    @views split_points = sizes_buffers.max_split_points[1:candidate_solution.n_dims, 1:n_cluster_points]
    @views split_points .= values_of_points[:, split_point_indices]

    sub_kmeans_result = flame_timed("kmeans") do
        return kmeans_in_rounds(split_points, 2; buffers = sizes_buffers.kmeans_buffers, rounds = kmeans_rounds, rng)
    end

    candidate_solution.k += 1
    new_cluster_index = candidate_solution.k
    @views candidate_solution.centers[1:candidate_solution.n_dims, old_cluster_index] .= sub_kmeans_result.centers[:, 1]
    @views candidate_solution.centers[1:candidate_solution.n_dims, new_cluster_index] .= sub_kmeans_result.centers[:, 2]
    sub_assignments = assignments(sub_kmeans_result)  # NOJET
    sub_counts = counts(sub_kmeans_result)  # NOJET
    old_weight = 0.0
    new_weight = 0.0
    @inbounds for cluster_point_position in 1:n_cluster_points
        point_index = split_point_indices[cluster_point_position]
        if sub_assignments[cluster_point_position] == 2
            candidate_solution.assignments[point_index] = new_cluster_index
            new_weight += weight_per_point[point_index]
        else
            old_weight += weight_per_point[point_index]
        end
    end
    candidate_solution.counts[old_cluster_index] = sub_counts[1]
    candidate_solution.counts[new_cluster_index] = sub_counts[2]
    candidate_solution.weight_per_cluster[old_cluster_index] = old_weight
    candidate_solution.weight_per_cluster[new_cluster_index] = new_weight
    compute_cluster_dispersion!(candidate_solution, old_cluster_index, sizes_buffers, dispersion_context)
    compute_cluster_dispersion!(candidate_solution, new_cluster_index, sizes_buffers, dispersion_context)
    return nothing
end

# Build a candidate solution by splitting the single most-deserving non-rejected cluster of `current_solution` (running
# a 2-cluster kmeans on its points; one sub-center replaces the original, the other is appended). Prefer the largest
# too-large cluster; if none qualifies, fall back to the widest too-wide cluster. A cluster is considered splittable
# only when its count and weight are at least twice the respective minimums (so that each half can plausibly clear
# them), and only when its bit in `is_rejected_for_split_per_cluster` is false. Returns the chosen cluster index in
# `current_solution`, or 0 if no eligible cluster remains.
function split_largest_cluster!(
    candidate_solution::KmeansSolution{T},
    current_solution::KmeansSolution{T},
    is_rejected_for_split_per_cluster::AbstractVector{Bool},
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,  # NOLINT
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Integer where {T <: AbstractFloat}
    k = current_solution.k
    max_cluster_size = 2 * current_solution.n_points / current_solution.k  # TODOX

    min_splittable_size = 2 * min_cluster_size
    min_splittable_weight = 2 * min_cluster_weight

    target_cluster = 0
    target_size = max_cluster_size
    @inbounds for cluster_index in 1:k
        if !is_rejected_for_split_per_cluster[cluster_index] &&
           current_solution.counts[cluster_index] > target_size &&
           current_solution.counts[cluster_index] >= min_splittable_size &&
           current_solution.weight_per_cluster[cluster_index] >= min_splittable_weight
            target_cluster = cluster_index
            target_size = current_solution.counts[cluster_index]
        end
    end
    if target_cluster == 0
        target_dispersion = dispersion_context.max_cells_dispersion
        @inbounds for cluster_index in 1:k
            if !is_rejected_for_split_per_cluster[cluster_index] &&
               current_solution.cells_dispersion_per_cluster[cluster_index] > target_dispersion &&
               current_solution.counts[cluster_index] >= min_splittable_size &&
               current_solution.weight_per_cluster[cluster_index] >= min_splittable_weight
                target_cluster = cluster_index
                target_dispersion = current_solution.cells_dispersion_per_cluster[cluster_index]
            end
        end
    end
    if target_cluster == 0
        return 0
    end

    candidate_solution.is_filled = true
    candidate_solution.n_dims = current_solution.n_dims
    candidate_solution.n_points = current_solution.n_points
    candidate_solution.k = k
    @views candidate_solution.assignments[1:current_solution.n_points] .=
        current_solution.assignments[1:current_solution.n_points]
    @views candidate_solution.centers[1:current_solution.n_dims, 1:k] .=
        current_solution.centers[1:current_solution.n_dims, 1:k]
    @views candidate_solution.counts[1:k] .= current_solution.counts[1:k]
    @views candidate_solution.weight_per_cluster[1:k] .= current_solution.weight_per_cluster[1:k]
    @views candidate_solution.cells_dispersion_per_cluster[1:k] .= current_solution.cells_dispersion_per_cluster[1:k]

    split_one_cluster!(
        candidate_solution,
        target_cluster,
        values_of_points,
        weight_per_point,
        kmeans_rounds,
        sizes_buffers,
        dispersion_context,
        rng,
    )

    update_size_statistics!(
        candidate_solution,
        min_cluster_size,
        min_cluster_weight,
        max_cluster_size,
        dispersion_context,
    )
    return target_cluster
end

# Run kmeans (in rounds), overwrite `rerun_solution` with the result, and if the new solution is perfect push a copy
# into `perfect_candidates`. Every K-means run by the Phase 1 walk flows through here, so this is the single point
# that collects perfect candidates for the block's `SolutionCandidates`.
function rerun_kmeans!(
    rerun_solution::KmeansSolution{T},
    perfect_candidates::Vector{SolutionCandidate},
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    initial_centers::Maybe{AbstractMatrix{<:AbstractFloat}},
    new_k::Integer,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Nothing where {T <: AbstractFloat}
    @assert (rerun_solution.n_dims, rerun_solution.n_points) == size(values_of_points)
    # When given centers, kmeans is essentially deterministic from those centers (rng is only used to repick a center
    # if a cluster empties out during iteration), so multiple rounds typically produce the same result.
    rounds = initial_centers === nothing ? kmeans_rounds : 1
    kmeans_result = flame_timed("kmeans_in_rounds") do
        return kmeans_in_rounds(
            values_of_points,
            new_k;
            centers = initial_centers,
            buffers = sizes_buffers.kmeans_buffers,
            rounds,
            rng,
        )
    end
    rerun_solution.is_filled = true
    rerun_solution.k = new_k
    @views rerun_solution.centers[1:rerun_solution.n_dims, 1:new_k] .= kmeans_result.centers
    copyto!(@view(rerun_solution.assignments[1:rerun_solution.n_points]), assignments(kmeans_result))  # NOJET
    copyto!(@view(rerun_solution.counts[1:new_k]), counts(kmeans_result))  # NOJET
    compute_weight_per_cluster!(rerun_solution, weight_per_point)
    compute_all_dispersions!(rerun_solution, sizes_buffers, dispersion_context)
    update_size_statistics!(rerun_solution, min_cluster_size, min_cluster_weight, max_cluster_size, dispersion_context)
    if is_perfect(rerun_solution)
        push!(perfect_candidates, build_candidate_from_solution(rerun_solution))
    end
    return nothing
end

@inline function is_perfect(solution::KmeansSolution)::Bool
    return solution.n_too_tight + solution.n_too_small + solution.n_too_wide + solution.n_too_large == 0
end

# Walk K values across `directions` (in order) starting from `best_solution.k`. For each direction, walks K outward
# while the latest `current_solution` stays "pure" in the same direction (no opposite-direction size issues), updating
# `best_solution` whenever a current solution improves on `penalty`. Pushes a fresh `SolutionCandidate` to
# `perfect_candidates` for every perfect `current_solution` encountered.
function walk_directions!(
    best_solution::KmeansSolution{T},
    current_solution::KmeansSolution{T},
    perfect_candidates::Vector{SolutionCandidate},
    directions::AbstractVector{<:Integer},
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    max_k::Integer,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Nothing where {T <: AbstractFloat}
    start_k = best_solution.k
    for direction in directions
        walk_k = start_k
        while true
            walk_k += direction
            if walk_k < 2 || max_k < walk_k
                break
            end

            rerun_kmeans!(
                current_solution,
                perfect_candidates,
                values_of_points,
                weight_per_point,
                nothing,
                walk_k,
                min_cluster_size,
                min_cluster_weight,
                max_cluster_size,
                kmeans_rounds,
                sizes_buffers,
                dispersion_context,
                rng,
            )

            (cells_range, UMIs_range, dispersion_range) = todox_solution_extrema(current_solution, min_cluster_size)
            todox_ranges_suffix =
                "C: $(cells_range[1])..$(cells_range[2]) " *
                "U: $(round(UMIs_range[1] / 1000; digits=1))K..$(round(UMIs_range[2] / 1000; digits=1))K " *
                "D: $(round(dispersion_range[1]; digits=3))..$(round(dispersion_range[2]; digits=3))"
            if penalty(current_solution) >= penalty(best_solution)
                @warn "TODOX - NEXT K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large) " *
                      todox_ranges_suffix
            else
                copy_solution!(best_solution, current_solution)
                if is_perfect(best_solution)
                    @warn "TODOX - PERFECT K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large) " *
                          todox_ranges_suffix
                else
                    @warn "TODOX - BETTER K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large) " *
                          todox_ranges_suffix
                end
            end
            if !(dispersion_range[1] > 0)
                todox_describe_solution_clusters(
                    current_solution,
                    min_cluster_size,
                    sizes_buffers,
                    dispersion_context,
                    "K $(current_solution.k)",
                )
            end

            # Continue walking the same direction only while the latest result is "pure" - it carries no
            # opposite-direction size issue. The opposite-direction issue means walking further only worsens things.
            current_is_pure_in_direction = if direction > 0
                current_solution.n_too_tight + current_solution.n_too_small == 0
            else
                current_solution.n_too_wide + current_solution.n_too_large == 0
            end
            if !current_is_pure_in_direction
                break
            end
        end
    end
    return nothing
end

# Allocate and return a `SolutionCandidate` whose arrays are sized exactly to `solution`'s active `k` (no aliasing of
# the per-thread max-k buffers). Called by Phase 1 every time it encounters a perfect K so that the candidate survives
# subsequent `rerun_kmeans!` overwrites of the working buffers.
function build_candidate_from_solution(solution::KmeansSolution)::SolutionCandidate
    return SolutionCandidate(
        solution.k,
        copy(@view(solution.counts[1:solution.k])),
        copy(@view(solution.weight_per_cluster[1:solution.k])),
        copy(@view(solution.assignments[1:solution.n_points])),
        is_perfect(solution),
    )
end

# Repeatedly split the largest non-rejected splittable cluster of `best_solution`, accepting the split if it improves
# the penalty (and then running an incremental kmeans seeded from the split centers, accepting the kmeans result while
# it improves further). Terminates naturally once `best_solution` becomes perfect (the next call to
# `split_largest_cluster!` returns 0 because no cluster is too-large or too-wide) or all candidate clusters have been
# rejected.
function walk_split!(
    best_solution::KmeansSolution{T},
    current_solution::KmeansSolution{T},
    candidate_solution::KmeansSolution{T},
    perfect_candidates::Vector{SolutionCandidate},
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    max_k::Integer,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Nothing where {T <: AbstractFloat}
    is_rejected_for_split_per_max_cluster = sizes_buffers.is_rejected_for_split_per_max_cluster
    @views fill!(is_rejected_for_split_per_max_cluster[1:best_solution.k], false)
    while best_solution.k < max_k
        # Split the largest non-rejected splittable cluster (too-large first, then too-wide).
        @views is_rejected_for_split_per_cluster = is_rejected_for_split_per_max_cluster[1:best_solution.k]
        split_cluster_index = split_largest_cluster!(
            candidate_solution,
            best_solution,
            is_rejected_for_split_per_cluster,
            values_of_points,
            weight_per_point,
            min_cluster_size,
            min_cluster_weight,
            max_cluster_size,
            kmeans_rounds,
            sizes_buffers,
            dispersion_context,
            rng,
        )
        if split_cluster_index == 0
            break  # No more candidate clusters to split.
        end
        @warn "TODOX - SPLIT K: $(candidate_solution.k) T: $(candidate_solution.n_too_tight) S: $(candidate_solution.n_too_small) W: $(candidate_solution.n_too_wide) L: $(candidate_solution.n_too_large)"

        if penalty(candidate_solution) < penalty(best_solution)
            @warn "TODOX   SPLIT IS BETTER"
            copy_solution!(best_solution, candidate_solution)
            if is_perfect(best_solution)
                push!(perfect_candidates, build_candidate_from_solution(best_solution))
            end

            # We can re-K-means based on the split clusters to try and improve the solution.
            rerun_kmeans!(
                current_solution,
                perfect_candidates,
                values_of_points,
                weight_per_point,
                @view(
                    candidate_solution.centers[1:candidate_solution.n_dims, 1:candidate_solution.k]
                )::AbstractMatrix{T},
                candidate_solution.k,
                min_cluster_size,
                min_cluster_weight,
                max_cluster_size,
                kmeans_rounds,
                sizes_buffers,
                dispersion_context,
                rng,
            )
            @warn "TODOX - INCR K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large)"

            if penalty(current_solution) < penalty(best_solution)
                @warn "TODOX   INCR IS BETTER"
                copy_solution!(best_solution, current_solution)
                # Cluster indices were reshuffled by K-means; rejected mask no longer corresponds.
                @views fill!(is_rejected_for_split_per_max_cluster[1:best_solution.k], false)
            end
        else
            # The split did not improve over the unsplit best: revert and mark this cluster rejected.
            is_rejected_for_split_per_max_cluster[split_cluster_index] = true
        end
    end
    return nothing
end

# Build an owned `LocalClusters` from `solution`, sharing `dispersion_context.block_cell_indices` and copying the
# active range of the assignment vector. `is_too_small_per_cluster` is derived from the per-cluster counts and weights
# vs. the size/weight thresholds.
function build_local_clusters_from_solution(
    solution::KmeansSolution,
    dispersion_context::DispersionContext,
    min_cluster_size::Real,
    min_cluster_weight::Real,
)::LocalClusters
    @views cluster_sizes = solution.counts[1:solution.k]
    @views cluster_weights = solution.weight_per_cluster[1:solution.k]
    is_too_small_per_cluster = (cluster_sizes .< min_cluster_size) .| (cluster_weights .< min_cluster_weight)
    cluster_index_per_block_cell = Vector{Int}(undef, solution.n_points)
    @views cluster_index_per_block_cell .= solution.assignments[1:solution.n_points]
    return LocalClusters(;
        n_clusters = solution.k,
        block_cell_indices = dispersion_context.block_cell_indices,
        cluster_index_per_block_cell,
        is_too_small_per_cluster,
    )
end

# Phase 1 of the two-phase walker: walk K outward in both directions (as constrained by the initial fit's size
# issues), collecting every perfect K result encountered (via `rerun_kmeans!`'s automatic push to
# `perfect_candidates`). At the end, sort perfects ascending by K so `candidates[1]` is the minimal-K baseline. If no
# perfect was found, the returned `SolutionCandidates` holds a single imperfect compromise (the best-by-penalty
# solution Phase 1 saw).
function kmeans_with_sizes_phase_1(
    values_of_points::AbstractMatrix{T},
    weight_per_point::AbstractVector{<:Real},
    initial_k::Integer;
    max_k::Integer,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::SolutionCandidates where {T <: AbstractFloat}
    @assert max_k >= initial_k
    @assert 2 <= min_cluster_size <= max_cluster_size
    @assert min_cluster_weight > 0
    n_dims, n_points = size(values_of_points)
    @assert length(weight_per_point) == n_points

    best_solution = KmeansSolution{T}(;
        is_filled = false,
        n_dims,
        n_points,
        k = 0,
        n_too_tight = 0,
        n_too_small = 0,
        n_too_wide = 0,
        n_too_large = 0,
        centers = sizes_buffers.max_best_centers,
        assignments = sizes_buffers.max_best_assignments,
        counts = sizes_buffers.max_best_counts,
        weight_per_cluster = sizes_buffers.max_best_weight_per_cluster,
        cells_dispersion_per_cluster = sizes_buffers.max_best_cells_dispersion,
    )
    current_solution = KmeansSolution{T}(;
        is_filled = false,
        n_dims,
        n_points,
        k = 0,
        n_too_tight = 0,
        n_too_small = 0,
        n_too_wide = 0,
        n_too_large = 0,
        centers = sizes_buffers.max_current_centers,
        assignments = sizes_buffers.max_current_assignments,
        counts = sizes_buffers.max_current_counts,
        weight_per_cluster = sizes_buffers.max_current_weight_per_cluster,
        cells_dispersion_per_cluster = sizes_buffers.max_current_cells_dispersion,
    )
    candidate_solution = KmeansSolution{T}(;
        is_filled = false,
        n_dims,
        n_points,
        k = 0,
        n_too_tight = 0,
        n_too_small = 0,
        n_too_wide = 0,
        n_too_large = 0,
        centers = sizes_buffers.max_candidate_centers,
        assignments = sizes_buffers.max_candidate_assignments,
        counts = sizes_buffers.max_candidate_counts,
        weight_per_cluster = sizes_buffers.max_candidate_weight_per_cluster,
        cells_dispersion_per_cluster = sizes_buffers.max_candidate_cells_dispersion,
    )

    perfect_candidates = SolutionCandidate[]

    rerun_kmeans!(
        best_solution,
        perfect_candidates,
        values_of_points,
        weight_per_point,
        nothing,
        initial_k,
        min_cluster_size,
        min_cluster_weight,
        max_cluster_size,
        kmeans_rounds,
        sizes_buffers,
        dispersion_context,
        rng,
    )

    @warn "TODOX ### KMEANS K: $(best_solution.k) T: $(best_solution.n_too_tight) S: $(best_solution.n_too_small) W: $(best_solution.n_too_wide) L: $(best_solution.n_too_large)"
    let (cells_range, UMIs_range, dispersion_range) = todox_solution_extrema(best_solution, min_cluster_size)
        initial_label = is_perfect(best_solution) ? "PERFECT" : "INITIAL"
        @warn "TODOX - $initial_label K: $(best_solution.k) " *
              "C: $(cells_range[1])..$(cells_range[2]) " *
              "U: $(round(UMIs_range[1] / 1000; digits=1))K..$(round(UMIs_range[2] / 1000; digits=1))K " *
              "D: $(round(dispersion_range[1]; digits=3))..$(round(dispersion_range[2]; digits=3))"
        if !(dispersion_range[1] > 0)
            todox_describe_solution_clusters(
                best_solution,
                min_cluster_size,
                sizes_buffers,
                dispersion_context,
                "$initial_label K $(best_solution.k)",
            )
        end
    end

    has_too_big = best_solution.n_too_wide + best_solution.n_too_large > 0
    has_too_tiny = best_solution.n_too_tight + best_solution.n_too_small > 0
    # No size issues -> both directions; only one kind of issue -> walk only the direction that fixes it; both at once
    # -> the direction walk is futile, leave it to the split loop.
    directions = if !has_too_big && !has_too_tiny
        [-1, +1]
    elseif has_too_big && !has_too_tiny
        [+1]
    elseif has_too_tiny && !has_too_big
        [-1]
    else
        Int[]
    end

    if !isempty(directions)
        walk_directions!(
            best_solution,
            current_solution,
            perfect_candidates,
            directions,
            values_of_points,
            weight_per_point,
            max_k,
            min_cluster_size,
            min_cluster_weight,
            max_cluster_size,
            kmeans_rounds,
            sizes_buffers,
            dispersion_context,
            rng,
        )
    end

    @warn "TODOX = INTER K: $(best_solution.k) T: $(best_solution.n_too_tight) S: $(best_solution.n_too_small) W: $(best_solution.n_too_wide) L: $(best_solution.n_too_large)"

    walk_split!(
        best_solution,
        current_solution,
        candidate_solution,
        perfect_candidates,
        values_of_points,
        weight_per_point,
        max_k,
        min_cluster_size,
        min_cluster_weight,
        max_cluster_size,
        kmeans_rounds,
        sizes_buffers,
        dispersion_context,
        rng,
    )

    @assert best_solution.is_filled
    @warn "TODOX - PHASE 1 BEST K: $(best_solution.k) T: $(best_solution.n_too_tight) S: $(best_solution.n_too_small) W: $(best_solution.n_too_wide) L: $(best_solution.n_too_large) PERFECT CANDIDATES: $(length(perfect_candidates))"

    if isempty(perfect_candidates)
        # No perfect found - take the best-by-penalty imperfect compromise as the lone candidate.
        return SolutionCandidates(n_points, [build_candidate_from_solution(best_solution)])
    end
    # Sort ascending by K so `candidates[1]` is the minimal-K baseline.
    sort!(perfect_candidates; by = candidate -> candidate.k)
    return SolutionCandidates(n_points, perfect_candidates)
end


# Build a `LocalClusters` from a `SolutionCandidate`. `block_cell_indices` is shared (not copied); the assignment
# vector and `is_too_small_per_cluster` are owned by the returned `LocalClusters`.
function build_local_clusters_from_candidate(
    candidate::SolutionCandidate,
    block_cell_indices::AbstractVector{<:Integer},
    min_cluster_size::Real,
    min_cluster_weight::Real,
)::LocalClusters
    is_too_small_per_cluster =
        (candidate.counts .< min_cluster_size) .| (candidate.weight_per_cluster .< min_cluster_weight)
    return LocalClusters(;
        n_clusters = candidate.k,
        block_cell_indices,
        cluster_index_per_block_cell = copy(candidate.assignments),
        is_too_small_per_cluster,
    )
end


end  # module
