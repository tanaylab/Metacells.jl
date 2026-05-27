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

using ..AnalyzeModules
using ..Contracts

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
    )
end

struct KmeansSizesResult
    assignments::AbstractVector{Int}
    counts::AbstractVector{Int}
    weight_per_cluster::AbstractVector{Float64}
end

Clustering.assignments(r::KmeansSizesResult) = r.assignments
Clustering.counts(r::KmeansSizesResult) = r.counts
Clustering.nclusters(r::KmeansSizesResult) = length(r.counts)

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
) function sharpen_metacells!(;
    sharp_daf::DafWriter,
    base_daf::DafReader,
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
                        max(n_neighborhood_cells_per_block[b], n_cells_per_block[b]) / mean_metacell_cells_per_block[b],
                    ),
                ),
                1,
            ) for b in 1:n_blocks
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
    post_migration_n_cells_per_block = [sum(block_index_per_cell .== b) for b in 1:n_blocks]
    post_max_n_kmeans_points = maximum(post_migration_n_cells_per_block)
    post_max_n_kmeans_clusters =
        2 * maximum(
            max(Int(round(post_migration_n_cells_per_block[b] / mean_metacell_cells_per_block[b])), 1) for
            b in 1:n_blocks
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
        max(Int(round(n_neighborhood_cells_per_block[b] / mean_metacell_cells_per_block[b])), 1) for b in 1:n_blocks
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

    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        policy = :static,
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
    kmeans_sizes_max_buffers_per_thread::AbstractVector{<:KmeansSizesBuffers},
    rng::AbstractRNG,
)::AbstractVector{Maybe{LocalClusters}}
    n_cells = length(total_UMIs_per_cell)
    n_blocks = size(is_found_per_module_per_block, 2)
    n_modules = axis_length(base_daf, "module")

    local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)

    n_cells_per_block = [sum(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    max_n_block_cells = maximum(n_cells_per_block)

    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)

    n_genes = size(module_index_per_gene_per_block, 1)

    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]
    is_in_block_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:maxthreadid()]
    z_score_per_max_module_per_max_block_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_block_cells) for _ in 1:maxthreadid()]
    cluster_index_per_block_cell_per_block =
        [Vector{Int}(undef, n_cells_per_block[block_index]) for block_index in 1:n_blocks]

    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        policy = :static,
        name = "compute_local_clusters",
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "local_clusters"),
    ) do block_index, rng
        is_in_block_per_cell = is_in_block_per_cell_per_thread[threadid()]
        is_in_block_per_cell .= block_index_per_cell .== block_index
        block_cell_indices = findall(is_in_block_per_cell)
        n_block_cells = n_cells_per_block[block_index]
        @assert n_block_cells == length(block_cell_indices)

        if n_block_cells == 0
            local_clusters_per_block[block_index] = nothing
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
        cluster_index_per_block_cell = cluster_index_per_block_cell_per_block[block_index]
        sizes_buffers = kmeans_sizes_max_buffers_per_thread[threadid()]

        @views weight_per_block_cell = total_UMIs_per_cell[block_cell_indices]
        if n_block_clusters == 1
            fill!(cluster_index_per_block_cell, 1)
            sizes_buffers.max_best_counts[1] = n_block_cells
            sizes_buffers.max_best_weight_per_cluster[1] = sum(weight_per_block_cell)
            cluster_sizes = @view(sizes_buffers.max_best_counts[1:1])
            cluster_weights = @view(sizes_buffers.max_best_weight_per_cluster[1:1])
        else
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
            kmeans_result = flame_timed("kmeans_with_sizes") do
                return kmeans_with_sizes(
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
                    best_assignments = cluster_index_per_block_cell,
                    rng,
                )
            end

            cluster_sizes = counts(kmeans_result)
            cluster_weights = kmeans_result.weight_per_cluster
        end
        local_clusters_per_block[block_index] = LocalClusters(;
            block_cell_indices,
            cluster_index_per_block_cell,
            is_too_small_per_cluster = (cluster_sizes .< min_cells_in_metacell) .|
                                       (cluster_weights .< min_metacell_total_UMIs),
        )

        return nothing
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
            n_total_new_metacells += length(local_clusters.is_too_small_per_cluster)
            @foreach_true_index local_clusters.is_too_small_per_cluster cluster_index begin  # NOLINT
                n_new_outlier_cells += sum(local_clusters.cluster_index_per_block_cell .== cluster_index)  # NOLINT
            end
        end
    end

    cells_of_new_metacells = Vector{AbstractVector{<:Integer}}()
    block_name_per_new_metacell = Vector{AbstractString}()

    for (block_index, local_clusters) in enumerate(local_clusters_per_block)
        block_name = name_per_block[block_index]
        if local_clusters !== nothing
            n_clusters = length(local_clusters.is_too_small_per_cluster)
            for cluster_index in 1:n_clusters
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

# Recompute n_too_tight (clusters whose maximal cells_dispersion is below dispersion_context.min_cells_dispersion),
# n_too_small (clusters whose count is below min_cluster_size OR whose weight is below min_cluster_weight),
# n_too_wide (clusters whose maximal cells_dispersion is above dispersion_context.max_cells_dispersion), and
# n_too_large (clusters whose count is above max_cluster_size). Assumes solution.cells_dispersion_per_cluster is
# already populated for the active clusters.
function update_size_statistics!(
    solution::KmeansSolution,
    min_cluster_size::Real,
    min_cluster_weight::Real,
    max_cluster_size::Real,
    dispersion_context::DispersionContext,
)::Nothing
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
    return (solution.n_too_tight + solution.n_too_small, solution.n_too_wide + solution.n_too_large, solution.k)
end

# Split a single cluster into two using a 2-cluster kmeans on its points. Replaces the old cluster's center with one
# sub-center, appends the other at the new last position (incrementing `candidate_solution.k`), reassigns the points,
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
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Integer where {T <: AbstractFloat}
    k = current_solution.k

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

# Run kmeans (in rounds) and overwrite the solution with the result.
function rerun_kmeans!(
    rerun_solution::KmeansSolution{T},
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
    return nothing
end

# Find a kmeans clustering where each cluster has count in [min_cluster_size, max_cluster_size] and total weight
# (sum of `weight_per_point` across the cluster's points) at least `min_cluster_weight`. Walks K from `initial_k`
# (incrementing to fix too-large clusters, decrementing to fix too-small) up to `max_k` / down to 2, tracking the
# best solution by `n_too_tight + n_too_small + n_too_wide + n_too_large` and stopping once K crosses out of its
# starting direction. If the resulting best still has too-large clusters, repeatedly splits them and runs an incremental
# kmeans seeded from the split centers, accepting the kmeans result while it improves and otherwise reverting to the
# split-only solution. Any too-small clusters that remain are left for the caller to dissolve as outliers.
function kmeans_with_sizes(
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
    best_assignments::AbstractVector{Int},
    rng::AbstractRNG,
)::KmeansSizesResult where {T <: AbstractFloat}
    @assert max_k >= initial_k
    @assert 2 <= min_cluster_size <= max_cluster_size
    @assert min_cluster_weight > 0
    n_dims, n_points = size(values_of_points)
    @assert length(weight_per_point) == n_points

    n_full_kmeans = 0
    n_incremental_kmeans = 0

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

    n_full_kmeans += 1
    rerun_kmeans!(
        best_solution,
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
    @warn "TODOX - PERFECT K: $(best_solution.k)"

    has_too_big = best_solution.n_too_wide + best_solution.n_too_large > 0
    has_too_tiny = best_solution.n_too_tight + best_solution.n_too_small > 0
    have_perfect = !has_too_big && !has_too_tiny
    if has_too_big != has_too_tiny || have_perfect
        # Modify K to avoid only too-small or only too-large clusters, or to have less clusters.
        direction = has_too_big ? +1 : -1
        walk_k = best_solution.k
        while true
            walk_k += direction
            if walk_k < 2 || max_k < walk_k
                break
            end

            n_full_kmeans += 1
            rerun_kmeans!(
                current_solution,
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

            if penalty(current_solution) >= penalty(best_solution)
                @warn "TODOX - NEXT K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large)"
            else
                copy_solution!(best_solution, current_solution)
                if best_solution.n_too_tight +
                   best_solution.n_too_small +
                   best_solution.n_too_wide +
                   best_solution.n_too_large > 0
                    @warn "TODOX - BETTER K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large)"
                else
                    if have_perfect
                        @warn "TODOX - REPERFECT K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large)"
                    else
                        @warn "TODOX - PERFECT K: $(current_solution.k) T: $(current_solution.n_too_tight) S: $(current_solution.n_too_small) W: $(current_solution.n_too_wide) L: $(current_solution.n_too_large)"
                    end
                    have_perfect = true
                    if direction > 0
                        break
                    end
                end
            end

            current_is_pure_in_same_direction = if direction > 0
                current_solution.n_too_tight + current_solution.n_too_small == 0
            else
                current_solution.n_too_wide + current_solution.n_too_large == 0
            end
            if !current_is_pure_in_same_direction
                break
            end
        end
    end

    @warn "TODOX = INTER K: $(best_solution.k) T: $(best_solution.n_too_tight) S: $(best_solution.n_too_small) W: $(best_solution.n_too_wide) L: $(best_solution.n_too_large)"

    best_is_split = false
    n_split_kmeans = 0
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
        n_split_kmeans += 1
        @warn "TODOX - SPLIT K: $(candidate_solution.k) T: $(candidate_solution.n_too_tight) S: $(candidate_solution.n_too_small) W: $(candidate_solution.n_too_wide) L: $(candidate_solution.n_too_large)"

        if penalty(candidate_solution) < penalty(best_solution)
            @warn "TODOX   SPLIT IS BETTER"
            copy_solution!(best_solution, candidate_solution)
            best_is_split = true

            # We can re-K-means based on the split clusters to try and improve the solution.
            n_incremental_kmeans += 1
            rerun_kmeans!(
                current_solution,
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
                best_is_split = false
                # Cluster indices were reshuffled by K-means; rejected mask no longer corresponds.
                @views fill!(is_rejected_for_split_per_max_cluster[1:best_solution.k], false)
            end
        else
            # The split did not improve over the unsplit best: revert and mark this cluster rejected.
            is_rejected_for_split_per_max_cluster[split_cluster_index] = true
        end
    end

    @assert best_solution.is_filled
    @warn "TODOX - BEST K: $(best_solution.k) T: $(best_solution.n_too_tight) S: $(best_solution.n_too_small) W: $(best_solution.n_too_wide) L: $(best_solution.n_too_large) SPLIT?: $(best_is_split) F: $(n_full_kmeans) I: $(n_incremental_kmeans) SP: $(n_split_kmeans)"

    @views best_assignments[1:n_points] .= best_solution.assignments[1:n_points]
    return KmeansSizesResult(
        best_assignments,
        @view(best_solution.counts[1:best_solution.k]),
        @view(best_solution.weight_per_cluster[1:best_solution.k]),
    )
end

end  # module
