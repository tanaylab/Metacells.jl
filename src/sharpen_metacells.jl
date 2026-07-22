"""
Compute "better" metacells based on the local linear model approximating the manifold.
"""
module SharpenMetacells

export compute_matrix_of_n_cells_per_prev_block_per_block!
export compute_matrix_of_n_cells_per_prev_block_type_per_block_type!
export compute_matrix_of_n_cells_per_prev_metacell_type_per_metacell_type!
export compute_vector_of_global_flow_order_per_type!
export sharpen_metacells!

using Base.Threads
using Clustering
using DataAxesFormats
using LoopVectorization
using Random
using SparseArrays
using StatsBase
using TanayLabUtilities

using ..AnalyzeBlocks
using ..AnalyzeModules
using ..Contracts
using ..Defaults

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
import Metacells.Contracts.matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_n_cells_per_prev_block_per_block
import Metacells.Contracts.matrix_of_n_cells_per_prev_block_type_per_block_type
import Metacells.Contracts.matrix_of_n_cells_per_prev_metacell_type_per_metacell_type
import Metacells.Contracts.matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.prev_block_axis
import Metacells.Contracts.type_axis
import Metacells.Contracts.vector_of_base_block_per_metacell
import Metacells.Contracts.vector_of_block_closest_by_pertinent_markers_per_cell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_global_flow_order_per_type
import Metacells.Contracts.vector_of_is_base_outlier_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_block
import Metacells.Contracts.vector_of_n_metacells_per_block
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_outlier_actual_UMIs_per_cell
import Metacells.Contracts.vector_of_outlier_by_prev_block_per_cell
import Metacells.Contracts.vector_of_outlier_by_prev_module_per_cell
import Metacells.Contracts.vector_of_outlier_expected_UMIs_per_cell
import Metacells.Contracts.vector_of_outlier_in_metacell_per_cell
import Metacells.Contracts.vector_of_outlier_in_prev_block_per_cell
import Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_type_per_block
import Metacells.Contracts.vector_of_type_per_metacell

struct KmeansSizesBuffers{T <: AbstractFloat}
    max_best_assignments::Vector{Int}
    max_best_counts::Vector{Int}
    max_best_weight_per_cluster::Vector{Float64}
    max_best_centers::Matrix{T}
    max_current_assignments::Vector{Int}
    max_current_counts::Vector{Int}
    max_current_weight_per_cluster::Vector{Float64}
    max_current_centers::Matrix{T}
    max_candidate_assignments::Vector{Int}
    max_candidate_counts::Vector{Int}
    max_candidate_weight_per_cluster::Vector{Float64}
    max_candidate_centers::Matrix{T}
    max_split_point_indices::Vector{Int}
    max_split_points::Matrix{T}
    is_rejected_for_split_per_max_cluster::BitVector
    max_best_cells_dispersion::Vector{Float32}
    max_current_cells_dispersion::Vector{Float32}
    max_candidate_cells_dispersion::Vector{Float32}
end

struct DispersionScratches
    cluster_cell_indices::Vector{Int}
    cells_dispersion_per_module::Vector{Float32}
    mean_normalized_per_module::Vector{Float64}
    total_UMIs_per_cluster_cell::Vector{Float32}
    normalized_factor_per_cluster_cell::Vector{Float32}
end

function DispersionScratches(; n_points::Integer, n_modules::Integer)
    return DispersionScratches(
        Vector{Int}(undef, n_points),
        Vector{Float32}(undef, n_modules),
        Vector{Float64}(undef, n_modules),
        Vector{Float32}(undef, n_points),
        Vector{Float32}(undef, n_points),
    )
end

function KmeansSizesBuffers{T}(; n_dims::Integer, max_k::Integer, n_points::Integer) where {T <: AbstractFloat}
    return KmeansSizesBuffers{T}(
        Vector{Int}(undef, n_points),         # max_best_assignments
        Vector{Int}(undef, max_k),            # max_best_counts
        Vector{Float64}(undef, max_k),        # max_best_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_best_centers
        Vector{Int}(undef, n_points),         # max_current_assignments
        Vector{Int}(undef, max_k),            # max_current_counts
        Vector{Float64}(undef, max_k),        # max_current_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_current_centers
        Vector{Int}(undef, n_points),         # max_candidate_assignments
        Vector{Int}(undef, max_k),            # max_candidate_counts
        Vector{Float64}(undef, max_k),        # max_candidate_weight_per_cluster
        Matrix{T}(undef, n_dims, max_k),      # max_candidate_centers
        Vector{Int}(undef, n_points),         # max_split_point_indices
        Matrix{T}(undef, n_dims, n_points),   # max_split_points
        BitVector(undef, max_k),              # is_rejected_for_split_per_max_cluster
        Vector{Float32}(undef, max_k),        # max_best_cells_dispersion
        Vector{Float32}(undef, max_k),        # max_current_cells_dispersion
        Vector{Float32}(undef, max_k),        # max_candidate_cells_dispersion
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
end

@kwdef struct DispersionContext
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer}
    total_UMIs_per_cell::AbstractVector{<:Integer}
    block_cell_indices::AbstractVector{<:Integer}
    is_found_per_module::BitVector
    gene_indices_per_module::Vector{Vector{Int}}
    normalized_UMIs_quantile::AbstractFloat
    min_module_UMIs::Integer
    min_cells_dispersion::AbstractFloat
    max_cells_dispersion::AbstractFloat
    min_cluster_size::Integer
end

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
    gene_index_per_friend::AbstractVector{<:Integer}
    friend_position_per_series_per_base_block::Dict{Int, Vector{Int}}
    block_cell_position_per_point_per_base_block::Dict{Int, Vector{Int}}
    # Per-thread scratch (sized to the walkable blocks' maxima).
    total_UMIs_per_max_cluster::Vector{Float64}
    UMIs_per_friend_per_max_cluster::Matrix{Float64}
    variable_per_block_cell_per_friend::Matrix{Float32}
    is_active_per_block_cell::Vector{Bool}
    punctuated_total_per_block_cell::Vector{Float64}
end

struct SolutionCandidate
    k::Int
    counts::Vector{Int}
    weight_per_cluster::Vector{Float64}
    assignments::Vector{Int}
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
        max_cells_in_metacell::Integer = $(DEFAULT.max_cells_in_metacell),
        target_cell_total_UMIs_quantile::AbstractFloat = $(DEFAULT.target_cell_total_UMIs_quantile),
        min_metacell_total_UMIs::Integer = $(DEFAULT.min_metacell_total_UMIs),
        min_migration_likelihood::AbstractFloat = $(DEFAULT.min_migration_likelihood),
        max_cells_dispersion_in_metacell::AbstractFloat = $(DEFAULT.max_cells_dispersion_in_metacell),
        min_cells_dispersion_in_metacell::AbstractFloat = $(DEFAULT.min_cells_dispersion_in_metacell),
        normalized_UMIs_quantile::AbstractFloat = $(DEFAULT.normalized_UMIs_quantile),
        min_module_UMIs::Integer = $(DEFAULT.min_module_UMIs),
        kmeans_rounds::Integer = $(DEFAULT.kmeans_rounds),
        sharpening_round::Integer,
        improvement_half_life::Maybe{Integer} = $(DEFAULT.improvement_half_life),
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        std_UMIs_regularization::AbstractFloat = $(DEFAULT.std_UMIs_regularization),
        min_outlier_fold::AbstractFloat = $(DEFAULT.min_outlier_fold),
        outlier_UMIs_regularization::AbstractFloat = $(DEFAULT.outlier_UMIs_regularization),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given an `base_daf` metacells repository with a blocks structure and local gene modules that describe the cell state
manifold, compute a `sharp_daf` metacells repository, which hopefully more faithfully captures this manifold.

 1. We cluster using K-means all the cells in each neighborhood using the z-score of the expression of the modules of
    that neighborhood. The z-score compares the module's linear fraction in the cell to the mean and the standard
    deviation of the module's linear fraction in the cells of the block's environment. It is regularized - we move the
    difference from the mean away from zero by `std_fraction_regularization`, and increase the standard deviation by
    it - so that modules with a tiny standard deviation will not turn single-UMI noise into a large z-score. This
    regularization is specified in UMIs (`std_UMIs_regularization`) and is converted to a linear fraction by dividing
    it by the mean total UMIs of a cell in the block's neighborhood. The number of clusters is the number of base
    metacells in that neighborhood.
 2. We assign to each cluster the base block which is most frequent in the cluster cells.
    Cells of the base base block of the neighborhood, which belong to a cluster that is assigned to a different
    block, and which also belong to a cluster of that block in the neighborhood of that block, are migrated to that
    block, but only if the enrichment of the cells of that block in the cluster is at least `min_migration_likelihood`
    times what would be expected assuming random clustering based on the relative sizes of the base and other blocks.
 3. Having finalized the block to which each cell belongs to, we cluster all the cells in each block using K-means
    using the modules of the neighborhood of that block. We start with the expected number of metacells in that block
    (based on the mean number of cells per metacell in the base block) and adjust the number of clusters to try and
    enforce the sizes of the clusters - not more than twice that mean, no more than `max_cells_in_metacell`, and no less
    than `min_cells_in_metacell`. A cluster is also considered too-large if it is larger than `max_cells_in_metacell`,
    if its maximal `cells_dispersion` (across the block's modules) is above `max_cells_dispersion_in_metacell`, and
    too-small if its maximal `cells_dispersion` is below `min_cells_dispersion_in_metacell`. In edge cases we dissolve
    too-small clusters, relocating each of their cells to the nearest surviving cluster, by the distance to the cluster's
    center in the same module z-score space K-means clustered them in. A block's last cluster is never dissolved, so
    every clustered cell always ends up in some cluster. We increase the number of target metacells for every input
    metacell whose maximal `cells_dispersion` is above `max_cells_dispersion_in_metacell`. This still leaves us with
    multiple possible numbers of clusters for the block; we break the tie using a prediction of the mean correlation
    with most-correlated marker genes in the (previous round's) neighborhoods. To stabilize the choice as sharpening
    proceeds, a candidate whose number of clusters differs from the expected number must improve this prediction by an
    extra `cooldown_margin` for each cluster of difference. The `sharpening_round` is 1-based: the first sharpening
    round (`sharpening_round == 1`) applies no cooldown (the plain correlation tie-break), and for each later round the
    margin's distance from its maximum halves every `improvement_half_life` rounds. Passing
    `improvement_half_life == nothing` disables the cooldown entirely.
 4. We detect the outlier cells and eject them from their metacell. For each cell in each metacell, and for every
    found module instance in the previous round's repository (across all blocks, not just the modules of the cell's own region), we
    compare the actual UMIs of the module's genes in the cell to the expected UMIs. The expectation is punctuated - the
    module's linear fraction in the cell's metacell excluding the cell itself, times the cell's total UMIs - so a strong
    outlier cannot inflate its own baseline. The comparison is a fold factor `log2((actual + outlier_UMIs_regularization) / (expected + outlier_UMIs_regularization))`, using a strong UMIs regularization so low-count modules do not trigger.
    A cell whose maximal fold (across all module instances) exceeds `min_outlier_fold` is an outlier: it loses its
    metacell, and we record the certificate (the metacell it was ejected from, the previous round's block and module that
    flagged it, and the expected and actual UMIs) of the module instance that most deviates. Ejecting the outliers can
    drop a metacell below `min_cells_in_metacell`; we dissolve such a metacell, relocate all its cells to the nearest
    surviving metacell of its block (by the same module z-score distance used above), and detect the outliers again - so
    a cell forced into a metacell it does not fit is caught. This repeats until no metacell is dissolved. A block's last
    metacell is never dissolved, so every outlier's certificate names a metacell that survives.
 5. The final clusters are the sharp metacells. We name them using the `prefix`, the convention is to advance the letter
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
    data = [
        vector_of_metacell_per_cell(CreatedOutput),
        vector_of_base_block_per_metacell(CreatedOutput),
        vector_of_outlier_in_prev_block_per_cell(CreatedOutput),
        vector_of_outlier_in_metacell_per_cell(CreatedOutput),
        vector_of_outlier_by_prev_block_per_cell(CreatedOutput),
        vector_of_outlier_by_prev_module_per_cell(CreatedOutput),
        vector_of_outlier_expected_UMIs_per_cell(CreatedOutput),
        vector_of_outlier_actual_UMIs_per_cell(CreatedOutput),
    ],
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
        vector_of_is_excluded_per_cell(RequiredInput),
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
        vector_of_total_neighborhood_UMIs_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
    ],
) function sharpen_metacells!(;
    sharp_daf::DafWriter,
    base_daf::DafReader,
    prefix::AbstractString = "M",
    min_cells_in_metacell::Integer = 12,
    max_cells_in_metacell::Integer = 800,
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
    sharpening_round::Integer,
    improvement_half_life::Maybe{Integer} = 2,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    std_UMIs_regularization::AbstractFloat = 2.0,
    min_outlier_fold::AbstractFloat = 3.0,
    outlier_UMIs_regularization::AbstractFloat = 2.0,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert min_cells_in_metacell >= 0
    @assert max_cells_in_metacell >= 2 * min_cells_in_metacell
    @assert 0 < target_cell_total_UMIs_quantile < 1
    @assert min_metacell_total_UMIs > 0
    @assert min_migration_likelihood > 0
    @assert 0 <= min_cells_dispersion_in_metacell <= max_cells_dispersion_in_metacell
    @assert 0 <= normalized_UMIs_quantile <= 1
    @assert min_module_UMIs >= 0
    @assert kmeans_rounds > 0
    @assert sharpening_round >= 1
    @assert improvement_half_life === nothing || improvement_half_life > 0
    @assert std_UMIs_regularization > 0
    @assert min_outlier_fold > 0
    @assert outlier_UMIs_regularization > 0

    n_cells = axis_length(base_daf, "cell")
    n_blocks = axis_length(base_daf, "block")
    n_base_metacells = axis_length(base_daf, "metacell")
    name_per_block = axis_vector(base_daf, "block")

    UMIs_per_cell_per_gene = get_matrix(base_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = mutable_array(densify(get_vector(base_daf, "cell", "total_UMIs").array))

    mean_linear_fraction_in_environment_cells_per_module_per_block =
        get_matrix(base_daf, "module", "block", "mean_linear_fraction_in_environment_cells").array
    std_linear_fraction_in_environment_cells_per_module_per_block =
        get_matrix(base_daf, "module", "block", "std_linear_fraction_in_environment_cells").array

    # The z-score regularization is specified in UMIs; convert it to a linear fraction using the mean total UMIs of a
    # cell in each block's neighborhood.
    mean_total_UMIs_in_neighborhood_cells_per_block =
        get_vector(base_daf, "block", "total_neighborhood_UMIs").array ./
        get_vector(base_daf, "block", "n_neighborhood_cells").array
    std_fraction_regularization_per_block =
        Float32.(std_UMIs_regularization ./ mean_total_UMIs_in_neighborhood_cells_per_block)

    is_in_neighborhood_per_other_block_per_base_block =
        get_matrix(base_daf, "block", "block", "is_in_neighborhood").array
    is_found_per_module_per_block = get_matrix(base_daf, "module", "block", "is_found").array
    module_index_per_gene_per_block = base_daf["@ gene @ block :: module ?? 0 : index"].array

    block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array
    is_base_outlier_per_cell = get_vector(base_daf, "cell", "is_base_outlier").array
    is_excluded_per_cell = get_vector(base_daf, "cell", "is_excluded").array
    closest_block_index_per_cell = base_daf["@ cell : block.closest_by_pertinent_markers : index"].array
    # Base outliers and excluded cells are never clustered. Every other cell gets a block: grouped cells use their
    # metacell's block; round outliers (lost their metacell in a previous sharpening round) use the closest block by
    # pertinent markers as a synthetic base.
    block_index_per_cell = ifelse.(
        is_base_outlier_per_cell .| is_excluded_per_cell,
        zero(eltype(block_index_per_cell)),
        ifelse.(block_index_per_cell .> 0, block_index_per_cell, closest_block_index_per_cell),
    )
    n_cells_per_block = [count(==(block_index), block_index_per_cell) for block_index in 1:n_blocks]

    # Floor on metacell size (in cells) needed for each block to clear `min_metacell_total_UMIs`, assuming each cell
    # contributes at the `target_cell_total_UMIs_quantile` of the block's per-cell UMI distribution.
    target_min_cells_per_block = Vector{Int}(undef, n_blocks)
    parallel_loop_wo_rng(
        1:n_blocks;
        name = "target_min_cells_per_block",
        weights = n_cells_per_block,
        progress = DebugProgress(sum(n_cells_per_block); group = :mcs_loops, desc = "target_min_cells_per_block"),
    ) do block_index
        @views block_total_UMIs_per_cell = total_UMIs_per_cell[block_index_per_cell .== block_index]
        if isempty(block_total_UMIs_per_cell)
            target_min_cells_per_block[block_index] = min_cells_in_metacell
        else
            target_total_UMIs_per_cell = quantile(block_total_UMIs_per_cell, target_cell_total_UMIs_quantile)  # NOLINT
            n_cells_for_min_UMIs = Int(ceil(min_metacell_total_UMIs / target_total_UMIs_per_cell))
            target_min_cells_per_block[block_index] =
                min(max(n_cells_for_min_UMIs, min_cells_in_metacell), div(max_cells_in_metacell, 2))
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
                        max(n_neighborhood_cells_per_block[block_index], n_cells_per_block[block_index]) /
                        mean_metacell_cells_per_block[block_index],
                    ),
                ),
                1,
            ) for block_index in 1:n_blocks
        )
    # Channel pool of `KMeansBuffers` for `kmeans_in_rounds`'s nested parallel rounds, shared across
    # `compute_preferred_block_index_per_cell_per_block` and `compute_local_clusters` (which run sequentially below).
    # One buffer per thread is sufficient because nested parallel rounds are bounded by Julia's thread pool.
    kmeans_buffer_pool = Channel{KMeansBuffers{Float32}}(maxthreadid())
    for _ in 1:maxthreadid()
        put!(
            kmeans_buffer_pool,
            KMeansBuffers{Float32}(;
                n_dims = max_n_block_modules,
                max_k = max_n_kmeans_clusters,
                n_points = max_n_kmeans_points,
            ),
        )
    end

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
        mean_linear_fraction_in_environment_cells_per_module_per_block,
        std_linear_fraction_in_environment_cells_per_module_per_block,
        std_fraction_regularization_per_block,
        min_migration_likelihood,
        kmeans_buffer_pool,
        rng,
    )

    block_index_per_cell =
        compute_preferred_block_index_of_cells(; block_index_per_cell, preferred_block_index_per_cell_per_block)

    # Post-migration blocks can have more cells than pre-migration; grow the pre-migration `kmeans_buffer_pool` if
    # needed (`kmeans_sizes_max_buffers_per_thread` is allocated once below, after the maxes are final).
    post_migration_n_cells_per_block = [count(==(block_index), block_index_per_cell) for block_index in 1:n_blocks]
    post_max_n_kmeans_points = maximum(post_migration_n_cells_per_block)
    post_max_n_kmeans_clusters =
        2 * maximum(
            max(
                Int(round(post_migration_n_cells_per_block[block_index] / mean_metacell_cells_per_block[block_index])),
                1,
            ) for block_index in 1:n_blocks
        )
    if post_max_n_kmeans_points > max_n_kmeans_points || post_max_n_kmeans_clusters > max_n_kmeans_clusters
        max_n_kmeans_points = max(max_n_kmeans_points, post_max_n_kmeans_points)
        max_n_kmeans_clusters = max(max_n_kmeans_clusters, post_max_n_kmeans_clusters)
        # The pool holds references to the old (smaller) buffers. Rebuild to match the new maxes. Only fires when at
        # least one of the maxes actually grew (gated by the surrounding `if`).
        close(kmeans_buffer_pool)
        kmeans_buffer_pool = Channel{KMeansBuffers{Float32}}(maxthreadid())
        for _ in 1:maxthreadid()
            put!(
                kmeans_buffer_pool,
                KMeansBuffers{Float32}(;
                    n_dims = max_n_block_modules,
                    max_k = max_n_kmeans_clusters,
                    n_points = max_n_kmeans_points,
                ),
            )
        end
    end

    # Sized to the final maxes; consumed only by `compute_local_clusters` below.
    kmeans_sizes_max_buffers_per_thread = [
        KmeansSizesBuffers{Float32}(;
            n_dims = max_n_block_modules,
            max_k = max_n_kmeans_clusters,
            n_points = max_n_kmeans_points,
        ) for _ in 1:maxthreadid()
    ]

    # The tie-break evaluates against the previous round (base_daf), so its blocks are the base blocks.
    base_block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array
    is_in_base_neighborhood_per_other_base_block_per_base_block = is_in_neighborhood_per_other_block_per_base_block
    most_correlated_gene_in_base_neighborhood_per_gene_per_base_block =
        get_matrix(base_daf, "gene", "block", "most_correlated_gene_in_neighborhood").array
    gene_name_to_index = axis_dict(base_daf, "gene")
    n_base_blocks = n_blocks

    # Per base block, the (gene, friend) correlation series as raw (global) gene indices - parallel vectors. The friend
    # is the gene's most-correlated-in-base-neighborhood partner; every gene with such a partner contributes a series.
    # Consumed by `build_grouped_for_base_block` (for reading `UMIs_per_cell_per_gene` columns) and by
    # `precompute_walkable_indirection` (to build each walkable block's friend-gene subspace).
    correlated_gene_indices_per_base_block = Vector{Vector{Int}}(undef, n_base_blocks)
    friend_gene_indices_per_base_block = Vector{Vector{Int}}(undef, n_base_blocks)
    for base_block_index in 1:n_base_blocks
        @views most_correlated_gene_in_base_neighborhood_per_gene =
            most_correlated_gene_in_base_neighborhood_per_gene_per_base_block[:, base_block_index]
        correlated_gene_indices = Int[]
        friend_gene_indices = Int[]
        for gene_index in 1:n_genes
            name = most_correlated_gene_in_base_neighborhood_per_gene[gene_index]
            if !isempty(name)
                friend_gene_index = gene_name_to_index[name]
                push!(correlated_gene_indices, gene_index)
                push!(friend_gene_indices, friend_gene_index)
            end
        end
        correlated_gene_indices_per_base_block[base_block_index] = correlated_gene_indices
        friend_gene_indices_per_base_block[base_block_index] = friend_gene_indices
    end

    local_clusters_per_block = compute_local_clusters(;
        base_daf,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        min_cells_in_metacell,
        max_cells_in_metacell,
        min_metacell_total_UMIs,
        min_cells_dispersion_in_metacell,
        max_cells_dispersion_in_metacell,
        normalized_UMIs_quantile,
        min_module_UMIs,
        kmeans_rounds,
        sharpening_round,
        improvement_half_life,
        is_found_per_module_per_block,
        module_index_per_gene_per_block,
        block_index_per_cell,
        mean_metacell_cells_per_block,
        mean_linear_fraction_in_environment_cells_per_module_per_block,
        std_linear_fraction_in_environment_cells_per_module_per_block,
        std_fraction_regularization_per_block,
        base_block_index_per_cell,
        is_in_base_neighborhood_per_other_base_block_per_base_block,
        correlated_gene_indices_per_base_block,
        friend_gene_indices_per_base_block,
        gene_fraction_regularization,
        max_n_kmeans_clusters,
        kmeans_sizes_max_buffers_per_thread,
        kmeans_buffer_pool,
        rng,
    )

    cells_of_sharp_metacells, block_index_per_sharp_metacell =
        combine_local_clusters(; local_clusters_per_block, n_base_metacells)

    outlier_certificates = eject_outliers!(;
        cells_of_sharp_metacells,
        block_index_per_sharp_metacell,
        min_cells_in_metacell,
        min_outlier_fold,
        outlier_UMIs_regularization,
        block_index_per_cell,
        is_found_per_module_per_block,
        module_index_per_gene_per_block,
        mean_linear_fraction_in_environment_cells_per_module_per_block,
        std_linear_fraction_in_environment_cells_per_module_per_block,
        std_fraction_regularization_per_block,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        name_per_block,
        name_per_module = axis_vector(base_daf, "module"),
        n_cells,
        n_genes,
        n_modules = axis_length(base_daf, "module"),
        n_blocks,
    )

    name_per_sharp_metacell = group_names(axis_vector(base_daf, "cell"), cells_of_sharp_metacells; prefix)  # NOJET
    sharp_metacell_name_per_cell = fill("", n_cells)
    for (sharp_metacell_name, cells_of_sharp_metacell) in zip(name_per_sharp_metacell, cells_of_sharp_metacells)
        sharp_metacell_name_per_cell[cells_of_sharp_metacell] .= sharp_metacell_name
    end

    # Each outlier cell names the block it was clustered in, and the metacell it was ejected from (the metacell whose
    # punctuated expectation flagged it, which always survives). Both are the empty string for a cell that kept its
    # metacell, and for base outliers and excluded cells (which are never clustered).
    outlier_in_block_name_per_cell = fill("", n_cells)
    outlier_in_metacell_name_per_cell = fill("", n_cells)
    for cell_index in 1:n_cells
        if outlier_certificates.is_outlier_per_cell[cell_index]
            outlier_in_block_name_per_cell[cell_index] = name_per_block[block_index_per_cell[cell_index]]
            outlier_in_metacell_name_per_cell[cell_index] =
                name_per_sharp_metacell[outlier_certificates.in_metacell_index_per_cell[cell_index]]
        end
    end

    add_axis!(sharp_daf, "metacell", name_per_sharp_metacell; overwrite)
    set_vector!(sharp_daf, "cell", "metacell", sharp_metacell_name_per_cell; overwrite)
    set_vector!(sharp_daf, "metacell", "base_block", name_per_block[block_index_per_sharp_metacell]; overwrite)
    set_vector!(sharp_daf, "cell", "prev_block.outlier_in", outlier_in_block_name_per_cell; overwrite)
    set_vector!(sharp_daf, "cell", "metacell.outlier_in", outlier_in_metacell_name_per_cell; overwrite)
    set_vector!(sharp_daf, "cell", "prev_block.outlier_by", outlier_certificates.by_prev_block_name_per_cell; overwrite)
    set_vector!(
        sharp_daf,
        "cell",
        "prev_module.outlier_by",
        outlier_certificates.by_prev_module_name_per_cell;
        overwrite,
    )
    set_vector!(sharp_daf, "cell", "outlier_expected_UMIs", outlier_certificates.expected_UMIs_per_cell; overwrite)
    set_vector!(sharp_daf, "cell", "outlier_actual_UMIs", outlier_certificates.actual_UMIs_per_cell; overwrite)

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
    mean_linear_fraction_in_environment_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_linear_fraction_in_environment_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    kmeans_buffer_pool::Channel{KMeansBuffers{Float32}},
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

    # Channel pool of per-neighborhood-cell Float32 accumulators for `compute_z_score_per_found_module_per_region_cell!`'s
    # nested parallel-over-found-modules loop. Each inner task takes one for the duration of its module's gene-outer
    # fill, then puts it back. Sized to `nthreads()` because each take/put is a single non-yielding scope.
    z_score_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(z_score_accumulator_pool, Vector{Float32}(undef, max_n_neighborhood_cells))
    end
    region_position_per_cell_per_thread = [zeros(Int, n_cells) for _ in 1:maxthreadid()]

    max_n_neighborhood_clusters = maximum(
        max(Int(round(n_neighborhood_cells_per_block[block_index] / mean_metacell_cells_per_block[block_index])), 1) for block_index in 1:n_blocks
    )
    preferred_block_index_per_cluster_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_clusters) for _ in 1:maxthreadid()]
    n_cluster_cells_per_block_scratch_per_thread = [zeros(Int, n_blocks) for _ in 1:maxthreadid()]
    preferred_block_index_per_max_neighborhood_cell_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_cells) for _ in 1:maxthreadid()]

    parallel_loop_with_rng(  # NOJET
        1:n_blocks;
        rng,
        weights = n_neighborhood_cells_per_block,
        name = "compute_preferred_block_index_per_cell_per_block",
        progress = DebugProgress(
            sum(n_neighborhood_cells_per_block);
            group = :mcs_loops,
            desc = "preferred_block_index_per_cell_per_block",
        ),
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
        n_block_modules = n_modules_per_block[block_index]

        z_score_per_max_module_per_max_neighborhood_cell =
            z_score_per_max_module_per_max_neighborhood_cell_per_thread[threadid()]
        @views z_score_per_found_module_per_neighborhood_cell =
            z_score_per_max_module_per_max_neighborhood_cell[1:n_block_modules, 1:n_neighborhood_cells]

        setup_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_neighborhood_cell,
            is_found_per_module,
            block_index,
            indices_of_region_cells = indices_of_neighborhood_cells,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            z_score_accumulator_pool,
            region_position_per_cell_per_thread,
            n_genes,
            n_modules,
        )
        @assert n_block_modules == sum(is_found_per_module)

        n_neighborhood_clusters = max(Int(round(n_neighborhood_cells / mean_metacell_cells_per_block[block_index])), 1)
        kmeans_result = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                z_score_per_found_module_per_neighborhood_cell,
                n_neighborhood_clusters;
                buffer_pool = kmeans_buffer_pool,
                rounds = kmeans_rounds,
                rng,
            )
        end
        cluster_index_per_neighborhood_cell = assignments(kmeans_result)
        n_cells_per_cluster = counts(kmeans_result)

        preferred_block_index_per_cluster = preferred_block_index_per_cluster_per_thread[threadid()]
        n_cluster_cells_per_block = n_cluster_cells_per_block_scratch_per_thread[threadid()]

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
            n_cluster_cells_per_block,
            preferred_block_index_per_max_neighborhood_cell,
        )

        preferred_block_index_per_cell_per_block[block_index] =
            SparseVector(n_cells, indices_of_neighborhood_cells, Vector(preferred_block_index_per_neighborhood_cell))
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
    n_cluster_cells_per_block::AbstractVector{<:Integer},
    preferred_block_index_per_max_neighborhood_cell::AbstractVector{<:Integer},
)::AbstractVector{<:Integer}
    n_block_cells = n_cells_per_block[block_index]
    n_clusters = length(n_cells_per_cluster)
    n_neighborhood_cells = length(indices_of_neighborhood_cells)

    @views preferred_block_index_per_cluster[1:n_clusters] .= 0

    for cluster_index in 1:n_clusters
        fill!(n_cluster_cells_per_block, 0)
        n_block_cells_in_cluster = 0

        for neighborhood_cell_position in 1:n_neighborhood_cells
            if cluster_index_per_neighborhood_cell[neighborhood_cell_position] == cluster_index
                cell_block_index = block_index_per_cell[indices_of_neighborhood_cells[neighborhood_cell_position]]
                n_cluster_cells_per_block[cell_block_index] += 1
                if cell_block_index == block_index
                    n_block_cells_in_cluster += 1
                end
            end
        end

        if n_block_cells_in_cluster > 0
            most_frequent_block_index = block_index
            most_frequent_count = n_block_cells_in_cluster
            for other_block_index in eachindex(n_cluster_cells_per_block)
                if n_cluster_cells_per_block[other_block_index] > most_frequent_count
                    most_frequent_count = n_cluster_cells_per_block[other_block_index]
                    most_frequent_block_index = other_block_index
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
    max_cells_in_metacell::Integer,
    min_metacell_total_UMIs::Integer,
    min_cells_dispersion_in_metacell::AbstractFloat,
    max_cells_dispersion_in_metacell::AbstractFloat,
    normalized_UMIs_quantile::AbstractFloat,
    min_module_UMIs::Integer,
    kmeans_rounds::Integer,
    sharpening_round::Integer,
    improvement_half_life::Maybe{Integer},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    mean_linear_fraction_in_environment_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_linear_fraction_in_environment_cells_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    base_block_index_per_cell::AbstractVector{<:Integer},
    is_in_base_neighborhood_per_other_base_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    correlated_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    friend_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    gene_fraction_regularization::AbstractFloat,
    max_n_kmeans_clusters::Integer,
    kmeans_sizes_max_buffers_per_thread::AbstractVector{<:KmeansSizesBuffers},
    kmeans_buffer_pool::Channel{KMeansBuffers{Float32}},
    rng::AbstractRNG,
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
    expected_n_clusters_per_block = zeros(Int, n_blocks)

    n_cells_per_block = [count(==(block_index), block_index_per_cell) for block_index in 1:n_blocks]
    max_n_block_cells = maximum(n_cells_per_block)

    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)

    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]
    is_in_block_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]
    z_score_per_max_module_per_max_block_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_block_cells) for _ in 1:maxthreadid()]

    # Channel pool of per-cell Float32 accumulators for `compute_z_score_per_found_module_per_region_cell!`'s nested
    # parallel-over-found-modules loop. Each inner task takes one for the duration of its module's gene-outer fill,
    # then puts it back. Sized to `nthreads()` because each take/put is a single non-yielding scope.
    z_score_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(z_score_accumulator_pool, Vector{Float32}(undef, max_n_block_cells))
    end
    region_position_per_cell_per_thread = [zeros(Int, n_cells) for _ in 1:maxthreadid()]

    # `DispersionScratches` channel pool for the nested parallel dispersion loops (compute_all_dispersions! /
    # compute_changed_dispersions!) and the serial split-then-disperse path in `split_one_cluster!`. Sized to
    # `nthreads()` because every take/put is a single non-yielding scope, so each running task holds at most one
    # scratch.
    dispersion_scratches_pool = Channel{DispersionScratches}(nthreads())
    for _ in 1:nthreads()
        put!(dispersion_scratches_pool, DispersionScratches(; n_points = max_n_block_cells, n_modules))
    end

    # ===== Phase 1: per-block walk to a perfect K (or exhaustion) =====
    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        weights = n_cells_per_block,
        name = "compute_local_clusters.clustering_options",
        progress = DebugProgress(sum(n_cells_per_block); group = :mcs_loops, desc = "clustering_options"),
    ) do block_index, rng
        is_in_block_per_cell = is_in_block_per_cell_per_thread[threadid()]
        is_in_block_per_cell .= block_index_per_cell .== block_index
        block_cell_indices = findall(is_in_block_per_cell)
        n_block_cells = n_cells_per_block[block_index]
        @assert n_block_cells == length(block_cell_indices)
        if n_block_cells == 0
            return nothing
        end

        is_found_per_module = is_found_per_module_per_thread[threadid()]
        n_block_modules = n_modules_per_block[block_index]

        z_score_per_max_module_per_max_block_cell = z_score_per_max_module_per_max_block_cell_per_thread[threadid()]
        @views z_score_per_found_module_per_block_cell =
            z_score_per_max_module_per_max_block_cell[1:n_block_modules, 1:n_block_cells]

        # The inverted `module_index_per_gene` returned here is reused below by the Phase 1 K-walk dispersion context.
        gene_indices_per_module = setup_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_block_cell,
            is_found_per_module,
            block_index,
            indices_of_region_cells = block_cell_indices,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            z_score_accumulator_pool,
            region_position_per_cell_per_thread,
            n_genes,
            n_modules,
        )
        @assert n_block_modules == sum(is_found_per_module)

        n_block_clusters = max(Int(round(n_block_cells / mean_metacell_cells_per_block[block_index])), 1)
        expected_n_clusters_per_block[block_index] = n_block_clusters
        sizes_buffers = kmeans_sizes_max_buffers_per_thread[threadid()]
        @views weight_per_block_cell = total_UMIs_per_cell[block_cell_indices]

        dispersion_context = DispersionContext(;
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            block_cell_indices,
            is_found_per_module,
            gene_indices_per_module,
            normalized_UMIs_quantile,
            min_module_UMIs,
            min_cells_dispersion = min_cells_dispersion_in_metacell,
            max_cells_dispersion = max_cells_dispersion_in_metacell,
            min_cluster_size = min_cells_in_metacell,
        )

        # The block's single cluster is also its last one, so it is kept even if it is too small.
        if n_block_clusters == 1
            local_clusters_per_block[block_index] = LocalClusters(;
                n_clusters = 1,
                block_cell_indices,
                cluster_index_per_block_cell = fill(1, n_block_cells),
            )
            return nothing
        end

        candidates = flame_timed("kmeans_with_sizes_candidates") do
            return kmeans_with_sizes_candidates(
                z_score_per_found_module_per_block_cell,
                weight_per_block_cell,
                n_block_clusters;
                max_k = 2 * n_block_clusters,
                min_cluster_size = min_cells_in_metacell,
                max_cluster_size = max_cells_in_metacell,
                min_cluster_weight = min_metacell_total_UMIs,
                kmeans_rounds,
                sizes_buffers,
                kmeans_buffer_pool,
                dispersion_scratches_pool,
                dispersion_context,
                rng,
            )
        end
        solution_candidates_per_block[block_index] = candidates
        # Tentative LocalClusters from the baseline (candidates[1]); Phase 2 may overwrite for walkable blocks.
        local_clusters_per_block[block_index] =
            build_local_clusters_from_candidate(candidates.candidates[1], block_cell_indices)
        return nothing
    end

    # `block_cell_indices_per_block` is constant across optimization passes (block membership is fixed by Phase 1;
    # only cluster assignments within a block change). Same for `n_points_per_base_block` (depends only on
    # `base_block_index_per_cell` and the neighborhood matrix).
    block_cell_indices_per_block = [
        local_clusters === nothing ? Int[] : Vector{Int}(local_clusters.block_cell_indices) for
        local_clusters in local_clusters_per_block
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

    (
        affected_base_block_indices_per_walkable_block,
        gene_index_per_friend_per_walkable_block,
        friend_position_per_series_per_base_block_per_walkable_block,
        block_cell_position_per_point_per_base_block_per_walkable_block,
    ) = precompute_walkable_indirection(
        walkable_block_indices,
        block_cell_indices_per_block,
        base_block_index_per_cell,
        is_in_base_neighborhood_per_other_base_block_per_base_block,
        friend_gene_indices_per_base_block,
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
    max_n_friends = 0
    max_n_walkable_block_cells = 0
    for walkable_position in 1:n_walkable_blocks
        block_index = walkable_block_indices[walkable_position]
        solution_candidates = solution_candidates_per_block[block_index]
        n_candidates = length(solution_candidates.candidates)
        n_block_cells = solution_candidates.n_points
        max_n_walkable_block_cells = max(max_n_walkable_block_cells, n_block_cells)
        max_n_friends = max(max_n_friends, length(gene_index_per_friend_per_walkable_block[walkable_position]))
        for candidate_index in 1:n_candidates
            push!(walkable_position_per_work_item, walkable_position)
            push!(candidate_index_per_work_item, candidate_index)
            push!(n_block_cells_per_work_item, n_block_cells)
        end
        work_item_start_per_walkable_block[walkable_position + 1] =
            work_item_start_per_walkable_block[walkable_position] + n_candidates
    end
    n_work_items = length(walkable_position_per_work_item)

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
            gene_index_per_friend = Int[],
            friend_position_per_series_per_base_block = Dict{Int, Vector{Int}}(),
            block_cell_position_per_point_per_base_block = Dict{Int, Vector{Int}}(),
            # Per-thread scratch.
            total_UMIs_per_max_cluster = Vector{Float64}(undef, max_n_kmeans_clusters),
            UMIs_per_friend_per_max_cluster = Matrix{Float64}(undef, max_n_friends, max_n_kmeans_clusters),
            variable_per_block_cell_per_friend = Matrix{Float32}(undef, max_n_walkable_block_cells, max_n_friends),
            is_active_per_block_cell = Vector{Bool}(undef, max_n_walkable_block_cells),
            punctuated_total_per_block_cell = Vector{Float64}(undef, max_n_walkable_block_cells),
        ) for _ in 1:maxthreadid()
    ]
    delta_per_work_item = Vector{Float64}(undef, n_work_items)

    # Per walkable block, the candidate index currently treated as the baseline. Initialized to `1` (the minimal-K
    # candidate, which is also what `local_clusters_per_block` was built from after Phase 1). Each optimization pass
    # updates this to whatever candidate it picks.
    chosen_candidate_index_per_walkable_block = fill(1, n_walkable_blocks)

    # Per-walkable-block scratches sized to that walkable block's own dimensions. Pass A (parallel argmax) populates
    # the variable + is_active scratches for every changed walkable block; Pass B (parallel `replace_group!`) reads them
    # per (changed_walkable, affected base block) without re-populating - the two-pass split eliminates the populate
    # redundancy that a parallel-by-base-block commit would otherwise pay. The `UMI_per_friend_per_block_cell`
    # cache (filled once below) lets every populate call read this walkable block's UMIs straight from a sequential
    # Float32 column instead of indexing the global `UMIs_per_cell_per_gene` matrix on every (cell, friend gene).
    variable_per_block_cell_per_friend_per_walkable_block = Vector{Matrix{Float32}}(undef, n_walkable_blocks)
    is_active_per_block_cell_per_walkable_block = Vector{Vector{Bool}}(undef, n_walkable_blocks)
    UMI_per_friend_per_block_cell_per_walkable_block = Vector{Matrix{Float32}}(undef, n_walkable_blocks)
    n_block_cells_per_walkable_block = Vector{Int}(undef, n_walkable_blocks)
    UMIs_fill_weight_per_walkable_block = Vector{Int}(undef, n_walkable_blocks)
    for walkable_position in 1:n_walkable_blocks
        block_index = walkable_block_indices[walkable_position]
        n_block_cells = solution_candidates_per_block[block_index].n_points
        n_friends = length(gene_index_per_friend_per_walkable_block[walkable_position])
        variable_per_block_cell_per_friend_per_walkable_block[walkable_position] =
            Matrix{Float32}(undef, n_block_cells, n_friends)
        is_active_per_block_cell_per_walkable_block[walkable_position] = Vector{Bool}(undef, n_block_cells)
        UMI_per_friend_per_block_cell_per_walkable_block[walkable_position] =
            Matrix{Float32}(undef, n_friends, n_block_cells)
        n_block_cells_per_walkable_block[walkable_position] = n_block_cells
        UMIs_fill_weight_per_walkable_block[walkable_position] = n_block_cells * n_friends
    end

    # One-time fill of the per-walkable UMI cache. Outer friend_position so each `UMIs_per_cell_per_gene[:, gene_index]`
    # column is touched contiguously; inner block_cell scatter reads only within that already-loaded column.
    parallel_loop_wo_rng(
        1:n_walkable_blocks;
        weights = UMIs_fill_weight_per_walkable_block,
        name = "compute_local_clusters.UMIs_cache_fill",
        progress = DebugProgress(
            sum(UMIs_fill_weight_per_walkable_block);
            group = :mcs_loops,
            desc = "UMIs_cache_fill",
        ),
    ) do walkable_position
        block_index = walkable_block_indices[walkable_position]
        block_cell_indices = block_cell_indices_per_block[block_index]
        gene_index_per_friend = gene_index_per_friend_per_walkable_block[walkable_position]
        UMI_per_friend_per_block_cell = UMI_per_friend_per_block_cell_per_walkable_block[walkable_position]
        n_block_cells = length(block_cell_indices)
        n_friends = length(gene_index_per_friend)
        @inbounds for friend_position in 1:n_friends
            gene_index = gene_index_per_friend[friend_position]
            for block_cell_position in 1:n_block_cells
                UMI_per_friend_per_block_cell[friend_position, block_cell_position] =
                    Float32(UMIs_per_cell_per_gene[block_cell_indices[block_cell_position], gene_index])
            end
        end
        return nothing
    end

    # Per-base-block locks serialize `replace_group!` updates to each base block's grouped correlations during the
    # commit phase. Allocated once and reused across all optimization passes.
    lock_per_base_block = [SpinLock() for _ in 1:n_base_blocks]

    # One-time setup: build the per-base-block grouped correlations from the Phase-1 `local_clusters_per_block`. All
    # later updates happen in place via the optimization loop's `replace_group!`.
    (total_UMIs_per_baseline_metacell, UMIs_per_gene_per_baseline_metacell, baseline_metacell_index_per_cell, _) =
        compute_baseline_metacell_aggregates(
            local_clusters_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            n_genes,
            n_cells,
        )

    parallel_loop_wo_rng(
        1:n_base_blocks;
        weights = n_points_per_base_block,
        name = "compute_local_clusters.build_grouped",
        progress = DebugProgress(sum(n_points_per_base_block); group = :mcs_loops, desc = "build_grouped"),
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
            UMIs_per_gene_per_baseline_metacell,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            correlated_gene_indices = correlated_gene_indices_per_base_block[base_block_index],
            friend_gene_indices = friend_gene_indices_per_base_block[base_block_index],
            gene_fraction_regularization,
        )
        grouped_per_base_block[base_block_index] = grouped
        @views group_index_per_block_per_base_block[:, base_block_index] .= group_index_per_block
        baseline_mean_correlation_per_base_block[base_block_index] = Float32(mean_correlation(grouped))
        return nothing
    end

    baseline_correlation = mean(baseline_mean_correlation_per_base_block)  # NOLINT
    @debug "Baseline mean correlation: $(baseline_correlation)" _group = :mcs_results

    epsilon = 1e-6
    if improvement_half_life === nothing
        cooldown_margin_per_k_distance = 0.0
    else
        max_cooldown_margin_per_k_distance = 100 * epsilon
        cooldown_margin_per_k_distance =
            max_cooldown_margin_per_k_distance * (1 - 2.0^(-(sharpening_round - 1) / improvement_half_life))
    end
    @debug "Cooldown margin per K distance: $(cooldown_margin_per_k_distance)" _group = :mcs_results
    pass_index = 0
    delta_correlation = 1
    while (delta_correlation < -epsilon || epsilon < delta_correlation) && pass_index < 50
        pass_index += 1

        parallel_loop_wo_rng(
            1:n_work_items;
            weights = n_block_cells_per_work_item,
            progress = DebugProgress(
                sum(n_block_cells_per_work_item);
                group = :mcs_loops,
                desc = "pass_$(pass_index).evaluate_options",
            ),
            name = "compute_local_clusters.pass_$(pass_index).evaluate_options",
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
                delta_context.gene_index_per_friend = gene_index_per_friend_per_walkable_block[walkable_position]
                delta_context.friend_position_per_series_per_base_block =
                    friend_position_per_series_per_base_block_per_walkable_block[walkable_position]
                delta_context.block_cell_position_per_point_per_base_block =
                    block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position]
                delta_per_work_item[work_index] = compute_delta_correlation(
                    candidate,
                    solution_candidates.n_points,
                    delta_context,
                    UMI_per_friend_per_block_cell_per_walkable_block[walkable_position],
                )
            end
            return nothing
        end

        n_changes = Threads.Atomic{Int}(0)
        parallel_loop_wo_rng(
            1:n_walkable_blocks;
            weights = n_block_cells_per_walkable_block,
            name = "compute_local_clusters.pass_$(pass_index).select_options",
            progress = DebugProgress(
                sum(n_block_cells_per_walkable_block);
                group = :mcs_loops,
                desc = "pass_$(pass_index).select_options",
            ),
        ) do walkable_position
            block_index = walkable_block_indices[walkable_position]
            solution_candidates = solution_candidates_per_block[block_index]
            candidates = solution_candidates.candidates
            previous_candidate_index = chosen_candidate_index_per_walkable_block[walkable_position]
            expected_n_clusters = expected_n_clusters_per_block[block_index]
            current_distance_from_expected = abs(candidates[previous_candidate_index].k - expected_n_clusters)
            best_candidate_index = previous_candidate_index
            # Seed with the current candidate's penalized score: its delta is 0, but it still owes its K-distance penalty.
            best_score = -cooldown_margin_per_k_distance * current_distance_from_expected
            slice_start = work_item_start_per_walkable_block[walkable_position]
            slice_end = work_item_start_per_walkable_block[walkable_position + 1] - 1
            for work_index in slice_start:slice_end
                candidate_index = candidate_index_per_work_item[work_index]
                candidate_distance_from_expected = abs(candidates[candidate_index].k - expected_n_clusters)
                score =
                    delta_per_work_item[work_index] - cooldown_margin_per_k_distance * candidate_distance_from_expected
                if score > best_score + epsilon
                    best_score = score
                    best_candidate_index = candidate_index
                end
            end
            if best_candidate_index == previous_candidate_index
                return nothing
            end
            chosen = candidates[best_candidate_index]
            Threads.atomic_add!(n_changes, 1)
            chosen_candidate_index_per_walkable_block[walkable_position] = best_candidate_index
            local_clusters_per_block[block_index] =
                build_local_clusters_from_candidate(chosen, block_cell_indices_per_block[block_index])

            # Populate this walkable's persistent variable + mask scratch under the chosen candidate.
            delta_context = delta_context_per_thread[threadid()]
            delta_context.walkable_block_index = block_index
            delta_context.block_cell_indices = block_cell_indices_per_block[block_index]
            delta_context.affected_base_block_indices =
                affected_base_block_indices_per_walkable_block[walkable_position]
            delta_context.gene_index_per_friend = gene_index_per_friend_per_walkable_block[walkable_position]
            delta_context.friend_position_per_series_per_base_block =
                friend_position_per_series_per_base_block_per_walkable_block[walkable_position]
            delta_context.block_cell_position_per_point_per_base_block =
                block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position]
            populate_candidate_scratches!(
                chosen,
                solution_candidates.n_points,
                delta_context,
                UMI_per_friend_per_block_cell_per_walkable_block[walkable_position],
                variable_per_block_cell_per_friend_per_walkable_block[walkable_position],
                is_active_per_block_cell_per_walkable_block[walkable_position],
            )

            # Apply the change directly to each affected base block's grouped correlations under its lock.
            for base_block_index in affected_base_block_indices_per_walkable_block[walkable_position]
                group_index_in_base_block = group_index_per_block_per_base_block[block_index, base_block_index]
                @assert group_index_in_base_block > 0
                grouped_for_base_block = grouped_per_base_block[base_block_index]
                block_cell_position_per_point =
                    block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position][base_block_index]
                friend_position_per_series =
                    friend_position_per_series_per_base_block_per_walkable_block[walkable_position][base_block_index]
                lock(lock_per_base_block[base_block_index]) do
                    replace_group!(
                        grouped_for_base_block,
                        group_index_in_base_block,
                        variable_per_block_cell_per_friend_per_walkable_block[walkable_position],
                        is_active_per_block_cell_per_walkable_block[walkable_position],
                        block_cell_position_per_point,
                        friend_position_per_series,
                    )
                    baseline_mean_correlation_per_base_block[base_block_index] =
                        Float32(mean_correlation(grouped_for_base_block))
                    return nothing
                end
            end
            return nothing
        end

        pass_correlation = mean(baseline_mean_correlation_per_base_block)  # NOLINT
        delta_correlation = pass_correlation - baseline_correlation
        baseline_correlation = pass_correlation

        @debug (
            "Pass $(pass_index) " *
            "changes: $(n_changes[]) " *
            "out of: $(n_blocks) " *
            "delta: $(delta_correlation) " *
            "updated: $(baseline_correlation)"
        ) _group = :mcs_results
    end

    return local_clusters_per_block
end

# Collect the clusters of all the blocks into a single flat list of the candidate sharp metacells (every clustered cell
# is in exactly one of them), and the block each was clustered in.
function combine_local_clusters(;
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    n_base_metacells::Integer,
)::Tuple{Vector{Vector{Int}}, Vector{Int}}
    cells_of_new_metacells = Vector{Vector{Int}}()
    block_index_per_new_metacell = Vector{Int}()

    for (block_index, local_clusters) in enumerate(local_clusters_per_block)
        if local_clusters !== nothing
            for cluster_index in 1:local_clusters.n_clusters
                cell_indices_of_new_metacell =
                    local_clusters.block_cell_indices[local_clusters.cluster_index_per_block_cell .== cluster_index]
                push!(cells_of_new_metacells, Vector{Int}(cell_indices_of_new_metacell))
                push!(block_index_per_new_metacell, block_index)
            end
        end
    end

    @debug ("Metacells Original: $(n_base_metacells) Sharpened: $(length(block_index_per_new_metacell))") _group =
        :mcs_results

    return (cells_of_new_metacells, block_index_per_new_metacell)
end

# The flattened enumeration of every found module instance across all blocks in a repository, indexed by
# `total_found_module_index`. The parallel vectors hold each instance's block and module names and the gene indices of
# the module in that block.
struct TotalFoundModules
    block_name_per_total_found_module::Vector{String}
    module_name_per_total_found_module::Vector{String}
    gene_indices_per_total_found_module::Vector{Vector{Int}}
end

# Collect every found module instance across all blocks into a single flat enumeration.
function collect_total_found_modules(;
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    name_per_block::AbstractVector{<:AbstractString},
    name_per_module::AbstractVector{<:AbstractString},
)::TotalFoundModules
    n_modules = size(is_found_per_module_per_block, 1)
    n_blocks = size(is_found_per_module_per_block, 2)
    n_genes = size(module_index_per_gene_per_block, 1)

    block_name_per_total_found_module = String[]
    module_name_per_total_found_module = String[]
    gene_indices_per_total_found_module = Vector{Vector{Int}}()
    for block_index in 1:n_blocks
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        @views is_found_per_module = is_found_per_module_per_block[:, block_index]
        gene_indices_per_module = [Int[] for _ in 1:n_modules]
        for gene_index in 1:n_genes
            module_index = module_index_per_gene[gene_index]
            if module_index > 0
                push!(gene_indices_per_module[module_index], gene_index)
            end
        end
        for module_index in 1:n_modules
            if is_found_per_module[module_index]
                push!(block_name_per_total_found_module, name_per_block[block_index])
                push!(module_name_per_total_found_module, name_per_module[module_index])
                push!(gene_indices_per_total_found_module, gene_indices_per_module[module_index])
            end
        end
    end

    return TotalFoundModules(
        block_name_per_total_found_module,
        module_name_per_total_found_module,
        gene_indices_per_total_found_module,
    )
end

# The per-cell certificates of the outlier detection. For each outlier cell, the metacell it was ejected from, the name of
# the most-deviant previous round's block and module that flagged it, and the expected and actual UMIs of that module's
# genes in the cell. All entries are the "not an outlier" defaults (zero / empty string) for cells that are not outliers.
struct OutlierCertificates
    is_outlier_per_cell::BitVector
    in_metacell_index_per_cell::Vector{Int}
    by_prev_block_name_per_cell::Vector{String}
    by_prev_module_name_per_cell::Vector{String}
    expected_UMIs_per_cell::Vector{Float32}
    actual_UMIs_per_cell::Vector{Float32}
end

# A single outlier fold factor found for a cell against one found module instance. Collected per thread during the
# parallel sweep over module instances, then reduced to the per-cell maximal fold in `compute_outlier_certificates`.
struct OutlierCandidate
    cell_index::Int
    total_found_module_index::Int
    fold::Float32
    expected_UMIs::Float32
    actual_UMIs::Float32
end

# Detect the outlier cells: cells whose expression of some base found module wildly exceeds what their metacell
# predicts. For every found module instance in the previous round's repository (across all blocks, not just the modules of a cell's
# own region), and for each grouped cell that expresses the module, compare the cell's actual UMIs of the module's genes
# to the expected UMIs. The expectation is punctuated - the module's linear fraction in the cell's metacell excluding the
# cell itself, times the cell's total UMIs - so a strong outlier cannot inflate its own baseline. A cell whose maximal
# fold (across all module instances) exceeds `min_outlier_fold` is an outlier, with the certificate of the most-deviant
# module instance. A metacell is never emptied - every outlier names the metacell it was ejected from, so if all the
# cells of a metacell are outliers, none of them is.
#
# The loop is parallel over module instances. `UMIs_per_cell_per_gene` is a CSC sparse matrix whose columns are genes,
# so each instance walks only the nonzeros of its module's gene columns (the cells that actually express it); cells with
# no UMIs of the module cannot be deviant for it (their fold is non-positive), so they are never visited. Per thread we
# reuse a dense `actual_per_cell` accumulator and a `is_touched_per_cell` mask, tracking the touched cells (and
# metacells) in preallocated index vectors so both the per-instance reduction and the reset are O(touched), not
# O(n_cells). Each instance emits only its fold-exceeding candidates (rare) to a per-thread buffer, reduced serially to
# the per-cell maximum at the end.
function compute_outlier_certificates(;
    cells_of_sharp_metacells::AbstractVector{<:AbstractVector{<:Integer}},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    name_per_block::AbstractVector{<:AbstractString},
    name_per_module::AbstractVector{<:AbstractString},
    min_outlier_fold::AbstractFloat,
    outlier_UMIs_regularization::AbstractFloat,
    n_cells::Integer,
)::OutlierCertificates
    n_sharp_metacells = length(cells_of_sharp_metacells)

    # The sparse-column sweep requires a CSC matrix whose columns are genes; unwrap the read-only Daf matrix to its
    # underlying `SparseMatrixCSC` (no copy).
    @assert issparse(UMIs_per_cell_per_gene)
    sparse_UMIs_per_cell_per_gene = mutable_array(UMIs_per_cell_per_gene)::SparseMatrixCSC

    # Allocated once with the "not an outlier" defaults; the candidate reduction below fills in the outlier cells.
    is_outlier_per_cell = falses(n_cells)
    in_metacell_index_per_cell = zeros(Int, n_cells)
    by_prev_block_name_per_cell = fill("", n_cells)
    by_prev_module_name_per_cell = fill("", n_cells)
    expected_UMIs_per_cell = zeros(Float32, n_cells)
    actual_UMIs_per_cell = zeros(Float32, n_cells)

    # Per-cell metacell membership (0 = ungrouped) and per-metacell total UMIs, both read-only inside the parallel loop.
    # The single cell of a single-cell metacell has no punctuated baseline to compare against, so it is left ungrouped -
    # it can never be an outlier.
    metacell_index_per_cell = zeros(Int32, n_cells)
    metacell_total_UMIs_per_metacell = zeros(Float64, n_sharp_metacells)
    for (sharp_metacell_index, cell_indices) in enumerate(cells_of_sharp_metacells)
        if length(cell_indices) > 1
            for cell_index in cell_indices
                metacell_index_per_cell[cell_index] = sharp_metacell_index
                metacell_total_UMIs_per_metacell[sharp_metacell_index] += total_UMIs_per_cell[cell_index]
            end
        end
    end

    base_total_modules = collect_total_found_modules(;
        is_found_per_module_per_block,
        module_index_per_gene_per_block,
        name_per_block,
        name_per_module,
    )
    n_total_found_modules = length(base_total_modules.gene_indices_per_total_found_module)

    # Weight each instance by the total nonzeros of its module's gene columns, so the heaviest instances are dispatched
    # first for load balance.
    weight_per_total_found_module = [
        sum(gene_index -> length(nzrange(sparse_UMIs_per_cell_per_gene, gene_index)), gene_indices_in_module; init = 0) for gene_indices_in_module in base_total_modules.gene_indices_per_total_found_module
    ]

    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)

    actual_per_cell_per_thread = [zeros(Float32, n_cells) for _ in 1:maxthreadid()]
    is_touched_per_cell_per_thread = [falses(n_cells) for _ in 1:maxthreadid()]
    touched_cells_per_thread = [Vector{Int32}(undef, n_cells) for _ in 1:maxthreadid()]
    module_total_UMIs_per_metacell_per_thread = [zeros(Float64, n_sharp_metacells) for _ in 1:maxthreadid()]
    touched_metacells_per_thread = [Vector{Int32}(undef, n_sharp_metacells) for _ in 1:maxthreadid()]
    candidates_per_thread = [OutlierCandidate[] for _ in 1:maxthreadid()]

    regularization = Float64(outlier_UMIs_regularization)

    parallel_loop_wo_rng(
        1:n_total_found_modules;
        name = "compute_outlier_certificates",
        weights = weight_per_total_found_module,
        progress = DebugProgress(sum(weight_per_total_found_module); group = :mcs_loops, desc = "outlier_certificates"),
    ) do total_found_module_index
        actual_per_cell = actual_per_cell_per_thread[threadid()]
        is_touched_per_cell = is_touched_per_cell_per_thread[threadid()]
        touched_cells = touched_cells_per_thread[threadid()]
        module_total_UMIs_per_metacell = module_total_UMIs_per_metacell_per_thread[threadid()]
        touched_metacells = touched_metacells_per_thread[threadid()]
        candidates = candidates_per_thread[threadid()]

        gene_indices_in_module = base_total_modules.gene_indices_per_total_found_module[total_found_module_index]

        # Gather the module's UMIs into `actual_per_cell` for every grouped cell that expresses it, walking only the
        # nonzeros of the module's gene columns and recording each cell once in `touched_cells`.
        n_touched_cells = 0
        @inbounds for gene_index in gene_indices_in_module
            for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
                cell_index = row_index_per_stored[stored_index]
                if metacell_index_per_cell[cell_index] > 0
                    if !is_touched_per_cell[cell_index]
                        is_touched_per_cell[cell_index] = true
                        n_touched_cells += 1
                        touched_cells[n_touched_cells] = cell_index
                    end
                    actual_per_cell[cell_index] += UMIs_per_stored[stored_index]
                end
            end
        end

        # Sum the module's UMIs per metacell (each touched cell counted once), recording each metacell once.
        n_touched_metacells = 0
        @inbounds for touched_index in 1:n_touched_cells
            cell_index = touched_cells[touched_index]
            metacell_index = metacell_index_per_cell[cell_index]
            if module_total_UMIs_per_metacell[metacell_index] == 0.0
                n_touched_metacells += 1
                touched_metacells[n_touched_metacells] = metacell_index
            end
            module_total_UMIs_per_metacell[metacell_index] += actual_per_cell[cell_index]
        end

        # Punctuated fold per touched cell; emit the fold-exceeding ones as candidates.
        @inbounds for touched_index in 1:n_touched_cells
            cell_index = touched_cells[touched_index]
            metacell_index = metacell_index_per_cell[cell_index]
            actual_UMIs = actual_per_cell[cell_index]
            cell_total_UMIs = total_UMIs_per_cell[cell_index]
            punctuated_total_UMIs = metacell_total_UMIs_per_metacell[metacell_index] - cell_total_UMIs
            @assert punctuated_total_UMIs > 0
            punctuated_module_UMIs = module_total_UMIs_per_metacell[metacell_index] - actual_UMIs
            module_linear_fraction = punctuated_module_UMIs / punctuated_total_UMIs
            expected_UMIs = module_linear_fraction * cell_total_UMIs
            fold = log2((actual_UMIs + regularization) / (expected_UMIs + regularization))
            if fold > min_outlier_fold
                push!(
                    candidates,
                    OutlierCandidate(
                        cell_index,
                        total_found_module_index,
                        Float32(fold),
                        Float32(expected_UMIs),
                        Float32(actual_UMIs),
                    ),
                )
            end
        end

        # Reset only the touched scratch entries for the next instance on this thread.
        @inbounds for touched_index in 1:n_touched_cells
            cell_index = touched_cells[touched_index]
            actual_per_cell[cell_index] = 0.0f0
            is_touched_per_cell[cell_index] = false
        end
        @inbounds for touched_index in 1:n_touched_metacells
            module_total_UMIs_per_metacell[touched_metacells[touched_index]] = 0.0
        end
        return nothing
    end

    # Reduce the per-thread candidates to the per-cell maximal fold. Only fold-exceeding candidates exist, so every cell
    # that appears is an outlier, recording the certificate of the module instance that flagged it most strongly.
    max_fold_per_cell = fill(-Inf32, n_cells)
    for candidates in candidates_per_thread
        for candidate in candidates
            cell_index = candidate.cell_index
            if candidate.fold > max_fold_per_cell[cell_index]
                max_fold_per_cell[cell_index] = candidate.fold
                is_outlier_per_cell[cell_index] = true
                in_metacell_index_per_cell[cell_index] = metacell_index_per_cell[cell_index]
                by_prev_block_name_per_cell[cell_index] =
                    base_total_modules.block_name_per_total_found_module[candidate.total_found_module_index]
                by_prev_module_name_per_cell[cell_index] =
                    base_total_modules.module_name_per_total_found_module[candidate.total_found_module_index]
                expected_UMIs_per_cell[cell_index] = candidate.expected_UMIs
                actual_UMIs_per_cell[cell_index] = candidate.actual_UMIs
            end
        end
    end

    # Every outlier names the metacell it was ejected from, so a metacell must keep at least one cell. If all the cells of
    # a metacell are outliers, none of them is.
    for cell_indices in cells_of_sharp_metacells
        if all(cell_index -> is_outlier_per_cell[cell_index], cell_indices)
            for cell_index in cell_indices
                is_outlier_per_cell[cell_index] = false
                in_metacell_index_per_cell[cell_index] = 0
                by_prev_block_name_per_cell[cell_index] = ""
                by_prev_module_name_per_cell[cell_index] = ""
                expected_UMIs_per_cell[cell_index] = 0
                actual_UMIs_per_cell[cell_index] = 0
            end
        end
    end

    n_outlier_cells = sum(is_outlier_per_cell)
    n_grouped_cells = sum(length(cell_indices) for cell_indices in cells_of_sharp_metacells; init = 0)
    @debug (
        "Outlier cells: $(n_outlier_cells)" *
        " ($(percent(n_outlier_cells, n_grouped_cells)) of grouped, $(percent(n_outlier_cells, n_cells)) of all)"
    ) _group = :mcs_results

    return OutlierCertificates(
        is_outlier_per_cell,
        in_metacell_index_per_cell,
        by_prev_block_name_per_cell,
        by_prev_module_name_per_cell,
        expected_UMIs_per_cell,
        actual_UMIs_per_cell,
    )
end

# Eject the outlier cells from the sharp metacells. Ejecting the outliers can drop a metacell below
# `min_cells_in_metacell` cells; such a metacell is dissolved, all its cells (outliers included) are relocated to the
# nearest surviving metacell of its block, and the outliers are detected again - so a cell relocated into a metacell it
# does not fit is caught. This repeats until no metacell is dissolved. Mutates `cells_of_sharp_metacells` (and the
# parallel `block_index_per_sharp_metacell`) into the final metacells, and returns the certificates of the cells ejected
# from them.
function eject_outliers!(;
    cells_of_sharp_metacells::Vector{Vector{Int}},
    block_index_per_sharp_metacell::Vector{Int},
    min_cells_in_metacell::Integer,
    min_outlier_fold::AbstractFloat,
    outlier_UMIs_regularization::AbstractFloat,
    block_index_per_cell::AbstractVector{<:Integer},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    mean_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    name_per_block::AbstractVector{<:AbstractString},
    name_per_module::AbstractVector{<:AbstractString},
    n_cells::Integer,
    n_genes::Integer,
    n_modules::Integer,
    n_blocks::Integer,
)::OutlierCertificates
    while true
        outlier_certificates = compute_outlier_certificates(;
            cells_of_sharp_metacells,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            name_per_block,
            name_per_module,
            min_outlier_fold,
            outlier_UMIs_regularization,
            n_cells,
        )
        is_outlier_per_cell = outlier_certificates.is_outlier_per_cell

        is_dissolved_per_sharp_metacell = pick_dissolved_metacells(;
            cells_of_sharp_metacells,
            block_index_per_sharp_metacell,
            is_outlier_per_cell,
            min_cells_in_metacell,
            n_blocks,
        )

        if !any(is_dissolved_per_sharp_metacell)
            n_ejected_cells = 0
            for (sharp_metacell_index, cell_indices) in enumerate(cells_of_sharp_metacells)
                kept_cell_indices = filter(cell_index -> !is_outlier_per_cell[cell_index], cell_indices)
                n_ejected_cells += length(cell_indices) - length(kept_cell_indices)
                cells_of_sharp_metacells[sharp_metacell_index] = kept_cell_indices
            end
            @debug (
                "Ejected outlier cells: $(n_ejected_cells)" *
                " Remaining metacells: $(length(cells_of_sharp_metacells))"
            ) _group = :mcs_results
            return outlier_certificates
        end

        relocate_dissolved_metacell_cells!(;
            cells_of_sharp_metacells,
            block_index_per_sharp_metacell,
            is_dissolved_per_sharp_metacell,
            is_outlier_per_cell,
            block_index_per_cell,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            n_cells,
            n_genes,
            n_modules,
            n_blocks,
        )
    end
end

# The metacells to dissolve: those left with fewer than `min_cells_in_metacell` cells once their outliers are ejected. A
# block's last metacell is never dissolved - there would be nowhere to relocate its cells to - so if all the metacells of
# a block would be dissolved, the one that keeps the most cells survives to absorb the rest.
function pick_dissolved_metacells(;
    cells_of_sharp_metacells::AbstractVector{<:AbstractVector{<:Integer}},
    block_index_per_sharp_metacell::AbstractVector{<:Integer},
    is_outlier_per_cell::Union{AbstractVector{Bool}, BitVector},
    min_cells_in_metacell::Integer,
    n_blocks::Integer,
)::BitVector
    n_kept_cells_per_sharp_metacell = [
        count(cell_index -> !is_outlier_per_cell[cell_index], cell_indices) for cell_indices in cells_of_sharp_metacells
    ]
    is_dissolved_per_sharp_metacell = BitVector(n_kept_cells_per_sharp_metacell .< min_cells_in_metacell)

    sharp_metacell_indices_per_block = [Int[] for _ in 1:n_blocks]
    for (sharp_metacell_index, block_index) in enumerate(block_index_per_sharp_metacell)
        push!(sharp_metacell_indices_per_block[block_index], sharp_metacell_index)
    end

    for sharp_metacell_indices in sharp_metacell_indices_per_block
        if !isempty(sharp_metacell_indices) &&
           all(sharp_metacell_index -> is_dissolved_per_sharp_metacell[sharp_metacell_index], sharp_metacell_indices)
            surviving_sharp_metacell_index =
                sharp_metacell_indices[argmax(n_kept_cells_per_sharp_metacell[sharp_metacell_indices])]
            is_dissolved_per_sharp_metacell[surviving_sharp_metacell_index] = false
        end
    end

    return is_dissolved_per_sharp_metacell
end

# Relocate all the cells of each dissolved metacell to the nearest surviving metacell of its block, in the
# environment-normalized module z-score space K-means clustered the block's cells in (recomputed here via the shared
# `setup_z_score_per_found_module_per_region_cell!`). A surviving metacell's center is the mean z-score of the cells it
# keeps (its outliers are exactly the cells that do not belong in it, so they do not contribute), and it is updated as it
# absorbs cells. Compacts `cells_of_sharp_metacells` and the parallel `block_index_per_sharp_metacell` to the surviving
# metacells.
function relocate_dissolved_metacell_cells!(;
    cells_of_sharp_metacells::Vector{Vector{Int}},
    block_index_per_sharp_metacell::Vector{Int},
    is_dissolved_per_sharp_metacell::BitVector,
    is_outlier_per_cell::Union{AbstractVector{Bool}, BitVector},
    block_index_per_cell::AbstractVector{<:Integer},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    mean_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    n_cells::Integer,
    n_genes::Integer,
    n_modules::Integer,
    n_blocks::Integer,
)::Nothing
    n_sharp_metacells = length(cells_of_sharp_metacells)

    sharp_metacell_index_per_cell = zeros(Int, n_cells)
    for (sharp_metacell_index, cell_indices) in enumerate(cells_of_sharp_metacells)
        sharp_metacell_index_per_cell[cell_indices] .= sharp_metacell_index
    end

    sharp_metacell_indices_per_block = [Int[] for _ in 1:n_blocks]
    for (sharp_metacell_index, block_index) in enumerate(block_index_per_sharp_metacell)
        push!(sharp_metacell_indices_per_block[block_index], sharp_metacell_index)
    end
    has_dissolved_per_block = [
        any(sharp_metacell_index -> is_dissolved_per_sharp_metacell[sharp_metacell_index], sharp_metacell_indices)
        for sharp_metacell_indices in sharp_metacell_indices_per_block
    ]

    n_block_cells_per_block = [count(==(block_index), block_index_per_cell) for block_index in 1:n_blocks]
    weight_per_block =
        [has_dissolved_per_block[block_index] ? n_block_cells_per_block[block_index] : 0 for block_index in 1:n_blocks]
    max_n_block_cells = maximum(weight_per_block)

    z_score_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(z_score_accumulator_pool, Vector{Float32}(undef, max_n_block_cells))
    end
    region_position_per_cell_per_thread = [zeros(Int, n_cells) for _ in 1:maxthreadid()]

    # Parallel over the blocks that dissolved a metacell (each writes only its own cells, so the shared
    # `sharp_metacell_index_per_cell` is not contended). The inner z-score fill is nested-parallel over the block's
    # found modules.
    parallel_loop_wo_rng(
        1:n_blocks;
        weights = weight_per_block,
        name = "relocate_dissolved_metacell_cells",
        progress = DebugProgress(sum(weight_per_block); group = :mcs_loops, desc = "relocate_dissolved_metacell_cells"),
    ) do block_index
        if !has_dissolved_per_block[block_index]
            return nothing
        end

        block_cell_indices = findall(==(block_index), block_index_per_cell)
        n_block_modules = count(@view is_found_per_module_per_block[:, block_index])

        is_found_per_module = falses(n_modules)
        z_score_per_found_module_per_block_cell = Matrix{Float32}(undef, n_block_modules, length(block_cell_indices))
        setup_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_block_cell,
            is_found_per_module,
            block_index,
            indices_of_region_cells = block_cell_indices,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            z_score_accumulator_pool,
            region_position_per_cell_per_thread,
            n_genes,
            n_modules,
        )

        # The block's metacells, as clusters of its cells - `relocate_dissolved_points!` works in the block's local
        # numbering, so map back to the global metacell indices when writing the results.
        sharp_metacell_indices = sharp_metacell_indices_per_block[block_index]
        local_index_per_sharp_metacell = Dict(
            sharp_metacell_index => local_index for
            (local_index, sharp_metacell_index) in enumerate(sharp_metacell_indices)
        )
        local_index_per_block_cell = [
            local_index_per_sharp_metacell[sharp_metacell_index_per_cell[cell_index]] for
            cell_index in block_cell_indices
        ]

        relocate_dissolved_points!(;
            cluster_index_per_point = local_index_per_block_cell,
            value_per_dim_per_point = z_score_per_found_module_per_block_cell,
            is_dissolved_per_cluster = is_dissolved_per_sharp_metacell[sharp_metacell_indices],
            is_center_per_point = [!is_outlier_per_cell[cell_index] for cell_index in block_cell_indices],
        )

        for (block_cell_position, cell_index) in enumerate(block_cell_indices)
            sharp_metacell_index_per_cell[cell_index] =
                sharp_metacell_indices[local_index_per_block_cell[block_cell_position]]
        end
        return nothing
    end

    # Compact the metacells, keeping the cells of each in ascending order.
    new_sharp_metacell_index_per_sharp_metacell = zeros(Int, n_sharp_metacells)
    n_new_sharp_metacells = 0
    for sharp_metacell_index in 1:n_sharp_metacells
        if !is_dissolved_per_sharp_metacell[sharp_metacell_index]
            n_new_sharp_metacells += 1
            new_sharp_metacell_index_per_sharp_metacell[sharp_metacell_index] = n_new_sharp_metacells
        end
    end

    cells_of_new_sharp_metacells = [Int[] for _ in 1:n_new_sharp_metacells]
    for cell_index in 1:n_cells
        sharp_metacell_index = sharp_metacell_index_per_cell[cell_index]
        if sharp_metacell_index > 0
            push!(
                cells_of_new_sharp_metacells[new_sharp_metacell_index_per_sharp_metacell[sharp_metacell_index]],
                cell_index,
            )
        end
    end
    block_index_per_new_sharp_metacell = block_index_per_sharp_metacell[.!is_dissolved_per_sharp_metacell]

    resize!(cells_of_sharp_metacells, n_new_sharp_metacells)
    cells_of_sharp_metacells .= cells_of_new_sharp_metacells
    resize!(block_index_per_sharp_metacell, n_new_sharp_metacells)
    block_index_per_sharp_metacell .= block_index_per_new_sharp_metacell

    @debug (
        "Dissolved metacells: $(n_sharp_metacells - n_new_sharp_metacells)" *
        " Remaining metacells: $(n_new_sharp_metacells)"
    ) _group = :mcs_results
    return nothing
end

# Build one base block's baseline `GroupedSeriesCorrelations{Float32}`. Each group corresponds to one block whose
# cells contribute to the base block's neighborhood; the points within a group are the contributing cells laid out
# contiguously, preserving the per-block cell order. Per (point, base-block series) we store the cell's
# `log2(correlated_UMIs / total_UMIs + reg)` as the fixed side and the punctuated baseline-metacell
# `log2((baseline_metacell_friend_UMIs - cell_friend_UMIs) / (baseline_metacell_total - cell_total) + reg)` as the
# variable side - with the mask `is_active_per_point` set false (and the variable value zeroed) when the cell's
# baseline metacell is missing or the punctuation denominator is non-positive. `baseline_metacell_index_per_cell` and
# the per-baseline-metacell aggregates (`total_UMIs_per_baseline_metacell`,
# `UMIs_per_gene_per_baseline_metacell`) encode "Phase 1's choice of K per block" - the baseline candidate's cluster
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
    UMIs_per_gene_per_baseline_metacell::AbstractMatrix{<:Integer},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    correlated_gene_indices::AbstractVector{<:Integer},
    friend_gene_indices::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat,
)::Tuple{GroupedSeriesCorrelations{Float32}, Vector{Int}}
    n_series = length(correlated_gene_indices)
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
        baseline_metacell_total =
            baseline_metacell_index > 0 ? Int(total_UMIs_per_baseline_metacell[baseline_metacell_index]) : 0
        is_active = baseline_metacell_index > 0 && baseline_metacell_total > cell_total
        is_active_per_point[point_index] = is_active

        for series_index in 1:n_series
            correlated_gene_index = correlated_gene_indices[series_index]
            friend_gene_index = friend_gene_indices[series_index]

            cell_correlated_UMIs = Int(UMIs_per_cell_per_gene[cell_index, correlated_gene_index])
            fixed_per_point_per_series[point_index, series_index] =
                log2(Float32(cell_correlated_UMIs) / Float32(cell_total) + regularization)

            if is_active
                baseline_metacell_friend_UMIs =
                    Int(UMIs_per_gene_per_baseline_metacell[baseline_metacell_index, friend_gene_index])
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
# (`total_UMIs_per_baseline_metacell` and `UMIs_per_gene_per_baseline_metacell`), the per-cell
# `baseline_metacell_index_per_cell` (the cell's baseline metacell index in the global array, or 0 if the cell was not
# clustered), and `first_baseline_metacell_per_block` indicating where each block's baseline
# metacells begin in the global array. For a cell in some block whose Phase 1 cluster is `cluster_index`, the
# baseline metacell index is `first_baseline_metacell_per_block[block_index] + cluster_index - 1`. Blocks with
# `local_clusters === nothing` contribute no baseline metacells.
function compute_baseline_metacell_aggregates(
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    n_genes::Integer,
    n_cells::Integer,
)::Tuple{Vector{UInt32}, Matrix{UInt32}, Vector{UInt32}, Vector{Int}}
    n_blocks = length(local_clusters_per_block)

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
    UMIs_per_gene_per_baseline_metacell = Matrix{UInt32}(undef, n_baseline_metacells, n_genes)
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
            baseline_metacell_index_per_cell[cell_index] = UInt32(baseline_metacell_index)
        end
    end

    n_cells_per_baseline_metacell = length.(cell_indices_per_baseline_metacell)
    parallel_loop_wo_rng(
        1:n_baseline_metacells;
        policy = :greedy,
        weights = n_cells_per_baseline_metacell,
        name = "compute_baseline_metacell_aggregates",
        progress = DebugProgress(
            sum(n_cells_per_baseline_metacell);
            group = :mcs_loops,
            desc = "baseline_metacell_aggregates",
        ),
    ) do baseline_metacell_index
        cell_indices = cell_indices_per_baseline_metacell[baseline_metacell_index]
        baseline_metacell_total = zero(UInt32)
        @inbounds for cell_index in cell_indices
            baseline_metacell_total += UInt32(total_UMIs_per_cell[cell_index])
        end
        total_UMIs_per_baseline_metacell[baseline_metacell_index] = baseline_metacell_total
        @inbounds for gene_index in 1:n_genes
            gene_sum = zero(UInt32)
            for cell_index in cell_indices
                gene_sum += UInt32(UMIs_per_cell_per_gene[cell_index, gene_index])
            end
            UMIs_per_gene_per_baseline_metacell[baseline_metacell_index, gene_index] = gene_sum
        end
        return nothing
    end

    return (
        total_UMIs_per_baseline_metacell,
        UMIs_per_gene_per_baseline_metacell,
        baseline_metacell_index_per_cell,
        first_baseline_metacell_per_block,
    )
end

# Per walkable block, precompute the indirection vectors the delta-correlation scoring will need: which base blocks'
# correlations the block can affect, the block's friend-gene subspace (column axis for the per-(block, K-candidate)
# log-fill cache), and per (block, affected base block) the (column position in the block's friend subspace, the
# block's `block_cell_position` per point) maps.
#   * `affected_base_block_indices_per_walkable_block[walkable_position]`: list of base blocks whose neighborhoods overlap any of the
#     block's cells.
#   * `gene_index_per_friend_per_walkable_block[walkable_position]`: the block's friend-gene subspace (global gene index per
#     column).
#   * `friend_position_per_series_per_base_block_per_walkable_block[walkable_position][base_block_index]`: per (block, base block), the
#     block's friend column for each of the base block's correlated-friend series.
#   * `block_cell_position_per_point_per_base_block_per_walkable_block[walkable_position][base_block_index]`: per (block, base block), the block's
#     `block_cell_position`s of the cells in the base block's neighborhood, in the block's canonical block-cell order.
#     Row indirection for the indirect-gather query.
function precompute_walkable_indirection(
    walkable_block_indices::AbstractVector{<:Integer},
    block_cell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    base_block_index_per_cell::AbstractVector{<:Integer},
    is_in_base_neighborhood_per_other_base_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    friend_gene_indices_per_base_block::AbstractVector{<:AbstractVector{<:Integer}},
    n_base_blocks::Integer,
    n_genes::Integer,
)::Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Dict{Int, Vector{Int}}}, Vector{Dict{Int, Vector{Int}}}}
    n_walkable_blocks = length(walkable_block_indices)

    affected_base_block_indices_per_walkable_block = Vector{Vector{Int}}(undef, n_walkable_blocks)
    gene_index_per_friend_per_walkable_block = Vector{Vector{Int}}(undef, n_walkable_blocks)
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
        weights = n_block_cells_per_walkable_block,
        name = "precompute_walkable_indirection",
        progress = DebugProgress(
            sum(n_block_cells_per_walkable_block);
            group = :mcs_loops,
            desc = "precompute_walkable_indirection",
        ),
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
        gene_index_per_friend = Int[]
        for gene_index in 1:n_genes
            if is_friend_gene_for_walkable_block[gene_index]
                push!(gene_index_per_friend, gene_index)
                friend_position_per_gene[gene_index] = length(gene_index_per_friend)
            end
        end
        gene_index_per_friend_per_walkable_block[walkable_position] = gene_index_per_friend

        friend_position_per_series_per_base_block = Dict{Int, Vector{Int}}()
        block_cell_position_per_point_per_base_block = Dict{Int, Vector{Int}}()
        for base_block_index in affected_base_block_indices
            friend_position_per_series_per_base_block[base_block_index] = [
                friend_position_per_gene[friend_gene_index] for
                friend_gene_index in friend_gene_indices_per_base_block[base_block_index]
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
        gene_index_per_friend_per_walkable_block,
        friend_position_per_series_per_base_block_per_walkable_block,
        block_cell_position_per_point_per_base_block_per_walkable_block,
    )
end

# Fill `z_score_per_found_module_per_region_cell` (rows = the block's found modules in `findall(is_found_per_module)`
# order, columns = the given region cells) with the environment-normalized z-score of each found module's linear
# fraction - the space cells are clustered in. Fills the caller's `is_found_per_module` buffer from the block, and
# returns the per-module gene index lists (the inverted `module_index_per_gene`, reused by callers). This is the shared
# setup around [`compute_z_score_per_found_module_per_region_cell!`](@ref); the caller owns (and sizes) the reused
# buffers.
function setup_z_score_per_found_module_per_region_cell!(;
    z_score_per_found_module_per_region_cell::AbstractMatrix{<:AbstractFloat},
    is_found_per_module::BitVector,
    block_index::Integer,
    indices_of_region_cells::AbstractVector{<:Integer},
    is_found_per_module_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    mean_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    z_score_accumulator_pool::Channel{Vector{Float32}},
    region_position_per_cell_per_thread::AbstractVector{<:AbstractVector{<:Integer}},
    n_genes::Integer,
    n_modules::Integer,
)::Vector{Vector{Int}}
    is_found_per_module .= @view is_found_per_module_per_block[:, block_index]

    # Invert `module_index_per_gene` once per block: parallel-over-found-modules in compute_z_score reads the list
    # directly instead of rebuilding a per-call `is_gene_in_module` BitVector.
    @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
    gene_indices_per_module = [Int[] for _ in 1:n_modules]
    for gene_index in 1:n_genes
        module_index = module_index_per_gene[gene_index]
        if module_index > 0
            push!(gene_indices_per_module[module_index], gene_index)
        end
    end

    @views mean_linear_fraction_in_environment_cells_per_module =
        mean_linear_fraction_in_environment_cells_per_module_per_block[:, block_index]
    @views std_linear_fraction_in_environment_cells_per_module =
        std_linear_fraction_in_environment_cells_per_module_per_block[:, block_index]

    compute_z_score_per_found_module_per_region_cell!(;
        z_score_per_found_module_per_region_cell,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        is_found_per_module,
        indices_of_region_cells,
        gene_indices_per_module,
        z_score_accumulator_pool,
        region_position_per_cell_per_thread,
        mean_linear_fraction_in_environment_cells_per_module,
        std_linear_fraction_in_environment_cells_per_module,
        std_fraction_regularization = std_fraction_regularization_per_block[block_index],
    )

    return gene_indices_per_module
end

function compute_z_score_per_found_module_per_region_cell!(;
    z_score_per_found_module_per_region_cell::AbstractMatrix{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    is_found_per_module::BitVector,
    indices_of_region_cells::AbstractVector{<:Integer},
    gene_indices_per_module::AbstractVector{<:AbstractVector{<:Integer}},
    z_score_accumulator_pool::Channel{Vector{Float32}},
    region_position_per_cell_per_thread::AbstractVector{<:AbstractVector{<:Integer}},
    mean_linear_fraction_in_environment_cells_per_module::AbstractVector{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module::AbstractVector{<:AbstractFloat},
    std_fraction_regularization::AbstractFloat,
)::Nothing
    n_region_cells = length(indices_of_region_cells)
    found_module_indices = findall(is_found_per_module)
    n_found_modules = length(found_module_indices)

    # A module's UMIs are accumulated per region cell by walking the nonzeros of each of its gene columns (no random
    # `[cell, gene]` access); `region_position_per_cell` maps each cell to its position among the region cells (0
    # otherwise). Owned by this outer (block) task and captured read-only by the nested found-module loop, so it is
    # indexed by the outer task's thread and reset once that loop is done.
    @assert issparse(UMIs_per_cell_per_gene)
    sparse_UMIs_per_cell_per_gene = mutable_array(UMIs_per_cell_per_gene)::SparseMatrixCSC
    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)
    region_position_per_cell = region_position_per_cell_per_thread[threadid()]
    for (region_cell_position, cell_index) in enumerate(indices_of_region_cells)
        region_position_per_cell[cell_index] = region_cell_position
    end

    # Parallel over found modules - each fills its own row of the z-score matrix (disjoint writes). Per task, take an
    # accumulator vector from the pool (contiguous Float32, sized to max region cells), do gene-outer / cell-inner
    # `@turbo` accumulation into it (each UMI column visited once and amortized across all region cells), then
    # normalize + write the strided row into the z-score matrix and return the accumulator to the pool. Weight
    # dispatch by each found module's gene count so the biggest modules go first.
    n_genes_per_found_module = [length(gene_indices_per_module[found_module_indices[k]]) for k in 1:n_found_modules]
    parallel_loop_wo_rng(
        1:n_found_modules;
        nested = true,
        policy = :greedy,
        weights = n_genes_per_found_module,
        name = "compute_z_score_per_found_module_per_region_cell",
    ) do found_module_position
        module_index = found_module_indices[found_module_position]
        gene_indices_in_module = gene_indices_per_module[module_index]
        accumulator = take!(z_score_accumulator_pool)
        try
            @views fill!(accumulator[1:n_region_cells], 0.0f0)
            # Walk only the nonzeros of the module's gene columns, mapping each to its region cell; the downstream
            # z-score normalization below operates on dense buffers and stays `@turbo`.
            @inbounds for gene_index in gene_indices_in_module
                for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
                    region_position = region_position_per_cell[row_index_per_stored[stored_index]]
                    if region_position > 0
                        accumulator[region_position] += UMIs_per_stored[stored_index]
                    end
                end
            end
            mean_linear_fraction = mean_linear_fraction_in_environment_cells_per_module[module_index]
            std_linear_fraction = std_linear_fraction_in_environment_cells_per_module[module_index]
            @assert std_linear_fraction > 0
            @check_turbo_vector(accumulator)
            @check_turbo_vector(indices_of_region_cells)
            @check_turbo_vector(total_UMIs_per_cell)
            @check_turbo_matrix(z_score_per_found_module_per_region_cell)
            @turbo for region_cell_position in 1:n_region_cells
                cell_index = indices_of_region_cells[region_cell_position]
                delta_fraction =
                    accumulator[region_cell_position] / total_UMIs_per_cell[cell_index] - mean_linear_fraction
                z_score_per_found_module_per_region_cell[found_module_position, region_cell_position] =
                    (delta_fraction + copysign(std_fraction_regularization, delta_fraction)) /
                    (std_linear_fraction + std_fraction_regularization)
            end
        finally
            put!(z_score_accumulator_pool, accumulator)
        end
        return nothing
    end

    for cell_index in indices_of_region_cells
        region_position_per_cell[cell_index] = 0
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
# it in solution.cells_dispersion_per_cluster. `dispersion_scratches` is the per-call scratch the caller obtained from
# the `Channel{DispersionScratches}` pool (held for the duration of this call and released by the caller).
function compute_cluster_dispersion!(
    solution::KmeansSolution,
    cluster_index::Integer,
    dispersion_scratches::DispersionScratches,
    dispersion_context::DispersionContext,
)::Nothing
    n_block_cells = length(dispersion_context.block_cell_indices)
    n_cluster_cells = 0
    @inbounds for point_index in 1:n_block_cells
        if solution.assignments[point_index] == cluster_index
            n_cluster_cells += 1
            dispersion_scratches.cluster_cell_indices[n_cluster_cells] =
                dispersion_context.block_cell_indices[point_index]
        end
    end

    @views indices_of_cluster_cells = dispersion_scratches.cluster_cell_indices[1:n_cluster_cells]
    dispersion = maximal_cells_dispersion_of_modules!(;
        cells_dispersion_per_module = dispersion_scratches.cells_dispersion_per_module,
        mean_normalized_per_module = dispersion_scratches.mean_normalized_per_module,
        indices_of_cells = indices_of_cluster_cells,
        UMIs_per_cell_per_gene = dispersion_context.UMIs_per_cell_per_gene,
        total_UMIs_per_cell = dispersion_context.total_UMIs_per_cell,
        is_found_per_module = dispersion_context.is_found_per_module,
        gene_indices_per_module = dispersion_context.gene_indices_per_module,
        total_UMIs_per_max_cells = dispersion_scratches.total_UMIs_per_cluster_cell,
        normalized_factor_per_max_cells = dispersion_scratches.normalized_factor_per_cluster_cell,
        normalized_UMIs_quantile = dispersion_context.normalized_UMIs_quantile,
        min_module_UMIs = dispersion_context.min_module_UMIs,
    )
    # `dispersion` is the maximal spread across the block's found modules, and can legitimately be zero: for a cluster
    # too small to estimate spread, or for a block with a single found module where a large cluster of cells that do
    # not express it has no other module to carry a positive value. Enforce positivity only outside these cases.
    n_found_modules = count(dispersion_context.is_found_per_module)
    @assert dispersion > 0 || n_cluster_cells < dispersion_context.min_cluster_size || n_found_modules <= 1
    solution.cells_dispersion_per_cluster[cluster_index] = dispersion
    return nothing
end

# Recompute the maximal cells_dispersion for every active cluster in the solution. Parallel over clusters with a
# `:greedy` nested loop (the caller - Phase 1 - is itself running under the top-level `parallel_loop_with_rng`); each
# inner task `take!`s a `DispersionScratches` from the pool, runs the per-cluster compute (no yield), and `put!`s it
# back. The pool size of `nthreads()` is sufficient because every running task holds at most one scratch.
function compute_all_dispersions!(
    solution::KmeansSolution,
    dispersion_scratches_pool::Channel{DispersionScratches},
    dispersion_context::DispersionContext,
)::Nothing
    parallel_loop_wo_rng(
        1:solution.k;
        nested = true,
        policy = :greedy,
        name = "compute_all_dispersions",
    ) do cluster_index
        dispersion_scratches = take!(dispersion_scratches_pool)
        try
            compute_cluster_dispersion!(solution, cluster_index, dispersion_scratches, dispersion_context)
        finally
            put!(dispersion_scratches_pool, dispersion_scratches)
        end
        return nothing
    end
    return nothing
end

# Incremental dispersion update against a reference solution that shares the same K (typical caller: post-split rerun
# of K-means with the split solution's centers as initial centers, so cluster identity is preserved by Clustering.jl
# and most cells stay in their original cluster). Per point, if the new assignment differs from the reference, both
# cluster indices involved are flagged as `changed`. Changed clusters are re-dispersed (parallel over them, take/put
# from the `DispersionScratches` pool); unchanged clusters copy their dispersion from the reference. The
# change-detection sweep is serial and writes to a locally-allocated `is_changed_per_cluster` (small `BitVector`); the
# `BitVector` is read-only inside the parallel inner tasks.
function compute_changed_dispersions!(
    solution::KmeansSolution,
    reference_assignments::AbstractVector{<:Integer},
    reference_cells_dispersion_per_cluster::AbstractVector{<:AbstractFloat},
    dispersion_scratches_pool::Channel{DispersionScratches},
    dispersion_context::DispersionContext,
)::Nothing
    is_changed_per_cluster = falses(solution.k)
    @inbounds for point_index in 1:solution.n_points
        new_cluster = solution.assignments[point_index]
        old_cluster = reference_assignments[point_index]
        if new_cluster != old_cluster
            is_changed_per_cluster[new_cluster] = true
            is_changed_per_cluster[old_cluster] = true
        end
    end
    parallel_loop_wo_rng(
        1:solution.k;
        nested = true,
        policy = :greedy,
        name = "compute_changed_dispersions",
    ) do cluster_index
        if is_changed_per_cluster[cluster_index]
            dispersion_scratches = take!(dispersion_scratches_pool)
            try
                compute_cluster_dispersion!(solution, cluster_index, dispersion_scratches, dispersion_context)
            finally
                put!(dispersion_scratches_pool, dispersion_scratches)
            end
        else
            solution.cells_dispersion_per_cluster[cluster_index] = reference_cells_dispersion_per_cluster[cluster_index]
        end
        return nothing
    end
    return nothing
end

# Fill `variable_per_block_cell_per_friend` and `is_active_per_block_cell` with `candidate`'s punctuated
# cluster-mean log-fractions over the walkable block's cells × friend genes; cells in single-cell or zero-punctuated-
# total clusters are masked off and their variable values zeroed. Uses `delta_context`'s per-thread cluster-aggregate
# scratches and reads the walkable block's UMIs out of the pre-filled per-walkable
# `UMI_per_friend_per_block_cell` cache (Float32, friend-position-major) - no scattered indexing into the
# global UMI matrix on each call. Allocation-free; the supplied output buffers may be per-thread (for Phase 2
# evaluation, overwritten per work item) or per-walkable-block (for Pass A, kept alive into Pass B).
function populate_candidate_scratches!(
    candidate::SolutionCandidate,
    n_block_cells::Integer,
    delta_context::DeltaCorrelationContext,
    UMI_per_friend_per_block_cell::AbstractMatrix{Float32},
    variable_per_block_cell_per_friend::AbstractMatrix{Float32},
    is_active_per_block_cell::AbstractVector{Bool},
)::Nothing
    n_friends = length(delta_context.gene_index_per_friend)
    k = candidate.k

    total_UMIs_per_max_cluster = delta_context.total_UMIs_per_max_cluster
    UMIs_per_friend_per_max_cluster = delta_context.UMIs_per_friend_per_max_cluster
    punctuated_total_per_block_cell = delta_context.punctuated_total_per_block_cell
    @views fill!(total_UMIs_per_max_cluster[1:k], 0.0)
    @views fill!(UMIs_per_friend_per_max_cluster[1:n_friends, 1:k], 0.0)

    @check_turbo_matrix(UMI_per_friend_per_block_cell)
    @check_turbo_matrix(UMIs_per_friend_per_max_cluster)
    @check_turbo_matrix(variable_per_block_cell_per_friend)

    # Per K-candidate cluster: total UMIs of the block's cells in the cluster + per-friend-gene UMI sums (only over
    # the walkable block's friend-gene subspace, not the full relevant gene set). Inner `friend` loop is contiguous on
    # both `UMI_per_friend_per_block_cell` (friend-major Float32 column read) and `UMIs_per_friend_per_max_cluster`
    # (friend-major Float64 column accumulate) - SIMD via `@turbo`.
    @inbounds for block_cell_position in 1:n_block_cells
        cluster_index = candidate.assignments[block_cell_position]
        cell_index = delta_context.block_cell_indices[block_cell_position]
        total_UMIs_per_max_cluster[cluster_index] += Float64(delta_context.total_UMIs_per_cell[cell_index])
        @turbo for friend_position in 1:n_friends
            UMIs_per_friend_per_max_cluster[friend_position, cluster_index] +=
                Float64(UMI_per_friend_per_block_cell[friend_position, block_cell_position])
        end
    end

    # Pre-pass per block cell: is_active + a safe `punctuated_total` (=1 when inactive, so the inner @turbo loop can
    # divide unconditionally without producing NaN/Inf for inactive cells - the `ifelse` then zeroes those entries).
    @inbounds for block_cell_position in 1:n_block_cells
        cluster_index = candidate.assignments[block_cell_position]
        cell_index = delta_context.block_cell_indices[block_cell_position]
        cell_total = Float64(delta_context.total_UMIs_per_cell[cell_index])
        cluster_total_UMIs = total_UMIs_per_max_cluster[cluster_index]
        punctuated_total = cluster_total_UMIs - cell_total
        is_active = candidate.counts[cluster_index] >= 2 && punctuated_total > 0
        is_active_per_block_cell[block_cell_position] = is_active
        punctuated_total_per_block_cell[block_cell_position] = is_active ? punctuated_total : 1.0
    end

    # Punctuated cluster-mean log-fraction per (block cell, friend gene). Inner `friend` loop reads contiguous friend-
    # major columns and writes strided into `variable_per_block_cell_per_friend` (scatter store) - @turbo SIMDs the
    # branchless arithmetic with `ifelse` selecting 0 for inactive cells.
    regularization = delta_context.gene_fraction_regularization
    @inbounds for block_cell_position in 1:n_block_cells
        cluster_index = candidate.assignments[block_cell_position]
        is_active = is_active_per_block_cell[block_cell_position]
        punctuated_total = punctuated_total_per_block_cell[block_cell_position]
        @turbo for friend_position in 1:n_friends
            cell_friend_UMIs = Float64(UMI_per_friend_per_block_cell[friend_position, block_cell_position])
            cluster_friend_UMIs = UMIs_per_friend_per_max_cluster[friend_position, cluster_index]
            punctuated_friend = cluster_friend_UMIs - cell_friend_UMIs
            result = log2(Float32(punctuated_friend / punctuated_total) + regularization)
            variable_per_block_cell_per_friend[block_cell_position, friend_position] = ifelse(is_active, result, 0.0f0)
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
    UMI_per_friend_per_block_cell::AbstractMatrix{Float32},
)::Float64
    variable_per_block_cell_per_friend = delta_context.variable_per_block_cell_per_friend
    is_active_per_block_cell = delta_context.is_active_per_block_cell
    populate_candidate_scratches!(
        candidate,
        n_block_cells,
        delta_context,
        UMI_per_friend_per_block_cell,
        variable_per_block_cell_per_friend,
        is_active_per_block_cell,
    )
    return query_delta_correlation_from_scratches(delta_context)
end

# Query half of `compute_delta_correlation`: reads `delta_context`'s per-thread variable + mask scratches (assumed
# already populated for some candidate via `populate_candidate_scratches!`) plus the per-walkable indirections, sums
# the per-affected-base-block delta against the cached baseline mean correlations. Allocation-free.
function query_delta_correlation_from_scratches(delta_context::DeltaCorrelationContext)::Float64
    variable_per_block_cell_per_friend = delta_context.variable_per_block_cell_per_friend
    is_active_per_block_cell = delta_context.is_active_per_block_cell

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
            variable_per_block_cell_per_friend,
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

function update_size_statistics!(
    solution::KmeansSolution,
    min_cluster_size::Real,
    max_cluster_size::Real,
    min_cluster_weight::Real,
    dispersion_context::DispersionContext,
)::Nothing
    max_cluster_size = min(max_cluster_size, 2 * solution.n_points / solution.k)
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
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
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
        return kmeans_in_rounds(split_points, 2; buffer_pool = kmeans_buffer_pool, rounds = kmeans_rounds, rng)
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
    # Take a `DispersionScratches` from the pool for both per-cluster dispersion calls (serial, no yield in between).
    dispersion_scratches = take!(dispersion_scratches_pool)
    try
        compute_cluster_dispersion!(candidate_solution, old_cluster_index, dispersion_scratches, dispersion_context)
        compute_cluster_dispersion!(candidate_solution, new_cluster_index, dispersion_scratches, dispersion_context)
    finally
        put!(dispersion_scratches_pool, dispersion_scratches)
    end
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
    max_cluster_size::Real,
    min_cluster_weight::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::Integer where {T <: AbstractFloat}
    k = current_solution.k

    min_splittable_size = min(2 * min_cluster_size, max_cluster_size)
    min_splittable_weight = 2 * min_cluster_weight

    target_cluster = 0
    target_size = 2 * current_solution.n_points / current_solution.k
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
        kmeans_buffer_pool,
        dispersion_scratches_pool,
        dispersion_context,
        rng,
    )

    update_size_statistics!(
        candidate_solution,
        min_cluster_size,
        max_cluster_size,
        min_cluster_weight,
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
    max_cluster_size::Real,
    min_cluster_weight::Real,
    kmeans_rounds::Integer,
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
    dispersion_context::DispersionContext,
    rng::AbstractRNG;
    reference_assignments::Maybe{AbstractVector{<:Integer}} = nothing,
    reference_cells_dispersion_per_cluster::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
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
            buffer_pool = kmeans_buffer_pool,
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
    # Incremental dispersion if a same-K reference is supplied (post-split rerun: initial centers came from that
    # reference, so cluster identity is preserved and most cells should stay put). Otherwise full recompute.
    can_incremental =
        reference_assignments !== nothing &&
        reference_cells_dispersion_per_cluster !== nothing &&
        length(reference_assignments) >= rerun_solution.n_points &&
        length(reference_cells_dispersion_per_cluster) >= new_k
    if can_incremental
        compute_changed_dispersions!(
            rerun_solution,
            reference_assignments::AbstractVector{<:Integer},
            reference_cells_dispersion_per_cluster::AbstractVector{<:AbstractFloat},
            dispersion_scratches_pool,
            dispersion_context,
        )
    else
        compute_all_dispersions!(rerun_solution, dispersion_scratches_pool, dispersion_context)
    end
    update_size_statistics!(rerun_solution, min_cluster_size, max_cluster_size, min_cluster_weight, dispersion_context)
    if is_perfect(rerun_solution)
        push!(
            perfect_candidates,
            build_candidate_from_solution(
                rerun_solution,
                values_of_points,
                weight_per_point,
                min_cluster_size,
                min_cluster_weight,
            ),
        )
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
    max_cluster_size::Real,
    min_cluster_weight::Real,
    kmeans_rounds::Integer,
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
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
                max_cluster_size,
                min_cluster_weight,
                kmeans_rounds,
                kmeans_buffer_pool,
                dispersion_scratches_pool,
                dispersion_context,
                rng,
            )

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
# the per-thread max-k buffers), so that the candidate survives subsequent `rerun_kmeans!` overwrites of the working
# buffers. The solution's too-small clusters are dissolved into the surviving ones, so a candidate never contains a
# too-small cluster. A perfect solution has none to dissolve; only the imperfect compromise solution does.
function build_candidate_from_solution(
    solution::KmeansSolution,
    values_of_points::AbstractMatrix{<:AbstractFloat},
    weight_per_point::AbstractVector{<:Real},
    min_cluster_size::Real,
    min_cluster_weight::Real,
)::SolutionCandidate
    k = solution.k
    n_points = solution.n_points

    @views is_dissolved_per_cluster = pick_dissolved_clusters(
        solution.counts[1:k],
        solution.weight_per_cluster[1:k],
        min_cluster_size,
        min_cluster_weight,
    )
    if !any(is_dissolved_per_cluster)
        return SolutionCandidate(
            k,
            copy(@view(solution.counts[1:k])),
            copy(@view(solution.weight_per_cluster[1:k])),
            copy(@view(solution.assignments[1:n_points])),
        )
    end

    cluster_index_per_point = copy(@view(solution.assignments[1:n_points]))
    relocate_dissolved_points!(;
        cluster_index_per_point,
        value_per_dim_per_point = values_of_points,
        is_dissolved_per_cluster,
        is_center_per_point = trues(n_points),
    )

    # Renumber the assignments to the surviving clusters, and recompute their sizes and weights.
    new_cluster_index_per_cluster = zeros(Int, k)
    n_new_clusters = 0
    for cluster_index in 1:k
        if !is_dissolved_per_cluster[cluster_index]
            n_new_clusters += 1
            new_cluster_index_per_cluster[cluster_index] = n_new_clusters
        end
    end

    counts = zeros(Int, n_new_clusters)
    weight_per_cluster = zeros(Float64, n_new_clusters)
    for point_index in 1:n_points
        new_cluster_index = new_cluster_index_per_cluster[cluster_index_per_point[point_index]]
        cluster_index_per_point[point_index] = new_cluster_index
        counts[new_cluster_index] += 1
        weight_per_cluster[new_cluster_index] += weight_per_point[point_index]
    end

    return SolutionCandidate(n_new_clusters, counts, weight_per_cluster, cluster_index_per_point)
end

# The clusters to dissolve: the too-small ones, with fewer than `min_cluster_size` points or less than
# `min_cluster_weight` weight. The last cluster is never dissolved - there would be nowhere to relocate its points to -
# so if all the clusters are too small, the largest one survives to absorb the rest.
function pick_dissolved_clusters(
    counts::AbstractVector{<:Integer},
    weight_per_cluster::AbstractVector{<:Real},
    min_cluster_size::Real,
    min_cluster_weight::Real,
)::BitVector
    is_dissolved_per_cluster = (counts .< min_cluster_size) .| (weight_per_cluster .< min_cluster_weight)
    if all(is_dissolved_per_cluster)
        is_dissolved_per_cluster[argmax(counts)] = false
    end
    return is_dissolved_per_cluster
end

# Relocate every point of the dissolved clusters into the nearest surviving cluster, in the space
# `value_per_dim_per_point` (the same space K-means clustered the points in). A cluster's center is the mean value of its
# center points, kept as a Float64 sum and a count, so absorbing a point immediately updates the center the next point is
# compared against. The points are visited in order, so the result does not depend on the number of threads. Updates
# `cluster_index_per_point` in place; the caller recomputes whatever per-cluster aggregates it needs.
function relocate_dissolved_points!(;
    cluster_index_per_point::AbstractVector{<:Integer},
    value_per_dim_per_point::AbstractMatrix{<:AbstractFloat},
    is_dissolved_per_cluster::Union{AbstractVector{Bool}, BitVector},
    is_center_per_point::Union{AbstractVector{Bool}, BitVector},
)::Nothing
    n_dims, n_points = size(value_per_dim_per_point)
    n_clusters = length(is_dissolved_per_cluster)

    sum_per_dim_per_cluster = zeros(Float64, n_dims, n_clusters)
    n_center_points_per_cluster = zeros(Int, n_clusters)
    for point_index in 1:n_points
        cluster_index = cluster_index_per_point[point_index]
        if is_center_per_point[point_index] && !is_dissolved_per_cluster[cluster_index]
            @views sum_per_dim_per_cluster[:, cluster_index] .+= value_per_dim_per_point[:, point_index]
            n_center_points_per_cluster[cluster_index] += 1
        end
    end

    for point_index in 1:n_points
        if !is_dissolved_per_cluster[cluster_index_per_point[point_index]]
            continue
        end

        nearest_cluster_index = 0
        nearest_distance = Inf
        for cluster_index in 1:n_clusters
            if is_dissolved_per_cluster[cluster_index] || n_center_points_per_cluster[cluster_index] == 0
                continue
            end
            scale = 1.0 / n_center_points_per_cluster[cluster_index]
            distance = 0.0
            @inbounds for dim_index in 1:n_dims
                difference =
                    value_per_dim_per_point[dim_index, point_index] -
                    sum_per_dim_per_cluster[dim_index, cluster_index] * scale
                distance += difference * difference
            end
            if distance < nearest_distance
                nearest_distance = distance
                nearest_cluster_index = cluster_index
            end
        end
        @assert nearest_cluster_index > 0

        cluster_index_per_point[point_index] = nearest_cluster_index
        @views sum_per_dim_per_cluster[:, nearest_cluster_index] .+= value_per_dim_per_point[:, point_index]
        n_center_points_per_cluster[nearest_cluster_index] += 1
    end

    return nothing
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
    max_cluster_size::Real,
    min_cluster_weight::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
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
            max_cluster_size,
            min_cluster_weight,
            kmeans_rounds,
            sizes_buffers,
            kmeans_buffer_pool,
            dispersion_scratches_pool,
            dispersion_context,
            rng,
        )
        if split_cluster_index == 0
            break  # No more candidate clusters to split.
        end

        if penalty(candidate_solution) < penalty(best_solution)
            copy_solution!(best_solution, candidate_solution)
            if is_perfect(best_solution)
                push!(
                    perfect_candidates,
                    build_candidate_from_solution(
                        best_solution,
                        values_of_points,
                        weight_per_point,
                        min_cluster_size,
                        min_cluster_weight,
                    ),
                )
            end

            # We can re-K-means based on the split clusters to try and improve the solution. Pass `candidate_solution`'s
            # assignments + cluster dispersions as the reference - K is identical and the initial centers come from the
            # same source, so `rerun_kmeans!` can short-circuit dispersion for clusters whose membership did not move.
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
                max_cluster_size,
                min_cluster_weight,
                kmeans_rounds,
                kmeans_buffer_pool,
                dispersion_scratches_pool,
                dispersion_context,
                rng;
                reference_assignments = @view(candidate_solution.assignments[1:candidate_solution.n_points]),
                reference_cells_dispersion_per_cluster = @view(
                    candidate_solution.cells_dispersion_per_cluster[1:candidate_solution.k]
                ),
            )

            if penalty(current_solution) < penalty(best_solution)
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

function kmeans_with_sizes_candidates(
    values_of_points::AbstractMatrix{T},
    weight_per_point::AbstractVector{<:Real},
    initial_k::Integer;
    max_k::Integer,
    min_cluster_size::Real,
    max_cluster_size::Real,
    min_cluster_weight::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers{T},
    kmeans_buffer_pool::Channel{KMeansBuffers{T}},
    dispersion_scratches_pool::Channel{DispersionScratches},
    dispersion_context::DispersionContext,
    rng::AbstractRNG,
)::SolutionCandidates where {T <: AbstractFloat}
    @assert max_k >= initial_k
    @assert min_cluster_size >= 2
    @assert max_cluster_size >= 2 * min_cluster_size
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
        max_cluster_size,
        min_cluster_weight,
        kmeans_rounds,
        kmeans_buffer_pool,
        dispersion_scratches_pool,
        dispersion_context,
        rng,
    )

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
            max_cluster_size,
            min_cluster_weight,
            kmeans_rounds,
            kmeans_buffer_pool,
            dispersion_scratches_pool,
            dispersion_context,
            rng,
        )
    end

    walk_split!(
        best_solution,
        current_solution,
        candidate_solution,
        perfect_candidates,
        values_of_points,
        weight_per_point,
        max_k,
        min_cluster_size,
        max_cluster_size,
        min_cluster_weight,
        kmeans_rounds,
        sizes_buffers,
        kmeans_buffer_pool,
        dispersion_scratches_pool,
        dispersion_context,
        rng,
    )

    @assert best_solution.is_filled

    if isempty(perfect_candidates)
        # No perfect found - take the best-by-penalty imperfect compromise as the lone candidate. This is the only
        # candidate that can contain too-small clusters, which are dissolved into the surviving ones.
        return SolutionCandidates(
            n_points,
            [
                build_candidate_from_solution(
                    best_solution,
                    values_of_points,
                    weight_per_point,
                    min_cluster_size,
                    min_cluster_weight,
                ),
            ],
        )
    end

    # Sort ascending by K so `candidates[1]` is the minimal-K baseline.
    sort!(perfect_candidates; by = candidate -> candidate.k)
    return SolutionCandidates(n_points, perfect_candidates)
end

# Build a `LocalClusters` from a `SolutionCandidate`. `block_cell_indices` is shared (not copied); the assignment vector
# is owned by the returned `LocalClusters`.
function build_local_clusters_from_candidate(
    candidate::SolutionCandidate,
    block_cell_indices::AbstractVector{<:Integer},
)::LocalClusters
    return LocalClusters(;
        n_clusters = candidate.k,
        block_cell_indices,
        cluster_index_per_block_cell = copy(candidate.assignments),
    )
end

"""
    compute_matrix_of_n_cells_per_prev_block_per_block!(;
        other_daf::DafWriter,
        base_daf::DafReader,
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_n_cells_per_prev_block_per_block`](@ref). This counts, for each pair of a `prev_daf` block
and an `other_daf` block, the cells that are grouped in both. This will also copy [`prev_block_axis`](@ref) from the
`prev_daf` into the `other_daf` if needed.

# Other

$(CONTRACT1)

# Previous

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "other_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        prev_block_axis(GuaranteedOutput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_n_cells_per_prev_block_per_block(CreatedOutput),
    ],
) Contract(;
    name = "prev_daf",
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [vector_of_metacell_per_cell(RequiredInput), vector_of_block_per_metacell(RequiredInput)],
) function compute_matrix_of_n_cells_per_prev_block_per_block!(;
    other_daf::DafWriter,
    prev_daf::DafReader,
    overwrite::Bool = false,
)::Nothing
    @assert axis_vector(other_daf, "cell") == axis_vector(prev_daf, "cell") "the cells differ between `other_daf` and `prev_daf`"

    if has_axis(other_daf, "prev_block")
        @assert axis_vector(other_daf, "prev_block") == axis_vector(prev_daf, "block")
    else
        copy_axis!(source = prev_daf, destination = other_daf, axis = "block", rename = "prev_block"; overwrite)
    end

    n_blocks = axis_length(other_daf, "block")
    n_prev_blocks = axis_length(prev_daf, "block")

    block_index_per_cell = other_daf["@ cell : metacell ?? 0 : block : index"].array
    prev_block_index_per_cell = prev_daf["@ cell : metacell ?? 0 : block : index"].array

    n_cells_per_prev_block_per_block = zeros(UInt32, n_prev_blocks, n_blocks)
    for (prev_block_index, block_index) in zip(prev_block_index_per_cell, block_index_per_cell)
        if prev_block_index > 0 && block_index > 0
            n_cells_per_prev_block_per_block[prev_block_index, block_index] += 1
        end
    end

    set_matrix!(other_daf, "prev_block", "block", "n_cells", bestify(n_cells_per_prev_block_per_block); overwrite)

    return nothing
end

"""
    compute_matrix_of_n_cells_per_prev_block_type_per_block_type!(;
        other_daf::DafWriter,
        prev_daf::DafReader,
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_n_cells_per_prev_block_type_per_block_type`](@ref). This counts, for each pair of a
previous-round block type and an `other_daf` block type, the cells that are of both. The block type of a cell is the type
of the block of the metacell of the cell, so it reflects this repository's blocks rather than any per-cell type
annotation. The (shared) [`type_axis`](@ref) must be identical between the `prev_daf` and the `other_daf`.

# Other

$(CONTRACT1)

# Previous

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "other_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        type_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_block(RequiredInput),
        matrix_of_n_cells_per_prev_block_type_per_block_type(CreatedOutput),
    ],
) Contract(;
    name = "prev_daf",
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        type_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_type_per_block(RequiredInput),
    ],
) function compute_matrix_of_n_cells_per_prev_block_type_per_block_type!(;
    other_daf::DafWriter,
    prev_daf::DafReader,
    overwrite::Bool = false,
)::Nothing
    @assert axis_vector(other_daf, "cell") == axis_vector(prev_daf, "cell") "the cells differ between `other_daf` and `prev_daf`"
    @assert axis_vector(other_daf, "type") == axis_vector(prev_daf, "type") "the types differ between `other_daf` and `prev_daf`"

    n_types = axis_length(other_daf, "type")

    block_type_index_per_cell = other_daf["@ cell : metacell ?? 0 : block : type ?? 0 : index"].array
    prev_block_type_index_per_cell = prev_daf["@ cell : metacell ?? 0 : block : type ?? 0 : index"].array

    n_cells_per_prev_block_type_per_block_type = zeros(UInt32, n_types, n_types)
    for (prev_block_type_index, block_type_index) in zip(prev_block_type_index_per_cell, block_type_index_per_cell)
        if prev_block_type_index > 0 && block_type_index > 0
            n_cells_per_prev_block_type_per_block_type[prev_block_type_index, block_type_index] += 1
        end
    end

    set_matrix!(
        other_daf,
        "type",
        "type",
        "n_cells_by_block",
        bestify(n_cells_per_prev_block_type_per_block_type);
        overwrite,
    )

    return nothing
end

"""
    compute_matrix_of_n_cells_per_prev_metacell_type_per_metacell_type!(;
        other_daf::DafWriter,
        prev_daf::DafReader,
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_n_cells_per_prev_metacell_type_per_metacell_type`](@ref). This counts, for each pair of a
previous-round metacell type and an `other_daf` metacell type, the cells that are of both. The metacell type of a cell is
the type of the metacell of the cell. The (shared) [`type_axis`](@ref) must be identical between the `prev_daf` and the
`other_daf`.

# Other

$(CONTRACT1)

# Previous

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "other_daf",
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), type_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_metacell(RequiredInput),
        matrix_of_n_cells_per_prev_metacell_type_per_metacell_type(CreatedOutput),
    ],
) Contract(;
    name = "prev_daf",
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), type_axis(RequiredInput)],
    data = [vector_of_metacell_per_cell(RequiredInput), vector_of_type_per_metacell(RequiredInput)],
) function compute_matrix_of_n_cells_per_prev_metacell_type_per_metacell_type!(;
    other_daf::DafWriter,
    prev_daf::DafReader,
    overwrite::Bool = false,
)::Nothing
    @assert axis_vector(other_daf, "cell") == axis_vector(prev_daf, "cell") "the cells differ between `other_daf` and `prev_daf`"
    @assert axis_vector(other_daf, "type") == axis_vector(prev_daf, "type") "the types differ between `other_daf` and `prev_daf`"

    n_types = axis_length(other_daf, "type")

    metacell_type_index_per_cell = other_daf["@ cell : metacell ?? 0 : type ?? 0 : index"].array
    prev_metacell_type_index_per_cell = prev_daf["@ cell : metacell ?? 0 : type ?? 0 : index"].array

    n_cells_per_prev_metacell_type_per_metacell_type = zeros(UInt32, n_types, n_types)
    for (prev_metacell_type_index, metacell_type_index) in
        zip(prev_metacell_type_index_per_cell, metacell_type_index_per_cell)
        if prev_metacell_type_index > 0 && metacell_type_index > 0
            n_cells_per_prev_metacell_type_per_metacell_type[prev_metacell_type_index, metacell_type_index] += 1
        end
    end

    set_matrix!(
        other_daf,
        "type",
        "type",
        "n_cells_by_metacell",
        bestify(n_cells_per_prev_metacell_type_per_metacell_type);
        overwrite,
    )

    return nothing
end

# The total weighted crossings of the type flow for a given order, summed over all the consecutive-round transitions.
# Each transition is a list of edges; an edge connects a base type (in the left column) to a type (in the right column)
# and weighs the number of cells with that base type and type. Two edges cross when their endpoints are in the opposite
# order in the two columns, and such a crossing weighs the product of the two edges' cell counts. Both columns of every
# transition use the same order, given by `position_per_type`.
function total_type_flow_crossings(
    position_per_type::AbstractVector{<:Integer},
    base_type_per_edge_per_transition::AbstractVector{Vector{Int}},
    type_per_edge_per_transition::AbstractVector{Vector{Int}},
    n_cells_per_edge_per_transition::AbstractVector{Vector{Float64}},
)::Float64
    total_crossings = 0.0
    @inbounds for transition_index in eachindex(base_type_per_edge_per_transition)
        base_type_per_edge = base_type_per_edge_per_transition[transition_index]
        type_per_edge = type_per_edge_per_transition[transition_index]
        n_cells_per_edge = n_cells_per_edge_per_transition[transition_index]
        n_edges = length(base_type_per_edge)
        for first_edge_index in 1:(n_edges - 1)
            first_left_position = position_per_type[base_type_per_edge[first_edge_index]]
            first_right_position = position_per_type[type_per_edge[first_edge_index]]
            first_n_cells = n_cells_per_edge[first_edge_index]
            for second_edge_index in (first_edge_index + 1):n_edges
                second_left_position = position_per_type[base_type_per_edge[second_edge_index]]
                second_right_position = position_per_type[type_per_edge[second_edge_index]]
                if (first_left_position - second_left_position) * (first_right_position - second_right_position) < 0
                    total_crossings += first_n_cells * n_cells_per_edge[second_edge_index]
                end
            end
        end
    end
    return total_crossings
end

# Fill `position_per_type` with the inverse of the `type_per_position` order (the position holding each type).
function fill_position_per_type!(
    position_per_type::AbstractVector{<:Integer},
    type_per_position::AbstractVector{<:Integer},
)::Nothing
    @inbounds for position in eachindex(type_per_position)
        position_per_type[type_per_position[position]] = position
    end
    return nothing
end

# Fill `target_per_position` with the order obtained from `source_per_position` by moving the type at `from_position` to
# `to_position`, shifting the types in between by one place.
function fill_type_flow_insertion!(
    target_per_position::AbstractVector{<:Integer},
    source_per_position::AbstractVector{<:Integer},
    from_position::Integer,
    to_position::Integer,
)::Nothing
    copyto!(target_per_position, source_per_position)
    moved_type = target_per_position[from_position]
    if from_position < to_position
        for position in from_position:(to_position - 1)
            target_per_position[position] = target_per_position[position + 1]
        end
    else
        for position in from_position:-1:(to_position + 1)
            target_per_position[position] = target_per_position[position - 1]
        end
    end
    target_per_position[to_position] = moved_type
    return nothing
end

# Improve the `type_per_position` order in place by repeatedly applying the best crossing-reducing move (exchanging the
# types at two positions, or moving one type to another position), until reaching a local optimum. Returns the local
# optimum's total weighted crossings. The other vectors are reused scratch buffers.
function hill_climb_type_flow_order!(
    type_per_position::AbstractVector{Int},
    position_per_type::AbstractVector{Int},
    candidate_per_position::AbstractVector{Int},
    best_per_position::AbstractVector{Int},
    base_type_per_edge_per_transition::AbstractVector{Vector{Int}},
    type_per_edge_per_transition::AbstractVector{Vector{Int}},
    n_cells_per_edge_per_transition::AbstractVector{Vector{Float64}},
)::Float64
    n_types = length(type_per_position)
    fill_position_per_type!(position_per_type, type_per_position)
    current_crossings = total_type_flow_crossings(
        position_per_type,
        base_type_per_edge_per_transition,
        type_per_edge_per_transition,
        n_cells_per_edge_per_transition,
    )

    while true
        best_crossings = current_crossings
        found_better = false

        for first_position in 1:(n_types - 1)
            for second_position in (first_position + 1):n_types
                copyto!(candidate_per_position, type_per_position)
                candidate_per_position[first_position], candidate_per_position[second_position] =
                    candidate_per_position[second_position], candidate_per_position[first_position]
                fill_position_per_type!(position_per_type, candidate_per_position)
                crossings = total_type_flow_crossings(
                    position_per_type,
                    base_type_per_edge_per_transition,
                    type_per_edge_per_transition,
                    n_cells_per_edge_per_transition,
                )
                if crossings < best_crossings
                    best_crossings = crossings
                    copyto!(best_per_position, candidate_per_position)
                    found_better = true
                end
            end
        end

        for from_position in 1:n_types
            for to_position in 1:n_types
                if from_position != to_position
                    fill_type_flow_insertion!(candidate_per_position, type_per_position, from_position, to_position)
                    fill_position_per_type!(position_per_type, candidate_per_position)
                    crossings = total_type_flow_crossings(
                        position_per_type,
                        base_type_per_edge_per_transition,
                        type_per_edge_per_transition,
                        n_cells_per_edge_per_transition,
                    )
                    if crossings < best_crossings
                        best_crossings = crossings
                        copyto!(best_per_position, candidate_per_position)
                        found_better = true
                    end
                end
            end
        end

        if found_better
            copyto!(type_per_position, best_per_position)
            current_crossings = best_crossings
        else
            return current_crossings
        end
    end
end

# Compute a global order of the types (the 1-based position of each type) that minimizes the total weighted crossings of
# the type flow across all the consecutive-round transitions. Each transition is described by an `n_cells` matrix whose
# `[base_type, type]` entry counts the cells of that base type (left column) and type (right column). Since minimizing
# crossings is NP-hard, we use random-restart hill climbing using the exact crossing count, seeding one restart with the
# types ordered by their total number of cells.
function compute_global_flow_order_per_type(
    n_cells_per_prev_type_per_type_per_transition::AbstractVector{<:AbstractMatrix{<:Real}};
    restarts::Integer,
    rng::AbstractRNG,
)::Vector{Int}
    @assert restarts >= 1
    @assert !isempty(n_cells_per_prev_type_per_type_per_transition)
    n_types = size(n_cells_per_prev_type_per_type_per_transition[1], 1)
    n_transitions = length(n_cells_per_prev_type_per_type_per_transition)

    base_type_per_edge_per_transition = Vector{Vector{Int}}(undef, n_transitions)
    type_per_edge_per_transition = Vector{Vector{Int}}(undef, n_transitions)
    n_cells_per_edge_per_transition = Vector{Vector{Float64}}(undef, n_transitions)
    total_cells_per_type = zeros(Float64, n_types)
    for transition_index in 1:n_transitions
        n_cells_per_prev_block_type_per_block_type = n_cells_per_prev_type_per_type_per_transition[transition_index]
        @assert size(n_cells_per_prev_block_type_per_block_type) == (n_types, n_types)
        base_type_per_edge = Int[]
        type_per_edge = Int[]
        n_cells_per_edge = Float64[]
        for type in 1:n_types
            for base_type in 1:n_types
                n_cells = Float64(n_cells_per_prev_block_type_per_block_type[base_type, type])
                if n_cells != 0
                    push!(base_type_per_edge, base_type)
                    push!(type_per_edge, type)
                    push!(n_cells_per_edge, n_cells)
                    total_cells_per_type[base_type] += n_cells
                    total_cells_per_type[type] += n_cells
                end
            end
        end
        base_type_per_edge_per_transition[transition_index] = base_type_per_edge
        type_per_edge_per_transition[transition_index] = type_per_edge
        n_cells_per_edge_per_transition[transition_index] = n_cells_per_edge
    end

    seed_type_per_position = sortperm(total_cells_per_type; rev = true)

    crossings_per_restart = Vector{Float64}(undef, restarts)
    type_per_position_per_restart = [Vector{Int}(undef, n_types) for _ in 1:restarts]
    position_per_type_per_restart = [Vector{Int}(undef, n_types) for _ in 1:restarts]
    candidate_per_position_per_restart = [Vector{Int}(undef, n_types) for _ in 1:restarts]
    best_per_position_per_restart = [Vector{Int}(undef, n_types) for _ in 1:restarts]
    parallel_loop_with_rng(1:restarts; rng, name = "compute_global_flow_order_per_type") do restart_index, rng
        type_per_position = type_per_position_per_restart[restart_index]
        if restart_index == 1
            copyto!(type_per_position, seed_type_per_position)
        else
            randperm!(rng, type_per_position)
        end
        crossings_per_restart[restart_index] = hill_climb_type_flow_order!(
            type_per_position,
            position_per_type_per_restart[restart_index],
            candidate_per_position_per_restart[restart_index],
            best_per_position_per_restart[restart_index],
            base_type_per_edge_per_transition,
            type_per_edge_per_transition,
            n_cells_per_edge_per_transition,
        )
        return nothing
    end

    best_type_per_position = type_per_position_per_restart[argmin(crossings_per_restart)]
    return invperm(best_type_per_position)
end

"""
    compute_vector_of_global_flow_order_per_type!(
        final_daf::DafWriter,
        base_daf_per_round::AbstractVector{<:DafReader};
        restarts::Integer = 20,
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = false,
    )::Nothing

Compute and set [`vector_of_global_flow_order_per_type`](@ref), a global order of the types minimizing the total weighted
crossings of the type flow across the sharpening rounds.

The `base_daf_per_round` are the repositories of the earlier rounds (rounds 0 to N-1) and `final_daf` is the last round
(round N); together they form the full sequence of rounds 0 to N, which must all share the same `type` axis. For each
consecutive pair of rounds we read both the [`matrix_of_n_cells_per_prev_block_type_per_block_type`](@ref) and the
[`matrix_of_n_cells_per_prev_metacell_type_per_metacell_type`](@ref) (the type flow from the previous round, by block type
and by metacell type) from the later round's repository, that is, from every repository except the first (round 0 has no
previous round). We minimize the combined (averaged) crossings of both using random-restart hill climbing (`restarts`
restarts, since the problem is NP-hard) and write the resulting order only into `final_daf`.

The contract below is that of `final_daf`. Each of the other repositories satisfies the same contract, except that the two
`n_cells` matrices are absent from the first (round 0), and the [`vector_of_global_flow_order_per_type`](@ref) output is
created only in `final_daf` (the last).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [type_axis(RequiredInput)],
    data = [
        matrix_of_n_cells_per_prev_block_type_per_block_type(RequiredInput),
        matrix_of_n_cells_per_prev_metacell_type_per_metacell_type(RequiredInput),
        vector_of_global_flow_order_per_type(CreatedOutput),
    ],
) function compute_vector_of_global_flow_order_per_type!(
    final_daf::DafWriter,
    base_daf_per_round::AbstractVector{<:DafReader};
    restarts::Integer = 20,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert restarts >= 1
    n_rounds = length(base_daf_per_round)
    @assert n_rounds > 0

    type_names = axis_vector(final_daf, "type")
    for base_daf in base_daf_per_round
        @assert axis_vector(base_daf, "type") == type_names "the types differ between the round repositories"
    end

    # The order minimizes the total crossings over both the by-block and by-metacell type-flow matrices of each
    # transition, which (up to a constant factor) is the average of the two.
    n_cells_per_prev_type_per_type_per_transition = AbstractMatrix{<:Real}[]
    for round_index in 2:n_rounds
        base_daf = base_daf_per_round[round_index]
        @assert has_matrix(base_daf, "type", "type", "n_cells_by_block") "missing the n_cells_by_block matrix"
        @assert has_matrix(base_daf, "type", "type", "n_cells_by_metacell") "missing the n_cells_by_metacell matrix"
        push!(
            n_cells_per_prev_type_per_type_per_transition,
            get_matrix(base_daf, "type", "type", "n_cells_by_block").array,
        )
        push!(
            n_cells_per_prev_type_per_type_per_transition,
            get_matrix(base_daf, "type", "type", "n_cells_by_metacell").array,
        )
    end
    push!(
        n_cells_per_prev_type_per_type_per_transition,
        get_matrix(final_daf, "type", "type", "n_cells_by_block").array,
    )
    push!(
        n_cells_per_prev_type_per_type_per_transition,
        get_matrix(final_daf, "type", "type", "n_cells_by_metacell").array,
    )

    global_flow_order_per_type =
        compute_global_flow_order_per_type(n_cells_per_prev_type_per_type_per_transition; restarts, rng)
    set_vector!(final_daf, "type", "global_flow_order", UInt32.(global_flow_order_per_type); overwrite)

    return nothing
end

end  # module
