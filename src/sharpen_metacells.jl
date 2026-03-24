"""
Compute "better" metacells based on the local linear model approximating the manifold.
"""
module SharpenMetacells

export sharpen_metacells!

using Base.Threads
using Clustering
using DataAxesFormats
using Random
using SparseArrays
using TanayLabUtilities

using ..Contracts

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
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_block
import Metacells.Contracts.vector_of_n_metacells_per_block
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell

struct KmeansSizesBuffers{T <: AbstractFloat}
    kmeans_buffers::Tuple{KMeansBuffers{T}, KMeansBuffers{T}}
    max_best_counts::Vector{Int}
    max_centers::Matrix{T}
    max_sorted_sizes::Vector{Int}
    max_split_cluster_indices::Vector{Int}
    max_split_point_indices::Vector{Int}
    max_split_points::Matrix{T}
    max_split_assignments::Vector{Int}
end

function KmeansSizesBuffers{T}(; n_dims::Integer, max_k::Integer, n_points::Integer) where {T <: AbstractFloat}
    return KmeansSizesBuffers{T}(
        (KMeansBuffers{T}(; n_dims, max_k, n_points), KMeansBuffers{T}(; n_dims, max_k, n_points)),
        Vector{Int}(undef, max_k),            # max_best_counts
        Matrix{T}(undef, n_dims, max_k),      # max_centers
        Vector{Int}(undef, max_k),            # max_sorted_sizes
        Vector{Int}(undef, max_k),            # max_split_cluster_indices
        Vector{Int}(undef, n_points),         # max_split_point_indices
        Matrix{T}(undef, n_dims, n_points),   # max_split_points
        Vector{Int}(undef, n_points),         # max_split_assignments
    )
end

struct KmeansSizesResult
    assignments::AbstractVector{Int}
    counts::AbstractVector{Int}
end

Clustering.assignments(r::KmeansSizesResult) = r.assignments
Clustering.counts(r::KmeansSizesResult) = r.counts
Clustering.nclusters(r::KmeansSizesResult) = length(r.counts)

@kwdef struct LocalClusters
    block_cell_indices::AbstractVector{<:Integer}
    cluster_index_per_block_cell::AbstractVector{<:Integer}
    is_too_small_per_cluster::BitVector
end

"""
    function sharpen_metacells!(;
        sharp_daf::DafWriter,
        base_daf::DafReader,
        prefix::AbstractString = $(DEFAULT.prefix),
        min_cells_in_metacell::Integer = $(DEFAULT.min_cells_in_metacell),
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
    enforce the sizes of the clusters - not more than twice that mean and no less than `min_cells_in_metacell`. In
    edge cases we dissolve too-small clusters, so this can create new outlier cells.
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
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_cells_per_block(RequiredInput),
        vector_of_n_metacells_per_block(RequiredInput),
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
    min_migration_likelihood::AbstractFloat = 1.5,
    kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert min_cells_in_metacell >= 0
    @assert min_migration_likelihood > 0
    @assert kmeans_rounds > 0

    n_cells = axis_length(base_daf, "cell")
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

    n_cells_per_block = copy_array(get_vector(base_daf, "block", "n_cells").array)
    n_metacells_per_block = get_vector(base_daf, "block", "n_metacells").array
    mean_metacell_cells_per_block = n_cells_per_block ./ max.(n_metacells_per_block, 1)

    block_index_per_cell = base_daf["@ cell : metacell ?? 0 : block : index"].array
    n_migrated = n_cells

    n_blocks = length(name_per_block)
    n_modules_per_block = get_vector(base_daf, "block", "n_modules").array
    max_n_block_modules = maximum(n_modules_per_block)
    n_neighborhood_cells_per_block = get_vector(base_daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

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
        ) for _ in 1:nthreads()
    ]

    preferred_block_index_per_cell_per_block = compute_preferred_block_index_per_cell_per_block(;
        base_daf,
        kmeans_rounds,
        name_per_block,
        n_cells_per_block,
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

    block_index_per_cell, n_migrated =
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
        @warn "TODOX post_max_n_kmeans_points: $(post_max_n_kmeans_points) ? max_n_kmeans_points: $(max_n_kmeans_points)"
        @warn "TODOX post_max_n_kmeans_clusters: $(post_max_n_kmeans_clusters) ? max_n_kmeans_clusters: $(max_n_kmeans_clusters)"
        max_n_kmeans_points = max(max_n_kmeans_points, post_max_n_kmeans_points)
        max_n_kmeans_clusters = max(max_n_kmeans_clusters, post_max_n_kmeans_clusters)
        kmeans_sizes_max_buffers_per_thread = [
            KmeansSizesBuffers{Float32}(;
                n_dims = max_n_block_modules,
                max_k = max_n_kmeans_clusters,
                n_points = max_n_kmeans_points,
            ) for _ in 1:nthreads()
        ]
    end

    local_clusters_per_block = compute_local_clusters(;
        base_daf,
        UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        min_cells_in_metacell,
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

    n_neighborhood_cells_per_block = get_vector(base_daf, "block", "n_neighborhood_cells").array
    max_n_neighborhood_cells = maximum(n_neighborhood_cells_per_block)

    preferred_block_index_per_cell_per_block = Vector{Maybe{SparseVector{<:Integer}}}(undef, n_blocks)
    preferred_block_index_per_cell_per_block .= nothing

    n_genes = size(module_index_per_gene_per_block, 1)

    z_score_per_max_module_per_max_neighborhood_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_neighborhood_cells) for _ in 1:nthreads()]
    is_in_neighborhood_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:nthreads()]
    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:nthreads()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:nthreads()]
    max_n_neighborhood_clusters = maximum(
        max(Int(round(n_neighborhood_cells_per_block[b] / mean_metacell_cells_per_block[b])), 1) for b in 1:n_blocks
    )
    preferred_block_index_per_cluster_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_clusters) for _ in 1:nthreads()]
    block_count_scratch_per_thread = [zeros(Int, n_blocks) for _ in 1:nthreads()]
    preferred_block_index_per_max_neighborhood_cell_per_thread =
        [Vector{Int}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]

    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        policy = :static,
        name = "compute_preferred_block_index_per_cell_per_block",
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "preferred_block_index_per_cell_per_block"),
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
)::Tuple{Vector{<:Integer}, Integer}
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

    return (new_block_index_per_cell, n_migrated[])
end

function compute_local_clusters(;
    base_daf::DafReader,
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    min_cells_in_metacell::Integer,
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

    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:nthreads()]
    is_in_block_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:nthreads()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:nthreads()]
    z_score_per_max_module_per_max_block_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_block_cells) for _ in 1:nthreads()]
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

        if n_block_clusters == 1
            fill!(cluster_index_per_block_cell, 1)
            sizes_buffers.max_best_counts[1] = n_block_cells
            cluster_sizes = @view(sizes_buffers.max_best_counts[1:1])
        else
            kmeans_result = flame_timed("kmeans_with_sizes") do
                return kmeans_with_sizes(
                    z_score_per_found_module_per_block_cell,
                    n_block_clusters;
                    max_k = 2 * n_block_clusters,
                    min_cluster_size = min_cells_in_metacell,
                    max_cluster_size = mean_metacell_cells_per_block[block_index] * 2,
                    kmeans_rounds,
                    sizes_buffers,
                    best_assignments = cluster_index_per_block_cell,
                    rng,
                )
            end

            cluster_sizes = counts(kmeans_result)
        end
        local_clusters_per_block[block_index] = LocalClusters(;
            block_cell_indices,
            cluster_index_per_block_cell,
            is_too_small_per_cluster = cluster_sizes .< min_cells_in_metacell,
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

function kmeans_with_sizes(
    values_of_points::AbstractMatrix{<:AbstractFloat},
    initial_k::Integer;
    max_k::Integer,
    min_cluster_size::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    sizes_buffers::KmeansSizesBuffers,
    best_assignments::AbstractVector{Int},
    rng::AbstractRNG,
)::KmeansSizesResult
    kmeans_result = nothing
    n_dims, n_points = size(values_of_points)
    k = min(initial_k, max_k)
    centers = nothing
    cluster_sizes = nothing

    do_next_phase = true
    retries = 0
    has_kmeans_result = false

    best_kmeans_result = nothing
    best_k = nothing
    best_n_too_small = nothing

    buffers = sizes_buffers.kmeans_buffers
    best_counts = sizes_buffers.max_best_counts
    centers_buffer = sizes_buffers.max_centers
    max_sorted_sizes = sizes_buffers.max_sorted_sizes
    max_split_cluster_indices = sizes_buffers.max_split_cluster_indices
    max_split_point_indices = sizes_buffers.max_split_point_indices
    max_split_points = sizes_buffers.max_split_points
    max_split_assignments = sizes_buffers.max_split_assignments

    while do_next_phase
        do_next_phase = false

        while true
            if has_kmeans_result
                has_kmeans_result = false
                cluster_sizes = counts(kmeans_result)  # NOJET
                @assert length(cluster_sizes) == k
            else
                kmeans_result = flame_timed("kmeans_in_rounds") do
                    return kmeans_in_rounds(values_of_points, k; centers, buffers, rounds = kmeans_rounds, rng)
                end
                cluster_sizes = counts(kmeans_result)
                @assert length(cluster_sizes) == k
            end

            if k < initial_k
                need_split = min(initial_k - k, k)
                sorted_sizes = @view(max_sorted_sizes[1:k])
                copyto!(sorted_sizes, cluster_sizes)
                partialsort!(sorted_sizes, need_split; rev = true)  # NOJET
                effective_max_cluster_size = min(max_cluster_size, sorted_sizes[need_split] - 1)
            else
                effective_max_cluster_size = max_cluster_size
            end
            effective_max_cluster_size = max(effective_max_cluster_size, 2 * min_cluster_size)

            largest_cluster_size = maximum(cluster_sizes)
            if largest_cluster_size <= effective_max_cluster_size
                break
            end

            n_large_clusters = 0
            @inbounds for i in 1:k
                if cluster_sizes[i] > effective_max_cluster_size
                    n_large_clusters += 1
                end
            end
            n_split_clusters = max(1, div(n_large_clusters + 1, 2))

            # Find the n_split_clusters largest clusters by insertion into split_cluster_indices
            split_cluster_indices = @view(max_split_cluster_indices[1:n_split_clusters])
            fill!(split_cluster_indices, 0)
            @inbounds for i in 1:k
                s = cluster_sizes[i]
                j = n_split_clusters
                while j >= 1 && (split_cluster_indices[j] == 0 || s > cluster_sizes[split_cluster_indices[j]])
                    if j < n_split_clusters
                        split_cluster_indices[j + 1] = split_cluster_indices[j]
                    end
                    split_cluster_indices[j] = i
                    j -= 1
                end
            end

            if k == max_k
                break
            end

            # Copy assignments and centers before reusing buffers for split
            copyto!(@view(max_split_assignments[1:n_points]), assignments(kmeans_result))  # NOJET
            clusters_of_points = @view(max_split_assignments[1:n_points])
            copyto!(@view(centers_buffer[1:n_dims, 1:k]), kmeans_result.centers)
            for si in 1:n_split_clusters
                split_cluster_index = split_cluster_indices[si]
                n_split_points = 0
                @inbounds for j in 1:n_points
                    if clusters_of_points[j] == split_cluster_index
                        n_split_points += 1
                        max_split_point_indices[n_split_points] = j
                    end
                end
                split_point_indices = @view(max_split_point_indices[1:n_split_points])
                @inbounds for j in 1:n_split_points
                    for i in 1:n_dims
                        max_split_points[i, j] = values_of_points[i, split_point_indices[j]]
                    end
                end
                split_points = @view(max_split_points[1:n_dims, 1:n_split_points])
                split_buf = KMeansBuffers(buffers[1]; n_dims, k = 2, n_points = n_split_points)
                split_result = flame_timed("kmeans") do
                    return kmeans_in_buffers(split_points, 2; buffers = split_buf, rng)
                end
                centers_buffer[1:n_dims, split_cluster_index] .= split_result.centers[:, 1]
                k += 1
                centers_buffer[1:n_dims, k] .= split_result.centers[:, 2]
                if k == max_k
                    break
                end
            end
            centers = @view(centers_buffer[1:n_dims, 1:k])
        end

        n_too_small = 0
        while k > 1
            @assert length(cluster_sizes) == k
            smallest_cluster_index = argmin(cluster_sizes)

            n_too_small = count(size -> size < min_cluster_size, cluster_sizes)
            n_too_large = count(size -> size > max_cluster_size, cluster_sizes)
            if n_too_large == 0 && (
                best_kmeans_result === nothing ||
                (best_k >= initial_k && n_too_small < best_n_too_small) ||
                (
                    best_k < initial_k && (
                        (k > best_k && n_too_small <= best_n_too_small) ||
                        (k == best_k && n_too_small < best_n_too_small)
                    )
                )
            )
                best_kmeans_result = kmeans_result
                best_k = k
                best_n_too_small = n_too_small
                copyto!(best_assignments, assignments(kmeans_result))  # NOJET
                copyto!(@view(best_counts[1:k]), counts(kmeans_result))  # NOJET
            end

            if n_too_small == 0
                break
            end

            j = 0
            for i in 1:k
                if i != smallest_cluster_index
                    j += 1
                    centers_buffer[1:n_dims, j] .= kmeans_result.centers[:, i]
                end
            end
            k -= 1
            centers = @view(centers_buffer[1:n_dims, 1:k])

            merged_kmeans_result = flame_timed("kmeans_in_rounds") do
                return kmeans_in_rounds(values_of_points, k; centers, buffers, rounds = kmeans_rounds, rng)
            end
            cluster_sizes = counts(merged_kmeans_result)
            @assert length(cluster_sizes) == k
            largest_cluster_size = maximum(cluster_sizes)

            if n_too_large == 0 && (
                best_kmeans_result === nothing ||
                (best_k >= initial_k && n_too_small < best_n_too_small) ||
                (
                    best_k < initial_k && (
                        (k > best_k && n_too_small <= best_n_too_small) ||
                        (k == best_k && n_too_small < best_n_too_small)
                    )
                )
            )
                best_kmeans_result = merged_kmeans_result
                best_k = k
                best_n_too_small = n_too_small
                copyto!(best_assignments, assignments(merged_kmeans_result))
                copyto!(@view(best_counts[1:k]), counts(merged_kmeans_result))
            end

            if k < initial_k || largest_cluster_size > max_cluster_size
                if retries < 10
                    do_next_phase = true
                    has_kmeans_result = true
                    kmeans_result = merged_kmeans_result
                    retries += 1
                    break
                end
            end

            kmeans_result = merged_kmeans_result
        end

        @assert length(counts(kmeans_result)) == k  # NOJET
    end

    @assert best_kmeans_result !== nothing
    return KmeansSizesResult(best_assignments, @view(best_counts[1:best_k]))
end

end  # module
