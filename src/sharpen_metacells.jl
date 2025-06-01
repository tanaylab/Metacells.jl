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
using StatsBase
using TanayLabUtilities

using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_gene_module_index_matrix
import Metacells.Contracts.block_module_is_used_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.module_axis
"""
    function sharpen_metacells!(
        original::DafReader,
        sharpened::DafWriter;
        scale::Bool = $(DEFAULT.scale),
        migrate::Bool = $(DEFAULT.migrate),
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        min_cells_in_metacell::Integer = $(DEFAULT.min_cells_in_metacell),
        min_migrated_cells_fraction::AbstractFloat = $(DEFAULT.min_migrated_cells_fraction),
        kmeans_rounds::Integer = $(DEFAULT.kmeans_rounds),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given an `original` repository with a local linear model approximating the manifold, compute new metacells into the
`sharpened` repository, which hopefully more faithfully capture this manifold.

TODOX

# Original Metacells

$(CONTRACT1)

# Sharpened Metacells

$(CONTRACT2)
"""
@logged @computation Contract(;
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_is_used_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(OptionalInput),
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), metacell_axis(GuaranteedOutput)],
    data = [cell_metacell_vector(GuaranteedOutput), metacell_block_vector(GuaranteedOutput)],
) function sharpen_metacells!(
    original::DafReader,
    sharpened::DafWriter;
    scale::Bool = false,
    migrate::Bool = true,
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    min_cells_in_metacell::Integer = 12,
    min_migrated_cells_fraction::AbstractFloat = 0.01,
    kmeans_rounds::Integer = 10,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert !scale  # TODOX
    n_cells = axis_length(original, "cell")
    n_blocks = axis_length(original, "block")
    name_per_block = axis_vector(original, "block")

    cells_of_new_metacells_per_block = Vector{Vector{Vector{UInt32}}}(undef, n_blocks)
    UMIs_per_covered_gene_per_cell = original["/ gene & is_covered / cell : UMIs"].array

    if migrate
        block_indices_per_neighborhood = [
            parent(original["/ block & is_in_neighborhood ;= $(block_name) : index"].array) for
            block_name in name_per_block
        ]

        mean_metacell_cells_per_block = [
            original["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : index %> Count"] /
            original["/ metacell & block => is_in_neighborhood ;= $(block_name) : index %> Count"] for
            block_name in name_per_block
        ]

        block_index_per_cell = original["/ cell : metacell ?? 0 => block => index"].array
        block_index_per_cell_per_round = AbstractVector{<:Integer}[block_index_per_cell]
        n_migrated = n_cells
        while n_migrated >= n_cells * min_migrated_cells_fraction
            preferred_block_index_per_cell_per_block = compute_preferred_block_index_per_cell_per_block(
                original;
                kmeans_rounds,
                name_per_block,
                block_index_per_cell,
                UMIs_per_covered_gene_per_cell,
                mean_metacell_cells_per_block,
                block_indices_per_neighborhood,
                min_downsamples,
                min_downsamples_quantile,
                max_downsamples_quantile,
                rng,
            )

            block_index_per_cell, n_migrated = compute_preferred_block_index_of_cells(;
                block_index_per_cell,
                preferred_block_index_per_cell_per_block,
                block_index_per_cell_per_round,
            )

            push!(block_index_per_cell_per_round, block_index_per_cell)
        end

        cells_of_new_metacells, block_name_per_new_metacell = compute_new_metacells(
            original;
            UMIs_per_covered_gene_per_cell,
            min_cells_in_metacell,
            kmeans_rounds,
            name_per_block,
            block_index_per_cell,
            mean_metacell_cells_per_block,
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
            rng,
        )

    else
        parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng
            block_name = name_per_block[block_index]
            cells_of_new_metacells_of_block = sharpen_block_without_migration(
                original;
                block_name,
                UMIs_per_covered_gene_per_cell,
                min_downsamples,
                min_downsamples_quantile,
                max_downsamples_quantile,
                min_cells_in_metacell,
                kmeans_rounds,
                rng,
            )
            cells_of_new_metacells_per_block[block_index] = cells_of_new_metacells_of_block
            return nothing
        end

        cells_of_new_metacells = vcat(cells_of_new_metacells_per_block...)
        block_name_per_new_metacell = vcat(
            [
                repeat([name_per_block[block_index]], length(cells_of_new_metacells_per_block[block_index])) for
                block_index in 1:n_blocks
            ]...,
        )
    end

    name_per_new_metacell = group_names(original, "cell", cells_of_new_metacells; prefix = "M")  # NOJET

    if overwrite
        delete_axis!(sharpened, "metacell"; must_exist = false)
    end
    add_axis!(sharpened, "metacell", name_per_new_metacell)

    new_metacell_name_per_cell = fill("", n_cells)
    for (new_metacell_name, cells_of_new_metacell) in zip(name_per_new_metacell, cells_of_new_metacells)
        new_metacell_name_per_cell[cells_of_new_metacell] .= new_metacell_name
    end
    set_vector!(sharpened, "cell", "metacell", new_metacell_name_per_cell)
    set_vector!(sharpened, "metacell", "block", block_name_per_new_metacell)

    return nothing
end

function compute_preferred_block_index_per_cell_per_block(
    daf::DafReader;
    kmeans_rounds::Integer,
    name_per_block::AbstractVector{<:AbstractString},
    block_index_per_cell::AbstractVector{<:Integer},
    UMIs_per_covered_gene_per_cell::AbstractMatrix{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    block_indices_per_neighborhood::AbstractVector{<:AbstractVector{<:Integer}},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
    rng::AbstractRNG,
)::Vector{Maybe{SparseVector{<:Integer}}}
    n_cells = length(block_index_per_cell)
    n_blocks = length(name_per_block)

    cell_indices_per_block = [findall(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    preferred_block_index_per_cell_per_block = Vector{Maybe{SparseVector{<:Integer}}}(undef, n_blocks)
    preferred_block_index_per_cell_per_block .= nothing

    # progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng
        block_name = name_per_block[block_index]

        neighborhood_cell_indices = vcat(cell_indices_per_block[block_indices_per_neighborhood[block_index]]...)
        n_neighborhood_cells = length(neighborhood_cell_indices)
        if n_neighborhood_cells == 0
            return nothing
        end

        fraction_per_used_module_per_neighborhood_cell = compute_fraction_per_used_module_per_region_cell(
            daf;
            block_name,
            UMIs_per_covered_gene_per_cell,
            region_cell_indices = neighborhood_cell_indices,
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        n_neighborhood_clusters = max(Int(round(n_neighborhood_cells / mean_metacell_cells_per_block[block_index])), 1)

        kmeans_result = kmeans_in_rounds(
            fraction_per_used_module_per_neighborhood_cell,
            n_neighborhood_clusters;
            kmeans_rounds,
            rng,
        )
        cluster_index_per_neighborhood_cell = assignments(kmeans_result)
        n_clusters = length(counts(kmeans_result))

        block_index_per_neighborhood_cell = vcat(
            [
                fill(block_index, length(cell_indices_per_block[block_index])) for
                block_index in block_indices_per_neighborhood[block_index]
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

        return preferred_block_index_per_cell_per_block[block_index] = SparseVector(
            n_cells,
            neighborhood_cell_indices[reorder],
            preferred_block_index_per_neighborhood_cell[reorder],
        )

        # counter = atomic_add!(progress_counter, 1)
        # n_stable_cells = sum(preferred_block_index_per_neighborhood_cell .== block_index)
        # @debug (
        #     "- Neighborhood: $(block_name) ($(percent(counter + 1, n_blocks)))" *
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

    new_block_index_per_cell = zeros(UInt32, n_cells)
    @threads :greedy for cell_index in 1:n_cells
        original_block_index_of_cell = block_index_per_cell[cell_index]
        if original_block_index_of_cell == 0
            continue
        end

        preferred_block_index_per_block = preferred_block_index_per_cell_per_block[original_block_index_of_cell]
        if preferred_block_index_per_block === nothing
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_quiescent, 1)
            continue
        end

        preferred_block_index_of_cell = preferred_block_index_per_block[cell_index]
        @assert preferred_block_index_of_cell > 0

        preferred_block_index_per_other_block = preferred_block_index_per_cell_per_block[preferred_block_index_of_cell]
        if preferred_block_index_per_other_block === nothing
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
            continue
        end

        back_preferred_block_index_of_cell = preferred_block_index_per_other_block[cell_index]
        @assert back_preferred_block_index_of_cell > 0

        if preferred_block_index_of_cell == original_block_index_of_cell ||
           back_preferred_block_index_of_cell != preferred_block_index_of_cell
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
        else
            new_block_index_per_cell[cell_index] = preferred_block_index_of_cell
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

    return (new_block_index_per_cell, n_migrated[])
end

function compute_new_metacells(
    daf::DafReader;
    UMIs_per_covered_gene_per_cell::AbstractMatrix{<:Integer},
    min_cells_in_metacell::Integer,
    kmeans_rounds::Integer,
    name_per_block::AbstractVector{<:AbstractString},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
    rng::AbstractRNG,
)::Tuple{AbstractVector{<:AbstractVector{<:Integer}}, AbstractVector{<:AbstractString}}
    local_clusters_per_block = compute_local_clusters(
        daf;
        UMIs_per_covered_gene_per_cell,
        min_cells_in_metacell,
        kmeans_rounds,
        name_per_block,
        block_index_per_cell,
        mean_metacell_cells_per_block,
        min_downsamples,
        min_downsamples_quantile,
        max_downsamples_quantile,
        rng,
    )

    n_cells = size(UMIs_per_covered_gene_per_cell, 2)
    n_old_metacells = axis_length(daf, "metacell")

    return combine_local_clusters(; local_clusters_per_block, name_per_block, n_old_metacells, n_cells)
end

@kwdef struct LocalClusters
    block_cell_indices::AbstractVector{<:Integer}
    cluster_index_per_block_cell::AbstractVector{<:Integer}
    is_too_small_per_cluster::Union{AbstractVector{Bool}, BitVector}
end

function compute_local_clusters(
    daf::DafReader;
    UMIs_per_covered_gene_per_cell::AbstractMatrix{<:Integer},
    min_cells_in_metacell::Integer,
    kmeans_rounds::Integer,
    name_per_block::AbstractVector{<:AbstractString},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
    rng::AbstractRNG,
)::AbstractVector{Maybe{LocalClusters}}
    n_blocks = length(name_per_block)

    local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)

    # progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng, policy = :serial) do block_index, rng
        block_name = name_per_block[block_index]

        block_cell_indices = findall(block_index_per_cell .== block_index)
        n_block_cells = length(block_cell_indices)
        if n_block_cells == 0
            local_clusters_per_block[block_index] = nothing
            # counter = atomic_add!(progress_counter, 1)
            # @debug "- Block: $(block_name) ($(percent(counter + 1, n_blocks))) Empty"
            return nothing
        end

        fraction_per_used_module_per_block_cell = compute_fraction_per_used_module_per_region_cell(
            daf;
            block_name,
            UMIs_per_covered_gene_per_cell,
            region_cell_indices = block_cell_indices,
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        n_block_clusters = max(Int(round(n_block_cells / mean_metacell_cells_per_block[block_index])), 1)

        kmeans_result = kmeans_with_sizes(
            fraction_per_used_module_per_block_cell,
            n_block_clusters;
            min_cluster_size = min_cells_in_metacell,
            max_cluster_size = mean_metacell_cells_per_block[block_index] * 2,
            kmeans_rounds,
            rng,
        )

        cluster_sizes = counts(kmeans_result)
        return local_clusters_per_block[block_index] = LocalClusters(;
            block_cell_indices,
            cluster_index_per_block_cell = assignments(kmeans_result),
            is_too_small_per_cluster = cluster_sizes .< min_cells_in_metacell,
        )

        # n_clusters = length(cluster_sizes)
        # n_large_clusters = sum(cluster_sizes .> 2 * mean_metacell_cells_per_block[block_index])
        # n_small_clusters = sum(cluster_sizes .< min_cells_in_metacell)
        # counter = atomic_add!(progress_counter, 1)
        # @debug (
        #     "- Block: $(block_name) ($(percent(counter + 1, n_blocks)))" *
        #     " Clusters: $(n_clusters)" *
        #     " Min: $(minimum(cluster_sizes))" *
        #     " Mean: $(mean(cluster_sizes))" *
        #     " (MCs: $(mean_metacell_cells_per_block[block_index]))" *
        #     " Max: $(maximum(cluster_sizes))" *
        #     " Too small: $(n_small_clusters)" *
        #     " Too large: $(n_large_clusters)"
        # )
    end

    return local_clusters_per_block
end

function combine_local_clusters(;
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    name_per_block::AbstractVector{<:AbstractString},
    n_old_metacells::Integer,
    n_cells::Integer,
)::Tuple{AbstractVector{<:AbstractVector{<:Integer}}, <:AbstractVector{<:AbstractString}}
    n_total_new_metacells = 0
    n_new_outlier_cells = 0
    for local_clusters in local_clusters_per_block
        if local_clusters !== nothing
            n_total_new_metacells += length(local_clusters.is_too_small_per_cluster)
            for too_small_cluster_index in findall(local_clusters.is_too_small_per_cluster)
                n_new_outlier_cells += sum(local_clusters.cluster_index_per_block_cell .== too_small_cluster_index)
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
                    push!(cells_of_new_metacells, cell_indices_of_new_metacell)
                    push!(block_name_per_new_metacell, block_name)
                end
            end
        end
    end

    @debug (
        "Metacells Old: $(n_old_metacells)" *
        " New: $(length(block_name_per_new_metacell))" *
        " Outliers: $(n_new_outlier_cells)" *
        " ($(percent(n_new_outlier_cells, n_cells)))"
    )

    return (cells_of_new_metacells, block_name_per_new_metacell)
end

function sharpen_block_without_migration(
    daf::DafReader;
    block_name::AbstractString,
    UMIs_per_covered_gene_per_cell::AbstractMatrix{<:Integer},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
    min_cells_in_metacell::Integer,
    kmeans_rounds::Integer = 10,
    rng::AbstractRNG,
)::Vector{Vector{UInt32}}
    block_cell_indices = daf["/ cell & metacell ?? => block = $(block_name) : index"].array

    fraction_per_used_module_per_block_cell = compute_fraction_per_used_module_per_region_cell(
        daf;
        block_name,
        UMIs_per_covered_gene_per_cell,
        region_cell_indices = block_cell_indices,
        min_downsamples,
        min_downsamples_quantile,
        max_downsamples_quantile,
    )

    n_block_metacells = daf["/ metacell & block = $(block_name) %> Count"]
    mean_cells_in_metacell = length(block_cell_indices) / n_block_metacells

    kmeans_result = kmeans_with_sizes(
        fraction_per_used_module_per_block_cell,
        n_block_metacells;
        min_cluster_size = min_cells_in_metacell,
        max_cluster_size = 2 * mean_cells_in_metacell,
        kmeans_rounds,
        rng,
    )

    cluster_index_per_block_cell = assignments(kmeans_result)
    cell_positions_per_cluster = collect_group_members(cluster_index_per_block_cell)
    cell_indices_per_cluster = [block_cell_indices[cell_positions] for cell_positions in cell_positions_per_cluster]
    @debug "- Block: $(block_name) Old metacells: $(n_block_metacells) New metacells: $(length(cell_indices_per_cluster))"
    return cell_indices_per_cluster
end

function compute_fraction_per_used_module_per_region_cell(
    daf::DafReader;
    block_name::AbstractString,
    UMIs_per_covered_gene_per_cell::AbstractMatrix{<:Integer},
    region_cell_indices::AbstractVector{<:Integer},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
)::Matrix{Float32}
    @views UMIs_per_covered_gene_per_region_cell = UMIs_per_covered_gene_per_cell[:, region_cell_indices]
    n_covered_genes, n_region_cells = size(UMIs_per_covered_gene_per_region_cell)

    covered_UMIs_per_region_cell = vec(sum(UMIs_per_covered_gene_per_region_cell; dims = 1))
    @assert_vector(covered_UMIs_per_region_cell, n_region_cells)

    downsamples_covered_UMIs =
        downsamples(covered_UMIs_per_region_cell; min_downsamples, min_downsamples_quantile, max_downsamples_quantile)

    downsampled_UMIs_per_covered_gene_per_region_cell =
        downsample(UMIs_per_covered_gene_per_region_cell, downsamples_covered_UMIs; dims = 1)

    total_downsampled_UMIs_per_region_cell = vec(sum(downsampled_UMIs_per_covered_gene_per_region_cell; dims = 1))
    @assert_vector(total_downsampled_UMIs_per_region_cell, n_region_cells)
    @assert all(total_downsampled_UMIs_per_region_cell .<= downsamples_covered_UMIs)

    fraction_per_covered_gene_per_region_cell =
        downsampled_UMIs_per_covered_gene_per_region_cell ./ transpose(total_downsampled_UMIs_per_region_cell)
    @assert_matrix(fraction_per_covered_gene_per_region_cell, n_covered_genes, n_region_cells)

    module_index_per_covered_gene = daf["/ gene & is_covered / block = $(block_name) : module_index"].array
    used_module_indices = daf["/ module & is_used ; block = $(block_name) : index"].array
    n_used_modules = length(used_module_indices)

    fraction_per_region_cell_per_used_module = Matrix{Float32}(undef, n_region_cells, n_used_modules)
    for (used_module_position, used_module_index) in enumerate(used_module_indices)
        covered_module_genes_positions = findall(module_index_per_covered_gene .== used_module_index)
        @assert length(covered_module_genes_positions) > 0

        @views fraction_per_covered_module_gene_per_region_cell =
            fraction_per_covered_gene_per_region_cell[covered_module_genes_positions, :]

        fraction_per_region_cell_per_used_module[:, used_module_position] =
            sum(fraction_per_covered_module_gene_per_region_cell; dims = 1)
    end

    return transposer(fraction_per_region_cell_per_used_module)
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
