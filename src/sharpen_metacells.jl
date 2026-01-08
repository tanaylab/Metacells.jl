"""
Compute "better" metacells based on the local linear model approximating the manifold.
"""
module SharpenMetacells

export sharpen_metacells!
export compute_metacells_radius!

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
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_is_strong_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.module_axis
"""
function sharpen_metacells!(
original::DafReader,
sharpened::DafWriter;
migrate::Bool = $(DEFAULT.migrate),
min_cells_in_metacell::Integer = $(DEFAULT.min_cells_in_metacell),
max_migrated_cells_fraction::AbstractFloat = $(DEFAULT.max_migrated_cells_fraction),
kmeans_rounds::Integer = $(DEFAULT.kmeans_rounds),
rng::AbstractRNG = default_rng(),
overwrite::Bool = $(DEFAULT.overwrite),
)::Nothing

Given an `original` metacells repository with a blocks structure and local gene modules that capture the essense of the
cell state manifold, compute a `sharpened` metacells repository, which hopefully more faithfully captures this manifold.

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
        cell_total_UMIs_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_sharpening_rounds_vector(OptionalInput),
        cell_gene_UMIs_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_is_strong_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(OptionalInput),
        block_module_neighborhood_mean_linear_fraction_matrix(OptionalInput),
        block_module_neighborhood_std_linear_fraction_matrix(OptionalInput),
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), metacell_axis(GuaranteedOutput)],
    data = [
        cell_metacell_vector(GuaranteedOutput),
        metacell_block_vector(GuaranteedOutput),
        metacell_sharpening_rounds_vector(GuaranteedOutput),
    ],
) function sharpen_metacells!(
    original::DafReader,
    sharpened::DafWriter;
    sharpening_rounds::Maybe{Integer} = nothing,
    migrate::Bool = false,
    min_cells_in_metacell::Integer = 12,
    min_migration_likelihood::AbstractFloat = 1.5,
    max_migrated_cells_fraction::AbstractFloat = 1.0,
    kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
    normalize_fractions::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert min_cells_in_metacell >= 0
    @assert min_migration_likelihood > 0
    @assert 0 <= max_migrated_cells_fraction <= 1
    @assert kmeans_rounds > 0

    n_cells = axis_length(original, "cell")
    n_blocks = axis_length(original, "block")
    n_modules = axis_length(original, "module")
    name_per_block = axis_vector(original, "block")

    cells_of_new_metacells_per_block = Vector{Vector{Vector{UInt32}}}(undef, n_blocks)
    UMIs_per_gene_per_cell = get_matrix(original, "gene", "cell", "UMIs").array
    total_UMIs_per_cell = get_vector(original, "cell", "total_UMIs").array

    mean_fraction_per_module_per_block = nothing
    std_fraction_per_module_per_block = nothing

    if migrate
        if normalize_fractions
            mean_fraction_per_module_per_block =
                get_matrix(original, "module", "block", "neighborhood_mean_linear_fraction").array
            std_fraction_per_module_per_block =
                get_matrix(original, "module", "block", "neighborhood_std_linear_fraction").array
        end

        block_indices_per_neighborhood = [
            parent(original["/ block & is_in_neighborhood ;= $(block_name) : index"].array) for
            block_name in name_per_block
        ]

        mean_metacell_cells_per_block = Vector{Float32}(undef, n_blocks)
        progress_counter = Atomic{Int}(0)
        @threads :greedy for block_index in 1:n_blocks
            block_name = name_per_block[block_index]
            mean_metacell_cells_per_block[block_index] =
                original["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : index %> Count"] /
                original["/ metacell & block => is_in_neighborhood ;= $(block_name) : index %> Count"]
            counter = atomic_add!(progress_counter, 1)
            print("\r$(counter) ($(percent(counter, n_blocks)))   ")
        end

        block_index_per_cell = original["/ cell : metacell ?? 0 => block => index"].array
        block_index_per_cell_per_round = AbstractVector{<:Integer}[block_index_per_cell]
        n_migrated = n_cells
        while n_migrated >= n_cells * max_migrated_cells_fraction
            preferred_block_index_per_cell_per_block = compute_preferred_block_index_per_cell_per_block(
                original;
                kmeans_rounds,
                normalize_fractions,
                name_per_block,
                block_index_per_cell,
                UMIs_per_gene_per_cell,
                total_UMIs_per_cell,
                mean_metacell_cells_per_block,
                block_indices_per_neighborhood,
                mean_fraction_per_module_per_block,
                std_fraction_per_module_per_block,
                min_migration_likelihood,
                rng,
            )

            block_index_per_cell, n_migrated = compute_preferred_block_index_of_cells(;
                block_index_per_cell,
                preferred_block_index_per_cell_per_block,
                block_index_per_cell_per_round,
            )

            push!(block_index_per_cell_per_round, block_index_per_cell)
        end

        local_clusters_per_block = compute_local_clusters(
            original;
            UMIs_per_gene_per_cell,
            total_UMIs_per_cell,
            min_cells_in_metacell,
            kmeans_rounds,
            normalize_fractions,
            name_per_block,
            block_index_per_cell,
            mean_metacell_cells_per_block,
            mean_fraction_per_module_per_block,
            std_fraction_per_module_per_block,
            rng,
        )

    else
        local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)
        parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng
            block_name = name_per_block[block_index]
            local_clusters_per_block[block_index] = sharpen_block_without_migration(
                original;
                block_name,
                UMIs_per_gene_per_cell,
                total_UMIs_per_cell,
                min_cells_in_metacell,
                kmeans_rounds,
                normalize_fractions,
                rng,
            )
            return nothing
        end
    end

    n_cells = axis_length(original, "cell")
    n_old_metacells = axis_length(original, "metacell")

    cells_of_new_metacells, block_name_per_new_metacell = combine_local_clusters(;
        min_cells_in_metacell,
        local_clusters_per_block,
        name_per_block,
        n_old_metacells,
        n_cells,
    )

    if sharpening_rounds === nothing
        sharpening_rounds = original["/ metacell : sharpening_rounds || 0 %> Max"] + 1
    end
    prefix = string('M' + sharpening_rounds)
    name_per_new_metacell = group_names(original, "cell", cells_of_new_metacells; prefix)  # NOJET

    new_metacell_name_per_cell = fill("", n_cells)
    for (new_metacell_name, cells_of_new_metacell) in zip(name_per_new_metacell, cells_of_new_metacells)
        new_metacell_name_per_cell[cells_of_new_metacell] .= new_metacell_name
    end

    add_axis!(sharpened, "metacell", name_per_new_metacell; overwrite)
    set_vector!(sharpened, "cell", "metacell", new_metacell_name_per_cell; overwrite)
    set_vector!(sharpened, "metacell", "block", block_name_per_new_metacell; overwrite)
    set_vector!(sharpened, "metacell", "sharpening_rounds", UInt32(sharpening_rounds); overwrite)

    return nothing
end

function compute_preferred_block_index_per_cell_per_block(
    daf::DafReader;
    kmeans_rounds::Integer,
    normalize_fractions::Bool,
    name_per_block::AbstractVector{<:AbstractString},
    block_index_per_cell::AbstractVector{<:Integer},
    UMIs_per_gene_per_cell::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    block_indices_per_neighborhood::AbstractVector{<:AbstractVector{<:Integer}},
    min_migration_likelihood::AbstractFloat,
    mean_fraction_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    std_fraction_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}},
    rng::AbstractRNG,
)::Vector{Maybe{SparseVector{<:Integer}}}
    n_cells = length(block_index_per_cell)
    n_blocks = length(name_per_block)

    cell_indices_per_block = [findall(block_index_per_cell .== block_index) for block_index in 1:n_blocks]
    preferred_block_index_per_cell_per_block = Vector{Maybe{SparseVector{<:Integer}}}(undef, n_blocks)
    preferred_block_index_per_cell_per_block .= nothing

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng  # TODOX
        block_name = name_per_block[block_index]

        neighborhood_cell_indices = vcat(cell_indices_per_block[block_indices_per_neighborhood[block_index]]...)
        n_neighborhood_cells = length(neighborhood_cell_indices)
        if n_neighborhood_cells == 0
            return nothing
        end

        indices_of_used_modules = daf["/ module & is_strong ; block = $(block_name) : index"].array
        n_used_modules = length(indices_of_used_modules)

        fraction_per_used_module_per_neighborhood_cell = compute_fraction_per_used_module_per_region_cell(
            daf;
            block_name,
            UMIs_per_gene_per_cell,
            total_UMIs_per_cell,
            region_cell_indices = neighborhood_cell_indices,
            indices_of_used_modules,
        )

        n_neighborhood_clusters = max(Int(round(n_neighborhood_cells / mean_metacell_cells_per_block[block_index])), 1)

        if normalize_fractions
            @views mean_fraction_per_module = mean_fraction_per_module_per_block[:, block_index]
            @views std_fraction_per_module = std_fraction_per_module_per_block[:, block_index]
            mean_fraction_per_used_module = mean_fraction_per_module[indices_of_used_modules]
            std_fraction_per_used_module = std_fraction_per_module[indices_of_used_modules]
            fraction_per_used_module_per_neighborhood_cell =
                (fraction_per_used_module_per_neighborhood_cell .- mean_fraction_per_used_module) ./
                std_fraction_per_used_module
        end

        kmeans_result = kmeans_in_rounds(
            fraction_per_used_module_per_neighborhood_cell,
            n_neighborhood_clusters;
            rounds = kmeans_rounds,
            rng,
        )
        cluster_index_per_neighborhood_cell = assignments(kmeans_result)
        n_cells_per_cluster = counts(kmeans_result)
        n_clusters = length(counts(kmeans_result))

        block_index_per_neighborhood_cell = vcat(
            [
                fill(block_index, length(cell_indices_per_block[block_index])) for
                block_index in block_indices_per_neighborhood[block_index]
            ]...,
        )
        @assert_vector(block_index_per_neighborhood_cell, n_neighborhood_cells)

        preferred_block_index_per_neighborhood_cell = pick_preferred_block_indices(;
            block_index,
            min_migration_likelihood,
            block_indices_in_neighborhood = block_indices_per_neighborhood[block_index],
            cell_indices_per_block,
            block_index_per_neighborhood_cell,
            cluster_index_per_neighborhood_cell,
            n_cells_per_cluster,
        )

        reorder = sortperm(neighborhood_cell_indices)

        preferred_block_index_per_cell_per_block[block_index] = SparseVector(
            n_cells,
            neighborhood_cell_indices[reorder],
            preferred_block_index_per_neighborhood_cell[reorder],
        )

        counter = atomic_add!(progress_counter, 1)
        print("\r$(counter) ($(percent(counter, n_blocks)))   ")
        return nothing
    end

    return preferred_block_index_per_cell_per_block
end

function pick_preferred_block_indices(;
    block_index::Integer,
    min_migration_likelihood::AbstractFloat,
    block_indices_in_neighborhood::AbstractVector{<:Integer},
    cell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    block_index_per_neighborhood_cell::AbstractVector{<:Integer},
    cluster_index_per_neighborhood_cell::AbstractVector{<:Integer},
    n_cells_per_cluster::AbstractVector{<:Integer},
)::AbstractVector{<:Integer}
    n_neighborhood_cells = length(block_index_per_neighborhood_cell)
    n_block_cells = length(cell_indices_per_block[block_index])

    n_clusters = length(n_cells_per_cluster)
    preferred_block_index_per_cluster = zeros(Int, n_clusters)

    for cluster_index in 1:n_clusters
        is_in_cluster_per_neighborhood_cell = cluster_index_per_neighborhood_cell .== cluster_index
        block_index_per_cluster_cell = block_index_per_neighborhood_cell[is_in_cluster_per_neighborhood_cell]
        n_block_cells_in_cluster = sum(block_index_per_cluster_cell .== block_index)

        if n_block_cells_in_cluster > 0
            most_frequent_block_index = mode(block_index_per_neighborhood_cell[is_in_cluster_per_neighborhood_cell])
            if most_frequent_block_index == block_index
                preferred_block_index_per_cluster[cluster_index] = most_frequent_block_index
            else
                n_cluster_cells = n_cells_per_cluster[cluster_index]

                n_most_frequent_block_cells_in_cluster = sum(block_index_per_cluster_cell .== most_frequent_block_index)
                most_frequent_block_fraction_in_cluster = n_most_frequent_block_cells_in_cluster / n_cluster_cells

                block_fraction_in_cluster = n_block_cells_in_cluster / n_cluster_cells

                observed_block_fraction_in_both =
                    block_fraction_in_cluster / (block_fraction_in_cluster + most_frequent_block_fraction_in_cluster)

                n_most_frequent_block_cells = length(cell_indices_per_block[most_frequent_block_index])
                expected_block_fraction_in_both = n_block_cells / (n_block_cells + n_most_frequent_block_cells)

                if observed_block_fraction_in_both * min_migration_likelihood <= expected_block_fraction_in_both
                    preferred_block_index_per_cluster[cluster_index] = most_frequent_block_index
                end
            end
        end
    end

    preferred_block_index_per_neighborhood_cell = preferred_block_index_per_cluster[cluster_index_per_neighborhood_cell]

    return preferred_block_index_per_neighborhood_cell
end

function compute_preferred_block_index_of_cells(;
    block_index_per_cell::AbstractVector{<:Integer},
    preferred_block_index_per_cell_per_block::Vector{Maybe{SparseVector{<:Integer}}},
    block_index_per_cell_per_round::AbstractVector{<:AbstractVector{<:Integer}},
)::Tuple{Vector{<:Integer}, Integer}
    n_cells = length(block_index_per_cell)

    n_stationary = Atomic{Int}(0)
    n_restless = Atomic{Int}(0)
    n_restless = Atomic{Int}(0)
    n_orbited = Atomic{Int}(0)
    n_migrated = Atomic{Int}(0)

    new_block_index_per_cell = zeros(UInt32, n_cells)
    progress_counter = Atomic{Int}(0)
    @threads :greedy for cell_index in 1:n_cells
        counter = atomic_add!(progress_counter, 1)
        print("\r$(counter) ($(percent(counter, n_cells)))   ")
        original_block_index_of_cell = block_index_per_cell[cell_index]
        if original_block_index_of_cell == 0
            continue
        end

        preferred_block_index_per_block = preferred_block_index_per_cell_per_block[original_block_index_of_cell]
        if preferred_block_index_per_block === nothing
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            continue
        end

        preferred_block_index_of_cell = preferred_block_index_per_block[cell_index]
        if preferred_block_index_of_cell == 0
            preferred_block_index_per_other_block = nothing
        else
            preferred_block_index_per_other_block =
                preferred_block_index_per_cell_per_block[preferred_block_index_of_cell]
        end

        if preferred_block_index_per_other_block === nothing
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
            continue
        end

        back_preferred_block_index_of_cell = preferred_block_index_per_other_block[cell_index]
        if preferred_block_index_of_cell == original_block_index_of_cell
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_stationary, 1)
        elseif back_preferred_block_index_of_cell != preferred_block_index_of_cell
            new_block_index_per_cell[cell_index] = original_block_index_of_cell
            atomic_add!(n_restless, 1)
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
        " Restless: $(n_restless[]) ($(percent(n_restless[], n_cells)))" *
        " Orbited: $(n_orbited[]) ($(percent(n_orbited[], n_cells)))" *
        " Migrated: $(n_migrated[]) ($(percent(n_migrated[], n_cells)))"
    )

    return (new_block_index_per_cell, n_migrated[])
end

@kwdef struct LocalClusters
    block_cell_indices::AbstractVector{<:Integer}
    cluster_index_per_block_cell::AbstractVector{<:Integer}
    is_too_small_per_cluster::Union{AbstractVector{Bool}, BitVector}
end

function compute_local_clusters(
    daf::DafReader;
    UMIs_per_gene_per_cell::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    min_cells_in_metacell::Integer,
    kmeans_rounds::Integer,
    normalize_fractions::Bool,
    name_per_block::AbstractVector{<:AbstractString},
    block_index_per_cell::AbstractVector{<:Integer},
    mean_metacell_cells_per_block::AbstractVector{<:AbstractFloat},
    mean_fraction_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing,
    std_fraction_per_module_per_block::Maybe{AbstractMatrix{<:AbstractFloat}} = nothing,
    rng::AbstractRNG,
)::AbstractVector{Maybe{LocalClusters}}
    n_blocks = length(name_per_block)

    local_clusters_per_block = Vector{Maybe{LocalClusters}}(undef, n_blocks)

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng, policy = :serial) do block_index, rng
        counter = atomic_add!(progress_counter, 1)
        print("\r$(counter) ($(percent(counter, n_blocks)))   ")
        block_name = name_per_block[block_index]

        block_cell_indices = findall(block_index_per_cell .== block_index)
        n_block_cells = length(block_cell_indices)
        if n_block_cells == 0
            local_clusters_per_block[block_index] = nothing
            counter = atomic_add!(progress_counter, 1)
            @debug "- Block: $(block_name) ($(percent(counter + 1, n_blocks))) Empty"
            return nothing
        end

        indices_of_used_modules = daf["/ module & is_strong ; block = $(block_name) : index"].array

        fraction_per_used_module_per_block_cell = compute_fraction_per_used_module_per_region_cell(
            daf;
            block_name,
            UMIs_per_gene_per_cell,
            total_UMIs_per_cell,
            region_cell_indices = block_cell_indices,
            indices_of_used_modules,
        )

        if normalize_fractions
            @views mean_fraction_per_module = mean_fraction_per_module_per_block[:, block_index]
            @views std_fraction_per_module = std_fraction_per_module_per_block[:, block_index]
            mean_fraction_per_used_module = mean_fraction_per_module[indices_of_used_modules]
            std_fraction_per_used_module = std_fraction_per_module[indices_of_used_modules]
            fraction_per_used_module_per_block_cell =
                (fraction_per_used_module_per_block_cell .- mean_fraction_per_used_module) ./
                std_fraction_per_used_module
        end

        n_block_clusters = max(Int(round(n_block_cells / mean_metacell_cells_per_block[block_index])), 1)
        if n_block_clusters == 1
            n_block_cells = length(block_cell_indices)
            cluster_index_per_block_cell = fill(1, n_block_cells)
            cluster_sizes = [n_block_cells]

        else
            kmeans_result = kmeans_with_sizes(
                fraction_per_used_module_per_block_cell,
                n_block_clusters;
                min_cluster_size = min_cells_in_metacell,
                max_cluster_size = mean_metacell_cells_per_block[block_index] * 2,
                kmeans_rounds,
                block_name,
                rng,
            )

            cluster_index_per_block_cell = assignments(kmeans_result)
            cluster_sizes = counts(kmeans_result)
        end

        return local_clusters_per_block[block_index] = LocalClusters(;
            block_cell_indices,
            cluster_index_per_block_cell,
            is_too_small_per_cluster = cluster_sizes .< min_cells_in_metacell,
        )
    end

    return local_clusters_per_block
end

function combine_local_clusters(;
    local_clusters_per_block::AbstractVector{Maybe{LocalClusters}},
    name_per_block::AbstractVector{<:AbstractString},
    min_cells_in_metacell::Integer,
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
                    @assert length(cell_indices_of_new_metacell) >= min_cells_in_metacell
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
    UMIs_per_gene_per_cell::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    min_cells_in_metacell::Integer,
    kmeans_rounds::Integer,
    normalize_fractions::Bool,
    rng::AbstractRNG,
)::LocalClusters
    block_cell_indices = daf["/ cell & metacell ?? => block = $(block_name) : index"].array
    indices_of_used_modules = daf["/ module & is_strong ; block = $(block_name) : index"].array

    fraction_per_used_module_per_block_cell = compute_fraction_per_used_module_per_region_cell(
        daf;
        block_name,
        UMIs_per_gene_per_cell,
        total_UMIs_per_cell,
        region_cell_indices = block_cell_indices,
        indices_of_used_modules,
    )

    if normalize_fractions  # TODOX REFACTOR
        n_used_modules = size(fraction_per_used_module_per_block_cell, 1)
        mean_fraction_per_used_module = vec(mean(fraction_per_used_module_per_block_cell; dims = 2))
        @assert_vector(mean_fraction_per_used_module, n_used_modules)
        std_fraction_per_used_module = vec(std(fraction_per_used_module_per_block_cell; dims = 2))
        @assert_vector(std_fraction_per_used_module, n_used_modules)
        std_fraction_per_used_module[std_fraction_per_used_module .== 0] .= 1.0
        fraction_per_used_module_per_block_cell =
            (fraction_per_used_module_per_block_cell .- mean_fraction_per_used_module) ./ std_fraction_per_used_module
    end

    n_block_metacells = daf["/ metacell & block = $(block_name) %> Count"]
    if n_block_metacells == 1
        n_block_cells = length(block_cell_indices)
        cluster_index_per_block_cell = fill(1, n_block_cells)
        cluster_sizes = [n_block_cells]

    else
        mean_cells_in_metacell = length(block_cell_indices) / n_block_metacells

        kmeans_result = kmeans_with_sizes(
            fraction_per_used_module_per_block_cell,
            n_block_metacells;
            min_cluster_size = min_cells_in_metacell,
            max_cluster_size = 2 * mean_cells_in_metacell,
            kmeans_rounds,
            block_name,
            rng,
        )

        cluster_index_per_block_cell = assignments(kmeans_result)
        cluster_sizes = counts(kmeans_result)
    end

    is_too_small_per_cluster = cluster_sizes .< min_cells_in_metacell
    @debug "- Block: $(block_name) Old metacells: $(n_block_metacells) New metacells: $(length(cluster_sizes)) Too small: $(sum(cluster_sizes .< min_cells_in_metacell))"
    return LocalClusters(; block_cell_indices, cluster_index_per_block_cell, is_too_small_per_cluster)
end

function compute_fraction_per_used_module_per_region_cell(
    daf::DafReader;
    block_name::AbstractString,
    UMIs_per_gene_per_cell::AbstractMatrix{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    region_cell_indices::AbstractVector{<:Integer},
    indices_of_used_modules::AbstractVector{<:Integer},
)::Matrix{Float32}
    n_region_cells = length(region_cell_indices)
    total_UMIs_per_region_cell = total_UMIs_per_cell[region_cell_indices]

    n_used_modules = length(indices_of_used_modules)

    is_strong_per_module = daf["/ module / block = $(block_name) : is_strong"].array
    module_index_per_gene = daf["/ gene / block = $(block_name) : module ?? 0 => index"].array
    is_strong_per_gene =
        [module_index > 0 && is_strong_per_module[module_index] for module_index in module_index_per_gene]
    n_strong_genes = sum(is_strong_per_gene)
    module_index_per_strong_gene = module_index_per_gene[is_strong_per_gene]

    @views UMIs_per_gene_per_region_cell = UMIs_per_gene_per_cell[:, region_cell_indices]
    @views UMIs_per_strong_gene_per_region_cell = UMIs_per_gene_per_region_cell[is_strong_per_gene, :]
    @assert_matrix(UMIs_per_strong_gene_per_region_cell, n_strong_genes, n_region_cells)

    fraction_per_strong_gene_per_region_cell =
        UMIs_per_strong_gene_per_region_cell ./ transpose(total_UMIs_per_region_cell)
    @assert_matrix(fraction_per_strong_gene_per_region_cell, n_strong_genes, n_region_cells)
    fraction_per_region_cell_per_strong_gene = flipped(fraction_per_strong_gene_per_region_cell)
    @assert_matrix(fraction_per_region_cell_per_strong_gene, n_region_cells, n_strong_genes)

    fraction_per_region_cell_per_used_module = Matrix{Float32}(undef, n_region_cells, n_used_modules)
    for (used_module_position, used_module_index) in enumerate(indices_of_used_modules)
        module_strong_gene_positions = findall(module_index_per_strong_gene .== used_module_index)
        n_module_strong_genes = length(module_strong_gene_positions)
        @assert n_module_strong_genes > 0

        @views fraction_per_region_cell_per_module_strong_gene =
            fraction_per_region_cell_per_strong_gene[:, module_strong_gene_positions]
        @assert_matrix(fraction_per_region_cell_per_module_strong_gene, n_region_cells, n_module_strong_genes)

        fraction_per_region_cell_per_used_module[:, used_module_position] =
            vec(sum(fraction_per_region_cell_per_module_strong_gene; dims = 2))
    end

    return flipped(fraction_per_region_cell_per_used_module)
end

function kmeans_with_sizes(
    values_of_points::AbstractMatrix{<:AbstractFloat},
    initial_k::Integer;
    min_cluster_size::Real,
    max_cluster_size::Real,
    kmeans_rounds::Integer,
    block_name::AbstractString,
    rng::AbstractRNG,
)::KmeansResult
    kmeans_result = nothing
    max_k = size(values_of_points, 2)
    k = min(initial_k, max_k)
    centers = nothing
    cluster_sizes = nothing

    do_next_phase = true
    retries = 0
    has_kmeans_result = false

    best_kmeans_result = nothing
    best_k = nothing
    best_n_too_small = nothing

    while do_next_phase
        #@debug "TODOX KMEANS $(block_name) P OLD: $(initial_k) K: $(k)"
        do_next_phase = false

        while true
            if has_kmeans_result
                has_kmeans_result = false
                cluster_sizes = counts(kmeans_result)
                #todox_n_too_small = sum(cluster_sizes .< min_cluster_size)
                #todox_n_too_large = sum(cluster_sizes .> max_cluster_size)
                @assert length(cluster_sizes) == k
                #@debug "TODOX KMEANS $(block_name) - A OLD: $(initial_k) K: $(k) small: $(todox_n_too_small) large: $(todox_n_too_large)"
            else
                kmeans_result = kmeans_in_rounds(values_of_points, k; centers, rounds = kmeans_rounds, rng)
                cluster_sizes = counts(kmeans_result)
                @assert length(cluster_sizes) == k
                #todox_n_too_small = sum(cluster_sizes .< min_cluster_size)
                #todox_n_too_large = sum(cluster_sizes .> max_cluster_size)
                #@debug "TODOX KMEANS $(block_name) - B OLD: $(initial_k) K: $(k) small: $(todox_n_too_small) large: $(todox_n_too_large)"
            end

            if k < initial_k
                need_split = min(initial_k - k, k)
                descending_sizes = sort(cluster_sizes; rev = true) # TODOX: quantile
                effective_max_cluster_size = min(max_cluster_size, descending_sizes[need_split] - 1)
                #@debug "TODOX KMEANS $(block_name)   descending_sizes: $(descending_sizes)"
                #@debug "TODOX KMEANS $(block_name)   max_cluster_size: $(max_cluster_size) effective_max_cluster_size: $(effective_max_cluster_size)"
            else
                effective_max_cluster_size = max_cluster_size
            end
            effective_max_cluster_size = max(effective_max_cluster_size, 2 * min_cluster_size)

            #todox_n_too_large = sum(cluster_sizes .> effective_max_cluster_size)
            #@debug "TODOX KMEANS $(block_name)   LARGE: $(todox_n_too_large)"

            largest_cluster_size = maximum(cluster_sizes)
            if largest_cluster_size <= effective_max_cluster_size
                break
            end

            n_large_clusters = sum(cluster_sizes .> effective_max_cluster_size)
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

        n_too_small = 0
        while k > 1
            @assert length(cluster_sizes) == k
            smallest_cluster_index = argmin(cluster_sizes)
            smallest_cluster_size = cluster_sizes[smallest_cluster_index]

            n_too_small = sum(cluster_sizes .< min_cluster_size)
            n_too_large = sum(cluster_sizes .> max_cluster_size)
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
                #@debug "TODOX KMEANS $(block_name)   !X OLD: $(initial_k) K: $(k) small: $(n_too_small)"
                best_kmeans_result = kmeans_result
                best_k = k
                best_n_too_small = n_too_small
            end

            if n_too_small == 0
                break
            end

            centers = kmeans_result.centers[:, 1:k .!= smallest_cluster_index]
            k -= 1

            merged_kmeans_result = kmeans_in_rounds(values_of_points, k; centers, rounds = kmeans_rounds, rng)
            cluster_sizes = counts(merged_kmeans_result)
            @assert length(cluster_sizes) == k
            #n_too_small = sum(cluster_sizes .< min_cluster_size)
            #n_too_large = sum(cluster_sizes .> max_cluster_size)
            #@debug "TODOX KMEANS $(block_name) - C OLD: $(initial_k) K: $(k) small: $(n_too_small) large: $(n_too_large)"
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
                #@debug "TODOX KMEANS $(block_name)   !Y OLD: $(initial_k) K: $(k) small: $(n_too_small)"
                best_kmeans_result = merged_kmeans_result
                best_k = k
                best_n_too_small = n_too_small
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

        @assert length(counts(kmeans_result)) == k
    end

    @assert best_kmeans_result !== nothing
    return best_kmeans_result
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
    ],
    data = [
        metacell_block_vector(RequiredInput),
        block_metacell_module_linear_fraction_tensor(RequiredInput),
        block_module_is_strong_matrix(RequiredInput),
        block_module_neighborhood_mean_linear_fraction_matrix(RequiredInput),
        block_module_neighborhood_std_linear_fraction_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_radius_vector(GuaranteedOutput),
    ],
) function compute_metacells_radius!(
    daf::DafWriter;
    metacell_radius_quantile::AbstractFloat = 0.95,
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_metacells = axis_length(daf, "metacell")
    n_modules = axis_length(daf, "module")

    name_per_block = axis_vector(daf, "block")
    block_index_per_metacell = daf["/ metacell : block => index"].array
    metacell_index_per_cell = daf["/ cell : metacell ?? 0 => index"].array
    module_index_per_gene_per_block = daf["/ gene / block : module ?? 0 => index"].array

    is_strong_per_module_per_block = get_matrix(daf, "module", "block", "is_strong").array
    mean_linear_fraction_per_module_per_block =
        get_matrix(daf, "module", "block", "neighborhood_mean_linear_fraction").array
    std_linear_fraction_per_module_per_block =
        get_matrix(daf, "module", "block", "neighborhood_std_linear_fraction").array

    radius_per_metacell = Vector{Float32}(undef, n_metacells)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    progress_counter = Atomic{Int}(0)
    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]

        indices_of_block_metacells = findall(block_index_per_metacell .== block_index)
        n_block_metacells = length(indices_of_block_metacells)
        @assert n_block_metacells > 0

        @views mean_linear_fraction_per_module = mean_linear_fraction_per_module_per_block[:, block_index]
        @views std_linear_fraction_per_module = std_linear_fraction_per_module_per_block[:, block_index]
        linear_fraction_per_module_per_metacell =
            get_matrix(daf, "module", "metacell", "$(block_name)_linear_fraction").array

        @views linear_fraction_per_module_per_block_metacell =
            linear_fraction_per_module_per_metacell[:, indices_of_block_metacells]
        normalized_fraction_per_module_per_block_metacell =
            (linear_fraction_per_module_per_block_metacell .- mean_linear_fraction_per_module) ./
            std_linear_fraction_per_module_per_block[:, block_index]

        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        @views is_strong_per_module = is_strong_per_module_per_block[:, block_index]
        indices_of_strong_modules = findall(is_strong_per_module)
        n_strong_modules = length(indices_of_strong_modules)
        @assert n_strong_modules > 0

        indices_of_genes_per_module = Vector{AbstractVector{<:Integer}}(undef, n_modules)
        for module_index in indices_of_strong_modules
            indices_of_genes_per_module[module_index] = findall(module_index_per_gene .== module_index)
            n_module_genes = length(indices_of_genes_per_module[module_index])
            @assert n_module_genes > 0
        end

        for (metacell_position, metacell_index) in enumerate(indices_of_block_metacells)
            indices_of_metacell_cells = findall(metacell_index_per_cell .== metacell_index)
            n_metacell_cells = length(indices_of_metacell_cells)
            @assert n_metacell_cells > 0

            @views normalized_fraction_per_module =
                normalized_fraction_per_module_per_block_metacell[:, metacell_position]

            distance_per_metacell_cell = zeros(Float32, n_metacell_cells)

            for module_index in indices_of_strong_modules
                indices_of_module_genes = indices_of_genes_per_module[module_index]

                UMIs_per_metacell_cell_per_module_gene =
                    UMIs_per_cell_per_gene[indices_of_metacell_cells, indices_of_module_genes]

                normalized_fraction_per_metacell_cell =
                    (
                        (
                            sum(UMIs_per_metacell_cell_per_module_gene; dims = 2) ./
                            total_UMIs_per_cell[indices_of_metacell_cells]
                        ) .- mean_linear_fraction_per_module[module_index]
                    ) ./ std_linear_fraction_per_module[module_index]

                normalized_fraction_of_metacell = normalized_fraction_per_module[module_index]
                @assert !isnan(normalized_fraction_of_metacell)

                distance_per_metacell_cell .+=
                    (normalized_fraction_per_metacell_cell .- normalized_fraction_of_metacell) .^ 2
            end

            distance_per_metacell_cell .= sqrt.(distance_per_metacell_cell)
            radius_per_metacell[metacell_index] = quantile(distance_per_metacell_cell, metacell_radius_quantile)
        end

        counter = atomic_add!(progress_counter, 1)
        @debug "- Block: $(block_name) ($(percent(counter + 1, n_blocks))) mean: $(mean(radius_per_metacell[indices_of_block_metacells]))"
    end

    @debug "Mean metacell radius: $(mean(radius_per_metacell))"
    set_vector!(daf, "metacell", "radius", radius_per_metacell)
    return nothing
end

end  # module
