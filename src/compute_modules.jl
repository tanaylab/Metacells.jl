"""
.
Group genes in each environment to modules, each associated with some anchor gene.
"""
module ComputeModules

export compute_blocks_modules!

using Base.Threads
using Clustering
using DataAxesFormats
using Distances
using TanayLabUtilities
using Random
using StatsBase
using Serialization

using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_distinct_per_gene_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_marker_per_gene_per_block
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_module_status_per_gene_per_block
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.vector_of_anchor_per_module
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_metacell_per_cell

"""
    compute_blocks_modules!(
        daf::DafWriter;
        max_clusters::Integer = 24,
        min_member_correlation::AbstractFloat = 0.5,
        min_orphan_correlation::AbstractFloat = 0.6,
        min_strong_UMIs::Integer = 8,
        min_strong_cells::Integer = 12,
        kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
        module_status::Bool = false,
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = false,
    )::Nothing

Compute and set [`vector_of_anchor_per_module`](@ref), [`matrix_of_is_found_per_module_per_block`](@ref),
[`matrix_of_module_per_gene_per_block`](@ref), and (if `module_status` is specified),
[`matrix_of_module_status_per_gene_per_block`](@ref). To do this, for each neighborhood:

 1. We cluster all the neighborhood lateral marker genes into `max_clusters` using K-means of the normalized
    gene expression levels (z-score using the mean and standard deviation in the neighborhood).
 2. We similarly cluster all the non-lateral genes which are correlated with a skeleton in the neighborhood.
 3. We identify as lateral-ish any of these genes which are more correlated to a center of some lateral cluster than to
    the center of their own cluster.
 4. We then only keep in each cluster the non-lateral-ish genes whose correlation with the center of the cluster is at
    least `min_member_correlation`. Some clusters may be dissolved as a result.
 5. We pick an anchor for each cluster - the skeleton gene closest to its center (if any).
 6. We look at the minimal number of UMIs in the `min_strong_cells` of each cluster. If the weakest cluster
    does not have at least `min_strong_UMIs`, we dissolve it.
 7. If any gene became unclustered above, and it has a correlation of at least `min_orphan_correlation` with
    any of the anchor genes, we add it to that anchor gene's cluster.
"""
@logged :mcs_ops @computation Contract(
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(CreatedOutput),
    ],
    data = [
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(RequiredInput),
        matrix_of_is_neighborhood_distinct_per_gene_per_block(RequiredInput),
        matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_anchor_per_module(CreatedOutput),
        matrix_of_is_found_per_module_per_block(CreatedOutput),
        matrix_of_module_per_gene_per_block(CreatedOutput),
        matrix_of_module_status_per_gene_per_block(OptionalInput),
    ],
) function compute_blocks_modules!(
    daf::DafWriter;
    max_clusters::Integer = 24,
    min_member_correlation::AbstractFloat = 0.5,
    min_orphan_correlation::AbstractFloat = 0.6,
    min_strong_UMIs::Integer = 8,
    min_strong_cells::Integer = 12,
    kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
    module_status::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= min_member_correlation <= min_orphan_correlation <= 1
    @assert min_strong_UMIs >= 0
    @assert min_strong_cells >= 0

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    name_per_block = axis_vector(daf, "block")
    name_per_gene = axis_vector(daf, "gene")

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array
    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array
    is_neighborhood_marker_per_gene_per_block = get_matrix(daf, "gene", "block", "is_neighborhood_marker").array
    is_correlated_with_skeletons_in_neighborhood_per_gene_per_block =
        get_matrix(daf, "gene", "block", "is_correlated_with_skeleton_in_neighborhood").array
    is_neighborhood_distinct_per_gene_per_block = get_matrix(daf, "gene", "block", "is_neighborhood_distinct").array
    block_index_per_metacell = daf["@ metacell : block : index"].array
    block_index_per_cell = daf["@ cell : metacell ?? 0 : block : index"].array
    log_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "log_linear_fraction").array

    genes_indices_of_anchor_index_per_block = Vector{Dict{Int, Vector{Int}}}(undef, n_blocks)

    if module_status
        module_status_per_gene_per_block = fill("", n_genes, n_blocks)
    else
        module_status_per_gene_per_block = nothing
    end

    parallel_loop_with_rng(
        1:n_blocks;
        rng,
        progress = DebugProgress(n_blocks; group = :mcs_details, desc = "compute_blocks_modules"),
    ) do block_index, rng
        if module_status_per_gene_per_block === nothing
            module_status_per_gene = nothing
        else
            @views module_status_per_gene = module_status_per_gene_per_block[:, block_index]
        end
        genes_indices_of_anchor_index_per_block[block_index] = compute_block_modules!(
            daf;
            block_index,
            max_clusters,
            min_member_correlation,
            min_orphan_correlation,
            min_strong_UMIs,
            min_strong_cells,
            kmeans_rounds,
            rng,
            module_status_per_gene,
            UMIs_per_cell_per_gene,
            is_lateral_per_gene,
            is_skeleton_per_gene,
            is_in_neighborhood_per_other_block_per_base_block,
            is_neighborhood_marker_per_gene_per_block,
            is_correlated_with_skeletons_in_neighborhood_per_gene_per_block,
            is_neighborhood_distinct_per_gene_per_block,
            block_index_per_metacell,
            log_fraction_per_metacell_per_gene,
            block_index_per_cell,
        )
        return nothing
    end

    is_anchor_per_gene = zeros(Bool, n_genes)
    for genes_indices_of_anchor_index in genes_indices_of_anchor_index_per_block
        for gene_index in keys(genes_indices_of_anchor_index)
            is_anchor_per_gene[gene_index] = true
        end
    end

    anchor_index_per_module = findall(is_anchor_per_gene)
    anchor_per_module = name_per_gene[anchor_index_per_module]
    name_per_module = anchor_per_module .* ".MOD"
    n_modules = length(name_per_module)

    add_axis!(daf, "module", name_per_module; overwrite)
    set_vector!(daf, "module", "gene.anchor", anchor_per_module; overwrite)

    is_found_per_module_per_block = zeros(Bool, n_modules, n_blocks)
    module_per_gene_per_block = fill("", n_genes, n_blocks)

    for (block_index, genes_indices_of_anchor_index) in enumerate(genes_indices_of_anchor_index_per_block)
        block_name = name_per_block[block_index]
        for (module_index, anchor_index) in enumerate(anchor_index_per_module)
            genes_indices_of_anchor = get(genes_indices_of_anchor_index, anchor_index, nothing)
            if genes_indices_of_anchor !== nothing
                module_name = name_per_module[module_index]
                module_per_gene_per_block[genes_indices_of_anchor, block_index] .= module_name
                is_found_per_module_per_block[module_index, block_index] = true
                @debug "Block: $(block_name) anchor: $(name_per_gene[anchor_index]) genes: [ $(join(sort(name_per_gene[genes_indices_of_anchor]), ", ")) ]" _group =
                    :mcs_details
            end
        end
    end

    set_matrix!(daf, "module", "block", "is_found", bestify(is_found_per_module_per_block); overwrite)
    set_matrix!(daf, "gene", "block", "module", module_per_gene_per_block; overwrite)
    if module_status_per_gene_per_block !== nothing
        set_matrix!(daf, "gene", "block", "module_status", module_status_per_gene_per_block; overwrite)
    end

    return nothing
end

function compute_block_modules!(
    daf::DafWriter;
    block_index::Integer,
    max_clusters::Integer,
    min_member_correlation::AbstractFloat,
    min_orphan_correlation::AbstractFloat,
    min_strong_UMIs::Integer,
    min_strong_cells::Integer,
    kmeans_rounds::Integer,
    rng::AbstractRNG,
    module_status_per_gene::Maybe{AbstractVector{<:AbstractString}},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Integer},
    is_lateral_per_gene::Union{AbstractVector{Bool}, BitVector},
    is_skeleton_per_gene::Union{AbstractVector{Bool}, BitVector},
    is_in_neighborhood_per_other_block_per_base_block::Union{AbstractMatrix{Bool}, BitMatrix},
    is_neighborhood_marker_per_gene_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    is_correlated_with_skeletons_in_neighborhood_per_gene_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    is_neighborhood_distinct_per_gene_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    block_index_per_metacell::AbstractVector{<:Integer},
    log_fraction_per_metacell_per_gene::AbstractMatrix{<:AbstractFloat},
    block_index_per_cell::AbstractVector{<:Integer},
)::Dict{Int, Vector{Int}}
    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    lateral_cluster_per_gene = zeros(UInt32, n_genes)

    local z_score_per_neighborhood_metacell_per_lateral_cluster
    local is_in_neighborhood_per_other_block
    local indices_of_neighborhood_metacells

    flame_timed("lateral_clusters") do
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        indices_of_neighborhood_metacells = findall(is_in_neighborhood_per_other_block[block_index_per_metacell])

        @views is_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene_per_block[:, block_index]
        is_lateral_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene .& is_lateral_per_gene
        indices_of_lateral_neighborhood_markers = findall(is_lateral_neighborhood_marker_per_gene)
        n_lateral_neighborhood_marker_genes = length(indices_of_lateral_neighborhood_markers)

        @views log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene =
            log_fraction_per_metacell_per_gene[
                indices_of_neighborhood_metacells,
                indices_of_lateral_neighborhood_markers,
            ]

        mean_log_fraction_per_lateral_neighborhood_marker_gene =
            vec(mean(log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene; dims = 1))  # NOLINT
        @assert_vector(mean_log_fraction_per_lateral_neighborhood_marker_gene, n_lateral_neighborhood_marker_genes)

        std_log_fraction_per_lateral_neighborhood_marker_gene = vec(
            std(  # NOLINT
                log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene;
                mean = transpose(mean_log_fraction_per_lateral_neighborhood_marker_gene),
                dims = 1,
            ),
        )
        @assert_vector(std_log_fraction_per_lateral_neighborhood_marker_gene, n_lateral_neighborhood_marker_genes)

        z_score_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene =
            (
                log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene .-
                transpose(mean_log_fraction_per_lateral_neighborhood_marker_gene)
            ) ./ transpose(std_log_fraction_per_lateral_neighborhood_marker_gene)

        kmeans_results = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                z_score_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene,
                max_clusters;
                rounds = kmeans_rounds,
                rng,
            )
        end

        lateral_cluster_per_gene[indices_of_lateral_neighborhood_markers] .= kmeans_results.assignments

        z_score_per_neighborhood_metacell_per_lateral_cluster = kmeans_results.centers  # NOLINT
        return nothing
    end

    local z_score_per_neighborhood_metacell_per_cluster
    local cluster_index_per_local_gene
    local indices_of_local_genes
    local indices_of_neighborhood_cells
    local is_skeleton_per_local_gene
    local n_local_genes
    local n_neighborhood_cells
    local z_score_per_neighborhood_metacell_per_local_gene

    flame_timed("local_clusters") do
        indices_of_neighborhood_cells = findall(
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1)),
        )
        n_neighborhood_cells = length(indices_of_neighborhood_cells)
        @assert n_neighborhood_cells > min_strong_cells

        @views is_correlated_with_skeletons_in_neighborhood_per_gene =
            is_correlated_with_skeletons_in_neighborhood_per_gene_per_block[:, block_index]
        is_local_per_gene = is_correlated_with_skeletons_in_neighborhood_per_gene .& .!is_lateral_per_gene
        indices_of_local_genes = findall(is_local_per_gene)
        n_local_genes = length(indices_of_local_genes)

        is_skeleton_per_local_gene = is_skeleton_per_gene[indices_of_local_genes]  # NOLINT

        log_fraction_per_neighborhood_metacell_per_local_gene =
            log_fraction_per_metacell_per_gene[indices_of_neighborhood_metacells, indices_of_local_genes]

        mean_log_fraction_per_local_gene = vec(mean(log_fraction_per_neighborhood_metacell_per_local_gene; dims = 1))  # NOLINT
        @assert_vector(mean_log_fraction_per_local_gene, n_local_genes)

        std_log_fraction_per_local_gene = vec(
            std(  # NOLINT
                log_fraction_per_neighborhood_metacell_per_local_gene;
                mean = transpose(mean_log_fraction_per_local_gene),
                dims = 1,
            ),
        )
        @assert_vector(std_log_fraction_per_local_gene, n_local_genes)

        z_score_per_neighborhood_metacell_per_local_gene =
            (log_fraction_per_neighborhood_metacell_per_local_gene .- transpose(mean_log_fraction_per_local_gene)) ./
            transpose(std_log_fraction_per_local_gene)

        kmeans_results = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                z_score_per_neighborhood_metacell_per_local_gene,
                max_clusters;
                rounds = kmeans_rounds,
                rng,
            )
        end

        cluster_index_per_local_gene = kmeans_results.assignments
        @assert_vector(cluster_index_per_local_gene, n_local_genes)

        z_score_per_neighborhood_metacell_per_cluster = kmeans_results.centers  # NOLINT
        if module_status_per_gene !== nothing
            module_status_per_gene[indices_of_local_genes] .= "cluster(" .* string.(cluster_index_per_local_gene) .* ")"
        end
        return nothing
    end

    local is_lateral_ish_per_local_gene

    flame_timed("identify_lateral_ish_genes") do
        is_lateral_ish_per_local_gene = zeros(Bool, n_local_genes)

        is_neighborhood_distinct_per_local_gene =
            is_neighborhood_distinct_per_gene_per_block[indices_of_local_genes, block_index]

        for local_gene_position in 1:n_local_genes
            @views z_score_per_neighborhood_metacell_of_local_gene =
                z_score_per_neighborhood_metacell_per_local_gene[:, local_gene_position]
            cluster_index = cluster_index_per_local_gene[local_gene_position]
            @views cluster_log_fraction_per_neighborhood_metacell =
                z_score_per_neighborhood_metacell_per_cluster[:, cluster_index]
            cluster_correlation = zero_cor_between_vectors(
                z_score_per_neighborhood_metacell_of_local_gene,
                cluster_log_fraction_per_neighborhood_metacell,
            )
            correlation_per_lateral_cluster = zero_cor_between_vector_and_matrix_columns(
                z_score_per_neighborhood_metacell_of_local_gene,
                z_score_per_neighborhood_metacell_per_lateral_cluster,
            )
            lateral_cluster_index = argmax(correlation_per_lateral_cluster)
            lateral_correlation = correlation_per_lateral_cluster[lateral_cluster_index]
            qualifier = ""
            if is_neighborhood_distinct_per_local_gene[local_gene_position]
                qualifier *= " distinct"
            end
            if is_skeleton_per_local_gene[local_gene_position]
                qualifier *= " skeleton"
            end
            if lateral_correlation > cluster_correlation
                if qualifier == ""
                    is_lateral_ish_per_local_gene[local_gene_position] = true
                else
                    qualifier *= " kept"
                end
            end
        end

        if module_status_per_gene !== nothing
            module_status_per_gene[indices_of_local_genes[is_lateral_ish_per_local_gene]] .*= ",lateral_ish"
        end

        return nothing
    end

    local anchor_local_position_per_cluster
    local cluster_index_of_anchor_local_position

    flame_timed("correlated_pertinent_gene_modules") do
        n_clusters = size(z_score_per_neighborhood_metacell_per_cluster, 2)
        @assert 0 < n_clusters <= max_clusters

        anchor_local_position_per_cluster = fill(UInt32(0), n_clusters)
        cluster_index_of_anchor_local_position = Dict{UInt32, UInt32}()

        for cluster_index in 1:n_clusters
            @views z_score_per_neighborhood_metacell_of_cluster =
                z_score_per_neighborhood_metacell_per_cluster[:, cluster_index]
            is_in_cluster_per_local_gene = cluster_index_per_local_gene .== cluster_index
            n_cluster_genes = sum(is_in_cluster_per_local_gene)
            @assert n_cluster_genes > 0
            is_lateral_ish_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& is_lateral_ish_per_local_gene

            cluster_index_per_local_gene[is_lateral_ish_in_cluster_per_local_gene] .= -1
            is_pertinent_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& .!is_lateral_ish_per_local_gene
            n_pertinent_in_cluster_local_genes = sum(is_pertinent_in_cluster_per_local_gene)

            is_skeleton_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& is_skeleton_per_local_gene
            local_positions_of_skeletons_in_cluster = findall(is_skeleton_in_cluster_per_local_gene)
            n_cluster_skeletons = length(local_positions_of_skeletons_in_cluster)

            if n_cluster_skeletons == 0
                cluster_index_per_local_gene[is_pertinent_in_cluster_per_local_gene] .= 0
                if module_status_per_gene !== nothing
                    module_status_per_gene[indices_of_local_genes[is_pertinent_in_cluster_per_local_gene]] .*= ",orphan_cluster"
                end
                continue
            end

            z_score_per_neighborhood_metacell_per_pertinent_local_gene =
                z_score_per_neighborhood_metacell_per_local_gene[:, is_pertinent_in_cluster_per_local_gene]
            correlation_with_cluster_per_pertinent_local_gene = zero_cor_between_vector_and_matrix_columns(
                z_score_per_neighborhood_metacell_of_cluster,
                z_score_per_neighborhood_metacell_per_pertinent_local_gene,
            )
            @assert_vector(correlation_with_cluster_per_pertinent_local_gene, n_pertinent_in_cluster_local_genes)

            is_correlated_pertinent_in_cluster_per_local_gene = zeros(Bool, n_local_genes)
            is_correlated_pertinent_in_cluster_per_local_gene[is_pertinent_in_cluster_per_local_gene] .=
                correlation_with_cluster_per_pertinent_local_gene .>= min_member_correlation

            is_uncorrelated_pertinent_in_cluster_per_local_gene =
                is_pertinent_in_cluster_per_local_gene .& .!is_correlated_pertinent_in_cluster_per_local_gene
            cluster_index_per_local_gene[is_uncorrelated_pertinent_in_cluster_per_local_gene] .= 0

            n_correlated_pertinent_local_genes = sum(is_correlated_pertinent_in_cluster_per_local_gene)
            n_uncorrelated_pertinent_cluster_genes = sum(is_uncorrelated_pertinent_in_cluster_per_local_gene)

            @assert n_correlated_pertinent_local_genes == sum(cluster_index_per_local_gene .== cluster_index)
            if n_correlated_pertinent_local_genes == 0
                if module_status_per_gene !== nothing
                    module_status_per_gene[indices_of_local_genes[is_uncorrelated_pertinent_in_cluster_per_local_gene]] .*= ",uncorrelated_gene,uncorrelated_cluster"
                end
                if n_uncorrelated_pertinent_cluster_genes == 0
                    reason = "all-lateral"
                else
                    reason = "all-uncorrelated"
                end
                continue
            end

            if module_status_per_gene !== nothing
                module_status_per_gene[indices_of_local_genes[is_uncorrelated_pertinent_in_cluster_per_local_gene]] .*= ",uncorrelated_gene"
            end

            if n_cluster_skeletons == 1
                anchor_local_position = local_positions_of_skeletons_in_cluster[1]

            else
                z_score_per_neighborhood_metacell_per_skeleton =
                    @views z_score_per_neighborhood_metacell_per_local_gene[:, local_positions_of_skeletons_in_cluster]
                distance_from_cluster_per_skeleton = flame_timed("colwise.Euclidean") do
                    return colwise(
                        Euclidean(),
                        z_score_per_neighborhood_metacell_of_cluster,
                        z_score_per_neighborhood_metacell_per_skeleton,
                    )
                end
                anchor_local_position =
                    local_positions_of_skeletons_in_cluster[argmin(distance_from_cluster_per_skeleton)]
            end

            if module_status_per_gene !== nothing
                module_status_per_gene[indices_of_local_genes[is_correlated_pertinent_in_cluster_per_local_gene]] .*= ",anchor($(name_per_gene[indices_of_local_genes[anchor_local_position]]))"
            end

            anchor_local_position_per_cluster[cluster_index] = anchor_local_position
            cluster_index_of_anchor_local_position[anchor_local_position] = cluster_index
        end

        @assert length(cluster_index_of_anchor_local_position) > 0
    end

    flame_timed("migrate_between_gene_modules") do
        @views UMIs_per_neighborhood_cell_per_gene = UMIs_per_cell_per_gene[indices_of_neighborhood_cells, :]

        while true
            local_positions_of_anchors = collect(keys(cluster_index_of_anchor_local_position))
            cluster_indices_of_anchors = [
                cluster_index_of_anchor_local_position[anchor_local_position] for
                anchor_local_position in local_positions_of_anchors
            ]
            @views z_score_per_neighborhood_metacell_per_anchor =
                z_score_per_neighborhood_metacell_per_cluster[:, cluster_indices_of_anchors]

            for local_gene_index in 1:n_local_genes
                if cluster_index_per_local_gene[local_gene_index] == 0
                    @views z_score_per_neighborhood_metacell_of_local_gene =
                        z_score_per_neighborhood_metacell_per_local_gene[:, local_gene_index]
                    correlation_per_anchor = zero_cor_between_vector_and_matrix_columns(
                        z_score_per_neighborhood_metacell_of_local_gene,
                        z_score_per_neighborhood_metacell_per_anchor,
                    )
                    @assert_vector(correlation_per_anchor, length(local_positions_of_anchors))
                    anchor_position = argmax(correlation_per_anchor)
                    correlation = correlation_per_anchor[anchor_position]
                    anchor_local_position = local_positions_of_anchors[anchor_position]
                    cluster_index = cluster_index_of_anchor_local_position[anchor_local_position]
                    if correlation >= min_orphan_correlation
                        cluster_index_per_local_gene[local_gene_index] = cluster_index
                        if module_status_per_gene !== nothing
                            module_status_per_gene[indices_of_local_genes[local_gene_index]] *= ",joined($(name_per_gene[indices_of_local_genes[anchor_local_position]]))"
                        end
                    end
                end
            end

            if length(cluster_index_of_anchor_local_position) == 1
                break
            end

            weakest_cluster_index = nothing
            weakest_cluster_strong_UMIs = nothing
            for (anchor_local_position, cluster_index) in cluster_index_of_anchor_local_position  # NOLINT
                is_in_cluster_per_local_gene = cluster_index_per_local_gene .== cluster_index
                @assert any(is_in_cluster_per_local_gene)
                @views UMIs_per_neighborhood_cell_per_cluster_local_gene =
                    UMIs_per_neighborhood_cell_per_gene[:, indices_of_local_genes[is_in_cluster_per_local_gene]]
                cluster_UMIs_per_neighborhood_cell =
                    vec(sum(UMIs_per_neighborhood_cell_per_cluster_local_gene; dims = 2))
                cluster_strong_UMIs = quantile(  # NOLINT
                    cluster_UMIs_per_neighborhood_cell,
                    1 - (min_strong_cells - 1) / (n_neighborhood_cells - 1),
                )
                if weakest_cluster_strong_UMIs === nothing || cluster_strong_UMIs < weakest_cluster_strong_UMIs
                    weakest_cluster_index = cluster_index
                    weakest_cluster_strong_UMIs = cluster_strong_UMIs
                end
            end

            @assert weakest_cluster_index !== nothing
            anchor_local_position = anchor_local_position_per_cluster[weakest_cluster_index]

            if weakest_cluster_strong_UMIs >= min_strong_UMIs
                break
            end

            is_in_weakest_cluster_per_local_gene = cluster_index_per_local_gene .== weakest_cluster_index
            if module_status_per_gene !== nothing
                module_status_per_gene[indices_of_local_genes[is_in_weakest_cluster_per_local_gene]] .*= ",weak_cluster"
            end
            cluster_index_per_local_gene[is_in_weakest_cluster_per_local_gene] .= 0
            anchor_local_position_per_cluster[weakest_cluster_index] = 0
            delete!(cluster_index_of_anchor_local_position, anchor_local_position)
        end
    end

    local genes_indices_of_anchor_index
    flame_timed("collect_gene_modules") do
        genes_indices_of_anchor_index = Dict{Int, Vector{Int}}()

        for (anchor_local_position, cluster_index) in cluster_index_of_anchor_local_position
            anchor_index = indices_of_local_genes[anchor_local_position]
            is_in_cluster_per_local_gene = cluster_index_per_local_gene .== cluster_index
            indices_of_cluster_genes = indices_of_local_genes[is_in_cluster_per_local_gene]
            genes_indices_of_anchor_index[anchor_index] = indices_of_cluster_genes
        end
    end

    return genes_indices_of_anchor_index
end

end  # module
