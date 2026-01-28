"""
.
Group genes in each environment to modules, each associated with some anchor gene.
"""
module ComputeModules

export compute_blocks_modules!
export compute_blocks_modules_k!

using Base.Threads
using Clustering
using DataAxesFormats
using Distances
using TanayLabUtilities
using Random
using StatsBase
using Serialization

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_gene_is_neighborhood_marker_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.module_anchor_gene_vector
import Metacells.Contracts.module_axis
import Random.default_rng

"""
    function compute_blocks_modules!(
        daf::DafWriter;
        mean_genes_per_cluster::Integer = $(DEFAULT.mean_genes_per_cluster),
        anchor_correlation_genes::Integer = $(DEFAULT.anchor_correlation_genes),
        min_anchor_correlation_genes_fraction::AbstractFloat = $(DEFAULT.min_anchor_correlation_genes_fraction),
        max_anchor_correlation_genes_fraction::AbstractFloat = $(DEFAULT.max_anchor_correlation_genes_fraction),
        min_anchor_merge_correlation::AbstractFloat = $(DEFAULT.min_anchor_merge_correlation), TODOX
        min_anchor_merge_migration_fraction::AbstractFloat = $(DEFAULT.min_anchor_merge_migration_fraction), TODOX
        min_gene_migration_correlation::AbstractFloat = $(DEFAULT.min_gene_migration_correlation),
        max_anchor_kept_cluster_genes_baseline::Integer = $(DEFAULT.max_anchor_kept_cluster_genes_baseline),
        min_anchor_kept_cluster_genes_likelihood::AbstractFloat = $(DEFAULT.min_anchor_kept_cluster_genes_likelihood),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group genes in each neighborhood to modules, each associated with some anchor gene.

For each block:

 1. Compute the correlation of the expression of the genes across all metacells in the block's neighborhood.
 2. Cluster the genes based on the Euclidean distances between their correlations, using Ward. The number of clusters is
    such that the mean number of genes per cluster is at least `mean_genes_per_cluster`.
 3. For each cluster, pick an anchor gene out of the skeleton genes. For each skeleton, look at its worst correlation
    in the top `anchor_correlation_genes` (but at least `min_anchor_correlation_genes_fraction` and at most
    `max_anchor_correlation_genes_fraction`). Pick as anchor the skeleton gene with the highest result.
 4. For each gene in the cluster, keep it in the anchor's module if it is more correlated with the anchor than
    with any other anchor. Any gene with a correlation of more than `min_gene_migration_correlation` with some
    anchor is moved to the module of its most correlated anchor, regardless of the original cluster it was in.
 5. For each cluster, consider the up to `max_anchor_kept_cluster_genes_baseline` genes most correlated with the
    cluster's anchor. If the fraction of the most correlated `max_anchor_kept_cluster_genes_baseline` genes that were
    kept in the anchor module is less than `min_anchor_kept_cluster_genes_likelihood` times the base fraction of
    1/number of anchors, then discard the anchor.
 6. The result is one module per anchor with the genes that were kept from its original cluster or were migrated to it.
    Some (most) genes will not participate in any module.

TODOX 0.15 parameter
TODOX Selection (100 or migrated)

We name the resulting modules after their anchor gene for interpretability, but add a `.MOD` suffix to clarify this is
the name of a module and not of a gene. It is possible that an anchor in one block will end up being a member of a
different anchor's module in a neighbouring block, or even not to be used in any of the neighbouring block gene modules.

Due to `Daf` limitations, the modules axis must therefore be the union of all the modules from all the blocks, so we
provide an `is_found` mask to specify which module is used for each block.

TODOX status

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(GuaranteedOutput),
    ],
    data = [
        gene_is_lateral_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_is_correlated_with_skeleton_in_neighborhood_matrix(RequiredInput),
        block_gene_is_neighborhood_varied_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        module_anchor_gene_vector(GuaranteedOutput),
        block_module_is_found_matrix(GuaranteedOutput),
        block_gene_module_matrix(GuaranteedOutput),
        block_gene_module_status_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules!(
    daf::DafWriter;
    mean_genes_per_cluster::Integer = 50,
    anchor_correlation_genes::Integer = 30,
    split_correlation_genes::Integer = 10,
    min_anchor_correlation_genes_fraction::AbstractFloat = 0.05,
    max_anchor_correlation_genes_fraction::AbstractFloat = 0.5,
    min_anchor_merge_correlation::AbstractFloat = 0.8,
    min_anchor_merge_migration_fraction::AbstractFloat = 0.5,
    min_orphan_cluster_correlation::AbstractFloat = 0.8,
    min_gene_migration_correlation::AbstractFloat = 0.6,
    max_anchor_kept_cluster_genes_baseline::Integer = 100,
    min_anchor_kept_genes_fraction::AbstractFloat = 0.15,
    min_anchor_kept_cluster_genes_likelihood::AbstractFloat = 1.5,
    min_split_anchor_correlation::AbstractFloat = 0.6,
    min_split_extra_correlation::AbstractFloat = 0.2,
    max_cluster_laterals_fraction::AbstractFloat = 0.25,
    overwrite::Bool = false,
)::Nothing
    @assert mean_genes_per_cluster > 0
    @assert 0 <= min_anchor_correlation_genes_fraction <= max_anchor_correlation_genes_fraction <= 1
    @assert 0 <= min_anchor_merge_correlation::AbstractFloat <= 1
    @assert 0 <= min_anchor_merge_migration_fraction::AbstractFloat <= 1
    @assert 0 <= min_orphan_cluster_correlation::AbstractFloat <= 1
    @assert 0 < split_correlation_genes < anchor_correlation_genes
    @assert 0 <= min_gene_migration_correlation <= 1
    @assert max_anchor_kept_cluster_genes_baseline > 0
    @assert min_anchor_kept_cluster_genes_likelihood >= 0
    @assert 0 <= max_cluster_laterals_fraction <= 1
    @assert 0 <= min_split_extra_correlation <= min_split_anchor_correlation <= 1

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    name_per_block = axis_vector(daf, "block")
    name_per_gene = axis_vector(daf, "gene")

    genes_indices_of_anchor_index_per_block = Vector{Dict{Int, Vector{Int}}}(undef, n_blocks)

    module_status_per_gene_per_block = fill("", n_genes, n_blocks)

    parallel_loop_wo_rng(1:n_blocks) do block_index
        block_name = name_per_block[block_index]
        @views module_status_per_gene = module_status_per_gene_per_block[:, block_index]
        genes_indices_of_anchor_index_per_block[block_index] = compute_block_modules!(
            daf;
            block_name,
            mean_genes_per_cluster,
            anchor_correlation_genes,
            split_correlation_genes,
            min_anchor_correlation_genes_fraction,
            max_anchor_correlation_genes_fraction,
            min_anchor_merge_correlation,
            min_anchor_merge_migration_fraction,
            min_orphan_cluster_correlation,
            min_gene_migration_correlation,
            max_anchor_kept_cluster_genes_baseline,
            min_anchor_kept_genes_fraction,
            min_anchor_kept_cluster_genes_likelihood,
            min_split_anchor_correlation,
            min_split_extra_correlation,
            max_cluster_laterals_fraction,
            module_status_per_gene,
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
        for (module_index, anchor_index) in enumerate(anchor_index_per_module)
            genes_indices_of_anchor = get(genes_indices_of_anchor_index, anchor_index, nothing)
            if genes_indices_of_anchor !== nothing
                module_name = name_per_module[module_index]
                module_per_gene_per_block[genes_indices_of_anchor, block_index] .= module_name
                is_found_per_module_per_block[module_index, block_index] = true
            end
        end
    end

    set_matrix!(daf, "module", "block", "is_found", bestify(is_found_per_module_per_block); overwrite)
    set_matrix!(daf, "gene", "block", "module", module_per_gene_per_block; overwrite)
    set_matrix!(daf, "gene", "block", "module_status", module_status_per_gene_per_block; overwrite)

    return nothing
end

function compute_block_modules!(
    daf::DafWriter;
    block_name::AbstractString,
    mean_genes_per_cluster::Integer,
    anchor_correlation_genes::Integer,
    split_correlation_genes::Integer,
    min_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_correlation_genes_fraction::AbstractFloat,
    min_anchor_merge_correlation::AbstractFloat,
    min_anchor_merge_migration_fraction::AbstractFloat,
    min_orphan_cluster_correlation::AbstractFloat,
    min_gene_migration_correlation::AbstractFloat,
    max_anchor_kept_cluster_genes_baseline::Integer,
    min_anchor_kept_genes_fraction::AbstractFloat,
    min_anchor_kept_cluster_genes_likelihood::AbstractFloat,
    min_split_anchor_correlation::AbstractFloat,
    min_split_extra_correlation::AbstractFloat,
    max_cluster_laterals_fraction::AbstractFloat,
    module_status_per_gene::AbstractVector{<:AbstractString},
)::Dict{Int, Vector{Int}}
    n_full_genes = axis_length(daf, "gene")
    name_per_full_gene = axis_vector(daf, "gene")
    return adapter(  # NOJET
        daf;
        input_axes = [
            "block" => "=",
            "metacell" => "/ metacell & block => is_in_neighborhood ;= $(block_name)",
            "gene" => "/ gene & is_correlated_with_skeleton_in_neighborhood ; block = $(block_name)",
        ],
        input_data = [
            ("gene", "block", "is_neighborhood_varied") => "=",
            ("metacell", "gene", "log_linear_fraction") => "=",
            ("gene", "is_lateral") => "=",
            ("gene", "is_skeleton") => "=",
            ("gene", "block", "is_neighborhood_varied") => "=",
            ("gene", "full_index") => "index",
        ],
        output_axes = [],
        output_data = [],
    ) do adapted
        n_genes = axis_length(adapted, "gene")
        name_per_gene = axis_vector(adapted, "gene")
        full_index_per_gene = get_vector(adapted, "gene", "full_index").array
        log_fraction_per_metacell_per_gene = get_matrix(adapted, "metacell", "gene", "log_linear_fraction").array
        todox_max_log_fraction_per_per_gene = vec(maximum(log_fraction_per_metacell_per_gene; dims = 1))
        @assert all(todox_max_log_fraction_per_per_gene .>= log2(1e-5))
        is_skeleton_per_gene = get_vector(adapted, "gene", "is_skeleton").array
        is_lateral_per_gene = get_vector(adapted, "gene", "is_lateral").array
        is_varied_per_gene = adapted["/ gene / block = $(block_name) : is_neighborhood_varied"].array
        skeleton_indices = findall(is_skeleton_per_gene)
        @assert length(skeleton_indices) > 0

        correlations_between_genes = flame_timed("cor") do
            return cor(log_fraction_per_metacell_per_gene)  # NOLINT
        end
        @assert_matrix(correlations_between_genes, n_genes, n_genes, Columns)

        correlations_between_genes[isnan.(correlations_between_genes)] .= 0
        distances_between_genes = parallel_pairwise(Euclidean(), correlations_between_genes; dims = 2)
        @assert_matrix(distances_between_genes, n_genes, n_genes, Columns)
        genes_tree = hclust(distances_between_genes; linkage = :ward)

        # n_clusters = Int(ceil(n_genes / mean_genes_per_cluster))
        @debug "TODOX $(block_name) N_GENES: $(n_genes)"
        n_clusters = 20  # TODOX

        cluster_index_per_gene = cutree(genes_tree; k = n_clusters)
        @assert_vector(cluster_index_per_gene, n_genes)

        group_index_per_gene = zeros(Int, n_genes)
        anchor_index_per_group = Int[]
        most_correlated_genes_per_group = Vector{Int}[]
        mean_correlation_with_anchor_per_group = Float32[]
        migrated_genes_per_group = Maybe{Vector{Int}}[]

        for cluster_index in 1:n_clusters
            collect_cluster_groups!(;
                cluster_index,
                correlations_between_genes,
                skeleton_indices,
                cluster_index_per_gene,
                anchor_correlation_genes,
                split_correlation_genes,
                min_anchor_correlation_genes_fraction,
                max_anchor_correlation_genes_fraction,
                min_anchor_merge_correlation,
                min_anchor_merge_migration_fraction,
                max_anchor_kept_cluster_genes_baseline,
                min_split_anchor_correlation,
                min_split_extra_correlation,
                block_name,
                name_per_gene,
                group_index_per_gene,
                module_status_per_gene,
                full_index_per_gene,
                anchor_index_per_group,
                most_correlated_genes_per_group,
                mean_correlation_with_anchor_per_group,
                migrated_genes_per_group,
            )
        end

        groups_order = sortperm(mean_correlation_with_anchor_per_group; rev = true)
        n_groups = length(anchor_index_per_group)
        if min_anchor_merge_correlation != 1 && n_groups > 1
            for group_position in 1:(n_groups - 1)
                group_index = groups_order[group_position]
                other_group_indices = groups_order[(group_position + 1):n_groups]
                merge_group!(;
                    group_index,
                    skeleton_indices,
                    other_group_indices,
                    group_index_per_gene,
                    anchor_index_per_group,
                    anchor_correlation_genes,
                    min_anchor_correlation_genes_fraction,
                    max_anchor_correlation_genes_fraction,
                    min_anchor_merge_correlation,
                    min_anchor_merge_migration_fraction,
                    min_gene_migration_correlation,
                    max_anchor_kept_cluster_genes_baseline,
                    correlations_between_genes,
                    name_per_gene,
                    block_name,
                    most_correlated_genes_per_group,
                )
            end
        end

        anchor_indices = anchor_index_per_group[anchor_index_per_group .!= 0]
        n_anchors = length(anchor_indices)
        @assert n_anchors > 0

        for gene_index in 1:n_genes
            patch_group_index_per_gene!(;
                group_index_per_gene,
                module_status_per_gene,
                full_index_per_gene,
                anchor_index_per_group,
                anchor_indices,
                gene_index,
                correlations_between_genes,
                min_gene_migration_correlation,
                name_per_gene,
                migrated_genes_per_group,
                block_name,
            )
        end

        n_original_groups = n_groups

        if min_orphan_cluster_correlation < 1
            is_orphan_per_gene = (group_index_per_gene .== 0) .& is_varied_per_gene
            indices_of_orphans = findall(is_orphan_per_gene)
            if length(indices_of_orphans) > 0
                correlations_between_orphans = correlations_between_genes[indices_of_orphans, indices_of_orphans]
                distances_between_orphans = clamp.(1.0 .- correlations_between_orphans, 0.0, 1.0)
                orphans_tree = hclust(distances_between_orphans; linkage = :complete)
                cluster_index_per_orphan = cutree(orphans_tree; h = 1.0 - min_orphan_cluster_correlation)
                n_clusters = maximum(cluster_index_per_orphan)
                for cluster_index in 1:n_clusters
                    indices_of_cluster_genes = indices_of_orphans[cluster_index_per_orphan .== cluster_index]
                    n_cluster_genes = length(indices_of_cluster_genes)
                    @assert n_cluster_genes > 0
                    n_cluster_lateral_genes = sum(is_lateral_per_gene[indices_of_cluster_genes])
                    if n_cluster_lateral_genes <= n_cluster_genes * max_cluster_laterals_fraction
                        anchor_index, _ = find_anchor(;
                            skeleton_indices = indices_of_cluster_genes,
                            gene_indices_of_cluster = indices_of_cluster_genes,
                            correlations_between_genes,
                            anchor_correlation_genes,
                            min_anchor_correlation_genes_fraction,
                            max_anchor_correlation_genes_fraction,
                        )
                        push!(anchor_indices, anchor_index)
                        n_anchors += 1
                        n_groups += 1
                        group_index_per_gene[indices_of_cluster_genes] .= n_groups
                        module_status_per_gene[full_index_per_gene[indices_of_cluster_genes]] .= "found(orphan)"
                        push!(anchor_index_per_group, anchor_index)
                        push!(most_correlated_genes_per_group, indices_of_cluster_genes)
                        push!(migrated_genes_per_group, nothing)
                        @debug "TODOX Block: $(block_name) Group: $(n_groups) Orphan anchor: $(name_per_gene[anchor_index]) Genes: $(length(indices_of_cluster_genes))"
                    else
                        module_status_per_gene[full_index_per_gene[indices_of_cluster_genes]] .= "lateral(orphan)"
                    end
                end
            end

            for gene_index in 1:n_genes
                patch_group_index_per_gene!(;
                    group_index_per_gene,
                    module_status_per_gene,
                    full_index_per_gene,
                    anchor_index_per_group,
                    anchor_indices,
                    gene_index,
                    correlations_between_genes,
                    min_gene_migration_correlation,
                    name_per_gene,
                    migrated_genes_per_group,
                    block_name,
                )
            end
        end

        min_kept_genes_fraction =
            max(min_anchor_kept_cluster_genes_likelihood / n_anchors, min_anchor_kept_genes_fraction)

        save_module_status_per_gene = copy_array(module_status_per_gene)
        genes_full_indices_of_anchor_full_index = Dict{Int, Vector{Int}}()
        final_anchor_indices = UInt32[]
        while true
            for group_index in 1:n_groups
                @debug "TODOX GROUP: $(group_index)"
                anchor_index = anchor_index_per_group[group_index]
                if anchor_index == 0
                    @debug "TODOX Block: $(block_name) Group: $(group_index) was merged"
                    continue
                end
                most_correlated_genes_of_anchor = most_correlated_genes_per_group[group_index]
                n_most_correlated_genes = length(most_correlated_genes_of_anchor)
                n_kept_most_correlated_genes =
                    sum(group_index_per_gene[most_correlated_genes_of_anchor] .== group_index)
                indices_of_group_genes = findall(group_index_per_gene .== group_index)

                if n_kept_most_correlated_genes >= n_most_correlated_genes * min_kept_genes_fraction
                    anchor_full_index = full_index_per_gene[anchor_index]
                    n_group_genes = length(indices_of_group_genes)
                    is_lateral_per_group_gene = is_lateral_per_gene[indices_of_group_genes]
                    n_lateral_group_genes = sum(is_lateral_per_group_gene)
                    if n_lateral_group_genes <= n_group_genes * max_cluster_laterals_fraction
                        indices_of_lateral_group_genes = indices_of_group_genes[is_lateral_per_group_gene]
                        module_status_per_gene[full_index_per_gene[indices_of_lateral_group_genes]] .=
                            "lateral;" .* module_status_per_gene[full_index_per_gene[indices_of_lateral_group_genes]]
                        indices_of_final_group_genes = indices_of_group_genes[.!is_lateral_per_group_gene]
                        genes_full_indices = full_index_per_gene[indices_of_final_group_genes]
                        genes_full_indices_of_anchor_full_index[anchor_full_index] = genes_full_indices
                        push!(final_anchor_indices, anchor_index)
                        @debug "TODOX Block: $(block_name) Group: $(group_index) n_kept_most_correlated_genes: $(n_kept_most_correlated_genes) >= n_most_correlated_genes: $(n_most_correlated_genes) * min_kept_genes_fraction: $(min_kept_genes_fraction) = $(n_most_correlated_genes * min_kept_genes_fraction), n_lateral_group_genes: $(n_lateral_group_genes) <= n_group_genes: $(n_group_genes) * max_cluster_laterals_fraction: $(max_cluster_laterals_fraction) = $(n_group_genes * max_cluster_laterals_fraction)"
                    else
                        @debug "TODOX Block: $(block_name) Group: $(group_index) n_kept_most_correlated_genes: $(n_kept_most_correlated_genes) >= n_most_correlated_genes: $(n_most_correlated_genes) * min_kept_genes_fraction: $(min_kept_genes_fraction) = $(n_most_correlated_genes * min_kept_genes_fraction), n_lateral_group_genes: $(n_lateral_group_genes) > n_group_genes: $(n_group_genes) * max_cluster_laterals_fraction: $(max_cluster_laterals_fraction) = $(n_group_genes * max_cluster_laterals_fraction)"
                        module_status_per_gene[full_index_per_gene[indices_of_group_genes]] .*= ";lateral"
                    end
                else
                    @debug "TODOX Block: $(block_name) Group: $(group_index) n_kept_most_correlated_genes: $(n_kept_most_correlated_genes) < n_most_correlated_genes: $(n_most_correlated_genes) * min_kept_genes_fraction: $(min_kept_genes_fraction) = $(n_most_correlated_genes * min_kept_genes_fraction)"
                    module_status_per_gene[full_index_per_gene[indices_of_group_genes]] .*= ";vacated"
                end
            end

            n_final_anchors = length(final_anchor_indices)
            @debug "- Block: $(block_name) Anchors: $(n_final_anchors) [ $(join(name_per_gene[final_anchor_indices], ", ")) ]"
            if n_final_anchors > 0
                break
            end

            @debug "TODOX NO ANCHORS IN BLOCK $(block_name) WITH $(max_cluster_laterals_fraction)"
            max_cluster_laterals_fraction *= 1.1
            module_status_per_gene .= save_module_status_per_gene
            empty!(genes_full_indices_of_anchor_full_index)
            empty!(final_anchor_indices)
        end

        todox_is_anchored = zeros(Bool, n_full_genes)
        for (_, gene_full_indices) in genes_full_indices_of_anchor_full_index
            for gene_full_index in gene_full_indices
                @assert !todox_is_anchored[gene_full_index]
                todox_is_anchored[gene_full_index] = true
                @assert module_status_per_gene[gene_full_index] != ""
            end
        end
        for gene_full_index in 1:n_full_genes
            if todox_is_anchored[gene_full_index]
                if module_status_per_gene[gene_full_index] == "" ||
                   contains(module_status_per_gene[gene_full_index], ";")
                    @debug "TODOX ANCHORED GENE $(name_per_full_gene[gene_full_index]) STATUS: $(module_status_per_gene[gene_full_index])"
                    @assert false
                end
            else
                if module_status_per_gene[gene_full_index] != "" &&
                   !contains(module_status_per_gene[gene_full_index], ";")
                    @debug "TODOX UNANCHORED GENE $(name_per_full_gene[gene_full_index]) STATUS: $(module_status_per_gene[gene_full_index])"
                    @assert false
                end
            end
        end

        return genes_full_indices_of_anchor_full_index
    end
end

function collect_cluster_groups!(;
    cluster_index::Integer,
    correlations_between_genes::AbstractMatrix{<:AbstractFloat},
    skeleton_indices::AbstractVector{<:Integer},
    cluster_index_per_gene::AbstractVector{<:Integer},
    anchor_correlation_genes::Integer,
    split_correlation_genes::Integer,
    min_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_correlation_genes_fraction::AbstractFloat,
    min_anchor_merge_correlation::AbstractFloat,
    min_anchor_merge_migration_fraction::AbstractFloat,
    max_anchor_kept_cluster_genes_baseline::Integer,
    min_split_anchor_correlation::AbstractFloat,
    min_split_extra_correlation::AbstractFloat,
    block_name::AbstractString,
    name_per_gene::AbstractVector{<:AbstractString},
    group_index_per_gene::Vector{Int},
    module_status_per_gene::AbstractVector{<:AbstractString},
    full_index_per_gene::AbstractVector{<:Integer},
    anchor_index_per_group::Vector{Int},
    most_correlated_genes_per_group::Vector{Vector{Int}},
    migrated_genes_per_group::Vector{Maybe{Vector{Int}}},
    mean_correlation_with_anchor_per_group::AbstractVector{<:AbstractFloat},
)::Nothing
    is_cluster_per_gene = cluster_index_per_gene .== cluster_index
    gene_indices_of_cluster = findall(is_cluster_per_gene)
    n_cluster_genes = length(gene_indices_of_cluster)
    @assert n_cluster_genes > 0

    anchor_index, anchor_correlation = find_anchor(;
        skeleton_indices,
        gene_indices_of_cluster,
        correlations_between_genes,
        anchor_correlation_genes,
        min_anchor_correlation_genes_fraction,
        max_anchor_correlation_genes_fraction,
    )

    if anchor_index == 0
        @debug "TODOX Block: $(block_name) Cluster: $(cluster_index) no anchors"
        return
    end

    split_correlation_quantile = clamp(
        split_correlation_genes / n_cluster_genes,
        min_anchor_correlation_genes_fraction,
        max_anchor_correlation_genes_fraction,
    )

    correlation_with_anchor_per_cluster_gene = correlations_between_genes[gene_indices_of_cluster, anchor_index]

    best_split_anchor_index = 0
    best_split_anchor_split_correlation = -2
    best_anchor_split_correlation = -2

    for skeleton_index in skeleton_indices
        if is_cluster_per_gene[skeleton_index] && skeleton_index != anchor_index
            correlation_with_skeleton_per_cluster_gene =
                correlations_between_genes[gene_indices_of_cluster, skeleton_index]
            skeleton_split_correlation =
                quantile(correlation_with_skeleton_per_cluster_gene, 1 - split_correlation_quantile)  # NOLINTA
            correlated_cluster_gene_positions =
                findall(correlation_with_skeleton_per_cluster_gene .>= skeleton_split_correlation)
            skeleton_split_correlation =
                minimum(correlation_with_skeleton_per_cluster_gene[correlated_cluster_gene_positions])
            anchor_split_correlation =
                minimum(correlation_with_anchor_per_cluster_gene[correlated_cluster_gene_positions])
            if skeleton_split_correlation - anchor_split_correlation >
               best_split_anchor_split_correlation - best_anchor_split_correlation
                best_split_anchor_index = skeleton_index
                best_split_anchor_split_correlation = skeleton_split_correlation
                best_anchor_split_correlation = anchor_split_correlation
            end
        end
    end

    if best_split_anchor_index != 0 &&
       best_split_anchor_split_correlation >= min_split_anchor_correlation &&
       best_split_anchor_split_correlation - best_anchor_split_correlation >= min_split_extra_correlation
        @debug "TODOX Block: $(block_name) Group: $(length(anchor_index_per_group)) Cluster: $(cluster_index) anchor: $(name_per_gene[anchor_index]) split correlation: $(best_anchor_split_correlation) split: $(name_per_gene[best_split_anchor_index]) correlation: $(best_split_anchor_split_correlation)"
        correlation_with_split_per_cluster_gene =
            correlations_between_genes[gene_indices_of_cluster, best_split_anchor_index]
        is_split_per_cluster_gene = correlation_with_split_per_cluster_gene .> correlation_with_anchor_per_cluster_gene

        gene_indices_of_anchor = gene_indices_of_cluster[.!is_split_per_cluster_gene]
        gene_indices_of_split = gene_indices_of_cluster[is_split_per_cluster_gene]
        mean_correlation_with_anchor = mean(correlation_with_anchor_per_cluster_gene[.!is_split_per_cluster_gene])
        mean_correlation_with_split = mean(correlation_with_anchor_per_cluster_gene[is_split_per_cluster_gene])
        groups = (
            (mean_correlation_with_anchor, anchor_index, gene_indices_of_anchor),
            (mean_correlation_with_split, best_split_anchor_index, gene_indices_of_split),
        )
    else
        mean_correlation_with_anchor = mean(correlation_with_anchor_per_cluster_gene)
        groups = ((mean_correlation_with_anchor, anchor_index, gene_indices_of_cluster),)
    end

    for (mean_correlation_with_anchor, selected_anchor_index, gene_indices_of_anchor) in groups
        most_correlated_genes_of_anchor = most_correlated_genes_of_cluster(;
            anchor_index = selected_anchor_index,
            gene_indices_of_anchor,
            max_anchor_kept_cluster_genes_baseline,
            correlations_between_genes,
        )

        push!(anchor_index_per_group, selected_anchor_index)
        push!(most_correlated_genes_per_group, most_correlated_genes_of_anchor)
        push!(mean_correlation_with_anchor_per_group, mean_correlation_with_anchor)
        push!(migrated_genes_per_group, nothing)
        group_index_per_gene[gene_indices_of_anchor] .= length(anchor_index_per_group)
        module_status_per_gene[full_index_per_gene[gene_indices_of_anchor]] .= "found(cluster)"

        @debug "TODOX Block: $(block_name) Group: $(length(anchor_index_per_group)) Cluster: $(cluster_index) anchor: $(name_per_gene[selected_anchor_index]) mean_correlation_with_anchor: $(mean_correlation_with_anchor) most_correlated_genes: $(length(most_correlated_genes_of_anchor))"
    end

    return nothing
end

function find_anchor(;
    skeleton_indices::AbstractVector{<:Integer},
    gene_indices_of_cluster::AbstractVector{<:Integer},
    correlations_between_genes::AbstractMatrix{<:AbstractFloat},
    anchor_correlation_genes::Integer,
    min_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_correlation_genes_fraction::AbstractFloat,
)::Tuple{Integer, AbstractFloat}
    n_cluster_genes = length(gene_indices_of_cluster)

    anchor_correlation_quantile = clamp(
        anchor_correlation_genes / n_cluster_genes,
        min_anchor_correlation_genes_fraction,
        max_anchor_correlation_genes_fraction,
    )

    anchor_index = 0
    anchor_correlation = -2.0

    for skeleton_index in skeleton_indices
        if in(skeleton_index, gene_indices_of_cluster)
            correlation_with_skeleton_per_cluster_gene =
                vec(correlations_between_genes[gene_indices_of_cluster, skeleton_index])
            skeleton_correlation = quantile(correlation_with_skeleton_per_cluster_gene, 1 - anchor_correlation_quantile)  # NOLINT
            if skeleton_correlation > anchor_correlation
                anchor_correlation = skeleton_correlation
                anchor_index = skeleton_index
            end
        end
    end

    return (anchor_index, anchor_correlation)
end

function merge_group!(;
    group_index::Integer,
    skeleton_indices::AbstractVector{<:Integer},
    other_group_indices::AbstractVector{<:Integer},
    group_index_per_gene::AbstractVector{<:Integer},
    anchor_index_per_group::AbstractVector{<:Integer},
    anchor_correlation_genes::Integer,
    min_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_correlation_genes_fraction::AbstractFloat,
    min_anchor_merge_correlation::AbstractFloat,
    min_anchor_merge_migration_fraction::AbstractFloat,
    min_gene_migration_correlation::AbstractFloat,
    max_anchor_kept_cluster_genes_baseline::Integer,
    correlations_between_genes::AbstractMatrix{<:AbstractFloat},
    name_per_gene::AbstractVector{<:AbstractString},  # NOLINT TODOX
    block_name::AbstractString,
    most_correlated_genes_per_group::Vector{Vector{Int}},
)::Nothing
    anchor_index = anchor_index_per_group[group_index]
    if anchor_index == 0
        return nothing
    end

    while true
        most_correlated_genes_of_group = most_correlated_genes_per_group[group_index]
        best_other_group_index = 0
        best_other_group_score = 0.0
        best_other_anchor_index = 0.0

        for other_group_index in other_group_indices
            other_anchor_index = anchor_index_per_group[other_group_index]
            if other_anchor_index != 0
                @assert other_anchor_index != anchor_index
                correlation_between_anchors = correlations_between_genes[anchor_index, other_anchor_index]
                if correlation_between_anchors >= min_anchor_merge_correlation
                    most_correlated_genes_of_other = most_correlated_genes_per_group[other_group_index]
                    @views correlation_of_group =
                        correlations_between_genes[other_anchor_index, most_correlated_genes_of_group]
                    @views correlation_of_others =
                        correlations_between_genes[anchor_index, most_correlated_genes_of_other]
                    fraction_of_correlated_group =
                        sum(correlation_of_group .>= min_gene_migration_correlation) / length(correlation_of_group)
                    fraction_of_correlated_others =
                        sum(correlation_of_others .>= min_gene_migration_correlation) / length(correlation_of_others)
                    if fraction_of_correlated_group >= min_anchor_merge_migration_fraction &&
                       fraction_of_correlated_others >= min_anchor_merge_migration_fraction
                        other_group_score =
                            correlation_between_anchors / (min_anchor_merge_correlation + 1e-6) +
                            min(fraction_of_correlated_group, fraction_of_correlated_others) /
                            (min_anchor_merge_migration_fraction + 1e-6)
                        @debug "TODOX Block: $(block_name) group: $(group_index) anchor: $(name_per_gene[anchor_index]) other group: $(other_group_index) anchor: $(name_per_gene[other_anchor_index]) correlation_between_anchors: $(correlation_between_anchors) fraction_of_correlated_genes: $(fraction_of_correlated_group) $(fraction_of_correlated_others) score: $(other_group_score)"
                        if other_group_score > best_other_group_score
                            best_other_group_index = other_group_index
                            best_other_group_score = other_group_score
                            best_other_anchor_index = other_anchor_index
                        end
                    else
                        @debug "TODOX Block: $(block_name) group: $(group_index) anchor: $(name_per_gene[anchor_index]) other group: $(other_group_index) anchor: $(name_per_gene[other_anchor_index]) correlation_between_anchors: $(correlation_between_anchors) LOW fraction_of_correlated_genes: $(fraction_of_correlated_group) $(fraction_of_correlated_others)"
                    end
                else
                    @debug "TODOX Block: $(block_name) group: $(group_index) anchor: $(name_per_gene[anchor_index]) other group: $(other_group_index) anchor: $(name_per_gene[other_anchor_index]) LOW correlation_between_anchors: $(correlation_between_anchors)"
                end
            end
        end

        if best_other_group_index == 0
            @assert best_other_anchor_index == 0
            return nothing
        else
            @assert best_other_anchor_index != 0
        end

        anchor_index_per_group[best_other_group_index] = 0
        empty!(most_correlated_genes_per_group[best_other_group_index])

        indices_of_genes_of_merged =
            findall((group_index_per_gene .== group_index) .| (group_index_per_gene .== best_other_group_index))
        group_index_per_gene[indices_of_genes_of_merged] .= group_index

        anchor_index_per_group[group_index], anchor_correlation = find_anchor(;
            skeleton_indices,
            gene_indices_of_cluster = indices_of_genes_of_merged,
            correlations_between_genes,
            anchor_correlation_genes,
            min_anchor_correlation_genes_fraction,
            max_anchor_correlation_genes_fraction,
        )

        @assert anchor_index_per_group[group_index] != 0
        @assert anchor_index != 0  # TODOX
        @debug "TODOX Block: $(block_name) merge group: $(group_index) anchor: $(name_per_gene[anchor_index]) new: $(name_per_gene[anchor_index_per_group[group_index]]) correlation: $(anchor_correlation) other group: $(best_other_group_index) anchor: $(name_per_gene[best_other_anchor_index]) score: $(best_other_group_score)"

        most_correlated_genes_per_group[group_index] = most_correlated_genes_of_cluster(;
            anchor_index,
            gene_indices_of_anchor = indices_of_genes_of_merged,
            max_anchor_kept_cluster_genes_baseline,
            correlations_between_genes,
        )
    end
    @assert false
end

function most_correlated_genes_of_cluster(;
    anchor_index::Integer,
    gene_indices_of_anchor::AbstractVector{<:Integer},
    max_anchor_kept_cluster_genes_baseline::Integer,
    correlations_between_genes::AbstractMatrix{<:AbstractFloat},
)::AbstractVector{<:Integer}
    n_cluster_genes = length(gene_indices_of_anchor)
    if n_cluster_genes <= max_anchor_kept_cluster_genes_baseline
        return gene_indices_of_anchor
    else
        correlations_of_cluster_genes_with_anchor =
            vec(correlations_between_genes[gene_indices_of_anchor, anchor_index])
        threshold = quantile(  # NOLINT
            correlations_of_cluster_genes_with_anchor,
            1 - max_anchor_kept_cluster_genes_baseline / n_cluster_genes,
        )
        return gene_indices_of_anchor[correlations_of_cluster_genes_with_anchor .>= threshold]
    end
end

function patch_group_index_per_gene!(;
    group_index_per_gene::AbstractVector{<:Integer},
    module_status_per_gene::AbstractVector{<:AbstractString},
    full_index_per_gene::AbstractVector{<:Integer},
    anchor_index_per_group::AbstractVector{<:Integer},
    anchor_indices::AbstractVector{<:Integer},
    gene_index::Integer,
    correlations_between_genes::AbstractMatrix{<:AbstractFloat},
    min_gene_migration_correlation::AbstractFloat,
    name_per_gene::AbstractVector{<:AbstractString},  # NOLINT TODOX
    migrated_genes_per_group::Vector{Maybe{Vector{Int}}},
    block_name::AbstractString,
)::Nothing
    current_group = group_index_per_gene[gene_index]
    anchors_correlations_with_gene = correlations_between_genes[anchor_indices, gene_index]
    most_correlated_anchor_index = anchor_indices[argmax(anchors_correlations_with_gene)]
    group_of_most_correlated_anchor = group_index_per_gene[most_correlated_anchor_index]

    if group_of_most_correlated_anchor != current_group
        if current_group < 0
            current_anchor_index = anchor_index_per_group[current_group]
        else
            current_anchor_index = 0
        end
        if current_anchor_index > 0
            current_anchor_name = name_per_gene[current_anchor_index]
        else
            current_anchor_name = "None"
        end
        most_correlated_anchor_name = name_per_gene[most_correlated_anchor_index]
        if correlations_between_genes[most_correlated_anchor_index, gene_index] < min_gene_migration_correlation
            # @debug "TODOX Block: $(block_name) gene: $(name_per_gene[gene_index]) group: $(current_group) != most_correlated_anchor_group: $(group_of_most_correlated_anchor) of anchor: $(name_per_gene[most_correlated_anchor_index]) correlation $(correlations_between_genes[most_correlated_anchor_index, gene_index]) < min_gene_migration_correlation: $(min_gene_migration_correlation)"
            group_index_per_gene[gene_index] = 0
            module_status_per_gene[full_index_per_gene[gene_index]] = "confused($(current_anchor_name);$(most_correlated_anchor_name))"
        else
            # @debug "TODOX Block: $(block_name) gene: $(name_per_gene[gene_index]) group: $(current_group) != most_correlated_anchor_group: $(group_of_most_correlated_anchor) of anchor: $(name_per_gene[most_correlated_anchor_index]) correlation $(correlations_between_genes[most_correlated_anchor_index, gene_index]) >= min_gene_migration_correlation: $(min_gene_migration_correlation)"
            group_index_per_gene[gene_index] = group_of_most_correlated_anchor
            module_status_per_gene[full_index_per_gene[gene_index]] = "migrated($(current_anchor_name),$(most_correlated_anchor_name))"
            migrated_genes_of_group_of_most_correlated_anchor =
                migrated_genes_per_group[group_of_most_correlated_anchor]
            if migrated_genes_of_group_of_most_correlated_anchor === nothing
                migrated_genes_per_group[group_of_most_correlated_anchor] = [gene_index]
            else
                push!(migrated_genes_of_group_of_most_correlated_anchor, gene_index)
            end
        end
    end

    return nothing
end

@logged @computation Contract(
    axes = [
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(GuaranteedOutput),
    ],
    data = [
        gene_is_lateral_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_is_neighborhood_marker_matrix(OptionalInput),  # TODOX
        block_gene_is_neighborhood_distinct_matrix(OptionalInput),  # TODOX
        block_gene_is_correlated_with_skeleton_in_neighborhood_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        module_anchor_gene_vector(GuaranteedOutput),
        block_module_is_found_matrix(GuaranteedOutput),
        block_gene_module_matrix(GuaranteedOutput),
        block_gene_module_status_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_k!(
    daf::DafWriter;
    max_clusters::Integer = 24,
    max_cluster_laterals_fraction::Maybe{AbstractFloat} = 0.4,
    min_member_correlation::AbstractFloat = 0.5,
    min_orphan_correlation::AbstractFloat = 0.6,
    min_strong_UMIs::Integer = 8,
    min_strong_cells::Integer = 12,
    kmeans_rounds::Integer = function_default(kmeans_in_rounds, :rounds),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert max_cluster_laterals_fraction === nothing || 0 <= max_cluster_laterals_fraction <= 1
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
    block_index_per_metacell = daf["/ metacell : block => index"].array
    block_index_per_cell = daf["/ cell : metacell ?? 0 => block => index"].array
    log_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "log_linear_fraction").array

    genes_indices_of_anchor_index_per_block = Vector{Dict{Int, Vector{Int}}}(undef, n_blocks)

    module_status_per_gene_per_block = fill("", n_genes, n_blocks)

    parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng
        block_name = name_per_block[block_index]
        @views module_status_per_gene = module_status_per_gene_per_block[:, block_index]
        return genes_indices_of_anchor_index_per_block[block_index] = compute_block_modules_k!(
            daf;
            block_index,
            block_name,
            max_clusters,
            max_cluster_laterals_fraction,
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
        for (module_index, anchor_index) in enumerate(anchor_index_per_module)
            genes_indices_of_anchor = get(genes_indices_of_anchor_index, anchor_index, nothing)
            if genes_indices_of_anchor !== nothing
                module_name = name_per_module[module_index]
                module_per_gene_per_block[genes_indices_of_anchor, block_index] .= module_name
                is_found_per_module_per_block[module_index, block_index] = true
            end
        end
    end

    set_matrix!(daf, "module", "block", "is_found", bestify(is_found_per_module_per_block); overwrite)
    set_matrix!(daf, "gene", "block", "module", module_per_gene_per_block; overwrite)
    set_matrix!(daf, "gene", "block", "module_status", module_status_per_gene_per_block; overwrite)

    return nothing
end

function compute_block_modules_k!(
    daf::DafWriter;
    block_index::Integer,
    block_name::AbstractString,
    max_clusters::Integer,
    max_cluster_laterals_fraction::Maybe{AbstractFloat},
    min_member_correlation::AbstractFloat,
    min_orphan_correlation::AbstractFloat,
    min_strong_UMIs::Integer,
    min_strong_cells::Integer,
    kmeans_rounds::Integer,
    rng::AbstractRNG,
    module_status_per_gene::AbstractVector{<:AbstractString},
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

    local center_per_neighborhood_metacell_per_lateral_cluster
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
            vec(mean(log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene; dims = 1))
        @assert_vector(mean_log_fraction_per_lateral_neighborhood_marker_gene, n_lateral_neighborhood_marker_genes)

        std_log_fraction_per_lateral_neighborhood_marker_gene = vec(
            std(
                log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene;
                mean = transpose(mean_log_fraction_per_lateral_neighborhood_marker_gene),
                dims = 1,
            ),
        )
        @assert_vector(std_log_fraction_per_lateral_neighborhood_marker_gene, n_lateral_neighborhood_marker_genes)

        normalized_log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene =
            (
                log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene .-
                transpose(mean_log_fraction_per_lateral_neighborhood_marker_gene)
            ) ./ transpose(std_log_fraction_per_lateral_neighborhood_marker_gene)

        kmeans_results = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                normalized_log_fraction_per_neighborhood_metacell_per_lateral_neighborhood_marker_gene,
                max_clusters;
                rounds = kmeans_rounds,
                rng,
            )
        end

        lateral_cluster_per_gene[indices_of_lateral_neighborhood_markers] .= kmeans_results.assignments

        n_lateral_clusters = size(kmeans_results.centers, 2)
        for lateral_cluster_index in 1:n_lateral_clusters
            @debug "Block: $(block_name) lateral_cluster: $(lateral_cluster_index) genes: [ \"$(join(name_per_gene[indices_of_lateral_neighborhood_markers[kmeans_results.assignments .== lateral_cluster_index]], "\", \""))\" ]"
        end

        center_per_neighborhood_metacell_per_lateral_cluster = kmeans_results.centers
        return nothing
    end

    local center_per_neighborhood_metacell_per_cluster
    local cluster_index_per_local_gene
    local indices_of_local_genes
    local indices_of_neighborhood_cells
    local is_skeleton_per_local_gene
    local is_skeleton_per_local_gene
    local n_local_genes
    local n_neighborhood_cells
    local normalized_log_fraction_per_neighborhood_metacell_per_local_gene

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

        is_skeleton_per_local_gene = is_skeleton_per_gene[indices_of_local_genes]

        log_fraction_per_neighborhood_metacell_per_local_gene =
            log_fraction_per_metacell_per_gene[indices_of_neighborhood_metacells, indices_of_local_genes]

        mean_log_fraction_per_local_gene = vec(mean(log_fraction_per_neighborhood_metacell_per_local_gene; dims = 1))
        @assert_vector(mean_log_fraction_per_local_gene, n_local_genes)

        std_log_fraction_per_local_gene = vec(
            std(
                log_fraction_per_neighborhood_metacell_per_local_gene;
                mean = transpose(mean_log_fraction_per_local_gene),
                dims = 1,
            ),
        )
        @assert_vector(std_log_fraction_per_local_gene, n_local_genes)

        normalized_log_fraction_per_neighborhood_metacell_per_local_gene =
            (log_fraction_per_neighborhood_metacell_per_local_gene .- transpose(mean_log_fraction_per_local_gene)) ./
            transpose(std_log_fraction_per_local_gene)

        kmeans_results = flame_timed("kmeans_in_rounds") do
            return kmeans_in_rounds(
                normalized_log_fraction_per_neighborhood_metacell_per_local_gene,
                max_clusters;
                rounds = kmeans_rounds,
                rng,
            )
        end

        cluster_index_per_local_gene = kmeans_results.assignments
        @assert_vector(cluster_index_per_local_gene, n_local_genes)

        center_per_neighborhood_metacell_per_cluster = kmeans_results.centers
        return module_status_per_gene[indices_of_local_genes] .=
            "cluster(" .* string.(cluster_index_per_local_gene) .* ")"
    end

    local is_lateral_ish_per_local_gene

    flame_timed("identify_lateral_ish_genes") do
        is_lateral_ish_per_local_gene = zeros(Bool, n_local_genes)

        is_neighborhood_distinct_per_local_gene =
            is_neighborhood_marker_per_gene_per_block[indices_of_local_genes, block_index]
        @debug "TODOX Block: $(block_name) distinct genes: $(sum(is_neighborhood_distinct_per_local_gene)) names: [ \"$(join(name_per_gene[indices_of_local_genes[is_neighborhood_distinct_per_local_gene]] , "\", \""))\" ]"

        for local_gene_position in 1:n_local_genes
            @views normalized_log_fraction_per_neighborhood_metacell =
                normalized_log_fraction_per_neighborhood_metacell_per_local_gene[:, local_gene_position]
            cluster_index = cluster_index_per_local_gene[local_gene_position]
            @views cluster_log_fraction_per_neighborhood_metacell =
                center_per_neighborhood_metacell_per_cluster[:, cluster_index]
            cluster_correlation = flame_timed("cor") do
                return cor(
                    normalized_log_fraction_per_neighborhood_metacell,
                    cluster_log_fraction_per_neighborhood_metacell,
                )
            end
            correlation_per_lateral_cluster = flame_timed("cor") do
                return vec(
                    cor(
                        normalized_log_fraction_per_neighborhood_metacell,
                        center_per_neighborhood_metacell_per_lateral_cluster,
                    ),
                )
            end
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
                @debug "Block $(block_name)$(qualifier) lateral-ish gene: $(name_per_gene[indices_of_local_genes[local_gene_position]]) lateral_correlation: $(lateral_correlation) > cluster_correlation: $(cluster_correlation) lateral cluster: $(lateral_cluster_index)"
            else
                @debug "Block $(block_name)$(qualifier) kept gene: $(name_per_gene[indices_of_local_genes[local_gene_position]]) lateral_correlation: $(lateral_correlation) <= cluster_correlation: $(cluster_correlation) lateral cluster: $(lateral_cluster_index)"
            end
        end

        return module_status_per_gene[indices_of_local_genes[is_lateral_ish_per_local_gene]] .*= ",lateral_ish"
    end

    local anchor_local_position_per_cluster
    local cluster_index_of_anchor_local_position

    flame_timed("correlated_pertinent_gene_modules") do
        n_clusters = size(center_per_neighborhood_metacell_per_cluster, 2)
        @debug "TODOX Block: $(block_name) raw clusters: $(n_clusters)"
        @assert 0 < n_clusters <= max_clusters

        anchor_local_position_per_cluster = fill(UInt32(0), n_clusters)
        cluster_index_of_anchor_local_position = Dict{UInt32, UInt32}()

        for cluster_index in 1:n_clusters
            @views cluster_center_per_neighborhood_metacell =
                center_per_neighborhood_metacell_per_cluster[:, cluster_index]
            is_in_cluster_per_local_gene = cluster_index_per_local_gene .== cluster_index
            n_cluster_genes = sum(is_in_cluster_per_local_gene)
            @assert n_cluster_genes > 0
            is_lateral_ish_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& is_lateral_ish_per_local_gene
            n_lateral_cluster_genes = sum(is_lateral_ish_in_cluster_per_local_gene)

            cluster_index_per_local_gene[is_lateral_ish_in_cluster_per_local_gene] .= -1
            is_pertinent_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& .!is_lateral_ish_per_local_gene

            is_skeleton_in_cluster_per_local_gene = is_in_cluster_per_local_gene .& is_skeleton_per_local_gene
            local_positions_of_skeletons_in_cluster = findall(is_skeleton_in_cluster_per_local_gene)
            n_cluster_skeletons = length(local_positions_of_skeletons_in_cluster)

            if n_cluster_skeletons == 0
                cluster_index_per_local_gene[is_pertinent_in_cluster_per_local_gene] .= 0
                module_status_per_gene[indices_of_local_genes[is_pertinent_in_cluster_per_local_gene]] .*= ",orphan_cluster"

                @debug "TODOX Block: $(block_name) orphan cluster: $(cluster_index) genes: $(n_cluster_genes) lateral: $(n_lateral_cluster_genes) ($(percent(n_lateral_cluster_genes, n_cluster_genes)))"
                continue
            end

            normalized_log_fraction_per_neighborhood_metacell_per_pertinent_local_gene =
                normalized_log_fraction_per_neighborhood_metacell_per_local_gene[
                    :,
                    is_pertinent_in_cluster_per_local_gene,
                ]
            correlation_between_center_and_pertinent_local_genes = flame_timed("cor") do
                return vec(
                    cor(
                        cluster_center_per_neighborhood_metacell,
                        normalized_log_fraction_per_neighborhood_metacell_per_pertinent_local_gene,
                    ),
                )
            end
            @assert_vector(
                correlation_between_center_and_pertinent_local_genes,
                sum(is_pertinent_in_cluster_per_local_gene)
            )

            is_correlated_pertinent_in_cluster_per_local_gene = zeros(Bool, n_local_genes)
            is_correlated_pertinent_in_cluster_per_local_gene[is_pertinent_in_cluster_per_local_gene] .=
                correlation_between_center_and_pertinent_local_genes .>= min_member_correlation

            is_uncorrelated_pertinent_in_cluster_per_local_gene =
                is_pertinent_in_cluster_per_local_gene .& .!is_correlated_pertinent_in_cluster_per_local_gene
            cluster_index_per_local_gene[is_uncorrelated_pertinent_in_cluster_per_local_gene] .= 0

            n_correlated_pertinent_local_genes = sum(is_correlated_pertinent_in_cluster_per_local_gene)
            n_uncorrelated_pertinent_cluster_genes = sum(is_uncorrelated_pertinent_in_cluster_per_local_gene)

            @assert n_correlated_pertinent_local_genes == sum(cluster_index_per_local_gene .== cluster_index)
            if n_correlated_pertinent_local_genes == 0
                module_status_per_gene[indices_of_local_genes[is_uncorrelated_pertinent_in_cluster_per_local_gene]] .*= ",uncorrelated_gene,uncorrelated_cluster"
                if n_uncorrelated_pertinent_cluster_genes == 0
                    reason = "all-lateral"
                else
                    reason = "all-uncorrelated"
                end
                @debug "TODOX Block: $(block_name) $(reason) cluster: $(cluster_index) genes: $(n_cluster_genes) lateral: $(n_lateral_cluster_genes) ($(percent(n_lateral_cluster_genes, n_cluster_genes))) uncorrelated pertinent: $(n_uncorrelated_pertinent_cluster_genes) ($(percent(n_uncorrelated_pertinent_cluster_genes, n_cluster_genes))) correlated pertinent: $(n_correlated_pertinent_local_genes) ($(percent(n_correlated_pertinent_local_genes, n_cluster_genes)))"
                #@debug "TODOX Block: $(block_name) $(reason) $(reason) cluster: $(cluster_index) gene names: [ \"$(join(name_per_gene[indices_of_local_genes[is_in_cluster_per_local_gene]] , "\", \""))\" ]"
                continue
            end

            module_status_per_gene[indices_of_local_genes[is_uncorrelated_pertinent_in_cluster_per_local_gene]] .*= ",uncorrelated_gene"

            if n_cluster_skeletons == 1
                anchor_local_position = local_positions_of_skeletons_in_cluster[1]

            else
                normalized_log_fraction_per_neighborhood_metacell_per_skeleton =
                    @views normalized_log_fraction_per_neighborhood_metacell_per_local_gene[
                        :,
                        local_positions_of_skeletons_in_cluster,
                    ]
                distance_from_center_per_skeleton = flame_timed("colwise.Euclidean") do
                    return colwise(
                        Euclidean(),
                        cluster_center_per_neighborhood_metacell,
                        normalized_log_fraction_per_neighborhood_metacell_per_skeleton,
                    )
                end
                anchor_local_position =
                    local_positions_of_skeletons_in_cluster[argmin(distance_from_center_per_skeleton)]
            end

            module_status_per_gene[indices_of_local_genes[is_correlated_pertinent_in_cluster_per_local_gene]] .*= ",anchor($(name_per_gene[indices_of_local_genes[anchor_local_position]]))"

            anchor_local_position_per_cluster[cluster_index] = anchor_local_position
            cluster_index_of_anchor_local_position[anchor_local_position] = cluster_index

            @debug "TODOX Block: $(block_name) anchor: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(cluster_index) genes: $(n_cluster_genes) lateral: $(n_lateral_cluster_genes) ($(percent(n_lateral_cluster_genes, n_cluster_genes))) uncorrelated pertinent: $(n_uncorrelated_pertinent_cluster_genes) ($(percent(n_uncorrelated_pertinent_cluster_genes, n_cluster_genes))) correlated pertinent: $(n_correlated_pertinent_local_genes) ($(percent(n_correlated_pertinent_local_genes, n_cluster_genes)))"
        end

        @assert length(cluster_index_of_anchor_local_position) > 0
    end

    flame_timed("migrate_between_gene_modules") do
        @views UMIs_per_neighborhood_cell_per_gene = UMIs_per_cell_per_gene[indices_of_neighborhood_cells, :]

        while true
            @debug "TODOX Block $(block_name) MIGRATE..."
            local_positions_of_anchors = collect(keys(cluster_index_of_anchor_local_position))
            cluster_indices_of_anchors = [
                cluster_index_of_anchor_local_position[anchor_local_position] for
                anchor_local_position in local_positions_of_anchors
            ]
            @views center_per_neighborhood_metacell_per_anchor =
                center_per_neighborhood_metacell_per_cluster[:, cluster_indices_of_anchors]

            for local_gene_index in 1:n_local_genes
                if cluster_index_per_local_gene[local_gene_index] == 0
                    @views normalized_log_fraction_per_neighborhood_metacell =
                        normalized_log_fraction_per_neighborhood_metacell_per_local_gene[:, local_gene_index]
                    correlation_per_anchor = flame_timed("cor") do
                        return vec(
                            cor(
                                normalized_log_fraction_per_neighborhood_metacell,
                                center_per_neighborhood_metacell_per_anchor,
                            ),
                        )
                    end
                    @assert_vector(correlation_per_anchor, length(local_positions_of_anchors))
                    anchor_position = argmax(correlation_per_anchor)
                    correlation = correlation_per_anchor[anchor_position]
                    anchor_local_position = local_positions_of_anchors[anchor_position]
                    cluster_index = cluster_index_of_anchor_local_position[anchor_local_position]
                    if correlation >= min_orphan_correlation
                        cluster_index_per_local_gene[local_gene_index] = cluster_index
                        module_status_per_gene[indices_of_local_genes[local_gene_index]] *= ",joined($(name_per_gene[indices_of_local_genes[anchor_local_position]]))"
                        @debug "TODOX Block $(block_name) gene: $(name_per_gene[indices_of_local_genes[local_gene_index]]) is migrated to: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(cluster_index) correlation: $(correlation) status: $(module_status_per_gene[indices_of_local_genes[local_gene_index]])"
                    else
                        @debug "TODOX Block $(block_name) gene: $(name_per_gene[indices_of_local_genes[local_gene_index]]) not migrated to: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(cluster_index) correlation: $(correlation) status: $(module_status_per_gene[indices_of_local_genes[local_gene_index]])"
                    end
                end
            end

            if length(cluster_index_of_anchor_local_position) == 1
                break
            end

            @debug "TODOX Block $(block_name) WEAK..."
            weakest_cluster_index = nothing
            weakest_cluster_strong_UMIs = nothing
            for (anchor_local_position, cluster_index) in cluster_index_of_anchor_local_position
                is_in_cluster_per_local_gene = cluster_index_per_local_gene .== cluster_index
                @assert any(is_in_cluster_per_local_gene)
                @views UMIs_per_neighborhood_cell_per_cluster_local_gene =
                    UMIs_per_neighborhood_cell_per_gene[:, indices_of_local_genes[is_in_cluster_per_local_gene]]
                cluster_UMIs_per_neighborhood_cell =
                    vec(sum(UMIs_per_neighborhood_cell_per_cluster_local_gene; dims = 2))
                cluster_strong_UMIs = quantile(
                    cluster_UMIs_per_neighborhood_cell,
                    1 - (min_strong_cells - 1) / (n_neighborhood_cells - 1),
                )
                @debug "TODOX Block $(block_name) anchor: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(cluster_index) strong UMIs: $(cluster_strong_UMIs)"
                if weakest_cluster_strong_UMIs === nothing || cluster_strong_UMIs < weakest_cluster_strong_UMIs
                    weakest_cluster_index = cluster_index
                    weakest_cluster_strong_UMIs = cluster_strong_UMIs
                end
            end

            @assert weakest_cluster_index !== nothing
            anchor_local_position = anchor_local_position_per_cluster[weakest_cluster_index]

            if weakest_cluster_strong_UMIs >= min_strong_UMIs
                @debug "TODOX Block $(block_name) strong anchor: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(weakest_cluster_index) strong UMIs: $(weakest_cluster_strong_UMIs)"
                break
            end

            @debug "TODOX Block $(block_name) weak anchor: $(name_per_gene[indices_of_local_genes[anchor_local_position]]) cluster: $(weakest_cluster_index) strong UMIs: $(weakest_cluster_strong_UMIs)"
            is_in_weakest_cluster_per_local_gene = cluster_index_per_local_gene .== weakest_cluster_index
            module_status_per_gene[indices_of_local_genes[is_in_weakest_cluster_per_local_gene]] .*= ",weak_cluster"
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

        @debug "TODOX Block: $(block_name) final clusters: $(length(cluster_index_of_anchor_local_position)) anchors: [ $(join(sort(name_per_gene[indices_of_local_genes[collect(keys(cluster_index_of_anchor_local_position))]]), ", ")) ]"
    end

    return genes_indices_of_anchor_index
end

function find_anchor_k(;
    normalized_log_fraction_per_metacell_per_local_gene::AbstractMatrix{<:AbstractFloat},
    cluster_center_per_metacell::AbstractVector{<:AbstractFloat},
    local_positions_of_skeletons_in_cluster::Union{AbstractVector{Bool}, BitVector},
)::Integer
    normalized_log_fraction_per_metacell_per_skeleton =
        @views normalized_log_fraction_per_metacell_per_local_gene[:, local_positions_of_skeletons_in_cluster]
    distance_from_center_per_skeleton = flame_timed("colwise.Euclidean") do
        return colwise(Euclidean(), cluster_center_per_metacell, normalized_log_fraction_per_metacell_per_skeleton)
    end
    return local_positions_of_skeletons_in_cluster[argmin(distance_from_center_per_skeleton)]
end

end  # module
