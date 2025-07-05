"""
Group covered genes in each environment to modules, each associated with some anchor gene.
"""
module ComputeModules

export compute_blocks_modules!

using Base.Threads
using Clustering
using DataAxesFormats
using TanayLabUtilities
using StatsBase

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_log_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.module_anchor_gene_vector
import Metacells.Contracts.module_axis

"""
    function compute_blocks_modules!(
        daf::DafWriter;
        mean_genes_per_cluster::Integer = $(DEFAULT.mean_genes_per_cluster),
        anchor_correlation_genes::Integer = $(DEFAULT.anchor_correlation_genes),
        min_anchor_correlation_genes_fraction::AbstractFloat = $(DEFAULT.min_anchor_correlation_genes_fraction),
        max_anchor_correlation_genes_fraction::AbstractFloat = $(DEFAULT.max_anchor_correlation_genes_fraction),
        min_gene_migration_correlation::AbstractFloat = $(DEFAULT.min_gene_migration_correlation),
        max_anchor_kept_cluster_genes_baseline::Integer = $(DEFAULT.max_anchor_kept_cluster_genes_baseline),
        min_anchor_kept_cluster_genes_likelihood::AbstractFloat = $(DEFAULT.min_anchor_kept_cluster_genes_likelihood),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group covered genes in each environment to modules, each associated with some anchor gene.

For each block:

 1. Compute the correlation of the expression of the covered genes across all metacells in the block's environment.
 2. Cluster the genes based on this correlation, using Ward. The number of clusters is such that the mean number of
    genes per cluster is at least `mean_genes_per_cluster`.
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

We name the resulting modules after their anchor gene for interpretability, but add a `.MOD` suffix to clarify this is
the name of a module and not of a gene. It is possible that an anchor in one block will end up being a member of a
different anchor's module in a neighbouring block, or even not to be used in any of the neighbouring block gene modules.

Due to `Daf` limitations, the modules axis must therefore be the union of all the modules from all the blocks, so we
provide an `is_found` mask to specify which module is used for each block.

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
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_is_environment_marker_matrix(RequiredInput),
        metacell_gene_log_covered_fraction_matrix(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        module_anchor_gene_vector(GuaranteedOutput),
        block_module_is_found_matrix(GuaranteedOutput),
        block_gene_module_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules!(
    daf::DafWriter;
    mean_genes_per_cluster::Integer = 200,
    anchor_correlation_genes::Integer = 30,
    min_anchor_correlation_genes_fraction::AbstractFloat = 0.05,
    max_anchor_correlation_genes_fraction::AbstractFloat = 0.5,
    min_gene_migration_correlation::AbstractFloat = 0.7,
    max_anchor_kept_cluster_genes_baseline::Integer = 100,
    min_anchor_kept_cluster_genes_likelihood::AbstractFloat = 1.5,
    overwrite::Bool = false,
)::Nothing
    @assert mean_genes_per_cluster > 0
    @assert 0 <= min_anchor_correlation_genes_fraction <= max_anchor_correlation_genes_fraction <= 1
    @assert anchor_correlation_genes > 0
    @assert 0 <= min_gene_migration_correlation <= 1
    @assert max_anchor_kept_cluster_genes_baseline > 0
    @assert min_anchor_kept_cluster_genes_likelihood >= 0

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    name_per_block = axis_vector(daf, "block")
    name_per_gene = axis_vector(daf, "gene")

    genes_indices_of_anchor_index_per_block = Vector{Dict{Int, Vector{Int}}}(undef, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        genes_indices_of_anchor_index_per_block[block_index] = compute_block_modules!(
            daf;
            block_name,
            mean_genes_per_cluster,
            anchor_correlation_genes,
            min_anchor_correlation_genes_fraction,
            max_anchor_correlation_genes_fraction,
            min_gene_migration_correlation,
            max_anchor_kept_cluster_genes_baseline,
            min_anchor_kept_cluster_genes_likelihood,
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

    if overwrite && has_axis(daf, "module")
        delete_axis!(daf, "module")
    end

    add_axis!(daf, "module", name_per_module)
    set_vector!(daf, "module", "gene.anchor", anchor_per_module)

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
    set_matrix!(daf, "gene", "block", "module", module_per_gene_per_block)
    return nothing
end

function compute_block_modules!(
    daf::DafWriter;
    block_name::AbstractString,
    mean_genes_per_cluster,
    anchor_correlation_genes,
    min_anchor_correlation_genes_fraction,
    max_anchor_correlation_genes_fraction,
    min_gene_migration_correlation,
    max_anchor_kept_cluster_genes_baseline,
    min_anchor_kept_cluster_genes_likelihood,
)::Dict{Int, Vector{Int}}
    #@debug "- Block: $(block_name)"
    return adapter(  # NOJET
        daf;
        input_axes = [
            "metacell" => "/ metacell & block => is_in_environment ;= $(block_name)",
            "gene" => "/ gene & is_environment_marker ; block = $(block_name)",
        ],
        input_data = [
            ("metacell", "gene", "log_covered_fraction") => "=",
            ("gene", "is_skeleton") => "=",
            ("gene", "full_index") => "index",
        ],
        output_axes = [],
        output_data = [],
    ) do adapted
        genes_full_indices_of_anchor_full_index = Dict{Int, Vector{Int}}()

        n_genes = axis_length(adapted, "gene")

        log_fraction_per_metacell_per_gene = get_matrix(adapted, "metacell", "gene", "log_covered_fraction").array
        is_skeleton_per_gene = get_vector(adapted, "gene", "is_skeleton")
        skeleton_indices = findall(is_skeleton_per_gene)
        @assert length(skeleton_indices) > 0

        distances_between_genes = cor(log_fraction_per_metacell_per_gene)  # NOLINT
        @assert_matrix(distances_between_genes, n_genes, n_genes, Columns)
        distances_between_genes .= 1 .- distances_between_genes
        distances_between_genes[isnan.(distances_between_genes)] .= 2

        genes_tree = hclust(distances_between_genes; linkage = :ward)

        n_clusters = Int(ceil(n_genes / mean_genes_per_cluster))
        #@debug "  Clusters: $(n_clusters)"

        cluster_index_per_gene = cutree(genes_tree; k = n_clusters)
        @assert_vector(cluster_index_per_gene, n_genes)

        anchor_index_per_cluster = Vector{Int}(undef, n_clusters)
        anchor_indices = Int[]
        most_correlated_genes_per_cluster = Vector{Maybe{Vector{Int}}}(undef, n_clusters)
        most_correlated_genes_per_cluster .= nothing

        name_per_gene = axis_vector(adapted, "gene")
        for cluster_index in 1:n_clusters
            #@debug "  - Cluster: $(cluster_index)"
            anchor_index, n_original_genes, most_correlated_genes_of_anchor = find_cluster_anchor(;  # NOLINT
                cluster_index,
                distances_between_genes,
                skeleton_indices,
                cluster_index_per_gene,
                anchor_correlation_genes,
                min_anchor_correlation_genes_fraction,
                max_anchor_correlation_genes_fraction,
                max_anchor_kept_cluster_genes_baseline,
            )
            #@debug "    Anchor: $(anchor_index == 0 ? nothing : name_per_gene[anchor_index]) Genes: $(n_original_genes)"
            anchor_index_per_cluster[cluster_index] = anchor_index
            if anchor_index > 0
                push!(anchor_indices, anchor_index)
                most_correlated_genes_per_cluster[cluster_index] = most_correlated_genes_of_anchor
            end
        end

        #@debug "  Patch:"
        for gene_index in 1:n_genes
            patch_cluster_index_per_gene!(;
                cluster_index_per_gene,
                anchor_index_per_cluster,
                gene_index,
                distances_between_genes,
                min_gene_migration_correlation,
                anchor_indices,
                name_per_gene,
            )
        end

        n_anchors = length(anchor_indices)
        @assert n_anchors > 0
        min_kept_genes_fraction = min_anchor_kept_cluster_genes_likelihood / n_anchors

        full_index_per_gene = get_vector(adapted, "gene", "full_index")

        #@debug "  Purge:"
        for cluster_index in 1:n_clusters
            anchor_index = anchor_index_per_cluster[cluster_index]
            if anchor_index > 0
                most_correlated_genes_of_anchor = most_correlated_genes_per_cluster[cluster_index]
                n_most_correlated_genes = length(most_correlated_genes_of_anchor)
                n_kept_genes = sum(cluster_index_per_gene[most_correlated_genes_of_anchor] .== cluster_index)
                if n_kept_genes >= n_most_correlated_genes * min_kept_genes_fraction
                    #@debug "  - Keep Anchor: $(name_per_gene[anchor_index]) Original: $(n_most_correlated_genes) Kept: $(n_kept_genes) Min: $(n_most_correlated_genes * min_kept_genes_fraction)"
                    anchor_full_index = full_index_per_gene[anchor_index]
                    genes_full_indices = full_index_per_gene[cluster_index_per_gene .== cluster_index]
                    genes_full_indices_of_anchor_full_index[anchor_full_index] = genes_full_indices
                else
                    #@debug "  - Discard Anchor: $(name_per_gene[anchor_index]) Original: $(n_most_correlated_genes) Kept: $(n_kept_genes) Min: $(n_most_correlated_genes * min_kept_genes_fraction)"
                    n_anchors -= 1
                end
            end
        end
        @debug "- Block: $(block_name) Anchors: $(n_anchors) [ $(join(name_per_gene[anchor_indices], ", ")) ]"
        return genes_full_indices_of_anchor_full_index
    end
end

function find_cluster_anchor(;
    cluster_index::Integer,
    cluster_index_per_gene::AbstractVector{<:Integer},
    distances_between_genes::AbstractMatrix{<:AbstractFloat},
    skeleton_indices::AbstractVector{<:Integer},
    anchor_correlation_genes::Integer,
    min_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_correlation_genes_fraction::AbstractFloat,
    max_anchor_kept_cluster_genes_baseline::Integer,
)::Tuple{Int, Int, Vector{Int}}
    is_cluster_gene = cluster_index_per_gene .== cluster_index
    n_cluster_genes = sum(is_cluster_gene)
    @assert n_cluster_genes > 0

    anchor_correlation_quantile = clamp(
        anchor_correlation_genes / n_cluster_genes,
        min_anchor_correlation_genes_fraction,
        max_anchor_correlation_genes_fraction,
    )

    #@debug "    Anchor quantile: $(anchor_correlation_quantile) Genes: $(anchor_correlation_quantile * n_cluster_genes) Out of: $(n_cluster_genes)"

    anchor_index = 0
    anchor_distance = 3

    for skeleton_index in skeleton_indices
        if is_cluster_gene[skeleton_index]
            @views distances_from_skeleton = distances_between_genes[:, skeleton_index]
            skeleton_distance = quantile(distances_from_skeleton, anchor_correlation_quantile)  # NOLINT
            if skeleton_distance < anchor_distance
                anchor_distance = skeleton_distance
                anchor_index = skeleton_index
            end
        end
    end

    if anchor_index == 0
        most_correlated_genes_of_anchor = Int[]
    else
        gene_indices_of_cluster = findall(cluster_index_per_gene .== cluster_index)
        if n_cluster_genes <= max_anchor_kept_cluster_genes_baseline
            most_correlated_genes_of_anchor = gene_indices_of_cluster
        else
            distances_of_cluster_genes_from_anchor = vec(distances_between_genes[gene_indices_of_cluster, anchor_index])
            threshold = quantile(  # NOLINT
                distances_of_cluster_genes_from_anchor,
                max_anchor_kept_cluster_genes_baseline / n_cluster_genes,
            )
            most_correlated_genes_of_anchor =
                gene_indices_of_cluster[distances_of_cluster_genes_from_anchor .<= threshold]
        end
    end

    #@debug "    Anchor Distance: $(anchor_distance) Correlation: $(1 - anchor_distance)"

    return (anchor_index, n_cluster_genes, most_correlated_genes_of_anchor)
end

function patch_cluster_index_per_gene!(;
    cluster_index_per_gene::AbstractVector{<:Integer},
    anchor_index_per_cluster::AbstractVector{<:Integer},
    gene_index::Integer,
    distances_between_genes::AbstractMatrix{<:AbstractFloat},
    min_gene_migration_correlation::AbstractFloat,
    anchor_indices::AbstractVector{<:Integer},
    name_per_gene::AbstractVector{<:AbstractString},  # NOLINT TODOX
)::Nothing
    cluster_index = cluster_index_per_gene[gene_index]
    anchor_index = anchor_index_per_cluster[cluster_index]  # NOLINT TODOX

    anchors_distances_from_gene = distances_between_genes[anchor_indices, gene_index]
    closest_anchor_index = anchor_indices[argmin(anchors_distances_from_gene)]
    closest_cluster_index = cluster_index_per_gene[closest_anchor_index]

    if closest_cluster_index != cluster_index
        if distances_between_genes[closest_anchor_index, gene_index] <= 1 - min_gene_migration_correlation
            #@debug "    - Move Gene: $(name_per_gene[gene_index]) from Cluster: $(cluster_index) Anchor: $(anchor_index == 0 ? nothing : name_per_gene[anchor_index]) Correlation: $(anchor_index == 0 ? -2 : 1 - distances_between_genes[anchor_index, gene_index]) to  Cluster: $(closest_cluster_index) Anchor: $(name_per_gene[closest_anchor_index]) Correlation: $(1 - distances_between_genes[closest_anchor_index, gene_index])"
            cluster_index_per_gene[gene_index] = closest_cluster_index
        else
            #@debug "    - Drop Gene: $(name_per_gene[gene_index]) from Cluster: $(cluster_index) Anchor: $(anchor_index == 0 ? nothing : name_per_gene[anchor_index]) Correlation: $(anchor_index == 0 ? -2 : 1 - distances_between_genes[anchor_index, gene_index]) because Cluster: $(closest_cluster_index) Anchor: $(name_per_gene[closest_anchor_index]) Correlation: $(1 - distances_between_genes[closest_anchor_index, gene_index])"
            cluster_index_per_gene[gene_index] = 0
        end
    end

    return nothing
end

end  # module
