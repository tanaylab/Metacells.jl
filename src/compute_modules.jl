"""
Group "very correlated" genes modules in each local environment.
"""
module ComputeModules

export compute_blocks_modules!

using Base.Threads
using Clustering
using DataAxesFormats
using Distributions
using StatsBase
using TanayLabUtilities

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_module_min_gene_correlation_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.module_axis

"""
    function compute_blocks_modules!(
        daf::DafWriter;
        min_module_downsampled_UMIs::Integer = $(DEFAULT.min_module_downsampled_UMIs),
        min_module_strong_cells::Integer = $(DEFAULT.min_module_strong_cells),
        min_merge_fragments_correlation::Real = $(DEFAULT.min_merge_fragments_correlation),
        min_merge_modules_correlation::Real = $(DEFAULT.min_merge_modules_correlation),
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Group the genes into modules based on their correlation.

For each block:

 1. Compute the correlation of the expression of the genes across all metacells in the block's environment.
 2. Cluster the genes based on this correlation, using the `:complete` distance measure.
 3. Bottom-up group genes into fragments and then into complete modules. Fragments have at least
    `min_merge_fragments_correlation` between all their genes. At some point, the fragment will contain enough genes so
    there will be at least `min_module_strong_cells`. A "strong" cell is expected to have at least
    `min_module_downsampled_UMIs`, if we downsample the cells in the environment using `min_downsamples`,
    `min_downsamples_quantile`, and `max_downsamples_quantile`. If the `divergence` is specified per gene, then it is
    applied to the UMIs. When this happens, we mark the fragment (which might be a single gene!) as a valid module. We
    will only merge two such modules if the correlation between all their genes is at least
    `min_merge_modules_correlation`.
 4. Having merged some of the genes into modules, we repeat steps 2 and 3 for the unmerged genes until no new modules
    are found.
 5. Most of the genes will not end up in a valid module, and are given the module index 0. The other genes are given a
    1-based module index. We need to store a matrix of module indices instead of module names due to `Daf` limitations
    (no matrices of strings).

!!! note

    This uses the virtual [`metacell_gene_fraction_matrix`](@ref). You will need an `adapter` to map these to concrete
    fractions (geomean, linear, scaled, ...). In addition, you might want to restrict the set of genes that are
    candidates to be included in the module, for example just to `is_covered` genes.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(GuaranteedOutput),
    ],
    data = [
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_fraction_matrix(RequiredInput),
        gene_divergence_vector(OptionalInput),
        block_gene_is_environment_marker_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_matrix(GuaranteedOutput),
        block_module_anchor_gene_matrix(GuaranteedOutput),
        block_module_is_found_matrix(GuaranteedOutput),
        block_module_min_gene_correlation_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules!(
    daf::DafWriter;
    min_module_downsampled_UMIs::Integer = 8,
    min_module_strong_cells::Integer = 12,
    min_merge_fragments_correlation::Real = 0.5,
    min_merge_modules_correlation::Real = 0.8,
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    overwrite::Bool = false,
)::Nothing
    @assert min_module_downsampled_UMIs >= 0
    @assert min_module_strong_cells >= 0
    @assert 0 < min_merge_fragments_correlation < 1
    @assert 0 < min_merge_modules_correlation < 1

    n_gene = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    gene_indices_per_module_per_block = Vector{Vector{Vector{UInt32}}}(undef, n_blocks)
    min_gene_correlation_per_module_vector_per_block = Vector{Vector{Float32}}(undef, n_blocks)

    n_modules = 0
    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        gene_names_per_module, min_gene_correlation_per_module = compute_block_modules!(
            daf;
            block_name,
            min_module_downsampled_UMIs,
            min_module_strong_cells,
            min_merge_fragments_correlation,
            min_merge_modules_correlation,
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )
        n_modules = max(n_modules, length(gene_names_per_module))
        gene_indices_per_module_per_block[block_index] =
            [axis_indices(daf, "gene", gene_names_of_module) for gene_names_of_module in gene_names_per_module]
        min_gene_correlation_per_module_vector_per_block[block_index] = min_gene_correlation_per_module
    end

    if overwrite
        delete_axis!(daf, "module"; must_exist = false)
    end
    name_per_module = "GM" .* string.(collect(1:n_modules))
    add_axis!(daf, "module", name_per_module)

    module_per_gene_per_block = fill("", n_gene, n_blocks)
    is_found_per_module_per_block = zeros(Bool, n_modules, n_blocks)
    min_gene_correlation_per_module_per_block = zeros(Float32, n_modules, n_blocks)
    for (block_index, gene_indices_per_module_of_block) in enumerate(gene_indices_per_module_per_block)
        for (module_index, gene_indices_of_module_of_block) in enumerate(gene_indices_per_module_of_block)
            module_per_gene_per_block[gene_indices_of_module_of_block, block_index] .= name_per_module[module_index]
            is_found_per_module_per_block[module_index, block_index] = true
        end
        min_gene_correlation_per_module = min_gene_correlation_per_module_vector_per_block[block_index]
        min_gene_correlation_per_module_per_block[1:length(min_gene_correlation_per_module), block_index] .=
            min_gene_correlation_per_module
    end

    set_matrix!(daf, "gene", "block", "module", module_per_gene_per_block)
    set_matrix!(daf, "module", "block", "is_found", bestify(is_found_per_module_per_block))  # NOJET
    set_matrix!(daf, "module", "block", "min_gene_correlation", bestify(min_gene_correlation_per_module_per_block))  # NOJET
    return nothing
end

function compute_block_modules!(
    daf;
    block_name::AbstractString,
    min_module_downsampled_UMIs::Integer,
    min_module_strong_cells::Integer,
    min_merge_fragments_correlation::Real,
    min_merge_modules_correlation::Real,
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
)::Tuple{Vector{Vector{AbstractString}}, Vector{Float32}}
    use_divergence = has_vector(daf, "gene", "divergence")
    return adapter(  # NOJET
        daf;
        input_axes = [
            "metacell" => "/ metacell & block => is_in_environment ;= $(block_name)",
            "cell" => "/ cell & metacell ?? => block => is_in_environment ;= $(block_name)",
            "gene" => "/ gene & is_environment_marker ; block = $(block_name)",
        ],
        input_data = [
            ("cell", "gene", "UMIs") => "=",
            ("metacell", "gene", "fraction") => "=",
            ("gene", "divergence") => use_divergence ? "=" : nothing,
        ],
        output_axes = [],
        output_data = [],
    ) do adapted
        n_cells = axis_length(adapted, "cell")
        n_genes = axis_length(adapted, "gene")
        name_per_gene = axis_vector(adapted, "gene")

        fraction_per_metacell_per_gene = get_matrix(adapted, "metacell", "gene", "fraction").array
        distances_between_genes = cor(fraction_per_metacell_per_gene)
        @assert_matrix(distances_between_genes, n_genes, n_genes, Columns)
        distances_between_genes .= 1 .- distances_between_genes
        distances_between_genes[isnan.(distances_between_genes)] .= 2

        UMIs_per_cell_per_gene = adapted["/ cell / gene : UMIs"].array
        if use_divergence
            divergence_per_gene = get_vector(adapted, "gene", "divergence").array
            scale_per_gene = 1.0 .- divergence_per_gene
            UMIs_per_cell_per_gene = UMIs_per_cell_per_gene .* transpose(scale_per_gene)
        else
            scale_per_gene = nothing
        end
        total_UMIs_per_cell = UInt32.(round.(vec(sum(UMIs_per_cell_per_gene; dims = 2))))
        @assert_vector(total_UMIs_per_cell, n_cells)
        downsampled_UMIs =
            downsamples(total_UMIs_per_cell; min_downsamples, min_downsamples_quantile, max_downsamples_quantile)
        @debug "Downsampled UMIs: $(downsampled_UMIs)"
        scale_UMIs_per_cell = total_UMIs_per_cell ./ downsampled_UMIs

        gene_names_per_module = Vector{AbstractVector{<:AbstractString}}()
        min_gene_correlation_per_module = Float32[]

        remaining_genes_mask = ones(Bool, n_genes)
        remaining_genes_indices = collect(1:n_genes)

        round_index = 0
        while length(remaining_genes_indices) > 0
            round_index += 1
            @debug "- Block: $(block_name) Round: $(round_index)"
            fragments = cluster_fragments(
                # name_per_gene,
                remaining_genes_indices,
                distances_between_genes,
                UMIs_per_cell_per_gene,
                scale_per_gene,
                scale_UMIs_per_cell,
                min_module_downsampled_UMIs,
                min_module_strong_cells,
                min_merge_fragments_correlation,
                min_merge_modules_correlation,
            )

            if collect_modules!(;
                block_name,
                name_per_gene,
                gene_names_per_module,
                min_gene_correlation_per_module,
                remaining_genes_mask,
                fragments,
            )
                remaining_genes_indices = findall(remaining_genes_mask)
            else
                break
            end
        end

        return (gene_names_per_module, min_gene_correlation_per_module)
    end
end

@kwdef mutable struct Fragment
    genes_indices::AbstractVector{<:Integer}
    total_scaled_UMIs_per_cell::AbstractVector{<:AbstractFloat}
    min_correlation_between_genes::AbstractFloat
    is_module::Bool
end

function collect_modules!(;
    block_name::AbstractString,
    name_per_gene::AbstractVector{<:AbstractString},
    gene_names_per_module::Vector{AbstractVector{<:AbstractString}},
    min_gene_correlation_per_module::Vector{Float32},
    remaining_genes_mask::Union{AbstractVector{Bool}, BitVector},
    fragments::AbstractVector{Maybe{Fragment}},
)::Bool
    did_collect = false
    for fragment in fragments
        if fragment !== nothing && fragment.is_module
            @debug (
                "  - Block: $(block_name)" *
                " Module: $(length(gene_names_per_module) + 1)" *
                " Genes: ($(length(fragment.genes_indices)))" *
                " Corr: $(fragment.min_correlation_between_genes)" *
                " Names: $(join(sort(name_per_gene[fragment.genes_indices]), " "))"
            )
            @assert all(remaining_genes_mask[fragment.genes_indices])
            remaining_genes_mask[fragment.genes_indices] .= false
            @views gene_names_of_module = name_per_gene[fragment.genes_indices]
            push!(gene_names_per_module, gene_names_of_module)
            push!(min_gene_correlation_per_module, fragment.min_correlation_between_genes)
            did_collect = true
        end
    end
    return did_collect
end

function cluster_fragments(
    # name_per_gene::AbstractVector{<:AbstractString},
    remaining_genes_indices::AbstractVector{<:Integer},
    distances_between_genes::AbstractMatrix{<:AbstractFloat},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Real},
    scale_per_gene::Maybe{AbstractVector{<:AbstractFloat}},
    scale_UMIs_per_cell::AbstractVector{<:AbstractFloat},
    min_module_downsampled_UMIs::Integer,
    min_module_strong_cells::Integer,
    min_merge_fragments_correlation::AbstractFloat,
    min_merge_modules_correlation::AbstractFloat,
)::AbstractVector{Maybe{Fragment}}
    n_remaining_genes = length(remaining_genes_indices)
    n_merges = n_remaining_genes - 1
    @views distances_between_remaining_genes = distances_between_genes[remaining_genes_indices, remaining_genes_indices]
    tree = hclust(distances_between_remaining_genes; linkage = :complete)
    @assert_matrix(tree.merges, n_merges, 2)

    fragments = Vector{Maybe{Fragment}}(undef, n_remaining_genes + n_merges)
    fragments .= nothing

    for merge_index in 1:n_merges
        left_index, right_index = tree.merges[merge_index, :]

        # @debug "  - Left: $(left_index) $(left_index < 0 ? name_per_gene[remaining_genes_indices[-left_index]] : "Module")"
        left_fragment = get_fragment(;
            fragments,
            index = left_index,
            remaining_genes_indices,
            UMIs_per_cell_per_gene,
            scale_per_gene,
            scale_UMIs_per_cell,
            min_module_downsampled_UMIs,
            min_module_strong_cells,
        )
        @assert left_fragment === nothing || length(left_fragment.genes_indices) > 0

        # @debug "    Right: $(right_index) $(right_index < 0 ? name_per_gene[remaining_genes_indices[-right_index]] : "Module")"
        right_fragment = get_fragment(;
            fragments,
            index = right_index,
            remaining_genes_indices,
            UMIs_per_cell_per_gene,
            scale_per_gene,
            scale_UMIs_per_cell,
            min_module_downsampled_UMIs,
            min_module_strong_cells,
        )
        @assert right_fragment === nothing || length(right_fragment.genes_indices) > 0

        min_correlation_between_genes = 1 - tree.heights[merge_index]

        merged_fragment = merge_fragments(
            left_fragment,
            right_fragment;
            min_correlation_between_genes,
            min_merge_fragments_correlation,
            min_merge_modules_correlation,
            min_module_downsampled_UMIs,
            min_module_strong_cells,
        )
        @assert merged_fragment === nothing || length(merged_fragment.genes_indices) > 0
        if merged_fragment === nothing
            # @debug "    Apart: $(merge_index) min_correlation_between_genes: $(min_correlation_between_genes)"
        else
            # @debug "    Merge: $(merge_index) min_correlation_between_genes: $(min_correlation_between_genes)"
        end
        fragments[merge_index + n_remaining_genes] = merged_fragment
    end

    return fragments
end

function get_fragment(;
    fragments::AbstractVector{Maybe{Fragment}},
    index::Integer,
    remaining_genes_indices::AbstractVector{<:Integer},
    UMIs_per_cell_per_gene::AbstractMatrix{<:Real},
    scale_per_gene::Maybe{AbstractVector{<:AbstractFloat}},
    scale_UMIs_per_cell::AbstractVector{<:AbstractFloat},
    min_module_downsampled_UMIs::Integer,
    min_module_strong_cells::Integer,
)::Maybe{Fragment}
    n_remaining_genes = length(remaining_genes_indices)
    @assert index != 0

    if index > 0
        return fragments[index + n_remaining_genes]
    end

    gene_index = remaining_genes_indices[-index]
    genes_indices = [gene_index]
    @views total_UMIs_per_cell = UMIs_per_cell_per_gene[:, gene_index]
    total_scaled_UMIs_per_cell = total_UMIs_per_cell ./ scale_UMIs_per_cell
    if scale_per_gene !== nothing
        total_scaled_UMIs_per_cell .*= scale_per_gene[gene_index]
    end

    is_module = sum(total_scaled_UMIs_per_cell .>= min_module_downsampled_UMIs) >= min_module_strong_cells
    fragment = Fragment(; genes_indices, total_scaled_UMIs_per_cell, is_module, min_correlation_between_genes = 1.0)

    fragments[-index] = fragment
    return fragment
end

function merge_fragments(
    ::Nothing,
    ::Nothing;
    min_correlation_between_genes::AbstractFloat,  # NOLINT
    min_merge_fragments_correlation::AbstractFloat,  # NOLINT
    min_merge_modules_correlation::AbstractFloat,  # NOLINT
    min_module_downsampled_UMIs::Integer,  # NOLINT
    min_module_strong_cells::Integer,  # NOLINT
)::Nothing
    return nothing
end

function merge_fragments(
    ::Fragment,
    ::Nothing;
    min_correlation_between_genes::AbstractFloat,  # NOLINT
    min_merge_fragments_correlation::AbstractFloat,  # NOLINT
    min_merge_modules_correlation::AbstractFloat,  # NOLINT
    min_module_downsampled_UMIs::Integer,  # NOLINT
    min_module_strong_cells::Integer,  # NOLINT
)::Nothing
    return nothing
end

function merge_fragments(
    ::Nothing,
    ::Fragment;
    min_correlation_between_genes::AbstractFloat,  # NOLINT
    min_merge_fragments_correlation::AbstractFloat,  # NOLINT
    min_merge_modules_correlation::AbstractFloat,  # NOLINT
    min_module_downsampled_UMIs::Integer,  # NOLINT
    min_module_strong_cells::Integer,  # NOLINT
)::Nothing
    return nothing
end

function merge_fragments(
    left_fragment::Fragment,
    right_fragment::Fragment;
    min_merge_fragments_correlation::AbstractFloat,
    min_merge_modules_correlation::AbstractFloat,
    min_correlation_between_genes::AbstractFloat,
    min_module_downsampled_UMIs::Integer,
    min_module_strong_cells::Integer,
)::Maybe{Fragment}
    if left_fragment.is_module != right_fragment.is_module
        return nothing
    end

    min_merge_correlation = left_fragment.is_module ? min_merge_modules_correlation : min_merge_fragments_correlation
    if min_correlation_between_genes < min_merge_correlation
        return nothing
    end

    total_scaled_UMIs_per_cell = left_fragment.total_scaled_UMIs_per_cell + right_fragment.total_scaled_UMIs_per_cell
    genes_indices = vcat(left_fragment.genes_indices, right_fragment.genes_indices)

    is_module = sum(total_scaled_UMIs_per_cell .>= min_module_downsampled_UMIs) >= min_module_strong_cells
    fragment = Fragment(; genes_indices, total_scaled_UMIs_per_cell, is_module, min_correlation_between_genes)

    if is_module
        left_fragment.is_module = false
        right_fragment.is_module = false
    end

    return fragment
end

end  # module
