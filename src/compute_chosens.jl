"""
Condense the local gene modules to global gene chosens.
"""
module ComputeChosens

export compute_chosens!

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
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.module_anchor_gene_vector
import Metacells.Contracts.module_axis

"""
    function compute_chosens!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        gene_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        chosen_axis(GuaranteedOutput),
    ],
    data = [
        gene_is_regulator_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_is_strong_matrix(RequiredInput),
        chosen_block_vector(GuaranteedOutput),
        chosen_module_vector(GuaranteedOutput),
        block_module_chosen_matrix(GuaranteedOutput),
        chosen_gene_is_member_matrix(GuaranteedOutput),
    ],
) function compute_chosens!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")
    n_modules = axis_length(daf, "module")

    name_per_gene = axis_vector(daf, "gene")
    name_per_block = axis_vector(daf, "block")
    name_per_module = axis_vector(daf, "module")

    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module").array

    is_strong_per_module_per_block = get_matrix(daf, "module", "block", "is_strong").array
    n_strong_modules = sum(is_strong_per_module_per_block)
    @debug "Strong blocks modules: $(n_strong_modules)"

    is_member_of_strong_modules_per_gene = zeros(Bool, n_genes)
    is_member_per_gene_per_strong_module = zeros(Bool, n_genes, n_strong_modules)
    block_index_per_strong_module = Vector{UInt32}(undef, n_strong_modules)
    module_index_per_strong_module = Vector{UInt32}(undef, n_strong_modules)

    m_strong_modules = 0
    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        @views module_per_gene = module_per_gene_per_block[:, block_index]
        @views is_strong_per_module = is_strong_per_module_per_block[:, block_index]
        @debug "- Block $(block_name) Strong modules: $(sum(is_strong_per_module))"
        @assert sum(is_strong_per_module) > 0
        for module_index in findall(is_strong_per_module)
            m_strong_modules += 1
            module_name = name_per_module[module_index]
            is_in_module_per_gene = module_per_gene .== module_name
            is_member_per_gene_per_strong_module[is_in_module_per_gene, m_strong_modules] .= true
            block_index_per_strong_module[m_strong_modules] = block_index
            module_index_per_strong_module[m_strong_modules] = module_index
            @debug "  - Module $(module_name) genes: $(sum(is_in_module_per_gene))"
        end
    end
    @assert m_strong_modules == n_strong_modules

    is_member_of_strong_modules_per_gene = vec(sum(is_member_per_gene_per_strong_module; dims = 2) .> 0)
    @assert_vector(is_member_of_strong_modules_per_gene, n_genes)

    n_strong_module_genes = sum(is_member_of_strong_modules_per_gene)
    @debug "n_strong_modules: $(n_strong_modules)"
    is_in_any_module_per_gene = vec(any(module_per_gene_per_block .!= ""; dims = 2))
    @assert_vector(is_in_any_module_per_gene, n_genes)

    distances_between_strong_modules = Matrix{Float32}(undef, n_strong_modules, n_strong_modules)
    distances_between_strong_modules[1, 1] = 0
    parallel_loop_wo_rng(
        reverse(2:n_strong_modules);
        progress = DebugProgress(n_strong_modules - 1),
    ) do base_strong_module_index
        distances_between_strong_modules[base_strong_module_index, base_strong_module_index] = 0
        @views is_member_per_gene_of_base_strong_module =
            is_member_per_gene_per_strong_module[:, base_strong_module_index]

        for other_strong_module_index in 1:base_strong_module_index
            @views is_member_per_gene_of_other_strong_module =
                is_member_per_gene_per_strong_module[:, other_strong_module_index]

            n_either_genes = sum(is_member_per_gene_of_other_strong_module .| is_member_per_gene_of_base_strong_module)
            @assert n_either_genes <= n_strong_module_genes

            n_both_genes = sum(is_member_per_gene_of_other_strong_module .& is_member_per_gene_of_base_strong_module)
            @assert n_both_genes <= n_either_genes

            n_different_genes = n_either_genes - n_both_genes
            distance = n_different_genes / n_either_genes

            distances_between_strong_modules[base_strong_module_index, other_strong_module_index] = distance
            distances_between_strong_modules[other_strong_module_index, base_strong_module_index] = distance
        end

        return nothing
    end

    chosens = hclust(distances_between_strong_modules; linkage = :ward)  # NOJET
    n_chosens = n_modules

    chosen_index_per_strong_module = cutree(chosens; k = n_chosens)

    block_per_chosen = Vector{AbstractString}(undef, n_chosens)
    module_per_chosen = Vector{AbstractString}(undef, n_chosens)
    is_member_per_gene_per_chosen = zeros(Bool, n_genes, n_chosens)

    fraction_used_per_gene_per_chosen = Matrix{Float32}(undef, n_genes, n_chosens)

    for chosen_index in 1:n_chosens
        indices_of_strong_modules = findall(chosen_index_per_strong_module .== chosen_index)
        distances_between_chosen_modules =
            distances_between_strong_modules[indices_of_strong_modules, indices_of_strong_modules]
        max_distance_per_chosen_module = vec(maximum(distances_between_chosen_modules; dims = 1))
        center_strong_module_index = indices_of_strong_modules[argmin(max_distance_per_chosen_module)]

        @views is_member_per_gene_per_chosen_module = is_member_per_gene_per_strong_module[:, indices_of_strong_modules]
        n_used_per_gene = vec(sum(is_member_per_gene_per_chosen_module; dims = 2))
        @assert_vector(n_used_per_gene, n_genes)
        fraction_used_per_gene_per_chosen[:, chosen_index] = n_used_per_gene ./ length(indices_of_strong_modules)

        block_per_chosen[chosen_index] = name_per_block[block_index_per_strong_module[center_strong_module_index]]
        module_per_chosen[chosen_index] = name_per_module[module_index_per_strong_module[center_strong_module_index]]
        is_member_per_gene_per_chosen[:, chosen_index] =
            is_member_per_gene_per_strong_module[:, center_strong_module_index]
    end

    indices_of_regulators = daf["/ gene & is_regulator : index"].array
    fraction_used_per_gene_per_chosen[indices_of_regulators, :] .*= 2.0

    name_per_chosen = fill("", n_chosens)
    is_used_as_name_per_gene = zeros(Bool, n_genes)

    order = sortperm(vec(fraction_used_per_gene_per_chosen); rev = true)
    n_named_chosens = 0
    for entry_index in order
        chosen_index, gene_index = divrem(entry_index - 1, n_genes)
        chosen_index += 1
        gene_index += 1
        if name_per_chosen[chosen_index] == "" && !is_used_as_name_per_gene[gene_index]
            name_per_chosen[chosen_index] = name_per_gene[gene_index] * ".CO"
            is_used_as_name_per_gene[gene_index] = true
            n_named_chosens += 1
            if n_named_chosens == n_chosens
                break
            end
        end
    end
    @assert !any(name_per_chosen .== "")

    chosen_per_module_per_block = fill("", n_modules, n_blocks)
    for chosen_index in 1:n_chosens
        indices_of_strong_modules = findall(chosen_index_per_strong_module .== chosen_index)
        for strong_module_index in indices_of_strong_modules
            chosen_per_module_per_block[
                module_index_per_strong_module[strong_module_index],
                block_index_per_strong_module[strong_module_index],
            ] = name_per_chosen[chosen_index]
        end
    end

    add_axis!(daf, "chosen", name_per_chosen; overwrite)
    set_vector!(daf, "chosen", "block", block_per_chosen; overwrite)
    set_vector!(daf, "chosen", "module", module_per_chosen; overwrite)
    set_matrix!(daf, "module", "block", "chosen", chosen_per_module_per_block; overwrite)
    set_matrix!(daf, "gene", "chosen", "is_member", bestify(is_member_per_gene_per_chosen); overwrite)
    return nothing
end

end  # module
