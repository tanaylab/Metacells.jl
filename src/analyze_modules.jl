"""
Do simple gene module analysis.
"""
module AnalyzeModules

export compute_blocks_modules_covered_UMIs!
export compute_blocks_modules_is_used!
export compute_blocks_modules_n_covered!
export compute_blocks_modules_n_genes!
export compute_blocks_modules_scaled_covered_UMIs!
export compute_blocks_modules_scaled_total_UMIs!
export compute_blocks_modules_total_UMIs!
export compute_blocks_n_found_modules!
export compute_blocks_n_used_modules!
export compute_metacells_modules_covered_UMIs!
export compute_metacells_modules_linear_covered_fractions!
export compute_metacells_modules_linear_fractions!
export compute_metacells_modules_log_linear_covered_fractions!
export compute_metacells_modules_log_linear_fractions!
export compute_metacells_modules_log_scaled_linear_covered_fractions!
export compute_metacells_modules_log_scaled_linear_fractions!
export compute_metacells_modules_scaled_covered_UMIs!
export compute_metacells_modules_scaled_linear_covered_fractions!
export compute_metacells_modules_scaled_linear_fractions!
export compute_metacells_modules_scaled_total_UMIs!
export compute_metacells_modules_total_UMIs!

using Base.Threads
using DataAxesFormats
using TanayLabUtilities

using ..AnalyzeMetacells
using ..Defaults
using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_gene_base_covered_fraction_matrix
import Metacells.Contracts.block_gene_module_index_matrix
import Metacells.Contracts.block_module_covered_UMIs_matrix
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_module_is_used_matrix
import Metacells.Contracts.block_module_n_covered_matrix
import Metacells.Contracts.block_module_n_genes_matrix
import Metacells.Contracts.block_module_scaled_covered_UMIs_matrix
import Metacells.Contracts.block_module_scaled_total_UMIs_matrix
import Metacells.Contracts.block_module_total_UMIs_matrix
import Metacells.Contracts.block_n_found_modules_vector
import Metacells.Contracts.block_n_used_modules_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_covered_UMIs_vector
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_module_covered_UMIs_matrix
import Metacells.Contracts.metacell_module_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_module_linear_fraction_matrix
import Metacells.Contracts.metacell_module_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_module_scaled_linear_fraction_matrix
import Metacells.Contracts.metacell_module_log_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_module_log_linear_fraction_matrix
import Metacells.Contracts.metacell_module_log_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_module_log_scaled_linear_fraction_matrix
import Metacells.Contracts.metacell_module_scaled_covered_UMIs_matrix
import Metacells.Contracts.metacell_module_scaled_total_UMIs_matrix
import Metacells.Contracts.metacell_module_total_UMIs_matrix
import Metacells.Contracts.metacell_scaled_covered_UMIs_vector
import Metacells.Contracts.metacell_scaled_total_UMIs_vector
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.module_axis

"""
    function compute_blocks_n_found_modules!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of found gene modules for each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput)],
    data = [block_module_is_found_matrix(RequiredInput), block_n_found_modules_vector(GuaranteedOutput)],
) function compute_blocks_n_found_modules!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_found_modules_per_block = daf["/ module / block : is_found %> Sum"].array
    set_vector!(daf, "block", "n_found_modules", n_found_modules_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_used_modules!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of used gene modules for each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput)],
    data = [block_module_is_used_matrix(RequiredInput), block_n_used_modules_vector(GuaranteedOutput)],
) function compute_blocks_n_used_modules!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_used_modules_per_block = daf["/ module / block : is_used %> Sum"].array
    set_vector!(daf, "block", "n_used_modules", n_used_modules_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_modules_n_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of genes in each gene module in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [block_gene_module_index_matrix(RequiredInput), block_module_n_genes_matrix(GuaranteedOutput)],
) function compute_blocks_modules_n_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_module_n_genes(daf; name = "n_genes", overwrite)
    return nothing
end

"""
    function compute_blocks_modules_n_covered!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of covered genes in each gene module in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_n_covered_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_n_covered!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_covered_per_gene = get_vector(daf, "gene", "is_covered").array
    do_compute_blocks_module_n_genes(daf; name = "n_covered", genes_mask = is_covered_per_gene, overwrite)
    return nothing
end

"""
    function compute_blocks_modules_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of the genes in each gene module in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_total_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_module_UMIs(daf; name = "total_UMIs", genes_mask = "", overwrite)
    return nothing
end

"""
    function compute_blocks_modules_scaled_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of the genes in each gene module in each block, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_divergence_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_scaled_total_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_scaled_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_gene = get_vector(daf, "gene", "divergence").array
    do_compute_blocks_module_UMIs(
        daf;
        name = "scaled_total_UMIs",
        scale_per_mask_gene = 1.0 .- divergence_per_gene,
        genes_mask = "",
        overwrite,
    )
    return nothing
end

"""
    function compute_blocks_modules_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of the covered genes in each gene module in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_covered_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_module_UMIs(daf; name = "covered_UMIs", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_blocks_modules_scaled_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of the covered genes in each gene module in each block, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_divergence_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        block_module_scaled_covered_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_scaled_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_covered_gene = daf["/ gene & is_covered : divergence"].array
    do_compute_blocks_module_UMIs(
        daf;
        name = "scaled_covered_UMIs",
        genes_mask = "is_covered",
        scale_per_mask_gene = 1.0 .- divergence_per_covered_gene,
        overwrite,
    )
    return nothing
end

"""
    function compute_blocks_modules_is_used!(
        daf::DafWriter;
        min_used_module_covered_UMIs_fraction::Real = $(DEFAULT.min_used_module_covered_UMIs_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Whether each gene module is used by each block. We only use modules that have a high fraction of covered gene UMIs
(specifically, at least `min_used_module_covered_UMIs_fraction`).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        block_module_covered_UMIs_matrix(RequiredInput),
        block_module_total_UMIs_matrix(RequiredInput),
        block_module_is_used_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_is_used!(  # UNTESTED
    daf::DafWriter;
    min_used_module_covered_UMIs_fraction::Real = 0.67,
    overwrite::Bool = false,
)::Nothing
    @assert min_used_module_covered_UMIs_fraction >= 0
    total_UMIs_per_module_per_block = get_matrix(daf, "module", "block", "total_UMIs").array
    covered_UMIs_per_module_per_block = get_matrix(daf, "module", "block", "covered_UMIs").array
    is_used_per_module_per_block =
        covered_UMIs_per_module_per_block .>=
        max.(1, total_UMIs_per_module_per_block .* min_used_module_covered_UMIs_fraction)
    set_matrix!(daf, "module", "block", "is_used", bestify(is_used_per_module_per_block); overwrite)
    return nothing
end

function do_compute_blocks_module_n_genes(
    daf::DafWriter;
    name::AbstractString,
    genes_mask::Maybe{Union{AbstractVector{Bool}, BitVector}} = nothing,
    overwrite::Bool,
)::Nothing
    n_modules = axis_length(daf, "module")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    n_genes_per_module_per_block = zeros(UInt32, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        module_index_per_gene = daf["/ gene / block = $(block_name) : module_index"].array
        for module_index in 1:n_modules
            module_genes_mask = module_index_per_gene .== module_index
            if genes_mask !== nothing
                module_genes_mask .&= module_genes_mask
            end
            n_genes_per_module_per_block[module_index, block_index] = sum(module_genes_mask; init = 0)
        end
    end

    set_matrix!(daf, "module", "block", name, bestify(n_genes_per_module_per_block); overwrite)
    return nothing
end

function do_compute_blocks_module_UMIs(
    daf::DafWriter;
    name::AbstractString,
    genes_mask::AbstractString,
    scale_per_mask_gene::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
    overwrite::Bool,
)::Nothing
    n_modules = axis_length(daf, "module")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    total_UMIs_per_module_per_block = zeros(UInt32, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        if genes_mask == ""
            module_index_per_mask_gene = daf["/ gene / block = $(block_name) : module_index"].array
            UMIs_per_metacell_per_mask_gene =
                daf["/ metacell & block => is_in_environment ;= $(block_name) / gene : UMIs"].array
        else
            module_index_per_mask_gene = daf["/ gene &$(genes_mask) / block = $(block_name) : module_index"].array
            UMIs_per_metacell_per_mask_gene =
                daf["/ metacell & block => is_in_environment ;= $(block_name) / gene &$(genes_mask) : UMIs"].array
        end

        if scale_per_mask_gene !== nothing
            UMIs_per_metacell_per_mask_gene = UMIs_per_metacell_per_mask_gene .* transpose(scale_per_mask_gene)
        end

        for module_index in 1:n_modules
            module_genes_mask = module_index_per_mask_gene .== module_index
            total_UMIs_per_module_per_block[module_index, block_index] =
                round(sum(UMIs_per_metacell_per_mask_gene[:, module_genes_mask]; init = 0))
        end
    end

    set_matrix!(daf, "module", "block", name, bestify(total_UMIs_per_module_per_block); overwrite)
    return nothing
end

"""
    function compute_metacells_modules_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of the genes of each module in each metacell.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        metacell_module_total_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    do_compute_metacells_modules_UMIs(daf; qualifier = "total", genes_mask = "", overwrite)
    return nothing
end

"""
    function compute_metacells_modules_scaled_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of the genes of each module in each metacell, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_divergence_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        metacell_module_scaled_total_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_scaled_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    divergence_per_gene = get_vector(daf, "gene", "divergence").array
    do_compute_metacells_modules_UMIs(
        daf;
        qualifier = "scaled_total",
        genes_mask = "",
        scale_per_mask_gene = 1.0 .- divergence_per_gene,
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_covered_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of the covered genes of each module in each metacell.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        metacell_module_covered_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_covered_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    do_compute_metacells_modules_UMIs(daf; qualifier = "covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_metacells_modules_scaled_covered_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of the covered genes of each module in each metacell, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_divergence_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_gene_module_index_matrix(RequiredInput),
        metacell_module_scaled_covered_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_scaled_covered_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    divergence_per_covered_gene = daf["/ gene & is_covered : divergence"].array
    do_compute_metacells_modules_UMIs(
        daf;
        qualifier = "scaled_covered",
        genes_mask = "is_covered",
        scale_per_mask_gene = 1.0 .- divergence_per_covered_gene,
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_total_UMIs_vector(RequiredInput),
        metacell_module_total_UMIs_matrix(RequiredInput),
        metacell_module_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array;
    total_UMIs_per_metacell_per_module = get_matrix(daf, "metacell", "module", "total_UMIs").array
    linear_fraction_per_metacell_per_module = Float32.(total_UMIs_per_metacell_per_module ./ total_UMIs_per_metacell)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "linear_fraction",
        bestify(linear_fraction_per_metacell_per_module);
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_scaled_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs, scaled by divergence.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_scaled_total_UMIs_vector(RequiredInput),
        metacell_module_scaled_total_UMIs_matrix(RequiredInput),
        metacell_module_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_scaled_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    scaled_total_UMIs_per_metacell = get_vector(daf, "metacell", "scaled_total_UMIs").array;
    scaled_total_UMIs_per_metacell_per_module = get_matrix(daf, "metacell", "module", "scaled_total_UMIs").array
    scaled_linear_fraction_per_metacell_per_module =
        Float32.(scaled_total_UMIs_per_metacell_per_module ./ scaled_total_UMIs_per_metacell)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "scaled_linear_fraction",
        bestify(scaled_linear_fraction_per_metacell_per_module);
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_covered_UMIs_vector(RequiredInput),
        metacell_module_covered_UMIs_matrix(RequiredInput),
        metacell_module_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_linear_covered_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    covered_UMIs_per_metacell = get_vector(daf, "metacell", "covered_UMIs").array;
    covered_UMIs_per_metacell_per_module = get_matrix(daf, "metacell", "module", "covered_UMIs").array
    linear_covered_fraction_per_metacell_per_module =
        Float32.(covered_UMIs_per_metacell_per_module ./ covered_UMIs_per_metacell)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "linear_covered_fraction",
        bestify(linear_covered_fraction_per_metacell_per_module);
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_scaled_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered
UMIs, scaled by divergence.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_scaled_covered_UMIs_vector(RequiredInput),
        metacell_module_scaled_covered_UMIs_matrix(RequiredInput),
        metacell_module_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_scaled_linear_covered_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    scaled_covered_UMIs_per_metacell = get_vector(daf, "metacell", "scaled_covered_UMIs").array;
    scaled_covered_UMIs_per_metacell_per_module = get_matrix(daf, "metacell", "module", "scaled_covered_UMIs").array
    scaled_linear_covered_fraction_per_metacell_per_module =
        Float32.(scaled_covered_UMIs_per_metacell_per_module ./ scaled_covered_UMIs_per_metacell)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "scaled_linear_covered_fraction",
        bestify(scaled_linear_covered_fraction_per_metacell_per_module);
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_log_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = false,
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the
total UMIs. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_module_linear_fraction_matrix(RequiredInput),
        metacell_module_log_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_log_linear_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    linear_fraction_per_metacell_per_module = get_matrix(daf, "metacell", "module", "linear_fraction").array
    log_linear_fraction_per_metacell_per_module =
        log2.(linear_fraction_per_metacell_per_module .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "log_linear_fraction",
        log_linear_fraction_per_metacell_per_module;
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_log_scaled_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = false,
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the
total UMIs, scaled by divergence. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_module_scaled_linear_fraction_matrix(RequiredInput),
        metacell_module_log_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_log_scaled_linear_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    scaled_linear_fraction_per_metacell_per_module =
        get_matrix(daf, "metacell", "module", "scaled_linear_fraction").array
    log_scaled_linear_fraction_per_metacell_per_module =
        log2.(scaled_linear_fraction_per_metacell_per_module .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "log_scaled_linear_fraction",
        log_scaled_linear_fraction_per_metacell_per_module;
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_log_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = false,
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of
the total covered UMIs. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_module_linear_covered_fraction_matrix(RequiredInput),
        metacell_module_log_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_log_linear_covered_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    linear_covered_fraction_per_metacell_per_module =
        get_matrix(daf, "metacell", "module", "linear_covered_fraction").array
    log_linear_covered_fraction_per_metacell_per_module =
        log2.(linear_covered_fraction_per_metacell_per_module .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "log_linear_covered_fraction",
        log_linear_covered_fraction_per_metacell_per_module;
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_modules_log_scaled_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = false,
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of
the total covered UMIs, scaled by divergence. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_module_scaled_linear_covered_fraction_matrix(RequiredInput),
        metacell_module_log_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_log_scaled_linear_covered_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    scaled_linear_covered_fraction_per_metacell_per_module =
        get_matrix(daf, "metacell", "module", "scaled_linear_covered_fraction").array
    log_scaled_linear_covered_fraction_per_metacell_per_module =
        log2.(scaled_linear_covered_fraction_per_metacell_per_module .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "metacell",
        "module",
        "log_scaled_linear_covered_fraction",
        log_scaled_linear_covered_fraction_per_metacell_per_module;
        overwrite,
    )
    return nothing
end

function do_compute_metacells_modules_UMIs(
    daf::DafWriter;
    qualifier::AbstractString,
    genes_mask::AbstractString,
    scale_per_mask_gene::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
    overwrite::Bool,
)::Nothing
    n_modules = axis_length(daf, "module")
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    if genes_mask == ""
        genes_axis_query = "/ gene"
    else
        genes_axis_query = "/ gene & $(genes_mask)"
    end

    UMIs_per_metacell_per_mask_gene = daf["/ metacell $(genes_axis_query) : UMIs"].array
    if scale_per_mask_gene !== nothing
        UMIs_per_metacell_per_mask_gene = UMIs_per_metacell_per_mask_gene .* transpose(scale_per_mask_gene)
    end

    total_UMIs_per_metacell_per_module = Matrix{UInt32}(undef, n_metacells, n_modules)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        block_metacell_indices = daf["/ metacell & block = $(block_name) : index"].array
        module_index_per_mask_gene = daf["/ block = $(block_name) $(genes_axis_query) : module_index"].array
        for module_index in 1:n_modules
            module_gene_indices = findall(module_index_per_mask_gene .== module_index)
            @views UMIs_per_block_metacell_per_module_gene =
                UMIs_per_metacell_per_mask_gene[block_metacell_indices, module_gene_indices]
            total_UMIs_per_metacell_per_module[block_metacell_indices, module_index] .=
                round.(sum(UMIs_per_block_metacell_per_module_gene; dims = 2, init = 0.0))
        end
    end

    set_matrix!(daf, "metacell", "module", "$(qualifier)_UMIs", bestify(total_UMIs_per_metacell_per_module); overwrite)
    return nothing
end

end  # module
