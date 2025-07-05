"""
Do simple gene module analysis.
"""
module AnalyzeModules

export compute_block_metacell_module_environment_covered_fractions!
export compute_block_metacell_module_environment_linear_fractions!
export compute_block_metacell_module_environment_log_covered_fractions!
export compute_block_metacell_module_environment_log_linear_fractions!
export compute_block_metacell_module_environment_total_UMIs!
export compute_blocks_modules_is_strong!
export compute_blocks_modules_n_genes!
export compute_blocks_modules_n_skeletons!

using Base.Threads
using DataAxesFormats
using TanayLabUtilities

using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_metacell_module_environment_covered_fraction_tensor
import Metacells.Contracts.block_metacell_module_environment_linear_fraction_tensor
import Metacells.Contracts.block_metacell_module_environment_log_covered_fraction_tensor
import Metacells.Contracts.block_metacell_module_environment_log_linear_fraction_tensor
import Metacells.Contracts.block_metacell_module_environment_total_UMIs_tensor
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_module_is_strong_matrix
import Metacells.Contracts.block_module_n_genes_matrix
import Metacells.Contracts.block_module_n_skeletons_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_covered_UMIs_vector
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.module_axis

"""
    function compute_blocks_modules_is_strong!(
        daf::DafWriter;
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        min_strong_UMIs::Integer = 8,
        min_strong_cells::Integer = $(DEFAULT.min_strong_cells),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

A mask of the strong modules that have enough UMIs in enough cells. Enough UMIs means that, if we downsample the cells
to a common total covered UMIs using `downsamples` with the `min_downsamples`, `min_downsamples_quantile`,
`max_downsamples_quantile`, the expected total UMIs of the module will be at least `min_strong_UMIs`. Enough cells mean
that this will happen in at least `min_strong_cells` in the environment. Modules which achieve this are marked as strong
enough to be meaningfully applicable to single cells.

$(CONTRACT)
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
        block_block_is_in_environment_matrix(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_covered_UMIs_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_module_is_strong_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_is_strong!(  # UNTESTED
    daf::DafWriter;
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    min_strong_UMIs::Integer = 8,
    min_strong_cells::Integer = 12,
    overwrite::Bool = false,
)::Nothing
    @assert min_downsamples >= 0
    @assert 0 <= min_downsamples_quantile <= max_downsamples_quantile <= 1
    @assert min_strong_UMIs >= 0
    @assert min_strong_cells >= 0

    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")

    name_per_block = axis_vector(daf, "block")
    name_per_module = axis_vector(daf, "module")

    is_found_per_module_per_block = get_matrix(daf, "module", "block", "is_found")
    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module")
    is_strong_per_module_per_block = zeros(Bool, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        compute_block_modules_is_strong!(
            daf;
            block_index,
            block_name,
            name_per_module,
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
            min_strong_UMIs,
            min_strong_cells,
            module_per_gene_per_block,
            is_found_per_module_per_block,
            is_strong_per_module_per_block,
        )
    end

    set_matrix!(daf, "module", "block", "is_strong", is_strong_per_module_per_block; overwrite)
    return nothing
end

function compute_block_modules_is_strong!(  # UNTESTED
    daf::DafWriter;
    block_index::Integer,
    block_name::AbstractString,
    name_per_module::AbstractVector{<:AbstractString},
    min_downsamples::Integer,
    min_downsamples_quantile::AbstractFloat,
    max_downsamples_quantile::AbstractFloat,
    min_strong_UMIs::Integer,
    min_strong_cells::Integer,
    module_per_gene_per_block::AbstractMatrix{<:AbstractString},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    is_strong_per_module_per_block::AbstractMatrix{Bool},
)::Nothing
    covered_UMIs_per_environment_cell =
        daf["/ cell & metacell ?? => block => is_in_environment ;= $(block_name) : covered_UMIs"].array
    n_environment_cells = length(covered_UMIs_per_environment_cell)

    UMIs_per_gene_per_environment_cell =
        daf["/ gene / cell & metacell ?? => block => is_in_environment ;= $(block_name) : UMIs"].array

    downsampled_UMIs = downsamples(
        covered_UMIs_per_environment_cell;
        min_downsamples,
        min_downsamples_quantile,
        max_downsamples_quantile,
    )
    downscale_per_environment_cell = downsampled_UMIs ./ covered_UMIs_per_environment_cell

    @views module_per_gene_of_block = module_per_gene_per_block[:, block_index]
    @views is_found_per_module_of_block = is_found_per_module_per_block[:, block_index]
    @views is_strong_per_module_of_block = is_strong_per_module_per_block[:, block_index]

    n_modules = length(name_per_module)
    for module_index in 1:n_modules
        if is_found_per_module_of_block[module_index]
            module_name = name_per_module[module_index]
            module_genes_mask = module_per_gene_of_block .== module_name
            @views UMIs_per_module_gene_per_environment_cell = UMIs_per_gene_per_environment_cell[module_genes_mask, :]
            scaled_module_UMIs_per_environment_cell = vec(sum(UMIs_per_module_gene_per_environment_cell; dims = 1))
            scaled_module_UMIs_per_environment_cell =
                scaled_module_UMIs_per_environment_cell .* downscale_per_environment_cell
            @assert_vector(scaled_module_UMIs_per_environment_cell, n_environment_cells)
            n_strong_cells = sum(scaled_module_UMIs_per_environment_cell .>= min_strong_UMIs)
            if n_strong_cells >= min_strong_cells
                is_strong_per_module_of_block[module_index] = true
                @debug "Block: $(block_name) Module: $(module_name) Strong: $(n_strong_cells) Out of: $(n_environment_cells)"
            else
                @debug "Block: $(block_name) Module: $(module_name) Weak: $(n_strong_cells) Out of: $(n_environment_cells)"
            end
        end
    end

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
    data = [block_gene_module_matrix(RequiredInput), block_module_n_genes_matrix(GuaranteedOutput)],
) function compute_blocks_modules_n_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_module_n_genes(daf; name = "n_genes", overwrite)
    return nothing
end

"""
    function compute_blocks_modules_n_skeletons!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of skeleton genes in each gene module in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_n_skeletons_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_n_skeletons!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array
    do_compute_blocks_module_n_genes(daf; name = "n_skeletons", genes_mask = is_skeleton_per_gene, overwrite)
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
    name_per_module = axis_vector(daf, "module")

    n_genes_per_module_per_block = zeros(UInt32, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        module_per_gene = daf["/ gene / block = $(block_name) : module"].array
        for module_index in 1:n_modules
            module_name = name_per_module[module_index]
            module_genes_mask = module_per_gene .== module_name
            if genes_mask !== nothing
                module_genes_mask .&= module_genes_mask
            end
            n_genes_per_module_per_block[module_index, block_index] = sum(module_genes_mask; init = 0)
        end
    end

    set_matrix!(daf, "module", "block", name, bestify(n_genes_per_module_per_block); overwrite)
    return nothing
end

"""
    compute_block_metacell_module_environment_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of each environment module in each environment metacell.
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
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_metacell_module_environment_total_UMIs_tensor(GuaranteedOutput),
    ],
) function compute_block_metacell_module_environment_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_metacells = axis_length(daf, "metacell")
    name_per_block = axis_vector(daf, "block")
    name_per_module = axis_vector(daf, "module")

    UMIs_per_metacell_per_gene = daf["/ metacell / gene : UMIs"].array

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        environment_module_per_gene = daf["/ block = $(block_name) / gene : module"]
        environment_metacell_indices = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"].array
        n_environment_metacells = length(environment_metacell_indices)

        total_UMIs_per_metacell_per_module = zeros(UInt32, n_metacells, n_modules)
        for module_index in n_modules
            module_name = name_per_module[module_index]
            module_genes_mask = environment_module_per_gene .== module_name
            if any(module_genes_mask)
                UMIs_per_environment_metacell_per_module_gene =
                    UMIs_per_metacell_per_gene[environment_metacell_indices, module_genes_mask]
                total_module_UMIs_per_environment_metacell =
                    vec(sum(UMIs_per_environment_metacell_per_module_gene; dims = 2))
                @assert_vector(total_module_UMIs_per_environment_metacell, n_environment_metacells)
                total_UMIs_per_metacell_per_module[environment_metacell_indices, module_index] .=
                    total_module_UMIs_per_environment_metacell
            end
        end

        set_matrix!(
            daf,
            "metacell",
            "module",
            "$(block_name)_total_UMIs",
            bestify(total_UMIs_per_metacell_per_module);
            overwrite,
        )
    end

    return nothing
end

"""
    compute_block_metacell_module_environment_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite)
    )::Nothing

The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total UMIs.
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_metacell_module_environment_total_UMIs_tensor(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        block_metacell_module_environment_linear_fraction_tensor(GuaranteedOutput),
    ],
) function compute_block_metacell_module_environment_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    return do_compute_block_metacell_module_linear_fractions!(
        daf;
        name = "linear",
        denominator = "total_UMIs",
        overwrite,
    )
end

"""
    compute_block_metacell_module_environment_log_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of
the total UMIs.
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_metacell_module_environment_linear_fraction_tensor(RequiredInput),
        block_metacell_module_environment_log_linear_fraction_tensor(GuaranteedOutput),
    ],
) function compute_block_metacell_module_environment_log_linear_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    return do_compute_block_metacell_module_log_linear_fraction(
        daf;
        name = "linear",
        gene_fraction_regularization,
        overwrite,
    )
end

"""
    compute_block_metacell_module_environment_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total covered
UMIs.
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_metacell_module_environment_total_UMIs_tensor(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        block_metacell_module_environment_covered_fraction_tensor(GuaranteedOutput),
    ],
) function compute_block_metacell_module_environment_covered_fractions!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    return do_compute_block_metacell_module_linear_fractions!(
        daf;
        name = "covered",
        denominator = "total_UMIs",
        overwrite,
    )
end

"""
    compute_block_metacell_module_environment_log_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of
the total covered UMIs.
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), module_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_metacell_module_environment_covered_fraction_tensor(RequiredInput),
        block_metacell_module_environment_log_covered_fraction_tensor(GuaranteedOutput),
    ],
) function compute_block_metacell_module_environment_log_covered_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_covered_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    return do_compute_block_metacell_module_log_linear_fraction(
        daf;
        name = "covered",
        gene_fraction_regularization,
        overwrite,
    )
end

function do_compute_block_metacell_module_linear_fractions!(
    daf::DafWriter;
    name::AbstractString,
    denominator::AbstractString,
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_metacells = axis_length(daf, "metacell")
    name_per_block = axis_vector(daf, "block")

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        environment_metacell_indices = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"].array
        total_UMIs_per_module_per_environment_metacell =
            daf["/ module / metacell & block => is_in_environment ;= $(block_name) : $(block_name)_total_UMIs"].array
        total_UMIs_per_environment_metacell =
            daf["/ metacell & block => is_in_environment ;= $(block_name) : $(denominator)"].array
        fraction_per_module_per_environment_metacell =
            total_UMIs_per_module_per_environment_metacell ./ transpose(total_UMIs_per_environment_metacell)
        linear_fraction_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        linear_fraction_per_module_per_metacell[:, environment_metacell_indices] .=
            fraction_per_module_per_environment_metacell
        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_$(name)_fraction",
            bestify(linear_fraction_per_module_per_metacell);
            overwrite,
        )
    end
end

function do_compute_block_metacell_module_log_linear_fraction(
    daf::DafWriter;
    name::AbstractString,
    gene_fraction_regularization::AbstractFloat,
    overwrite::Bool,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_metacells = axis_length(daf, "metacell")
    name_per_block = axis_vector(daf, "block")

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        environment_metacell_indices = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"].array
        is_found_per_module = daf["/ block = $(block_name) / module : is_found"].array
        linear_fraction_per_found_module_per_environment_metacell =
            daf["/ module & is_found ; block = $(block_name) / metacell & block => is_in_environment ;= $(block_name) : $(block_name)_$(name)_fraction"].array
        log_linear_fraction_per_found_module_per_environment_metacell =
            log2.(linear_fraction_per_found_module_per_environment_metacell .+ gene_fraction_regularization)
        log_linear_fraction_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        log_linear_fraction_per_module_per_metacell[is_found_per_module, environment_metacell_indices] .=
            log_linear_fraction_per_found_module_per_environment_metacell
        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_log_$(name)_fraction",
            bestify(log_linear_fraction_per_module_per_metacell);
            overwrite,
        )
    end

    return nothing
end

end  # module
