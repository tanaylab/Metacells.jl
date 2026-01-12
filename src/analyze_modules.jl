"""
Do simple gene module analysis.
"""
module AnalyzeModules

export compute_blocks_modules_is_strong!
export compute_blocks_modules_n_genes!
export compute_blocks_modules_n_skeletons!
export compute_blocks_modules_neighborhood_metacells_linear_fraction_todox!
export compute_blocks_modules_neighborhood_metacells_linear_fraction!
export compute_blocks_modules_neighborhood_metacells_log_linear_fraction!
export compute_metacells_modules_variance_over_mean!
export compute_metacells_modules_significance!

using Base.Threads
using DataAxesFormats
using Distances
using Random
using StatsBase
using TanayLabUtilities

using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_module_is_marker_in_environment_matrix
import Metacells.Contracts.block_module_is_strong_matrix
import Metacells.Contracts.block_module_n_genes_matrix
import Metacells.Contracts.block_module_n_skeletons_matrix
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_module_variance_over_mean_matrix
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.metacell_total_UMIs_vector
import Metacells.Contracts.module_axis

"""
    function compute_blocks_modules_is_strong(
        daf::DafWriter;
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        min_strong_UMIs::Integer = 8,
        min_strong_cells::Integer = $(DEFAULT.min_strong_cells),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

A mask of the strong modules that have enough UMIs in enough cells. Enough UMIs means that, if we downsample the cells
to a common total UMIs using `downsamples` with the `min_downsamples`, `min_downsamples_quantile`,
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
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_module_is_strong_matrix(GuaranteedOutput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
    ],
) function compute_blocks_modules_is_strong!(  # UNTESTED
    daf::DafWriter;
    min_strong_fraction::AbstractFloat = 2e-3,
    min_strong_UMIs::Integer = 7,
    min_strong_cells::Integer = 12,
    overwrite::Bool = false,
)::Nothing
    @assert min_strong_fraction >= 0
    @assert min_strong_UMIs >= 0
    @assert min_strong_cells >= 0

    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")

    name_per_block = axis_vector(daf, "block")
    name_per_module = axis_vector(daf, "module")

    is_found_per_module_per_block = get_matrix(daf, "module", "block", "is_found").array
    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module").array
    is_strong_per_module_per_block = zeros(Bool, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        compute_block_modules_is_strong!(
            daf;
            block_index,
            block_name,
            name_per_module,
            min_strong_fraction,
            min_strong_UMIs,
            min_strong_cells,
            module_per_gene_per_block,
            is_found_per_module_per_block,
            is_strong_per_module_per_block,
        )
    end

    set_matrix!(daf, "module", "block", "is_strong", bestify(is_strong_per_module_per_block); overwrite)
    return nothing
end

function compute_block_modules_is_strong!(  # UNTESTED
    daf::DafWriter;
    block_index::Integer,
    block_name::AbstractString,
    name_per_module::AbstractVector{<:AbstractString},
    min_strong_fraction::AbstractFloat,
    min_strong_UMIs::Integer,
    min_strong_cells::Integer,
    module_per_gene_per_block::AbstractMatrix{<:AbstractString},
    is_found_per_module_per_block::AbstractMatrix{Bool},
    is_strong_per_module_per_block::AbstractMatrix{Bool},
)::Nothing
    total_UMIs_per_neighborhood_cell =
        daf["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : total_UMIs"].array
    n_neighborhood_cells = length(total_UMIs_per_neighborhood_cell)

    UMIs_per_gene_per_neighborhood_cell =
        daf["/ gene / cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : UMIs"].array

    @views module_per_gene_of_block = module_per_gene_per_block[:, block_index]
    @views is_found_per_module_of_block = is_found_per_module_per_block[:, block_index]
    @views is_strong_per_module_of_block = is_strong_per_module_per_block[:, block_index]

    n_modules = length(name_per_module)

    strong_UMIs_per_module = Vector{UInt32}(undef, n_modules)
    strong_fraction_per_module = Vector{Float32}(undef, n_modules)
    for module_index in 1:n_modules
        if !is_found_per_module_of_block[module_index]
            strong_UMIs_per_module[module_index] = 0
            strong_fraction_per_module[module_index] = 0
        else
            module_name = name_per_module[module_index]
            module_genes_mask = module_per_gene_of_block .== module_name
            @views UMIs_per_module_gene_per_neighborhood_cell =
                UMIs_per_gene_per_neighborhood_cell[module_genes_mask, :]
            total_module_UMIs_per_neighborhood_cell = vec(sum(UMIs_per_module_gene_per_neighborhood_cell; dims = 1))
            @assert_vector(total_module_UMIs_per_neighborhood_cell, n_neighborhood_cells)
            total_module_fraction_per_neighborhood_cell =
                total_module_UMIs_per_neighborhood_cell ./ total_UMIs_per_neighborhood_cell
            @assert n_neighborhood_cells > min_strong_cells
            strong_UMIs_per_module[module_index] = UInt32(
                round(
                    quantile(
                        total_module_UMIs_per_neighborhood_cell,
                        1 - (min_strong_cells - 1) / (n_neighborhood_cells - 1),
                    ),
                ),
            )
            strong_fraction_per_module[module_index] = quantile(
                total_module_fraction_per_neighborhood_cell,
                1 - (min_strong_cells - 1) / (n_neighborhood_cells - 1),
            )
            @debug "Block: $(block_name) Large: $(n_neighborhood_cells) Module: $(module_name) Strong: $(strong_fraction_per_module[module_index]) UMIs: $(strong_UMIs_per_module[module_index])"
        end
    end

    is_strong_per_module_of_block .=
        (strong_UMIs_per_module .>= min_strong_UMIs) .& (strong_fraction_per_module .>= min_strong_fraction)
    n_strong_modules = sum(is_strong_per_module_of_block)
    if n_strong_modules > 1
        @debug "Block: $(block_name) picked Strong modules: $(n_strong_modules)"
    else
        strong_score_per_module =
            strong_UMIs_per_module ./ Float32(min_strong_UMIs) .+ strong_fraction_per_module ./ min_strong_fraction
        is_score_per_module = strong_score_per_module .> 0
        n_score_modules = sum(is_score_per_module)
        if n_score_modules == 0
            @debug "TODOX A BLOCK $(block_name) HAS NO STRONG MODULES"
        end
        @assert n_score_modules > 0
        if n_score_modules == 1
            is_strong_per_module_of_block .= is_score_per_module
            n_strong_modules = sum(is_strong_per_module_of_block)
            @assert n_strong_modules == 1
        else
            strong_module_indices = partialsortperm(strong_score_per_module, 1:2; rev = true)
            is_strong_per_module_of_block .= false
            is_strong_per_module_of_block[strong_module_indices] .= true
            n_strong_modules = sum(is_strong_per_module_of_block)
            @assert n_strong_modules == 2
        end
        @assert all((is_found_per_module_of_block .& is_strong_per_module_of_block) .== is_strong_per_module_of_block)
        @debug "Block: $(block_name) forced Strong modules: $(n_strong_modules)"
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
    compute_metacells_modules_linear_fraction_todox!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_module_neighborhood_mean_linear_fraction_matrix(GuaranteedOutput),
        block_module_neighborhood_std_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_neighborhood_metacells_linear_fraction_todox!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(daf, "cell")
    n_metacells = axis_length(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module").array

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    mean_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    std_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)

    progress_counter = Atomic{Int}(0)  # TODOX Refactor this everywhere
    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        indices_of_neighborhood_cells =
            daf["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : index"].array
        n_neighborhood_cells = length(indices_of_neighborhood_cells)
        @views module_per_gene = module_per_gene_per_block[:, block_index]
        for module_index in 1:n_modules
            module_name = name_per_module[module_index]
            indices_of_module_genes = findall(module_per_gene .== module_name)
            if length(indices_of_module_genes) == 0
                mean_linear_fraction_per_module_per_block[module_index, block_index] = 0
                std_linear_fraction_per_module_per_block[module_index, block_index] = 0
            else
                module_UMIs_per_neighborhood_cell =
                    vec(sum(UMIs_per_cell_per_gene[indices_of_neighborhood_cells, indices_of_module_genes]; dims = 2))
                @assert_vector(module_UMIs_per_neighborhood_cell, n_neighborhood_cells)

                linear_fraction_per_neighborhood_cell =
                    module_UMIs_per_neighborhood_cell ./ total_UMIs_per_cell[indices_of_neighborhood_cells]
                mean_linear_fraction_per_module_per_block[module_index, block_index] =
                    mean(linear_fraction_per_neighborhood_cell)
                std_linear_fraction_per_module_per_block[module_index, block_index] =
                    std(linear_fraction_per_neighborhood_cell)
            end
        end
        counter = atomic_add!(progress_counter, 1)
        print("\r$(progress_counter[]) ($(percent(counter + 1, n_blocks))) ...")
    end

    set_matrix!(
        daf,
        "module",
        "block",
        "neighborhood_mean_linear_fraction",
        bestify(mean_linear_fraction_per_module_per_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "module",
        "block",
        "neighborhood_std_linear_fraction",
        bestify(std_linear_fraction_per_module_per_block);
        overwrite,
    )

    return nothing
end

"""
    compute_metacells_modules_linear_fraction!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        metacell_total_UMIs_vector(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_metacell_module_linear_fraction_tensor(GuaranteedOutput),
    ],
) function compute_blocks_modules_neighborhood_metacells_linear_fraction!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module").array

    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    @threads :greedy for block_index in 1:n_blocks
        linear_fraction_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        block_name = name_per_block[block_index]
        indices_of_neighborhood_metacells =
            daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array
        n_neighborhood_metacells = length(indices_of_neighborhood_metacells)
        @views module_per_gene = module_per_gene_per_block[:, block_index]
        n_found_modules = 0
        for module_index in 1:n_modules
            module_name = name_per_module[module_index]
            indices_of_module_genes = findall(module_per_gene .== module_name)
            if length(indices_of_module_genes) > 0
                n_found_modules += 1
                module_UMIs_per_neighborhood_metacell = vec(
                    sum(
                        UMIs_per_metacell_per_gene[indices_of_neighborhood_metacells, indices_of_module_genes];
                        dims = 2,
                    ),
                )
                @assert_vector(module_UMIs_per_neighborhood_metacell, n_neighborhood_metacells)

                linear_fraction_per_module_per_metacell[module_index, indices_of_neighborhood_metacells] .=
                    module_UMIs_per_neighborhood_metacell ./ total_UMIs_per_metacell[indices_of_neighborhood_metacells]
            end
        end
        @debug "Block: $(block_name) neighborhood metacells: $(n_neighborhood_metacells) found modules: $(n_found_modules)"
        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_linear_fraction",
            bestify(linear_fraction_per_module_per_metacell);
            overwrite,
        )
    end

    return nothing
end

"""
    compute_blocks_modules_neighborhood_metacells_log_linear_fraction!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        block_gene_module_matrix(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_metacell_module_linear_fraction_tensor(RequiredInput),
        block_metacell_module_log_linear_fraction_tensor(GuaranteedOutput),
    ],
) function compute_blocks_modules_neighborhood_metacells_log_linear_fraction!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    module_per_gene_per_block = get_matrix(daf, "gene", "block", "module").array
    is_found_per_module_per_block = get_matrix(daf, "module", "block", "is_found").array

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        linear_fraction_per_module_per_metacell = get_matrix(daf, "module", "metacell", "$(block_name)_linear_fraction")
        indices_of_neighborhood_metacells =
            daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array
        indices_of_found_modules = daf["/ module & is_found ; block = $(block_name) : index"].array

        @views linear_fraction_per_found_module_per_neighborhood_metacell =
            linear_fraction_per_module_per_metacell[indices_of_found_modules, indices_of_neighborhood_metacells]

        log_linear_fraction_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        log_linear_fraction_per_module_per_metacell[indices_of_found_modules, indices_of_neighborhood_metacells] .=
            log2.(linear_fraction_per_found_module_per_neighborhood_metacell .+ gene_fraction_regularization)

        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_log_linear_fraction",
            bestify(log_linear_fraction_per_module_per_metacell);
            overwrite,
        )
        @debug "Block: $(block_name) neighborhood metacells: $(length(indices_of_neighborhood_metacells)) found modules: $(length(indices_of_found_modules))"
    end

    return nothing
end

"""
    compute_metacells_modules_variance_over_mean!(
        daf::DafWriter;
        min_downsamples::Integer = $(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = $(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = $(DEFAULT.max_downsamples_quantile),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The variance over mean of the total downsampled UMIs of each found module in the neighborhood of each block. This
ignores cells with less than the target downsampled UMIs per cell, to make the estimate more robust.
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_excluded_vector(RequiredInput),
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_module_variance_over_mean_matrix(GuaranteedOutput),
    ],
) function compute_metacells_modules_variance_over_mean!(
    daf::DafWriter;
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = 0.25,
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    block_index_per_metacell = daf["/ metacell : block => index"].array
    module_per_included_gene_per_block = daf["/ gene &! is_excluded / block : module"].array

    variance_over_mean_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_metacells; rng, policy = :serial) do metacell_index, rng  # TODOX
        metacell_name = name_per_metacell[metacell_index]
        block_index = block_index_per_metacell[metacell_index]

        @views module_per_included_gene = module_per_included_gene_per_block[:, block_index]

        UMIs_per_included_gene_per_metacell_cell =
            daf["/ gene &! is_excluded / cell & metacell = $(metacell_name) : UMIs"].array
        total_UMIs_per_metacell_cell = daf["/ cell & metacell = $(metacell_name) : total_UMIs"].array

        downsamples_UMIs = downsamples(
            total_UMIs_per_metacell_cell;
            min_downsamples,
            min_downsamples_quantile,
            max_downsamples_quantile,
        )

        is_significant_per_metacell_cell = total_UMIs_per_metacell_cell .>= downsamples_UMIs
        n_significant_metacell_cells = sum(is_significant_per_metacell_cell)

        @views UMIs_per_included_per_significant_metacell_cell =
            UMIs_per_included_gene_per_metacell_cell[:, is_significant_per_metacell_cell]
        downsampled_UMIs_per_included_per_significant_metacell_cell =
            downsample(UMIs_per_included_per_significant_metacell_cell, downsamples_UMIs; dims = 2, rng)

        total_downsampled_UMIs_per_significant_metacell_cell =
            vec(sum(downsampled_UMIs_per_included_per_significant_metacell_cell; dims = 1))
        @assert_vector(total_downsampled_UMIs_per_significant_metacell_cell, n_significant_metacell_cells)
        @assert all(total_downsampled_UMIs_per_significant_metacell_cell .== downsamples_UMIs)

        for module_index in 1:n_modules
            module_name = name_per_module[module_index]
            is_of_module_per_included_gene = module_per_included_gene .== module_name
            if any(is_of_module_per_included_gene)
                @views downsampled_UMIs_per_module_gene_per_significant_metacell_cell =
                    downsampled_UMIs_per_included_per_significant_metacell_cell[is_of_module_per_included_gene, :]
                module_total_downsampled_UMIs_per_significant_metacell_cell =
                    vec(sum(downsampled_UMIs_per_module_gene_per_significant_metacell_cell; dims = 1))
                @assert_vector(
                    module_total_downsampled_UMIs_per_significant_metacell_cell,
                    n_significant_metacell_cells
                )
                mean_module_total_downsampled_UMIs = mean(module_total_downsampled_UMIs_per_significant_metacell_cell)  # NOLINT
                if mean_module_total_downsampled_UMIs == 0
                    variance_over_mean = 1
                else
                    variance_of_module_total_downsampled_UMIs =
                        var(module_total_downsampled_UMIs_per_significant_metacell_cell)  # NOLINT
                    variance_over_mean = variance_of_module_total_downsampled_UMIs ./ mean_module_total_downsampled_UMIs
                end
                variance_over_mean_per_module_per_metacell[module_index, metacell_index] = variance_over_mean
            end
        end

        counter = atomic_add!(progress_counter, 1)
        print("\r$(progress_counter[]) ($(percent(counter + 1, n_metacells))) ...")
        return nothing
    end

    set_matrix!(
        daf,
        "module",
        "metacell",
        "variance_over_mean",
        bestify(variance_over_mean_per_module_per_metacell);
        overwrite,
    )
    return nothing
end

"""
    compute_metacells_modules_significance!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_mean_modules_distance_vector(GuaranteedOutput),
        metacell_std_modules_distance_vector(GuaranteedOutput),
    ],
) function compute_metacells_modules_significance!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_metacells = axis_length(daf, "metacell")
    name_per_metacell = axis_vector(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    block_index_per_metacell = daf["/ metacell : block => index"].array
    module_per_gene_per_block = daf["/ gene / block : module"].array

    mean_distance_per_metacell = zeros(Float32, n_metacells)
    std_distance_per_meyacell = zeros(Float32, n_metacells)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    progress_counter = Atomic{Int}(0)
    @threads :greedy for metacell_index in 1:n_metacells
        metacell_name = name_per_metacell[metacell_index]
        block_index = block_index_per_metacell[metacell_index]

        indices_of_metacell_cells = daf["/ cell & metacell = $(metacell_name) : index"].array
        n_metacell_cells = length(indices_of_metacell_cells)
        @assert n_metacell_cells > 0
        total_UMIs_per_metacell_cell = total_UMIs_per_cell[indices_of_metacell_cells]
        total_UMIs_of_metacell = sum(total_UMIs_per_metacell_cell)

        fraction_per_module_per_cell = zeros(Float32, n_modules, n_metacell_cells)
        mean_fraction_per_module = zeros(Float32, n_modules)

        @views module_per_gene = module_per_gene_per_block[:, block_index]
        for module_index in 1:n_modules
            module_name = name_per_module[module_index]
            indices_of_module_genes = findall(module_per_gene .== module_name)
            n_module_genes = length(indices_of_module_genes)
            if n_module_genes > 0
                module_UMIs_per_metacell_cell =
                    vec(sum(UMIs_per_cell_per_gene[indices_of_metacell_cells, indices_of_module_genes]; dims = 2))
                @assert_vector(module_UMIs_per_metacell_cell, n_metacell_cells)
                fraction_per_module_per_cell[module_index, :] =
                    module_UMIs_per_metacell_cell ./ total_UMIs_per_metacell_cell
                module_total_UMIs_of_metacell = sum(module_UMIs_per_metacell_cell)
                mean_fraction_per_module[module_index] = module_total_UMIs_of_metacell ./ total_UMIs_of_metacell
            end
        end

        distance_per_cell = pairwise(Euclidean(), Ref(mean_fraction_per_module), eachcol(fraction_per_module_per_cell))
        mean_distance_per_metacell[metacell_index] = mean(distance_per_cell)
        std_distance_per_meyacell[metacell_index] = std(distance_per_cell)

        counter = atomic_add!(progress_counter, 1)
        print("\r$(progress_counter[]) ($(percent(counter + 1, n_metacells))) ...")
    end

    set_vector!(daf, "metacell", "mean_modules_distance", bestify(mean_distance_per_metacell); overwrite)
    set_vector!(daf, "metacell", "std_modules_distance", bestify(std_distance_per_meyacell); overwrite)
    return nothing
end

end  # module
