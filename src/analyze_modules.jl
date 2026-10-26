"""
Do simple gene module analysis.
"""
module AnalyzeModules

export compute_matrix_of_cells_dispersion_per_metacell_per_module!
export compute_matrix_of_n_genes_per_module_per_block!
export compute_projected_modules_z_score_cutoff
export compute_projected_modules_z_score_cutoff_per_metacell!
export compute_stats_of_euclidean_modules_cells_distance_per_metacell!
export compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!
export compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!
export compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!
export compute_vector_of_n_modules_per_block!

using Base.Threads
using DataAxesFormats
using Distances
using LoopVectorization
using Random
using SparseArrays
using StatsBase
using TanayLabUtilities

using ..AnalyzeCells
using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

import Base.Threads.maxthreadid
import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_cells_dispersion_per_metacell_per_module
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_environment_per_metacell_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_n_genes_per_module_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.tensor_of_mean_cells_z_score_per_block_per_module_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_base_outlier_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_mean_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_metacell_modules_z_score_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_projected_metacell_per_cell
import Metacells.Contracts.vector_of_projected_modules_z_score_cutoff_level_per_metacell
import Metacells.Contracts.vector_of_projected_modules_z_score_cutoff_n_cells_per_metacell
import Metacells.Contracts.vector_of_projected_modules_z_score_cutoff_per_metacell
import Metacells.Contracts.vector_of_projected_modules_z_score_per_cell
import Metacells.Contracts.vector_of_std_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell

"""
    function compute_vector_of_n_modules_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_modules_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput)],
    data = [matrix_of_is_found_per_module_per_block(RequiredInput), vector_of_n_modules_per_block(CreatedOutput)],
) function compute_vector_of_n_modules_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_modules_per_block = daf["@ module @ block :: is_found >- Sum"].array
    set_vector!(daf, "block", "n_modules", n_modules_per_block; overwrite)
    @debug "Mean found modules per block: $(mean(n_modules_per_block))" _group = :mcs_results  # NOLINT
    return nothing
end

"""
    function compute_matrix_of_n_genes_per_module_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_n_genes_per_module_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [matrix_of_module_per_gene_per_block(RequiredInput), matrix_of_n_genes_per_module_per_block(CreatedOutput)],
) function compute_matrix_of_n_genes_per_module_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_modules = axis_length(daf, "module")
    n_blocks = axis_length(daf, "block")

    n_genes_per_module_per_block = Matrix{UInt32}(undef, n_modules, n_blocks)
    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "n_genes_per_module_per_block"),
    ) do block_index
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        @views n_genes_per_module = n_genes_per_module_per_block[:, block_index]
        fill!(n_genes_per_module, 0)
        for module_index in module_index_per_gene
            if module_index > 0
                n_genes_per_module[module_index] += 1
            end
        end
        return nothing
    end

    set_matrix!(daf, "module", "block", "n_genes", bestify(n_genes_per_module_per_block); overwrite)
    @debug (
        "Mean genes per found module: $(mean(n_genes_per_module_per_block[n_genes_per_module_per_block .!= 0]))"  # NOLINT
    ) _group = :mcs_results

    return nothing
end

"""
    compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!(
        daf::DafWriter;
        std_UMIs_regularization::AbstractFloat = $(DEFAULT.std_UMIs_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`tensor_of_mean_cells_z_score_per_block_per_module_per_metacell`](@ref). The z-score of each cell is
computed the same way sharpening computes it. That is, it is relative to the mean and standard deviation of the block's
environment cells, regularized by `std_UMIs_regularization`. This regularization is given in UMIs, and is converted to a
linear fraction using the mean total UMIs of the cells in the block's neighborhood.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_is_found_per_module_per_block(RequiredInput),
        vector_of_n_modules_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_total_neighborhood_UMIs_per_block(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        tensor_of_mean_cells_z_score_per_block_per_module_per_metacell(CreatedOutput),
    ],
) function compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!(  # UNTESTED
    daf::DafWriter;
    std_UMIs_regularization::AbstractFloat = 2.0,
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_genes = axis_length(daf, "gene")
    name_per_block = axis_vector(daf, "block")

    block_index_per_metacell = daf["@ metacell : block : index"].array
    metacell_index_per_cell = daf["@ cell : metacell ?? 0 : index"].array
    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array
    is_found_per_module_per_block = plain_array(get_matrix(daf, "module", "block", "is_found"))
    n_modules_per_block = plain_array(get_vector(daf, "block", "n_modules"))
    mean_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(daf, "module", "block", "mean_linear_fraction_in_environment_cells"))
    std_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(daf, "module", "block", "std_linear_fraction_in_environment_cells"))
    is_in_neighborhood_per_other_block_per_base_block =
        plain_array(get_matrix(daf, "block", "block", "is_in_neighborhood"))
    sparse_UMIs_per_cell_per_gene = plain_array(get_matrix(daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_cell = plain_array(get_vector(daf, "cell", "total_UMIs"))
    std_fraction_regularization_per_block = compute_std_fraction_regularization_per_block(daf; std_UMIs_regularization)

    block_index_per_cell = [
        metacell_index == 0 ? 0 : block_index_per_metacell[metacell_index] for metacell_index in metacell_index_per_cell
    ]
    n_cells_per_block = zeros(Int, n_blocks)
    for block_index in block_index_per_cell
        if block_index > 0
            n_cells_per_block[block_index] += 1
        end
    end
    n_region_cells_per_block = vec(transpose(is_in_neighborhood_per_other_block_per_base_block) * n_cells_per_block)
    max_n_region_cells = maximum(n_region_cells_per_block)
    max_n_block_modules = maximum(n_modules_per_block)

    z_score_per_max_module_per_max_region_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_region_cells) for _ in 1:maxthreadid()]
    mean_cells_z_score_per_max_module_per_thread =
        [Vector{Float32}(undef, max_n_block_modules) for _ in 1:maxthreadid()]
    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]

    # Channel pool of per-region-cell Float32 accumulators for `compute_z_score_per_found_module_per_region_cell!`'s
    # nested parallel-over-found-modules loop, as in sharpening.
    z_score_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(z_score_accumulator_pool, Vector{Float32}(undef, max_n_region_cells))
    end

    parallel_loop_wo_rng(
        1:n_blocks;
        weights = n_region_cells_per_block,
        progress = DebugProgress(
            sum(n_region_cells_per_block);
            group = :mcs_loops,
            desc = "mean_cells_z_score_per_block_per_module_per_metacell",
        ),
    ) do block_index
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        indices_of_region_cells = findall(
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1)),
        )
        n_region_cells = length(indices_of_region_cells)
        @assert n_region_cells == n_region_cells_per_block[block_index]

        # Order the cells by their metacell, so the cells of each metacell are a contiguous range of columns below.
        indices_of_region_cells = indices_of_region_cells[sortperm(metacell_index_per_cell[indices_of_region_cells])]

        n_block_modules = n_modules_per_block[block_index]
        is_found_per_module = is_found_per_module_per_thread[threadid()]
        @views z_score_per_found_module_per_region_cell =
            z_score_per_max_module_per_max_region_cell_per_thread[threadid()][1:n_block_modules, 1:n_region_cells]

        setup_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell,
            is_found_per_module,
            block_index,
            indices_of_region_cells,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            sparse_UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            z_score_accumulator_pool,
            n_genes,
            n_modules,
        )
        @assert n_block_modules == sum(is_found_per_module)
        indices_of_found_modules = findall(is_found_per_module)

        mean_cells_z_score_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        @views mean_cells_z_score_per_found_module =
            mean_cells_z_score_per_max_module_per_thread[threadid()][1:n_block_modules]
        metacell_index_per_region_cell = metacell_index_per_cell[indices_of_region_cells]
        first_position = 1
        while first_position <= n_region_cells
            metacell_index = metacell_index_per_region_cell[first_position]
            last_position = first_position
            while last_position < n_region_cells && metacell_index_per_region_cell[last_position + 1] == metacell_index
                last_position += 1
            end

            n_metacell_cells = last_position - first_position + 1
            @check_turbo_matrix(z_score_per_found_module_per_region_cell)
            @check_turbo_vector(mean_cells_z_score_per_found_module)
            @turbo for found_module_position in 1:n_block_modules
                sum_of_z_scores = 0.0f0
                for position in first_position:last_position
                    sum_of_z_scores += z_score_per_found_module_per_region_cell[found_module_position, position]
                end
                mean_cells_z_score_per_found_module[found_module_position] = sum_of_z_scores / n_metacell_cells
            end
            mean_cells_z_score_per_module_per_metacell[indices_of_found_modules, metacell_index] .=
                mean_cells_z_score_per_found_module

            first_position = last_position + 1
        end

        block_name = name_per_block[block_index]
        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_mean_cells_z_score",
            bestify(mean_cells_z_score_per_module_per_metacell);
            overwrite,
        )
        return nothing
    end

    return nothing
end

"""
    compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block`](@ref) and
and set [`matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block(CreatedOutput),
        matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block(CreatedOutput),
    ],
) function compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")

    sparse_UMIs_per_cell_per_gene = plain_array(get_matrix(daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_cell = plain_array(densify(get_vector(daf, "cell", "total_UMIs")))

    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array

    is_in_neighborhood_per_other_block_per_base_block =
        plain_array(get_matrix(daf, "block", "block", "is_in_neighborhood"))
    block_index_per_cell = daf["@ cell : metacell ?? 0 : block : index"].array

    n_cells, n_genes = size(sparse_UMIs_per_cell_per_gene)
    mean_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    std_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    is_in_neighborhood_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]

    # Per-block work iterates the cells in the block's neighborhood; weight blocks heaviest-first by that count.
    n_cells_per_block_with_metacell = zeros(Int, n_blocks)
    for cell_index in 1:n_cells
        block_of_cell = block_index_per_cell[cell_index]
        if block_of_cell > 0
            n_cells_per_block_with_metacell[block_of_cell] += 1
        end
    end
    n_cells_in_neighborhood_per_block =
        vec(is_in_neighborhood_per_other_block_per_base_block' * n_cells_per_block_with_metacell)
    max_n_neighborhood_cells = maximum(n_cells_in_neighborhood_per_block)

    # Channel pool of per-neighborhood-cell Float32 accumulators for the nested parallel-over-modules loop. Each inner
    # task takes one for its module, then puts it back. Sized to `nthreads()` because each take/put is a single
    # non-yielding scope.
    module_UMIs_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(module_UMIs_accumulator_pool, Vector{Float32}(undef, max_n_neighborhood_cells))
    end

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(
            sum(n_cells_in_neighborhood_per_block);
            group = :mcs_loops,
            desc = "stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block",
        ),
        weights = n_cells_in_neighborhood_per_block,
    ) do block_index
        is_in_neighborhood_per_cell = is_in_neighborhood_per_cell_per_thread[threadid()]
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_cell .=
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1))
        indices_of_neighborhood_cells = findall(is_in_neighborhood_per_cell)
        n_neighborhood_cells = length(indices_of_neighborhood_cells)
        @assert n_neighborhood_cells > 0

        # Invert `module_index_per_gene` once per block: the per-module loop below reads the gene list directly.
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        gene_indices_per_module = [Int[] for _ in 1:n_modules]
        for gene_index in 1:n_genes
            module_index = module_index_per_gene[gene_index]
            if module_index > 0
                push!(gene_indices_per_module[module_index], gene_index)
            end
        end

        # Parallel over modules - each writes to its own (module_index, block_index) slot in the output matrices, so
        # the inner tasks are race-free. Per task: take an accumulator vector, gene-outer / cell-inner `@turbo` fill,
        # scalar reduction, write mean+std, return accumulator.
        parallel_loop_wo_rng(
            1:n_modules;
            nested = true,
            policy = :greedy,
            name = "stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block.modules",
        ) do module_index
            gene_indices_in_module = gene_indices_per_module[module_index]
            if isempty(gene_indices_in_module)
                mean_linear_fraction_per_module_per_block[module_index, block_index] = 0.0f0
                std_linear_fraction_per_module_per_block[module_index, block_index] = 0.0f0
                return nothing
            end
            module_UMIs_per_neighborhood_cell = take!(module_UMIs_accumulator_pool)
            try
                # The module's genes are a single group of columns, summed into the neighborhood cells.
                @views module_UMIs_per_neighborhood_cell_matrix =
                    reshape(module_UMIs_per_neighborhood_cell[1:n_neighborhood_cells], 1, n_neighborhood_cells)
                module_UMIs_per_neighborhood_cell_matrix .= 0.0f0
                sum_sparse_columns_per_group!(
                    module_UMIs_per_neighborhood_cell_matrix,
                    sparse_UMIs_per_cell_per_gene,
                    indices_of_neighborhood_cells,
                    view(gene_indices_per_module, module_index:module_index),
                )
                sum_fractions = 0.0
                sum_squared_fractions = 0.0
                @inbounds for neighborhood_cell_position in 1:n_neighborhood_cells
                    cell_index = indices_of_neighborhood_cells[neighborhood_cell_position]
                    fraction =
                        module_UMIs_per_neighborhood_cell[neighborhood_cell_position] / total_UMIs_per_cell[cell_index]
                    sum_fractions += fraction
                    sum_squared_fractions += fraction * fraction
                end
                mean_fraction = sum_fractions / n_neighborhood_cells
                mean_linear_fraction_per_module_per_block[module_index, block_index] = Float32(mean_fraction)
                std_linear_fraction_per_module_per_block[module_index, block_index] = if n_neighborhood_cells > 1
                    Float32(
                        sqrt(
                            max(
                                (sum_squared_fractions - sum_fractions * mean_fraction) / (n_neighborhood_cells - 1),
                                0.0,
                            ),
                        ),
                    )
                else
                    0.0f0
                end
            finally
                put!(module_UMIs_accumulator_pool, module_UMIs_per_neighborhood_cell)
            end
            return nothing
        end
        return nothing
    end

    set_matrix!(
        daf,
        "module",
        "block",
        "mean_linear_fraction_in_neighborhood_cells",
        bestify(mean_linear_fraction_per_module_per_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "module",
        "block",
        "std_linear_fraction_in_neighborhood_cells",
        bestify(std_linear_fraction_per_module_per_block);
        overwrite,
    )

    return nothing
end

"""
    compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block`](@ref) and
and set [`matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_is_in_environment_per_metacell_per_block(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(CreatedOutput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(CreatedOutput),
    ],
) function compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_metacells = axis_length(daf, "metacell")

    sparse_UMIs_per_cell_per_gene = plain_array(get_matrix(daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_cell = plain_array(densify(get_vector(daf, "cell", "total_UMIs")))

    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array

    is_in_environment_per_metacell_per_block = plain_array(get_matrix(daf, "metacell", "block", "is_in_environment"))
    metacell_index_per_cell = daf["@ cell : metacell ?? 0 : index"].array

    n_cells, n_genes = size(sparse_UMIs_per_cell_per_gene)
    mean_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    std_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    is_in_environment_per_metacell_per_thread = [BitVector(undef, n_metacells) for _ in 1:maxthreadid()]
    is_in_environment_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:maxthreadid()]

    # Per-block work iterates the cells in the block's environment; weight blocks heaviest-first by that count.
    n_cells_per_metacell = zeros(Int, n_metacells)
    for cell_index in 1:n_cells
        metacell_of_cell = metacell_index_per_cell[cell_index]
        if metacell_of_cell > 0
            n_cells_per_metacell[metacell_of_cell] += 1
        end
    end
    n_cells_in_environment_per_block = vec(is_in_environment_per_metacell_per_block' * n_cells_per_metacell)
    max_n_environment_cells = maximum(n_cells_in_environment_per_block)

    # Channel pool of per-environment-cell Float32 accumulators for the nested parallel-over-modules loop. Each inner
    # task takes one for its module, then puts it back. Sized to `nthreads()` because each take/put is a single
    # non-yielding scope.
    module_UMIs_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(module_UMIs_accumulator_pool, Vector{Float32}(undef, max_n_environment_cells))
    end

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(
            sum(n_cells_in_environment_per_block);
            group = :mcs_loops,
            desc = "stats_of_linear_fraction_in_environment_cells_per_module_per_block",
        ),
        weights = n_cells_in_environment_per_block,
    ) do block_index
        is_in_environment_per_metacell = is_in_environment_per_metacell_per_thread[threadid()]
        @views is_in_environment_per_metacell .= is_in_environment_per_metacell_per_block[:, block_index]
        is_in_environment_per_cell = is_in_environment_per_cell_per_thread[threadid()]
        is_in_environment_per_cell .=
            (metacell_index_per_cell .> 0) .&
            getindex.(Ref(is_in_environment_per_metacell), max.(metacell_index_per_cell, 1))
        indices_of_environment_cells = findall(is_in_environment_per_cell)
        n_environment_cells = length(indices_of_environment_cells)
        @assert n_environment_cells > 0

        # Invert `module_index_per_gene` once per block: the per-module loop below reads the gene list directly.
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        gene_indices_per_module = [Int[] for _ in 1:n_modules]
        for gene_index in 1:n_genes
            module_index = module_index_per_gene[gene_index]
            if module_index > 0
                push!(gene_indices_per_module[module_index], gene_index)
            end
        end

        # Parallel over modules - each writes to its own (module_index, block_index) slot in the output matrices, so
        # the inner tasks are race-free. Per task: take an accumulator vector, gene-outer / cell-inner `@turbo` fill,
        # scalar reduction, write mean+std, return accumulator.
        parallel_loop_wo_rng(
            1:n_modules;
            nested = true,
            policy = :greedy,
            name = "stats_of_linear_fraction_in_environment_cells_per_module_per_block.modules",
        ) do module_index
            gene_indices_in_module = gene_indices_per_module[module_index]
            if isempty(gene_indices_in_module)
                mean_linear_fraction_per_module_per_block[module_index, block_index] = 0.0f0
                std_linear_fraction_per_module_per_block[module_index, block_index] = 0.0f0
                return nothing
            end
            module_UMIs_per_environment_cell = take!(module_UMIs_accumulator_pool)
            try
                # The module's genes are a single group of columns, summed into the environment cells.
                @views module_UMIs_per_environment_cell_matrix =
                    reshape(module_UMIs_per_environment_cell[1:n_environment_cells], 1, n_environment_cells)
                module_UMIs_per_environment_cell_matrix .= 0.0f0
                sum_sparse_columns_per_group!(
                    module_UMIs_per_environment_cell_matrix,
                    sparse_UMIs_per_cell_per_gene,
                    indices_of_environment_cells,
                    view(gene_indices_per_module, module_index:module_index),
                )
                sum_fractions = 0.0
                sum_squared_fractions = 0.0
                @inbounds for environment_cell_position in 1:n_environment_cells
                    cell_index = indices_of_environment_cells[environment_cell_position]
                    fraction =
                        module_UMIs_per_environment_cell[environment_cell_position] / total_UMIs_per_cell[cell_index]
                    sum_fractions += fraction
                    sum_squared_fractions += fraction * fraction
                end
                mean_fraction = sum_fractions / n_environment_cells
                mean_linear_fraction_per_module_per_block[module_index, block_index] = Float32(mean_fraction)
                std_linear_fraction_per_module_per_block[module_index, block_index] = if n_environment_cells > 1
                    Float32(
                        sqrt(
                            max(
                                (sum_squared_fractions - sum_fractions * mean_fraction) / (n_environment_cells - 1),
                                0.0,
                            ),
                        ),
                    )
                else
                    0.0f0
                end
            finally
                put!(module_UMIs_accumulator_pool, module_UMIs_per_environment_cell)
            end
            return nothing
        end
        return nothing
    end

    set_matrix!(
        daf,
        "module",
        "block",
        "mean_linear_fraction_in_environment_cells",
        bestify(mean_linear_fraction_per_module_per_block);
        overwrite,
    )
    set_matrix!(
        daf,
        "module",
        "block",
        "std_linear_fraction_in_environment_cells",
        bestify(std_linear_fraction_per_module_per_block);
        overwrite,
    )

    return nothing
end

"""
    compute_stats_of_euclidean_modules_cells_distance_per_metacell!(
        daf::DafWriter;
        std_UMIs_regularization::AbstractFloat = $(DEFAULT.std_UMIs_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_mean_euclidean_modules_cells_distance_per_metacell`](@ref),
[`vector_of_std_euclidean_modules_cells_distance_per_metacell`](@ref) and
[`vector_of_metacell_modules_z_score_per_cell`](@ref). The distance between a cell and its metacell is measured in the
space sharpening clusters cells in. That is, it is the Euclidean distance between the cell's z-scores of the found
modules of the metacell's block and the metacell's
[`tensor_of_mean_cells_z_score_per_block_per_module_per_metacell`](@ref). The z-scores are computed as in
[`compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!`](@ref), using the same
`std_UMIs_regularization`. The z score of each cell is its distance minus the mean distance of its metacell's cells,
divided by their standard deviation. This is computed the same way
[`compute_cells_projection!`](@ref Metacells.ProjectCells.compute_cells_projection!) computes the z score of a query
cell relative to its projected metacell.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_is_found_per_module_per_block(RequiredInput),
        vector_of_n_modules_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        vector_of_total_neighborhood_UMIs_per_block(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        tensor_of_mean_cells_z_score_per_block_per_module_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_mean_euclidean_modules_cells_distance_per_metacell(CreatedOutput),
        vector_of_std_euclidean_modules_cells_distance_per_metacell(CreatedOutput),
        vector_of_metacell_modules_z_score_per_cell(CreatedOutput),
    ],
) function compute_stats_of_euclidean_modules_cells_distance_per_metacell!(  # UNTESTED
    daf::DafWriter;
    std_UMIs_regularization::AbstractFloat = function_default(
        compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!,
        :std_UMIs_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(daf, "cell")
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    n_genes = axis_length(daf, "gene")
    name_per_block = axis_vector(daf, "block")

    block_index_per_metacell = daf["@ metacell : block : index"].array
    metacell_index_per_cell = daf["@ cell : metacell ?? 0 : index"].array
    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array
    is_found_per_module_per_block = plain_array(get_matrix(daf, "module", "block", "is_found"))
    n_modules_per_block = plain_array(get_vector(daf, "block", "n_modules"))
    mean_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(daf, "module", "block", "mean_linear_fraction_in_environment_cells"))
    std_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(daf, "module", "block", "std_linear_fraction_in_environment_cells"))
    sparse_UMIs_per_cell_per_gene = plain_array(get_matrix(daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_cell = plain_array(get_vector(daf, "cell", "total_UMIs"))
    std_fraction_regularization_per_block = compute_std_fraction_regularization_per_block(daf; std_UMIs_regularization)

    block_index_per_cell = [
        metacell_index == 0 ? 0 : block_index_per_metacell[metacell_index] for metacell_index in metacell_index_per_cell
    ]
    n_cells_per_block = zeros(Int, n_blocks)
    for block_index in block_index_per_cell
        if block_index > 0
            n_cells_per_block[block_index] += 1
        end
    end
    max_n_block_cells = maximum(n_cells_per_block)
    max_n_block_modules = maximum(n_modules_per_block)

    mean_distance_per_metacell = zeros(Float32, n_metacells)
    std_distance_per_metacell = zeros(Float32, n_metacells)
    z_score_per_cell = fill(Inf32, n_cells)

    z_score_per_max_module_per_max_block_cell_per_thread =
        [Matrix{Float32}(undef, max_n_block_modules, max_n_block_cells) for _ in 1:maxthreadid()]
    distance_per_max_block_cell_per_thread = [Vector{Float32}(undef, max_n_block_cells) for _ in 1:maxthreadid()]
    is_found_per_module_per_thread = [BitVector(undef, n_modules) for _ in 1:maxthreadid()]

    # Channel pool of per-block-cell Float32 accumulators for `compute_z_score_per_found_module_per_region_cell!`'s
    # nested parallel-over-found-modules loop, as in sharpening.
    z_score_accumulator_pool = Channel{Vector{Float32}}(nthreads())
    for _ in 1:nthreads()
        put!(z_score_accumulator_pool, Vector{Float32}(undef, max_n_block_cells))
    end

    parallel_loop_wo_rng(
        1:n_blocks;
        weights = n_cells_per_block,
        progress = DebugProgress(
            sum(n_cells_per_block);
            group = :mcs_loops,
            desc = "stats_of_euclidean_modules_cells_distance_per_metacell",
        ),
    ) do block_index
        indices_of_block_cells = findall(block_index_per_cell .== block_index)
        n_block_cells = length(indices_of_block_cells)
        if n_block_cells == 0
            return nothing
        end

        # Order the cells by their metacell, so the cells of each metacell are a contiguous range of columns below.
        indices_of_block_cells = indices_of_block_cells[sortperm(metacell_index_per_cell[indices_of_block_cells])]

        n_block_modules = n_modules_per_block[block_index]
        is_found_per_module = is_found_per_module_per_thread[threadid()]
        @views z_score_per_found_module_per_block_cell =
            z_score_per_max_module_per_max_block_cell_per_thread[threadid()][1:n_block_modules, 1:n_block_cells]

        setup_z_score_per_found_module_per_region_cell!(;
            z_score_per_found_module_per_region_cell = z_score_per_found_module_per_block_cell,
            is_found_per_module,
            block_index,
            indices_of_region_cells = indices_of_block_cells,
            is_found_per_module_per_block,
            module_index_per_gene_per_block,
            mean_linear_fraction_in_environment_cells_per_module_per_block,
            std_linear_fraction_in_environment_cells_per_module_per_block,
            std_fraction_regularization_per_block,
            sparse_UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            z_score_accumulator_pool,
            n_genes,
            n_modules,
        )
        @assert n_block_modules == sum(is_found_per_module)
        indices_of_found_modules = findall(is_found_per_module)

        block_name = name_per_block[block_index]
        mean_cells_z_score_per_module_per_metacell =
            plain_array(get_matrix(daf, "module", "metacell", "$(block_name)_mean_cells_z_score"))

        distance_per_block_cell = distance_per_max_block_cell_per_thread[threadid()]
        metacell_index_per_block_cell = metacell_index_per_cell[indices_of_block_cells]
        first_position = 1
        while first_position <= n_block_cells
            metacell_index = metacell_index_per_block_cell[first_position]
            last_position = first_position
            while last_position < n_block_cells && metacell_index_per_block_cell[last_position + 1] == metacell_index
                last_position += 1
            end

            mean_cells_z_score_per_found_module =
                densify(mean_cells_z_score_per_module_per_metacell[indices_of_found_modules, metacell_index])
            @check_turbo_matrix(z_score_per_found_module_per_block_cell)
            @check_turbo_vector(mean_cells_z_score_per_found_module)
            @check_turbo_vector(distance_per_block_cell)
            @turbo for position in first_position:last_position
                squared_distance = 0.0f0
                for found_module_position in 1:n_block_modules
                    delta =
                        z_score_per_found_module_per_block_cell[found_module_position, position] -
                        mean_cells_z_score_per_found_module[found_module_position]
                    squared_distance += delta * delta
                end
                distance_per_block_cell[position] = sqrt(squared_distance)
            end

            @views distance_per_metacell_cell = distance_per_block_cell[first_position:last_position]
            mean_distance = mean(distance_per_metacell_cell)  # NOLINT
            std_distance = std(distance_per_metacell_cell)  # NOLINT
            mean_distance_per_metacell[metacell_index] = mean_distance
            std_distance_per_metacell[metacell_index] = std_distance
            z_score_per_cell[indices_of_block_cells[first_position:last_position]] .=
                (distance_per_metacell_cell .- mean_distance) ./ std_distance

            first_position = last_position + 1
        end

        return nothing
    end

    set_vector!(
        daf,
        "metacell",
        "mean_euclidean_modules_cells_distance",
        bestify(mean_distance_per_metacell);
        overwrite,
    )
    set_vector!(daf, "metacell", "std_euclidean_modules_cells_distance", bestify(std_distance_per_metacell); overwrite)
    set_vector!(daf, "cell", "metacell_modules_z_score", z_score_per_cell; overwrite)
    return nothing
end

"""
    compute_projected_modules_z_score_cutoff(
        daf::DafReader;
        rejected_fraction::AbstractFloat = $(DEFAULT.rejected_fraction),
    )::AbstractFloat

Compute a cutoff for the [`vector_of_projected_modules_z_score_per_cell`](@ref) of query cells projected onto this
(atlas) `daf`. This is the z score which only `rejected_fraction` of the atlas cells exceed, using their
[`vector_of_metacell_modules_z_score_per_cell`](@ref). That is, a query cell above the cutoff fits its projected
metacell worse than all but `rejected_fraction` of the atlas cells fit their own metacell. Atlas cells which do not
belong to any metacell are ignored.

This makes no assumption about the distribution of the z scores, which is not normal. It also makes no assumption about
the query cells. It only describes how well the atlas cells fit their own metacells.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [cell_axis(RequiredInput)],
    data = [vector_of_metacell_modules_z_score_per_cell(RequiredInput)],
) function compute_projected_modules_z_score_cutoff(  # UNTESTED
    daf::DafReader;
    rejected_fraction::AbstractFloat = 0.01,
)::AbstractFloat
    @assert 0 < rejected_fraction < 1
    z_score_per_cell = plain_array(get_vector(daf, "cell", "metacell_modules_z_score"))
    finite_z_scores = filter(isfinite, z_score_per_cell)
    @assert !isempty(finite_z_scores)
    return quantile(finite_z_scores, 1 - rejected_fraction)  # NOLINT
end

"""
    compute_projected_modules_z_score_cutoff_per_metacell!(;
        atlas_daf::DafWriter,
        calibration_daf::DafReader,
        rejected_fraction::AbstractFloat = $(DEFAULT.rejected_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_projected_modules_z_score_cutoff_per_metacell`](@ref),
[`vector_of_projected_modules_z_score_cutoff_level_per_metacell`](@ref) and
[`vector_of_projected_modules_z_score_cutoff_n_cells_per_metacell`](@ref) of the `atlas_daf`, from the cells of the
`calibration_daf` projected onto it (see
[`compute_cells_projection!`](@ref Metacells.ProjectCells.compute_cells_projection!)). The calibration cells are those
which are neither excluded nor base outliers. For a data set too large to compute its metacells directly, these are the
cells of the sample the atlas was computed from. Otherwise, they are the atlas's own cells.

A metacell's cutoff is the quantile at position `q * (n + 1)` (where `q = 1 - rejected_fraction`), interpolated, of the
projected modules z scores of the `n` calibration cells projected onto it. Query cells from the same distribution as the
calibration cells exceed it at the `rejected_fraction` on average, whatever `n` is, as long as `q * (n + 1) <= n`.

A metacell with fewer calibration cells than that takes the cutoff of the calibration cells projected onto its block,
computed the same way. A block with too few takes the cutoff of its neighborhood, and a neighborhood with too few takes
the cutoff of all the calibration cells. A metacell taking the cutoff of all the calibration cells is a strong signal
that something is wrong with the data set or the parameters, since its whole neighborhood got too few calibration cells.

# Atlas

$(CONTRACT1)

# Calibration

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "atlas_daf",
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        vector_of_projected_modules_z_score_cutoff_per_metacell(CreatedOutput),
        vector_of_projected_modules_z_score_cutoff_level_per_metacell(CreatedOutput),
        vector_of_projected_modules_z_score_cutoff_n_cells_per_metacell(CreatedOutput),
    ],
) Contract(;
    name = "calibration_daf",
    axes = [cell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_cell(RequiredInput),
        vector_of_is_base_outlier_per_cell(RequiredInput),
        vector_of_projected_metacell_per_cell(RequiredInput),
        vector_of_projected_modules_z_score_per_cell(RequiredInput),
    ],
) function compute_projected_modules_z_score_cutoff_per_metacell!(;  # UNTESTED
    atlas_daf::DafWriter,
    calibration_daf::DafReader,
    rejected_fraction::AbstractFloat = function_default(compute_projected_modules_z_score_cutoff, :rejected_fraction),
    overwrite::Bool = false,
)::Nothing
    @assert 0 < rejected_fraction < 1
    quantile_fraction = 1 - rejected_fraction

    n_metacells = axis_length(atlas_daf, "metacell")
    n_blocks = axis_length(atlas_daf, "block")
    block_index_per_metacell = atlas_daf["@ metacell : block : index"].array
    is_in_neighborhood_per_other_block_per_base_block =
        plain_array(get_matrix(atlas_daf, "block", "block", "is_in_neighborhood"))

    is_excluded_per_cell = plain_array(get_vector(calibration_daf, "cell", "is_excluded"))
    is_base_outlier_per_cell = plain_array(get_vector(calibration_daf, "cell", "is_base_outlier"))
    is_calibration_per_cell = .!is_excluded_per_cell .& .!is_base_outlier_per_cell
    metacell_index_per_cell = axis_indices(
        atlas_daf,
        "metacell",
        plain_array(get_vector(calibration_daf, "cell", "projected_metacell"));
        allow_empty = true,
    )
    z_score_per_cell = plain_array(get_vector(calibration_daf, "cell", "projected_metacell_modules_z_score"))

    z_scores_per_metacell = [Float32[] for _ in 1:n_metacells]
    for cell_index in findall(is_calibration_per_cell)
        metacell_index = metacell_index_per_cell[cell_index]
        @assert metacell_index > 0
        push!(z_scores_per_metacell[metacell_index], z_score_per_cell[cell_index])
    end
    z_scores_per_block = [
        reduce(vcat, z_scores_per_metacell[block_index_per_metacell .== block_index]; init = Float32[]) for
        block_index in 1:n_blocks
    ]

    all_z_scores = reduce(vcat, z_scores_per_block; init = Float32[])
    @assert !isempty(all_z_scores)
    all_cutoff = z_score_cutoff(all_z_scores, quantile_fraction)

    cutoff_per_block = Vector{Float32}(undef, n_blocks)
    level_per_block = Vector{String}(undef, n_blocks)
    n_cells_per_block = Vector{UInt32}(undef, n_blocks)
    for block_index in 1:n_blocks
        block_z_scores = z_scores_per_block[block_index]
        if is_enough_for_z_score_cutoff(block_z_scores, quantile_fraction)
            cutoff_per_block[block_index] = z_score_cutoff(block_z_scores, quantile_fraction)
            level_per_block[block_index] = "block"
            n_cells_per_block[block_index] = length(block_z_scores)
        else
            neighborhood_z_scores = reduce(
                vcat,
                z_scores_per_block[findall(view(is_in_neighborhood_per_other_block_per_base_block, :, block_index))];
                init = Float32[],
            )
            if is_enough_for_z_score_cutoff(neighborhood_z_scores, quantile_fraction)
                cutoff_per_block[block_index] = z_score_cutoff(neighborhood_z_scores, quantile_fraction)
                level_per_block[block_index] = "neighborhood"
                n_cells_per_block[block_index] = length(neighborhood_z_scores)
            else
                cutoff_per_block[block_index] = all_cutoff
                level_per_block[block_index] = "all"
                n_cells_per_block[block_index] = length(all_z_scores)
            end
        end
    end

    cutoff_per_metacell = Vector{Float32}(undef, n_metacells)
    level_per_metacell = Vector{String}(undef, n_metacells)
    n_cells_per_metacell = Vector{UInt32}(undef, n_metacells)
    for metacell_index in 1:n_metacells
        metacell_z_scores = z_scores_per_metacell[metacell_index]
        if is_enough_for_z_score_cutoff(metacell_z_scores, quantile_fraction)
            cutoff_per_metacell[metacell_index] = z_score_cutoff(metacell_z_scores, quantile_fraction)
            level_per_metacell[metacell_index] = "metacell"
            n_cells_per_metacell[metacell_index] = length(metacell_z_scores)
        else
            block_index = block_index_per_metacell[metacell_index]
            cutoff_per_metacell[metacell_index] = cutoff_per_block[block_index]
            level_per_metacell[metacell_index] = level_per_block[block_index]
            n_cells_per_metacell[metacell_index] = n_cells_per_block[block_index]
        end
    end

    for level in ("metacell", "block", "neighborhood", "all")
        @debug "Metacells with a $(level) cutoff: $(count(level_per_metacell .== level))" _group = :mcs_results
    end
    n_all_level_metacells = count(level_per_metacell .== "all")
    if n_all_level_metacells > 0
        @warn (
            "$(n_all_level_metacells) metacells take the cutoff of all the calibration cells, since their" *
            " neighborhood got too few of them. Something is likely wrong with the data set or the parameters."
        )
    end

    set_vector!(atlas_daf, "metacell", "projected_modules_z_score_cutoff", cutoff_per_metacell; overwrite)
    set_vector!(atlas_daf, "metacell", "projected_modules_z_score_cutoff_level", level_per_metacell; overwrite)
    set_vector!(atlas_daf, "metacell", "projected_modules_z_score_cutoff_n_cells", n_cells_per_metacell; overwrite)
    return nothing
end

# Whether there are enough z scores for their quantile at position `quantile_fraction * (n + 1)` to be at most `n`.
function is_enough_for_z_score_cutoff(z_scores::AbstractVector{Float32}, quantile_fraction::AbstractFloat)::Bool  # UNTESTED
    return quantile_fraction * (length(z_scores) + 1) <= length(z_scores)
end

# The quantile of the z scores at position `quantile_fraction * (n + 1)`, interpolated. New z scores from the same
# distribution exceed it at `1 - quantile_fraction` on average.
function z_score_cutoff(z_scores::AbstractVector{Float32}, quantile_fraction::AbstractFloat)::Float32  # UNTESTED
    return Float32(quantile(z_scores, quantile_fraction; alpha = 0, beta = 0))  # NOLINT
end

# Compute and return the maximum (across the block's found modules) cells_dispersion - the ratio between the actual
# standard deviation of normalized module UMIs in the given group of cells (e.g. a metacell or a candidate cluster)
# and the standard deviation expected from pure multinomial sampling noise. Also fills the per-module values into the
# `cells_dispersion_per_module` scratch (zero for non-found modules and for modules that don't pass the effective
# threshold). If no found module's mean reaches `min_module_UMIs`, the threshold is auto-reduced by integer steps until
# at least one module qualifies (concretely: `effective_threshold = min(min_module_UMIs, floor(max_mean))`). Returns 0
# when there are fewer than two cells, no found modules, or every found module has mean = 0. If
# `module_UMIs_per_found_module_per_region_cell` is given, the module UMIs of the cells are read from it, at the
# `region_positions_of_cells`, instead of being summed from the `sparse_UMIs_per_cell_per_gene`.
function maximal_cells_dispersion_of_modules!(;  # UNTESTED
    cells_dispersion_per_module::AbstractVector{<:AbstractFloat},
    mean_normalized_per_module::AbstractVector{<:AbstractFloat},
    indices_of_cells::AbstractVector{<:Integer},
    module_UMIs_per_found_module_per_region_cell::Maybe{AbstractMatrix{Float32}} = nothing,
    region_positions_of_cells::Maybe{AbstractVector{<:Integer}} = nothing,
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    is_found_per_module::Union{AbstractVector{Bool}, BitVector},
    gene_indices_per_module::AbstractVector{<:AbstractVector{<:Integer}},
    total_UMIs_per_max_cells::AbstractVector{<:AbstractFloat},
    normalized_factor_per_max_cells::AbstractVector{<:AbstractFloat},
    normalized_UMIs_quantile::AbstractFloat,
    min_module_UMIs::Integer,
)::AbstractFloat
    @assert any(is_found_per_module)
    cells_dispersion_per_module .= 0
    mean_normalized_per_module .= 0
    n_cells = length(indices_of_cells)
    if n_cells < 2
        return 0.0f0
    end

    @views total_UMIs_per_cells = total_UMIs_per_max_cells[1:n_cells]
    @views normalized_factor_per_cells = normalized_factor_per_max_cells[1:n_cells]

    normalized_factor_per_cells .= getindex.(Ref(total_UMIs_per_cell), indices_of_cells)
    total_UMIs_per_cells .= normalized_factor_per_cells
    normalized_total_UMIs = quantile!(total_UMIs_per_cells, normalized_UMIs_quantile)
    @. normalized_factor_per_cells = normalized_total_UMIs / normalized_factor_per_cells

    # Single pass: per found module with a positive mean, compute mean + dispersion in the same per-cell loop. Store
    # both in scratch. Track max mean so the effective threshold can be derived from the data after the loop instead of
    # iterating thresholds. The module UMIs of each cell are summed into the total UMIs buffer, whose content
    # `quantile!` has already reordered.
    module_UMIs_per_cells = reshape(total_UMIs_per_cells, 1, n_cells)
    n_modules = length(is_found_per_module)
    max_mean_normalized_module_UMIs = 0.0
    found_module_position = 0
    for module_index in 1:n_modules
        if !is_found_per_module[module_index]
            continue
        end
        found_module_position += 1

        if module_UMIs_per_found_module_per_region_cell === nothing
            module_UMIs_per_cells .= 0
            sum_sparse_columns_per_group!(
                module_UMIs_per_cells,
                sparse_UMIs_per_cell_per_gene,
                indices_of_cells,
                view(gene_indices_per_module, module_index:module_index),
            )
        else
            @assert region_positions_of_cells !== nothing
            for cell_position in 1:n_cells
                module_UMIs_per_cells[cell_position] = module_UMIs_per_found_module_per_region_cell[
                    found_module_position,
                    region_positions_of_cells[cell_position],
                ]
            end
        end

        sum_normalized_module_UMIs = 0.0
        sum_squared_normalized_module_UMIs = 0.0
        for cell_position in 1:n_cells
            normalized_module_UMIs = module_UMIs_per_cells[cell_position] * normalized_factor_per_cells[cell_position]
            sum_normalized_module_UMIs += normalized_module_UMIs
            sum_squared_normalized_module_UMIs += normalized_module_UMIs * normalized_module_UMIs
        end

        mean_normalized_module_UMIs = sum_normalized_module_UMIs / n_cells
        if mean_normalized_module_UMIs <= 0
            continue
        end
        mean_normalized_per_module[module_index] = mean_normalized_module_UMIs

        actual_var_normalized_module_UMIs = max(
            (sum_squared_normalized_module_UMIs - sum_normalized_module_UMIs * mean_normalized_module_UMIs) /
            (n_cells - 1),
            0.0,
        )
        cells_dispersion_per_module[module_index] =
            Float32(sqrt(actual_var_normalized_module_UMIs / mean_normalized_module_UMIs))
        max_mean_normalized_module_UMIs = max(max_mean_normalized_module_UMIs, mean_normalized_module_UMIs)
    end

    if max_mean_normalized_module_UMIs <= 0
        return 0.0f0
    end

    effective_threshold = min(Float64(min_module_UMIs), floor(max_mean_normalized_module_UMIs))

    max_cells_dispersion = 0.0f0
    for module_index in 1:n_modules
        if !is_found_per_module[module_index]
            continue
        end
        if mean_normalized_per_module[module_index] < effective_threshold
            cells_dispersion_per_module[module_index] = 0
        else
            max_cells_dispersion = max(max_cells_dispersion, cells_dispersion_per_module[module_index])
        end
    end

    return max_cells_dispersion
end

"""
    compute_matrix_of_cells_dispersion_per_metacell_per_module!(
        daf::DafWriter;
        normalized_UMIs_quantile::AbstractFloat = $(DEFAULT.normalized_UMIs_quantile),
        min_module_UMIs::Integer = $(DEFAULT.min_module_UMIs),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_cells_dispersion_per_metacell_per_module`](@ref). For each metacell and each gene module of
the metacell's block, set the value to the ratio between the actual standard deviation of the (normalized) UMIs of the
module in the cells of the metacell, and the expected standard deviation assuming all the noise is technical multinomial
sampling noise.

We normalize each cell total UMIs to the `normalized_UMIs_quantile` of the cells in each metacell. The expected standard
deviation is the square root of the mean normalized UMIs of the module in the cells (this reflect pure multinomial
sampling noise). We set the dispersion to 0 for non-found modules, and found modules whose mean normalized UMIs in the
cells is below `min_module_UMIs`.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_n_cells_per_metacell(RequiredInput),
        matrix_of_is_found_per_module_per_block(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_cells_dispersion_per_metacell_per_module(CreatedOutput),
    ],
) function compute_matrix_of_cells_dispersion_per_metacell_per_module!(  # UNTESTED
    daf::DafWriter;
    normalized_UMIs_quantile::AbstractFloat = 0.25,
    min_module_UMIs::Integer = 9,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= normalized_UMIs_quantile <= 1
    @assert min_module_UMIs >= 0

    n_cells = axis_length(daf, "cell")
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")
    n_modules = axis_length(daf, "module")

    block_index_per_metacell = daf["@ metacell : block : index"].array
    metacell_index_per_cell = daf["@ cell : metacell ?? 0 : index"].array
    module_index_per_gene_per_block = get_query(daf, "@ gene @ block :: module ?? 0 : index").array
    is_found_per_module_per_block = plain_array(get_matrix(daf, "module", "block", "is_found"))

    sparse_UMIs_per_cell_per_gene = plain_array(get_matrix(daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_cell = plain_array(densify(get_vector(daf, "cell", "total_UMIs")))
    n_cells_per_metacell = plain_array(get_vector(daf, "metacell", "n_cells"))
    max_n_metacell_cells = maximum(n_cells_per_metacell)

    max_cells_dispersion_per_per_metacell = zeros(Float32, n_metacells)
    cells_dispersion_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)

    indices_of_metacell_max_cells_per_thread = [Vector{Int}(undef, max_n_metacell_cells) for _ in 1:maxthreadid()]
    total_UMIs_per_max_metacell_cell_per_thread =
        [Vector{Float32}(undef, max_n_metacell_cells) for _ in 1:maxthreadid()]
    normalized_factor_per_max_metacell_cell_per_thread =
        [Vector{Float32}(undef, max_n_metacell_cells) for _ in 1:maxthreadid()]
    mean_normalized_per_module_per_thread = [Vector{Float64}(undef, n_modules) for _ in 1:maxthreadid()]

    # Per-block inverted index from `module_index_per_gene`: for every (block, module) the gene indices that belong to
    # the module. Computed once and reused across all the block's metacells in the parallel loop below; replaces the
    # per-call `is_gene_in_module` BitVector that every `maximal_cells_dispersion_of_modules!` used to rebuild.
    n_blocks = size(module_index_per_gene_per_block, 2)
    gene_indices_per_module_per_block = [[Int[] for _ in 1:n_modules] for _ in 1:n_blocks]
    for block_index in 1:n_blocks
        gene_indices_per_module = gene_indices_per_module_per_block[block_index]
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        for gene_index in 1:n_genes
            module_index = module_index_per_gene[gene_index]
            if module_index > 0
                push!(gene_indices_per_module[module_index], gene_index)
            end
        end
    end

    parallel_loop_wo_rng(
        1:n_metacells;
        progress = DebugProgress(
            sum(n_cells_per_metacell);
            group = :mcs_loops,
            desc = "cells_dispersion_per_metacell_per_module",
        ),
        weights = n_cells_per_metacell,
    ) do metacell_index
        indices_of_metacell_max_cells = indices_of_metacell_max_cells_per_thread[threadid()]
        total_UMIs_per_max_metacell_cell = total_UMIs_per_max_metacell_cell_per_thread[threadid()]
        normalized_factor_per_max_metacell_cell = normalized_factor_per_max_metacell_cell_per_thread[threadid()]

        n_metacell_cells = 0
        for cell_index in 1:n_cells
            if metacell_index_per_cell[cell_index] == metacell_index
                n_metacell_cells += 1
                indices_of_metacell_max_cells[n_metacell_cells] = cell_index
            end
        end
        @assert n_metacell_cells == n_cells_per_metacell[metacell_index]

        if n_metacell_cells == 1
            return nothing
        end

        block_index = block_index_per_metacell[metacell_index]

        @views indices_of_metacell_cells = indices_of_metacell_max_cells[1:n_metacell_cells]
        @views is_found_per_module = is_found_per_module_per_block[:, block_index]
        @views cells_dispersion_per_module = cells_dispersion_per_module_per_metacell[:, metacell_index]

        max_cells_dispersion_per_per_metacell[metacell_index] = maximal_cells_dispersion_of_modules!(;
            cells_dispersion_per_module,
            mean_normalized_per_module = mean_normalized_per_module_per_thread[threadid()],
            indices_of_cells = indices_of_metacell_cells,
            sparse_UMIs_per_cell_per_gene,
            total_UMIs_per_cell,
            is_found_per_module,
            gene_indices_per_module = gene_indices_per_module_per_block[block_index],
            total_UMIs_per_max_cells = total_UMIs_per_max_metacell_cell,
            normalized_factor_per_max_cells = normalized_factor_per_max_metacell_cell,
            normalized_UMIs_quantile,
            min_module_UMIs,
        )
        return nothing
    end

    set_matrix!(
        daf,
        "module",
        "metacell",
        "cells_dispersion",
        bestify(cells_dispersion_per_module_per_metacell);
        overwrite,
    )
    @debug (
        "Mean max cells dispersion per metacell per found module: $(mean(max_cells_dispersion_per_per_metacell))"  # NOLINT
    ) _group = :mcs_results

    return nothing
end

# The z-score regularization of each block. It is specified in UMIs; convert it to a linear fraction using the mean
# total UMIs of a cell in each block's neighborhood.
function compute_std_fraction_regularization_per_block(  # UNTESTED
    daf::DafReader;
    std_UMIs_regularization::AbstractFloat,
)::Vector{Float32}
    mean_total_UMIs_in_neighborhood_cells_per_block =
        plain_array(get_vector(daf, "block", "total_neighborhood_UMIs")) ./
        plain_array(get_vector(daf, "block", "n_neighborhood_cells"))
    return Float32.(std_UMIs_regularization ./ mean_total_UMIs_in_neighborhood_cells_per_block)
end

# Fill `z_score_per_found_module_per_region_cell` (rows = the block's found modules in `findall(is_found_per_module)`
# order, columns = the given region cells) with the environment-normalized z-score of each found module's linear
# fraction - the space cells are clustered in. Fills the caller's `is_found_per_module` buffer from the block, and
# returns the per-module gene index lists (the inverted `module_index_per_gene`, reused by callers). This is the shared
# setup around [`compute_z_score_per_found_module_per_region_cell!`](@ref); the caller owns (and sizes) the reused
# buffers. If `module_UMIs_per_found_module_per_region_cell` is given, it is filled (in the same layout) with the module
# UMIs the z-scores are computed from.
function setup_z_score_per_found_module_per_region_cell!(;  # UNTESTED
    z_score_per_found_module_per_region_cell::AbstractMatrix{<:AbstractFloat},
    module_UMIs_per_found_module_per_region_cell::Maybe{AbstractMatrix{Float32}} = nothing,
    is_found_per_module::BitVector,
    block_index::Integer,
    indices_of_region_cells::AbstractVector{<:Integer},
    is_found_per_module_per_block::Union{AbstractMatrix{Bool}, BitMatrix},
    module_index_per_gene_per_block::AbstractMatrix{<:Integer},
    mean_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module_per_block::AbstractMatrix{<:AbstractFloat},
    std_fraction_regularization_per_block::AbstractVector{<:AbstractFloat},
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    z_score_accumulator_pool::Channel{Vector{Float32}},
    n_genes::Integer,
    n_modules::Integer,
)::Vector{Vector{Int}}
    is_found_per_module .= @view is_found_per_module_per_block[:, block_index]

    # Invert `module_index_per_gene` once per block: parallel-over-found-modules in compute_z_score reads the list
    # directly instead of rebuilding a per-call `is_gene_in_module` BitVector.
    @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
    gene_indices_per_module = [Int[] for _ in 1:n_modules]
    for gene_index in 1:n_genes
        module_index = module_index_per_gene[gene_index]
        if module_index > 0
            push!(gene_indices_per_module[module_index], gene_index)
        end
    end

    @views mean_linear_fraction_in_environment_cells_per_module =
        mean_linear_fraction_in_environment_cells_per_module_per_block[:, block_index]
    @views std_linear_fraction_in_environment_cells_per_module =
        std_linear_fraction_in_environment_cells_per_module_per_block[:, block_index]

    compute_z_score_per_found_module_per_region_cell!(;
        z_score_per_found_module_per_region_cell,
        module_UMIs_per_found_module_per_region_cell,
        sparse_UMIs_per_cell_per_gene,
        total_UMIs_per_cell,
        is_found_per_module,
        indices_of_region_cells,
        gene_indices_per_module,
        z_score_accumulator_pool,
        mean_linear_fraction_in_environment_cells_per_module,
        std_linear_fraction_in_environment_cells_per_module,
        std_fraction_regularization = std_fraction_regularization_per_block[block_index],
    )

    return gene_indices_per_module
end

function compute_z_score_per_found_module_per_region_cell!(;  # UNTESTED
    z_score_per_found_module_per_region_cell::AbstractMatrix{<:AbstractFloat},
    module_UMIs_per_found_module_per_region_cell::Maybe{AbstractMatrix{Float32}} = nothing,
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC{<:Integer},
    total_UMIs_per_cell::AbstractVector{<:Integer},
    is_found_per_module::BitVector,
    indices_of_region_cells::AbstractVector{<:Integer},
    gene_indices_per_module::AbstractVector{<:AbstractVector{<:Integer}},
    z_score_accumulator_pool::Channel{Vector{Float32}},
    mean_linear_fraction_in_environment_cells_per_module::AbstractVector{<:AbstractFloat},
    std_linear_fraction_in_environment_cells_per_module::AbstractVector{<:AbstractFloat},
    std_fraction_regularization::AbstractFloat,
)::Nothing
    n_region_cells = length(indices_of_region_cells)
    found_module_indices = findall(is_found_per_module)
    n_found_modules = length(found_module_indices)

    # `sum_sparse_columns_per_group!` needs the cells sorted by index, while the region cells may be in any order.
    region_position_per_sorted_position = sortperm(indices_of_region_cells)
    sorted_indices_of_region_cells = indices_of_region_cells[region_position_per_sorted_position]
    sorted_position_per_region_position = invperm(region_position_per_sorted_position)

    # Parallel over found modules. Each fills its own row of the z-score matrix. Per task, take an accumulator from the
    # pool, sum the module's gene columns into it (in the sorted order of the cells), then write the normalized row
    # and return the accumulator to the pool. The modules with the most genes go first.
    n_genes_per_found_module = [length(gene_indices_per_module[found_module_indices[k]]) for k in 1:n_found_modules]
    parallel_loop_wo_rng(
        1:n_found_modules;
        nested = true,
        policy = :greedy,
        weights = n_genes_per_found_module,
        name = "compute_z_score_per_found_module_per_region_cell",
    ) do found_module_position
        module_index = found_module_indices[found_module_position]
        accumulator = take!(z_score_accumulator_pool)
        try
            @views module_UMIs_per_sorted_region_cell = reshape(accumulator[1:n_region_cells], 1, n_region_cells)
            module_UMIs_per_sorted_region_cell .= 0.0f0
            sum_sparse_columns_per_group!(
                module_UMIs_per_sorted_region_cell,
                sparse_UMIs_per_cell_per_gene,
                sorted_indices_of_region_cells,
                view(gene_indices_per_module, module_index:module_index),
            )
            mean_linear_fraction = mean_linear_fraction_in_environment_cells_per_module[module_index]
            std_linear_fraction = std_linear_fraction_in_environment_cells_per_module[module_index]
            @assert std_linear_fraction > 0
            @check_turbo_vector(accumulator)
            @check_turbo_vector(indices_of_region_cells)
            @check_turbo_vector(sorted_position_per_region_position)
            @check_turbo_vector(total_UMIs_per_cell)
            @check_turbo_matrix(z_score_per_found_module_per_region_cell)
            @turbo for region_cell_position in 1:n_region_cells
                cell_index = indices_of_region_cells[region_cell_position]
                module_UMIs = accumulator[sorted_position_per_region_position[region_cell_position]]
                delta_fraction = module_UMIs / total_UMIs_per_cell[cell_index] - mean_linear_fraction
                z_score_per_found_module_per_region_cell[found_module_position, region_cell_position] =
                    (delta_fraction + copysign(std_fraction_regularization, delta_fraction)) /
                    (std_linear_fraction + std_fraction_regularization)
            end
            if module_UMIs_per_found_module_per_region_cell !== nothing
                for region_cell_position in 1:n_region_cells
                    module_UMIs_per_found_module_per_region_cell[found_module_position, region_cell_position] =
                        accumulator[sorted_position_per_region_position[region_cell_position]]
                end
            end
        finally
            put!(z_score_accumulator_pool, accumulator)
        end
        return nothing
    end

    return nothing
end

end  # module
