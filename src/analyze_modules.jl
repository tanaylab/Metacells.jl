"""
Do simple gene module analysis.
"""
module AnalyzeModules

export compute_matrix_of_n_genes_per_module_per_block!
export compute_stats_of_euclidean_modules_cells_distance_per_metacell!
export compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!
export compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!
export compute_vector_of_n_modules_per_block!

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
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_n_genes_per_module_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.tensor_of_linear_fraction_per_block_per_module_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_mean_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_std_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell

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
    module_index_per_gene_per_block = daf["@ gene @ block :: module ?? 0 : index"].array

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(n_blocks; group = :mcs_loops, desc = "n_genes_per_module_per_block"),
    ) do block_index
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        @views n_genes_per_module = n_genes_per_module_per_block[:, block_index]
        fill!(n_genes_per_module, 0)
        for gene_index in 1:length(module_index_per_gene)
            module_index = module_index_per_gene[gene_index]
            if module_index > 0
                n_genes_per_module[module_index] += 1
            end
        end
        return nothing
    end

    set_matrix!(daf, "module", "block", "n_genes", bestify(n_genes_per_module_per_block); overwrite)  # NOJET
    @debug "Mean genes per found module: $(mean(n_genes_per_module_per_block[n_genes_per_module_per_block .!= 0]))" _group =  # NOLINT
        :mcs_results

    return nothing
end

"""
    compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        tensor_of_linear_fraction_per_block_per_module_per_metacell(CreatedOutput),
    ],
) function compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")
    name_per_block = axis_vector(daf, "block")

    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array
    block_index_per_metacell = daf["@ metacell : block : index"].array

    module_index_per_gene_per_block = daf["@ gene @ block :: module ?? 0 : index"].array

    # TODO: This does a lot of memory allocations in the parallel loop.
    # Given it also sets a matrix in each iteration, there's probably not much point in optimization?
    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(
            n_blocks;
            group = :mcs_loops,
            desc = "linear_fraction_per_block_per_module_per_metacell",
        ),
    ) do block_index
        linear_fraction_per_module_per_metacell = zeros(Float32, n_modules, n_metacells)
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        indices_of_neighborhood_metacells = findall(is_in_neighborhood_per_other_block[block_index_per_metacell])
        n_neighborhood_metacells = length(indices_of_neighborhood_metacells)
        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]

        for module_index in 1:n_modules
            indices_of_module_genes = findall(module_index_per_gene .== module_index)
            if length(indices_of_module_genes) > 0
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

        block_name = name_per_block[block_index]
        set_matrix!(
            daf,
            "module",
            "metacell",
            "$(block_name)_linear_fraction",
            bestify(linear_fraction_per_module_per_metacell);
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
) function compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    n_modules = axis_length(daf, "module")

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    module_index_per_gene_per_block = daf["@ gene @ block :: module ?? 0 : index"].array

    is_in_neighborhood_per_other_block_per_base_block = get_matrix(daf, "block", "block", "is_in_neighborhood").array
    block_index_per_cell = daf["@ cell : metacell ?? 0 : block : index"].array

    n_cells, n_genes = size(UMIs_per_cell_per_gene)
    mean_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    std_linear_fraction_per_module_per_block = Matrix{Float32}(undef, n_modules, n_blocks)
    is_in_neighborhood_per_cell_per_thread = [BitVector(undef, n_cells) for _ in 1:nthreads()]
    is_gene_in_module_per_thread = [BitVector(undef, n_genes) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_blocks;
        progress = DebugProgress(
            n_blocks;
            group = :mcs_loops,
            desc = "stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block",
        ),
        policy = :static,
    ) do block_index
        is_in_neighborhood_per_cell = is_in_neighborhood_per_cell_per_thread[threadid()]
        is_gene_in_module = is_gene_in_module_per_thread[threadid()]
        @views is_in_neighborhood_per_other_block = is_in_neighborhood_per_other_block_per_base_block[:, block_index]
        is_in_neighborhood_per_cell .=
            (block_index_per_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(block_index_per_cell, 1))
        n_neighborhood_cells = count(is_in_neighborhood_per_cell)
        @assert n_neighborhood_cells > 0

        @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
        for module_index in 1:n_modules
            @. is_gene_in_module = module_index_per_gene == module_index
            sum_fractions = 0.0
            sum_squared_fractions = 0.0
            @foreach_true_index is_in_neighborhood_per_cell cell_index begin
                module_UMIs = 0
                @foreach_true_index is_gene_in_module gene_index begin
                    module_UMIs += UMIs_per_cell_per_gene[cell_index, gene_index]
                end
                fraction = module_UMIs / Float64(total_UMIs_per_cell[cell_index])
                sum_fractions += fraction
                sum_squared_fractions += fraction * fraction
            end
            mean_fraction = sum_fractions / n_neighborhood_cells
            mean_linear_fraction_per_module_per_block[module_index, block_index] = Float32(mean_fraction)
            std_linear_fraction_per_module_per_block[module_index, block_index] = if n_neighborhood_cells > 1
                Float32(
                sqrt(max((sum_squared_fractions - sum_fractions * mean_fraction) / (n_neighborhood_cells - 1), 0.0)),
            )
            else
                0.0f0
            end
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
    compute_stats_of_euclidean_modules_cells_distance_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_mean_euclidean_modules_cells_distance_per_metacell`](@ref) and
and set [`vector_of_std_euclidean_modules_cells_distance_per_metacell`](@ref).

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
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_mean_euclidean_modules_cells_distance_per_metacell(CreatedOutput),
        vector_of_std_euclidean_modules_cells_distance_per_metacell(CreatedOutput),
    ],
) function compute_stats_of_euclidean_modules_cells_distance_per_metacell!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")

    n_modules = axis_length(daf, "module")
    name_per_module = axis_vector(daf, "module")

    block_index_per_metacell = daf["@ metacell : block : index"].array
    module_per_gene_per_block = daf["@ gene @ block :: module"].array
    metacell_index_per_cell = daf["@ cell : metacell ?? 0 : index"].array

    mean_distance_per_metacell = zeros(Float32, n_metacells)
    std_distance_per_meyacell = zeros(Float32, n_metacells)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array

    # TODO: This does a lot of memory allocations in the parallel loop.
    parallel_loop_wo_rng(
        1:n_metacells;
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "stats_of_euclidean_modules_cells_distance_per_metacell",
        ),
    ) do metacell_index
        block_index = block_index_per_metacell[metacell_index]

        indices_of_metacell_cells = findall(metacell_index_per_cell .== metacell_index)
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

        distance_per_cell = flame_timed("pairwise.Euclidean") do
            return pairwise(Euclidean(), Ref(mean_fraction_per_module), eachcol(fraction_per_module_per_cell))
        end
        mean_distance_per_metacell[metacell_index] = mean(distance_per_cell)  # NOLINT
        std_distance_per_meyacell[metacell_index] = std(distance_per_cell)  # NOLINT

        return nothing
    end

    set_vector!(
        daf,
        "metacell",
        "mean_euclidean_modules_cells_distance",
        bestify(mean_distance_per_metacell);
        overwrite,
    )
    set_vector!(daf, "metacell", "std_euclidean_modules_cells_distance", bestify(std_distance_per_meyacell); overwrite)
    return nothing
end

end  # module
