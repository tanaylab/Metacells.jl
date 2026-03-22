"""
Project cells on an atlas with modules.
"""
module ProjectCells

export compute_cells_projection!
export compute_vector_of_correlation_between_cells_and_projected_metacells!
export compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!

using Base.Threads
using DataAxesFormats
using Distances
using LinearAlgebra
using LoopVectorization
using ProgressMeter
using StatsBase
using TanayLabUtilities

using ..Contracts
using ..Defaults

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_marker_per_gene_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.projected_block_axis
import Metacells.Contracts.tensor_of_linear_fraction_per_block_per_module_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_correlation_between_cells_and_projected_metacells_per_gene
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_mean_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_projected_block_per_cell
import Metacells.Contracts.vector_of_projected_metacell_per_cell
import Metacells.Contracts.vector_of_projected_modules_z_score_per_cell
import Metacells.Contracts.vector_of_std_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell

"""
    compute_cells_projection(;
        query_daf::DafWriter,
        atlas_daf::DafReader,
        gene_fraction_regularization::Real = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_projected_block_per_cell`] and [`vector_of_projected_modules_z_score_per_cell`](@ref). To pick
the best metacell in the `atlas_daf` for each cell of the `query_daf`, we:

  - Pick a provisional block for each cell. This is the block with the minimal Euclidean distance of the non-lateral marker
    genes.
  - Then, for each cell, consider the expression level of the found modules of the provisional block's neighborhood. We use
    this to pick a metacell in the neighborhood which has the closest Euclidean distance.
  - This metacell might belong to a different block. If so, we repeat the process using that block's neighborhood modules.
  - If the resulting best match metacell is in the same (new) block, we accept it.
  - Otherwise, we just look for the closest metacell in the original block (using that block's modules) and settle for that.

Distances are computed on the log (base 2) of the gene expression using the `gene_fraction_regularization` to handle
zero fractions. We also compute the z-score (final distance between the cell and projected metacell, minus the mean
distance of the cells in the metacell, divided by their standard deviation). If this is too large, then the query cell
isn't a good match for anything in the atlas.

# Query

$(CONTRACT1)

# Atlas

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "query_daf",
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_projected_metacell_per_cell(CreatedOutput),
        vector_of_projected_modules_z_score_per_cell(CreatedOutput),
        vector_of_projected_block_per_cell(CreatedOutput),
    ],
) Contract(;
    name = "atlas_daf",
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        matrix_of_is_found_per_module_per_block(RequiredInput),
        matrix_of_module_per_gene_per_block(RequiredInput),
        matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block(RequiredInput),
        tensor_of_linear_fraction_per_block_per_module_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_mean_euclidean_modules_cells_distance_per_metacell(RequiredInput),
        vector_of_std_euclidean_modules_cells_distance_per_metacell(RequiredInput),
    ],
) function compute_cells_projection!(;
    query_daf::DafWriter,
    atlas_daf::DafReader,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    query_gene_index_per_atlas_pertinent_marker, query_gene_index_per_atlas_gene =
        map_query_atlas_genes(; query_daf, atlas_daf)

    UMIs_per_query_cell_per_query_gene = get_matrix(query_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_query_cell = get_vector(query_daf, "cell", "total_UMIs").array

    provisional_block_index_per_query_cell = compute_provisional_projection_per_query_cell!(;
        query_daf,
        atlas_daf,
        gene_fraction_regularization,
        UMIs_per_query_cell_per_query_gene,
        total_UMIs_per_query_cell,
        query_gene_index_per_atlas_pertinent_marker,
    )

    compute_final_projection_per_query_cell(;
        query_daf,
        atlas_daf,
        query_gene_index_per_atlas_gene,
        UMIs_per_query_cell_per_query_gene,
        total_UMIs_per_query_cell,
        provisional_block_index_per_query_cell,
        overwrite,
    )

    return nothing
end

function map_query_atlas_genes(;
    query_daf::DafReader,
    atlas_daf::DafReader,
)::Tuple{AbstractVector{<:Integer}, AbstractVector{<:Integer}}
    atlas_genes = axis_vector(atlas_daf, "gene")

    is_in_module_per_atlas_gene = atlas_daf["@ block @ gene :: module != '' >- Max"].array
    atlas_in_module_genes = atlas_genes[is_in_module_per_atlas_gene]
    axis_indices(query_daf, "gene", atlas_in_module_genes)

    is_pertinent_marker_per_atlas_gene =
        get_vector(atlas_daf, "gene", "is_marker").array .& .! get_vector(atlas_daf, "gene", "is_lateral").array
    atlas_pertinent_marker_genes = atlas_genes[is_pertinent_marker_per_atlas_gene]

    query_gene_index_per_atlas_pertinent_marker = axis_indices(query_daf, "gene", atlas_pertinent_marker_genes)
    query_gene_index_per_atlas_gene = axis_indices(query_daf, "gene", atlas_genes; allow_missing = true)

    return query_gene_index_per_atlas_pertinent_marker, query_gene_index_per_atlas_gene
end

function compute_provisional_projection_per_query_cell!(;
    query_daf::DafWriter,
    atlas_daf::DafReader,
    gene_fraction_regularization::AbstractFloat,
    UMIs_per_query_cell_per_query_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_query_cell::AbstractVector{<:Integer},
    query_gene_index_per_atlas_pertinent_marker::AbstractVector{<:Integer},
)::AbstractVector{<:Integer}
    n_blocks = axis_length(atlas_daf, "block")

    is_excluded_per_query_cell = get_vector(query_daf, "cell", "is_excluded").array

    linear_fraction_per_block_per_pertinent_marker = atlas_daf["@ block @ gene [ is_marker & ! is_lateral ] :: linear_fraction"].array
    linear_fraction_per_block_per_pertinent_marker =
    flame_timed("linear_fraction_per_block_per_pertinent_marker") do
        return mutable_array(densify(linear_fraction_per_block_per_pertinent_marker))
    end

    log_linear_fraction_per_pertinent_marker_per_block = Matrix{Float32}(
        undef,
        size(linear_fraction_per_block_per_pertinent_marker, 2),
        size(linear_fraction_per_block_per_pertinent_marker, 1),
    )

    flame_timed("log_linear_fraction_per_pertinent_marker_per_block") do
        @assert LoopVectorization.check_args(log_linear_fraction_per_pertinent_marker_per_block) "check_args failed in compute_provisional_projection_per_query_cell! (log_linear_fraction)\nfor log_linear_fraction_per_pertinent_marker_per_block: $(brief(log_linear_fraction_per_pertinent_marker_per_block))"
        @assert LoopVectorization.check_args(linear_fraction_per_block_per_pertinent_marker) "check_args failed in compute_provisional_projection_per_query_cell! (log_linear_fraction)\nfor linear_fraction_per_block_per_pertinent_marker: $(brief(linear_fraction_per_block_per_pertinent_marker))"
        @turbo for j in axes(linear_fraction_per_block_per_pertinent_marker, 1)
            for i in axes(linear_fraction_per_block_per_pertinent_marker, 2)
                log_linear_fraction_per_pertinent_marker_per_block[i, j] =
                    log2(linear_fraction_per_block_per_pertinent_marker[j, i] + gene_fraction_regularization)
            end
        end
    end

    n_query_cells = axis_length(query_daf, "cell")

    provisional_block_index_per_query_cell = Vector{UInt32}(undef, n_query_cells)

    n_pertinent_neighborhood_markers = length(query_gene_index_per_atlas_pertinent_marker)
    log_linear_fraction_per_pertinent_marker_per_thread =
        [Vector{Float32}(undef, n_pertinent_neighborhood_markers) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_query_cells;
        name = "provisional_projection_per_cell",
        policy = :static,
        progress = DebugProgress(n_query_cells; group = :mcs_loops, desc = "provisional_projection_per_cell"),
    ) do query_cell_index
        if is_excluded_per_query_cell[query_cell_index]
            provisional_block_index_per_query_cell[query_cell_index] = 0
        else
            total_query_cell_UMIs = total_UMIs_per_query_cell[query_cell_index]
            @views UMIs_per_pertinent_marker =
                UMIs_per_query_cell_per_query_gene[query_cell_index, query_gene_index_per_atlas_pertinent_marker]

            log_linear_fraction_per_pertinent_marker = log_linear_fraction_per_pertinent_marker_per_thread[threadid()]
            for gene_position in 1:n_pertinent_neighborhood_markers
                log_linear_fraction_per_pertinent_marker[gene_position] =
                    UMIs_per_pertinent_marker[gene_position] / total_query_cell_UMIs
            end
            @assert LoopVectorization.check_args(log_linear_fraction_per_pertinent_marker) "check_args failed in compute_provisional_projection_per_query_cell! (log2)\nfor log_linear_fraction_per_pertinent_marker: $(brief(log_linear_fraction_per_pertinent_marker))"
            @turbo for gene_position in 1:n_pertinent_neighborhood_markers
                log_linear_fraction_per_pertinent_marker[gene_position] =
                    log2(log_linear_fraction_per_pertinent_marker[gene_position] + gene_fraction_regularization)
            end

            distances_to_blocks = flame_timed("pairwise.Euclidean") do
                return vec(
                    pairwise(
                        Euclidean(),
                        Ref(log_linear_fraction_per_pertinent_marker),
                        eachcol(log_linear_fraction_per_pertinent_marker_per_block),
                    ),
                )
            end
            @assert_vector(distances_to_blocks, n_blocks)

            provisional_block_index_per_query_cell[query_cell_index] = argmin(distances_to_blocks)
        end
        return nothing
    end

    return provisional_block_index_per_query_cell
end

function compute_final_projection_per_query_cell(;
    query_daf::DafWriter,
    atlas_daf::DafReader,
    query_gene_index_per_atlas_gene::AbstractVector{<:Integer},
    UMIs_per_query_cell_per_query_gene::AbstractMatrix{<:Integer},
    total_UMIs_per_query_cell::AbstractVector{<:Integer},
    provisional_block_index_per_query_cell::AbstractVector{<:Integer},
    overwrite::Bool,
)::Nothing
    n_query_cells = axis_length(query_daf, "cell")

    n_blocks = axis_length(atlas_daf, "block")
    name_per_block = axis_vector(atlas_daf, "block")
    name_per_metacell = axis_vector(atlas_daf, "metacell")

    module_index_per_atlas_gene_per_block = atlas_daf["@ gene @ block :: module ?? 0 : index"].array
    is_found_per_module_per_block = get_matrix(atlas_daf, "module", "block", "is_found").array

    is_in_neighborhood_per_other_block_per_base_block =
        get_matrix(atlas_daf, "block", "block", "is_in_neighborhood").array
    block_index_per_metacell = atlas_daf["@ metacell : block : index"].array

    mean_linear_fraction_in_neighborhood_cells_per_module_per_block =
        get_matrix(atlas_daf, "module", "block", "mean_linear_fraction_in_neighborhood_cells").array
    std_linear_fraction_in_neighborhood_cells_per_module_per_block =
        get_matrix(atlas_daf, "module", "block", "std_linear_fraction_in_neighborhood_cells").array

    mean_euclidean_modules_cells_distance_per_metacell =
        get_vector(atlas_daf, "metacell", "mean_euclidean_modules_cells_distance").array
    std_euclidean_modules_cells_distance_per_metacell =
        get_vector(atlas_daf, "metacell", "std_euclidean_modules_cells_distance").array

    projected_metacell_per_query_cell = Vector{AbstractString}(undef, n_query_cells)
    projected_metacell_per_query_cell .= ""
    projected_modules_z_score_per_query_cell = zeros(Float32, n_query_cells)

    next_indices_of_undetermined_query_cells_per_block =
        [findall(provisional_block_index_per_query_cell .== block_index) for block_index in 1:n_blocks]

    spin_lock = SpinLock()

    n_determined_query_cells = Atomic{Int}(0)

    for phase in 1:3
        phase_name = "final_projection_per_cell.phase_$(phase)"
        flame_timed(phase_name) do
            indices_of_undetermined_query_cells_per_block = next_indices_of_undetermined_query_cells_per_block
            next_indices_of_undetermined_query_cells_per_block = [Int32[] for _ in 1:n_blocks]
            @debug "Determined: $(n_determined_query_cells[]) $(percent(n_determined_query_cells[], n_query_cells)) Phase $(phase)..." _group =
                :mcs_details

            # TODO: Many memory allocations inside the parallel loop.
            parallel_loop_wo_rng(
                1:n_blocks;
                name = "final_projection_per_cell.loop",
                progress = DebugProgress(n_blocks; group = :mcs_loops, desc = phase_name),
            ) do block_index
                block_name = name_per_block[block_index]

                indices_of_undetermined_query_cells = indices_of_undetermined_query_cells_per_block[block_index]
                n_undetermined_query_cells = length(indices_of_undetermined_query_cells)
                if n_undetermined_query_cells > 0
                    @views is_in_neighborhood_per_other_block =
                        is_in_neighborhood_per_other_block_per_base_block[:, block_index]
                    if 1 <= phase <= 2
                        indices_of_candidate_metacells =
                            findall(is_in_neighborhood_per_other_block[block_index_per_metacell])
                        block_index_per_candidate_metacell = block_index_per_metacell[indices_of_candidate_metacells]
                    elseif phase == 3
                        indices_of_candidate_metacells = findall(block_index_per_metacell .== block_index)
                        block_index_per_candidate_metacell = nothing
                    else
                        @assert false
                    end

                    @views module_index_per_atlas_gene = module_index_per_atlas_gene_per_block[:, block_index]
                    @views is_found_per_module = is_found_per_module_per_block[:, block_index]
                    indices_of_found_modules = findall(is_found_per_module)
                    n_found_modules = length(indices_of_found_modules)
                    @assert n_found_modules > 0

                    linear_fraction_per_module_per_metacell =
                        get_matrix(atlas_daf, "module", "metacell", "$(block_name)_linear_fraction").array

                    linear_fraction_per_found_module_per_candidate_metacell = linear_fraction_per_module_per_metacell[
                        indices_of_found_modules,
                        indices_of_candidate_metacells,
                    ]

                    @views module_index_per_atlas_gene = module_index_per_atlas_gene_per_block[:, block_index]
                    total_UMIs_per_undetermined_query_cell =
                        total_UMIs_per_query_cell[indices_of_undetermined_query_cells]
                    linear_fraction_per_found_module_per_undetermined_query_cell =
                        Matrix{Float32}(undef, n_found_modules, n_undetermined_query_cells)

                    for (found_module_position, found_module_index) in enumerate(indices_of_found_modules)
                        indices_of_found_module_atlas_genes =
                            findall(module_index_per_atlas_gene .== found_module_index)
                        indices_of_found_module_query_genes =
                            query_gene_index_per_atlas_gene[indices_of_found_module_atlas_genes]
                        @assert !any(indices_of_found_module_atlas_genes .== 0)

                        n_found_module_genes = length(indices_of_found_module_query_genes)
                        @assert n_found_module_genes > 0
                        found_module_UMIs_per_undetermined_query_cell = vec(
                            sum(
                                UMIs_per_query_cell_per_query_gene[
                                    indices_of_undetermined_query_cells,
                                    indices_of_found_module_query_genes,
                                ];
                                dims = 2,
                            ),
                        )
                        @assert_vector(found_module_UMIs_per_undetermined_query_cell, n_undetermined_query_cells)
                        linear_fraction_per_found_module_per_undetermined_query_cell[found_module_position, :] .=
                            found_module_UMIs_per_undetermined_query_cell ./ total_UMIs_per_undetermined_query_cell
                    end

                    mean_linear_fraction_in_neighborhood_cells_per_found_module =
                        mean_linear_fraction_in_neighborhood_cells_per_module_per_block[
                            indices_of_found_modules,
                            block_index,
                        ]
                    std_linear_fraction_in_neighborhood_cells_per_found_module =
                        std_linear_fraction_in_neighborhood_cells_per_module_per_block[
                            indices_of_found_modules,
                            block_index,
                        ]
                    z_score_per_found_module_per_candidate_metacell =
                        (
                            linear_fraction_per_found_module_per_candidate_metacell .-
                            mean_linear_fraction_in_neighborhood_cells_per_found_module
                        ) ./ std_linear_fraction_in_neighborhood_cells_per_found_module

                    z_score_per_found_module_per_undetermined_query_cell =
                        (
                            linear_fraction_per_found_module_per_undetermined_query_cell .-
                            mean_linear_fraction_in_neighborhood_cells_per_found_module
                        ) ./ std_linear_fraction_in_neighborhood_cells_per_found_module

                    distances_between_candidate_metacells_and_undetermined_query_cells =
                        flame_timed("pairwise.Euclidean") do
                            return pairwise(
                                Euclidean(),
                                z_score_per_found_module_per_candidate_metacell,
                                z_score_per_found_module_per_undetermined_query_cell,
                            )
                        end

                    minimal_position_per_undetermined_query_cell =
                        vec(argmin(distances_between_candidate_metacells_and_undetermined_query_cells; dims = 1))
                    minimal_distance_per_undetermined_query_cell =
                        distances_between_candidate_metacells_and_undetermined_query_cells[minimal_position_per_undetermined_query_cell]
                    position_of_nearest_candidate_metacell_per_undetermined_query_cell =
                        [position.I[1] for position in minimal_position_per_undetermined_query_cell]

                    if block_index_per_candidate_metacell === nothing
                        positions_of_stable_undetermined_query_cells = 1:n_undetermined_query_cells
                    else
                        nearest_block_index_per_undetermined_query_cell =
                            block_index_per_candidate_metacell[position_of_nearest_candidate_metacell_per_undetermined_query_cell]
                        positions_of_stable_undetermined_query_cells =
                            findall(nearest_block_index_per_undetermined_query_cell .== block_index)

                        try
                            lock(spin_lock)
                            for (undetermined_query_cell_position, nearest_block_index) in
                                enumerate(nearest_block_index_per_undetermined_query_cell)
                                if nearest_block_index != block_index
                                    query_cell_index =
                                        indices_of_undetermined_query_cells[undetermined_query_cell_position]
                                    if phase == 1
                                        push!(
                                            next_indices_of_undetermined_query_cells_per_block[nearest_block_index],
                                            query_cell_index,
                                        )
                                    elseif phase == 2
                                        provisional_block_index =
                                            provisional_block_index_per_query_cell[query_cell_index]
                                        @assert block_index != provisional_block_index
                                        push!(
                                            next_indices_of_undetermined_query_cells_per_block[provisional_block_index],
                                            query_cell_index,
                                        )
                                    else
                                        @assert false
                                    end
                                end
                            end
                        finally
                            unlock(spin_lock)
                        end
                    end

                    n_stable_query_cells = length(positions_of_stable_undetermined_query_cells)
                    if n_stable_query_cells > 0
                        position_of_nearest_candidate_metacell_per_stable_undetermined_query_cell =
                            position_of_nearest_candidate_metacell_per_undetermined_query_cell[positions_of_stable_undetermined_query_cells]
                        metacell_index_per_stable_undetermined_query_cell =
                            indices_of_candidate_metacells[position_of_nearest_candidate_metacell_per_stable_undetermined_query_cell]
                        metacell_per_stable_undetermined_query_cell =
                            name_per_metacell[metacell_index_per_stable_undetermined_query_cell]

                        indices_of_stable_query_cells =
                            indices_of_undetermined_query_cells[positions_of_stable_undetermined_query_cells]
                        projected_metacell_per_query_cell[indices_of_stable_query_cells] =
                            metacell_per_stable_undetermined_query_cell

                        projected_modules_z_score_per_query_cell[indices_of_stable_query_cells] .=
                            (
                                minimal_distance_per_undetermined_query_cell[positions_of_stable_undetermined_query_cells] .-
                                mean_euclidean_modules_cells_distance_per_metacell[metacell_index_per_stable_undetermined_query_cell]
                            ) ./
                            (std_euclidean_modules_cells_distance_per_metacell[metacell_index_per_stable_undetermined_query_cell])
                    end
                    atomic_add!(n_determined_query_cells, n_stable_query_cells)
                end
                return nothing
            end
        end
    end

    @assert n_determined_query_cells[] == n_query_cells

    set_vector!(query_daf, "cell", "projected_metacell", projected_metacell_per_query_cell; overwrite)
    set_vector!(
        query_daf,
        "cell",
        "projected_metacell_modules_z_score",
        projected_modules_z_score_per_query_cell;
        overwrite,
    )

    projected_metacell_index_per_query_cell = axis_indices(atlas_daf, "metacell", projected_metacell_per_query_cell)
    block_per_metacell = get_vector(atlas_daf, "metacell", "block").array
    projected_block_per_query_cell = block_per_metacell[projected_metacell_index_per_query_cell]
    set_vector!(query_daf, "cell", "projected_block", projected_block_per_query_cell; overwrite)

    return nothing
end

"""
    function compute_vector_of_correlation_between_cells_and_projected_metacells!(;
        atlas_daf::DafReader,
        query_daf::DafWriter,
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set (in query) [`vector_of_correlation_between_cells_and_projected_metacells_per_gene`](@ref).

# Query

$(CONTRACT1)

# Atlas

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "query_daf",
    axes = [gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_cell(RequiredInput),
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_projected_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_correlation_between_cells_and_projected_metacells_per_gene(CreatedOutput),
    ],
) Contract(
    name = "atlas_daf",
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
    ],
) function compute_vector_of_correlation_between_cells_and_projected_metacells!(;  # UNTESTED
    query_daf::DafWriter,
    atlas_daf::DafReader,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    gene_fraction_regularization = Float32(gene_fraction_regularization)
    n_genes_in_atlas = axis_length(atlas_daf, "gene")
    n_genes_in_query = axis_length(query_daf, "gene")

    query_included_genes = query_daf["@ gene [ ! is_excluded ]"]
    atlas_included_genes = atlas_daf["@ gene [ ! is_excluded ]"]
    common_included_genes = intersect(query_included_genes, atlas_included_genes)
    n_common_included_genes = length(common_included_genes)

    @debug "Common included genes: $(n_common_included_genes)" *
           " ($(percent(n_common_included_genes, n_genes_in_atlas)) out of $(n_genes_in_atlas) atlas genes" *
           ", $(percent(n_common_included_genes, n_genes_in_query)) out of $(n_genes_in_query) query genes)" _group =
        :mcs_results

    indices_of_common_included_atlas_genes = axis_indices(atlas_daf, "gene", common_included_genes)
    indices_of_common_included_query_genes = axis_indices(query_daf, "gene", common_included_genes)

    indices_of_included_query_cells = query_daf["@ cell [ ! is_excluded ] : index"].array
    n_included_query_cells = length(indices_of_included_query_cells)

    UMIs_per_query_cell_per_gene = get_matrix(query_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_query_cell = get_vector(query_daf, "cell", "total_UMIs").array

    projected_metacell_per_query_cell = get_vector(query_daf, "cell", "projected_metacell").array
    projected_metacell_per_included_query_cell = projected_metacell_per_query_cell[indices_of_included_query_cells]
    atlas_metacell_index_per_included_query_cell =
        axis_indices(atlas_daf, "metacell", projected_metacell_per_included_query_cell)

    UMIs_per_atlas_metacell_per_gene = get_matrix(atlas_daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_atlas_metacell = get_vector(atlas_daf, "metacell", "total_UMIs").array

    correlation_between_cells_and_projected_metacells_per_query_gene = zeros(Float32, n_genes_in_query)

    gene_cell_log_fraction_per_included_query_cell_per_thread =
        [Vector{Float32}(undef, n_included_query_cells) for _ in 1:nthreads()]
    gene_metacell_log_fraction_per_included_query_cell_per_thread =
        [Vector{Float32}(undef, n_included_query_cells) for _ in 1:nthreads()]

    todox_all_zero = Atomic{Int32}(0)
    parallel_loop_wo_rng(
        1:n_common_included_genes;
        policy = :static,
        progress = DebugProgress(
            n_common_included_genes;
            group = :mcs_loops,
            desc = "correlation_between_cells_and_projected_metacells",
        ),
        progress_chunk = 100,
    ) do common_included_gene_position
        atlas_gene_index = indices_of_common_included_atlas_genes[common_included_gene_position]
        query_gene_index = indices_of_common_included_query_genes[common_included_gene_position]

        gene_cell_log_fraction_per_included_query_cell =
            gene_cell_log_fraction_per_included_query_cell_per_thread[threadid()]
        gene_metacell_log_fraction_per_included_query_cell =
            gene_metacell_log_fraction_per_included_query_cell_per_thread[threadid()]

        all_zero_included_query_cell_UMIs = true
        for included_query_cell_position in 1:n_included_query_cells
            if UMIs_per_query_cell_per_gene[
                indices_of_included_query_cells[included_query_cell_position],
                query_gene_index,
            ] > 0
                all_zero_included_query_cell_UMIs = false
                break
            end
        end
        if all_zero_included_query_cell_UMIs
            atomic_add!(todox_all_zero, Int32(1))
            correlation_between_cells_and_projected_metacells_per_query_gene[query_gene_index] = 0.0f0
            return nothing
        end

        all_zero_included_atlas_metacell_UMIs = true
        for included_query_cell_position in 1:n_included_query_cells
            if UMIs_per_atlas_metacell_per_gene[
                atlas_metacell_index_per_included_query_cell[included_query_cell_position],
                atlas_gene_index,
            ] > 0
                all_zero_included_atlas_metacell_UMIs = false
                break
            end
        end
        if all_zero_included_atlas_metacell_UMIs
            atomic_add!(todox_all_zero, Int32(1))
            correlation_between_cells_and_projected_metacells_per_query_gene[query_gene_index] = 0.0f0
            return nothing
        end

        for included_query_cell_position in 1:n_included_query_cells
            atlas_metacell_index = atlas_metacell_index_per_included_query_cell[included_query_cell_position]
            query_cell_index = indices_of_included_query_cells[included_query_cell_position]
            atlas_metacell_UMIs = UMIs_per_atlas_metacell_per_gene[atlas_metacell_index, atlas_gene_index]
            query_cell_UMIs = UMIs_per_query_cell_per_gene[query_cell_index, query_gene_index]
            gene_metacell_log_fraction_per_included_query_cell[included_query_cell_position] =
                atlas_metacell_UMIs / total_UMIs_per_atlas_metacell[atlas_metacell_index] + gene_fraction_regularization
            gene_cell_log_fraction_per_included_query_cell[included_query_cell_position] =
                query_cell_UMIs / total_UMIs_per_query_cell[query_cell_index] + gene_fraction_regularization
        end
        @assert LoopVectorization.check_args(gene_metacell_log_fraction_per_included_query_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_projected_metacells!\nfor gene_metacell_log_fraction_per_included_query_cell: $(brief(gene_metacell_log_fraction_per_included_query_cell))"
        @assert LoopVectorization.check_args(gene_cell_log_fraction_per_included_query_cell) "check_args failed in compute_vector_of_correlation_between_cells_and_projected_metacells!\nfor gene_cell_log_fraction_per_included_query_cell: $(brief(gene_cell_log_fraction_per_included_query_cell))"
        @turbo for included_query_cell_position in 1:n_included_query_cells
            gene_metacell_log_fraction_per_included_query_cell[included_query_cell_position] =
                log2(gene_metacell_log_fraction_per_included_query_cell[included_query_cell_position])
            gene_cell_log_fraction_per_included_query_cell[included_query_cell_position] =
                log2(gene_cell_log_fraction_per_included_query_cell[included_query_cell_position])
        end

        correlation_between_cells_and_projected_metacells_per_query_gene[query_gene_index] = zero_cor_between_vectors(
            gene_cell_log_fraction_per_included_query_cell,
            gene_metacell_log_fraction_per_included_query_cell,
        )

        return nothing
    end

    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(n_common_included_genes) ($(percent(todox_all_zero[], n_common_included_genes)))" _group =
        :todox
    set_vector!(
        query_daf,
        "gene",
        "correlation_between_cells_and_projected_metacells",
        bestify(correlation_between_cells_and_projected_metacells_per_query_gene);
        overwrite,
    )

    atlas_pertinent_markers = atlas_daf["@ gene [ is_marker & ! is_lateral ]"].array
    common_pertinent_markers = intersect(atlas_pertinent_markers, common_included_genes)
    index_in_query_per_common_pertinent_marker = axis_indices(query_daf, "gene", common_pertinent_markers)

    if query_daf isa DataAxesFormats.Contracts.ContractDaf
        query_daf = query_daf.daf
    end
    if atlas_daf isa DataAxesFormats.Contracts.ContractDaf
        atlas_daf = atlas_daf.daf
    end
    if query_daf === atlas_daf
        qualifier = "self"
    else
        qualifier = "common"
    end
    @debug (
        "Mean correlation of $(qualifier) pertinent marker genes between cells and their projected metacells: " *
        "$(mean(correlation_between_cells_and_projected_metacells_per_query_gene[index_in_query_per_common_pertinent_marker]))"  # NOLINT
    ) _group = :mcs_results

    return nothing
end

"""
    compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!(
        query_daf::DafWriter,
        atlas_daf::DafReader,
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_neighborhood_query_cells::Integer = $(DEFAULT.min_neighborhood_query_cells),
        overwrite::Bool = false,
    )::Nothing

Compute and set [`matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block`](@ref).
If there are less than `min_neighborhood_query_cells`, we set this to zero. This will create [`projected_block_axis`](@ref)
in the `query_daf` if necessary.

# Query

$(CONTRACT1)

# Atlas

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(
    name = "query_daf",
    axes = [gene_axis(RequiredInput), cell_axis(RequiredInput), projected_block_axis(GuaranteedOutput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_projected_metacell_per_cell(RequiredInput),
        vector_of_projected_block_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block(
            CreatedOutput,
        ),
    ],
) Contract(
    name = "atlas_daf",
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(OptionalInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
    ],
) function compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!(;
    query_daf::DafWriter,
    atlas_daf::DafReader,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    min_neighborhood_query_cells::Integer = 30,
    overwrite::Bool = false,
    bin::Maybe{Integer} = nothing,
)::Nothing
    @assert 0 <= gene_fraction_regularization <= 1
    gene_fraction_regularization = Float32(gene_fraction_regularization)
    @assert min_neighborhood_query_cells > 1

    if has_axis(query_daf, "projected_block")
        @assert axis_vector(atlas_daf, "block") == axis_vector(query_daf, "projected_block")
    else
        copy_axis!(source = atlas_daf, destination = query_daf, axis = "block", rename = "projected_block"; overwrite)
    end

    n_query_cells = axis_length(query_daf, "cell")

    n_atlas_blocks = axis_length(atlas_daf, "block")
    n_genes_in_atlas = axis_length(atlas_daf, "gene")
    n_genes_in_query = axis_length(query_daf, "gene")

    name_per_atlas_block = axis_vector(atlas_daf, "block")
    query_included_genes = query_daf["@ gene [ ! is_excluded ]"]
    atlas_included_genes = atlas_daf["@ gene [ ! is_excluded ]"]
    common_included_genes = intersect(query_included_genes, atlas_included_genes)
    n_common_included_genes = length(common_included_genes)

    @debug "Common included genes: $(n_common_included_genes)" *
           " ($(percent(n_common_included_genes, n_genes_in_atlas)) out of $(n_genes_in_atlas) atlas genes" *
           ", $(percent(n_common_included_genes, n_genes_in_query)) out of $(n_genes_in_query) query genes)" _group =
        :mcs_results

    indices_of_common_included_atlas_genes = axis_indices(atlas_daf, "gene", common_included_genes)
    indices_of_common_included_query_genes = axis_indices(query_daf, "gene", common_included_genes)

    UMIs_per_query_cell_per_gene = get_matrix(query_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_query_cell = get_vector(query_daf, "cell", "total_UMIs").array

    projected_metacell_per_query_cell = get_vector(query_daf, "cell", "projected_metacell").array
    atlas_metacell_index_per_query_cell =
        axis_indices(atlas_daf, "metacell", projected_metacell_per_query_cell; allow_missing = true)

    projected_block_per_query_cell = get_vector(query_daf, "cell", "projected_block").array
    atlas_block_index_per_query_cell =
        axis_indices(atlas_daf, "block", projected_block_per_query_cell; allow_missing = true)

    UMIs_per_atlas_metacell_per_gene = get_matrix(atlas_daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_atlas_metacell = get_vector(atlas_daf, "metacell", "total_UMIs").array

    is_lateral_per_atlas_gene = get_vector(atlas_daf, "gene", "is_lateral").array
    is_lateral_per_common_included_gene = is_lateral_per_atlas_gene[indices_of_common_included_atlas_genes]

    is_marker_per_atlas_gene = get_vector(atlas_daf, "gene", "is_marker").array
    is_marker_per_common_included_gene = is_marker_per_atlas_gene[indices_of_common_included_atlas_genes]
    is_pertinent_marker_per_common_included_gene =
        is_marker_per_common_included_gene .& .! is_lateral_per_common_included_gene

    is_in_neighborhood_per_other_block_per_base_block =
        get_matrix(atlas_daf, "block", "block", "is_in_neighborhood").array
    is_neighborhood_marker_per_atlas_gene_per_atlas_block =
        get_matrix(atlas_daf, "gene", "block", "is_neighborhood_marker"; default = nothing)
    if is_neighborhood_marker_per_atlas_gene_per_atlas_block !== nothing
        is_neighborhood_marker_per_atlas_gene_per_atlas_block =
            is_neighborhood_marker_per_atlas_gene_per_atlas_block.array
    end

    correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block =
        zeros(Float32, n_genes_in_query, n_atlas_blocks)

    mean_pertinent_marker_genes_correlation_per_atlas_block = nothing

    mean_neighborhood_pertinent_markers_correlation_per_atlas_block = nothing

    is_in_bin_per_common_included_gene = nothing
    is_in_bin_pertinent_marker_per_common_included_gene = nothing
    is_out_bin_pertinent_marker_per_common_included_gene = nothing

    in_bin_mean_pertinent_markers_correlation_per_atlas_block = nothing
    out_bin_mean_pertinent_markers_correlation_per_atlas_block = nothing

    in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block = nothing
    out_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block = nothing

    if bin === nothing
        mean_pertinent_marker_genes_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
        if is_neighborhood_marker_per_atlas_gene_per_atlas_block !== nothing
            mean_neighborhood_pertinent_markers_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
        end

    else
        if atlas_daf isa DataAxesFormats.Contracts.ContractDaf
            daf = atlas_daf.daf
        else
            daf = atlas_daf
        end
        is_in_bin_per_atlas_gene = get_vector(daf, "gene", "bin").array .== bin
        is_in_bin_per_common_included_gene = is_in_bin_per_atlas_gene[indices_of_common_included_atlas_genes]
        @assert any(is_in_bin_per_common_included_gene)
        @assert !all(is_in_bin_per_common_included_gene)

        is_in_bin_pertinent_marker_per_common_included_gene =
            is_pertinent_marker_per_common_included_gene .& is_in_bin_per_common_included_gene
        is_out_bin_pertinent_marker_per_common_included_gene =
            is_pertinent_marker_per_common_included_gene .& .!is_in_bin_per_common_included_gene
        @assert any(is_in_bin_pertinent_marker_per_common_included_gene)
        @assert any(is_out_bin_pertinent_marker_per_common_included_gene)

        in_bin_mean_pertinent_markers_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
        out_bin_mean_pertinent_markers_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
        if is_neighborhood_marker_per_atlas_gene_per_atlas_block !== nothing
            in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
            out_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block = zeros(Float32, n_atlas_blocks)
        end
    end

    progress = DebugProgress(
        n_atlas_blocks * n_common_included_genes;
        group = :mcs_loops,
        desc = "correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block",
    )

    is_of_neighborhood_per_query_cell = BitVector(undef, n_query_cells)

    max_n_neighborhood_cells = 0
    for atlas_block_index in 1:n_atlas_blocks
        @views is_in_neighborhood_per_other_block =
            is_in_neighborhood_per_other_block_per_base_block[:, atlas_block_index]
        is_of_neighborhood_per_query_cell .=
            (atlas_block_index_per_query_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(atlas_block_index_per_query_cell, 1))
        max_n_neighborhood_cells = max(max_n_neighborhood_cells, sum(is_of_neighborhood_per_query_cell))
    end

    cell_log_fraction_per_max_neighborhood_query_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]
    projected_metacell_log_fraction_per_max_neighborhood_query_cell_per_thread =
        [Vector{Float32}(undef, max_n_neighborhood_cells) for _ in 1:nthreads()]

    if mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing &&
       in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
        is_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
    else
        is_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
            [BitVector(undef, n_common_included_genes) for _ in 1:nthreads()]
        if in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
            is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
            is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        else
            is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
                [BitVector(undef, n_common_included_genes) for _ in 1:nthreads()]
            is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
                [BitVector(undef, n_common_included_genes) for _ in 1:nthreads()]
        end
    end

    todox_all_zero = Atomic{Int32}(0)
    for atlas_block_index in 1:n_atlas_blocks
        @views is_in_neighborhood_per_other_block =
            is_in_neighborhood_per_other_block_per_base_block[:, atlas_block_index]
        is_of_neighborhood_per_query_cell .=
            (atlas_block_index_per_query_cell .> 0) .&
            getindex.(Ref(is_in_neighborhood_per_other_block), max.(atlas_block_index_per_query_cell, 1))

        indices_of_neighborhood_query_cells = findall(is_of_neighborhood_per_query_cell)

        n_neighborhood_query_cells = length(indices_of_neighborhood_query_cells)
        if n_neighborhood_query_cells < min_neighborhood_query_cells
            if n_neighborhood_query_cells > 0
                @warn "Ignoring too few query cells: $(n_neighborhood_query_cells) for the neighborhood of the atlas block: $(name_per_atlas_block[atlas_block_index])"
            end
            next!(progress; step = n_common_included_genes)  # NOJET

        else
            atlas_metacell_index_per_neighborhood_query_cell =
                atlas_metacell_index_per_query_cell[indices_of_neighborhood_query_cells]

            total_atlas_metacell_UMIs_per_neighborhood_query_cell =
                total_UMIs_per_atlas_metacell[atlas_metacell_index_per_neighborhood_query_cell]

            parallel_loop_wo_rng(
                1:n_common_included_genes;
                policy = :static,
                progress,
                progress_chunk = 100,
            ) do common_included_gene_position
                atlas_gene_index = indices_of_common_included_atlas_genes[common_included_gene_position]
                query_gene_index = indices_of_common_included_query_genes[common_included_gene_position]

                cell_log_fraction_per_max_neighborhood_query_cell =
                    cell_log_fraction_per_max_neighborhood_query_cell_per_thread[threadid()]
                @views cell_log_fraction_per_neighborhood_query_cell =
                    cell_log_fraction_per_max_neighborhood_query_cell[1:n_neighborhood_query_cells]

                projected_metacell_log_fraction_per_max_neighborhood_query_cell =
                    projected_metacell_log_fraction_per_max_neighborhood_query_cell_per_thread[threadid()]
                @views projected_metacell_log_fraction_per_neighborhood_query_cell =
                    projected_metacell_log_fraction_per_max_neighborhood_query_cell[1:n_neighborhood_query_cells]

                all_zero_neighborhood_query_cell_UMIs = true
                for neighborhood_query_cell_position in 1:n_neighborhood_query_cells
                    if UMIs_per_query_cell_per_gene[
                        indices_of_neighborhood_query_cells[neighborhood_query_cell_position],
                        query_gene_index,
                    ] > 0
                        all_zero_neighborhood_query_cell_UMIs = false
                        break
                    end
                end
                if all_zero_neighborhood_query_cell_UMIs
                    atomic_add!(todox_all_zero, Int32(1))
                    correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block[
                        query_gene_index,
                        atlas_block_index,
                    ] = 0.0f0
                    return nothing
                end

                all_zero_neighborhood_projected_metacell_UMIs = true
                for neighborhood_query_cell_position in 1:n_neighborhood_query_cells
                    if UMIs_per_atlas_metacell_per_gene[
                        atlas_metacell_index_per_neighborhood_query_cell[neighborhood_query_cell_position],
                        atlas_gene_index,
                    ] > 0
                        all_zero_neighborhood_projected_metacell_UMIs = false
                        break
                    end
                end
                if all_zero_neighborhood_projected_metacell_UMIs
                    atomic_add!(todox_all_zero, Int32(1))
                    correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block[
                        query_gene_index,
                        atlas_block_index,
                    ] = 0.0f0
                    return nothing
                end

                for (neighborhood_query_cell_position, query_cell_index) in
                    enumerate(indices_of_neighborhood_query_cells)
                    cell_UMIs = UMIs_per_query_cell_per_gene[query_cell_index, query_gene_index]
                    atlas_metacell_index =
                        atlas_metacell_index_per_neighborhood_query_cell[neighborhood_query_cell_position]
                    projected_metacell_UMIs = UMIs_per_atlas_metacell_per_gene[atlas_metacell_index, atlas_gene_index]
                    cell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position] =
                        cell_UMIs / total_UMIs_per_query_cell[query_cell_index] + gene_fraction_regularization
                    projected_metacell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position] =
                        projected_metacell_UMIs
                end
                @assert LoopVectorization.check_args(cell_log_fraction_per_neighborhood_query_cell) "check_args failed in compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!\nfor cell_log_fraction_per_neighborhood_query_cell: $(brief(cell_log_fraction_per_neighborhood_query_cell))"
                @assert LoopVectorization.check_args(projected_metacell_log_fraction_per_neighborhood_query_cell) "check_args failed in compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!\nfor projected_metacell_log_fraction_per_neighborhood_query_cell: $(brief(projected_metacell_log_fraction_per_neighborhood_query_cell))"
                @assert LoopVectorization.check_args(total_atlas_metacell_UMIs_per_neighborhood_query_cell) "check_args failed in compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!\nfor total_atlas_metacell_UMIs_per_neighborhood_query_cell: $(brief(total_atlas_metacell_UMIs_per_neighborhood_query_cell))"
                @turbo for neighborhood_query_cell_position in 1:n_neighborhood_query_cells
                    cell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position] =
                        log2(cell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position])
                    projected_metacell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position] =
                        log2(
                            projected_metacell_log_fraction_per_neighborhood_query_cell[neighborhood_query_cell_position] /
                            total_atlas_metacell_UMIs_per_neighborhood_query_cell[neighborhood_query_cell_position] +
                            gene_fraction_regularization,
                        )
                end

                correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block[
                    query_gene_index,
                    atlas_block_index,
                ] = zero_cor_between_vectors(
                    cell_log_fraction_per_neighborhood_query_cell,
                    projected_metacell_log_fraction_per_neighborhood_query_cell,
                )

                return nothing
            end

            if mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing &&
               in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
                is_neighborhood_pertinent_marker_per_common_included_gene = nothing
                is_in_bin_neighborhood_pertinent_marker_per_common_included_gene = nothing
                is_out_bin_neighborhood_pertinent_marker_per_common_included_gene = nothing
            else
                is_neighborhood_pertinent_marker_per_common_included_gene =
                    is_neighborhood_pertinent_marker_per_common_included_gene_per_thread[threadid()]
                is_neighborhood_pertinent_marker_per_common_included_gene .=
                    is_neighborhood_marker_per_atlas_gene_per_atlas_block[
                        indices_of_common_included_atlas_genes,
                        atlas_block_index,
                    ] .& .! is_lateral_per_common_included_gene
                @assert any(is_neighborhood_pertinent_marker_per_common_included_gene)

                if in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
                    is_in_bin_neighborhood_pertinent_marker_per_common_included_gene = nothing
                    is_out_bin_neighborhood_pertinent_marker_per_common_included_gene = nothing
                else
                    is_in_bin_neighborhood_pertinent_marker_per_common_included_gene =
                        is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread[threadid()]
                    is_out_bin_neighborhood_pertinent_marker_per_common_included_gene =
                        is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread[threadid()]

                    @. is_in_bin_neighborhood_pertinent_marker_per_common_included_gene =
                        is_neighborhood_pertinent_marker_per_common_included_gene & is_in_bin_per_common_included_gene
                    @. is_out_bin_neighborhood_pertinent_marker_per_common_included_gene =
                        is_neighborhood_pertinent_marker_per_common_included_gene & is_in_bin_per_common_included_gene
                    @assert any(is_in_bin_neighborhood_pertinent_marker_per_common_included_gene)
                    @assert any(is_out_bin_neighborhood_pertinent_marker_per_common_included_gene)
                end
            end

            for (mean_per_atlas_block, is_in_mask_per_common_included_genes) in (
                (mean_pertinent_marker_genes_correlation_per_atlas_block, is_pertinent_marker_per_common_included_gene),
                (
                    mean_neighborhood_pertinent_markers_correlation_per_atlas_block,
                    is_neighborhood_pertinent_marker_per_common_included_gene,
                ),
                (
                    in_bin_mean_pertinent_markers_correlation_per_atlas_block,
                    is_in_bin_pertinent_marker_per_common_included_gene,
                ),
                (
                    out_bin_mean_pertinent_markers_correlation_per_atlas_block,
                    is_out_bin_pertinent_marker_per_common_included_gene,
                ),
                (
                    in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block,
                    is_in_bin_neighborhood_pertinent_marker_per_common_included_gene,
                ),
                (
                    out_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block,
                    is_out_bin_neighborhood_pertinent_marker_per_common_included_gene,
                ),
            )
                if mean_per_atlas_block !== nothing
                    @assert is_in_mask_per_common_included_genes !== nothing
                    sum_correlations = 0.0
                    num_correlations = 0
                    for common_included_gene_position in 1:n_common_included_genes
                        if is_in_mask_per_common_included_genes[common_included_gene_position]
                            sum_correlations +=
                                correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block[
                                    indices_of_common_included_query_genes[common_included_gene_position],
                                    atlas_block_index,
                                ]
                            num_correlations += 1
                        end
                    end
                    @assert num_correlations > 0
                    mean_per_atlas_block[atlas_block_index] = sum_correlations / num_correlations
                end
            end
        end
    end
    @debug "TODOX ALL-ZERO: $(todox_all_zero[]) OUT OF: $(n_atlas_blocks * n_common_included_genes) ($(percent(todox_all_zero[], n_atlas_blocks * n_common_included_genes)))" _group =
        :todox

    set_matrix!(
        query_daf,
        "gene",
        "projected_block",
        "correlation_between_neighborhood_cells_and_projected_metacells",
        bestify(correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block);
        overwrite,
    )

    if query_daf isa DataAxesFormats.Contracts.ContractDaf
        query_daf = query_daf.daf
    end
    if atlas_daf isa DataAxesFormats.Contracts.ContractDaf
        atlas_daf = atlas_daf.daf
    end
    if query_daf === atlas_daf
        first_qualifier = "self"
        projected_qualifier = " self"
    else
        first_qualifier = "common"
        if bin === nothing
            projected_qualifier = ""
        else
            projected_qualifier = " cross"
        end
    end

    for (mean_per_atlas_block, second_qualifier) in (
        (mean_pertinent_marker_genes_correlation_per_atlas_block, ""),
        (mean_neighborhood_pertinent_markers_correlation_per_atlas_block, " neighborhood"),
        (in_bin_mean_pertinent_markers_correlation_per_atlas_block, " bin"),
        (out_bin_mean_pertinent_markers_correlation_per_atlas_block, " !bin"),
        (in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block, " bin neighborhood"),
        (out_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block, " !bin neighborhood"),
    )
        if mean_per_atlas_block !== nothing
            @debug (
                "Mean correlation of $(first_qualifier)$(second_qualifier) pertinent marker genes between neighborhood cells and their$(projected_qualifier) projected metacells: " *
                "$(mean(mean_per_atlas_block[mean_per_atlas_block .!= 0]))"  # NOLINT
            ) _group = :mcs_results
        end
    end

    return nothing
end

end  # module
