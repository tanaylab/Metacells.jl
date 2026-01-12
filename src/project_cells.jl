"""
Project cells on an atlas with modules.
"""
module ProjectCells

export compute_blocks_cells_pertinent_markers_distances!
export compute_blocks_cells_pertinent_markers_significance!
export compute_blocks_cells_skeleton_distances!
export compute_cell_closest_blocks_by_pertinent_markers!
export compute_cell_closest_blocks_by_skeletons!
export compute_cell_most_correlated_blocks_by_markers!
export compute_cell_most_correlated_blocks_by_skeletons!
export provisional_cells_projection!
export final_cells_projection!

using Base.Threads
using DataAxesFormats
using Distances
using LinearAlgebra
using StatsBase
using TanayLabUtilities

using ..Contracts

#TODOX import Base.Threads.Mutex

"""
TODOX
"""
@logged @computation Contract(;
    axes = [metacell_axis(OptionalInput), block_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(OptionalInput),
        block_gene_linear_fraction_matrix(RequiredInput),
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), block_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        block_cell_pertinent_markers_eucildean_distance_matrix(GuaranteedOutput),
        block_cell_pertinent_markers_correlation_matrix(GuaranteedOutput),
    ],
) function compute_blocks_cells_pertinent_markers_distances!(
    blocks_daf::DafReader,
    cells_daf::DafWriter;
    gene_UMIs_regularization::Real = 1 / 7,
    overwrite::Bool = false,
)::Nothing
    todox_distances(;
        blocks_daf,
        cells_daf,
        gene_UMIs_regularization,
        prefix = "pertinent_",
        which = "marker",
        overwrite,
    )
    return nothing
end

function todox_distances(;
    blocks_daf::DafReader,
    cells_daf::DafWriter,
    gene_UMIs_regularization::Real,
    which::AbstractString,
    prefix::AbstractString,
    overwrite::Bool,
)::Nothing
    @assert axis_vector(cells_daf, "block") == axis_vector(blocks_daf, "block")

    n_blocks = axis_length(blocks_daf, "block")
    n_cells = axis_length(cells_daf, "cell")
    n_cells_genes = axis_length(cells_daf, "gene")

    name_per_which = blocks_daf["/ gene & is_$(which) &! is_lateral : name"].array
    cells_gene_index_per_which = axis_indices(cells_daf, "gene", name_per_which)
    is_block_which_per_cells_gene = zeros(Bool, n_cells_genes)
    is_block_which_per_cells_gene[cells_gene_index_per_which] .= true

    total_UMIs_per_cell = get_vector(cells_daf, "cell", "total_UMIs").array
    UMIs_per_cell_per_gene = get_matrix(cells_daf, "cell", "gene", "UMIs").array
    @views UMIs_per_cell_per_which = UMIs_per_cell_per_gene[:, cells_gene_index_per_which]
    @debug "TODOX prepare..."
    UMIs_per_which_per_cell = densify(flip(UMIs_per_cell_per_which))

    @debug "TODOX log_fraction_per_which_per_cell..."
    log_fraction_per_which_per_cell = log2.(1e-4 .+ UMIs_per_which_per_cell ./ transpose(total_UMIs_per_cell))  # TODOX
    log_fraction_per_which_per_block =
        log2.(1e-4 .+ densify(blocks_daf["/ gene & is_$(which) &! is_lateral / block : linear_fraction"].array))

    @debug "TODOX distances..."
    distances_between_cells_and_blocks =
        parallel_pairwise(Euclidean(), log_fraction_per_which_per_cell, log_fraction_per_which_per_block; dims = 2)
    set_matrix!(cells_daf, "cell", "block", "$(prefix)$(which)s_euclidean_distance", distances_between_cells_and_blocks)

    @debug "TODOX correlations..."
    correlations_between_cells_and_blocks = cor(log_fraction_per_which_per_cell, log_fraction_per_which_per_block)
    set_matrix!(cells_daf, "cell", "block", "$(prefix)$(which)s_correlation", correlations_between_cells_and_blocks)

    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [metacell_axis(OptionalInput), block_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(OptionalInput),
        block_gene_linear_fraction_matrix(RequiredInput),
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        block_cell_skeleton_euclidean_distance_matrix(GuaranteedOutput),
        block_cell_skeleton_correlation_matrix(GuaranteedOutput),
    ],
) function compute_blocks_cells_skeleton_distances!(
    blocks_daf::DafReader,
    cells_daf::DafWriter;
    gene_UMIs_regularization::Real = 1 / 7,
    overwrite::Bool = false,
)::Nothing
    todox_distances(; blocks_daf, cells_daf, gene_UMIs_regularization, prefix = "", which = "skeleton", overwrite)
    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [cell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_cell_pertinent_markers_eucildean_distance_matrix(RequiredInput),
        block_mean_pertinent_markers_distance_vector(GuaranteedOutput),
        block_std_pertinent_markers_distance_vector(GuaranteedOutput),
    ],
) function compute_blocks_cells_pertinent_markers_significance!(daf::DafReader; overwrite::Bool = false)::Nothing
    n_blocks = axis_length(daf, "block")

    distance_per_cell_per_block = get_matrix(daf, "cell", "block", "pertinent_markers_euclidean_distance").array

    mean_distance_per_block = vec(mean(distance_per_cell_per_block; dims = 1))
    @assert_vector(mean_distance_per_block, n_blocks)
    set_vector!(daf, "block", "mean_pertinent_markers_distance", mean_distance_per_block)

    std_distance_per_block = vec(std(distance_per_cell_per_block; dims = 1))
    @assert_vector(std_distance_per_block, n_blocks)
    set_vector!(daf, "block", "std_pertinent_markers_distance", std_distance_per_block)

    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        block_cell_skeleton_euclidean_distance_matrix(RequiredInput),
        cell_closest_by_skeletons_block_vector(GuaranteedOutput),
    ],
) function compute_cell_closest_blocks_by_skeletons!(daf::DafWriter; overwrite::Bool = false)::Nothing
    todox_nearest(daf, "skeletons_euclidean_distance", "closest_by_skeletons", overwrite)
    return nothing
end

function todox_nearest(daf::DafWriter, by::AbstractString, property::AbstractString, overwrite::Bool)::Nothing
    n_cells = axis_length(daf, "cell")
    name_per_block = axis_vector(daf, "block")
    distances_between_blocks_and_cells = get_matrix(daf, "block", "cell", by).array
    nearest_block_index_per_cell = vec(argmin(distances_between_blocks_and_cells; dims = 1))
    @assert_vector(nearest_block_index_per_cell, n_cells)
    nearest_block_name_per_cell = name_per_block[first.(Tuple.(nearest_block_index_per_cell))]
    set_vector!(daf, "cell", "block.$(property)", nearest_block_name_per_cell; overwrite)
    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        block_cell_pertinent_markers_eucildean_distance_matrix(RequiredInput),
        cell_closest_by_pertinent_markers_block_vector(GuaranteedOutput),
    ],
) function compute_cell_closest_blocks_by_pertinent_markers!(daf::DafWriter; overwrite::Bool = false)::Nothing
    todox_nearest(daf, "pertinent_markers_euclidean_distance", "closest_by_pertinent_markers", overwrite)
    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        block_cell_skeleton_correlation_matrix(RequiredInput),
        cell_most_correlated_by_skeletons_block_vector(GuaranteedOutput),
    ],
) function compute_cell_most_correlated_blocks_by_skeletons!(daf::DafWriter; overwrite::Bool = false)::Nothing
    todox_nearest(daf, "skeletons_correlation", "most_correlated_by_skeletons", overwrite)
    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        block_cell_pertinent_markers_correlation_matrix(RequiredInput),
        cell_most_correlated_by_pertinent_markers_block_vector(GuaranteedOutput),
    ],
) function compute_cell_most_correlated_blocks_by_markers!(daf::DafWriter; overwrite::Bool = false)::Nothing
    todox_nearest(daf, "pertinent_markers_correlation", "most_correlated_by_pertinent_markers", overwrite)
    return nothing
end

@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        block_gene_linear_fraction_matrix(RequiredInput),
        block_mean_pertinent_markers_distance_vector(RequiredInput),
        block_std_pertinent_markers_distance_vector(RequiredInput),
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_provisional_block_vector(GuaranteedOutput),
        cell_provisional_block_pertinent_markers_z_score(GuaranteedOutput),
    ],
) function provisional_cells_projection!(
    modules_daf::DafReader,
    cells_daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX coherent regularization
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(modules_daf, "metacell")
    name_per_metacell = axis_vector(modules_daf, "metacell")

    n_blocks = axis_length(modules_daf, "block")
    name_per_block = axis_vector(modules_daf, "block")

    mean_pertinent_markers_distance_per_block =
        get_vector(modules_daf, "block", "mean_pertinent_markers_distance").array
    std_pertinent_markers_distance_per_block = get_vector(modules_daf, "block", "std_pertinent_markers_distance").array

    linear_fraction_per_pertinent_marker_per_block =
        modules_daf["/ gene & is_marker &! is_lateral / block : linear_fraction"].array
    log_linear_fraction_per_pertinent_marker_per_block =
        log2.(gene_cell_fraction_regularization .+ linear_fraction_per_pertinent_marker_per_block)

    name_per_pertinent_marker = modules_daf["/ gene & is_marker &! is_lateral"]
    index_per_pertinent_marker = axis_indices(cells_daf, "gene", name_per_pertinent_marker)

    n_cells = axis_length(cells_daf, "cell")

    total_UMIs_per_cell = get_vector(cells_daf, "cell", "total_UMIs").array
    UMIs_per_cell_per_gene = get_matrix(cells_daf, "cell", "gene", "UMIs").array

    provisional_block_per_cell = Vector{AbstractString}(undef, n_cells)
    z_score_per_cell = Vector{Float32}(undef, n_cells)

    progress_counter = Atomic{Int}(0)
    @threads :greedy for cell_index in 1:n_cells
        UMIs_per_pertinent_marker = densify(UMIs_per_cell_per_gene[cell_index, index_per_pertinent_marker])
        log_linear_fraction_per_pertinent_marker =
            log2.(gene_cell_fraction_regularization .+ UMIs_per_pertinent_marker ./ total_UMIs_per_cell[cell_index])

        distances_to_blocks = vec(
            pairwise(
                Euclidean(),
                Ref(log_linear_fraction_per_pertinent_marker),
                eachcol(log_linear_fraction_per_pertinent_marker_per_block),
            ),
        )
        @assert_vector(distances_to_blocks, n_blocks)

        provisional_block_index = argmin(distances_to_blocks)
        provisional_block_per_cell[cell_index] = name_per_block[provisional_block_index]

        minimal_distance = distances_to_blocks[provisional_block_index]
        mean_pertinent_markers_distance = mean_pertinent_markers_distance_per_block[provisional_block_index]
        std_pertinent_markers_distance = std_pertinent_markers_distance_per_block[provisional_block_index]
        z_score_per_cell[cell_index] =
            (minimal_distance - mean_pertinent_markers_distance) / std_pertinent_markers_distance
        counter = atomic_add!(progress_counter, 1)
        if cell_index % 100 == 1
            print("\r$(counter) ($(percent(counter, n_cells)))   ")
        end
    end

    set_vector!(cells_daf, "cell", "block.provisional", provisional_block_per_cell; overwrite)
    set_vector!(cells_daf, "cell", "provisional_block_pertinent_markers_z_score", z_score_per_cell; overwrite)

    return nothing
end

@logged @computation Contract(;
    axes = [
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_module_is_strong_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_metacell_module_linear_fraction_tensor(RequiredInput),
        metacell_mean_modules_distance_vector(RequiredInput),
        metacell_std_modules_distance_vector(RequiredInput),
        block_module_neighborhood_std_linear_fraction_matrix(RequiredInput),  # TODOX: Normalized only
        block_module_neighborhood_mean_linear_fraction_matrix(RequiredInput),  # TODOX: Normalized only
    ],
) Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_provisional_block_vector(RequiredInput),
        cell_provisional_block_pertinent_markers_z_score(RequiredInput),
        cell_projected_block_vector(GuaranteedOutput),
        cell_projected_metacell_vector(GuaranteedOutput),
        cell_projected_metacell_modules_z_score_vector(GuaranteedOutput),
    ],
) function final_cells_projection!(
    modules_daf::DafReader,
    cells_daf::DafWriter;
    normalize_fractions::Bool = false,
    overwrite::Bool = false,
)::Nothing
    n_cells = axis_length(cells_daf, "cell")
    provisional_block_per_cell = get_vector(cells_daf, "cell", "block.provisional").array
    provisional_markers_z_score_per_cell =
        get_vector(cells_daf, "cell", "provisional_block_pertinent_markers_z_score").array

    UMIs_per_cell_per_gene = get_matrix(cells_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(cells_daf, "cell", "total_UMIs").array

    n_blocks = axis_length(modules_daf, "block")
    name_per_block = axis_vector(modules_daf, "block")
    name_per_module = axis_vector(modules_daf, "module")
    name_per_metacell = axis_vector(modules_daf, "metacell")

    module_index_per_gene_per_block = modules_daf["/ gene / block : module ?? 0 => index"].array
    is_strong_per_module_per_block = get_matrix(modules_daf, "module", "block", "is_strong").array

    provisional_block_index_per_cell = axis_indices(modules_daf, "block", provisional_block_per_cell)

    projected_block_per_cell = copy_array(provisional_block_per_cell; eltype = AbstractString)

    mean_linear_fraction_per_module_per_block =
        get_matrix(modules_daf, "module", "block", "neighborhood_mean_linear_fraction").array
    std_linear_fraction_per_module_per_block =
        get_matrix(modules_daf, "module", "block", "neighborhood_std_linear_fraction").array
    module_per_gene_per_block = get_matrix(modules_daf, "gene", "block", "module").array

    mean_modules_distance_per_metacell = get_vector(modules_daf, "metacell", "mean_modules_distance").array
    std_modules_distance_per_metacell = get_vector(modules_daf, "metacell", "std_modules_distance").array

    projected_block_per_cell = fill("", n_cells)
    projected_metacell_per_cell = fill("", n_cells)
    modules_z_score_per_cell = zeros(Float32, n_cells)

    next_indices_of_cells_per_block =
        [findall(provisional_block_index_per_cell .== block_index) for block_index in 1:n_blocks]

    spin_lock = SpinLock()

    determined_counter = Atomic{Int}(0)

    for phase in 1:3
        indices_of_cells_per_block = next_indices_of_cells_per_block
        next_indices_of_cells_per_block = [Int32[] for _ in 1:n_blocks]

        @debug "TODOX Determined: $(determined_counter[]) $(percent(determined_counter[], n_cells)) Phase $(phase)..."
        @threads :greedy for block_index in 1:n_blocks
            block_name = name_per_block[block_index]

            indices_of_block_cells = indices_of_cells_per_block[block_index]
            n_block_cells = length(indices_of_block_cells)
            if n_block_cells > 0
                if 1 <= phase <= 2
                    indices_of_candidate_metacells =
                        modules_daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array
                    n_candidate_metacells = length(indices_of_candidate_metacells)
                    block_index_per_candidate_metacell =
                        modules_daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : block => index"].array
                elseif phase == 3
                    indices_of_candidate_metacells = modules_daf["/ metacell & block = $(block_name) : index"].array
                    n_candidate_metacells = length(indices_of_candidate_metacells)
                    block_index_per_candidate_metacell = nothing
                else
                    @assert false
                end

                @views module_index_per_gene = module_index_per_gene_per_block[:, block_index]
                @views is_strong_per_module = is_strong_per_module_per_block[:, block_index]
                indices_of_strong_modules = findall(is_strong_per_module)
                n_strong_modules = length(indices_of_strong_modules)
                @assert n_strong_modules > 0

                linear_fraction_per_module_per_metacell =
                    get_matrix(modules_daf, "module", "metacell", "$(block_name)_linear_fraction").array
                linear_fraction_per_strong_module_per_candidate_metacell =
                    linear_fraction_per_module_per_metacell[indices_of_strong_modules, indices_of_candidate_metacells]

                @views module_per_gene = module_per_gene_per_block[:, block_index]
                total_UMIs_per_block_cell = total_UMIs_per_cell[indices_of_block_cells]
                linear_fraction_per_strong_module_per_block_cell =
                    Matrix{Float32}(undef, n_strong_modules, n_block_cells)

                for (strong_module_position, strong_module_index) in enumerate(indices_of_strong_modules)
                    strong_module_name = name_per_module[strong_module_index]
                    indices_of_strong_module_genes = findall(module_per_gene .== strong_module_name)
                    n_strong_module_genes = length(indices_of_strong_module_genes)
                    @assert n_strong_module_genes > 0
                    strong_module_UMIs_per_block_cell = vec(
                        sum(UMIs_per_cell_per_gene[indices_of_block_cells, indices_of_strong_module_genes]; dims = 2),
                    )
                    @assert_vector(strong_module_UMIs_per_block_cell, n_block_cells)
                    linear_fraction_per_strong_module_per_block_cell[strong_module_position, :] .=
                        strong_module_UMIs_per_block_cell ./ total_UMIs_per_block_cell
                end

                if normalize_fractions
                    mean_linear_fraction_per_strong_module =
                        mean_linear_fraction_per_module_per_block[indices_of_strong_modules, block_index]
                    std_linear_fraction_per_strong_module =
                        std_linear_fraction_per_module_per_block[indices_of_strong_modules, block_index]
                    linear_fraction_per_strong_module_per_candidate_metacell .=
                        (
                            linear_fraction_per_strong_module_per_candidate_metacell .-
                            mean_linear_fraction_per_strong_module
                        ) ./ std_linear_fraction_per_strong_module
                    linear_fraction_per_strong_module_per_block_cell .=
                        (linear_fraction_per_strong_module_per_block_cell .- mean_linear_fraction_per_strong_module) ./
                        std_linear_fraction_per_strong_module
                end

                distances_between_candidate_metacells_and_block_cells = pairwise(
                    Euclidean(),
                    linear_fraction_per_strong_module_per_candidate_metacell,
                    linear_fraction_per_strong_module_per_block_cell,
                )

                minimal_position_per_block_cell =
                    vec(argmin(distances_between_candidate_metacells_and_block_cells; dims = 1))
                minimal_distance_per_block_cell =
                    distances_between_candidate_metacells_and_block_cells[minimal_position_per_block_cell]
                position_of_nearest_candidate_metacell_per_block_cell =
                    [position.I[1] for position in minimal_position_per_block_cell]

                if block_index_per_candidate_metacell === nothing
                    positions_of_stable_cells = 1:n_block_cells
                else
                    nearest_block_index_per_block_cell =
                        block_index_per_candidate_metacell[position_of_nearest_candidate_metacell_per_block_cell]
                    positions_of_stable_cells = findall(nearest_block_index_per_block_cell .== block_index)

                    try
                        lock(spin_lock)
                        for (cell_position, nearest_block_index) in enumerate(nearest_block_index_per_block_cell)
                            if nearest_block_index != block_index
                                cell_index = indices_of_block_cells[cell_position]
                                if phase == 1
                                    push!(next_indices_of_cells_per_block[nearest_block_index], cell_index)
                                elseif phase == 2
                                    provisional_block_index = provisional_block_index_per_cell[cell_index]
                                    @assert block_index != provisional_block_index
                                    push!(next_indices_of_cells_per_block[provisional_block_index], cell_index)
                                else
                                    @assert false
                                end
                            end
                        end
                    finally
                        unlock(spin_lock)
                    end
                end

                n_stable_cells = length(positions_of_stable_cells)
                if n_stable_cells > 0
                    metacell_index_per_stable_cell =
                        indices_of_candidate_metacells[position_of_nearest_candidate_metacell_per_block_cell[positions_of_stable_cells]]
                    metacell_per_stable_cell = name_per_metacell[metacell_index_per_stable_cell]

                    indices_of_stable_cells = indices_of_block_cells[positions_of_stable_cells]
                    projected_metacell_per_cell[indices_of_stable_cells] = metacell_per_stable_cell

                    modules_z_score_per_cell[indices_of_stable_cells] .=
                        (
                            minimal_distance_per_block_cell[positions_of_stable_cells] .-
                            mean_modules_distance_per_metacell[metacell_index_per_stable_cell]
                        ) ./ (std_modules_distance_per_metacell[metacell_index_per_stable_cell])
                end
                atomic_add!(determined_counter, n_stable_cells)
            end
        end
    end

    @assert determined_counter[] == n_cells

    set_vector!(cells_daf, "cell", "metacell.projected", projected_metacell_per_cell; overwrite)
    set_vector!(cells_daf, "cell", "projected_metacell_modules_z_score", modules_z_score_per_cell; overwrite)

    block_per_metacell = get_vector(modules_daf, "metacell", "block").array
    projected_metacell_index_per_cell = axis_indices(modules_daf, "metacell", projected_metacell_per_cell)
    projected_block_per_cell = block_per_metacell[projected_metacell_index_per_cell]
    set_vector!(cells_daf, "cell", "block.projected", projected_block_per_cell; overwrite)

    return nothing
end

"""
TODOX
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_cell_skeleton_euclidean_distance_matrix(RequiredInput),
    ],
) function compute_cell_projected_blocks_correlations!(daf::DafWriter; overwrite::Bool = false)::Nothing
    block_index_per_cell = daf["/ cell : metacell.projected ?? 0 => block => index"].array
    correlation_between_blocks_and_cells = get_matrix(daf, "block", "cell", "correlation").array
    self_correlation_per_cell = [
        block_index == 0 ? 0 : correlation_between_blocks_and_cells[block_index, cell_index] for
        (cell_index, block_index) in enumerate(block_index_per_cell)
    ]
    set_vector!(daf, "cell", "projected_block_correlation", self_correlation_per_cell; overwrite = true)
    return nothing
end

end  # module
