"""
Project cells on an atlas with modules.
"""
module ProjectCells

export compute_cells_projection!
export compute_vector_of_is_deviant_per_cell!
export compute_vector_of_correlation_between_cells_and_projected_metacells!
export compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!

using Base.Threads
using DataAxesFormats
using LinearAlgebra
using LoopVectorization
using ProgressMeter
using SparseArrays
using StatsBase
using TanayLabUtilities

using ..AnalyzeCells
using ..AnalyzeModules
using ..Contracts
using ..Defaults

import Base.Threads.maxthreadid
import Metacells.AnalyzeCells.cell_group_membership_from_block_neighborhoods
import Metacells.AnalyzeCells.sparse_cell_reference_correlation_per_gene_per_group!
import Metacells.AnalyzeModules.compute_std_fraction_regularization_per_block

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
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.module_axis
import Metacells.Contracts.projected_block_axis
import Metacells.Contracts.tensor_of_mean_cells_z_score_per_block_per_module_per_metacell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_correlation_between_cells_and_projected_metacells_per_gene
import Metacells.Contracts.vector_of_is_deviant_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_mean_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_projected_block_per_cell
import Metacells.Contracts.vector_of_projected_metacell_per_cell
import Metacells.Contracts.vector_of_projected_modules_z_score_per_cell
import Metacells.Contracts.vector_of_std_euclidean_modules_cells_distance_per_metacell
import Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell

"""
    compute_cells_projection!(;
        query_daf::DafWriter,
        atlas_daf::DafReader,
        gene_fraction_regularization::Real = $(DEFAULT.gene_fraction_regularization),
        std_UMIs_regularization::AbstractFloat = $(DEFAULT.std_UMIs_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_projected_metacell_per_cell`](@ref), [`vector_of_projected_block_per_cell`](@ref) and
[`vector_of_projected_modules_z_score_per_cell`](@ref). To pick the best metacell in the `atlas_daf` for each cell of
the `query_daf`, we:

  - Pick a provisional block for each cell. This is the block with the minimal Euclidean distance of the non-lateral
    marker genes. It is the cell's first current block.
  - Score the cell in the found modules of its current block, and find the metacell of the block's neighborhood with
    the closest Euclidean distance.
  - If this metacell is in the current block, we accept it.
  - If it is in a block the cell hasn't visited yet, that block becomes the current block, and we repeat.
  - Otherwise, it is in a block the cell already visited. We settle for the closest metacell of the current block.

Each cell is projected on its own, reading only its own column of the query's `UMIs`. The cells are read in a single
pass over the gene × cell layout.

The provisional block is picked using the log (base 2) of the expression of the marker genes, using the
`gene_fraction_regularization` to handle zero fractions. The metacells are compared using the z-scores of the found
modules, computed the same way sharpening computes them, using the same `std_UMIs_regularization` (see
[`compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!`](@ref)). Each metacell is represented by the
mean z-scores of its cells. We also compute the z-score of the final distance between the cell and its projected
metacell, relative to the mean and standard deviation of the distances of the metacell's own cells (see
[`compute_stats_of_euclidean_modules_cells_distance_per_metacell!`](@ref)). If this is too large, then the query cell
isn't a good match for anything in the atlas.

Excluded query cells are not projected. Their projected metacell and block are empty, and their z-score is infinite.

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
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(RequiredInput),
        tensor_of_mean_cells_z_score_per_block_per_module_per_metacell(RequiredInput),
        matrix_of_is_in_neighborhood_per_block_per_block(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_block(RequiredInput),
        vector_of_block_per_metacell(RequiredInput),
        vector_of_total_neighborhood_UMIs_per_block(RequiredInput),
        vector_of_n_neighborhood_cells_per_block(RequiredInput),
        vector_of_mean_euclidean_modules_cells_distance_per_metacell(RequiredInput),
        vector_of_std_euclidean_modules_cells_distance_per_metacell(RequiredInput),
    ],
) function compute_cells_projection!(;  # UNTESTED
    query_daf::DafWriter,
    atlas_daf::DafReader,
    gene_fraction_regularization::Real = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    std_UMIs_regularization::AbstractFloat = function_default(
        compute_tensor_of_mean_cells_z_score_per_block_per_module_per_metacell!,
        :std_UMIs_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    query_gene_index_per_atlas_pertinent_marker, query_gene_index_per_atlas_gene =
        map_query_atlas_genes(; query_daf, atlas_daf)
    projection_atlas = load_projection_atlas(;
        query_daf,
        atlas_daf,
        query_gene_index_per_atlas_pertinent_marker,
        query_gene_index_per_atlas_gene,
        gene_fraction_regularization,
        std_UMIs_regularization,
    )

    is_excluded_per_query_cell = plain_array(get_vector(query_daf, "cell", "is_excluded"))
    total_UMIs_per_query_cell = plain_array(get_vector(query_daf, "cell", "total_UMIs"))
    UMIs_per_query_gene_per_query_cell = get_matrix(query_daf, "gene", "cell", "UMIs")
    @assert issparse(UMIs_per_query_gene_per_query_cell)

    n_query_cells = length(is_excluded_per_query_cell)
    indices_of_included_query_cells = findall(.!is_excluded_per_query_cell)
    projected_metacell_index_per_query_cell = zeros(UInt32, n_query_cells)
    projected_modules_z_score_per_query_cell = fill(Inf32, n_query_cells)

    n_blocks = length(projection_atlas.squared_norm_per_block)
    max_n_found_modules = maximum(projection_atlas.n_found_modules_per_block)
    score_per_block_per_thread = [Vector{Float64}(undef, n_blocks) for _ in 1:maxthreadid()]
    module_UMIs_per_found_module_per_thread = [Vector{Float32}(undef, max_n_found_modules) for _ in 1:maxthreadid()]
    z_score_per_found_module_per_thread = [Vector{Float32}(undef, max_n_found_modules) for _ in 1:maxthreadid()]
    visited_block_indices_per_thread = [Int[] for _ in 1:maxthreadid()]

    parallel_foreach_column(
        UMIs_per_query_gene_per_query_cell;
        columns = indices_of_included_query_cells,
        name = "cells_projection",
        progress = DebugProgress(
            length(indices_of_included_query_cells);
            group = :mcs_loops,
            desc = "cells_projection",
        ),
        progress_chunk = 100,
    ) do _, query_cell_index, UMIs_per_query_gene
        thread = threadid()
        query_gene_index_per_stored = SparseArrays.nonzeroinds(UMIs_per_query_gene)
        UMIs_per_stored = nonzeros(UMIs_per_query_gene)
        total_UMIs = total_UMIs_per_query_cell[query_cell_index]

        provisional_block_index = provisional_block_index_of_query_cell(
            projection_atlas;
            query_gene_index_per_stored,
            UMIs_per_stored,
            total_UMIs,
            gene_fraction_regularization,
            score_per_block = score_per_block_per_thread[thread],
        )
        metacell_index, distance = projected_metacell_of_query_cell(
            projection_atlas;
            provisional_block_index,
            query_gene_index_per_stored,
            UMIs_per_stored,
            total_UMIs,
            module_UMIs_per_found_module = module_UMIs_per_found_module_per_thread[thread],
            z_score_per_found_module = z_score_per_found_module_per_thread[thread],
            visited_block_indices = visited_block_indices_per_thread[thread],
        )

        projected_metacell_index_per_query_cell[query_cell_index] = metacell_index
        projected_modules_z_score_per_query_cell[query_cell_index] =
            (distance - projection_atlas.mean_euclidean_modules_cells_distance_per_metacell[metacell_index]) /
            projection_atlas.std_euclidean_modules_cells_distance_per_metacell[metacell_index]
        return nothing
    end

    name_per_metacell = axis_vector(atlas_daf, "metacell")
    block_per_metacell = plain_array(get_vector(atlas_daf, "metacell", "block"))
    projected_metacell_per_query_cell = [
        metacell_index == 0 ? "" : name_per_metacell[metacell_index] for
        metacell_index in projected_metacell_index_per_query_cell
    ]
    projected_block_per_query_cell = [
        metacell_index == 0 ? "" : block_per_metacell[metacell_index] for
        metacell_index in projected_metacell_index_per_query_cell
    ]

    set_vector!(query_daf, "cell", "projected_metacell", projected_metacell_per_query_cell; overwrite)
    set_vector!(
        query_daf,
        "cell",
        "projected_metacell_modules_z_score",
        projected_modules_z_score_per_query_cell;
        overwrite,
    )
    set_vector!(query_daf, "cell", "projected_block", projected_block_per_query_cell; overwrite)

    return nothing
end

function map_query_atlas_genes(;  # UNTESTED
    query_daf::DafReader,
    atlas_daf::DafReader,
)::Tuple{AbstractVector{<:Integer}, AbstractVector{<:Integer}}
    atlas_genes = axis_vector(atlas_daf, "gene")

    is_in_module_per_atlas_gene = atlas_daf["@ block @ gene :: module != '' >- Max"].array
    atlas_in_module_genes = atlas_genes[is_in_module_per_atlas_gene]
    axis_indices(query_daf, "gene", atlas_in_module_genes)

    is_pertinent_marker_per_atlas_gene =
        plain_array(get_vector(atlas_daf, "gene", "is_marker")) .&
        .!plain_array(get_vector(atlas_daf, "gene", "is_lateral"))
    atlas_pertinent_marker_genes = atlas_genes[is_pertinent_marker_per_atlas_gene]

    query_gene_index_per_atlas_pertinent_marker = axis_indices(query_daf, "gene", atlas_pertinent_marker_genes)
    query_gene_index_per_atlas_gene = axis_indices(query_daf, "gene", atlas_genes; allow_missing = true)

    return query_gene_index_per_atlas_pertinent_marker, query_gene_index_per_atlas_gene
end

# What projecting a query cell needs from the atlas, indexed by the query's genes.
struct ProjectionAtlas
    # The position of each query gene among the atlas's pertinent markers (0 if it isn't one).
    pertinent_marker_position_per_query_gene::Vector{Int32}
    # The `log2(1 + fraction / regularization)` of each pertinent marker in each block.
    shifted_log_linear_fraction_per_block_per_pertinent_marker::Matrix{Float64}
    squared_norm_per_block::Vector{Float64}
    # The position of each query gene among the atlas genes in some module (0 if it is in none).
    module_gene_position_per_query_gene::Vector{Int32}
    # The position of each module gene's module among the found modules of each block (0 if it isn't found there).
    found_module_position_per_module_gene_per_block::Matrix{UInt16}
    n_found_modules_per_block::Vector{Int}
    mean_linear_fraction_in_environment_cells_per_found_module_per_block::Matrix{Float32}
    std_linear_fraction_in_environment_cells_per_found_module_per_block::Matrix{Float32}
    std_fraction_regularization_per_block::Vector{Float32}
    # The metacells of each block's neighborhood, and their mean z-scores in the block's found modules.
    indices_of_candidate_metacells_per_block::Vector{Vector{Int}}
    mean_cells_z_score_per_found_module_per_candidate_metacell_per_block::Vector{Matrix{Float32}}
    block_index_per_metacell::Vector{Int}
    mean_euclidean_modules_cells_distance_per_metacell::Vector{Float32}
    std_euclidean_modules_cells_distance_per_metacell::Vector{Float32}
end

function load_projection_atlas(;  # UNTESTED
    query_daf::DafReader,
    atlas_daf::DafReader,
    query_gene_index_per_atlas_pertinent_marker::AbstractVector{<:Integer},
    query_gene_index_per_atlas_gene::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat,
    std_UMIs_regularization::AbstractFloat,
)::ProjectionAtlas
    n_query_genes = axis_length(query_daf, "gene")
    n_blocks = axis_length(atlas_daf, "block")
    n_modules = axis_length(atlas_daf, "module")
    name_per_block = axis_vector(atlas_daf, "block")

    pertinent_marker_position_per_query_gene = zeros(Int32, n_query_genes)
    pertinent_marker_position_per_query_gene[query_gene_index_per_atlas_pertinent_marker] .=
        1:length(query_gene_index_per_atlas_pertinent_marker)

    # As in `compute_vector_of_block_closest_by_pertinent_markers_per_cell!`: the distances use
    # `log2(1 + fraction / regularization)`, so zero fractions stay zero and each query cell's work is proportional to
    # its non-zero UMIs, and the `|cell|²` term of the squared distance is skipped as it is the same for all blocks.
    linear_fraction_per_block_per_pertinent_marker =
        atlas_daf["@ block @ gene [ is_marker & ! is_lateral ] :: linear_fraction"].array
    shifted_log_linear_fraction_per_block_per_pertinent_marker = log2.(
        1 .+
        Float64.(Matrix(plain_array(densify(linear_fraction_per_block_per_pertinent_marker)))) ./
        gene_fraction_regularization,
    )
    @check_turbo_matrix(shifted_log_linear_fraction_per_block_per_pertinent_marker)
    squared_norm_per_block = vec(sum(abs2, shifted_log_linear_fraction_per_block_per_pertinent_marker; dims = 2))

    module_index_per_atlas_gene_per_block = get_query(atlas_daf, "@ gene @ block :: module ?? 0 : index").array
    is_found_per_module_per_block = plain_array(get_matrix(atlas_daf, "module", "block", "is_found"))
    is_in_neighborhood_per_other_block_per_base_block =
        plain_array(get_matrix(atlas_daf, "block", "block", "is_in_neighborhood"))
    block_index_per_metacell = Vector{Int}(atlas_daf["@ metacell : block : index"].array)
    mean_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(atlas_daf, "module", "block", "mean_linear_fraction_in_environment_cells"))
    std_linear_fraction_in_environment_cells_per_module_per_block =
        plain_array(get_matrix(atlas_daf, "module", "block", "std_linear_fraction_in_environment_cells"))

    indices_of_module_atlas_genes = findall(vec(any(module_index_per_atlas_gene_per_block .> 0; dims = 2)))
    n_module_genes = length(indices_of_module_atlas_genes)
    module_gene_position_per_query_gene = zeros(Int32, n_query_genes)
    module_gene_position_per_query_gene[query_gene_index_per_atlas_gene[indices_of_module_atlas_genes]] .=
        1:n_module_genes

    n_found_modules_per_block = vec(sum(is_found_per_module_per_block; dims = 1))
    max_n_found_modules = maximum(n_found_modules_per_block)
    found_module_position_per_module_gene_per_block = zeros(UInt16, n_module_genes, n_blocks)
    mean_linear_fraction_in_environment_cells_per_found_module_per_block = zeros(Float32, max_n_found_modules, n_blocks)
    std_linear_fraction_in_environment_cells_per_found_module_per_block = zeros(Float32, max_n_found_modules, n_blocks)
    indices_of_candidate_metacells_per_block = Vector{Vector{Int}}(undef, n_blocks)
    mean_cells_z_score_per_found_module_per_candidate_metacell_per_block = Vector{Matrix{Float32}}(undef, n_blocks)

    found_module_position_per_module = zeros(UInt16, n_modules)
    for block_index in 1:n_blocks
        indices_of_found_modules = findall(view(is_found_per_module_per_block, :, block_index))
        n_found_modules = length(indices_of_found_modules)
        @assert n_found_modules > 0

        found_module_position_per_module .= 0
        found_module_position_per_module[indices_of_found_modules] .= 1:n_found_modules
        for (module_gene_position, atlas_gene_index) in enumerate(indices_of_module_atlas_genes)
            module_index = module_index_per_atlas_gene_per_block[atlas_gene_index, block_index]
            if module_index > 0
                found_module_position_per_module_gene_per_block[module_gene_position, block_index] =
                    found_module_position_per_module[module_index]
            end
        end

        mean_linear_fraction_in_environment_cells_per_found_module_per_block[1:n_found_modules, block_index] .=
            view(mean_linear_fraction_in_environment_cells_per_module_per_block, indices_of_found_modules, block_index)
        std_linear_fraction_in_environment_cells_per_found_module_per_block[1:n_found_modules, block_index] .=
            view(std_linear_fraction_in_environment_cells_per_module_per_block, indices_of_found_modules, block_index)
        @assert all(
            std_linear_fraction_in_environment_cells_per_found_module_per_block[1:n_found_modules, block_index] .> 0,
        )

        # A block is always in its own neighborhood, so its own metacells are among the candidates.
        indices_of_candidate_metacells =
            findall(view(is_in_neighborhood_per_other_block_per_base_block, :, block_index)[block_index_per_metacell])
        indices_of_candidate_metacells_per_block[block_index] = indices_of_candidate_metacells
        mean_cells_z_score_per_module_per_metacell = plain_array(
            get_matrix(atlas_daf, "module", "metacell", "$(name_per_block[block_index])_mean_cells_z_score"),
        )
        mean_cells_z_score_per_found_module_per_candidate_metacell_per_block[block_index] = Matrix{Float32}(
            mean_cells_z_score_per_module_per_metacell[indices_of_found_modules, indices_of_candidate_metacells],
        )
    end

    return ProjectionAtlas(
        pertinent_marker_position_per_query_gene,
        shifted_log_linear_fraction_per_block_per_pertinent_marker,
        squared_norm_per_block,
        module_gene_position_per_query_gene,
        found_module_position_per_module_gene_per_block,
        n_found_modules_per_block,
        mean_linear_fraction_in_environment_cells_per_found_module_per_block,
        std_linear_fraction_in_environment_cells_per_found_module_per_block,
        compute_std_fraction_regularization_per_block(atlas_daf; std_UMIs_regularization),
        indices_of_candidate_metacells_per_block,
        mean_cells_z_score_per_found_module_per_candidate_metacell_per_block,
        block_index_per_metacell,
        Vector{Float32}(plain_array(get_vector(atlas_daf, "metacell", "mean_euclidean_modules_cells_distance"))),
        Vector{Float32}(plain_array(get_vector(atlas_daf, "metacell", "std_euclidean_modules_cells_distance"))),
    )
end

# The block whose pertinent markers are closest to the query cell's.
function provisional_block_index_of_query_cell(  # UNTESTED
    projection_atlas::ProjectionAtlas;
    query_gene_index_per_stored::AbstractVector{<:Integer},
    UMIs_per_stored::AbstractVector{<:Integer},
    total_UMIs::Integer,
    gene_fraction_regularization::AbstractFloat,
    score_per_block::Vector{Float64},
)::Int
    shifted_log_linear_fraction_per_block_per_pertinent_marker =
        projection_atlas.shifted_log_linear_fraction_per_block_per_pertinent_marker
    pertinent_marker_position_per_query_gene = projection_atlas.pertinent_marker_position_per_query_gene
    n_blocks = length(score_per_block)

    score_per_block .= projection_atlas.squared_norm_per_block
    scale = 1.0 / (total_UMIs * gene_fraction_regularization)
    @inbounds for stored_index in eachindex(query_gene_index_per_stored)
        pertinent_marker_position = pertinent_marker_position_per_query_gene[query_gene_index_per_stored[stored_index]]
        if pertinent_marker_position > 0
            twice_shifted_log_linear_fraction = 2 * log2(1 + UMIs_per_stored[stored_index] * scale)
            @turbo for block_index in 1:n_blocks
                score_per_block[block_index] -=
                    twice_shifted_log_linear_fraction *
                    shifted_log_linear_fraction_per_block_per_pertinent_marker[block_index, pertinent_marker_position]
            end
        end
    end

    return argmin(score_per_block)
end

# The z-scores of the query cell in the found modules of a block. These are the z-scores
# `compute_z_score_per_found_module_per_region_cell!` computes, for a single cell.
function compute_z_score_per_found_module_of_query_cell!(  # UNTESTED
    projection_atlas::ProjectionAtlas;
    block_index::Integer,
    query_gene_index_per_stored::AbstractVector{<:Integer},
    UMIs_per_stored::AbstractVector{<:Integer},
    total_UMIs::Integer,
    module_UMIs_per_found_module::AbstractVector{Float32},
    z_score_per_found_module::AbstractVector{Float32},
)::Nothing
    module_gene_position_per_query_gene = projection_atlas.module_gene_position_per_query_gene
    found_module_position_per_module_gene_per_block = projection_atlas.found_module_position_per_module_gene_per_block

    module_UMIs_per_found_module .= 0.0f0
    @inbounds for stored_index in eachindex(query_gene_index_per_stored)
        module_gene_position = module_gene_position_per_query_gene[query_gene_index_per_stored[stored_index]]
        if module_gene_position > 0
            found_module_position = found_module_position_per_module_gene_per_block[module_gene_position, block_index]
            if found_module_position > 0
                module_UMIs_per_found_module[found_module_position] += UMIs_per_stored[stored_index]
            end
        end
    end

    std_fraction_regularization = projection_atlas.std_fraction_regularization_per_block[block_index]
    @inbounds for found_module_position in eachindex(z_score_per_found_module)
        mean_linear_fraction = projection_atlas.mean_linear_fraction_in_environment_cells_per_found_module_per_block[
            found_module_position,
            block_index,
        ]
        std_linear_fraction = projection_atlas.std_linear_fraction_in_environment_cells_per_found_module_per_block[
            found_module_position,
            block_index,
        ]
        delta_fraction = module_UMIs_per_found_module[found_module_position] / total_UMIs - mean_linear_fraction
        z_score_per_found_module[found_module_position] =
            (delta_fraction + copysign(std_fraction_regularization, delta_fraction)) /
            (std_linear_fraction + std_fraction_regularization)
    end

    return nothing
end

# The metacell a query cell is projected on, and its distance from it. Starting at the provisional block, the cell is
# scored in the current block's found modules, against the metacells of the block's neighborhood. If the nearest is in
# the current block, the cell settles on it. If it is in a block the cell hasn't visited, the cell moves there.
# Otherwise, the cell settles on the nearest metacell of the current block.
function projected_metacell_of_query_cell(  # UNTESTED
    projection_atlas::ProjectionAtlas;
    provisional_block_index::Integer,
    query_gene_index_per_stored::AbstractVector{<:Integer},
    UMIs_per_stored::AbstractVector{<:Integer},
    total_UMIs::Integer,
    module_UMIs_per_found_module::Vector{Float32},
    z_score_per_found_module::Vector{Float32},
    visited_block_indices::Vector{Int},
)::Tuple{Int, Float32}
    block_index_per_metacell = projection_atlas.block_index_per_metacell
    block_index = provisional_block_index
    empty!(visited_block_indices)

    while true
        push!(visited_block_indices, block_index)
        n_found_modules = projection_atlas.n_found_modules_per_block[block_index]
        @views block_z_score_per_found_module = z_score_per_found_module[1:n_found_modules]
        compute_z_score_per_found_module_of_query_cell!(
            projection_atlas;
            block_index,
            query_gene_index_per_stored,
            UMIs_per_stored,
            total_UMIs,
            module_UMIs_per_found_module = view(module_UMIs_per_found_module, 1:n_found_modules),
            z_score_per_found_module = block_z_score_per_found_module,
        )

        indices_of_candidate_metacells = projection_atlas.indices_of_candidate_metacells_per_block[block_index]
        mean_cells_z_score_per_found_module_per_candidate_metacell =
            projection_atlas.mean_cells_z_score_per_found_module_per_candidate_metacell_per_block[block_index]
        nearest_position = 0
        nearest_squared_distance = Inf32
        nearest_own_position = 0
        nearest_own_squared_distance = Inf32
        @inbounds for candidate_position in eachindex(indices_of_candidate_metacells)
            squared_distance = 0.0f0
            for found_module_position in 1:n_found_modules
                delta =
                    block_z_score_per_found_module[found_module_position] -
                    mean_cells_z_score_per_found_module_per_candidate_metacell[
                        found_module_position,
                        candidate_position,
                    ]
                squared_distance += delta * delta
            end
            if squared_distance < nearest_squared_distance
                nearest_squared_distance = squared_distance
                nearest_position = candidate_position
            end
            if squared_distance < nearest_own_squared_distance &&
               block_index_per_metacell[indices_of_candidate_metacells[candidate_position]] == block_index
                nearest_own_squared_distance = squared_distance
                nearest_own_position = candidate_position
            end
        end
        @assert nearest_own_position > 0

        nearest_block_index = block_index_per_metacell[indices_of_candidate_metacells[nearest_position]]
        if nearest_block_index == block_index
            return indices_of_candidate_metacells[nearest_position], sqrt(nearest_squared_distance)
        elseif !(nearest_block_index in visited_block_indices)
            block_index = nearest_block_index
        else
            return indices_of_candidate_metacells[nearest_own_position], sqrt(nearest_own_squared_distance)
        end
    end
end

const IS_DEVIANT_ATLAS_CONTRACT =
    renamed_contract(function_contract(compute_projected_modules_z_score_cutoff), "atlas_daf")

"""
    compute_vector_of_is_deviant_per_cell!(;
        query_daf::DafWriter,
        atlas_daf::DafReader,
        rejected_fraction::AbstractFloat = $(DEFAULT.rejected_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_is_deviant_per_cell`](@ref). A query cell is deviant when its projected modules z score is
above the cutoff of the atlas ([`compute_projected_modules_z_score_cutoff`](@ref)). That is, it fits its projected
metacell worse than all but `rejected_fraction` of the atlas cells fit their own metacell. Excluded cells are never
deviant. If the query has a metacell for each cell, then the cells which belong to a metacell are never deviant either.

# Query

$(CONTRACT1)

# Atlas

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "query_daf",
    axes = [cell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_cell(RequiredInput),
        vector_of_metacell_per_cell(OptionalInput),
        vector_of_projected_modules_z_score_per_cell(RequiredInput),
        vector_of_is_deviant_per_cell(CreatedOutput),
    ],
) IS_DEVIANT_ATLAS_CONTRACT function compute_vector_of_is_deviant_per_cell!(;  # UNTESTED
    query_daf::DafWriter,
    atlas_daf::DafReader,
    rejected_fraction::AbstractFloat = function_default(compute_projected_modules_z_score_cutoff, :rejected_fraction),
    overwrite::Bool = false,
)::Nothing
    z_score_cutoff = compute_projected_modules_z_score_cutoff(atlas_daf; rejected_fraction)

    is_excluded_per_query_cell = plain_array(get_vector(query_daf, "cell", "is_excluded"))
    z_score_per_query_cell = plain_array(get_vector(query_daf, "cell", "projected_metacell_modules_z_score"))
    is_deviant_per_query_cell = .!is_excluded_per_query_cell .& (z_score_per_query_cell .> z_score_cutoff)

    metacell_per_query_cell = get_vector(query_daf, "cell", "metacell"; default = nothing)
    if metacell_per_query_cell !== nothing
        is_deviant_per_query_cell .&= plain_array(metacell_per_query_cell) .== ""
    end

    n_included_query_cells = sum(.!is_excluded_per_query_cell)
    n_deviant_query_cells = sum(is_deviant_per_query_cell)
    @debug (
        "Z score cutoff: $(z_score_cutoff) Deviant cells: $(n_deviant_query_cells)" *
        " $(percent(n_deviant_query_cells, n_included_query_cells)) of the included cells"
    ) _group = :mcs_results

    set_vector!(query_daf, "cell", "is_deviant", is_deviant_per_query_cell; overwrite)
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

    # Only the atlas marker genes are correlated (the consumed mean uses a subset); zero for non-marker genes.
    query_included_genes = query_daf["@ gene [ ! is_excluded ]"]
    atlas_included_genes = atlas_daf["@ gene [ is_marker ]"]
    common_included_genes = intersect(query_included_genes, atlas_included_genes)
    n_common_included_genes = length(common_included_genes)

    @debug "Common included genes: $(n_common_included_genes)" *
           " ($(percent(n_common_included_genes, n_genes_in_atlas)) out of $(n_genes_in_atlas) atlas genes" *
           ", $(percent(n_common_included_genes, n_genes_in_query)) out of $(n_genes_in_query) query genes)" _group =
        :mcs_results

    indices_of_common_included_atlas_genes = axis_indices(atlas_daf, "gene", common_included_genes)
    indices_of_common_included_query_genes = axis_indices(query_daf, "gene", common_included_genes)

    indices_of_included_query_cells = query_daf["@ cell [ ! is_excluded ] : index"].array

    sparse_UMIs_per_query_cell_per_gene = plain_array(get_matrix(query_daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_query_cell = plain_array(get_vector(query_daf, "cell", "total_UMIs"))

    projected_metacell_per_query_cell = plain_array(get_vector(query_daf, "cell", "projected_metacell"))
    atlas_metacell_index_per_query_cell =
        axis_indices(atlas_daf, "metacell", projected_metacell_per_query_cell; allow_empty = true, allow_missing = true)

    UMIs_per_atlas_metacell_per_gene = plain_array(densify(get_matrix(atlas_daf, "metacell", "gene", "UMIs")))
    total_UMIs_per_atlas_metacell = plain_array(get_vector(atlas_daf, "metacell", "total_UMIs"))

    # All included query cells form a single group, each mapped to its projected atlas metacell (cells without one do
    # not participate); the sparse kernel then correlates each gene over these cells against their projected metacells.
    n_query_cells = axis_length(query_daf, "cell")
    is_participating_per_query_cell = falses(n_query_cells)
    is_participating_per_query_cell[indices_of_included_query_cells] .= true
    is_participating_per_query_cell .&= atlas_metacell_index_per_query_cell .> 0
    participating_query_cell_indices = findall(is_participating_per_query_cell)
    is_member_per_group_per_query_cell = sparse(
        ones(Int32, length(participating_query_cell_indices)),
        participating_query_cell_indices,
        trues(length(participating_query_cell_indices)),
        1,
        n_query_cells,
    )

    correlation_between_cells_and_projected_metacells_per_query_gene_per_group = zeros(Float32, n_genes_in_query, 1)

    sparse_cell_reference_correlation_per_gene_per_group!(
        correlation_between_cells_and_projected_metacells_per_query_gene_per_group,
        is_member_per_group_per_query_cell,
        sparse_UMIs_per_query_cell_per_gene,
        total_UMIs_per_query_cell,
        atlas_metacell_index_per_query_cell,
        UMIs_per_atlas_metacell_per_gene,
        total_UMIs_per_atlas_metacell,
        indices_of_common_included_query_genes,
        indices_of_common_included_atlas_genes,
        gene_fraction_regularization;
        name = "correlation_between_cells_and_projected_metacells",
    )
    correlation_between_cells_and_projected_metacells_per_query_gene =
        vec(correlation_between_cells_and_projected_metacells_per_query_gene_per_group)

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

Compute and set
[`matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block`](@ref). If
there are less than `min_neighborhood_query_cells`, we set this to zero. This will create [`projected_block_axis`](@ref)
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
) function compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!(;  # UNTESTED
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

    n_atlas_blocks = axis_length(atlas_daf, "block")
    n_genes_in_atlas = axis_length(atlas_daf, "gene")
    n_genes_in_query = axis_length(query_daf, "gene")

    name_per_atlas_block = axis_vector(atlas_daf, "block")
    # Only the atlas marker genes are correlated (the consumed mean uses a subset); zero for non-marker genes.
    query_included_genes = query_daf["@ gene [ ! is_excluded ]"]
    atlas_included_genes = atlas_daf["@ gene [ is_marker ]"]
    common_included_genes = intersect(query_included_genes, atlas_included_genes)
    n_common_included_genes = length(common_included_genes)

    @debug "Common included genes: $(n_common_included_genes)" *
           " ($(percent(n_common_included_genes, n_genes_in_atlas)) out of $(n_genes_in_atlas) atlas genes" *
           ", $(percent(n_common_included_genes, n_genes_in_query)) out of $(n_genes_in_query) query genes)" _group =
        :mcs_results

    indices_of_common_included_atlas_genes = axis_indices(atlas_daf, "gene", common_included_genes)
    indices_of_common_included_query_genes = axis_indices(query_daf, "gene", common_included_genes)

    sparse_UMIs_per_query_cell_per_gene = plain_array(get_matrix(query_daf, "cell", "gene", "UMIs"))::SparseMatrixCSC
    total_UMIs_per_query_cell = plain_array(get_vector(query_daf, "cell", "total_UMIs"))

    projected_metacell_per_query_cell = plain_array(get_vector(query_daf, "cell", "projected_metacell"))
    atlas_metacell_index_per_query_cell =
        axis_indices(atlas_daf, "metacell", projected_metacell_per_query_cell; allow_empty = true, allow_missing = true)

    projected_block_per_query_cell = plain_array(get_vector(query_daf, "cell", "projected_block"))
    atlas_block_index_per_query_cell =
        axis_indices(atlas_daf, "block", projected_block_per_query_cell; allow_empty = true, allow_missing = true)

    UMIs_per_atlas_metacell_per_gene = plain_array(densify(get_matrix(atlas_daf, "metacell", "gene", "UMIs")))
    total_UMIs_per_atlas_metacell = plain_array(get_vector(atlas_daf, "metacell", "total_UMIs"))

    is_lateral_per_atlas_gene = plain_array(get_vector(atlas_daf, "gene", "is_lateral"))
    is_lateral_per_common_included_gene = is_lateral_per_atlas_gene[indices_of_common_included_atlas_genes]

    is_marker_per_atlas_gene = plain_array(get_vector(atlas_daf, "gene", "is_marker"))
    is_marker_per_common_included_gene = is_marker_per_atlas_gene[indices_of_common_included_atlas_genes]
    is_pertinent_marker_per_common_included_gene =
        is_marker_per_common_included_gene .& .! is_lateral_per_common_included_gene

    is_in_neighborhood_per_other_block_per_base_block =
        plain_array(get_matrix(atlas_daf, "block", "block", "is_in_neighborhood"))
    is_neighborhood_marker_per_atlas_gene_per_atlas_block =
        get_matrix(atlas_daf, "gene", "block", "is_neighborhood_marker"; default = nothing)
    if is_neighborhood_marker_per_atlas_gene_per_atlas_block !== nothing
        is_neighborhood_marker_per_atlas_gene_per_atlas_block =
            plain_array(is_neighborhood_marker_per_atlas_gene_per_atlas_block)
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
        is_in_bin_per_atlas_gene = plain_array(get_vector(daf, "gene", "bin")) .== bin
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

    # Correlate each query gene's expression across each atlas block's neighborhood query cells against those cells'
    # projected metacells, filling the (gene, block) matrix. Each query cell (in the atlas block it projected onto)
    # belongs to the neighborhoods (groups) of the atlas blocks whose neighborhood includes that block; a cell without a
    # projected metacell or block participates in no group.
    is_participating_per_query_cell =
        (atlas_block_index_per_query_cell .> 0) .& (atlas_metacell_index_per_query_cell .> 0)
    is_member_per_atlas_block_per_query_cell = cell_group_membership_from_block_neighborhoods(
        atlas_block_index_per_query_cell,
        is_in_neighborhood_per_other_block_per_base_block,
        is_participating_per_query_cell,
    )
    n_neighborhood_query_cells_per_atlas_block = vec(sum(is_member_per_atlas_block_per_query_cell; dims = 2))

    sparse_cell_reference_correlation_per_gene_per_group!(
        correlation_between_neighborhood_query_cells_and_projected_metacells_per_query_gene_per_atlas_block,
        is_member_per_atlas_block_per_query_cell,
        sparse_UMIs_per_query_cell_per_gene,
        total_UMIs_per_query_cell,
        atlas_metacell_index_per_query_cell,
        UMIs_per_atlas_metacell_per_gene,
        total_UMIs_per_atlas_metacell,
        indices_of_common_included_query_genes,
        indices_of_common_included_atlas_genes,
        gene_fraction_regularization;
        min_group_cells = min_neighborhood_query_cells,
        name = "correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block",
    )

    if mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing &&
       in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
        is_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
    else
        is_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
            [BitVector(undef, n_common_included_genes) for _ in 1:maxthreadid()]
        if in_bin_mean_neighborhood_pertinent_markers_correlation_per_atlas_block === nothing
            is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
            is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread = nothing
        else
            is_in_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
                [BitVector(undef, n_common_included_genes) for _ in 1:maxthreadid()]
            is_out_bin_neighborhood_pertinent_marker_per_common_included_gene_per_thread =
                [BitVector(undef, n_common_included_genes) for _ in 1:maxthreadid()]
        end
    end

    for atlas_block_index in 1:n_atlas_blocks
        n_neighborhood_query_cells = n_neighborhood_query_cells_per_atlas_block[atlas_block_index]
        if n_neighborhood_query_cells < min_neighborhood_query_cells
            if n_neighborhood_query_cells > 0
                @warn (
                    "Ignoring too few query cells: $(n_neighborhood_query_cells) for the neighborhood of the atlas" *
                    " block: $(name_per_atlas_block[atlas_block_index])"
                )
            end

        else
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
                if !any(is_neighborhood_pertinent_marker_per_common_included_gene)
                    n_block_neighborhood_markers =
                        sum(@view(is_neighborhood_marker_per_atlas_gene_per_atlas_block[:, atlas_block_index]))
                    n_common_included_block_neighborhood_markers = sum(
                        @view(
                            is_neighborhood_marker_per_atlas_gene_per_atlas_block[
                                indices_of_common_included_atlas_genes,
                                atlas_block_index,
                            ]
                        ),
                    )
                    n_common_included_lateral = sum(is_lateral_per_common_included_gene)
                    @warn (
                        "Skipping atlas block: $(name_per_atlas_block[atlas_block_index]): no non-lateral" *
                        " neighborhood markers\n" *
                        "  block neighborhood markers: $(n_block_neighborhood_markers)\n" *
                        "  common+included block neighborhood markers:" *
                        " $(n_common_included_block_neighborhood_markers)\n" *
                        "  common+included lateral genes: $(n_common_included_lateral)\n" *
                        "  common+included genes: $(length(indices_of_common_included_atlas_genes))"
                    )
                    continue
                end

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
                        is_neighborhood_pertinent_marker_per_common_included_gene & !is_in_bin_per_common_included_gene
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
                "Mean correlation of $(first_qualifier)$(second_qualifier) pertinent marker genes between" *
                " neighborhood cells and their$(projected_qualifier) projected metacells: " *
                "$(mean(mean_per_atlas_block[mean_per_atlas_block .!= 0]))"  # NOLINT
            ) _group = :mcs_results
        end
    end

    return nothing
end

end  # module
