"""
ro simple blocks analysis.
"""
module AnalyzeBlocks

export compute_blocks_cells_correlations!
export compute_blocks_confusion_by_closest_by_pertinent_markers!
export compute_blocks_fitting!
export compute_blocks_genes_UMIs!
export compute_blocks_genes_correlation_between_neighborhood_cells_and_metacells_matrix!
export compute_blocks_genes_is_correlated_with_skeletons_in_neighborhood!
export compute_blocks_genes_is_environment_markers!
export compute_blocks_genes_is_neighborhood_distincts!
export compute_blocks_genes_is_neighborhood_markers!
export compute_blocks_genes_is_neighborhood_varied!
export compute_blocks_genes_linear_fractions!
export compute_blocks_genes_log_linear_fractions!
export compute_blocks_is_in_neighborhood_by_confusion!
export compute_blocks_max_skeleton_fold_distances!
export compute_blocks_mean_euclidean_distances!
export compute_blocks_n_cells!
export compute_blocks_n_environment_blocks!
export compute_blocks_n_environment_cells!
export compute_blocks_n_environment_metacells!
export compute_blocks_n_metacells!
export compute_blocks_n_neighborhood_blocks!
export compute_blocks_n_neighborhood_cells!
export compute_blocks_n_neighborhood_metacells!
export compute_blocks_refitting!
export compute_blocks_total_UMIs!
export compute_blocks_types_by_metacells!
export compute_genes_most_changed_correlation_in_neighborhood_metacells_ranks!
export compute_suspect_genes!

using Base.Threads
using DataAxesFormats
using Distances
using LinearAlgebra
using StatsBase
using TanayLabUtilities

using ..AnalyzeGenes
using ..Defaults
using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_block_max_skeleton_fold_distance
import Metacells.Contracts.block_block_mean_euclidean_skeleton_distance_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_linear_fraction_matrix
import Metacells.Contracts.block_gene_log_linear_fraction_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_gene_correlation_between_neighborhood_cells_and_metacells_matrix
import Metacells.Contracts.block_gene_UMIs_matrix
import Metacells.Contracts.block_n_cells_vector
import Metacells.Contracts.block_n_environment_blocks_vector
import Metacells.Contracts.block_n_environment_cells_vector
import Metacells.Contracts.block_n_environment_metacells_vector
import Metacells.Contracts.block_n_metacells_vector
import Metacells.Contracts.block_n_neighborhood_blocks_vector
import Metacells.Contracts.block_n_neighborhood_cells_vector
import Metacells.Contracts.block_n_neighborhood_metacells_vector
import Metacells.Contracts.block_total_UMIs_vector
import Metacells.Contracts.block_type_vector
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.cell_total_UMIs_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_metacell_euclidean_skeleton_distance
import Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.metacell_type_vector

"""
    compute_blocks_max_skeleton_fold_distances!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The maximal significant skeleton genes fractions fold factor between metacells of the blocks.
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_metacell_max_skeleton_fold_distance(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_max_skeleton_fold_distance(GuaranteedOutput),
    ],
) function compute_blocks_max_skeleton_fold_distances!(daf::DafWriter; overwrite::Bool = false)::Nothing
    distances_between_metacells = get_matrix(daf, "metacell", "metacell", "max_skeleton_fold_distance").array
    block_index_per_metacell = daf["/ metacell : block => index"].array

    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    distances_between_blocks = compute_distances_between_blocks(;
        distances_between_metacells,
        metacell_indices_per_block,
        aggregation = maximum,
    )

    set_matrix!(daf, "block", "block", "max_skeleton_fold_distance", distances_between_blocks; overwrite)
    return nothing
end

"""
    compute_blocks_mean_euclidean_distances!(
        daf::DafWriter;
        overwrite::Bool = false,
    )::Nothing

The mean Euclidean skeleton genes fractions distance between the metacells of the blocks.
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_metacell_euclidean_skeleton_distance(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_mean_euclidean_skeleton_distance_matrix(GuaranteedOutput),
    ],
) function compute_blocks_mean_euclidean_distances!(daf::DafWriter; overwrite::Bool = false)::Nothing
    distances_between_metacells = get_matrix(daf, "metacell", "metacell", "euclidean_skeleton_distance").array
    block_index_per_metacell = daf["/ metacell : block => index"].array

    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    distances_between_blocks =
        compute_distances_between_blocks(; distances_between_metacells, metacell_indices_per_block, aggregation = mean)

    set_matrix!(daf, "block", "block", "mean_euclidean_skeleton_distance", distances_between_blocks; overwrite)
    return nothing
end

function compute_distances_between_blocks(;  # untested
    distances_between_metacells::AbstractMatrix{<:AbstractFloat},
    metacell_indices_per_block::AbstractVector{<:AbstractVector{<:Integer}},
    aggregation::Function,
)::AbstractMatrix{<:AbstractFloat}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacell_indices_per_block)
    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    @threads :greedy for base_block_index in reverse(1:n_blocks)
        indices_of_base_metacells = metacell_indices_per_block[base_block_index]

        @views distance_per_metacell_per_base_metacell = distances_between_metacells[:, indices_of_base_metacells]

        aggregated_distance_from_base_per_metacell = vec(aggregation(distance_per_metacell_per_base_metacell; dims = 2))
        @assert length(aggregated_distance_from_base_per_metacell) == n_metacells

        for other_block_index in 1:base_block_index
            indices_of_other_metacells = metacell_indices_per_block[other_block_index]

            aggregated_distance_from_base_per_other_metacell =
                aggregated_distance_from_base_per_metacell[indices_of_other_metacells]

            aggregated_distance_between_base_and_other_block =
                aggregation(aggregated_distance_from_base_per_other_metacell)

            distances_between_blocks[base_block_index, other_block_index] =
                aggregated_distance_between_base_and_other_block
            distances_between_blocks[other_block_index, base_block_index] =
                aggregated_distance_between_base_and_other_block
        end

        if base_block_index > 1
            @views distances_between_base_and_other_blocks =
                distances_between_blocks[1:(base_block_index - 1), base_block_index]
        end
    end

    return distances_between_blocks
end

"""
    function compute_blocks_genes_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of UMIs used to estimate the fraction of each gene in each block. This can used to estimate the
robustness of the fraction.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_gene_UMIs_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_gene_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    UMIs_per_block_per_gene = daf["/ metacell / gene : UMIs @ block ! %> Sum || 0"].array
    set_matrix!(daf, "block", "gene", "UMIs", bestify(UMIs_per_block_per_gene); overwrite)
    return nothing
end

"""
    function compute_blocks_genes_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

An estimated linear fraction of the UMIs of each non-excluded gene in each block. This is just the total UMIs of the
gene in the metacell divided by the total UMIs of the metacell, which is the "best" estimate assuming multinomial
sampling noise. However, this is sensitive to a few metacells with very high expression levels ("bursty" genes).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_is_excluded_vector(RequiredInput),
        block_gene_UMIs_matrix(RequiredInput),
        block_gene_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_linear_fractions(daf; qualifier = "linear", genes_mask = "!is_excluded", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_log_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of each gene in each block. This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [block_gene_linear_fraction_matrix(RequiredInput), block_gene_log_linear_fraction_matrix(GuaranteedOutput)],
) function compute_blocks_genes_log_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "linear", overwrite)
    return nothing
end

"""
    function compute_blocks_n_neighborhood_blocks!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of blocks in the neighborhood centered at each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [block_block_is_in_neighborhood_matrix(RequiredInput), block_n_neighborhood_blocks_vector(GuaranteedOutput)],
) function compute_blocks_n_neighborhood_blocks!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_neighborhood_blocks_per_block = daf["/ block / block : is_in_neighborhood %> Sum"].array
    set_vector!(daf, "block", "n_neighborhood_blocks", n_neighborhood_blocks_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_environment_blocks!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of blocks in the environment centered at each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [block_block_is_in_environment_matrix(RequiredInput), block_n_environment_blocks_vector(GuaranteedOutput)],
) function compute_blocks_n_environment_blocks!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_environment_blocks_per_block = daf["/ block / block : is_in_environment %> Sum"].array
    set_vector!(daf, "block", "n_environment_blocks", n_environment_blocks_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of metacells per block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [metacell_block_vector(RequiredInput), block_n_metacells_vector(GuaranteedOutput)],
) function compute_blocks_n_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells_per_block = daf["/ metacell @ block ! %> Count || 0"].array
    set_vector!(daf, "block", "n_metacells", n_metacells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_neighborhood_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of metacells in the blocks of the neighborhood centered at each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_n_neighborhood_metacells_vector(GuaranteedOutput),
    ],
) function compute_blocks_n_neighborhood_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")
    n_neighborhood_metacells_per_block = Vector{UInt32}(undef, n_blocks)
    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        n_neighborhood_metacells_per_block[block_index] =
            daf["/ metacell & block => is_in_neighborhood ;= $(block_name) %> Count || 0"]
    end
    set_vector!(daf, "block", "n_neighborhood_metacells", n_neighborhood_metacells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_environment_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of metacells in the blocks of the environment centered at each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_n_environment_metacells_vector(GuaranteedOutput),
    ],
) function compute_blocks_n_environment_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")
    n_environment_metacells_per_block = Vector{UInt32}(undef, n_blocks)
    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        n_environment_metacells_per_block[block_index] =
            daf["/ metacell & block => is_in_environment ;= $(block_name) %> Count || 0"]
    end
    set_vector!(daf, "block", "n_environment_metacells", n_environment_metacells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of cells per block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_n_cells_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_n_cells_vector(GuaranteedOutput),
    ],
) function compute_blocks_n_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells_per_block = daf["/ metacell : n_cells @ block ! %> Sum || 0"].array
    set_vector!(daf, "block", "n_cells", n_cells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_neighborhood_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of cells in the metacells of the blocks of the neighborhood centered at a block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_n_cells_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_n_neighborhood_cells_vector(GuaranteedOutput),
    ],
) function compute_blocks_n_neighborhood_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")
    n_neighborhood_cells_per_block = Vector{UInt32}(undef, n_blocks)
    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        n_neighborhood_cells_per_block[block_index] =
            daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : n_cells %> Sum"]
    end
    set_vector!(daf, "block", "n_neighborhood_cells", n_neighborhood_cells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_n_environment_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total number of cells in the metacells of the blocks of the environment centered at each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_n_cells_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_n_environment_cells_vector(GuaranteedOutput),
    ],
) function compute_blocks_n_environment_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")
    n_environment_cells_per_block = Vector{UInt32}(undef, n_blocks)
    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        n_environment_cells_per_block[block_index] =
            daf["/ metacell & block => is_in_environment ;= $(block_name) : n_cells %> Sum"]
    end
    set_vector!(daf, "block", "n_environment_cells", n_environment_cells_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of non-excluded genes per block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [block_gene_UMIs_matrix(RequiredInput), block_total_UMIs_vector(GuaranteedOutput)],
) function compute_blocks_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_UMIs(daf; qualifier = "total", genes_mask = "", overwrite)
    return nothing
end

"""
    compute_blocks_genes_is_environment_markers!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        overwrite::Bool = false,
    )::Nothing

A mask of genes that distinguish between cell states in each environment. Otherwise similar to
[`identify_marker_genes!`](@ref).
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_environment_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_environment_marker_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_environment_markers!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 2 ^ -14.5,
    min_marker_gene_range_fold::Real = 1.5,
    min_marker_quantile::Real = 0.1,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_environment_marker_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = ["metacell" => "/ metacell & block => is_in_environment ;= $(block_name)", "gene" => "="],
            input_data = [
                ("metacell", "gene", "linear_fraction") => "=",
                ("metacell", "gene", "log_linear_fraction") => "=",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            identify_marker_genes!(
                adapted;
                min_marker_gene_max_fraction,
                min_marker_gene_range_fold,
                min_marker_quantile,
            )
            is_environment_marker_per_gene = get_vector(adapted, "gene", "is_marker").array
            is_environment_marker_per_gene_per_block[:, block_index] = is_environment_marker_per_gene
            @debug "Block: $(block_name) environment markers: $(sum(is_environment_marker_per_gene))"
            return nothing
        end
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_environment_marker",
        bestify(is_environment_marker_per_gene_per_block);
        overwrite,
    )
    return nothing
end

@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_neighborhood_marker_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_neighborhood_markers!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 2 ^ -13.5,
    min_marker_gene_range_fold::Real = 1.0,
    min_marker_quantile::Real = 0.1,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_neighborhood_marker_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = ["metacell" => "/ metacell & block => is_in_neighborhood ;= $(block_name)", "gene" => "="],
            input_data = [
                ("metacell", "gene", "linear_fraction") => "=",
                ("metacell", "gene", "log_linear_fraction") => "=",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            identify_marker_genes!(
                adapted;
                min_marker_gene_max_fraction,
                min_marker_gene_range_fold,
                min_marker_quantile,
            )
            is_neighborhood_marker_per_gene = get_vector(adapted, "gene", "is_marker").array
            is_neighborhood_marker_per_gene_per_block[:, block_index] = is_neighborhood_marker_per_gene
            @debug "Block: $(block_name) neighborhood markers: $(sum(is_neighborhood_marker_per_gene))"
            return nothing
        end
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_neighborhood_marker",
        bestify(is_neighborhood_marker_per_gene_per_block);
        overwrite,
    )
    return nothing
end

"""
TODOX
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_is_skeleton_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_correlated_with_skeleton_in_neighborhood_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_correlated_with_skeletons_in_neighborhood!(
    daf::DafWriter;
    min_gene_max_fraction::AbstractFloat = function_default(
        identify_genes_correlated_with_skeletons!,
        :min_gene_max_fraction,
    ),
    min_gene_range_fold::Real = function_default(identify_genes_correlated_with_skeletons!, :min_gene_range_fold),
    min_gene_correlation::AbstractFloat = function_default(
        identify_genes_correlated_with_skeletons!,
        :min_gene_correlation,
    ),
    min_gene_correlation_quantile::AbstractFloat = function_default(
        identify_genes_correlated_with_skeletons!,
        :min_gene_correlation_quantile,
    ),
    genes_correlation_window::Integer = function_default(
        identify_genes_correlated_with_skeletons!,
        :genes_correlation_window,
    ),
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_correlated_with_skeleton_in_neighborhood_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = ["metacell" => "/ metacell & block => is_in_neighborhood ;= $(block_name)", "gene" => "="],
            input_data = [
                ("gene", "is_skeleton") => "=",
                ("metacell", "gene", "linear_fraction") => "=",
                ("metacell", "gene", "log_linear_fraction") => "=",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            identify_genes_correlated_with_skeletons!(
                adapted;
                min_gene_max_fraction,
                min_gene_range_fold,
                min_gene_correlation,
                min_gene_correlation_quantile,
                genes_correlation_window,
            )
            is_correlated_with_skeleton_in_neighborhood_per_gene =
                get_vector(adapted, "gene", "is_correlated_with_skeleton").array
            is_correlated_with_skeleton_in_neighborhood_per_gene_per_block[:, block_index] =
                is_correlated_with_skeleton_in_neighborhood_per_gene
            @debug "Block: $(block_name) genes correlated with skeleton in neighborhood: $(sum(is_correlated_with_skeleton_in_neighborhood_per_gene))"
            return nothing
        end
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_correlated_with_skeleton_in_neighborhood",
        bestify(is_correlated_with_skeleton_in_neighborhood_per_gene_per_block);
        overwrite,
    )
    return nothing
end

"""
    old_compute_blocks_genes_is_neighborhood_markers!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        overwrite::Bool = false,
    )::Nothing

A mask of genes that distinguish between cell states in each neighborhood. Otherwise similar to
[`identify_marker_genes!`](@ref).
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_neighborhood_marker_matrix(GuaranteedOutput),
    ],
) function old_compute_blocks_genes_is_neighborhood_markers!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 2 ^ -14.5,
    min_marker_gene_range_fold::Real = 1.5,
    min_marker_quantile::Real = 0.1,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_neighborhood_marker_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = ["metacell" => "/ metacell & block => is_in_neighborhood ;= $(block_name)", "gene" => "="],
            input_data = [
                ("metacell", "gene", "linear_fraction") => "=",
                ("metacell", "gene", "log_linear_fraction") => "=",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            identify_marker_genes!(
                adapted;
                min_marker_gene_max_fraction,
                min_marker_gene_range_fold,
                min_marker_quantile,
            )
            is_neighborhood_marker_per_gene_per_block[:, block_index] = get_vector(adapted, "gene", "is_marker").array
            return nothing
        end
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_neighborhood_marker",
        bestify(is_neighborhood_marker_per_gene_per_block);
        overwrite,
    )
    return nothing
end

"""
    compute_blocks_genes_is_neighborhood_varied!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        overwrite::Bool = false,
    )::Nothing

A mask of genes that "very strongly" distinguish between cell states in each neighborhood. Otherwise similar to
[`identify_marker_genes!`](@ref).
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_neighborhood_varied_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_neighborhood_varied!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 2e-3,
    min_marker_gene_range_fold::Real = 2.0,
    min_marker_quantile::Real = function_default(identify_marker_genes!, :min_marker_quantile),
    overwrite::Bool = false,
)::Nothing
    adapter(
        daf;
        input_axes = ["block" => "=", "metacell" => "=", "gene" => "="],
        input_data = [
            ("block", "block", "is_in_neighborhood") => "=",
            ("metacell", "block") => "=",
            ("metacell", "gene", "linear_fraction") => "=",
            ("metacell", "gene", "log_linear_fraction") => "=",
        ],
        output_axes = ["block" => "=", "gene" => "="],
        output_data = [("block", "gene", "is_neighborhood_varied") => "is_neighborhood_marker"],
        overwrite,
    ) do adapted
        return old_compute_blocks_genes_is_neighborhood_markers!(
            adapted;
            min_marker_gene_max_fraction,
            min_marker_gene_range_fold,
            min_marker_quantile,
            overwrite = true,
        )
    end
    return nothing
end

"""
    function compute_blocks_types_by_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The type of each block, based on is metacell types.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        metacell_type_vector(RequiredInput),
        block_type_vector(GuaranteedOutput),
    ],
) function compute_blocks_types_by_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_block = daf["/ metacell : type @ block ! %> Mode"].array
    set_vector!(daf, "block", "type", type_per_block; overwrite)
    return nothing
end

function do_compute_blocks_genes_linear_fractions(
    daf::DafWriter;
    qualifier::AbstractString,
    genes_mask::AbstractString,
    overwrite::Bool,
)::Nothing
    UMIs_per_block_per_masked_gene = daf["/ block / gene &$(genes_mask) : UMIs"].array
    total_masked_UMIs_per_block = vec(sum(UMIs_per_block_per_masked_gene; dims = 2))

    indices_of_mask_genes = daf["/ gene &$(genes_mask) : index"].array

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    linear_fraction_per_block_per_gene = zeros(Float32, n_blocks, n_genes)
    linear_fraction_per_block_per_gene[:, indices_of_mask_genes] .=
        UMIs_per_block_per_masked_gene ./ total_masked_UMIs_per_block

    set_matrix!(daf, "block", "gene", "$(qualifier)_fraction", bestify(linear_fraction_per_block_per_gene); overwrite)  # NOJET

    return nothing
end

function do_compute_blocks_genes_log_fractions(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat,
    qualifier::AbstractString,
    overwrite::Bool,
)::Nothing
    name = "$(qualifier)_fraction"
    fraction_per_gene_per_block = daf["/ gene / block : $(name)"].array

    empty_dense_matrix!(daf, "gene", "block", "log_$(name)", Float32; overwrite) do log_fraction_per_gene_per_block
        log_fraction_per_gene_per_block .= log2.(fraction_per_gene_per_block .+ gene_fraction_regularization)
        return nothing
    end

    return nothing
end

function do_compute_blocks_UMIs(
    daf::DafWriter;
    qualifier::AbstractString,
    genes_mask::AbstractString,
    scale_per_gene::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
    overwrite::Bool,
)::Nothing
    if genes_mask == ""
        query = "/ gene / block : UMIs"
    else
        query = "/ gene &$(genes_mask) / block : UMIs"
    end
    UMIs_per_gene_per_block = daf[query].array
    if scale_per_gene !== nothing
        UMIs_per_gene_per_block = round.(UMIs_per_gene_per_block .* scale_per_gene)
    end
    UMIs_per_block = UInt32.(vec(sum(UMIs_per_gene_per_block; dims = 1)))
    set_vector!(daf, "block", "$(qualifier)_UMIs", UMIs_per_block; overwrite)
    return nothing
end

"""
    compute_blocks_genes_correlation_between_neighborhood_cells_and_metacells_matrix!(
        daf::DafWriter;
        gene_cell_fraction_regularization::AbstractFloat = $(DEFAULT.gene_cell_fraction_regularization),
        overwrite::Bool,
    )::Nothing

The correlation between cells and metacells gene expression levels in each block's neighborhood. TODOX
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        gene_is_excluded_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        metacell_gene_UMIs_matrix(RequiredInput),
        block_gene_correlation_between_neighborhood_cells_and_metacells_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_correlation_between_neighborhood_cells_and_metacells_matrix!(
    daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    indices_of_included_genes = daf["/ gene &! is_excluded : index"].array
    n_included_genes = length(indices_of_included_genes)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array
    metacell_index_per_cell = daf["/ cell : metacell ?? 0 => index"].array

    correlation_per_gene_per_block = zeros(Float32, n_genes, n_blocks)

    mean_correlation_per_block = Vector{Float32}(undef, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]

        indices_of_neighborhood_cells =
            daf["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : index"].array

        if length(indices_of_neighborhood_cells) == 0
            @debug "Block: $(block_name) Mean correlation: NA"
        else
            metacell_index_per_neighborhood_cell = metacell_index_per_cell[indices_of_neighborhood_cells]

            total_UMIs_per_neighborhood_cell = total_UMIs_per_cell[indices_of_neighborhood_cells]
            total_metacell_UMIs_per_neighborhood_cell = total_UMIs_per_metacell[metacell_index_per_neighborhood_cell]
            total_punctuated_metacell_UMIs_per_neighborhood_cell =
                total_metacell_UMIs_per_neighborhood_cell .- total_UMIs_per_neighborhood_cell

            @views UMIs_per_neighborhood_cell_per_gene = UMIs_per_cell_per_gene[indices_of_neighborhood_cells, :]
            @views metacell_UMIs_per_neighborhood_cell_per_gene =
                UMIs_per_metacell_per_gene[metacell_index_per_neighborhood_cell, :]

            @threads :greedy for included_gene_position in 1:n_included_genes
                gene_index = indices_of_included_genes[included_gene_position]
                UMIs_per_neighborhood_cell = UMIs_per_neighborhood_cell_per_gene[:, gene_index]
                metacell_UMIs_per_neighborhood_cell = metacell_UMIs_per_neighborhood_cell_per_gene[:, gene_index]

                cell_log_fraction_per_neighborhood_cell = log2.(
                    UMIs_per_neighborhood_cell ./ total_UMIs_per_neighborhood_cell .+ gene_cell_fraction_regularization,
                )
                punctuated_metacell_log_fraction_per_neighborhood_cell = log2.(
                    (
                        (metacell_UMIs_per_neighborhood_cell .- UMIs_per_neighborhood_cell) ./
                        total_punctuated_metacell_UMIs_per_neighborhood_cell
                    ) .+ gene_cell_fraction_regularization,
                )

                correlation =
                    cor(cell_log_fraction_per_neighborhood_cell, punctuated_metacell_log_fraction_per_neighborhood_cell)
                if isnan(correlation)
                    correlation = 0
                end

                correlation_per_gene_per_block[gene_index, block_index] = correlation
            end

            mean_correlation_per_block[block_index] =
                mean(correlation_per_gene_per_block[indices_of_included_genes, block_index])
            @debug "Block: $(block_name) Mean correlation: $(mean_correlation_per_block[block_index])"
        end
    end

    @debug "Mean correlation in blocks: $(mean(mean_correlation_per_block))"

    set_matrix!(
        daf,
        "gene",
        "block",
        "correlation_between_neighborhood_cells_and_metacells",
        bestify(correlation_per_gene_per_block);
        overwrite,
    )
    return nothing
end

"""
    compute_blocks_cells_correlations!(daf::DafWriter; overwrite::Bool = $(DEFAULT.overwrite))::Nothing

TODOX

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        cell_genes_correlation_with_punctuated_metacells_vector(RequiredInput),
        cell_markers_correlation_with_punctuated_metacells_vector(RequiredInput),
        cell_pertinent_markers_correlation_with_punctuated_metacells_vector(RequiredInput),
        cell_regulators_correlation_with_punctuated_metacells_vector(RequiredInput),
        cell_pertinent_regulators_correlation_with_punctuated_metacells_vector(RequiredInput),
        block_mean_cells_genes_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        block_mean_cells_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        block_mean_cells_pertinent_markers_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        block_mean_cells_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
        block_mean_cells_pertinent_regulators_correlation_with_punctuated_metacells_vector(GuaranteedOutput),
    ],
) function compute_blocks_cells_correlations!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_cells_genes_correlation_with_punctuated_metacells_per_block =
        daf["/ cell : genes_correlation_with_punctuated_metacells @ metacell ?? => block ! %> Mean"].array
    set_vector!(
        daf,
        "block",
        "mean_cells_genes_correlation_with_punctuated_metacells",
        mean_cells_genes_correlation_with_punctuated_metacells_per_block;
        overwrite,
    )
    @debug "Blocks mean_cells_genes_correlation_with_punctuated_metacells: $(mean(mean_cells_genes_correlation_with_punctuated_metacells_per_block))"

    mean_cells_markers_correlation_with_punctuated_metacells_per_block =
        daf["/ cell : markers_correlation_with_punctuated_metacells @ metacell ?? => block ! %> Mean"].array
    set_vector!(
        daf,
        "block",
        "mean_cells_markers_correlation_with_punctuated_metacells",
        mean_cells_markers_correlation_with_punctuated_metacells_per_block;
        overwrite,
    )
    @debug "Blocks mean_cells_markers_correlation_with_punctuated_metacells: $(mean(mean_cells_markers_correlation_with_punctuated_metacells_per_block))"

    mean_cells_pertinent_markers_correlation_with_punctuated_metacells_per_block =
        daf["/ cell : pertinent_markers_correlation_with_punctuated_metacells @ metacell ?? => block ! %> Mean"].array
    set_vector!(
        daf,
        "block",
        "mean_cells_pertinent_markers_correlation_with_punctuated_metacells",
        mean_cells_pertinent_markers_correlation_with_punctuated_metacells_per_block;
        overwrite,
    )
    @debug "Blocks mean_cells_pertinent_markers_correlation_with_punctuated_metacells: $(mean(mean_cells_pertinent_markers_correlation_with_punctuated_metacells_per_block))"

    mean_cells_regulators_correlation_with_punctuated_metacells_per_block =
        daf["/ cell : regulators_correlation_with_punctuated_metacells @ metacell ?? => block ! %> Mean"].array
    set_vector!(
        daf,
        "block",
        "mean_cells_regulators_correlation_with_punctuated_metacells",
        mean_cells_regulators_correlation_with_punctuated_metacells_per_block;
        overwrite,
    )
    @debug "Blocks mean_cells_regulators_correlation_with_punctuated_metacells: $(mean(mean_cells_regulators_correlation_with_punctuated_metacells_per_block))"

    mean_cells_pertinent_regulators_correlation_with_punctuated_metacells_per_block =
        daf["/ cell : pertinent_regulators_correlation_with_punctuated_metacells @ metacell ?? => block ! %> Mean"].array
    set_vector!(
        daf,
        "block",
        "mean_cells_pertinent_regulators_correlation_with_punctuated_metacells",
        mean_cells_pertinent_regulators_correlation_with_punctuated_metacells_per_block;
        overwrite,
    )
    @debug "Blocks mean_cells_pertinent_regulators_correlation_with_punctuated_metacells: $(mean(mean_cells_pertinent_regulators_correlation_with_punctuated_metacells_per_block))"

    return nothing
end

"""
    compute_blocks_fitting!(
        daf::DafWriter;
        gene_cell_fraction_regularization::AbstractFloat = $(DEFAULT.gene_cell_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_lateral_vector(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_is_neighborhood_marker_matrix(RequiredInput),
        block_gene_gene_is_most_correlated_in_neighborhood_cells_tensor(GuaranteedOutput),
    ],
) function compute_blocks_fitting!(
    daf::DafWriter;
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= gene_cell_fraction_regularization <= 1

    log_linear_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "log_linear_fraction").array

    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array

    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    is_neighborhood_marker_per_gene_per_block = get_matrix(daf, "gene", "block", "is_neighborhood_marker").array

    is_most_correlated_per_gene_per_gene = Matrix{Bool}(undef, n_genes, n_genes)  # TODOX: Should be vector instead of matrix

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        fill!(is_most_correlated_per_gene_per_gene, false)

        @views is_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene_per_block[:, block_index]
        is_pertinent_neighborhood_marker_per_gene = is_neighborhood_marker_per_gene .& .!is_lateral_per_gene
        indices_of_pertinent_neighborhood_markers = findall(is_pertinent_neighborhood_marker_per_gene)
        n_pertinent_neighborhood_markers = length(indices_of_pertinent_neighborhood_markers)

        indices_of_neighborhood_cells =
            daf["/ cell & metacell ?? => block => is_in_neighborhood ;= $(block_name) : index"].array

        total_UMIs_per_neighborhood_cell = total_UMIs_per_cell[indices_of_neighborhood_cells]

        UMIs_per_neighborhood_cell_per_pertinent_neighborhood_marker =
            densify(UMIs_per_cell_per_gene[indices_of_neighborhood_cells, indices_of_pertinent_neighborhood_markers])
        log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker = log2.(
            UMIs_per_neighborhood_cell_per_pertinent_neighborhood_marker ./ total_UMIs_per_neighborhood_cell .+
            gene_cell_fraction_regularization,
        )

        correlation_between_pertinent_neighborhood_markers =
            Matrix{Float32}(undef, n_pertinent_neighborhood_markers, n_pertinent_neighborhood_markers)
        correlation_between_pertinent_neighborhood_markers[1, 1] = 0
        @threads :greedy for base_position in reverse(2:n_pertinent_neighborhood_markers)
            correlation_between_pertinent_neighborhood_markers[base_position, base_position] = 0
            @views base_column =
                log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker[:, base_position]
            for other_position in 1:(base_position - 1)
                @views other_column =
                    log_fraction_per_neighborhood_cell_per_pertinent_neighborhood_marker[:, other_position]
                correlation = cor(base_column, other_column)
                if isnan(correlation)
                    correlation = 0
                end
                correlation_between_pertinent_neighborhood_markers[base_position, other_position] = correlation
                correlation_between_pertinent_neighborhood_markers[other_position, base_position] = correlation
            end
        end

        correlation_between_pertinent_neighborhood_markers[isnan.(
            correlation_between_pertinent_neighborhood_markers,
        )] .= 0

        for pertinent_neighborhood_marker_position in 1:n_pertinent_neighborhood_markers
            pertinent_neighborhood_marker_index =
                indices_of_pertinent_neighborhood_markers[pertinent_neighborhood_marker_position]
            @views correlation_with_pertinent_neighborhood_marker =
                correlation_between_pertinent_neighborhood_markers[:, pertinent_neighborhood_marker_position]
            most_correlated_pertinent_neighborhood_marker_position =
                argmax(correlation_with_pertinent_neighborhood_marker)
            most_correlated_neighborhood_marker_index =
                indices_of_pertinent_neighborhood_markers[most_correlated_pertinent_neighborhood_marker_position]
            is_most_correlated_per_gene_per_gene[
                most_correlated_neighborhood_marker_index,
                pertinent_neighborhood_marker_index,
            ] = true
        end

        @debug "- Block: $(block_name) cells: $(length(indices_of_neighborhood_cells)) base genes: $(n_pertinent_neighborhood_markers)"

        set_matrix!(
            daf,
            "gene",
            "gene",
            "$(block_name)_is_most_correlated_in_neighborhood_cells",
            sparsify(is_most_correlated_per_gene_per_gene);
            overwrite,
        )
    end

    return nothing
end

"""
    compute_blocks_refitting!(
        base_daf::DafWriter,
        other_daf::DafWriter;
        gene_cell_fraction_regularization::AbstractFloat = $(DEFAULT.gene_cell_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

$(CONTRACT)
"""
@logged @computation Contract(
    is_relaxed = true,
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_gene_is_most_correlated_in_neighborhood_cells_tensor(RequiredInput),
    ],
) Contract(
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        cell_projected_metacell_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
    ],
) function compute_blocks_refitting!(
    base_daf::DafWriter,
    other_daf::DafReader;
    other_name::Maybe{AbstractString} = nothing,
    gene_cell_fraction_regularization::AbstractFloat = 1e-4,  # TODOX
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= gene_cell_fraction_regularization <= 1
    if other_name === nothing
        other_name = "." * other_daf.name
    end

    @assert axis_vector(base_daf, "gene") == axis_vector(other_daf, "gene")
    @assert axis_vector(base_daf, "cell") == axis_vector(other_daf, "cell")

    n_genes = axis_length(base_daf, "gene")
    name_per_gene = axis_vector(base_daf, "gene")

    linear_fraction_per_other_metacell_per_gene =
        get_matrix(other_daf, "metacell", "gene", "linear_fraction").array
    other_metacells_per_cell = get_vector(other_daf, "cell", "metacell").array
    projected_metacells_per_cell = get_vector(other_daf, "cell", "metacell.projected").array

    n_base_blocks = axis_length(base_daf, "block")
    name_per_base_block = axis_vector(base_daf, "block")

    UMIs_per_cell_per_gene = get_matrix(base_daf, "cell", "gene", "UMIs").array
    total_UMIs_per_cell = get_vector(base_daf, "cell", "total_UMIs").array

    most_correlation_per_gene_per_base_block = zeros(Float32, n_genes, n_base_blocks)
    correlation_per_gene_per_base_block = zeros(Float32, n_genes, n_base_blocks)
    projected_correlation_per_gene_per_base_block = zeros(Float32, n_genes, n_base_blocks)

    mean_most_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
    mean_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)
    mean_projected_correlation_per_base_block = Vector{Float32}(undef, n_base_blocks)

    @threads :greedy for base_block_index in 1:n_base_blocks
        base_block_name = name_per_base_block[base_block_index]

        is_most_correlated_per_gene_per_gene =
            get_matrix(base_daf, "gene", "gene", "$(base_block_name)_is_most_correlated_in_neighborhood_cells").array
        is_correlated_per_gene = vec(maximum(is_most_correlated_per_gene_per_gene; dims = 1))
        indices_of_correlated_genes = findall(is_correlated_per_gene)
        n_correlated_genes = length(indices_of_correlated_genes)

        indices_of_neighborhood_cells =
            base_daf["/ cell & metacell ?? => block => is_in_neighborhood ;= $(base_block_name) : index"].array
        other_metacell_per_neighborhood_cell = other_metacells_per_cell[indices_of_neighborhood_cells]
        projected_metacell_per_neighborhood_cell = projected_metacells_per_cell[indices_of_neighborhood_cells]
        is_grouped_per_neighborhood_cell = other_metacell_per_neighborhood_cell .!= ""

        indices_of_grouped_neighborhood_cells = indices_of_neighborhood_cells[is_grouped_per_neighborhood_cell]
        n_grouped_neighborhood_cells = length(indices_of_grouped_neighborhood_cells)
        other_metacell_per_grouped_neighborhood_cell =
            other_metacell_per_neighborhood_cell[is_grouped_per_neighborhood_cell]
        projected_metacell_per_grouped_neighborhood_cell =
            projected_metacell_per_neighborhood_cell[is_grouped_per_neighborhood_cell]
        other_metacell_index_per_grouped_neighborhood_cell =
            axis_indices(other_daf, "metacell", other_metacell_per_grouped_neighborhood_cell)
        projected_metacell_index_per_grouped_neighborhood_cell =
            axis_indices(other_daf, "metacell", projected_metacell_per_grouped_neighborhood_cell)

        total_UMIs_per_grouped_neighborhood_cell = total_UMIs_per_cell[indices_of_grouped_neighborhood_cells]

        UMIs_per_grouped_neighborhood_cell_per_correlated_gene =
            densify(UMIs_per_cell_per_gene[indices_of_grouped_neighborhood_cells, indices_of_correlated_genes])
        cell_log_fraction_per_grouped_neighborhood_cell_per_correlated_gene = log2.(
            UMIs_per_grouped_neighborhood_cell_per_correlated_gene ./ total_UMIs_per_grouped_neighborhood_cell .+
            gene_cell_fraction_regularization,
        )

        @views most_correlation_per_gene = most_correlation_per_gene_per_base_block[:, base_block_index]
        @views correlation_per_gene = correlation_per_gene_per_base_block[:, base_block_index]
        @views projected_correlation_per_gene = projected_correlation_per_gene_per_base_block[:, base_block_index]

        @threads :greedy for correlated_gene_position in 1:n_correlated_genes
            @views cell_log_fraction_per_grouped_neighborhood_cell =
                cell_log_fraction_per_grouped_neighborhood_cell_per_correlated_gene[:, correlated_gene_position]
            @assert_vector(cell_log_fraction_per_grouped_neighborhood_cell, n_grouped_neighborhood_cells)

            correlated_gene_index = indices_of_correlated_genes[correlated_gene_position]
            @views is_friend_per_gene = is_most_correlated_per_gene_per_gene[:, correlated_gene_index]
            @assert sum(is_friend_per_gene) == 1
            friend_gene_index = findfirst(is_friend_per_gene)

            @views metacell_fraction_per_other_metacell_of_gene =
                linear_fraction_per_other_metacell_per_gene[:, correlated_gene_index]
            @views metacell_fraction_per_other_metacell_of_friend_gene =
                linear_fraction_per_other_metacell_per_gene[:, friend_gene_index]

            @views other_metacell_log_fraction_per_grouped_neighborhood_cell_of_friend_gene =
                log2.(
                    metacell_fraction_per_other_metacell_of_friend_gene[other_metacell_index_per_grouped_neighborhood_cell] .+
                    gene_cell_fraction_regularization
                )
            @views other_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene =
                log2.(
                    metacell_fraction_per_other_metacell_of_gene[other_metacell_index_per_grouped_neighborhood_cell] .+
                    gene_cell_fraction_regularization
                )
            @views projected_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene =
                log2.(
                    metacell_fraction_per_other_metacell_of_gene[projected_metacell_index_per_grouped_neighborhood_cell] .+
                    gene_cell_fraction_regularization
                )

            if minimum(cell_log_fraction_per_grouped_neighborhood_cell) ==
               maximum(cell_log_fraction_per_grouped_neighborhood_cell) ||
               minimum(other_metacell_log_fraction_per_grouped_neighborhood_cell_of_friend_gene) ==
               maximum(other_metacell_log_fraction_per_grouped_neighborhood_cell_of_friend_gene)
                most_correlation = 0
            else
                most_correlation = cor(
                    cell_log_fraction_per_grouped_neighborhood_cell,
                    other_metacell_log_fraction_per_grouped_neighborhood_cell_of_friend_gene,
                )
                @assert !isnan(most_correlation) && most_correlation != 0.0
            end
            most_correlation_per_gene[correlated_gene_index] = most_correlation

            if minimum(cell_log_fraction_per_grouped_neighborhood_cell) ==
               maximum(cell_log_fraction_per_grouped_neighborhood_cell) ||
               minimum(other_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene) ==
               maximum(other_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene)
                correlation = 0
            else
                correlation = cor(
                    cell_log_fraction_per_grouped_neighborhood_cell,
                    other_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene,
                )
                @assert !isnan(correlation) && correlation != 0.0
            end
            correlation_per_gene[correlated_gene_index] = correlation

            if minimum(cell_log_fraction_per_grouped_neighborhood_cell) ==
               maximum(cell_log_fraction_per_grouped_neighborhood_cell) ||
               minimum(projected_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene) ==
               maximum(projected_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene)
                projected_correlation = 0
            else
                projected_correlation = cor(
                    cell_log_fraction_per_grouped_neighborhood_cell,
                    projected_metacell_log_fraction_per_grouped_neighborhood_cell_of_gene,
                )
                @assert !isnan(projected_correlation) && projected_correlation != 0.0
            end
            projected_correlation_per_gene[correlated_gene_index] = projected_correlation
        end

        mean_most_correlation = mean(most_correlation_per_gene[is_correlated_per_gene])
        @assert mean_most_correlation > 0

        mean_correlation = mean(correlation_per_gene[is_correlated_per_gene])
        @assert mean_correlation > 0

        mean_projected_correlation = mean(projected_correlation_per_gene[is_correlated_per_gene])
        @assert mean_projected_correlation > 0

        mean_most_correlation_per_base_block[base_block_index] = mean_most_correlation
        mean_correlation_per_base_block[base_block_index] = mean_correlation
        mean_projected_correlation_per_base_block[base_block_index] = mean_projected_correlation

        @debug "- Block: $(base_block_name) cells: $(percent(length(indices_of_grouped_neighborhood_cells), length(indices_of_neighborhood_cells))) out of $(length(indices_of_neighborhood_cells)) genes: $(n_correlated_genes) mean most: $(mean_most_correlation) self: $(mean_correlation) projected: $(mean_projected_correlation)"
    end

    set_matrix!(
        base_daf,
        "gene",
        "block",
        "most_correlation_in_neighborhood_metacells$(other_name)",
        bestify(most_correlation_per_gene_per_base_block);
        overwrite,
    )

    set_matrix!(
        base_daf,
        "gene",
        "block",
        "correlation_in_neighborhood_metacells$(other_name)",
        bestify(correlation_per_gene_per_base_block);
        overwrite,
    )

    set_matrix!(
        base_daf,
        "gene",
        "block",
        "projected_correlation_in_neighborhood_metacells$(other_name)",
        bestify(projected_correlation_per_gene_per_base_block);
        overwrite,
    )

    @debug "Mean correlation of pertinent marker genes with friends (between cells and their metacells in neighborhoods): $(mean(mean_most_correlation_per_base_block))"
    @debug "Mean correlation of pertinent marker genes with self (between cells and their metacells in neighborhoods): $(mean(mean_correlation_per_base_block))"
    @debug "Mean correlation of pertinent marker genes with self (between cells and projected metacells in neighborhoods): $(mean(mean_projected_correlation_per_base_block))"

    return nothing
end

@logged function compute_genes_most_changed_correlation_in_neighborhood_metacells_ranks!(
    daf::DafWriter;
    base_name::AbstractString,
    other_name::AbstractString,
    rank_name::AbstractString,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")

    base_most_correlation_per_gene_per_block =
        get_matrix(daf, "gene", "block", "most_correlation_in_neighborhood_metacells$(base_name)").array

    other_most_correlation_per_gene_per_block =
        get_matrix(daf, "gene", "block", "most_correlation_in_neighborhood_metacells$(other_name)").array

    is_relevant_per_gene = vec(maximum(base_most_correlation_per_gene_per_block; dims = 2) .!= 0)
    @assert_vector(is_relevant_per_gene, n_genes)

    @views relevant_base_most_correlation_per_gene_per_block =
        base_most_correlation_per_gene_per_block[is_relevant_per_gene, :]
    @views relevant_other_most_correlation_per_gene_per_block =
        other_most_correlation_per_gene_per_block[is_relevant_per_gene, :]

    diff_most_correlation_per_relevant_gene_per_block =
        relevant_other_most_correlation_per_gene_per_block .- relevant_base_most_correlation_per_gene_per_block

    improved_rank_per_relevant_gene = rank_markers(diff_most_correlation_per_relevant_gene_per_block)
    diff_most_correlation_per_relevant_gene_per_block .*= -1
    degraded_rank_per_relevant_gene = rank_markers(diff_most_correlation_per_relevant_gene_per_block)

    improved_rank_per_gene = fill(typemax(UInt32), n_genes)
    improved_rank_per_gene[is_relevant_per_gene] .= improved_rank_per_relevant_gene

    degraded_rank_per_gene = fill(typemax(UInt32), n_genes)
    degraded_rank_per_gene[is_relevant_per_gene] .= degraded_rank_per_relevant_gene

    set_vector!(
        daf,
        "gene",
        "improved_correlation_in_neighborhood_metacells_rank$(rank_name)",
        improved_rank_per_gene;
        overwrite,
    )
    set_vector!(
        daf,
        "gene",
        "degraded_correlation_in_neighborhood_metacells_rank$(rank_name)",
        degraded_rank_per_gene;
        overwrite,
    )

    return nothing
end

@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_gene_is_neighborhood_distinct_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_neighborhood_distincts!(
    daf::DafWriter;
    min_distinct_gene_max_fraction::AbstractFloat = 2 ^ -14.5,
    min_distinct_gene_mean_fold::Real = 2,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_neighborhood_distinct_per_gene_per_block = zeros(Bool, n_genes, n_blocks)
    linear_fraction_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", "linear_fraction").array
    log_linear_fraction_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", "log_linear_fraction").array
    median_log_linear_fraction_per_gene = daf["/ metacell / gene : log_linear_fraction %> Median"].array
    @assert_vector(median_log_linear_fraction_per_gene, n_genes)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        @views is_neighborhood_distinct_per_gene = is_neighborhood_distinct_per_gene_per_block[:, block_index]
        compute_block_genes_is_neighborhood_distincts!(
            daf;
            min_distinct_gene_max_fraction,
            min_distinct_gene_mean_fold,
            block_name,
            linear_fraction_per_gene_per_metacell,
            log_linear_fraction_per_gene_per_metacell,
            median_log_linear_fraction_per_gene,
            is_neighborhood_distinct_per_gene,
        )
    end
    set_matrix!(
        daf,
        "gene",
        "block",
        "is_neighborhood_distinct",
        bestify(is_neighborhood_distinct_per_gene_per_block);
        overwrite,
    )
    return nothing
end

function compute_block_genes_is_neighborhood_distincts!(
    daf::DafWriter;
    min_distinct_gene_max_fraction::AbstractFloat,
    min_distinct_gene_mean_fold::Real,
    linear_fraction_per_gene_per_metacell::AbstractMatrix{<:AbstractFloat},
    log_linear_fraction_per_gene_per_metacell::AbstractMatrix{<:AbstractFloat},
    median_log_linear_fraction_per_gene::AbstractVector{<:AbstractFloat},
    block_name::AbstractString,
    is_neighborhood_distinct_per_gene::Union{AbstractVector{Bool}, BitVector},
)::Nothing
    indices_of_neighborhood_metacells = daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array

    @views linear_fraction_per_gene_per_neighborhood_metacell =
        linear_fraction_per_gene_per_metacell[:, indices_of_neighborhood_metacells]
    @views log_linear_fraction_per_gene_per_neighborhood_metacell =
        log_linear_fraction_per_gene_per_metacell[:, indices_of_neighborhood_metacells]
    max_linear_fraction_per_gene = vec(maximum(linear_fraction_per_gene_per_neighborhood_metacell; dims = 2))

    fold_per_gene_per_neighborhood_metacell =
        log_linear_fraction_per_gene_per_neighborhood_metacell .- median_log_linear_fraction_per_gene
    mean_fold_per_gene = vec(mean(fold_per_gene_per_neighborhood_metacell; dims = 2))

    is_neighborhood_distinct_per_gene .=
        (mean_fold_per_gene .>= min_distinct_gene_mean_fold) .&
        (max_linear_fraction_per_gene .>= min_distinct_gene_max_fraction)

    @debug "Block: $(block_name) distinct: $(sum(is_neighborhood_distinct_per_gene))"

    return nothing
end

"""
TODOX

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [cell_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        cell_metacell_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        cell_closest_by_pertinent_markers_block_vector(RequiredInput),
        block_n_cells_vector(RequiredInput),
        block_block_confusion_by_closest_by_pertinent_markers_matrix(GuaranteedOutput),
    ],
) function compute_blocks_confusion_by_closest_by_pertinent_markers!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_blocks = axis_length(daf, "block")

    n_cells_per_block = get_vector(daf, "block", "n_cells").array
    closest_block_index_per_cell = daf["/ cell & metacell : block.closest_by_pertinent_markers => index"].array
    original_block_index_per_cell = daf["/ cell : metacell ?? => block => index"].array

    confusion_per_other_block_per_original_block = zeros(UInt32, n_blocks, n_blocks)
    for (closest_block_index, original_block_index) in zip(closest_block_index_per_cell, original_block_index_per_cell)
        confusion_per_other_block_per_original_block[closest_block_index, original_block_index] += 1
    end

    set_matrix!(daf, "block", "block", "confusion_by_closest_by_pertinent_markers", confusion_per_other_block_per_original_block)

    n_stable_cells = sum(confusion_per_other_block_per_original_block[diagind(confusion_per_other_block_per_original_block)])
    @debug "Stable cells: $(n_stable_cells) ($(percent(n_stable_cells, sum(n_cells_per_block))))"

    return nothing
end

"""
TODOX

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput)],
    data = [
        block_n_cells_vector(RequiredInput),
        block_block_confusion_by_closest_by_pertinent_markers_matrix(RequiredInput),
        block_total_UMIs_vector(RequiredInput),
        block_n_metacells_vector(RequiredInput),
        block_block_mean_euclidean_skeleton_distance_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
    ],
) function compute_blocks_is_in_neighborhood_by_confusion!(
    daf::DafWriter;
    min_neighbour_confusion_fractions::AbstractFloat = 0.01,
    min_blocks_in_neighborhood::Integer = 5,
    min_metacells_in_neighborhood::Integer = 20,
    min_total_UMIs_in_neighborhood::Integer = 2_000_000,
    max_blocks_in_neighborhood::Integer = 100,
    max_metacells_in_neighborhood::Integer = 400,
    max_total_UMIs_in_neighborhood::Integer = 4_000_000_0,
    overwrite::Bool = false,
)::Nothing
    @assert 0 < min_blocks_in_neighborhood < max_blocks_in_neighborhood
    @assert 0 < min_metacells_in_neighborhood < max_metacells_in_neighborhood
    @assert 0 < min_total_UMIs_in_neighborhood < max_total_UMIs_in_neighborhood

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    n_cells_per_block = get_vector(daf, "block", "n_cells").array

    mean_euclidean_skeleton_distance_per_block_per_block =
        get_matrix(daf, "block", "block", "mean_euclidean_skeleton_distance")
    n_metacells_per_block = get_vector(daf, "block", "n_metacells")
    total_UMIs_per_block = get_vector(daf, "block", "total_UMIs")
    confusion_per_other_block_per_original_block = get_matrix(daf, "block", "block", "confusion_by_closest_by_pertinent_markers").array

    confusion_fraction_per_block_per_block = Float32.(confusion_per_other_block_per_original_block ./ n_cells_per_block)
    confusion_fraction_per_block_per_block .+= transpose(confusion_fraction_per_block_per_block)

    is_in_neighborhood_per_other_block_per_base_block = zeros(Bool, n_blocks, n_blocks)

    for base_block_index in 1:n_blocks
        @views confusion_fraction_per_other_block = confusion_fraction_per_block_per_block[:, base_block_index]
        @views distance_per_other_block = mean_euclidean_skeleton_distance_per_block_per_block[:, base_block_index]
        priority_per_other_block = collect(zip(.-confusion_fraction_per_other_block, distance_per_other_block))
        priority_per_other_block[base_block_index] = (-Inf32, -Inf32)
        ordered_block_indices = sortperm(priority_per_other_block)

        neighborhood_n_confused_blocks = 0
        neighborhood_n_disjoint_blocks = 0
        neighborhood_n_overflow_blocks = 0

        neighborhood_n_blocks = 0
        neighborhood_n_metacells = 0
        neighborhood_total_UMIs = 0

        while true
            next_block_index = ordered_block_indices[neighborhood_n_blocks + 1]
            include_next_block = false

            if confusion_fraction_per_other_block[next_block_index] >= min_neighbour_confusion_fractions
                if neighborhood_n_blocks <= max_blocks_in_neighborhood &&
                        neighborhood_n_metacells <= max_metacells_in_neighborhood &&
                        neighborhood_total_UMIs <= max_total_UMIs_in_neighborhood
                    include_next_block = true
                    neighborhood_n_confused_blocks += 1
                end
            else
                if neighborhood_n_blocks < min_blocks_in_neighborhood ||
                        neighborhood_n_metacells < min_metacells_in_neighborhood ||
                        neighborhood_total_UMIs < min_total_UMIs_in_neighborhood
                    include_next_block = true
                    neighborhood_n_disjoint_blocks += 1
                end
            end

            if !include_next_block
                break
            end

            neighborhood_n_blocks += 1
            is_in_neighborhood_per_other_block_per_base_block[next_block_index, base_block_index] = true
            neighborhood_n_metacells += n_metacells_per_block[next_block_index]
            neighborhood_total_UMIs += total_UMIs_per_block[next_block_index]
        end

        while neighborhood_n_blocks + neighborhood_n_overflow_blocks < n_blocks
            next_block_index = ordered_block_indices[neighborhood_n_blocks + neighborhood_n_overflow_blocks + 1]
            if confusion_fraction_per_other_block[next_block_index] >= min_neighbour_confusion_fractions
                neighborhood_n_overflow_blocks += 1
            else
                break
            end
        end

        @debug (
            "Neighborhood: $(name_per_block[base_block_index])" *
            " Blocks: $(neighborhood_n_blocks)" *
            " Confused: $(neighborhood_n_confused_blocks)" *
            " Disjoint: $(neighborhood_n_disjoint_blocks)" *
            " Overflow: $(neighborhood_n_overflow_blocks)" *
            " Metacells: $(neighborhood_n_metacells)" *
            " Covered M-UMIs: $(neighborhood_total_UMIs / 1e6)"
        )
        @assert neighborhood_n_disjoint_blocks == 0 || neighborhood_n_overflow_blocks == 0
    end

    set_matrix!(daf, "block", "block", "is_in_neighborhood", is_in_neighborhood_per_other_block_per_base_block)
    return nothing
end

mutable struct IncrementalCorrelation
    sum_x::Float64
    sum_x2::Float64
    sum_y::Float64
    sum_y2::Float64
    sum_xy::Float64
    num::UInt64
end

function IncrementalCorrelation()::IncrementalCorrelation
    return IncrementalCorrelation(0, 0, 0, 0, 0, 0)
end

function add_sample!(incremental_correlation::IncrementalCorrelation, x::Float64, y::Float64)::Nothing
    @assert !isnan(x)
    @assert !isnan(y)
    incremental_correlation.sum_x += x
    incremental_correlation.sum_x2 += x ^ 2
    incremental_correlation.sum_y += y
    incremental_correlation.sum_y2 += y ^ 2
    incremental_correlation.sum_xy += x * y
    incremental_correlation.num += 1
    return nothing
end

function combine_into!(into_incremental_correlation::IncrementalCorrelation, from_incremental_correlation::IncrementalCorrelation)::IncrementalCorrelation
    @assert !isnan(into_incremental_correlation.sum_x)
    @assert !isnan(into_incremental_correlation.sum_x2)
    @assert !isnan(into_incremental_correlation.sum_y)
    @assert !isnan(into_incremental_correlation.sum_y2)
    @assert !isnan(into_incremental_correlation.sum_xy)
    @assert !isnan(into_incremental_correlation.num)
    @assert !isnan(from_incremental_correlation.sum_x)
    @assert !isnan(from_incremental_correlation.sum_x2)
    @assert !isnan(from_incremental_correlation.sum_y)
    @assert !isnan(from_incremental_correlation.sum_y2)
    @assert !isnan(from_incremental_correlation.sum_xy)
    @assert !isnan(from_incremental_correlation.num)
    into_incremental_correlation.sum_x += from_incremental_correlation.sum_x
    into_incremental_correlation.sum_x2 += from_incremental_correlation.sum_x2
    into_incremental_correlation.sum_y += from_incremental_correlation.sum_y
    into_incremental_correlation.sum_y2 += from_incremental_correlation.sum_y2
    into_incremental_correlation.sum_xy += from_incremental_correlation.sum_xy
    into_incremental_correlation.num += from_incremental_correlation.num
    return into_incremental_correlation
end

function correlation(incremental_correlation::IncrementalCorrelation)::Float64
    if incremental_correlation.num == 0
        return 0
    else
        result = sqrt(
            max(
                0.0,
                (
                    incremental_correlation.num * incremental_correlation.sum_xy -
                    incremental_correlation.sum_x * incremental_correlation.sum_y
                ) / (
                    (
                        incremental_correlation.num * incremental_correlation.sum_x2 -
                        incremental_correlation.sum_x ^ 2
                    ) * (
                        incremental_correlation.num * incremental_correlation.sum_y2 -
                        incremental_correlation.sum_y ^ 2
                    )
                )
            )
        )
        if isnan(result)
            return 0.0
        else
            return result
        end
    end
end

@logged function compute_suspect_genes!(
    daf::DafWriter;
    method::Symbol,
    overwrite::Bool = false
)::Nothing
    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    indices_of_pertinent_markers = daf["/ gene & is_marker &! is_lateral : index"].array
    n_pertinent_markers = length(indices_of_pertinent_markers)

    confusion_per_other_block_per_original_block = get_matrix(daf, "block", "block", "confusion_by_closest_by_pertinent_markers").array

    into_per_other_block = vec(sum(confusion_per_other_block_per_original_block; dims = 2))

    log_fraction_per_block_per_pertinent_marker = daf["/ block / gene & is_marker &! is_lateral : log_linear_fraction"].array
    log_fraction_per_pertinent_marker_per_block = daf["/ gene & is_marker &! is_lateral / block : log_linear_fraction"].array

    UMIs_per_cell_per_pertinent_marker = daf["/ cell / gene & is_marker &! is_lateral : UMIs"].array
    total_UMIs_per_cell = daf["/ cell : total_UMIs"].array

    if method == :a
        correlation_per_pertinent_marker = vec(cor(into_per_other_block, log_fraction_per_block_per_pertinent_marker))
        correlation_per_pertinent_marker[isnan.(correlation_per_pertinent_marker)] .= 0

        @debug "TODOX QUANTILES: $(quantile(correlation_per_pertinent_marker, (0:10)./10))"
        @debug "TODOX MEDIAN: $(median(correlation_per_pertinent_marker))"
        @debug "TODOX MEAN: $(mean(correlation_per_pertinent_marker))"
        @debug "TODOX STD: $(std(correlation_per_pertinent_marker))"
        stdev = std(correlation_per_pertinent_marker)
        low = -2 * stdev
        high = 2 * stdev
        @debug "TODOX NON-CONFUSING RANGE: $(low) .. $(high)"

        z_score_per_pertinent_marker = correlation_per_pertinent_marker ./ stdev
        is_suspect_per_pertinent_marker = (correlation_per_pertinent_marker .< low) .| (high .< correlation_per_pertinent_marker)

    elseif method == :b
        # TODOX
        log_fraction_per_pertinent_marker_per_block =
            log2.(1e-4 .+ daf["/ gene & is_marker &! is_lateral / block : linear_fraction"].array)

        n_cells_per_block = get_vector(daf, "block", "n_cells").array
        original_block_index_per_cell = daf["/ cell : metacell ?? 0 => block => index"].array
        closest_block_index_per_cell = daf["/ cell : block.closest_by_pertinent_markers => index"].array

        incremental_correlation_per_pertinent_marker_per_block = Matrix{IncrementalCorrelation}(undef, n_pertinent_markers, n_blocks)
        for block_index in 1:n_blocks
            for pertinent_marker_position in 1:n_pertinent_markers
                incremental_correlation_per_pertinent_marker_per_block[pertinent_marker_position, block_index] =
                    IncrementalCorrelation()
            end
        end

        progress_counter = Atomic{Int}(0)
        for into_block_index in 1:n_blocks
            log_fraction_of_into_per_pertinent_marker = log_fraction_per_pertinent_marker_per_block[:, into_block_index]
            @assert_vector(log_fraction_of_into_per_pertinent_marker, n_pertinent_markers)

            @threads :greedy for from_block_index in 1:n_blocks
                if from_block_index != into_block_index
                    indices_of_migrated_cells = findall(
                        (original_block_index_per_cell .== from_block_index) .&
                        (closest_block_index_per_cell .== into_block_index)
                    )
                    n_migrated_cells = length(indices_of_migrated_cells)
                    @assert n_migrated_cells == confusion_per_other_block_per_original_block[into_block_index, from_block_index]
                    if n_migrated_cells > 0
                        fraction_from_to_into =
                            confusion_per_other_block_per_original_block[into_block_index, from_block_index] ./
                            n_cells_per_block[from_block_index]
                        fraction_into_to_from =
                            confusion_per_other_block_per_original_block[from_block_index, into_block_index] ./
                            n_cells_per_block[into_block_index]
                        into_more_than_from_fraction = fraction_from_to_into - fraction_into_to_from

                        log_fraction_per_migrated_cell_per_pertinent_marker =
                            log2.(
                                1e-4 .+  # TODOx
                                UMIs_per_cell_per_pertinent_marker[indices_of_migrated_cells, :] ./
                                total_UMIs_per_cell[indices_of_migrated_cells]
                            )
                        log_fraction_per_pertinent_marker_per_migrated_cell =
                            flip(log_fraction_per_migrated_cell_per_pertinent_marker)
                        @assert_matrix(log_fraction_per_pertinent_marker_per_migrated_cell, n_pertinent_markers, n_migrated_cells)

                        into_distance_per_pertinent_marker_per_migrated_cell =
                            abs.(
                                log_fraction_per_pertinent_marker_per_migrated_cell .-
                                log_fraction_of_into_per_pertinent_marker
                            )
                        @assert_matrix(into_distance_per_pertinent_marker_per_migrated_cell, n_pertinent_markers, n_migrated_cells)

                        log_fraction_of_from_per_pertinent_marker = log_fraction_per_pertinent_marker_per_block[:, from_block_index]
                        @assert_vector(log_fraction_of_from_per_pertinent_marker, n_pertinent_markers)
                        from_distance_per_pertinent_marker_per_migrated_cell =
                            abs.(
                                log_fraction_per_pertinent_marker_per_migrated_cell .-
                                log_fraction_of_from_per_pertinent_marker
                            )
                        @assert_matrix(from_distance_per_pertinent_marker_per_migrated_cell, n_pertinent_markers, n_migrated_cells)

                        into_better_than_from_per_pertinent_marker_per_migrated_cell =
                            from_distance_per_pertinent_marker_per_migrated_cell .-
                            into_distance_per_pertinent_marker_per_migrated_cell

                        if method == :c
                            into_better_than_from_per_pertinent_marker_per_migrated_cell .*=
                                into_better_than_from_per_pertinent_marker_per_migrated_cell
                        end

                        for migrated_cell_position in 1:n_migrated_cells
                            for pertinent_marker_position in 1:n_pertinent_markers
                                add_sample!(
                                    incremental_correlation_per_pertinent_marker_per_block[
                                        pertinent_marker_position,
                                        from_block_index
                                    ],
                                    into_better_than_from_per_pertinent_marker_per_migrated_cell[
                                        pertinent_marker_position,
                                        migrated_cell_position,
                                    ],
                                    into_more_than_from_fraction,
                                )
                            end
                        end
                    end
                end
            end
            counter = atomic_add!(progress_counter, 1)
            print("\r$(progress_counter[]) ($(percent(counter + 1, n_blocks))) ...")
        end

        correlation_per_pertinent_marker = Vector{Float32}(undef, n_pertinent_markers)
        @threads :greedy for pertinent_marker_position in 1:n_pertinent_markers
            pertinent_marker_incremental_correlation = IncrementalCorrelation()
            for block_index in 1:n_blocks
                combine_into!(
                    pertinent_marker_incremental_correlation,
                    incremental_correlation_per_pertinent_marker_per_block[pertinent_marker_position, block_index]
                )
            end
            correlation_per_pertinent_marker[pertinent_marker_position] = correlation(pertinent_marker_incremental_correlation)
        end

        @debug "TODOX QUANTILES: $(quantile(correlation_per_pertinent_marker, (0:10)./10))"
        @debug "TODOX MEDIAN: $(median(correlation_per_pertinent_marker))"
        @debug "TODOX MEAN: $(mean(correlation_per_pertinent_marker))"
        @debug "TODOX STD: $(std(correlation_per_pertinent_marker))"
        high = mean(correlation_per_pertinent_marker) + 2 * std(correlation_per_pertinent_marker)
        @debug "TODOX NON-CONFUSING RANGE: 0 .. $(high)"

        z_score_per_pertinent_marker = (correlation_per_pertinent_marker .- mean(correlation_per_pertinent_marker)) ./ std(correlation_per_pertinent_marker)
        is_suspect_per_pertinent_marker = correlation_per_pertinent_marker .> high
    else
        @assert false
    end

    @debug "Confusing total genes: $(sum(is_suspect_per_pertinent_marker))"

    is_suspect_per_gene = zeros(Bool, n_genes)
    is_suspect_per_gene[indices_of_pertinent_markers] .= is_suspect_per_pertinent_marker
    set_vector!(daf, "gene", "is_suspect", bestify(is_suspect_per_gene); overwrite)

    z_score_per_gene = zeros(Float32, n_genes)
    z_score_per_gene[indices_of_pertinent_markers] .= z_score_per_pertinent_marker
    set_vector!(daf, "gene", "suspicion_z_score", bestify(z_score_per_gene); overwrite)

    return nothing
end

end  # module
