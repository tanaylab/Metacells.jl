"""
Do simple blocks analysis.
"""
module AnalyzeBlocks

export compute_blocks_covered_UMIs!
export compute_blocks_genes_covered_fractions!
export compute_blocks_genes_is_environment_markers!
export compute_blocks_genes_linear_fractions!
export compute_blocks_genes_log_covered_fractions!
export compute_blocks_genes_log_linear_fractions!
export compute_blocks_genes_UMIs!
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
export compute_blocks_total_UMIs!
export compute_blocks_types!

using Base.Threads
using DataAxesFormats
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
import Metacells.Contracts.block_block_mean_euclidean_skeleton_distance
import Metacells.Contracts.block_covered_UMIs_vector
import Metacells.Contracts.block_gene_covered_fraction_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_linear_fraction_matrix
import Metacells.Contracts.block_gene_log_covered_fraction_matrix
import Metacells.Contracts.block_gene_log_linear_fraction_matrix
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
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_log_covered_fraction_matrix
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

The maximal significant skeleton genes covered fractions fold factor between metacells of the blocks.
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

The mean Euclidean skeleton genes covered fractions distance between the metacells of the blocks.
"""
@logged @computation Contract(;
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        metacell_metacell_euclidean_skeleton_distance(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_mean_euclidean_skeleton_distance(GuaranteedOutput),
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
        base_metacells_indices = metacell_indices_per_block[base_block_index]

        @views distance_per_metacell_per_base_metacell = distances_between_metacells[:, base_metacells_indices]

        aggregated_distance_from_base_per_metacell = vec(aggregation(distance_per_metacell_per_base_metacell; dims = 2))
        @assert length(aggregated_distance_from_base_per_metacell) == n_metacells

        for other_block_index in 1:base_block_index
            other_metacells_indices = metacell_indices_per_block[other_block_index]

            aggregated_distance_from_base_per_other_metacell =
                aggregated_distance_from_base_per_metacell[other_metacells_indices]

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
    UMIs_per_block_per_gene = daf["/ metacell / gene : UMIs @ block ! %> Sum"].array
    set_matrix!(daf, "block", "gene", "UMIs", bestify(UMIs_per_block_per_gene); overwrite)
    return nothing
end

"""
    function compute_blocks_genes_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

An estimated linear fraction of the UMIs of each gene in each block. This is just the total UMIs of the gene in the
metacell divided by the total UMIs of the metacell, which is the "best" estimate assuming multinomial sampling noise.
However, this is sensitive to a few metacells with very high expression levels ("bursty" genes).

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
    function compute_blocks_genes_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

An estimated linear fraction of the UMIs of each covered gene in each block. By considering only the covered genes
this avoid the impact of highly-expressed lateral genes (e.g., cell cycle). Otherwise is similar to
[`compute_blocks_genes_linear_fractions!`](@ref).

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_gene_UMIs_matrix(RequiredInput),
        block_gene_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_linear_fractions(daf; qualifier = "covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_log_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each block. This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_gene_covered_fraction_matrix(RequiredInput),
        block_gene_log_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_log_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "covered", overwrite)
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
    n_metacells_per_block = daf["/ metacell @ block ! %> Count"].array
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
            daf["/ metacell & block => is_in_neighborhood ;= $(block_name) %> Count"]
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
            daf["/ metacell & block => is_in_environment ;= $(block_name) %> Count"]
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
    n_cells_per_block = daf["/ metacell : n_cells @ block ! %> Sum"].array
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

The total UMIs of genes per block.

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
    function compute_blocks_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of covered genes per block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_gene_UMIs_matrix(RequiredInput),
        block_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_blocks_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_UMIs(daf; qualifier = "covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    compute_blocks_genes_is_environment_markers!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
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
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_gene_log_covered_fraction_matrix(RequiredInput),
        block_gene_is_environment_marker_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_is_environment_markers!(
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = function_default(
        identify_marker_genes!,
        :min_marker_gene_max_fraction,
    ),
    min_marker_gene_range_fold::Real = function_default(identify_marker_genes!, :min_marker_gene_range_fold),
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    is_environment_marker_per_gene_per_block = Matrix{Bool}(undef, n_genes, n_blocks)

    for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        adapter(  # NOJET
            daf;
            input_axes = ["metacell" => "/ metacell & block => is_in_environment ;= $(block_name)", "gene" => "="],
            input_data = [
                ("metacell", "gene", "linear_fraction") => "covered_fraction",
                ("metacell", "gene", "log_linear_fraction") => "log_covered_fraction",
            ],
            output_axes = [],
            output_data = [],
        ) do adapted
            identify_marker_genes!(adapted; min_marker_gene_max_fraction, min_marker_gene_range_fold)
            is_environment_marker_per_gene_per_block[:, block_index] = get_vector(adapted, "gene", "is_marker").array
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

"""
    function compute_blocks_types!(
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
) function compute_blocks_types!(  # UNTESTED
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

    mask_gene_indices = daf["/ gene &$(genes_mask) : index"].array

    n_blocks = axis_length(daf, "block")
    n_genes = axis_length(daf, "gene")

    linear_fraction_per_block_per_gene = zeros(Float32, n_blocks, n_genes)
    linear_fraction_per_block_per_gene[:, mask_gene_indices] .=
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

end  # module
