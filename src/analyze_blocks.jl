"""
Do simple blocks analysis.
"""
module AnalyzeBlocks

export compute_blocks_covered_UMIs!
export compute_blocks_genes_is_environment_markers!
export compute_blocks_genes_UMIs!
export compute_blocks_genes_linear_covered_fractions!
export compute_blocks_genes_linear_fractions!
export compute_blocks_genes_log_linear_covered_fractions!
export compute_blocks_genes_log_linear_fractions!
export compute_blocks_genes_log_scaled_linear_covered_fractions!
export compute_blocks_genes_log_scaled_linear_fractions!
export compute_blocks_genes_scaled_linear_covered_fractions!
export compute_blocks_genes_scaled_linear_fractions!
export compute_blocks_n_cells!
export compute_blocks_n_environment_blocks!
export compute_blocks_n_environment_cells!
export compute_blocks_n_environment_metacells!
export compute_blocks_n_metacells!
export compute_blocks_n_neighborhood_blocks!
export compute_blocks_n_neighborhood_cells!
export compute_blocks_n_neighborhood_metacells!
export compute_blocks_scaled_covered_UMIs!
export compute_blocks_scaled_total_UMIs!
export compute_blocks_total_UMIs!
export compute_blocks_types!

using DataAxesFormats
using TanayLabUtilities

using ..AnalyzeGenes
using ..Defaults
using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_covered_UMIs_vector
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_linear_covered_fraction_matrix
import Metacells.Contracts.block_gene_linear_fraction_matrix
import Metacells.Contracts.block_gene_log_linear_covered_fraction_matrix
import Metacells.Contracts.block_gene_log_linear_fraction_matrix
import Metacells.Contracts.block_gene_log_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.block_gene_log_scaled_linear_fraction_matrix
import Metacells.Contracts.block_gene_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.block_gene_scaled_linear_fraction_matrix
import Metacells.Contracts.block_gene_UMIs_matrix
import Metacells.Contracts.block_n_cells_vector
import Metacells.Contracts.block_n_environment_blocks_vector
import Metacells.Contracts.block_n_environment_cells_vector
import Metacells.Contracts.block_n_environment_metacells_vector
import Metacells.Contracts.block_n_metacells_vector
import Metacells.Contracts.block_n_neighborhood_blocks_vector
import Metacells.Contracts.block_n_neighborhood_cells_vector
import Metacells.Contracts.block_n_neighborhood_metacells_vector
import Metacells.Contracts.block_scaled_covered_UMIs_vector
import Metacells.Contracts.block_scaled_total_UMIs_vector
import Metacells.Contracts.block_total_UMIs_vector
import Metacells.Contracts.block_type_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_fraction_matrix
import Metacells.Contracts.metacell_gene_log_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.metacell_type_vector

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
    function compute_blocks_genes_scaled_linear_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each gene in each block, scaled by divergence. We apply this scaling to
reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        block_gene_linear_fraction_matrix(RequiredInput),
        block_gene_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_scaled_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_scaled_fractions(daf; qualifier = "linear", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_log_scaled_linear_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each gene in each block, scaled by divergence. We apply this scaling to
reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_gene_scaled_linear_fraction_matrix(RequiredInput),
        block_gene_log_scaled_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_log_scaled_linear_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_blocks_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    do_compute_blocks_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "scaled_linear", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_linear_covered_fractions!(
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
        block_gene_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_linear_fractions(daf; qualifier = "linear_covered", genes_mask = "is_covered", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_log_linear_covered_fractions!(
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
        block_gene_linear_covered_fraction_matrix(RequiredInput),
        block_gene_log_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_log_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_log_fractions(daf; gene_fraction_regularization, qualifier = "linear_covered", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_scaled_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each covered gene in each block, scaled by divergence. We apply this
scaling to reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        block_gene_linear_covered_fraction_matrix(RequiredInput),
        block_gene_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_scaled_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_blocks_genes_scaled_fractions(daf; qualifier = "linear_covered", overwrite)
    return nothing
end

"""
    function compute_blocks_genes_log_scaled_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The estimated linear fraction of the UMIs of each covered gene in each block, scaled by divergence. We apply this
scaling to reduce the disproportionate impact of highly variable ("bursty") genes when using square-error methods.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_gene_scaled_linear_covered_fraction_matrix(RequiredInput),
        block_gene_log_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_genes_log_scaled_linear_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_blocks_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    do_compute_blocks_genes_log_fractions(
        daf;
        gene_fraction_regularization,
        qualifier = "scaled_linear_covered",
        overwrite,
    )
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
    function compute_blocks_scaled_total_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of genes per block, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        block_gene_UMIs_matrix(RequiredInput),
        block_scaled_total_UMIs_vector(GuaranteedOutput),
    ],
) function compute_blocks_scaled_total_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_gene = get_vector(daf, "gene", "divergence").array
    do_compute_blocks_UMIs(
        daf;
        qualifier = "scaled_total",
        genes_mask = "",
        scale_per_gene = 1.0 .- divergence_per_gene,
        overwrite,
    )
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
    function compute_blocks_scaled_covered_UMIs!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The total UMIs of covered genes per block, scaled by divergence.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [gene_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_divergence_vector(RequiredInput),
        gene_is_covered_vector(RequiredInput),
        block_gene_UMIs_matrix(RequiredInput),
        block_scaled_covered_UMIs_vector(GuaranteedOutput),
    ],
) function compute_blocks_scaled_covered_UMIs!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    divergence_per_scaled_gene = daf["/ gene & is_covered : divergence"].array
    do_compute_blocks_UMIs(
        daf;
        qualifier = "scaled_covered",
        genes_mask = "is_covered",
        scale_per_gene = 1.0 .- divergence_per_scaled_gene,
        overwrite,
    )
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

!!! note

    This uses the virtual [`metacell_gene_fraction_matrix`](@ref) and [`metacell_gene_log_fraction_matrix`](@ref). You
    will need an `adapter` to map these to concrete fractions (geomean, linear, scaled, ...). Make sure you are
    consistent when mapping the fractions and log-fraction matrices.
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        block_block_is_in_environment_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_fraction_matrix(RequiredInput),
        metacell_gene_log_fraction_matrix(RequiredInput),
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
            input_data = [("metacell", "gene", "fraction") => "=", ("metacell", "gene", "log_fraction") => "="],
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

function do_compute_blocks_genes_scaled_fractions(daf::DafWriter; qualifier::AbstractString, overwrite::Bool)::Nothing
    name = "$(qualifier)_fraction"
    fraction_per_gene_per_block = get_matrix(daf, "gene", "block", name).array

    divergence_per_gene = get_vector(daf, "gene", "divergence").array
    scale_per_gene = 1.0 .- divergence_per_gene

    scaled_fraction_per_gene_per_block = fraction_per_gene_per_block .* scale_per_gene

    set_matrix!(
        daf,
        "gene",
        "block",
        "scaled_$(name)",
        bestify(scaled_fraction_per_gene_per_block; eltype = Float32);
        overwrite,
    )

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
