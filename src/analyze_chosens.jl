"""
Do simple gene chosen analysis.
"""
module AnalyzeChosens

export compute_chosens_n_genes!
export compute_metacells_chosens_total_UMIs!
export compute_metacells_chosens_linear_fractions!
export compute_metacells_chosens_log_linear_fractions!
export compute_metacells_chosens_variance_over_means!
export compute_blocks_chosens_mean_variance_over_means!

using Base.Threads
using DataAxesFormats
using Random
using TanayLabUtilities
using StatsBase

using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.chosen_axis
import Metacells.Contracts.chosen_gene_is_member_matrix
import Metacells.Contracts.chosen_n_genes_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.metacell_chosen_total_UMIs_matrix
import Metacells.Contracts.metacell_total_UMIs_vector

"""
    function compute_chosens_n_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The number of genes in each chosen in each block.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [chosen_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [chosen_gene_is_member_matrix(RequiredInput), chosen_n_genes_vector(GuaranteedOutput)],
) function compute_chosens_n_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_chosens = axis_length(daf, "chosen")
    n_genes = axis_length(daf, "gene")
    is_member_per_gene_per_chosen = get_matrix(daf, "gene", "chosen", "is_member").array
    n_genes_per_chosen = vec(sum(is_member_per_gene_per_chosen; dims = 1))
    @assert_vector(n_genes_per_chosen, n_chosens)
    set_vector!(daf, "chosen", "n_genes", n_genes_per_chosen; eltype = UInt32, overwrite)
    return nothing
end

"""
    compute_metacells_chosens_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing

The total UMIs of each chosen each metacell.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), chosen_axis(RequiredInput)],
    data = [
        metacell_gene_UMIs_matrix(RequiredInput),
        chosen_gene_is_member_matrix(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        metacell_chosen_total_UMIs_matrix(GuaranteedOutput),
    ],
) function compute_metacells_chosens_total_UMIs!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_chosen = axis_length(daf, "chosen")

    is_member_per_gene_per_chosen = get_matrix(daf, "gene", "chosen", "is_member").array
    UMIs_per_gene_per_metacell = get_matrix(daf, "gene", "metacell", "UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array
    total_UMIs_per_metacell_per_chosen = Matrix{UInt32}(undef, n_metacells, n_chosen)

    parallel_loop_wo_rng(1:n_chosen; progress = DebugProgress(n_chosen)) do chosen_index
        is_member_per_gene = is_member_per_gene_per_chosen[:, chosen_index]
        chosen_UMIs_per_gene_per_metacell = UMIs_per_gene_per_metacell[is_member_per_gene, :]
        total_chosen_UMIs_per_metacell = vec(sum(chosen_UMIs_per_gene_per_metacell; dims = 1))
        @assert_vector(total_chosen_UMIs_per_metacell, n_metacells)
        total_UMIs_per_metacell_per_chosen[:, chosen_index] .= total_chosen_UMIs_per_metacell
        return nothing
    end

    set_matrix!(daf, "metacell", "chosen", "total_UMIs", bestify(total_UMIs_per_metacell_per_chosen); overwrite)
    return nothing
end

"""
    compute_metacells_chosens_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), chosen_axis(RequiredInput)],
    data = [
        metacell_total_UMIs_vector(RequiredInput),
        metacell_chosen_total_UMIs_matrix(RequiredInput),
        metacell_chosen_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_chosens_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing
    total_UMIs_per_metacell_per_chosen = get_matrix(daf, "metacell", "chosen", "total_UMIs").array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array
    linear_fraction_per_metacell_per_chosen = total_UMIs_per_metacell_per_chosen ./ total_UMIs_per_metacell
    set_matrix!(daf, "metacell", "chosen", "linear_fraction", linear_fraction_per_metacell_per_chosen; overwrite)
    return nothing
end

"""
    compute_metacells_chosens_log_linear_fractions!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), chosen_axis(RequiredInput)],
    data = [
        metacell_chosen_linear_fraction_matrix(RequiredInput),
        metacell_chosen_log_linear_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_chosens_log_linear_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    linear_fraction_per_metacell_per_chosen = get_matrix(daf, "metacell", "chosen", "linear_fraction").array
    log_linear_fraction_per_metacell_per_chosen =
        log2.(linear_fraction_per_metacell_per_chosen .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "metacell",
        "chosen",
        "log_linear_fraction",
        log_linear_fraction_per_metacell_per_chosen;
        overwrite,
    )
    return nothing
end

"""
    compute_metacells_chosens_v!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [
        cell_axis(RequiredInput),
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        chosen_axis(RequiredInput),
    ],
    data = [
        gene_is_excluded_vector(RequiredInput),
        chosen_gene_is_member_matrix(RequiredInput),
        cell_metacell_vector(RequiredInput),
        cell_gene_UMIs_matrix(RequiredInput),
        cell_total_UMIs_vector(RequiredInput),
        metacell_chosen_variance_over_mean_matrix(GuaranteedOutput),
    ],
) function compute_metacells_chosens_variance_over_means!(
    daf::DafWriter;
    min_downsamples::Integer = function_default(downsamples, :min_downsamples),
    min_downsamples_quantile::AbstractFloat = function_default(downsamples, :min_downsamples_quantile),
    max_downsamples_quantile::AbstractFloat = function_default(downsamples, :max_downsamples_quantile),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_chosen = axis_length(daf, "chosen")

    name_per_metacell = axis_vector(daf, "metacell")

    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    is_member_per_included_gene_per_chosen = daf["/ gene &! is_excluded / chosen : is_member"].array

    variance_over_mean_per_chosen_per_metacell = Matrix{Float32}(undef, n_chosen, n_metacells)

    parallel_loop_wo_rng(1:n_metacells; progress = DebugProgress(n_metacells)) do metacell_index
        metacell_name = name_per_metacell[metacell_index]
        indices_of_metacell_cells = daf["/ cell & metacell = $(metacell_name) : index"].array
        n_metacell_cells = length(indices_of_metacell_cells)

        UMIs_per_included_gene_per_metacell_cell =
            daf["/ gene &! is_excluded / cell & metacell = $(metacell_name) : UMIs"].array

        total_UMIs_per_metacell_cell = total_UMIs_per_cell[indices_of_metacell_cells]

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

        for chosen_index in 1:n_chosen
            @views is_of_chosen_per_included_gene = is_member_per_included_gene_per_chosen[:, chosen_index]

            @views downsampled_UMIs_per_chosen_gene_per_significant_metacell_cell =
                downsampled_UMIs_per_included_per_significant_metacell_cell[is_of_chosen_per_included_gene, :]

            chosen_total_downsampled_UMIs_per_significant_metacell_cell =
                vec(sum(downsampled_UMIs_per_chosen_gene_per_significant_metacell_cell; dims = 1))
            @assert_vector(chosen_total_downsampled_UMIs_per_significant_metacell_cell, n_significant_metacell_cells)

            mean_chosen_total_downsampled_UMIs = mean(chosen_total_downsampled_UMIs_per_significant_metacell_cell)
            if mean_chosen_total_downsampled_UMIs == 0
                variance_over_mean = 0
            else
                variance_of_chosen_total_downsampled_UMIs =
                    var(chosen_total_downsampled_UMIs_per_significant_metacell_cell)  # NOLINT
                variance_over_mean = variance_of_chosen_total_downsampled_UMIs ./ mean_chosen_total_downsampled_UMIs
                @assert variance_over_mean >= 0
            end

            variance_over_mean_per_chosen_per_metacell[chosen_index, metacell_index] = variance_over_mean
        end

        return nothing
    end

    set_matrix!(daf, "chosen", "metacell", "variance_over_mean", variance_over_mean_per_chosen_per_metacell; overwrite)
    return nothing
end

"""
    compute_blocks_chosens_mean_variance_over_means!(daf::DafWriter; overwrite::Bool = false)::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), block_axis(RequiredInput), chosen_axis(RequiredInput)],
    data = [
        metacell_block_vector(RequiredInput),
        metacell_chosen_variance_over_mean_matrix(RequiredInput),
        block_chosen_mean_variance_over_mean_matrix(GuaranteedOutput),
    ],
) function compute_blocks_chosens_mean_variance_over_means!(daf::DafWriter; overwrite::Bool = false)::Nothing
    mean_variance_over_mean_per_block_per_chosen =
        daf["/ metacell / chosen : variance_over_mean @ block ! %> Mean"].array
    set_matrix!(
        daf,
        "block",
        "chosen",
        "mean_variance_over_mean",
        mean_variance_over_mean_per_block_per_chosen;
        overwrite,
    )
    return nothing
end

end  # module
