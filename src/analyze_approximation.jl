"""
Do simple blocks analysis.
"""
module AnalyzeApproximation

export compute_metacells_genes_approximated_linear_covered_fractions!
export compute_metacells_genes_approximated_scaled_linear_covered_fractions!
export compute_metacells_genes_approximated_log_linear_covered_fractions!
export compute_metacells_genes_approximated_log_scaled_linear_covered_fractions!
export compute_blocks_modules_block_RMSE!
export compute_blocks_modules_neighborhood_RMSE!

using Base.Threads
using DataAxesFormats
using TanayLabUtilities
using Statistics

using ..AnalyzeMetacells
using ..Contracts

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_gene_base_covered_fraction_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_base_covered_fraction_matrix
import Metacells.Contracts.block_module_gene_covered_coefficient_tensor
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_modules_block_RMSE_vector
import Metacells.Contracts.block_modules_neighborhood_RMSE_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_approximated_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_approximated_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_approximated_log_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_approximated_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_module_linear_covered_fraction_matrix
import Metacells.Contracts.metacell_module_scaled_linear_covered_fraction_matrix
import Metacells.Contracts.module_axis

"""
    function compute_metacells_genes_approximated_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model. This
adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_module_gene_covered_coefficient_tensor(RequiredInput),
        block_module_base_covered_fraction_matrix(RequiredInput),
        block_gene_base_covered_fraction_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_module_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_approximated_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_approximated_linear_covered_fractions!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_approximated_linear_covered_fractions!(
        daf;
        covered_fraction = "linear_covered_fraction",
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_genes_approximated_scaled_linear_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The linear fraction of each covered gene UMIs out of the total covered UMIs, scaled by divergence, as approximated by
the linear model. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
        gene_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_module_gene_covered_coefficient_tensor(RequiredInput),
        block_module_base_covered_fraction_matrix(RequiredInput),
        block_gene_base_covered_fraction_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_module_scaled_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_approximated_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_approximated_scaled_linear_covered_fractions!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    do_compute_metacells_genes_approximated_linear_covered_fractions!(
        daf;
        covered_fraction = "scaled_linear_covered_fraction",
        overwrite,
    )
    return nothing
end

function do_compute_metacells_genes_approximated_linear_covered_fractions!(
    daf::DafWriter;
    covered_fraction::AbstractString,
    overwrite::Bool,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    covered_genes_indices = daf["/ gene & is_covered : index"].array

    approximated_linear_covered_fraction_per_gene_per_metacell = zeros(Float32, n_genes, n_metacells)

    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]

        block_metacell_indices = daf["/ metacell & block = $(block_name) : index"].array

        coefficient_per_covered_gene_per_found_module =
            daf["/ gene & is_covered / module & is_found ; block = $(block_name) : $(block_name)_covered_coefficient"].array

        @views base_covered_fraction_per_found_module_of_block =
            daf["/ module & is_found ; block = $(block_name) / block = $(block_name) : base_covered_fraction"].array

        @views base_covered_fraction_per_covered_gene_of_block =
            daf["/ gene & is_covered / block = $(block_name) : base_covered_fraction"].array

        @views linear_covered_fraction_per_found_module_per_block_metacell =
            daf["/ module & is_found ; block = $(block_name) / metacell & block = $(block_name) : $(covered_fraction)"].array

        offset_linear_covered_fraction_per_found_module_per_block_metacell =
            linear_covered_fraction_per_found_module_per_block_metacell .-
            base_covered_fraction_per_found_module_of_block

        approximated_linear_covered_fraction_per_gene_per_metacell[covered_genes_indices, block_metacell_indices] .=
            max.(
                (
                    coefficient_per_covered_gene_per_found_module *
                    offset_linear_covered_fraction_per_found_module_per_block_metacell
                ) .+ base_covered_fraction_per_covered_gene_of_block,
                0.0,
            )
    end

    set_matrix!(
        daf,
        "gene",
        "metacell",
        "approximated_$(covered_fraction)",
        bestify(approximated_linear_covered_fraction_per_gene_per_metacell);
        overwrite,
    )
    return nothing
end

"""
    function compute_metacells_genes_approximated_log_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the
linear model. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        metacell_gene_approximated_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_approximated_log_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_approximated_log_linear_covered_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    approximated_linear_covered_fraction_per_gene_per_metacell =
        get_matrix(daf, "gene", "metacell", "approximated_linear_covered_fraction")
    approximated_log_linear_covered_fraction_per_gene_per_metacell =
        log2.(approximated_linear_covered_fraction_per_gene_per_metacell .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "gene",
        "metacell",
        "approximated_log_linear_covered_fraction",
        approximated_log_linear_covered_fraction_per_gene_per_metacell;
        overwrite,
    );
    return nothing
end

"""
    function compute_metacells_genes_approximated_log_scaled_linear_covered_fractions!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the
linear model, scaled by divergence. This adds the `gene_fraction_regularization` to deal with zero fractions.
"""
@logged @computation Contract(
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        metacell_gene_approximated_scaled_linear_covered_fraction_matrix(RequiredInput),
        metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_metacells_genes_approximated_log_scaled_linear_covered_fractions!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_metacells_genes_log_linear_fractions!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    approximated_scaled_linear_covered_fraction_per_gene_per_metacell =
        get_matrix(daf, "gene", "metacell", "approximated_scaled_linear_covered_fraction")
    approximated_log_scaled_linear_covered_fraction_per_gene_per_metacell =
        log2.(approximated_scaled_linear_covered_fraction_per_gene_per_metacell .+ gene_fraction_regularization)
    set_matrix!(
        daf,
        "gene",
        "metacell",
        "approximated_log_scaled_linear_covered_fraction",
        approximated_log_scaled_linear_covered_fraction_per_gene_per_metacell;
        overwrite,
    );
    return nothing
end

"""
    function compute_blocks_modules_block_RMSE!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block,
evaluated at the block metacells.

!!! note

    This compares the virtual [`metacell_gene_approximated_covered_fraction_matrix`](@ref) and
    [`metacell_gene_covered_fraction_matrix`](@ref). You will need an `adapter` to map these to concrete fractions
    (geomean, linear, scaled, ...).
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_gene_approximated_covered_fraction_matrix(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_modules_block_RMSE_vector(GuaranteedOutput),
    ],
) function compute_blocks_modules_block_RMSE!(daf::DafWriter; overwrite::Bool = false)::Nothing
    estimated_covered_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "covered_fraction")
    approximated_covered_fraction_per_metacell_per_gene =
        get_matrix(daf, "metacell", "gene", "approximated_covered_fraction")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    covered_gene_indices = findall(get_vector(daf, "gene", "is_covered"))

    RMSE_per_block = Vector{Float32}(undef, n_blocks)
    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        block_metacell_indices = daf["/ metacell & block = $(block_name) : index"].array
        n_metacells_in_block = length(block_metacell_indices)

        @views estimated_covered_fraction_per_block_metacell_per_covered =
            estimated_covered_fraction_per_metacell_per_gene[block_metacell_indices, covered_gene_indices]

        @views approximated_covered_fraction_per_block_metacell_per_covered =
            approximated_covered_fraction_per_metacell_per_gene[block_metacell_indices, covered_gene_indices]

        RMSE = sqrt.(
            mean(
                (
                    approximated_covered_fraction_per_block_metacell_per_covered .-
                    estimated_covered_fraction_per_block_metacell_per_covered
                ) .^ 2,
            ),
        )

        @debug "- $(block_name): Metacell: $(n_metacells_in_block) RMSE: $(RMSE)"
        RMSE_per_block[block_index] = RMSE
    end

    @debug "Mean Block RMSE: $(mean(RMSE_per_block[.!isnan.(RMSE_per_block)]))"
    set_vector!(daf, "block", "modules_block_RMSE", RMSE_per_block; overwrite)
    return nothing
end

"""
    function compute_blocks_modules_neighborhood_RMSE!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block,
evaluated at the neighborhood metacells. This internally computes the approximation of the neighborhood metacells based
on the block's model.

!!! note

    This uses the [`metacell_gene_covered_fraction_matrix`](@ref). You will need an `adapter` to map it to a concrete
    fraction (geomean, linear, scaled, ...).
"""
@logged @computation Contract(
    axes = [
        block_axis(RequiredInput),
        metacell_axis(RequiredInput),
        gene_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_module_is_found_matrix(RequiredInput),
        block_module_gene_covered_coefficient_tensor(RequiredInput),
        block_module_base_covered_fraction_matrix(RequiredInput),
        block_gene_base_covered_fraction_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        block_module_gene_covered_coefficient_tensor(RequiredInput),
        block_modules_neighborhood_RMSE_vector(GuaranteedOutput),
    ],
) function compute_blocks_modules_neighborhood_RMSE!(daf::DafWriter; overwrite::Bool = false)::Nothing
    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    covered_fraction_per_covered_per_metacell = daf["/ gene & is_covered / metacell : covered_fraction"].array

    RMSE_per_block = Vector{Float32}(undef, n_blocks)
    @threads :greedy for block_index in 1:n_blocks
        block_name = name_per_block[block_index]
        neighborhood_metacell_indices = daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array
        n_neighborhood_metacells = length(neighborhood_metacell_indices)

        is_found_per_module = daf["/ block = $(block_name) / module : is_found"].array
        found_module_indices = findall(is_found_per_module)
        n_found_modules = length(found_module_indices)
        @assert n_found_modules > 0

        base_covered_fraction_per_module = daf["/ block = $(block_name) / module : base_covered_fraction"].array

        base_covered_fraction_per_covered =
            daf["/ block = $(block_name) / gene & is_covered : base_covered_fraction"].array

        module_index_per_covered = daf["/ block = $(block_name) / gene & is_covered : module => index"].array

        offset_covered_fraction_per_neighborhood_metacell_per_found_module =
            Matrix{Float32}(undef, n_neighborhood_metacells, n_found_modules)
        for (found_module_position, found_module_index) in enumerate(found_module_indices)
            covered_module_genes_positions = findall(module_index_per_covered .== found_module_index)
            @assert length(covered_module_genes_positions) > 0

            @views covered_fraction_per_covered_module_gene_per_neighborhood_metacell =
                covered_fraction_per_covered_per_metacell[covered_module_genes_positions, neighborhood_metacell_indices]

            offset_covered_fraction_per_neighborhood_metacell_per_found_module[:, found_module_position] .=
                vec(sum(covered_fraction_per_covered_module_gene_per_neighborhood_metacell; dims = 1)) .-
                base_covered_fraction_per_module[found_module_index]
        end

        coefficient_per_module_per_covered =
            daf["/ module / gene & is_covered : $(block_name)_covered_coefficient"].array
        coefficient_per_found_module_per_covered =
            coefficient_per_module_per_covered = coefficient_per_module_per_covered[found_module_indices, :]

        predicted_covered_fraction_per_covered_per_neighborhood_metacell = max.(
            (
                transpose(coefficient_per_found_module_per_covered) *
                transpose(offset_covered_fraction_per_neighborhood_metacell_per_found_module)
            ) .+ base_covered_fraction_per_covered,
            0.0,
        )

        estimated_covered_fraction_per_covered_per_neighborhood_metacell =
            daf["/ gene & is_covered / metacell & block => is_in_neighborhood ;= $(block_name) : covered_fraction"].array

        RMSE_per_block[block_index] = sqrt(
            mean(
                (
                    predicted_covered_fraction_per_covered_per_neighborhood_metacell .-
                    estimated_covered_fraction_per_covered_per_neighborhood_metacell
                ) .^ 2,
            ),
        )
    end

    @debug "Mean Neighborhood RMSE: $(mean(RMSE_per_block[.!isnan.(RMSE_per_block)]))"
    set_vector!(daf, "block", "modules_neighborhood_RMSE", RMSE_per_block; overwrite)
    return nothing
end

end  # module
