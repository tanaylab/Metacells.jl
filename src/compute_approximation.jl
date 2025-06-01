"""
Compute a local linear approximation to the manifold based on gene modules.
"""
module ComputeApproximation

export compute_blocks_modules_base_covered_fractions!
export compute_approximation!

using Base.Threads
using DataAxesFormats
using Statistics
using TanayLabUtilities
using Random

using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.block_block_is_in_environment_matrix
import Metacells.Contracts.block_block_is_in_neighborhood_matrix
import Metacells.Contracts.block_gene_base_covered_fraction_matrix
import Metacells.Contracts.block_gene_is_environment_marker_matrix
import Metacells.Contracts.block_gene_module_matrix
import Metacells.Contracts.block_module_base_covered_fraction_matrix
import Metacells.Contracts.block_module_gene_covered_coefficient_tensor
import Metacells.Contracts.block_module_is_found_matrix
import Metacells.Contracts.block_module_min_gene_correlation_matrix
import Metacells.Contracts.block_modules_neighborhood_RMSE_vector
import Metacells.Contracts.block_modules_neighborhood_XRMSE_vector
import Metacells.Contracts.cell_axis
import Metacells.Contracts.cell_gene_UMIs_matrix
import Metacells.Contracts.cell_metacell_vector
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_block_vector
import Metacells.Contracts.metacell_gene_covered_fraction_matrix
import Metacells.Contracts.metacell_gene_UMIs_matrix
import Metacells.Contracts.metacell_n_cells_vector
import Metacells.Contracts.module_axis

@kwdef struct ReusableMatrices
    coefficient_per_module_per_gene::Matrix{Float32}
end

function TanayLabUtilities.reset_reusable_storage!(reusable_matrices::ReusableMatrices)::Nothing
    reusable_matrices.coefficient_per_module_per_gene .= 0
    return nothing
end

"""
    function compute_blocks_modules_base_covered_fractions!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

The base fractions used by the gene modules analysis. This is based on the per-gene `base_covered_fraction`.

$(CONTRACT)
"""
@logged @computation Contract(;
    axes = [block_axis(RequiredInput), module_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_gene_base_covered_fraction_matrix(RequiredInput),
        block_module_base_covered_fraction_matrix(GuaranteedOutput),
    ],
) function compute_blocks_modules_base_covered_fractions!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_modules = axis_length(daf, "module")
    n_blocks = axis_length(daf, "block")

    is_covered_per_gene = get_vector(daf, "gene", "is_covered")
    module_index_per_gene_per_block = daf["/ gene / block : module => index"].array
    base_covered_fraction_per_gene_per_block = get_matrix(daf, "gene", "block", "base_covered_fraction")

    base_covered_fraction_per_module_per_block = zeros(Float32, n_modules, n_blocks)

    @threads :greedy for block_index in 1:n_blocks
        @views module_index_per_gene_of_block = module_index_per_gene_per_block[:, block_index]
        @views base_covered_fraction_per_gene_of_block = base_covered_fraction_per_gene_per_block[:, block_index]

        for module_index in 1:n_modules
            module_genes_mask = (module_index_per_gene_of_block .== module_index) .& is_covered_per_gene
            base_covered_fraction_per_module_per_block[module_index, block_index] =
                sum(base_covered_fraction_per_gene_of_block[module_genes_mask]; init = Float32(0))
        end
    end

    set_matrix!(
        daf,
        "module",
        "block",
        "base_covered_fraction",
        bestify(base_covered_fraction_per_module_per_block);
        overwrite,
    )
    return nothing
end

"""
    function compute_approximation!(
        daf::DafWriter;
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

For each block, compute a local linear model for its environment, based on the found gene modules. This is a simple
least-squares approximation predicting the fraction of the covered genes based on the fraction of the gene modules
(using the `base_covered_fraction` of the modules in the environment). This is assumed to over-fit the solution, so we use
cross-validation (in `cross_validation_parts`) using subsets of the metacells to compute both the RMSE and XRMSE to
estimate this.

!!! note

    This uses the virtual [`metacell_gene_covered_fraction_matrix`](@ref). You will need an `adapter` to map these to
    concrete fractions (geomean, linear, scaled, ...).

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [
        gene_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(RequiredInput),
        module_axis(RequiredInput),
    ],
    data = [
        gene_is_covered_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_covered_fraction_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_gene_module_matrix(RequiredInput),
        block_gene_base_covered_fraction_matrix(RequiredInput),
        block_module_base_covered_fraction_matrix(RequiredInput),
        block_module_gene_covered_coefficient_tensor(GuaranteedOutput),
        block_modules_neighborhood_RMSE_vector(GuaranteedOutput),
        block_modules_neighborhood_XRMSE_vector(GuaranteedOutput),
    ],
) function compute_approximation!(
    daf::DafWriter;
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert cross_validation_parts > 1
    n_genes = axis_length(daf, "gene")
    n_modules = axis_length(daf, "module")
    n_blocks = axis_length(daf, "block")

    name_per_block = axis_vector(daf, "block")
    covered_indices = daf["/ gene & is_covered : index"].array
    covered_fraction_per_covered_per_metacell = daf["/ gene & is_covered / metacell : covered_fraction"].array

    modules_neighborhood_RMSE_per_block = Vector{Float32}(undef, n_blocks)
    modules_neighborhood_XRMSE_per_block = Vector{Float32}(undef, n_blocks)

    reusable_storage = ReusableStorage() do
        return ReusableMatrices(; coefficient_per_module_per_gene = zeros(Float32, n_modules, n_genes))
    end

    progress_counter = Atomic{Int}(0)
    parallel_loop_with_rng(1:n_blocks; rng) do block_index, rng
        block_name = name_per_block[block_index]

        n_found_modules, rmse, xrmse = with_reusable(reusable_storage) do reusable_matrices
            return compute_block_approximation(
                daf;
                block_name,
                covered_fraction_per_covered_per_metacell,
                covered_indices,
                cross_validation_parts,
                rng,
                overwrite,
                reusable_matrices.coefficient_per_module_per_gene,
            )
        end

        modules_neighborhood_RMSE_per_block[block_index] = rmse
        modules_neighborhood_XRMSE_per_block[block_index] = xrmse

        counter = atomic_add!(progress_counter, 1)
        @debug (
            "- Environment: $(block_name) ($(percent(counter + 1, n_blocks)))" *
            " Gene Modules: $(n_found_modules)" *
            " RMSE: $(rmse)" *
            " XRMSE: $(xrmse)"
        )

        return nothing
    end

    @debug "Mean RMSE: $(mean(modules_neighborhood_RMSE_per_block)) XRMSE: $(mean(modules_neighborhood_XRMSE_per_block))"
    set_vector!(daf, "block", "modules_neighborhood_RMSE", modules_neighborhood_RMSE_per_block; overwrite)
    set_vector!(daf, "block", "modules_neighborhood_XRMSE", modules_neighborhood_XRMSE_per_block; overwrite)
    return nothing
end

function compute_block_approximation(
    daf::DafWriter;
    block_name::AbstractString,
    covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    cross_validation_parts::Integer,
    covered_indices::AbstractVector{<:Integer},
    rng::AbstractRNG,
    overwrite::Bool,
    coefficient_per_module_per_gene::Matrix{Float32},
)::Tuple{Integer, AbstractFloat, AbstractFloat}
    metacell_indices_of_neighborhood = daf["/ metacell & block => is_in_neighborhood ;= $(block_name) : index"].array
    metacell_indices_of_environment = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"].array
    additional_indices_of_environment = setdiff(metacell_indices_of_environment, metacell_indices_of_neighborhood)

    is_found_per_module = daf["/ block = $(block_name) / module : is_found"].array
    found_module_indices = findall(is_found_per_module)
    n_found_modules = length(found_module_indices)
    @assert n_found_modules > 0

    cross_validation_indices =
        pick_cross_validation_indices(; full_indices = metacell_indices_of_neighborhood, cross_validation_parts, rng)

    rmse_per_part = Vector{Float32}(undef, cross_validation_parts)
    for part_index in 1:cross_validation_parts
        rmse_of_part, _ = compute_environment_model(
            daf;
            block_name,
            covered_fraction_per_covered_per_metacell,
            found_module_indices,
            train_metacell_indices = vcat(
                cross_validation_indices.train_indices_per_part[part_index],
                additional_indices_of_environment,
            ),
            test_metacell_indices = cross_validation_indices.test_indices_per_part[part_index],
        )
        rmse_per_part[part_index] = rmse_of_part
    end
    xrmse = mean(rmse_per_part)

    rmse, coefficient_per_found_module_per_covered = compute_environment_model(
        daf;
        block_name,
        covered_fraction_per_covered_per_metacell,
        found_module_indices,
        train_metacell_indices = metacell_indices_of_environment,
        test_metacell_indices = metacell_indices_of_neighborhood,
    )

    coefficient_per_module_per_gene[found_module_indices, covered_indices] .= coefficient_per_found_module_per_covered
    set_matrix!(
        daf,
        "module",
        "gene",
        "$(block_name)_covered_coefficient",
        bestify(coefficient_per_module_per_gene);
        overwrite,
    )

    return (n_found_modules, rmse, xrmse)
end

function compute_environment_model(
    daf::DafWriter;
    block_name::AbstractString,
    covered_fraction_per_covered_per_metacell::AbstractMatrix{<:AbstractFloat},
    found_module_indices::AbstractVector{<:Integer},
    train_metacell_indices::AbstractVector{<:Integer},
    test_metacell_indices::AbstractVector{<:Integer},
)::Tuple{AbstractFloat, AbstractMatrix{<:AbstractFloat}}
    base_covered_fraction_per_module = daf["/ block = $(block_name) / module : base_covered_fraction"].array
    module_index_per_covered = daf["/ block = $(block_name) / gene & is_covered : module => index"].array
    base_covered_fraction_per_covered = daf["/ block = $(block_name) / gene & is_covered : base_covered_fraction"].array

    n_found_modules = length(found_module_indices)
    n_train_metacells = length(train_metacell_indices)
    n_test_metacells = length(test_metacell_indices)
    n_covered_genes = size(covered_fraction_per_covered_per_metacell, 1)
    @assert n_found_modules > 0

    @views covered_fraction_per_covered_per_train_metacell =
        covered_fraction_per_covered_per_metacell[:, train_metacell_indices]
    offset_covered_fraction_per_covered_per_train_metacell =
        covered_fraction_per_covered_per_train_metacell .- base_covered_fraction_per_covered

    offset_covered_fraction_per_train_metacell_per_found_module =
        Matrix{Float32}(undef, n_train_metacells, n_found_modules)
    if train_metacell_indices === test_metacell_indices
        offset_covered_fraction_per_test_metacell_per_found_module =
            offset_covered_fraction_per_train_metacell_per_found_module
        offset_covered_fraction_per_covered_per_test_metacell = offset_covered_fraction_per_covered_per_train_metacell
    else
        offset_covered_fraction_per_test_metacell_per_found_module =
            Matrix{Float32}(undef, n_test_metacells, n_found_modules)
        @views covered_fraction_per_covered_per_test_metacell =
            covered_fraction_per_covered_per_metacell[:, test_metacell_indices]
        offset_covered_fraction_per_covered_per_test_metacell =
            covered_fraction_per_covered_per_test_metacell .- base_covered_fraction_per_covered
    end

    for (found_module_position, found_module_index) in enumerate(found_module_indices)
        covered_module_genes_positions = findall(module_index_per_covered .== found_module_index)
        @assert length(covered_module_genes_positions) > 0

        @views covered_fraction_per_covered_module_gene_per_train_metacell =
            covered_fraction_per_covered_per_train_metacell[covered_module_genes_positions, :]
        offset_covered_fraction_per_train_metacell_per_found_module[:, found_module_position] .=
            vec(sum(covered_fraction_per_covered_module_gene_per_train_metacell; dims = 1)) .-
            base_covered_fraction_per_module[found_module_index]

        if offset_covered_fraction_per_test_metacell_per_found_module !==
           offset_covered_fraction_per_train_metacell_per_found_module
            @views covered_fraction_per_covered_module_gene_per_test_metacell =
                covered_fraction_per_covered_per_test_metacell[covered_module_genes_positions, :]
            offset_covered_fraction_per_test_metacell_per_found_module[:, found_module_position] .=
                vec(sum(covered_fraction_per_covered_module_gene_per_test_metacell; dims = 1)) .-
                base_covered_fraction_per_module[found_module_index]
        end
    end

    coefficient_per_found_module_per_covered =
        offset_covered_fraction_per_train_metacell_per_found_module \
        transpose(offset_covered_fraction_per_covered_per_train_metacell)
    @assert_matrix(coefficient_per_found_module_per_covered, n_found_modules, n_covered_genes)

    predicted_covered_fraction_per_covered_per_test_metacell = max.(
        (
            transpose(coefficient_per_found_module_per_covered) *
            transpose(offset_covered_fraction_per_test_metacell_per_found_module)
        ) .+ base_covered_fraction_per_covered,
        0.0,
    )
    @assert_matrix(predicted_covered_fraction_per_covered_per_test_metacell, n_covered_genes, n_test_metacells)

    rmse = sqrt(
        mean(
            (
                predicted_covered_fraction_per_covered_per_test_metacell .-
                covered_fraction_per_covered_per_test_metacell
            ) .^ 2,
        ),
    )

    return (rmse, coefficient_per_found_module_per_covered)
end

end  # module
