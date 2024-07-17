"""
Compute gene programs to model the biological cell state.
"""
module Programs

using ..Contracts
using ..Defaults
using Base.Threads
using Convex
using Daf
using Daf.GenericTypes
using DataStructures
using Random
using SCS
using SparseArrays
using StatsBase
using VectorizedStatistics

import Daf

struct Problem
    factors_mask::Union{AbstractVector{Bool}, BitVector}
    genes_mask::Union{AbstractVector{Bool}, BitVector}
    metacells_mask::Union{AbstractVector{Bool}, BitVector}
    values_of_genes_of_metacells::AbstractMatrix{<:AbstractFloat}
end

struct Solution
    selected_factors_indices::AbstractVector{<:Integer}
    selected_factors_mask::Union{AbstractVector{Bool}, BitVector}
    coefficients_of_genes_of_selected_factors::AbstractMatrix{<:AbstractFloat}
    left_out_r2::AbstractFloat
end

function solve_without_overfitting(
    problem::Problem;
    min_abs_coefficient::AbstractFloat,
    max_factors::Integer,
    names_of_genes::AbstractVector{<:AbstractString},
)::Solution
    @assert min_abs_coefficient >= 0
    @assert max_factors > 0

    n_metacells, n_genes = size(problem.values_of_genes_of_metacells)
    @assert length(problem.metacells_mask) == n_metacells
    @assert length(problem.factors_mask) == n_genes
    @assert length(problem.genes_mask) == n_genes
    n_factors = sum(problem.factors_mask)
    @assert n_factors > 0

    @debug "solve_without_overfitting n_metacells: $(sum(problem.metacells_mask)) n_genes: $(sum(problem.genes_mask)) n_factors: $(sum(problem.factors_mask))"

    masked_values_of_genes_of_metacells =
        problem.values_of_genes_of_metacells[problem.metacells_mask, problem.genes_mask]
    n_masked_metacells, n_masked_genes = size(masked_values_of_genes_of_metacells)
    @assert n_masked_genes <= n_genes

    masked_values_of_factors_of_metacells =
        problem.values_of_genes_of_metacells[problem.metacells_mask, problem.factors_mask]
    m_masked_metacells, n_masked_factors = size(masked_values_of_factors_of_metacells)
    @assert m_masked_metacells == n_masked_metacells
    @assert n_masked_factors <= n_genes

    masked_correlation_between_factors_and_genes =
        abs.(cor(masked_values_of_factors_of_metacells, masked_values_of_genes_of_metacells))
    masked_correlation_between_factors_and_genes[isnan.(masked_correlation_between_factors_and_genes)] .= 0
    @assert size(masked_correlation_between_factors_and_genes) == (n_factors, n_genes)

    masked_factors_indices = findall(problem.factors_mask)

    masked_selected_factors_mask = zeros(Bool, n_masked_factors)
    selected_factors_mask = zeros(Bool, n_genes)
    selected_solution = nothing
    n_selected_factors = 0

    metacell_indices = collect(1:n_metacells)
    masked_metacell_indices = metacell_indices[problem.metacells_mask]

    max_factors = min(max_factors, n_factors)

    selected_factors_indices = Vector{Int32}()

    # TODOX Certificates for neighborhoods
    # TODOX Dump few neighborhoods
    # TODOX Dump all-type

    best_left_out_r2 = -1.0
    best_solution = nothing

    while length(selected_factors_indices) < max_factors
        if n_genes < 5
            masked_correlation_of_factors = vec(sum(masked_correlation_between_factors_and_genes; dims = 2))
        else
            quantile = (n_genes - 4) / (n_genes - 1)  # TODOX - 3rd best (not including the factor itself).
            masked_correlation_of_factors = vec(vquantile(masked_correlation_between_factors_and_genes, quantile; dims = 2))
        end
        masked_correlation_of_factors[masked_selected_factors_mask] .= -1
        @assert length(masked_correlation_of_factors) == n_masked_factors
        masked_correlation_between_factors_and_genes = nothing

        best_factor_index = nothing
        best_masked_factor_index = nothing

        n_tested = 0
        while best_factor_index === nothing || n_tested < 10  # TODOX
            masked_factor_index = argmax(masked_correlation_of_factors)
            if masked_correlation_of_factors[masked_factor_index] < 0
                break
            end

            factor_index = masked_factors_indices[masked_factor_index]

            @assert !selected_factors_mask[factor_index]
            selected_factors_mask[factor_index] = true

            left_out_r2 = estimate_leave_one_out_r2(
                problem;
                min_abs_coefficient = min_abs_coefficient,
                selected_factors_mask = selected_factors_mask,
                masked_metacell_indices = masked_metacell_indices,
            )

            selected_factors_mask[factor_index] = false
            masked_correlation_of_factors[masked_factor_index] = -1

            if left_out_r2 > best_left_out_r2
                @debug "factor number: $(length(selected_factors_indices)) factor: $(names_of_genes[factor_index]) => left_out_r2: $(left_out_r2)"
                best_left_out_r2 = left_out_r2
                best_factor_index = factor_index
                best_masked_factor_index = masked_factor_index
                n_tested = 0
            else
                @debug "factor number: $(length(selected_factors_indices)) factor: $(names_of_genes[factor_index]) ... left_out_r2: $(left_out_r2)"
                n_tested += 1
            end
        end

        if best_factor_index === nothing
            @assert best_solution !== nothing
            return best_solution
        end

        @assert best_left_out_r2 !== nothing
        @assert best_masked_factor_index !== nothing
        push!(selected_factors_indices, best_factor_index)

        @assert !selected_factors_mask[best_factor_index]
        selected_factors_mask[best_factor_index] = true

        @assert !masked_selected_factors_mask[best_masked_factor_index]
        masked_selected_factors_mask[best_masked_factor_index] = true

        selected_problem = Problem(
            selected_factors_mask,
            problem.genes_mask,
            problem.metacells_mask,
            problem.values_of_genes_of_metacells,
        )

        coefficients_of_genes_of_selected_factors,
        masked_measured_values_of_genes_of_metacells,
        masked_computed_values_of_genes_of_metacells =
            compute_least_squares(selected_problem; min_abs_coefficient = min_abs_coefficient)

        best_solution = Solution(
            copy_array(selected_factors_indices),
            copy_array(selected_factors_mask),
            coefficients_of_genes_of_selected_factors,
            best_left_out_r2
        )

        masked_residual_values_of_genes_of_metacells =
            masked_measured_values_of_genes_of_metacells .- masked_computed_values_of_genes_of_metacells

        masked_correlation_between_factors_and_genes =
            abs.(cor(masked_values_of_factors_of_metacells, masked_residual_values_of_genes_of_metacells))
        masked_correlation_between_factors_and_genes[isnan.(masked_correlation_between_factors_and_genes)] .= 0
    @assert size(masked_correlation_between_factors_and_genes) == (n_factors, n_genes)
    end

    @assert best_solution !== nothing
    return best_solution
end

function estimate_leave_one_out_r2(
    problem::Problem;
    min_abs_coefficient::AbstractFloat,
    selected_factors_mask::Union{AbstractVector{Bool}, BitVector},
    masked_metacell_indices::AbstractVector{<:Integer},
)::AbstractFloat
    n_metacells, n_genes = size(problem.values_of_genes_of_metacells)
    @assert length(problem.metacells_mask) == n_metacells
    @assert length(problem.factors_mask) == n_genes
    @assert length(problem.genes_mask) == n_genes
    n_factors = sum(problem.factors_mask)
    @assert n_factors > 0
    @assert length(selected_factors_mask) == n_genes
    n_selected_factors = sum(selected_factors_mask)
    @assert n_selected_factors <= n_factors

    n_masked_genes = sum(problem.genes_mask)
    n_masked_metacells = length(masked_metacell_indices)

    r2s_of_left_outs = Vector{Float64}(undef, n_masked_metacells)
    @threads for left_out_position in 1:n_masked_metacells
        left_out_index = masked_metacell_indices[left_out_position]
        metacells_mask = copy_array(problem.metacells_mask)
        @assert metacells_mask[left_out_index]
        metacells_mask[left_out_index] = false
        genes_mask = problem.genes_mask .& .! selected_factors_mask
        @assert sum(genes_mask) == n_masked_genes - n_selected_factors

        left_one_out_problem = Problem(
            selected_factors_mask,
            genes_mask,
            metacells_mask,
            problem.values_of_genes_of_metacells,
        )

        coefficients_of_genes_of_selected_factors,
        masked_measured_values_of_genes_of_metacells,
        masked_computed_values_of_genes_of_metacells =
            compute_least_squares(left_one_out_problem; min_abs_coefficient = min_abs_coefficient)

        @assert size(coefficients_of_genes_of_selected_factors) == (n_selected_factors, n_masked_genes - n_selected_factors)
        @assert size(masked_measured_values_of_genes_of_metacells) == (n_masked_metacells - 1, n_masked_genes - n_selected_factors)
        @assert size(masked_computed_values_of_genes_of_metacells) == (n_masked_metacells - 1, n_masked_genes - n_selected_factors)

        values_of_selected_factors_of_left_out = vec(problem.values_of_genes_of_metacells[left_out_index, selected_factors_mask])
        @assert length(values_of_selected_factors_of_left_out) == n_selected_factors

        measured_of_genes_of_left_out = vec(problem.values_of_genes_of_metacells[left_out_index, genes_mask])
        @assert length(measured_of_genes_of_left_out) == n_masked_genes - n_selected_factors

        computed_of_genes_of_left_out = transpose(values_of_selected_factors_of_left_out) * coefficients_of_genes_of_selected_factors
        @assert length(computed_of_genes_of_left_out) == n_masked_genes - n_selected_factors

        correlation = vcor(vec(computed_of_genes_of_left_out), vec(measured_of_genes_of_left_out))
        if isnan(correlation)
            r2s_of_left_outs[left_out_position] = -1
        else
            r2s_of_left_outs[left_out_position] = correlation * correlation
        end
    end

    left_out_r2 = vmedian(r2s_of_left_outs[r2s_of_left_outs .> 0])
    return left_out_r2
end

function compute_least_squares(  # untested
    problem::Problem;
    min_abs_coefficient::AbstractFloat,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}}
    n_metacells, n_genes = size(problem.values_of_genes_of_metacells)
    @assert length(problem.metacells_mask) == n_metacells
    @assert length(problem.factors_mask) == n_genes
    @assert length(problem.genes_mask) == n_genes

    masked_values_of_factors_of_metacells =
        problem.values_of_genes_of_metacells[problem.metacells_mask, problem.factors_mask]
    masked_measured_values_of_genes_of_metacells =
        problem.values_of_genes_of_metacells[problem.metacells_mask, problem.genes_mask]

    masked_coefficients_of_genes_of_factors =
        masked_values_of_factors_of_metacells \ masked_measured_values_of_genes_of_metacells
    masked_coefficients_of_genes_of_factors[abs.(masked_coefficients_of_genes_of_factors) .< min_abs_coefficient] .= 0.0

    masked_computed_values_of_genes_of_metacells =
        masked_values_of_factors_of_metacells * masked_coefficients_of_genes_of_factors
    @assert size(masked_computed_values_of_genes_of_metacells) == (sum(problem.metacells_mask), sum(problem.genes_mask))

    return (
        masked_coefficients_of_genes_of_factors,
        masked_measured_values_of_genes_of_metacells,
        masked_computed_values_of_genes_of_metacells,
    )
end

end
