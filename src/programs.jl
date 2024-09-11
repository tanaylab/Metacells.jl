"""
Approximate the manifold of actual cell states (captured by metacells) using linear programs in each local region.
"""
module Programs

export compute_factor_priority_of_genes!
export compute_global_predictive_factors!
export compute_blocks!
export compute_blocks_vicinities!
export compute_local_predictive_factors!
export compute_programs!

using ..Contracts
using ..Defaults
using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using Distributions
using MultivariateStats
using NonNegLeastSquares
using Printf
using Random
using SparseArrays
using Statistics
using VectorizedStatistics

import Random.default_rng

@kwdef struct Context
    daf::DafReader
    n_genes::Integer
    n_metacells::Integer
    names_of_genes::AbstractVector{<:AbstractString}
    divergence_of_genes::AbstractVector{<:AbstractFloat}
    fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
    log_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_factor_priority_of_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        max_principal_components::Integer = $(DEFAULT.max_principal_components),
        factors_per_principal_component::Real = $(DEFAULT.factors_per_principal_component),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given the transcription factors in the data, identify an ordered subset of these to use as candidates for predicting the
rest of the genes. To achieve this, we run PCA analysis with up to `max_principal_components` on the log (base 2) of the
gene expressions (using the `gene_fraction_regularization`). We then rank each gene in each principal component (based
on the absolute coefficient value) and compute the minimal rank of each gene (giving lower weight to genes in later
components). We prioritize the transcription factors according to this minimal rank, choosing the top
`factors_per_principal_component` times the number of principal components.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_transcription_factor_vector(RequiredInput),
        gene_is_forbidden_factor_vector(OptionalInput),
        gene_factor_priority_vector(GuaranteedOutput),
    ],
) function compute_factor_priority_of_genes!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 2 * GENE_FRACTION_REGULARIZATION,
    max_principal_components::Integer = 40,
    factors_per_principal_component::Real = 2,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert max_principal_components > 0
    @assert factors_per_principal_component > 0

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)

    factor_genes_mask =
        get_vector(daf, "gene", "is_transcription_factor") .&
        .!get_vector(daf, "gene", "is_forbidden_factor"; default = false)
    factor_genes = findall(factor_genes_mask)
    @assert length(factor_genes) > 0

    log_fractions_of_genes_in_metacells = context.log_fractions_of_genes_in_metacells
    divergence_of_genes = context.divergence_of_genes
    log_fractions_of_genes_in_metacells .*= 1.0 .- divergence_of_genes

    pca = fit(PCA, log_fractions_of_genes_in_metacells; maxoutdim = max_principal_components)

    n_principal_components = size(pca, 2)

    coefficients_of_genes_of_principal_components = loadings(pca)
    @assert_matrix(coefficients_of_genes_of_principal_components, context.n_genes, n_principal_components, Columns)

    coefficients_of_factor_genes_of_principal_components =
        coefficients_of_genes_of_principal_components[factor_genes, :]
    @assert_matrix(
        coefficients_of_factor_genes_of_principal_components,
        length(factor_genes),
        n_principal_components,
        Columns,
    )

    minimal_rank_of_factor_genes = fill(length(factor_genes) * 100.0, length(factor_genes))
    for principal_component_index in 1:n_principal_components
        abs_coefficients_of_factors_in_principal_component =
            abs.(coefficients_of_factor_genes_of_principal_components[:, principal_component_index])

        factor_genes_order = sortperm(abs_coefficients_of_factors_in_principal_component; rev = true)
        rank_of_factor_genes = invperm(factor_genes_order) .* (1 + principal_component_index / n_principal_components)
        minimal_rank_of_factor_genes .= min.(minimal_rank_of_factor_genes, rank_of_factor_genes)
    end

    factor_genes_order = sortperm(minimal_rank_of_factor_genes)
    rank_of_factor_genes = invperm(factor_genes_order)

    n_top_genes = min(Int(round(n_principal_components * factors_per_principal_component)), length(factor_genes))
    top_genes = factor_genes[factor_genes_order[1:n_top_genes]]
    @debug "top $(n_top_genes) ranked transcription factors: $(join(context.names_of_genes[top_genes], ", "))"

    factor_priority_of_genes = zeros(UInt16, context.n_genes)
    factor_priority_of_genes[factor_genes] .= max.((n_top_genes + 1) .- rank_of_factor_genes, 0)
    set_vector!(daf, "gene", "factor_priority", factor_priority_of_genes; overwrite = overwrite)

    return nothing
end

function load_context(daf::DafReader; gene_fraction_regularization::AbstractFloat)::Context  # untested
    n_genes = axis_length(daf, "gene")
    @assert n_genes > 0

    n_metacells = axis_length(daf, "metacell")
    @assert n_metacells > 0

    names_of_genes = axis_array(daf, "gene")

    divergence_of_genes = get_vector(daf, "gene", "divergence").array

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    log_fractions_of_genes_in_metacells = log2.(fractions_of_genes_in_metacells .+ gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    return Context(;
        daf = daf,
        n_genes = n_genes,
        n_metacells = n_metacells,
        names_of_genes = names_of_genes,
        divergence_of_genes = divergence_of_genes,
        fractions_of_genes_in_metacells = fractions_of_genes_in_metacells,
        log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells,
    )
end

"""
    function compute_global_predictive_factors!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given a prioritized subset of the factor genes likely to be predictive, identify a subset of these genes which best
predict the values of the rest of the genes, using non-negative least-squares approximation using the log (base 2) of
the genes expression (using the `gene_fraction_regularization`). To avoid over-fitting, we use `cross_validation_parts`.

We require each additional used factor to improve (reduce) the RMS of the cross validation error by at least
`min_rms_improvement` (on average). Factors with higher prioriry are allowed a lower improvement and factors with a
lower priority require higher improvement; the range of this is `rms_priority_improvement` times the
`min_rms_improvement`. That is, if this range is 1, then the highest priority factor requires 0.5 times the minimal
improvement, and the lowest priority factor requires 1.5 times the minimal improvement.

Since this prediction is global for the whole data set, we do not expect it to be useful as of itself. We therefore
intentionally do not store the coefficients. Instead we only store the mask of chosen predictive factors.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_factor_priority_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(GuaranteedOutput),
        gene_global_rms_vector(GuaranteedOutput),
        gene_global_r2_vector(GuaranteedOutput),
    ],
) function compute_global_predictive_factors!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_factor_priority_of_genes!,
        :gene_fraction_regularization,
    ),
    min_rms_improvement::AbstractFloat = 1e-2,
    rms_priority_improvement::Real = 1,
    cross_validation_parts::Integer = 5,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert cross_validation_parts > 1
    @assert min_rms_improvement >= 0
    @assert 0 <= rms_priority_improvement <= 2

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)
    candidate_factors = load_candidate_factors(
        daf;
        min_rms_improvement = min_rms_improvement,
        rms_priority_improvement = rms_priority_improvement,
    )

    global_predictive_genes, cross_validation = identify_predictive_genes!(;
        context = context,
        cross_validation_parts = cross_validation_parts,
        included_genes = 1:(context.n_genes),
        included_metacells = 1:(context.n_metacells),
        core_metacells = 1:(context.n_metacells),
        candidate_factors = candidate_factors,
        rng = rng,
    )

    is_global_predictive_factor_of_genes = zeros(Bool, context.n_genes)
    is_global_predictive_factor_of_genes[global_predictive_genes] .= true

    set_vector!(daf, "gene", "global_rms", cross_validation.rms_of_genes; overwrite = overwrite)
    set_vector!(daf, "gene", "global_r2", cross_validation.r2_of_genes; overwrite = overwrite)
    set_vector!(daf, "gene", "is_global_predictive_factor", is_global_predictive_factor_of_genes; overwrite = overwrite)

    @debug "global predictive factors: [ $(join(context.names_of_genes[global_predictive_genes], ", ")) ]"
    @debug "mean cross-validation genes RMS: $(mean(cross_validation.rms_of_genes))"
    @debug "mean cross-validation genes R^2: $(mean(cross_validation.r2_of_genes))"

    return nothing
end

@kwdef struct CandidateFactors
    ordered_factor_genes::AbstractVector{<:Integer}
    factor_priority_of_genes::AbstractVector{<:Unsigned}
    rms_cost_factor_of_genes::AbstractVector{<:AbstractFloat}
end

function load_candidate_factors(  # untested
    daf::DafReader;
    min_rms_improvement::AbstractFloat,
    rms_priority_improvement::Real,
)::CandidateFactors
    factor_priority_of_genes = get_vector(daf, "gene", "factor_priority")
    ordered_genes_indices = sortperm(factor_priority_of_genes; rev = true)

    n_genes = axis_length(daf, "gene")
    n_ordered_factor_genes = sum(factor_priority_of_genes .> 0)

    ordered_factor_genes = ordered_genes_indices[1:n_ordered_factor_genes]
    rms_cost_factor_of_genes = compute_rms_cost_factor_of_genes(;
        min_rms_improvement = min_rms_improvement,
        rms_priority_improvement = rms_priority_improvement,
        ordered_factor_genes = ordered_factor_genes,
        n_genes = n_genes,
    )

    return CandidateFactors(;
        ordered_factor_genes = ordered_factor_genes,
        factor_priority_of_genes = factor_priority_of_genes,
        rms_cost_factor_of_genes = rms_cost_factor_of_genes,
    )
end

function compute_rms_cost_factor_of_genes(;  # untested
    min_rms_improvement::AbstractFloat,
    rms_priority_improvement::Real,
    ordered_factor_genes::AbstractVector{<:Integer},
    n_genes::Integer,
)::AbstractVector{<:AbstractFloat}
    rms_cost_factor_of_genes = zeros(Float32, n_genes)
    rms_cost_factor_of_genes[ordered_factor_genes] .= (
        1 .+
        min_rms_improvement .* (
            1 .+ (
                rms_priority_improvement .*
                ((0:(length(ordered_factor_genes) - 1)) ./ length(ordered_factor_genes) .- 0.5)
            )
        )
    )
    return rms_cost_factor_of_genes
end

@kwdef struct CrossValidation
    mean_rms::AbstractFloat
    rms_of_genes::AbstractVector{<:AbstractFloat}
    r2_of_genes::AbstractVector{<:AbstractFloat}
end

function identify_predictive_genes!(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    candidate_factors::CandidateFactors,
    rng::AbstractRNG,
)::Tuple{AbstractVector{<:Integer}, CrossValidation}
    predictive_genes = Int32[]
    predictive_cost = Maybe{Float32}[nothing, nothing, nothing]

    op_index = 2
    op_success = [true, true]
    last_added_gene_index = nothing
    last_removed_gene_index = nothing
    cross_validation = nothing
    while sum(op_success) > 0
        op_index = 1 + op_index % 2
        if op_index == 1
            last_added_gene_index, new_cross_validation = try_add_factors!(;
                context = context,
                cross_validation_parts = cross_validation_parts,
                included_genes = included_genes,
                included_metacells = included_metacells,
                core_metacells = core_metacells,
                candidate_factors = candidate_factors,
                rng = rng,
                predictive_genes = predictive_genes,
                predictive_cost = predictive_cost,
                last_removed_gene_index = last_removed_gene_index,
            )
            op_success[op_index] = last_added_gene_index !== nothing
        else
            last_removed_gene_index, new_cross_validation = try_remove_factors!(;
                context = context,
                cross_validation_parts = cross_validation_parts,
                included_genes = included_genes,
                included_metacells = included_metacells,
                core_metacells = core_metacells,
                candidate_factors = candidate_factors,
                rng = rng,
                predictive_genes = predictive_genes,
                predictive_cost = predictive_cost,
                last_added_gene_index = last_added_gene_index,
            )
            op_success[op_index] = last_removed_gene_index !== nothing
        end
        if new_cross_validation !== nothing
            cross_validation = new_cross_validation
        end
    end

    @assert cross_validation !== nothing
    return predictive_genes, cross_validation
end  # NOJET

function try_add_factors!(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    candidate_factors::CandidateFactors,
    rng::AbstractRNG,
    predictive_genes::AbstractVector{<:Integer},
    predictive_cost::Vector{Maybe{Float32}},
    last_removed_gene_index::Maybe{<:Integer},
)::Tuple{Maybe{<:Integer}, Maybe{CrossValidation}}
    last_added_gene_index = nothing
    cross_validation = nothing
    for (test_index, factor_gene_index) in enumerate(candidate_factors.ordered_factor_genes)
        if factor_gene_index != last_removed_gene_index && !(factor_gene_index in predictive_genes)
            new_cross_validation = try_add_factor!(;
                context = context,
                cross_validation_parts = cross_validation_parts,
                included_genes = included_genes,
                included_metacells = included_metacells,
                core_metacells = core_metacells,
                candidate_factors = candidate_factors,
                rng = rng,
                predictive_genes = predictive_genes,
                predictive_cost = predictive_cost,
                last_added_gene_index = last_added_gene_index,
                test_index = test_index,
                factor_gene_index = factor_gene_index,
            )
            if new_cross_validation !== nothing
                cross_validation = new_cross_validation
                last_added_gene_index = factor_gene_index
            end
        end
    end
    return (last_added_gene_index, cross_validation)
end

function try_add_factor!(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    candidate_factors::CandidateFactors,
    rng::AbstractRNG,
    predictive_genes::AbstractVector{<:Integer},
    predictive_cost::Vector{Maybe{Float32}},
    last_added_gene_index::Maybe{<:Integer},
    test_index::Integer,
    factor_gene_index::Integer,
)::Maybe{CrossValidation}
    push!(predictive_genes, factor_gene_index)
    order_predictive_genes!(; predictive_genes = predictive_genes, candidate_factors = candidate_factors)

    print("\n\e[1A")
    print(
        "\e7$(length(predictive_genes) - 1)" *
        (last_added_gene_index === nothing ? "" : "+") *
        " > $(test_index)" *
        " / $(length(candidate_factors.ordered_factor_genes)) " *
        join(
            [
                (gene_index == last_added_gene_index ? "+" : "") *
                (gene_index == factor_gene_index ? "?+" : "") *
                context.names_of_genes[gene_index] for gene_index in predictive_genes
            ],
            " ",
        ) *
        " ~ $(predictive_cost[1] === nothing ? "NA" : @sprintf("%.5f", predictive_cost[1]))" *
        " RMS $(predictive_cost[2] === nothing ? "NA" : @sprintf("%.5f", predictive_cost[2]))" *
        " R^2 $(predictive_cost[3] === nothing ? "NA" : @sprintf("%.5f", predictive_cost[3]))" *
        " ...\e[J\e8",
    )

    cross_validation = compute_cross_validation(;
        context = context,
        cross_validation_parts = cross_validation_parts,
        predictive_genes = predictive_genes,
        included_genes = included_genes,
        included_metacells = included_metacells,
        core_metacells = core_metacells,
        rng = rng,
    )

    rms_cost_factor = reduce(*, candidate_factors.rms_cost_factor_of_genes[predictive_genes])
    cost = cross_validation.mean_rms * rms_cost_factor

    if predictive_cost[1] === nothing
        improvement = 1.0
    else
        improvement = predictive_cost[1] - cost
    end

    if improvement > 0.0
        predictive_cost[1] = cost
        predictive_cost[2] = cross_validation.mean_rms
        predictive_cost[3] = mean(cross_validation.r2_of_genes[included_genes])
        return cross_validation
    else
        filter!(predictive_genes) do predictive_gene
            return predictive_gene != factor_gene_index
        end
        return nothing
    end
end

function try_remove_factors!(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    candidate_factors::CandidateFactors,
    rng::AbstractRNG,
    predictive_genes::AbstractVector{<:Integer},
    predictive_cost::Vector{Maybe{Float32}},
    last_added_gene_index::Maybe{<:Integer},
)::Tuple{Maybe{<:Integer}, Maybe{CrossValidation}}
    last_removed_gene_index = nothing
    cross_validation = nothing
    predictive_index = 0
    while predictive_index < length(predictive_genes)
        predictive_index += 1
        factor_gene_index = predictive_genes[length(predictive_genes) + 1 - predictive_index]
        if factor_gene_index !== last_added_gene_index
            new_cross_validation = try_remove_factor!(;
                context = context,
                cross_validation_parts = cross_validation_parts,
                included_genes = included_genes,
                included_metacells = included_metacells,
                core_metacells = core_metacells,
                candidate_factors = candidate_factors,
                rng = rng,
                predictive_genes = predictive_genes,
                predictive_cost = predictive_cost,
                last_removed_gene_index = last_removed_gene_index,
                predictive_index = predictive_index,
                factor_gene_index = factor_gene_index,
            )
            if new_cross_validation !== nothing
                last_removed_gene_index = factor_gene_index
                cross_validation = new_cross_validation
            end
        end
    end
    return last_removed_gene_index, cross_validation
end

function try_remove_factor!(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    candidate_factors::CandidateFactors,
    rng::AbstractRNG,
    predictive_genes::AbstractVector{<:Integer},
    predictive_cost::Vector{Maybe{Float32}},
    last_removed_gene_index::Maybe{<:Integer},
    predictive_index::Integer,
    factor_gene_index::Integer,
)::Maybe{CrossValidation}
    print("\n\e[1A")
    print(
        "\e7$(length(predictive_genes))" *
        (last_removed_gene_index === nothing ? "" : "-") *
        " < $(predictive_index)" *
        " / $(length(predictive_genes)) " *
        join(
            [
                (gene_index == last_removed_gene_index ? "-" : "") *
                (gene_index == factor_gene_index ? "?-" : "") *
                context.names_of_genes[gene_index] for gene_index in candidate_factors.ordered_factor_genes if
                gene_index in predictive_genes || gene_index == last_removed_gene_index
            ],
            " ",
        ) *
        " ~ $(@sprintf("%.5f", predictive_cost[1]))" *
        " RMS $(@sprintf("%.5f", predictive_cost[2]))" *
        " R^2 $(@sprintf("%.5f", predictive_cost[3]))" *
        " ...\e[J\e8",
    )

    filter!(predictive_genes) do predictive_gene
        return predictive_gene != factor_gene_index
    end

    cross_validation = compute_cross_validation(;
        context = context,
        cross_validation_parts = cross_validation_parts,
        predictive_genes = predictive_genes,
        included_genes = included_genes,
        included_metacells = included_metacells,
        core_metacells = core_metacells,
        rng = rng,
    )

    rms_cost_factor = reduce(*, candidate_factors.rms_cost_factor_of_genes[predictive_genes])
    cost = cross_validation.mean_rms * rms_cost_factor

    improvement = predictive_cost[1] - cost
    if improvement >= 0
        predictive_cost[1] = cost
        predictive_cost[2] = cross_validation.mean_rms
        predictive_cost[3] = mean(cross_validation.r2_of_genes[included_genes])
        return cross_validation
    else
        push!(predictive_genes, factor_gene_index)
        order_predictive_genes!(; predictive_genes = predictive_genes, candidate_factors = candidate_factors)
        return nothing
    end
end

function order_predictive_genes!(; predictive_genes::Vector{<:Integer}, candidate_factors::CandidateFactors)::Nothing  # untested
    factor_priority_of_predictive_genes = candidate_factors.factor_priority_of_genes[predictive_genes]
    predictive_genes_order = sortperm(factor_priority_of_predictive_genes; rev = true)
    predictive_genes .= predictive_genes[predictive_genes_order]
    return nothing
end

function compute_cross_validation(;  # untested
    context::Context,
    cross_validation_parts::Integer,
    predictive_genes::AbstractVector{<:Integer},
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    core_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    rng::AbstractRNG,
)::CrossValidation
    core_metacells_indices = collect(core_metacells)
    shuffle!(rng, core_metacells_indices)

    floor_chunk_size = max(div(length(core_metacells_indices), cross_validation_parts), 1)
    n_chunks = ceil(length(core_metacells_indices) / floor_chunk_size)
    chunk_size = length(core_metacells_indices) / n_chunks

    included_metacells_mask = zeros(Bool, context.n_metacells)
    included_metacells_mask[included_metacells] .= true

    rms_of_genes = zeros(Float32, context.n_genes)
    r2_of_genes = zeros(Float32, context.n_genes)

    for chunk_index in 1:n_chunks
        first_left_out_position = Int(round((chunk_index - 1) * chunk_size)) + 1
        last_left_out_position = Int(round(chunk_index * chunk_size))
        left_out_positions = first_left_out_position:last_left_out_position

        chunk_metacells_indices = core_metacells_indices[left_out_positions]

        @assert all(included_metacells_mask[chunk_metacells_indices])
        included_metacells_mask[chunk_metacells_indices] .= false
        included_metacells_indices = findall(included_metacells_mask)
        included_metacells_mask[chunk_metacells_indices] .= true

        least_squares = solve_least_squares(;
            context = context,
            predictive_genes = predictive_genes,
            included_genes = included_genes,
            included_metacells = included_metacells_indices,
        )

        evaluate_genes_by_least_squares(;
            context = context,
            least_squares = least_squares,
            chunk_metacells = chunk_metacells_indices,
            rms_of_genes = rms_of_genes,
            r2_of_genes = r2_of_genes,
        )
    end

    rms_of_genes ./= n_chunks
    r2_of_genes ./= n_chunks
    return CrossValidation(;
        mean_rms = mean(rms_of_genes[included_genes]),
        rms_of_genes = rms_of_genes,
        r2_of_genes = r2_of_genes,
    )
end

@kwdef struct LeastSquares
    predictive_genes::AbstractVector{<:Integer}
    included_genes::AbstractVector{<:Integer}
    divergence_of_predictive_genes::AbstractVector{<:AbstractFloat}
    divergence_of_included_genes::AbstractVector{<:AbstractFloat}
    mean_log_fractions_of_predictive_genes::AbstractVector{<:AbstractFloat}
    mean_log_fractions_of_included_genes::AbstractVector{<:AbstractFloat}
    coefficients_of_included_genes_of_predictive_genes::AbstractMatrix{<:AbstractFloat}
end

function solve_least_squares(;  # untested
    context::Context,
    predictive_genes::AbstractVector{<:Integer},
    included_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    included_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
)::LeastSquares
    relative_log_fractions_in_included_metacells_of_included_genes,
    divergence_of_included_genes,
    mean_log_fractions_of_included_genes =
        prepare_data(; context = context, selected_genes = included_genes, selected_metacells = included_metacells)

    relative_log_fractions_in_included_metacells_of_predictive_genes,
    divergence_of_predictive_genes,
    mean_log_fractions_of_predictive_genes =
        prepare_data(; context = context, selected_genes = predictive_genes, selected_metacells = included_metacells)

    coefficients_of_predictive_genes_of_included_genes = nonneg_lsq(
        relative_log_fractions_in_included_metacells_of_predictive_genes,
        relative_log_fractions_in_included_metacells_of_included_genes;
        alg = :fnnls,
    )
    @assert_matrix(
        coefficients_of_predictive_genes_of_included_genes,
        length(predictive_genes),
        length(included_genes),
        Columns
    )

    return LeastSquares(;
        predictive_genes = predictive_genes,
        included_genes = included_genes,
        divergence_of_predictive_genes = divergence_of_predictive_genes,
        divergence_of_included_genes = divergence_of_included_genes,
        mean_log_fractions_of_predictive_genes = mean_log_fractions_of_predictive_genes,
        mean_log_fractions_of_included_genes = mean_log_fractions_of_included_genes,
        coefficients_of_included_genes_of_predictive_genes = transposer(
            coefficients_of_predictive_genes_of_included_genes,
        ),
    )
end

function prepare_data(;  # untested
    context::Context,
    selected_genes::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    selected_metacells::Union{UnitRange{<:Integer}, AbstractVector{<:Integer}},
    mean_log_fractions_of_selected_genes::Maybe{AbstractVector{<:AbstractFloat}} = nothing,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractVector{<:AbstractFloat}, AbstractVector{<:AbstractFloat}}
    log_fractions_of_selected_genes_in_selected_metacells =
        context.log_fractions_of_genes_in_metacells[selected_genes, selected_metacells]
    @assert_matrix(
        log_fractions_of_selected_genes_in_selected_metacells,
        length(selected_genes),
        length(selected_metacells),
        Columns,
    )

    relative_log_fractions_in_selected_metacells_of_selected_genes =
        transpose(log_fractions_of_selected_genes_in_selected_metacells)
    @assert_matrix(
        relative_log_fractions_in_selected_metacells_of_selected_genes,
        length(selected_metacells),
        length(selected_genes),
        Rows,
    )

    if mean_log_fractions_of_selected_genes === nothing
        mean_log_fractions_of_selected_genes =
            vec(mean(relative_log_fractions_in_selected_metacells_of_selected_genes; dims = 1))  # NOJET
    end
    @assert_vector(mean_log_fractions_of_selected_genes, length(selected_genes))

    relative_log_fractions_in_selected_metacells_of_selected_genes .-= transpose(mean_log_fractions_of_selected_genes)

    divergence_of_selected_genes = context.divergence_of_genes[selected_genes]

    relative_log_fractions_in_selected_metacells_of_selected_genes .*= transpose(1.0 .- divergence_of_selected_genes)

    return (
        relative_log_fractions_in_selected_metacells_of_selected_genes,
        divergence_of_selected_genes,
        mean_log_fractions_of_selected_genes,
    )
end

function evaluate_genes_by_least_squares(;  # untested
    context::Context,
    least_squares::LeastSquares,
    chunk_metacells::AbstractVector{<:Integer},
    rms_of_genes::AbstractVector{<:AbstractFloat},
    r2_of_genes::AbstractVector{<:AbstractFloat},
)::Nothing
    n_included_genes = length(least_squares.included_genes)

    divergence_of_included_genes = context.divergence_of_genes[least_squares.included_genes]
    log_fractions_in_chunk_metacells_of_included_genes =
        transposer(context.log_fractions_of_genes_in_metacells[least_squares.included_genes, chunk_metacells])
    predicted_log_fractions_in_chunk_metacells_of_included_genes = predict_by_least_squares(;
        context = context,
        least_squares = least_squares,
        selected_metacells = chunk_metacells,
        divergence_of_included_genes = divergence_of_included_genes,
    )

    residual_log_fractions_in_chunk_metacells_of_included_genes =
        predicted_log_fractions_in_chunk_metacells_of_included_genes .-
        log_fractions_in_chunk_metacells_of_included_genes
    residual_log_fractions_in_chunk_metacells_of_included_genes .*= transpose(1.0 .- divergence_of_included_genes)
    residual_log_fractions_in_chunk_metacells_of_included_genes .*=
        residual_log_fractions_in_chunk_metacells_of_included_genes

    rms_of_included_genes = vec(sqrt.(mean(residual_log_fractions_in_chunk_metacells_of_included_genes; dims = 1)))
    @assert_vector(rms_of_included_genes, n_included_genes)
    rms_of_genes[least_squares.included_genes] .+= rms_of_included_genes

    @threads for included_gene_position in 1:n_included_genes
        included_gene_index = least_squares.included_genes[included_gene_position]
        @views log_fractions_in_chunk_metacells_of_included_gene =
            log_fractions_in_chunk_metacells_of_included_genes[:, included_gene_position]
        @views predicted_log_fractions_in_chunk_metacells_of_included_gene =
            predicted_log_fractions_in_chunk_metacells_of_included_genes[:, included_gene_position]
        correlation_of_gene = vcor(
            log_fractions_in_chunk_metacells_of_included_gene,
            predicted_log_fractions_in_chunk_metacells_of_included_gene,
        )
        if isnan(correlation_of_gene)
            if minimum(log_fractions_in_chunk_metacells_of_included_gene) ==
               maximum(log_fractions_in_chunk_metacells_of_included_gene) &&
               minimum(predicted_log_fractions_in_chunk_metacells_of_included_gene) ==
               maximum(predicted_log_fractions_in_chunk_metacells_of_included_gene)
                correlation_of_gene = 1
            else
                correlation_of_gene = 0
            end
        end
        r2_of_genes[included_gene_index] += correlation_of_gene * correlation_of_gene
    end

    return nothing
end

function predict_by_least_squares(;  # untested
    context::Context,
    least_squares::LeastSquares,
    selected_metacells::AbstractVector{<:Integer},
    divergence_of_included_genes::AbstractVector{<:AbstractFloat},
)::AbstractMatrix{<:AbstractFloat}
    relative_log_fractions_in_selected_metacells_of_predictive_genes, _, _ = prepare_data(;
        context = context,
        selected_metacells = selected_metacells,
        selected_genes = least_squares.predictive_genes,
        mean_log_fractions_of_selected_genes = least_squares.mean_log_fractions_of_predictive_genes,
    )

    predicted_log_fractions_in_selected_metacells_of_included_genes =
        relative_log_fractions_in_selected_metacells_of_predictive_genes *
        transpose(least_squares.coefficients_of_included_genes_of_predictive_genes)
    @assert_matrix(
        predicted_log_fractions_in_selected_metacells_of_included_genes,
        length(selected_metacells),
        length(least_squares.included_genes),
        Columns,
    )

    predicted_log_fractions_in_selected_metacells_of_included_genes ./= transpose(1.0 .- divergence_of_included_genes)
    predicted_log_fractions_in_selected_metacells_of_included_genes .+=
        transpose(least_squares.mean_log_fractions_of_included_genes)

    return predicted_log_fractions_in_selected_metacells_of_included_genes
end

@kwdef struct Confidence
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned}
    log_decreased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
    log_increased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_blocks!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        max_block_span::Real = $(DEFAULT.max_block_span),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given a set of transcription factors that can be used to predict the rest of the genes across the whole manifold, group
the metacells into distinct blocks, such that within each block, the expression level of all these factors differ by a
fold factor of more than `max_block_span` between any of the metacells of the block. When evaluating the difference in
expression between the genes, we reduce the distance using the `fold_confidence` based on the number of UMIs used to
estimate the expression in the metacells; we also ignore any difference between the expression level of genes whose
total UMIs (in both compared metacells) is less thabn `min_significant_gene_UMIs`.

We expect all the metacells within each block to be "essentially identical", so blocks serve as a useful way to group
metacells together for the purpose of type annotations. We also assume that the same local linear program applies for
all metacells of each block. The distance matrix between the blocks expresses the overall structure of the manifold.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        gene_metacell_total_UMIs_matrix(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        block_block_distance_matrix(GuaranteedOutput),
    ],
) function compute_blocks!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_factor_priority_of_genes!,
        :gene_fraction_regularization,
    ),
    fold_confidence::AbstractFloat = 0.9,
    max_block_span::Real = max(function_default(identify_marker_genes!, :min_marker_gene_range_fold) - 1, 1),
    overwrite::Bool = false,
)::Nothing
    @assert min_significant_gene_UMIs >= 0
    @assert gene_fraction_regularization >= 0
    @assert 0 <= fold_confidence <= 1
    @assert max_block_span > 0

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)
    confidence = compute_confidence(;
        context = context,
        gene_fraction_regularization = gene_fraction_regularization,
        fold_confidence = fold_confidence,
    )
    blocks = compute_blocks_by_confidence(;
        context = context,
        confidence = confidence,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        max_block_span = max_block_span,
    )

    block_names = group_names(daf, "metacell", blocks.metacells_of_blocks; prefix = "B")
    add_axis!(daf, "block", block_names)
    set_vector!(daf, "metacell", "block", block_names[blocks.blocks_of_metacells]; overwrite = overwrite)
    distances_between_blocks = blocks.distances_between_blocks
    @assert distances_between_blocks !== nothing
    set_matrix!(daf, "block", "block", "distance", distances_between_blocks; overwrite = overwrite)

    return nothing
end

function compute_confidence(;  # untested
    context::Context,
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
)::Confidence
    total_UMIs_of_metacells = get_vector(context.daf, "metacell", "total_UMIs").array
    total_UMIs_of_genes_in_metacells = get_matrix(context.daf, "gene", "metacell", "total_UMIs").array
    @assert_matrix(total_UMIs_of_genes_in_metacells, context.n_genes, context.n_metacells, Columns)

    log_decreased_fractions_of_genes_in_metacells, log_increased_fractions_of_genes_in_metacells =
        compute_confidence_log_fraction_of_genes_in_metacells(;
            gene_fraction_regularization = gene_fraction_regularization,
            fractions_of_genes_in_metacells = context.fractions_of_genes_in_metacells,
            total_UMIs_of_metacells = total_UMIs_of_metacells,
            fold_confidence = fold_confidence,
        )

    @assert_matrix(log_decreased_fractions_of_genes_in_metacells, context.n_genes, context.n_metacells, Columns,)

    @assert_matrix(log_increased_fractions_of_genes_in_metacells, context.n_genes, context.n_metacells, Columns,)

    return Confidence(;
        total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
        log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_genes_in_metacells,
        log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_genes_in_metacells,
    )
end

function compute_confidence_log_fraction_of_genes_in_metacells(;  # untested
    gene_fraction_regularization::AbstractFloat,
    fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_metacells::AbstractVector{<:Unsigned},
    fold_confidence::AbstractFloat,
)::Tuple{AbstractMatrix{<:AbstractFloat}, AbstractMatrix{<:AbstractFloat}}
    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_fractions_of_genes_in_metacells =  # NOJET
        confidence_stdevs .* sqrt.(transpose(total_UMIs_of_metacells) .* fractions_of_genes_in_metacells) ./
        transpose(total_UMIs_of_metacells)

    log_decreased_fractions_of_genes_in_metacells =
        log2.(
            max.(fractions_of_genes_in_metacells .- confidence_fractions_of_genes_in_metacells, 0.0) .+
            gene_fraction_regularization
        )

    log_increased_fractions_of_genes_in_metacells =
        log2.(
            fractions_of_genes_in_metacells .+ confidence_fractions_of_genes_in_metacells .+
            gene_fraction_regularization
        )

    return (log_decreased_fractions_of_genes_in_metacells, log_increased_fractions_of_genes_in_metacells)
end

@kwdef struct Blocks
    blocks_of_metacells::AbstractVector{<:Integer}
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}}
    n_blocks::Integer
    distances_between_blocks::Maybe{AbstractMatrix{<:AbstractFloat}}
end

function compute_blocks_by_confidence(;  # untested
    context::Context,
    confidence::Confidence,
    min_significant_gene_UMIs::Integer,
    max_block_span::Real,
)::Blocks
    global_predictive_genes = findall(get_vector(context.daf, "gene", "is_global_predictive_factor").array)
    distances_between_metacells = compute_distances_between_metacells(;
        context = context,
        confidence = confidence,
        global_predictive_genes = global_predictive_genes,
        min_significant_gene_UMIs = min_significant_gene_UMIs,
    )

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    blocks_of_metacells = Vector{UInt32}(cutree(clusters; h = max_block_span))
    metacells_of_blocks = collect_group_members(blocks_of_metacells)
    n_blocks = length(metacells_of_blocks)

    distances_between_blocks = compute_distances_between_blocks(;
        distances_between_metacells = distances_between_metacells,
        metacells_of_blocks = metacells_of_blocks,
    )

    return Blocks(;
        blocks_of_metacells = blocks_of_metacells,
        metacells_of_blocks = metacells_of_blocks,
        n_blocks = n_blocks,
        distances_between_blocks = distances_between_blocks,
    )
end

function compute_distances_between_metacells(;  # untested
    context::Context,
    confidence::Confidence,
    global_predictive_genes::AbstractVector{<:Integer},
    min_significant_gene_UMIs::Integer,
)::AbstractMatrix{<:AbstractFloat}
    distances_between_metacells = Matrix{Float32}(undef, context.n_metacells, context.n_metacells)

    divergence_of_global_predictive_genes = context.divergence_of_genes[global_predictive_genes]

    log_decreased_fractions_of_predictive_genes_in_metacells =
        confidence.log_decreased_fractions_of_genes_in_metacells[global_predictive_genes, :]
    log_increased_fractions_of_predictive_genes_in_metacells =
        confidence.log_increased_fractions_of_genes_in_metacells[global_predictive_genes, :]
    total_UMIs_of_predictive_genes_in_metacells =
        confidence.total_UMIs_of_genes_in_metacells[global_predictive_genes, :]

    distances_between_metacells[1, 1] = 0.0
    @threads for base_metacell_index in reverse(2:(context.n_metacells))
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views total_UMIs_of_predictive_genes_in_base_metacell =
            vec(total_UMIs_of_predictive_genes_in_metacells[:, base_metacell_index])
        @views log_decreased_fractions_of_predictive_genes_in_base_metacell =
            vec(log_decreased_fractions_of_predictive_genes_in_metacells[:, base_metacell_index])
        @views log_increased_fractions_of_predictive_genes_in_base_metacell =
            vec(log_increased_fractions_of_predictive_genes_in_metacells[:, base_metacell_index])

        @views total_UMIs_of_predictive_genes_in_other_metacells =
            total_UMIs_of_predictive_genes_in_metacells[:, 1:(base_metacell_index - 1)]
        @views log_decreased_fractions_of_predictive_genes_in_other_metacells =
            log_decreased_fractions_of_predictive_genes_in_metacells[:, 1:(base_metacell_index - 1)]
        @views log_increased_fractions_of_predictive_genes_in_other_metacells =
            log_increased_fractions_of_predictive_genes_in_metacells[:, 1:(base_metacell_index - 1)]

        significant_fold_factors_of_predictive_genes_in_other_metacells =
            gene_distance.(
                min_significant_gene_UMIs,
                total_UMIs_of_predictive_genes_in_base_metacell,
                log_decreased_fractions_of_predictive_genes_in_base_metacell,
                log_increased_fractions_of_predictive_genes_in_base_metacell,
                total_UMIs_of_predictive_genes_in_other_metacells,
                log_decreased_fractions_of_predictive_genes_in_other_metacells,
                log_increased_fractions_of_predictive_genes_in_other_metacells,
                divergence_of_global_predictive_genes,
            )
        @assert_matrix(
            significant_fold_factors_of_predictive_genes_in_other_metacells,
            length(global_predictive_genes),
            base_metacell_index - 1,
            Columns,
        )

        distances_between_base_and_other_metacells =
            vec(maximum(significant_fold_factors_of_predictive_genes_in_other_metacells; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, base_metacell_index - 1)

        distances_between_metacells[1:(base_metacell_index - 1), base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, 1:(base_metacell_index - 1)] .=
            distances_between_base_and_other_metacells
    end

    return distances_between_metacells
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Integer,
    total_UMIs_of_gene_in_base_metacell::Integer,
    log_decreased_fractions_of_gene_in_base_metacell::AbstractFloat,
    log_increased_fractions_of_gene_in_base_metacell::AbstractFloat,
    total_UMIs_of_genes_in_other_metacell::Integer,
    log_decreased_fractions_of_gene_in_other_metacell::AbstractFloat,
    log_increased_fractions_of_gene_in_other_metacell::AbstractFloat,
    divergence::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = total_UMIs_of_gene_in_base_metacell + total_UMIs_of_genes_in_other_metacell
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = log_increased_fractions_of_gene_in_base_metacell < log_increased_fractions_of_gene_in_other_metacell

    log_increased_low_fractions_of_gene =
        is_base_low * log_increased_fractions_of_gene_in_base_metacell +
        !is_base_low * log_increased_fractions_of_gene_in_other_metacell

    log_decreased_high_fractions_of_gene =
        is_base_low * log_decreased_fractions_of_gene_in_other_metacell +
        !is_base_low * log_decreased_fractions_of_gene_in_base_metacell

    return (
        is_significant *
        max.((log_decreased_high_fractions_of_gene - log_increased_low_fractions_of_gene) * (1.0 - divergence), 0.0)
    )
end

function compute_distances_between_blocks(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}},
)::AbstractMatrix{<:AbstractFloat}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacells_of_blocks)
    mean_distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    mean_distances_between_blocks[1, 1] = 0.0
    @threads for base_block_index in reverse(2:n_blocks)
        mean_distances_between_blocks[base_block_index, base_block_index] = 0.0
        metacells_of_base_block = metacells_of_blocks[base_block_index]

        distance_of_metacells_from_base_block_metacells = distances_between_metacells[:, metacells_of_base_block]

        distance_of_metacells_from_base_block = vec(mean(distance_of_metacells_from_base_block_metacells; dims = 2))
        @assert length(distance_of_metacells_from_base_block) == n_metacells

        for other_block_index in 1:(base_block_index - 1)
            metacells_of_other_block = metacells_of_blocks[other_block_index]

            @views distance_of_other_block_metacells_from_base_block =
                distance_of_metacells_from_base_block[metacells_of_other_block]

            distance_between_other_and_base_block = mean(distance_of_other_block_metacells_from_base_block)

            mean_distances_between_blocks[base_block_index, other_block_index] = distance_between_other_and_base_block
            mean_distances_between_blocks[other_block_index, base_block_index] = distance_between_other_and_base_block
        end
    end

    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)
    @threads for base_block_index in 1:n_blocks
        @views mean_distance_between_others_and_base_block = mean_distances_between_blocks[:, base_block_index]
        rank_of_others_for_base_block = invperm(sortperm(mean_distance_between_others_and_base_block))
        @assert rank_of_others_for_base_block[base_block_index] == 1
        distances_between_blocks[:, base_block_index] .= rank_of_others_for_base_block
    end

    distances_between_blocks .*= transpose(distances_between_blocks)
    distances_between_blocks .-= 1
    distances_between_blocks ./= maximum(distances_between_blocks)
    @assert minimum(distances_between_blocks) == 0

    return distances_between_blocks
end

"""
    function compute_blocks_vicinities!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_rms_improvement::AbstractFloat = $(DEFAULT.min_rms_improvement),
        rms_priority_improvement::Real = $(DEFAULT.rms_priority_improvement),
        block_rms_bonus::AbstractFloat = $(DEFAULT.block_rms_bonus),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Given a partition of the metacells into distinct blocks, compute for each block its immediate neighborhood (of at least
`min_blocks_in_neighborhood` blocks and at least `min_metacells_in_neighborhood` metacells). Then expact this to an
environment such that computing a linear model for the environment (based on the precomputed set of global predictive
genes) minimizes the RMS of the error in the neighborhood.

The motivation is that the small neighborhood is large enough to evaluate the linear model, but is too small by itself
for detecting the linear model. Expanding the neighborhood to a larger environment allows us to detect weak linear
programs that exist in the neighborhood, therefore reducing the RMS. Since the manifold is not linear, expanding the
environment beyond a certain point starts to increase the RMS again. We give a slight bonus of `block_rms_bonus` for per
each block added to the environment to reduce the chances of stopping in a local minimum.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_distance_matrix(RequiredInput),
        block_block_is_in_neighborhood_matrix(GuaranteedOutput),
        block_block_is_in_environment_matrix(GuaranteedOutput),
    ],
) function compute_blocks_vicinities!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_factor_priority_of_genes!,
        :gene_fraction_regularization,
    ),
    min_rms_improvement::AbstractFloat = function_default(compute_global_predictive_factors!, :min_rms_improvement),
    rms_priority_improvement::Real = function_default(compute_global_predictive_factors!, :rms_priority_improvement),
    block_rms_bonus::AbstractFloat = 5e-4,
    cross_validation_parts::Integer = function_default(compute_global_predictive_factors!, :cross_validation_parts),
    min_blocks_in_neighborhood::Integer = 4,
    min_metacells_in_neighborhood::Integer = 100,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert cross_validation_parts > 1
    @assert min_rms_improvement >= 0
    @assert 0 <= rms_priority_improvement <= 2
    @assert block_rms_bonus >= 0
    @assert min_blocks_in_neighborhood > 0
    @assert min_metacells_in_neighborhood > 0

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)
    @assert min_metacells_in_neighborhood <= context.n_metacells

    global_predictive_genes = findall(get_vector(daf, "gene", "is_global_predictive_factor").array)

    blocks = load_blocks(daf)
    @assert min_blocks_in_neighborhood <= blocks.n_blocks

    block_block_is_in_neighborhood = zeros(Bool, blocks.n_blocks, blocks.n_blocks)
    block_block_is_in_environment = zeros(Bool, blocks.n_blocks, blocks.n_blocks)

    for block_index in 1:(blocks.n_blocks)
        @views block_is_in_neighborhood = block_block_is_in_neighborhood[:, block_index]
        @views block_is_in_environment = block_block_is_in_environment[:, block_index]
        compute_vicinity_of_block!(;
            context = context,
            cross_validation_parts = cross_validation_parts,
            block_rms_bonus = block_rms_bonus,
            min_blocks_in_neighborhood = min_blocks_in_neighborhood,
            min_metacells_in_neighborhood = min_metacells_in_neighborhood,
            rng = rng,
            blocks = blocks,
            global_predictive_genes = global_predictive_genes,
            block_index = block_index,
            block_is_in_neighborhood = block_is_in_neighborhood,
            block_is_in_environment = block_is_in_environment,
        )
    end

    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_neighborhood",
        SparseMatrixCSC(block_block_is_in_neighborhood);
        overwrite = overwrite,
    )
    set_matrix!(
        daf,
        "block",
        "block",
        "is_in_environment",
        SparseMatrixCSC(block_block_is_in_environment);
        overwrite = overwrite,
    )

    return nothing
end

function load_blocks(daf::DafReader)::Blocks  # untested
    n_blocks = axis_length(daf, "block")
    blocks_of_metacells = axis_indices(daf, "block", get_vector(daf, "metacell", "block"))
    metacells_of_blocks = [findall(blocks_of_metacells .== block_index) for block_index in 1:n_blocks]
    distances_between_blocks = get_matrix(daf, "block", "block", "distance"; default = nothing)

    return Blocks(;
        blocks_of_metacells = blocks_of_metacells,
        metacells_of_blocks = metacells_of_blocks,
        n_blocks = n_blocks,
        distances_between_blocks = distances_between_blocks,
    )
end

function compute_vicinity_of_block!(;  # untested
    context::Context,
    block_rms_bonus::AbstractFloat,
    cross_validation_parts::Integer,
    min_blocks_in_neighborhood::Integer,
    min_metacells_in_neighborhood::Integer,
    global_predictive_genes::AbstractVector{<:Integer},
    rng::AbstractRNG,
    blocks::Blocks,
    block_index::Integer,
    block_is_in_neighborhood::AbstractVector{Bool},
    block_is_in_environment::AbstractVector{Bool},
)::Nothing
    distances_between_blocks = blocks.distances_between_blocks
    @assert distances_between_blocks !== nothing
    distances_between_others_and_block = distances_between_blocks[:, block_index]
    @assert_vector(distances_between_others_and_block, blocks.n_blocks)

    ordered_block_indices = sortperm(distances_between_others_and_block)
    @assert ordered_block_indices[1] == block_index

    n_blocks_in_neighborhood = compute_neighborhood_of_block(;
        blocks = blocks,
        min_blocks_in_neighborhood = min_blocks_in_neighborhood,
        min_metacells_in_neighborhood = min_metacells_in_neighborhood,
        ordered_block_indices = ordered_block_indices,
    )
    block_is_in_neighborhood[ordered_block_indices[1:n_blocks_in_neighborhood]] .= true

    n_blocks_in_environment = compute_environment_of_block(;
        context = context,
        blocks = blocks,
        block_rms_bonus = block_rms_bonus,
        cross_validation_parts = cross_validation_parts,
        global_predictive_genes = global_predictive_genes,
        rng = rng,
        block_index = block_index,
        ordered_block_indices = ordered_block_indices,
        n_blocks_in_neighborhood = n_blocks_in_neighborhood,
    )
    block_is_in_environment[ordered_block_indices[1:n_blocks_in_environment]] .= true

    @debug "$(sum(region_metacells_mask(blocks, ordered_block_indices[1:1]))) mcs" *
           " @ block $(block_index)" *
           " in $(sum(region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood]))) mcs" *
           " @ $(n_blocks_in_neighborhood) neighborhood" *
           " in $(sum(region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_environment]))) mcs" *
           " @ $(n_blocks_in_environment) environment"

    return nothing
end

function compute_neighborhood_of_block(;  # untested
    blocks::Blocks,
    min_blocks_in_neighborhood::Integer,
    min_metacells_in_neighborhood::Integer,
    ordered_block_indices::AbstractVector{<:Integer},
)::Integer
    n_metacells_in_neighborhood = 0
    n_blocks_in_neighborhood = 0
    neighborhood_metacells_mask = nothing

    while n_blocks_in_neighborhood < min_blocks_in_neighborhood ||
        n_metacells_in_neighborhood < min_metacells_in_neighborhood
        n_blocks_in_neighborhood += 1
        neighborhood_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood])
        @assert neighborhood_metacells_mask !== nothing
        n_metacells_in_neighborhood = sum(neighborhood_metacells_mask)
    end

    return n_blocks_in_neighborhood
end

function region_metacells_mask(blocks::Blocks, block_indices::AbstractVector{<:Integer})::Union{Vector{Bool}, BitVector}  # untested
    metacells_mask = zeros(Bool, length(blocks.blocks_of_metacells))
    for block_index in block_indices
        metacells_mask[blocks.metacells_of_blocks[block_index]] .= true
    end
    return metacells_mask
end

function compute_environment_of_block(;  # untested
    context::Context,
    blocks::Blocks,
    block_rms_bonus::AbstractFloat,
    cross_validation_parts::Integer,
    global_predictive_genes::AbstractVector{<:Integer},
    rng::AbstractRNG,
    block_index::Integer,
    ordered_block_indices::AbstractVector{<:Integer},
    n_blocks_in_neighborhood::Integer,
)::Integer
    cost_of_environments = Vector{Maybe{Float32}}(undef, blocks.n_blocks)
    cost_of_environments .= nothing

    n_blocks_of_results = Vector{Int32}()
    cost_of_results = Vector{Float32}()
    n_results = 0

    best_n_blocks = next_n_blocks = n_blocks_in_neighborhood
    while true
        @assert cost_of_environments[next_n_blocks] === nothing

        cost_of_environments[next_n_blocks] =
            cost = cost_of_environment(;
                context = context,
                blocks = blocks,
                block_rms_bonus = block_rms_bonus,
                cross_validation_parts = cross_validation_parts,
                global_predictive_genes = global_predictive_genes,
                rng = rng,
                block_index = block_index,
                ordered_block_indices = ordered_block_indices,
                n_blocks_in_neighborhood = n_blocks_in_neighborhood,
                n_blocks_in_environment = next_n_blocks,
            )

        push!(n_blocks_of_results, next_n_blocks)
        push!(cost_of_results, cost)
        n_results += 1

        order = sortperm(n_blocks_of_results)
        n_blocks_of_results .= n_blocks_of_results[order]
        cost_of_results .= cost_of_results[order]

        best_result_index = argmin(cost_of_results)
        best_n_blocks = n_blocks_of_results[best_result_index]

        if best_result_index == n_results
            high_n_blocks = div(best_n_blocks + blocks.n_blocks + 1, 2)
        else
            high_n_blocks = n_blocks_of_results[best_result_index + 1]
            if high_n_blocks > best_n_blocks + 1
                high_n_blocks = div(best_n_blocks + high_n_blocks, 2)
            end
        end

        if best_result_index == 1
            @assert best_n_blocks == n_blocks_in_neighborhood
            low_n_blocks = best_n_blocks
        else
            low_n_blocks = n_blocks_of_results[best_result_index - 1]
            if low_n_blocks < best_n_blocks - 1
                low_n_blocks = div(low_n_blocks + best_n_blocks, 2)
            end
        end

        high_diff = high_n_blocks - best_n_blocks
        low_diff = best_n_blocks - low_n_blocks

        @assert high_diff >= 0
        @assert low_diff >= 0

        if (
            high_diff <= 1 &&
            low_diff <= 1 &&
            cost_of_environments[high_n_blocks] !== nothing &&
            cost_of_environments[low_n_blocks] !== nothing
        )
            @assert best_n_blocks >= n_blocks_in_neighborhood
            return best_n_blocks

        elseif cost_of_environments[low_n_blocks] !== nothing
            @assert cost_of_environments[high_n_blocks] === nothing
            next_n_blocks = high_n_blocks

        elseif cost_of_environments[high_n_blocks] !== nothing
            @assert cost_of_environments[low_n_blocks] === nothing
            next_n_blocks = low_n_blocks

        elseif high_diff > low_diff
            @assert cost_of_environments[high_n_blocks] === nothing
            next_n_blocks = high_n_blocks

        else
            @assert cost_of_environments[low_n_blocks] === nothing
            next_n_blocks = low_n_blocks
        end
    end

    @assert false
end

function cost_of_environment(;  # untested
    context::Context,
    blocks::Blocks,
    block_rms_bonus::AbstractFloat,
    cross_validation_parts::Integer,
    global_predictive_genes::AbstractVector{<:Integer},
    rng::AbstractRNG,
    block_index::Integer,
    ordered_block_indices::AbstractVector{<:Integer},
    n_blocks_in_neighborhood::Integer,
    n_blocks_in_environment::Integer,
)::AbstractFloat
    @assert 0 < n_blocks_in_neighborhood <= n_blocks_in_environment

    neighborhoods_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood])
    neighborhoods_metacells = findall(neighborhoods_metacells_mask)

    environment_metacells_mask = region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_environment])
    environment_metacells = findall(environment_metacells_mask)

    cross_validation = compute_cross_validation(;
        context = context,
        cross_validation_parts = cross_validation_parts,
        predictive_genes = global_predictive_genes,
        included_genes = 1:(context.n_genes),
        included_metacells = environment_metacells,
        core_metacells = neighborhoods_metacells,
        rng = rng,
    )

    cost =
        cross_validation.mean_rms *
        (1 - block_rms_bonus * sum(1 ./ (1:(n_blocks_in_environment - n_blocks_in_neighborhood))))

    print("\n\e[1A")
    print(
        "\e7$(sum(region_metacells_mask(blocks, ordered_block_indices[1:1]))) mcs" *
        " @ block $(block_index)" *
        " in $(sum(region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_neighborhood]))) mcs" *
        " @ $(n_blocks_in_neighborhood) neighborhood" *
        " in $(sum(region_metacells_mask(blocks, ordered_block_indices[1:n_blocks_in_environment]))) mcs" *
        " @ $(n_blocks_in_environment) environment" *
        " : $(@sprintf("%.5f", cost))" *
        " RMS $(@sprintf("%.5f", cross_validation.mean_rms))" *
        " R^2 $(@sprintf("%.5f", mean(cross_validation.r2_of_genes)))" *
        " ...\e[J\e8",
    )

    return cost
end

"""
    function compute_local_predictive_factors!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_rms_improvement::AbstractFloat = $(DEFAULT.min_rms_improvement),
        rms_priority_improvement::Real = $(DEFAULT.rms_priority_improvement),
        cross_validation_parts::Integer = $(DEFAULT.cross_validation_parts),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Having computed the neighborhoods and environments, then for each block, figure out the set of transcription factors for
best approximating the gene expression of the metacells in each neighborhood based on the environment. This set of local
predictive factors will be different from the set of global predictive factors (and will be different for each block).

When computing this set, we only consider the genes that are marker genes within the environment (as per
[`identify_marker_genes!`](@ref), using `min_marker_gene_max_fraction` and a tighter `min_marker_gene_range_fold`.

$(CONTRACT)
"""
@logged @computation Contract(
    is_relaxed = true,
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_factor_priority_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_neighborhood_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        gene_block_is_local_marker_matrix(GuaranteedOutput),
        gene_block_is_local_predictive_factor_matrix(GuaranteedOutput),
        gene_block_local_rms_matrix(GuaranteedOutput),
        gene_block_local_r2_matrix(GuaranteedOutput),
    ],
) function compute_local_predictive_factors!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_factor_priority_of_genes!,
        :gene_fraction_regularization,
    ),
    min_rms_improvement::AbstractFloat = function_default(compute_global_predictive_factors!, :min_rms_improvement),
    rms_priority_improvement::Real = function_default(compute_global_predictive_factors!, :rms_priority_improvement),
    cross_validation_parts::Integer = function_default(compute_global_predictive_factors!, :cross_validation_parts),
    min_marker_gene_range_fold::Real = max(
        function_default(identify_marker_genes!, :min_marker_gene_range_fold) - 1,
        1,
    ),
    min_marker_gene_max_fraction::AbstractFloat = function_default(
        identify_marker_genes!,
        :min_marker_gene_max_fraction,
    ),
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_rms_improvement >= 0
    @assert cross_validation_parts > 0
    @assert min_marker_gene_range_fold >= 0
    @assert min_marker_gene_max_fraction >= 0

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)

    candidate_factors = load_candidate_factors(
        daf;
        min_rms_improvement = min_rms_improvement,
        rms_priority_improvement = rms_priority_improvement,
    )

    global_predictive_genes = findall(get_vector(daf, "gene", "is_global_predictive_factor").array)

    blocks = load_blocks(daf)
    block_block_is_in_neighborhood = get_matrix(daf, "block", "block", "is_in_neighborhood")
    block_block_is_in_environment = get_matrix(daf, "block", "block", "is_in_environment")

    gene_block_is_local_predictive_factor = zeros(Bool, context.n_genes, blocks.n_blocks)
    gene_block_is_local_marker = zeros(Bool, context.n_genes, blocks.n_blocks)
    gene_block_local_rms = zeros(Float32, context.n_genes, blocks.n_blocks)
    gene_block_local_r2 = zeros(Float32, context.n_genes, blocks.n_blocks)

    for block_index in 1:(blocks.n_blocks)
        blocks_in_neighborhood = findall(block_block_is_in_neighborhood[:, block_index])
        blocks_in_environment = findall(block_block_is_in_environment[:, block_index])
        @views block_local_predictive_factors_mask = gene_block_is_local_predictive_factor[:, block_index]
        @views block_local_markers_mask = gene_block_is_local_marker[:, block_index]
        @views block_local_rms_of_genes = gene_block_local_rms[:, block_index]
        @views block_local_r2_of_genes = gene_block_local_r2[:, block_index]
        compute_local_predictive_factors_of_block(;
            context = context,
            global_predictive_genes = global_predictive_genes,
            candidate_factors = candidate_factors,
            blocks = blocks,
            gene_fraction_regularization = gene_fraction_regularization,
            cross_validation_parts = cross_validation_parts,
            min_marker_gene_range_fold = min_marker_gene_range_fold,
            min_marker_gene_max_fraction = min_marker_gene_max_fraction,
            rng = rng,
            block_index = block_index,
            blocks_in_neighborhood = blocks_in_neighborhood,
            blocks_in_environment = blocks_in_environment,
            block_local_predictive_factors_mask = block_local_predictive_factors_mask,
            block_local_markers_mask = block_local_markers_mask,
            block_local_rms_of_genes = block_local_rms_of_genes,
            block_local_r2_of_genes = block_local_r2_of_genes,
            overwrite = overwrite,
        )
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "is_local_predictive_factor",
        SparseMatrixCSC(gene_block_is_local_predictive_factor);
        overwrite = overwrite,
    )

    set_matrix!(
        daf,
        "gene",
        "block",
        "is_local_marker",
        SparseMatrixCSC(gene_block_is_local_marker);
        overwrite = overwrite,
    )

    set_matrix!(daf, "gene", "block", "local_rms", SparseMatrixCSC(gene_block_local_rms); overwrite = overwrite)

    set_matrix!(daf, "gene", "block", "local_r2", SparseMatrixCSC(gene_block_local_r2); overwrite = overwrite)

    return nothing
end

function compute_local_predictive_factors_of_block(;  # untested
    context::Context,
    global_predictive_genes::AbstractVector{<:Integer},
    candidate_factors::CandidateFactors,
    blocks::Blocks,
    gene_fraction_regularization::AbstractFloat,
    cross_validation_parts::Integer,
    min_marker_gene_range_fold::Real,
    min_marker_gene_max_fraction::AbstractFloat,
    rng::AbstractRNG,
    block_index::Integer,
    blocks_in_neighborhood::AbstractVector{<:Integer},
    blocks_in_environment::AbstractVector{<:Integer},
    block_local_predictive_factors_mask::Union{AbstractVector{Bool}, BitVector},
    block_local_markers_mask::Union{AbstractVector{Bool}, BitVector},
    block_local_rms_of_genes::AbstractVector{<:AbstractFloat},
    block_local_r2_of_genes::AbstractVector{<:AbstractFloat},
    overwrite::Bool,
)::Nothing
    neighborhood_metacells_mask = region_metacells_mask(blocks, blocks_in_neighborhood)
    neighborhood_metacells = findall(neighborhood_metacells_mask)

    environment_metacells_mask = region_metacells_mask(blocks, blocks_in_environment)
    environment_metacells = findall(environment_metacells_mask)

    environment_genes_mask = copy_array(
        marker_genes_of_environment(;
            context = context,
            gene_fraction_regularization = gene_fraction_regularization,
            min_marker_gene_range_fold = min_marker_gene_range_fold,
            min_marker_gene_max_fraction = min_marker_gene_max_fraction,
            environment_metacells_mask = environment_metacells_mask,
            overwrite = overwrite,
        ),
    )
    block_local_markers_mask .= environment_genes_mask
    environment_genes_mask[global_predictive_genes] .= true
    environment_genes = findall(environment_genes_mask)

    local_predictive_genes, cross_validation = identify_predictive_genes!(;
        context = context,
        cross_validation_parts = cross_validation_parts,
        included_genes = environment_genes,
        included_metacells = environment_metacells,
        core_metacells = neighborhood_metacells,
        candidate_factors = candidate_factors,
        rng = rng,
    )

    @debug "Block $(block_index) Predictive : [ $(join(context.names_of_genes[local_predictive_genes], ", ")) ]"
    @debug "Block $(block_index) Markers : $(length(environment_genes))"
    @debug "Block $(block_index) RMS : $(mean(cross_validation.rms_of_genes[environment_genes])) R^2 : $(mean(cross_validation.r2_of_genes[environment_genes]))"

    block_local_predictive_factors_mask[local_predictive_genes] .= true
    block_local_rms_of_genes .= cross_validation.rms_of_genes
    block_local_r2_of_genes .= cross_validation.r2_of_genes

    reused = 0
    added = 0
    removed = 0

    for predictive_gene_index in global_predictive_genes
        if predictive_gene_index in local_predictive_genes
            reused += 1
        else
            removed += 1
        end
    end

    for predictive_gene_index in local_predictive_genes
        if !(predictive_gene_index in global_predictive_genes)
            added += 1
        end
    end

    @debug "- markers: $(length(environment_genes)) local: $(length(local_predictive_genes)) reused: $(reused) added: $(added) removed: $(removed)"

    return nothing
end

function marker_genes_of_environment(;  # untested
    context::Context,
    gene_fraction_regularization::AbstractFloat,
    min_marker_gene_range_fold::Real,
    min_marker_gene_max_fraction::AbstractFloat,
    environment_metacells_mask::Union{AbstractVector{Bool}, BitVector},
    overwrite::Bool,
)::Union{AbstractVector{Bool}, BitVector}
    chain = chain_writer([context.daf, MemoryDaf(; name = "environment")]; name = "mask_chain")
    set_vector!(chain, "metacell", "is_in_environment", environment_metacells_mask; overwrite = overwrite)
    adapter(  # NOJET
        chain;
        input_axes = ["metacell" => "/metacell & is_in_environment", "gene" => "="],
        input_data = [("gene", "divergence") => "=", ("metacell", "gene", "fraction") => "="],
        output_axes = ["gene" => "="],
        output_data = [("gene", "is_marker") => "="],
        overwrite = true,
    ) do adapted
        return identify_marker_genes!(
            adapted;
            gene_fraction_regularization = gene_fraction_regularization,
            min_marker_gene_range_fold = min_marker_gene_range_fold,
            min_marker_gene_max_fraction = min_marker_gene_max_fraction,
            overwrite = true,
        )
    end
    return get_vector(chain, "gene", "is_marker").array
end

"""
    function compute_programs!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Having computed, for each block, the set of transcription factors for best approximating the gene expression of the
metacells in each neighborhood based on the environment, then compute the actual coefficients for doing this prediction.
If `overwrite`, will overwrite existing data.

The program for predicting the value of a gene `i` based on the values of the predictive transcription factors `j` is
expressed as `G_i = MG_i + Sum P_ij (F_j - MF_j)`, where `MG_i` is the mean of the (log base 2 of the) expression of the
gene `i` in the environment, and similarly `MF_j` is the mean of the (log base 2 of the) expression of the factor `j` in
the same environment; and the `P_ij` is the coefficient for computing the (log base 2 of the) expression of gene `i`
using the (log base 2 of the) expression of the factor `F_j`.

We store the `gene_fraction_regularization` used to compute the log base 2 of the fractions as a scalar called
`program_gene_fraction_regularization`, and the means as the per-block-per-gene matrix `program_mean_log_fraction`. The
coefficients themselves are stored as a set of matrices `<block_name>_program_coefficient` which hold the `P_ij` matrix
(columns for transcription factors, rows for genes).

$(CONTRACT)

**gene, gene @ <block_name>_program_coefficient**::AbstractFloat (guaranteed): The (sparse) coefficient for predicting
the (log base 2) expression of each (row) gene based on the (log base 2) expression of each (column) transcription
factor.
"""
@logged @computation Contract(
    is_relaxed = true,
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        metacell_block_vector(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        gene_block_is_local_predictive_factor_matrix(RequiredInput),
        program_gene_fraction_regularization_scalar(GuaranteedOutput),
        gene_block_program_mean_log_fraction_matrix(GuaranteedOutput),
    ],
) function compute_programs!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = function_default(
        compute_factor_priority_of_genes!,
        :gene_fraction_regularization,
    ),
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0

    context = load_context(daf; gene_fraction_regularization = gene_fraction_regularization)

    blocks = load_blocks(daf)
    block_block_is_in_environment = get_matrix(daf, "block", "block", "is_in_environment").array
    gene_block_is_local_predictive_factor = get_matrix(daf, "gene", "block", "is_local_predictive_factor").array
    names_of_blocks = axis_array(daf, "block")

    set_scalar!(daf, "program_gene_fraction_regularization", gene_fraction_regularization; overwrite = overwrite)

    n_genes = axis_length(daf, "gene")
    n_blocks = axis_length(daf, "block")
    mean_log_fraction_of_genes_in_blocks = zeros(Float32, n_genes, n_blocks)

    for block_index in 1:(blocks.n_blocks)
        block_name = names_of_blocks[block_index]
        blocks_in_environment = findall(block_block_is_in_environment[:, block_index])
        environment_metacells_mask = region_metacells_mask(blocks, blocks_in_environment)

        @debug "Block $(block_name) ..."

        @views block_local_predictive_factors_mask = gene_block_is_local_predictive_factor[:, block_index]
        predictive_genes = findall(block_local_predictive_factors_mask)
        least_squares = solve_least_squares(;
            context = context,
            predictive_genes = predictive_genes,
            included_genes = 1:n_genes,
            included_metacells = findall(environment_metacells_mask),
        )

        mean_log_fraction_of_genes_in_blocks[:, block_index] .= least_squares.mean_log_fractions_of_included_genes

        coefficients_of_included_genes_of_predictive_genes =
            SparseMatrixCSC(least_squares.coefficients_of_included_genes_of_predictive_genes)

        empty_sparse_matrix!(
            daf,
            "gene",
            "gene",
            "$(block_name)_program_coefficient",
            eltype(coefficients_of_included_genes_of_predictive_genes),
            nnz(coefficients_of_included_genes_of_predictive_genes);
            overwrite = overwrite,
        ) do empty_colptr, empty_rowval, empty_nzval
            empty_rowval .= coefficients_of_included_genes_of_predictive_genes.rowval
            empty_nzval .= coefficients_of_included_genes_of_predictive_genes.nzval

            next_predictive_gene_position = 1
            next_predictive_gene_index = predictive_genes[1]
            empty_colptr[1] = 1
            for gene_index in 1:n_genes
                if gene_index < next_predictive_gene_index
                    empty_colptr[gene_index + 1] = empty_colptr[gene_index]
                else
                    @assert gene_index == next_predictive_gene_index
                    empty_colptr[gene_index + 1] =
                        coefficients_of_included_genes_of_predictive_genes.colptr[next_predictive_gene_position + 1]
                    next_predictive_gene_position += 1
                    if next_predictive_gene_position <= length(predictive_genes)
                        next_predictive_gene_index = predictive_genes[next_predictive_gene_position]
                    else
                        next_predictive_gene_index = n_genes + 1
                    end
                end
            end
            @assert next_predictive_gene_position == length(predictive_genes) + 1
            @assert next_predictive_gene_index == n_genes + 1
            @assert empty_colptr[n_genes + 1] == length(empty_nzval) + 1
        end
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "program_mean_log_fraction",
        mean_log_fraction_of_genes_in_blocks;
        overwrite = overwrite,
    )

    return nothing
end

end  # module
