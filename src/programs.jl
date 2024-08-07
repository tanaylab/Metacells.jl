"""
Given a set of raw metacells, partition them into blocks such that all metacells in the same block are within some
(fold factor) radius of each other. The centroids of these blocks can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states. Group these blocks in overlapping neighborhoods
of "similar" blocks for further analysis.
"""
module Programs

export compute_global_predictive_factors!
export compute_local_predictive_factors!

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

@kwdef struct ParametersContext
    max_principal_components::Integer
    min_significant_gene_UMIs::Integer = -1
    fold_confidence::AbstractFloat = -1
    max_block_span::AbstractFloat = -1
    min_blocks_in_neighborhood::Integer = -1
    min_metacells_in_neighborhood::Integer = -1
    gene_fraction_regularization::AbstractFloat
    rng::AbstractRNG
end

@kwdef struct SomeIndices
    indices::AbstractVector{<:Integer}
    n_entries::Integer
end

function Daf.depict(some_indices::SomeIndices)::String
    return "$(some_indices.n_entries)"
end

@kwdef struct SomeMask
    mask::Union{AbstractVector{Bool}, BitVector}
    n_entries::Integer
end

function Daf.depict(some_mask::SomeMask)::String
    return "$(some_mask.n_entries)"
end

@kwdef struct DataContext
    n_genes::Integer
    n_metacells::Integer
    names_of_genes::AbstractVector{<:AbstractString}
    factor_genes::SomeIndices
    divergence_of_genes::AbstractVector{<:AbstractFloat}
    fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
    log_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
end

@kwdef struct Confidence
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned}
    log_decreased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
    log_increased_fractions_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat}
end

@kwdef struct LeftOutBuffers
    log_fractions_of_genes_in_left_out_metacells::AbstractMatrix{<:AbstractFloat}
    predicted_log_fractions_of_genes_in_left_out_metacells::AbstractMatrix{<:AbstractFloat}
end

@kwdef struct Variability
    mean_rms::AbstractFloat
end

function Daf.depict(variability::Variability)::String
    return "{ RMS: $(variability.mean_rms) }"
end

@kwdef struct CrossValidation
    mean_r2::AbstractFloat
    mean_shuffled_r2::AbstractFloat
    mean_extra_r2::AbstractFloat
    mean_rms::AbstractFloat
end

function Daf.depict(cross_validation::CrossValidation)::String
    return "{ RMS: $(cross_validation.mean_rms), R2: $(cross_validation.mean_extra_r2) }"
end

@kwdef struct LeastSquares
    mean_log_fractions_of_genes::AbstractVector{<:AbstractFloat}
    coefficients_of_genes_of_predictive_genes::AbstractMatrix{<:AbstractFloat}
end

@kwdef struct Blocks
    blocks_of_metacells::AbstractVector{<:Integer}
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}}
    n_blocks::Integer
    distances_between_blocks::AbstractMatrix{<:AbstractFloat}
end

function Daf.depict(blocks::Blocks)::String
    return "$(blocks.n_blocks)"
end

@kwdef struct Analysis
    core_metacells::SomeMask
    included_metacells::SomeMask
    predictive_genes::SomeIndices
    least_squares::LeastSquares
    variability::Variability
    cross_validation::CrossValidation
end

function Daf.depict(analysis::Analysis)::String
    return "{ variability: $(depict(analysis.variability)), cross_validation: $(depict(analysis.cross_validation)) }"
end

@kwdef struct Region
    blocks::SomeIndices
    metacells::SomeMask
end

function Daf.depict(region::Region)::String
    return "$(region.metacells.n_entries) metacells in $(region.blocks.n_entries) blocks"
end

@kwdef struct Environment
    region::Region
    neighborhood::Region
    local_analysis::Analysis
    by_global_analysis::Analysis
    global_analysis::Analysis
end

function Daf.depict(environment::Environment)::String
    return "$(depict(environment.neighborhood)) inside $(depict(environment.region))"
end

@kwdef mutable struct Context
    parameters::ParametersContext
    data::DataContext
    confidence::Maybe{Confidence} = nothing
    ordered_factor_gene_indices::Maybe{AbstractVector{<:Integer}} = nothing
    minimal_rank_of_genes::Maybe{AbstractVector{<:AbstractFloat}} = nothing
    core_metacells::Maybe{SomeMask} = nothing
    included_metacells::Maybe{SomeMask} = nothing
    predictive_genes::Maybe{SomeIndices} = nothing
    left_out_metacells::Maybe{SomeIndices} = nothing
    chunk_metacells::Maybe{SomeIndices} = nothing
    least_squares::Maybe{LeastSquares} = nothing
    left_out_buffers::Maybe{LeftOutBuffers} = nothing
    cross_validation::Maybe{CrossValidation} = nothing
    included_variability::Maybe{Variability} = nothing
    blocks::Maybe{Blocks} = nothing
    environments::Maybe{AbstractVector{Environment}} = nothing
end

function Daf.depict(context::Context)::String
    fields = String[]
    if context.core_metacells !== nothing
        push!(fields, "core_metacells: $(depict(context.core_metacells))")
    end
    if context.included_metacells !== nothing
        push!(fields, "included_metacells: $(depict(context.included_metacells))")
    end
    if context.predictive_genes !== nothing
        push!(fields, "predictive_genes: $(depict(context.predictive_genes))")
    end
    if context.left_out_metacells !== nothing
        push!(fields, "left_out_metacells: $(depict(context.left_out_metacells))")
    end
    if context.chunk_metacells !== nothing
        push!(fields, "chunk_metacells: $(depict(context.chunk_metacells))")
    end
    if context.least_squares !== nothing
        push!(fields, depict(context.least_squares))
    end
    if context.left_out_buffers !== nothing
        push!(fields, depict(context.left_out_buffers))
    end
    if context.cross_validation !== nothing
        push!(fields, "cross_validation: $(depict(context.cross_validation))")
    end
    if context.included_variability !== nothing
        push!(fields, "included_variability: $(depict(context.included_variability))")
    end
    if context.blocks !== nothing
        push!(fields, "blocks: $(depict(context.blocks))")
    end
    return join(["Context: {", join(fields, ", "), "}"], " ")
end

@logged function DataContext(daf::DafReader, parameters::ParametersContext)::DataContext
    n_genes = axis_length(daf, "gene")
    @assert n_genes > 0

    n_metacells = axis_length(daf, "metacell")
    @assert n_metacells > 0

    names_of_genes = axis_array(daf, "gene")

    factor_genes_mask = get_vector(daf, "gene", "is_transcription_factor")
    factor_genes_indices = findall(factor_genes_mask)
    factor_genes = SomeIndices(factor_genes_indices, length(factor_genes_indices))
    @assert factor_genes.n_entries > 0

    divergence_of_genes = get_vector(daf, "gene", "divergence").array

    fractions_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    @assert_matrix(fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    log_fractions_of_genes_in_metacells =
        log2.(fractions_of_genes_in_metacells .+ parameters.gene_fraction_regularization)
    @assert_matrix(log_fractions_of_genes_in_metacells, n_genes, n_metacells, Columns)

    return DataContext(;
        n_genes = n_genes,
        n_metacells = n_metacells,
        names_of_genes = names_of_genes,
        factor_genes = factor_genes,
        divergence_of_genes = divergence_of_genes,
        fractions_of_genes_in_metacells = fractions_of_genes_in_metacells,
        log_fractions_of_genes_in_metacells = log_fractions_of_genes_in_metacells,
    )
end

function reset_core!(
    context::Context,
    core_metacells_mask::Union{AbstractVector{Bool}, BitVector};
    keep_predictive::Bool = false,
)::SomeMask
    @assert_vector(core_metacells_mask, context.data.n_metacells)

    context.core_metacells = core_metacells = SomeMask(core_metacells_mask, sum(core_metacells_mask))
    context.included_metacells = nothing

    if keep_predictive
        clear_left_out_buffers!(context)
    else
        clear_predictive!(context)
    end

    context.cross_validation = nothing
    return core_metacells
end

function reset_included!(
    context::Context,
    included_metacells_mask::Union{AbstractVector{Bool}, BitVector};
    keep_predictive::Bool = false,
)::SomeMask
    @assert_vector(included_metacells_mask, context.data.n_metacells)

    context.included_metacells = included_metacells = SomeMask(included_metacells_mask, sum(included_metacells_mask))

    if keep_predictive
        clear_left_out_buffers!(context)
    else
        clear_predictive!(context)
    end

    context.cross_validation = nothing
    return included_metacells
end

function clear_predictive!(context::Context)::Nothing
    context.predictive_genes = nothing
    clear_left_out_buffers!(context)
    return nothing
end

function reset_predictive!(context::Context, predictive_genes_indices::AbstractVector{<:Integer})::SomeIndices
    context.predictive_genes =
        predictive_genes = SomeIndices(predictive_genes_indices, length(predictive_genes_indices))

    clear_left_out_buffers!(context)

    return predictive_genes
end

function push_predictive!(context::Context, factor_gene_index::Integer)::SomeIndices
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    push!(predictive_genes.indices, factor_gene_index)
    context.predictive_genes =
        predictive_genes = SomeIndices(predictive_genes.indices, length(predictive_genes.indices))

    clear_left_out_buffers!(context)

    return predictive_genes
end

function pop_predictive!(context::Context, factor_gene_index::Integer)::SomeIndices
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    popped_factor_index = pop!(predictive_genes.indices)
    @assert popped_factor_index == factor_gene_index
    context.predictive_genes =
        predictive_genes = SomeIndices(predictive_genes.indices, length(predictive_genes.indices))

    clear_left_out_buffers!(context)

    return predictive_genes
end

function remove_predictive!(context::Context, factor_gene_index::Integer)::SomeIndices
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    filtered_gene_indices = [
        current_factor_index for
        current_factor_index in predictive_genes.indices if current_factor_index != factor_gene_index
    ]
    @assert length(filtered_gene_indices) == predictive_genes.n_entries - 1
    context.predictive_genes = predictive_genes = SomeIndices(filtered_gene_indices, predictive_genes.n_entries - 1)

    clear_left_out_buffers!(context)

    return predictive_genes
end

function clear_left_out_buffers!(context::Context)::Nothing
    context.left_out_metacells = nothing
    context.left_out_buffers = nothing

    clear_chunk!(context)

    return nothing
end

function reset_left_out_buffers!(
    context::Context,
    left_out_metacells_indices::AbstractVector{<:Integer},
)::Tuple{SomeIndices, LeftOutBuffers}
    @assert context.included_metacells !== nothing
    @assert context.predictive_genes !== nothing

    context.left_out_metacells =
        left_out_metacells =
            SomeIndices(; indices = left_out_metacells_indices, n_entries = length(left_out_metacells_indices))

    context.left_out_buffers =
        left_out_buffers = LeftOutBuffers(;
            log_fractions_of_genes_in_left_out_metacells = Matrix{Float32}(
                undef,
                context.data.n_genes,
                left_out_metacells.n_entries,
            ),
            predicted_log_fractions_of_genes_in_left_out_metacells = Matrix{Float32}(
                undef,
                context.data.n_genes,
                left_out_metacells.n_entries,
            ),
        )

    clear_chunk!(context)

    return left_out_metacells, left_out_buffers
end

function clear_chunk!(context::Context)::Nothing
    context.chunk_metacells = nothing
    return clear_least_squares!(context)
end

function reset_chunk!(context::Context, chunk_metacells_indices::AbstractVector{<:Integer})::SomeIndices
    @assert context.included_metacells !== nothing
    @assert context.predictive_genes !== nothing
    @assert context.left_out_metacells !== nothing

    context.chunk_metacells =
        chunk_metacells = SomeIndices(; indices = chunk_metacells_indices, n_entries = length(chunk_metacells_indices))

    clear_least_squares!(context)

    return chunk_metacells
end

function clear_least_squares!(context::Context)::Nothing
    context.least_squares = nothing
    return nothing
end

function clear_cross_validation!(context::Context)::Nothing
    context.cross_validation = nothing
    return nothing
end

"""
    function compute_global_predictive_factors!(  # untested
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization)
        max_principal_components::Integer = $(DEFAULT.max_principal_components),
        rng::AbstractRNG = default_rng(),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_transcription_factor_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(GuaranteedOutput),
    ],
) function compute_global_predictive_factors!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    max_principal_components::Integer = 30,
    rng::AbstractRNG = default_rng(),
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert max_principal_components > 0

    parameters = ParametersContext(;
        max_principal_components = max_principal_components,
        gene_fraction_regularization = gene_fraction_regularization,
        rng = rng,
    )
    context = Context(; parameters = parameters, data = DataContext(daf, parameters))

    order_factors!(context)
    reset_included!(context, ones(Bool, context.data.n_metacells))
    reset_predictive!(
        context,
        context.ordered_factor_gene_indices[1:div(length(context.ordered_factor_gene_indices), 4)],
    )
    global_analysis = analyze_without_overfitting!(context)

    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    @debug "global predictive factors: $(join(context.data.names_of_genes[predictive_genes.indices], ", "))"
    @debug "global rms: $(global_analysis.variability.mean_rms)"
    @debug "left out rms: $(global_analysis.cross_validation.mean_rms)"
    @debug "left out extra r2: $(global_analysis.cross_validation.mean_extra_r2)"
    @debug "left out shuffled r2: $(global_analysis.cross_validation.mean_shuffled_r2)"
    @debug "left out r2: $(global_analysis.cross_validation.mean_r2)"

    predictive_genes_mask = zeros(Bool, context.data.n_genes)
    predictive_genes_mask[predictive_genes.indices] .= true
    set_vector!(daf, "gene", "is_global_predictive_factor", predictive_genes_mask)
    return nothing
end

@logged function order_factors!(context::Context)::Nothing
    @assert context.ordered_factor_gene_indices === nothing

    log_fractions_of_genes_in_metacells = context.data.log_fractions_of_genes_in_metacells

    pca = fit(PCA, log_fractions_of_genes_in_metacells; maxoutdim = context.parameters.max_principal_components)

    n_principal_components = size(pca, 2)

    coefficients_of_genes_of_principal_components = loadings(pca)
    @assert_matrix(coefficients_of_genes_of_principal_components, context.data.n_genes, n_principal_components, Columns)

    coefficients_of_factor_genes_of_principal_components =
        coefficients_of_genes_of_principal_components[context.data.factor_genes.indices, :]
    @assert_matrix(
        coefficients_of_factor_genes_of_principal_components,
        context.data.factor_genes.n_entries,
        n_principal_components,
        Columns
    )

    minimal_rank_of_factor_genes = fill(context.data.factor_genes.n_entries * 3.0, context.data.factor_genes.n_entries)

    for principal_component_index in 1:n_principal_components
        abs_coefficients_of_factors_in_principal_component =
            abs.(coefficients_of_factor_genes_of_principal_components[:, principal_component_index])

        factor_genes_order = sortperm(abs_coefficients_of_factors_in_principal_component; rev = true)
        rank_of_factor_genes = invperm(factor_genes_order) .* (1 + principal_component_index / n_principal_components)
        minimal_rank_of_factor_genes .= min.(minimal_rank_of_factor_genes, rank_of_factor_genes)
    end

    factor_genes_order = sortperm(minimal_rank_of_factor_genes)
    ordered_factor_gene_indices = context.data.factor_genes.indices[factor_genes_order[1:(n_principal_components * 2)]]
    @debug "top $(length(ordered_factor_gene_indices)) ranked transcription factors: $(join(context.data.names_of_genes[ordered_factor_gene_indices], ", "))"
    open("ordered_factors.csv", "w") do file
        println(file, "gene_index,gene_name")
        for gene_index in ordered_factor_gene_indices
            println(file, "$(gene_index),$(context.data.names_of_genes[gene_index])")
        end
    end

    minimal_rank_of_genes = fill(context.data.factor_genes.n_entries * 3.0, context.data.n_genes)
    minimal_rank_of_genes[context.data.factor_genes.indices] .= minimal_rank_of_factor_genes
    context.minimal_rank_of_genes = minimal_rank_of_genes
    context.ordered_factor_gene_indices = ordered_factor_gene_indices

    return nothing
end

function analyze_without_overfitting!(context::Context)::Analysis
    @assert context.predictive_genes !== nothing
    analyze_left_outs!(context)

    op_index = 1
    op_success = [true, true]
    last_added_gene_index = nothing
    last_removed_gene_index = nothing
    while sum(op_success) > 0
        op_index = 1 + op_index % 2
        if op_index == 1
            last_added_gene_index = try_add_factors!(context, last_removed_gene_index)
            op_success[op_index] = last_added_gene_index !== nothing
        else
            last_removed_gene_index = try_remove_factors!(context, last_added_gene_index)
            op_success[op_index] = last_removed_gene_index !== nothing
        end
    end

    analysis = final_analysis!(context)

    clear_left_out_buffers!(context)

    return analysis
end

function try_add_factors!(context::Context, last_removed_gene_index::Maybe{<:Integer})::Maybe{<:Integer}
    ordered_factor_gene_indices = context.ordered_factor_gene_indices
    @assert ordered_factor_gene_indices !== nothing

    last_added_gene_index = nothing
    for (test_index, factor_gene_index) in enumerate(ordered_factor_gene_indices)
        if try_add_factor!(context, test_index, factor_gene_index, last_removed_gene_index)
            order_predictive!(context)
            last_added_gene_index = factor_gene_index
        end
    end
    return last_added_gene_index
end

function try_add_factor!(
    context::Context,
    test_index::Integer,
    factor_gene_index::Integer,
    last_removed_gene_index::Maybe{<:Integer},
)::Bool
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    if factor_gene_index == last_removed_gene_index || factor_gene_index in predictive_genes.indices
        return false
    end

    current_cross_validation = context.cross_validation
    @assert current_cross_validation !== nothing
    current_n_predictive_genes = predictive_genes.n_entries

    predictive_genes = push_predictive!(context, factor_gene_index)
    n_predictive_genes = predictive_genes.n_entries
    @assert n_predictive_genes == current_n_predictive_genes + 1

    analyze_left_outs!(context)

    cross_validation = context.cross_validation
    @assert context.cross_validation !== nothing

    improvement =
        current_cross_validation.mean_rms * (1 + 1e-2)^current_n_predictive_genes -
        cross_validation.mean_rms * (1 + 1e-2)^n_predictive_genes

    if improvement > 0
        print(
            stderr,
            "TODOX # $(current_n_predictive_genes) ++ $(test_index) / $(length(context.ordered_factor_gene_indices)) / $(context.data.factor_genes.n_entries) $(context.data.names_of_genes[factor_gene_index]) -> $(context.cross_validation.mean_rms) > $(current_cross_validation.mean_rms) ...          \r",
        )
        return true
    end

    print(
        stderr,
        "TODOX # $(current_n_predictive_genes) ?+ $(test_index) / $(length(context.ordered_factor_gene_indices)) / $(context.data.factor_genes.n_entries) $(context.data.names_of_genes[factor_gene_index]) -> $(context.cross_validation.mean_rms) <~ $(current_cross_validation.mean_rms) ...          \r",
    )

    pop_predictive!(context, factor_gene_index)
    context.cross_validation = current_cross_validation
    return false
end

function try_remove_factors!(context::Context, last_added_gene_index::Maybe{<:Integer})::Maybe{<:Integer}
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing
    predictive_gene_indices = copy_array(predictive_genes.indices)

    last_removed_gene_index = nothing
    for (predictive_index, factor_gene_index) in Iterators.reverse(enumerate(predictive_gene_indices))
        if try_remove_factor!(context, predictive_index, factor_gene_index, last_added_gene_index)
            order_predictive!(context)
            last_removed_gene_index = factor_gene_index
        end
    end
    return last_removed_gene_index
end

function try_remove_factor!(
    context,
    predictive_index::Integer,
    factor_gene_index::Integer,
    last_added_gene_index::Maybe{<:Integer},
)::Bool
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    if factor_gene_index == last_added_gene_index
        return false
    end

    current_predictive_genes_indices = copy_array(predictive_genes.indices)
    current_n_predictive_genes = predictive_genes.n_entries

    current_cross_validation = context.cross_validation
    @assert current_cross_validation !== nothing

    remove_predictive!(context, factor_gene_index)

    n_predictive_genes = context.predictive_genes.n_entries
    @assert n_predictive_genes == current_n_predictive_genes - 1

    analyze_left_outs!(context)

    cross_validation = context.cross_validation
    @assert context.cross_validation !== nothing

    improvement =
        current_cross_validation.mean_rms * (1 + 1e-2)^current_n_predictive_genes -
        cross_validation.mean_rms * (1 + 1e-2)^n_predictive_genes
    if improvement >= 0 # -1e-3
        print(
            stderr,
            "TODOX # $(current_n_predictive_genes) -- $(predictive_index) / $(current_n_predictive_genes) / $(context.data.factor_genes.n_entries) $(context.data.names_of_genes[factor_gene_index]) -> $(context.cross_validation.mean_rms) <= $(current_cross_validation.mean_rms) ...          \r",
        )
        return true
    end

    print(
        stderr,
        "TODOX # $(current_n_predictive_genes) ?- $(predictive_index) / $(current_n_predictive_genes) / $(context.data.factor_genes.n_entries) $(context.data.names_of_genes[factor_gene_index]) -> $(context.cross_validation.mean_rms) > $(current_cross_validation.mean_rms) ...          \r",
    )

    reset_predictive!(context, current_predictive_genes_indices)
    context.cross_validation = current_cross_validation
    return false
end

function analyze_left_outs!(context::Context; use_current_least_squares::Bool = false)::Nothing
    if use_current_least_squares
        current_least_squares = context.least_squares
        @assert current_least_squares !== nothing
    else
        current_least_squares = nothing
    end

    included_metacells = context.included_metacells
    @assert included_metacells !== nothing

    core_metacells = context.core_metacells
    if core_metacells === nothing
        core_metacells = included_metacells
    end

    core_metacells_indices = findall(core_metacells.mask)
    shuffle!(context.parameters.rng, core_metacells_indices)
    reset_left_out_buffers!(context, core_metacells_indices)

    chunk_size = max(div(core_metacells.n_entries, 5), 1)
    n_chunks = ceil(core_metacells.n_entries / chunk_size)
    chunk_size = core_metacells.n_entries / n_chunks

    selected_metacells_mask = similar(included_metacells.mask)

    for chunk_index in 1:n_chunks
        first_left_out_position = Int(round((chunk_index - 1) * chunk_size)) + 1
        last_left_out_position = Int(round(chunk_index * chunk_size))
        left_out_positions = first_left_out_position:last_left_out_position

        chunk_metacells_indices = core_metacells_indices[left_out_positions]

        chunk_metacells = reset_chunk!(context, chunk_metacells_indices)

        selected_metacells_mask .= included_metacells.mask
        selected_metacells_mask[chunk_metacells.indices] .= false
        selected_metacells = SomeMask(selected_metacells_mask, included_metacells.n_entries - chunk_metacells.n_entries)
        if current_least_squares === nothing
            solve_least_squares!(context, selected_metacells)
        else
            context.least_squares = current_least_squares
        end

        collect_chunk!(context, left_out_positions)

        clear_chunk!(context)
    end

    evaluate_left_out!(context)

    clear_left_out_buffers!(context)

    context.least_squares = current_least_squares

    return nothing
end

function solve_least_squares!(context::Context, selected_metacells::Union{SomeMask, SomeIndices})::LeastSquares
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    if selected_metacells isa SomeMask
        log_fractions_of_genes_in_selected_metacells =
            context.data.log_fractions_of_genes_in_metacells[:, selected_metacells.mask]
    else
        log_fractions_of_genes_in_selected_metacells =
            context.data.log_fractions_of_genes_in_metacells[:, selected_metacells.indices]
    end
    @assert_matrix(
        log_fractions_of_genes_in_selected_metacells,
        context.data.n_genes,
        selected_metacells.n_entries,
        Columns
    )

    log_fractions_in_selected_metacells_of_genes = transposer(log_fractions_of_genes_in_selected_metacells)
    @assert_matrix(
        log_fractions_in_selected_metacells_of_genes,
        selected_metacells.n_entries,
        context.data.n_genes,
        Columns
    )

    log_fractions_in_selected_metacells_of_predictive_genes =
        log_fractions_in_selected_metacells_of_genes[:, predictive_genes.indices]
    @assert_matrix(
        log_fractions_in_selected_metacells_of_predictive_genes,
        selected_metacells.n_entries,
        predictive_genes.n_entries,
        Columns
    )

    mean_log_fractions_of_genes = vec(vmean(log_fractions_in_selected_metacells_of_genes; dims = 1))  # NOJET
    @assert_vector(mean_log_fractions_of_genes, context.data.n_genes)

    mean_log_fractions_of_predictive_genes = mean_log_fractions_of_genes[predictive_genes.indices]
    @assert_vector(mean_log_fractions_of_predictive_genes, predictive_genes.n_entries)

    log_fractions_in_selected_metacells_of_genes .-= transpose(mean_log_fractions_of_genes)
    log_fractions_in_selected_metacells_of_predictive_genes .-= transpose(mean_log_fractions_of_predictive_genes)

    divergence_of_genes = context.data.divergence_of_genes
    divergence_of_predictive_genes = divergence_of_genes[predictive_genes.indices]

    log_fractions_in_selected_metacells_of_genes .*= transpose(1.0 .- divergence_of_genes)
    log_fractions_in_selected_metacells_of_predictive_genes .*= transpose(1.0 .- divergence_of_predictive_genes)

    # Least squares: log_fractions_in_selected_metacells_of_predictive_genes \ log_fractions_in_selected_metacells_of_genes
    coefficients_of_predictive_genes_of_genes = nonneg_lsq(
        log_fractions_in_selected_metacells_of_predictive_genes,
        log_fractions_in_selected_metacells_of_genes;
        alg = :fnnls,
    )
    @assert_matrix(coefficients_of_predictive_genes_of_genes, predictive_genes.n_entries, context.data.n_genes, Columns)

    context.least_squares =
        least_squares = LeastSquares(;
            mean_log_fractions_of_genes = mean_log_fractions_of_genes,
            coefficients_of_genes_of_predictive_genes = transposer(coefficients_of_predictive_genes_of_genes),
        )

    return least_squares
end

function predict_by_least_squares(
    context::Context,
    selected_metacells::Union{SomeMask, SomeIndices},
)::AbstractMatrix{<:AbstractFloat}
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    least_squares = context.least_squares
    @assert least_squares !== nothing

    if selected_metacells isa SomeMask
        log_fractions_of_predictive_genes_in_selected_metacells =
            context.data.log_fractions_of_genes_in_metacells[predictive_genes.indices, selected_metacells.mask]
    else
        log_fractions_of_predictive_genes_in_selected_metacells =
            context.data.log_fractions_of_genes_in_metacells[predictive_genes.indices, selected_metacells.indices]
    end
    @assert_matrix(
        log_fractions_of_predictive_genes_in_selected_metacells,
        predictive_genes.n_entries,
        selected_metacells.n_entries,
        Columns
    )

    mean_log_fractions_of_predictive_genes = least_squares.mean_log_fractions_of_genes[predictive_genes.indices]
    @assert_vector(mean_log_fractions_of_predictive_genes, predictive_genes.n_entries)

    log_fractions_of_predictive_genes_in_selected_metacells .-= mean_log_fractions_of_predictive_genes

    divergence_of_genes = context.data.divergence_of_genes
    divergence_of_predictive_genes = divergence_of_genes[predictive_genes.indices]

    log_fractions_of_predictive_genes_in_selected_metacells .*= 1.0 .- divergence_of_predictive_genes

    predicted_log_fractions_of_genes_in_selected_metacells =
        least_squares.coefficients_of_genes_of_predictive_genes *
        log_fractions_of_predictive_genes_in_selected_metacells
    @assert_matrix(
        predicted_log_fractions_of_genes_in_selected_metacells,
        context.data.n_genes,
        selected_metacells.n_entries
    )

    scale_of_genes = 1.0 ./ (1.0 .- divergence_of_genes)
    predicted_log_fractions_of_genes_in_selected_metacells .*= scale_of_genes

    predicted_log_fractions_of_genes_in_selected_metacells .+= least_squares.mean_log_fractions_of_genes

    return predicted_log_fractions_of_genes_in_selected_metacells
end

function collect_chunk!(context::Context, left_out_positions::UnitRange{<:Integer})::Nothing
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    chunk_metacells = context.chunk_metacells
    @assert chunk_metacells !== nothing

    left_out_buffers = context.left_out_buffers
    @assert left_out_buffers !== nothing

    log_fractions_of_genes_in_chunk_metacells =
        context.data.log_fractions_of_genes_in_metacells[:, chunk_metacells.indices]
    @assert_matrix(log_fractions_of_genes_in_chunk_metacells, context.data.n_genes, chunk_metacells.n_entries, Columns)

    predicted_log_fractions_of_genes_in_chunk_metacells = predict_by_least_squares(context, chunk_metacells)
    @assert_matrix(
        predicted_log_fractions_of_genes_in_chunk_metacells,
        context.data.n_genes,
        chunk_metacells.n_entries,
        Columns
    )

    left_out_buffers.log_fractions_of_genes_in_left_out_metacells[:, left_out_positions] .=
        log_fractions_of_genes_in_chunk_metacells
    left_out_buffers.predicted_log_fractions_of_genes_in_left_out_metacells[:, left_out_positions] .=
        predicted_log_fractions_of_genes_in_chunk_metacells

    return nothing
end

function evaluate_left_out!(context::Context)::Nothing
    left_out_metacells = context.left_out_metacells
    @assert left_out_metacells !== nothing

    left_out_buffers = context.left_out_buffers
    @assert left_out_buffers !== nothing

    log_fractions_in_left_out_metacells_of_genes =
        transposer(left_out_buffers.log_fractions_of_genes_in_left_out_metacells)
    @assert_matrix(
        log_fractions_in_left_out_metacells_of_genes,
        left_out_metacells.n_entries,
        context.data.n_genes,
        Columns
    )

    predicted_log_fractions_in_left_out_metacells_of_genes =
        transposer(left_out_buffers.predicted_log_fractions_of_genes_in_left_out_metacells)
    @assert_matrix(
        predicted_log_fractions_in_left_out_metacells_of_genes,
        left_out_metacells.n_entries,
        context.data.n_genes,
        Columns
    )

    r2_of_genes = Vector{Float32}(undef, context.data.n_genes)
    shuffled_r2_of_genes = Vector{Float32}(undef, context.data.n_genes)
    rms_of_genes = Vector{Float32}(undef, context.data.n_genes)

    @threads for gene_index in 1:(context.data.n_genes)
        @views log_fractions_in_left_out_metacells_of_gene = log_fractions_in_left_out_metacells_of_genes[:, gene_index]
        @views predicted_log_fractions_in_left_out_metacells_of_gene =
            predicted_log_fractions_in_left_out_metacells_of_genes[:, gene_index]

        residual_log_fractions_in_left_out_metacells_of_gene =
            predicted_log_fractions_in_left_out_metacells_of_gene .- log_fractions_in_left_out_metacells_of_gene
        residual_log_fractions_in_left_out_metacells_of_gene .*= residual_log_fractions_in_left_out_metacells_of_gene
        rms_of_genes[gene_index] = sqrt(vmean(residual_log_fractions_in_left_out_metacells_of_gene))

        correlation_of_gene =
            vcor(log_fractions_in_left_out_metacells_of_gene, predicted_log_fractions_in_left_out_metacells_of_gene)
        if isnan(correlation_of_gene)
            correlation_of_gene = 0
        end

        r2_of_genes[gene_index] = correlation_of_gene * correlation_of_gene

        n_shuffles = 6
        shuffled_r2s_of_gene = Vector{Float32}(undef, n_shuffles)
        for shuffle_index in 1:n_shuffles
            shuffle!(context.parameters.rng, predicted_log_fractions_in_left_out_metacells_of_gene)
            shuffled_correlation_of_gene =
                vcor(log_fractions_in_left_out_metacells_of_gene, predicted_log_fractions_in_left_out_metacells_of_gene)
            if isnan(shuffled_correlation_of_gene)
                shuffled_correlation_of_gene = 0
            end
            shuffled_r2s_of_gene[shuffle_index] = shuffled_correlation_of_gene * shuffled_correlation_of_gene
        end
        shuffled_r2_of_gene = max(maximum(shuffled_r2s_of_gene), 0.0)
        shuffled_r2_of_genes[gene_index] = shuffled_r2_of_gene
    end

    mean_r2 = vmean(r2_of_genes)
    mean_shuffled_r2 = vmean(shuffled_r2_of_genes)
    mean_extra_r2 = mean_r2 - mean_shuffled_r2

    mean_rms = vmean(rms_of_genes)

    context.cross_validation = CrossValidation(;
        mean_r2 = mean_r2,
        mean_shuffled_r2 = mean_shuffled_r2,
        mean_extra_r2 = mean_extra_r2,
        mean_rms = mean_rms,
    )
    return nothing
end

function final_analysis!(context; use_current_least_squares::Bool = false)::Analysis
    included_metacells = context.included_metacells
    @assert included_metacells !== nothing

    core_metacells = context.core_metacells
    if core_metacells === nothing
        core_metacells = included_metacells
    end

    if use_current_least_squares
        current_least_squares = context.least_squares
        @assert current_least_squares !== nothing
    else
        current_least_squares = nothing
    end

    order_predictive!(context)

    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    cross_validation = context.cross_validation
    @assert cross_validation !== nothing

    if current_least_squares !== nothing
        context.least_squares = least_squares = current_least_squares
    else
        least_squares = solve_least_squares!(context, included_metacells)
    end

    variability = compute_variability!(context)

    context.least_squares = current_least_squares

    return Analysis(;
        core_metacells = core_metacells,
        included_metacells = included_metacells,
        predictive_genes = predictive_genes,
        least_squares = least_squares,
        variability = variability,
        cross_validation = cross_validation,
    )
end

function compute_variability!(context::Context)::Variability
    included_metacells = context.included_metacells
    @assert included_metacells !== nothing

    core_metacells = context.core_metacells
    if core_metacells === nothing
        core_metacells = included_metacells
    end

    least_squares = context.least_squares
    @assert least_squares !== nothing

    log_fractions_in_core_metacells_of_genes =
        transposer(context.data.log_fractions_of_genes_in_metacells[:, core_metacells.mask])
    @assert_matrix(log_fractions_in_core_metacells_of_genes, core_metacells.n_entries, context.data.n_genes, Columns)

    rms_of_genes = Vector{Float32}(undef, context.data.n_genes)

    @threads for gene_index in 1:(context.data.n_genes)
        @views log_fractions_in_core_metacells_of_gene = log_fractions_in_core_metacells_of_genes[:, gene_index]
        residual_log_fractions_in_core_metacells_of_gene =
            log_fractions_in_core_metacells_of_gene .- least_squares.mean_log_fractions_of_genes[gene_index]
        residual_log_fractions_in_core_metacells_of_gene .*= residual_log_fractions_in_core_metacells_of_gene
        rms_of_genes[gene_index] = sqrt(vmean(residual_log_fractions_in_core_metacells_of_gene))
    end

    context.included_variability = variability = Variability(; mean_rms = vmean(rms_of_genes))

    return variability
end

function order_predictive!(context::Context)::Nothing
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    minimal_rank_of_genes = context.minimal_rank_of_genes
    @assert minimal_rank_of_genes !== nothing

    minimal_rank_of_predictive_genes = minimal_rank_of_genes[predictive_genes.indices]
    predictive_genes_order = sortperm(minimal_rank_of_predictive_genes)
    reset_predictive!(context, predictive_genes.indices[predictive_genes_order])
    return nothing
end

"""
    function compute_local_predictive_factors!(  # untested
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        max_principal_components::Integer = $(DEFAULT.max_principal_components),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        max_block_span::AbstractFloat = $(DEFAULT.max_block_span),
        min_blocks_in_neighborhood::Integer = $(DEFAULT.min_blocks_in_neighborhood),
        min_metacells_in_neighborhood::Integer = $(DEFAULT.min_metacells_in_neighborhood),
        rng::AbstractRNG = default_rng(),
    )::Nothing

TODOX
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), block_axis(GuaranteedOutput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_transcription_factor_vector(RequiredInput),
        gene_is_global_predictive_factor_vector(RequiredInput),
        metacell_total_UMIs_vector(RequiredInput),
        gene_metacell_total_UMIs_matrix(RequiredInput),
        metacell_block_vector(GuaranteedOutput),
        block_block_distance_matrix(GuaranteedOutput),
    ],
) function compute_local_predictive_factors!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    max_principal_components::Integer = 30,
    fold_confidence::AbstractFloat = 0.9,
    max_block_span::AbstractFloat = function_default(identify_marker_genes!, :min_marker_gene_range_fold),
    min_blocks_in_neighborhood::Integer = 4,
    min_metacells_in_neighborhood::Integer = 100,
    rng::AbstractRNG = default_rng(),
)::Nothing
    @assert min_significant_gene_UMIs >= 0
    @assert gene_fraction_regularization >= 0
    @assert 0 <= fold_confidence <= 1
    @assert max_block_span > 0
    @assert min_blocks_in_neighborhood > 0
    @assert min_metacells_in_neighborhood > 0

    parameters = ParametersContext(;
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        gene_fraction_regularization = gene_fraction_regularization,
        max_principal_components = max_principal_components,
        fold_confidence = fold_confidence,
        max_block_span = max_block_span,
        min_blocks_in_neighborhood = min_blocks_in_neighborhood,
        min_metacells_in_neighborhood = min_metacells_in_neighborhood,
        rng = rng,
    )
    context = Context(; parameters = parameters, data = DataContext(daf, parameters))
    order_factors!(context)
    load_predictive!(context, daf)
    compute_confidence!(context, daf)

    blocks = compute_blocks!(context)
    block_names = group_names(daf, "metacell", blocks.metacells_of_blocks; prefix = "B")
    blocks_of_metacells = block_names[blocks.blocks_of_metacells]
    add_axis!(daf, "block", block_names)
    set_vector!(daf, "metacell", "block", blocks_of_metacells)
    set_matrix!(daf, "block", "block", "distance", blocks.distances_between_blocks)

    compute_environments!(context)

    return nothing
end

@logged function load_predictive!(context::Context, daf::DafReader)::SomeIndices
    predictive_genes_mask = get_vector(daf, "gene", "is_global_predictive_factor")
    predictive_genes_indices = findall(predictive_genes_mask)
    predictive_genes = SomeIndices(predictive_genes_indices, length(predictive_genes_indices))
    @assert predictive_genes.n_entries > 0
    context.predictive_genes = predictive_genes
    order_predictive!(context)
    @debug "global predictive genes: [ $(join(context.data.names_of_genes[context.predictive_genes.indices], ", ")) ]"
    open("global_predictive_factors.csv", "w") do file
        println(file, "gene_index,gene_name")
        for gene_index in context.predictive_genes.indices
            println(file, "$(gene_index),$(context.data.names_of_genes[gene_index])")
        end
    end
    return predictive_genes
end

@logged function compute_confidence!(context::Context, daf::DafReader)::Confidence
    total_UMIs_of_metacells = get_vector(daf, "metacell", "total_UMIs").array
    total_UMIs_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    @assert_matrix(total_UMIs_of_genes_in_metacells, context.data.n_genes, context.data.n_metacells, Columns)

    log_decreased_fractions_of_genes_in_metacells, log_increased_fractions_of_genes_in_metacells =
        compute_confidence_log_fraction_of_genes_in_metacells(;
            gene_fraction_regularization = context.parameters.gene_fraction_regularization,
            fractions_of_genes_in_metacells = context.data.fractions_of_genes_in_metacells,
            total_UMIs_of_metacells = total_UMIs_of_metacells,
            fold_confidence = context.parameters.fold_confidence,
        )

    @assert_matrix(
        log_decreased_fractions_of_genes_in_metacells,
        context.data.n_genes,
        context.data.n_metacells,
        Columns
    )

    @assert_matrix(
        log_increased_fractions_of_genes_in_metacells,
        context.data.n_genes,
        context.data.n_metacells,
        Columns
    )

    context.confidence =
        confidence = Confidence(;
            total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
            log_decreased_fractions_of_genes_in_metacells = log_decreased_fractions_of_genes_in_metacells,
            log_increased_fractions_of_genes_in_metacells = log_increased_fractions_of_genes_in_metacells,
        )

    return confidence
end

function compute_confidence_log_fraction_of_genes_in_metacells(;
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

@logged function compute_blocks!(context::Context)::Blocks
    distances_between_metacells = compute_distances_between_metacells(context)

    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    blocks_of_metacells = Vector{UInt32}(cutree(clusters; h = 0))
    metacells_of_blocks = collect_group_members(blocks_of_metacells)
    n_blocks = length(metacells_of_blocks)

    distances_between_blocks = compute_distances_between_blocks(distances_between_metacells, metacells_of_blocks)

    context.blocks =
        blocks = Blocks(;
            blocks_of_metacells = blocks_of_metacells,
            metacells_of_blocks = metacells_of_blocks,
            n_blocks = n_blocks,
            distances_between_blocks = distances_between_blocks,
        )

    open("blocks_of_metacells.csv", "w") do file
        println(file, "metacell_index,block_index")
        for (metacell_index, block_index) in enumerate(blocks_of_metacells)
            println(file, "$(metacell_index),$(block_index)")
        end
    end

    open("distances_between_blocks.csv", "w") do file
        println(file, "left_block_index,right_block_index,distance")
        for left_block_index in 1:n_blocks
            for right_block_index in 1:n_blocks
                println(
                    file,
                    "$(left_block_index),$(right_block_index),$(distances_between_blocks[left_block_index, right_block_index])",
                )
            end
        end
    end

    return blocks
end

function compute_distances_between_metacells(context::Context)::AbstractMatrix{<:AbstractFloat}
    predictive_genes = context.predictive_genes
    @assert predictive_genes !== nothing

    distances_between_metacells = Matrix{Float32}(undef, context.data.n_metacells, context.data.n_metacells)

    divergence_of_predictive_genes = context.data.divergence_of_genes[predictive_genes.indices]

    log_decreased_fractions_of_predictive_genes_in_metacells =
        context.confidence.log_decreased_fractions_of_genes_in_metacells[predictive_genes.indices, :]
    log_increased_fractions_of_predictive_genes_in_metacells =
        context.confidence.log_increased_fractions_of_genes_in_metacells[predictive_genes.indices, :]
    total_UMIs_of_predictive_genes_in_metacells =
        context.confidence.total_UMIs_of_genes_in_metacells[predictive_genes.indices, :]

    distances_between_metacells[1, 1] = 0.0
    @threads for base_metacell_index in reverse(2:(context.data.n_metacells))
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
                context.parameters.min_significant_gene_UMIs,
                total_UMIs_of_predictive_genes_in_base_metacell,
                log_decreased_fractions_of_predictive_genes_in_base_metacell,
                log_increased_fractions_of_predictive_genes_in_base_metacell,
                total_UMIs_of_predictive_genes_in_other_metacells,
                log_decreased_fractions_of_predictive_genes_in_other_metacells,
                log_increased_fractions_of_predictive_genes_in_other_metacells,
                divergence_of_predictive_genes,
            )
        @assert_matrix(
            significant_fold_factors_of_predictive_genes_in_other_metacells,
            predictive_genes.n_entries,
            base_metacell_index - 1,
            Columns
        )

        distances_between_base_and_other_metacells = vec(
            vmean(
                max.(
                    significant_fold_factors_of_predictive_genes_in_other_metacells .-
                    context.parameters.max_block_span,
                    0.0,
                );
                dims = 1,
            ),
        )
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

function compute_distances_between_blocks(
    distances_between_metacells::Matrix{Float32},
    metacells_of_blocks::AbstractVector{<:AbstractVector{<:Integer}},
)::Matrix{Float32}
    n_metacells = size(distances_between_metacells, 1)
    @assert_matrix(distances_between_metacells, n_metacells, n_metacells, Columns)

    n_blocks = length(metacells_of_blocks)
    distances_between_blocks = Matrix{Float32}(undef, n_blocks, n_blocks)

    distances_between_blocks[1, 1] = 0.0
    @threads for base_block_index in reverse(2:n_blocks)
        distances_between_blocks[base_block_index, base_block_index] = 0.0
        metacells_of_base_block = metacells_of_blocks[base_block_index]

        distance_of_metacells_from_base_block_metacells = distances_between_metacells[:, metacells_of_base_block]

        distance_of_metacells_from_base_block = vec(vmean(distance_of_metacells_from_base_block_metacells; dims = 2))
        @assert length(distance_of_metacells_from_base_block) == n_metacells

        for other_block_index in 1:(base_block_index - 1)
            metacells_of_other_block = metacells_of_blocks[other_block_index]

            @views distance_of_other_block_metacells_from_base_block =
                distance_of_metacells_from_base_block[metacells_of_other_block]

            distance_between_other_and_base_block = mean(distance_of_other_block_metacells_from_base_block)

            distances_between_blocks[base_block_index, other_block_index] = distance_between_other_and_base_block
            distances_between_blocks[other_block_index, base_block_index] = distance_between_other_and_base_block
        end
    end

    return distances_between_blocks
end

@logged function compute_environments!(context::Context)::AbstractVector{Environment}
    blocks = context.blocks
    @assert blocks !== nothing

    @assert context.parameters.min_blocks_in_neighborhood <= blocks.n_blocks
    @assert context.parameters.min_metacells_in_neighborhood <= context.data.n_metacells

    included_metacells = reset_included!(context, ones(Bool, context.data.n_metacells); keep_predictive = true)
    global_least_squares = solve_least_squares!(context, included_metacells)
    global_predictive_genes = context.predictive_genes
    @assert global_predictive_genes !== nothing

    global_predictive_set = Set(global_predictive_genes.indices)
    total_predictive_set = Set(global_predictive_genes.indices)

    environments = Vector{Environment}(undef, blocks.n_blocks)

    uses_of_predictive_genes = Dict{UInt32, UInt32}()
    for predictive_gene_index in global_predictive_genes.indices
        uses_of_predictive_genes[predictive_gene_index] = 0
    end

    for block_index in 1:(blocks.n_blocks)
        environments[block_index] =
            environment =
                compute_environment_of_block(context, global_predictive_genes, global_least_squares, block_index)

        local_predictive_set = Set(environment.local_analysis.predictive_genes.indices)
        for predictive_gene_index in environment.local_analysis.predictive_genes.indices
            if haskey(uses_of_predictive_genes, predictive_gene_index)
                uses_of_predictive_genes[predictive_gene_index] += 1
            else
                uses_of_predictive_genes[predictive_gene_index] = 1
            end
        end

        total_predictive_set = union(total_predictive_set, local_predictive_set)

        total = 0
        reused = 0
        added = 0
        removed = 0
        for predictive_gene_index in total_predictive_set
            is_in_global = predictive_gene_index in global_predictive_set
            is_in_local = predictive_gene_index in local_predictive_set
            if is_in_local
                total += 1
                if is_in_global
                    reused += 1
                else
                    added += 1
                end
            else
                if is_in_global
                    removed += 1
                end
            end
        end

        n_block_metacells = sum(blocks.blocks_of_metacells .== block_index)
        @debug "TODOX $(block_index): $(n_block_metacells) metacells => $(depict(environments[block_index]))"
        @debug "TODOX   Predictive genes: total: $(total) reused: $(reused) added: $(added) removed: $(removed)"
        @debug "TODOX   Variability RMS: global: $(@sprintf("%.5f", environment.global_analysis.variability.mean_rms)) by_global: $(@sprintf("%.5f", environment.by_global_analysis.variability.mean_rms)) local: $(@sprintf("%.5f", environment.local_analysis.variability.mean_rms))"
        @debug "TODOX   Analysis RMS: global: $(@sprintf("%.5f", environment.global_analysis.cross_validation.mean_rms)) by_global: $(@sprintf("%.5f", environment.by_global_analysis.cross_validation.mean_rms)) local: $(@sprintf("%.5f", environment.local_analysis.cross_validation.mean_rms))"
        @debug "TODOX   Analysis R2 global: $(@sprintf("%.5f", environment.global_analysis.cross_validation.mean_r2)) by_global: $(@sprintf("%.5f", environment.by_global_analysis.cross_validation.mean_r2)) local: $(@sprintf("%.5f", environment.local_analysis.cross_validation.mean_r2))"
    end

    open("blocks_qc.csv", "w") do file
        println(
            file,
            "block_index,environment_blocks,block_metacells,neighborhood_metacells,environment_metacells,global_rms,by_global_rms,local_rms,global_r2,by_global_r2,local_r2",
        )
        for block_index in 1:(blocks.n_blocks)
            n_block_metacells = sum(blocks.blocks_of_metacells .== block_index)
            environment = environments[block_index]
            n_environment_blocks = environment.region.blocks.n_entries
            n_neighborhood_metacells = 0
            for block_index in environment.neighborhood.blocks.indices
                n_neighborhood_metacells += sum(blocks.blocks_of_metacells .== block_index)
            end
            n_environment_metacells = 0
            for block_index in environment.region.blocks.indices
                n_environment_metacells += sum(blocks.blocks_of_metacells .== block_index)
            end
            println(
                file,
                "$(block_index),$(n_environment_blocks),$(n_block_metacells),$(n_neighborhood_metacells),$(n_environment_metacells),$(environment.global_analysis.cross_validation.mean_rms),$(environment.by_global_analysis.cross_validation.mean_rms),$(environment.local_analysis.cross_validation.mean_rms),$(environment.global_analysis.cross_validation.mean_r2),$(environment.by_global_analysis.cross_validation.mean_r2),$(environment.local_analysis.cross_validation.mean_r2)",
            )
        end
    end

    open("blocks_coefficients.csv", "w") do file
        println(file, "block_index,predictive_gene_name,predicted_gene_name,coefficient")
        for block_index in 1:(blocks.n_blocks)
            local_analysis = environments[block_index].local_analysis
            predictive_genes = local_analysis.predictive_genes
            coefficients_of_genes_of_predictive_genes =
                local_analysis.least_squares.coefficients_of_genes_of_predictive_genes
            for (predictive_index, predictive_gene_index) in enumerate(predictive_genes.indices)
                @views coefficients_of_genes_of_predictive_gene =
                    coefficients_of_genes_of_predictive_genes[:, predictive_index]
                genes_order = sortperm(abs.(coefficients_of_genes_of_predictive_gene); rev = true)
                for gene_index in genes_order
                    coefficient = coefficients_of_genes_of_predictive_gene[gene_index]
                    if abs(coefficient) == 0
                        break
                    end
                    println(
                        file,
                        "$(block_index),$(context.data.names_of_genes[predictive_gene_index]),$(context.data.names_of_genes[gene_index]),$(coefficient)",
                    )
                end
            end
        end
    end

    open("blocks_neighborhoods.csv", "w") do file
        println(file, "seed_block_index,near_block_index")
        for block_index in 1:(blocks.n_blocks)
            for near_index in environments[block_index].neighborhood.blocks.indices
                println(file, "$(block_index),$(near_index)")
            end
        end
    end

    open("blocks_environments.csv", "w") do file
        println(file, "seed_block_index,near_block_index")
        for block_index in 1:(blocks.n_blocks)
            for near_index in environments[block_index].region.blocks.indices
                println(file, "$(block_index),$(near_index)")
            end
        end
    end

    open("blocks_predictive_factors.csv", "w") do file
        println(file, "block_index,gene_name")
        for block_index in 1:(blocks.n_blocks)
            for gene_index in environments[block_index].local_analysis.predictive_genes.indices
                println(file, "$(block_index),$(context.data.names_of_genes[gene_index])")
            end
        end
    end

    predictive_gene_indices = collect(total_predictive_set)
    reset_predictive!(context, predictive_gene_indices)
    order_predictive!(context)

    open("local_predictive_factors.csv", "w") do file
        println(file, "gene_name,is_global,used_in_blocks")
        for gene_index in context.predictive_genes.indices
            is_global = gene_index in global_predictive_set
            println(
                file,
                "$(context.data.names_of_genes[gene_index]),$(is_global),$(uses_of_predictive_genes[gene_index])",
            )
        end
    end

    return environments
end

function compute_environment_of_block(
    context::Context,
    global_predictive_genes::SomeIndices,
    global_least_squares::LeastSquares,
    block_index::Integer,
)::Environment
    blocks = context.blocks
    @assert blocks !== nothing

    distances_between_other_and_base_blocks = blocks.distances_between_blocks[:, block_index]
    @assert_vector(distances_between_other_and_base_blocks, blocks.n_blocks)
    ordered_block_indices = sortperm(distances_between_other_and_base_blocks)
    @assert ordered_block_indices[1] == block_index

    n_metacells_in_neighborhood = 0
    n_blocks_in_neighborhood = 0
    neighborhood_metacells_mask = nothing
    while n_blocks_in_neighborhood < context.parameters.min_blocks_in_neighborhood ||
        n_metacells_in_neighborhood < context.parameters.min_metacells_in_neighborhood
        n_blocks_in_neighborhood += 1
        neighborhood_metacells_mask = region_metacells_mask(context, n_blocks_in_neighborhood, ordered_block_indices)
        @assert neighborhood_metacells_mask !== nothing
        reset_core!(context, neighborhood_metacells_mask; keep_predictive = true)
        n_metacells_in_neighborhood = context.core_metacells.n_entries
    end

    analysis_of_environments = Vector{Maybe{Analysis}}(undef, blocks.n_blocks)
    analysis_of_environments .= nothing

    n_blocks_of_results = Vector{Int32}()
    mean_rms_of_results = Vector{Float32}()
    n_results = 0

    best_n_blocks = next_n_blocks = n_blocks_in_neighborhood
    while true
        @assert analysis_of_environments[next_n_blocks] === nothing
        analysis_of_environments[next_n_blocks] =
            analysis = analyze_environment(context, n_blocks_in_neighborhood, next_n_blocks, ordered_block_indices)
        push!(n_blocks_of_results, next_n_blocks)
        push!(mean_rms_of_results, analysis.cross_validation.mean_rms)
        n_results += 1
        order = sortperm(n_blocks_of_results)
        n_blocks_of_results .= n_blocks_of_results[order]
        mean_rms_of_results .= mean_rms_of_results[order]

        best_result_index = argmin(mean_rms_of_results)
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
            analysis_of_environments[high_n_blocks] !== nothing &&
            analysis_of_environments[low_n_blocks] !== nothing
        )
            high_n_blocks = min(blocks.n_blocks, best_n_blocks + 10)
            has_all = true
            for test_n_blocks in reverse(best_n_blocks:high_n_blocks)
                if analysis_of_environments[test_n_blocks] === nothing
                    next_n_blocks = test_n_blocks
                    has_all = false
                    break
                end
            end
            if has_all
                break
            end

        elseif analysis_of_environments[low_n_blocks] !== nothing
            @assert analysis_of_environments[high_n_blocks] === nothing
            next_n_blocks = high_n_blocks

        elseif analysis_of_environments[high_n_blocks] !== nothing
            @assert analysis_of_environments[low_n_blocks] === nothing
            next_n_blocks = low_n_blocks

        elseif high_diff > low_diff
            @assert analysis_of_environments[high_n_blocks] === nothing
            next_n_blocks = high_n_blocks

        else
            @assert analysis_of_environments[low_n_blocks] === nothing
            next_n_blocks = low_n_blocks
        end
    end

    @assert best_n_blocks >= n_blocks_in_neighborhood
    by_global_analysis = analysis_of_environments[best_n_blocks]

    neighborhood_blocks = SomeIndices(ordered_block_indices[1:n_blocks_in_neighborhood], n_blocks_in_neighborhood)
    neighborhood_metacells = SomeMask(neighborhood_metacells_mask, n_metacells_in_neighborhood)

    environment_blocks = SomeIndices(ordered_block_indices[1:best_n_blocks], best_n_blocks)
    environment_metacells_mask = region_metacells_mask(context, best_n_blocks, ordered_block_indices)
    environment_metacells = SomeMask(environment_metacells_mask, sum(environment_metacells_mask))

    @assert neighborhood_metacells_mask !== nothing
    reset_core!(context, neighborhood_metacells_mask; keep_predictive = true)
    reset_included!(context, environment_metacells_mask; keep_predictive = true)

    local_analysis = analyze_without_overfitting!(context)

    context.predictive_genes = global_predictive_genes
    context.least_squares = global_least_squares

    analyze_left_outs!(context; use_current_least_squares = true)
    global_analysis = final_analysis!(context; use_current_least_squares = true)

    return Environment(;
        region = Region(; blocks = environment_blocks, metacells = environment_metacells),
        neighborhood = Region(; blocks = neighborhood_blocks, metacells = neighborhood_metacells),
        local_analysis = local_analysis,
        by_global_analysis = by_global_analysis,
        global_analysis = global_analysis,
    )
end

function analyze_environment(
    context::Context,
    n_blocks_in_neighborhood::Integer,
    next_n_blocks::Integer,
    ordered_block_indices::Vector{<:Integer},
)::Analysis
    blocks = context.blocks
    @assert blocks !== nothing

    environment_metacells_mask = region_metacells_mask(context, next_n_blocks, ordered_block_indices)

    reset_included!(context, environment_metacells_mask; keep_predictive = true)

    analyze_left_outs!(context)
    analysis = final_analysis!(context)

    print(
        stderr,
        "TODOX $(context.included_metacells.n_entries) metacells in $(n_blocks_in_neighborhood) <= $(next_n_blocks) => RMS: $(analysis.cross_validation.mean_rms) ...        \r",
    )

    return analysis
end

function region_metacells_mask(
    context::Context,
    n_blocks::Integer,
    ordered_block_indices::Vector{<:Integer},
)::Union{Vector{Bool}, BitVector}
    blocks = context.blocks
    @assert blocks !== nothing

    metacells_mask = zeros(Bool, context.data.n_metacells)
    for block_index in ordered_block_indices[1:n_blocks]
        metacells_mask[blocks.metacells_of_blocks[block_index]] .= true
    end

    return metacells_mask
end

end  # module
