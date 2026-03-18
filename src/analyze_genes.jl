"""
Do simple per-gene analysis.
"""
module AnalyzeGenes

export compute_vector_of_is_correlated_with_skeleton_per_gene!
export compute_vector_of_is_marker_per_gene!
export compute_vector_of_is_skeleton_per_gene!
export compute_vector_of_marker_rank_per_gene!
export fetch_gmara_vector_of_is_regulator_per_gene!
export fetch_gmara_vector_of_is_transcription_factor_per_gene!
export rank_variables

using DataAxesFormats
using LinearAlgebra
using TanayLabUtilities
using StatsBase
using Random

using ..AnalyzeMetacells
using ..Contracts
using ..Defaults
using ..Gmara

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_is_correlated_with_skeleton_per_gene
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_forbidden_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_regulator_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_is_transcription_factor_per_gene
import Metacells.Contracts.vector_of_marker_rank_per_gene

"""
    function compute_vector_of_is_marker_per_gene!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set `vector_of_is_marker_per_gene`. Such genes:

  - Have a maximal gene expression level of at least `min_marker_gene_max_fraction`, and
  - Have a range of expression of at least `min_marker_gene_range_fold` between the maximal and minimal
    (log base 2) of the expression level. For the minimal expression, we take the `min_marker_quantile`
    of the expression levels in the population.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_marker_per_gene(CreatedOutput),
    ],
) function compute_vector_of_is_marker_per_gene!(  # UNTESTED
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 1e-4,
    min_marker_gene_range_fold::Real = 2,
    min_marker_quantile::Real = 0,
    overwrite::Bool = false,
)::Nothing
    @assert min_marker_gene_range_fold >= 0
    @assert 0 <= min_marker_gene_max_fraction <= 1
    @assert 0 <= min_marker_quantile < 1

    maximal_fraction_per_gene = daf["@ metacell @ gene :: linear_fraction >- Max"].array
    if min_marker_quantile == 0
        log_minimal_fraction_per_gene = daf["@ metacell @ gene :: log_linear_fraction >- Min"].array
    else
        log_minimal_fraction_per_gene =
            daf["@ metacell @ gene :: log_linear_fraction >- Quantile p $(min_marker_quantile)"].array
    end
    log_maximal_fraction_per_gene = daf["@ metacell @ gene :: log_linear_fraction >- Max"].array
    range_fold_per_gene = log_maximal_fraction_per_gene .- log_minimal_fraction_per_gene

    is_marker_per_gene =
        (maximal_fraction_per_gene .>= min_marker_gene_max_fraction) .&
        (range_fold_per_gene .>= min_marker_gene_range_fold)

    set_vector!(daf, "gene", "is_marker", is_marker_per_gene; overwrite)
    @debug "Marker genes: $(sum(is_marker_per_gene))" _group = :mcs_results

    return nothing
end

"""
    function compute_vector_of_is_correlated_with_skeleton_per_gene!(
        daf::DafWriter;
        min_gene_correlation::AbstractFloat = $(DEFAULT.min_gene_correlation),
        min_gene_correlation_quantile::AbstractFloat = $(DEFAULT.min_gene_correlation_quantile),
        genes_correlation_window::Integer = $(DEFAULT.genes_correlation_window),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set `vector_of_is_correlated_with_skeleton_per_gene`. Such genes:

  - Are markers (see [`vector_of_is_marker_per_gene`](@ref).
  - Have a correlation of at least `min_gene_correlation` with some skeleton gene (see
    [`vector_of_is_skeleton_per_gene`](@ref)), OR
  - Have a correlation with a skeleton gene which is in the `min_gene_correlation_quantile` out of the
    `genes_correlation_window` genes with a similar correlation with some skeleton.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_is_correlated_with_skeleton_per_gene(CreatedOutput),
    ],
) function compute_vector_of_is_correlated_with_skeleton_per_gene!(  # UNTESTED
    daf::DafWriter;
    min_gene_correlation::AbstractFloat = 0.60,
    min_gene_correlation_quantile::AbstractFloat = 0.90,
    genes_correlation_window::Integer = 100,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= min_gene_correlation_quantile < 1
    @assert genes_correlation_window > 0

    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array
    indices_of_skeletons = findall(is_skeleton_per_gene)
    n_skeletons = length(indices_of_skeletons)

    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)
    n_markers = length(indices_of_markers)

    log_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "log_linear_fraction").array
    @views log_fraction_per_metacell_per_marker = log_fraction_per_metacell_per_gene[:, indices_of_markers]
    @views log_fraction_per_metacell_per_skeleton = log_fraction_per_metacell_per_gene[:, indices_of_skeletons]

    is_correlated_with_skeleton_per_gene = zeros(Bool, n_genes)
    correlation_per_skeleton_per_marker = Matrix{Float32}(undef, n_skeletons, n_markers)
    max_correlation_per_marker = Vector{Float32}(undef, n_markers)
    scratch_per_metacell_per_skeleton = Matrix{Float32}(undef, n_metacells, n_skeletons)
    scratch_per_metacell_per_marker = Matrix{Float32}(undef, n_metacells, n_markers)
    scratch_per_skeleton = Vector{Float32}(undef, n_skeletons)
    scratch_per_marker = Vector{Float32}(undef, n_markers)
    scratch_marker_window = Vector{Float32}(undef, genes_correlation_window)

    fill_vector_of_is_correlated_with_skeleton_per_gene!(;
        min_gene_correlation,
        min_gene_correlation_quantile,
        genes_correlation_window,
        indices_of_markers,
        log_fraction_per_metacell_per_skeleton,
        log_fraction_per_metacell_per_marker,
        is_correlated_with_skeleton_per_gene,
        correlation_per_skeleton_per_marker,
        max_correlation_per_marker,
        scratch_per_metacell_per_skeleton,
        scratch_per_metacell_per_marker,
        scratch_per_skeleton,
        scratch_per_marker,
        scratch_marker_window,
    )

    set_vector!(daf, "gene", "is_correlated_with_skeleton", is_correlated_with_skeleton_per_gene; overwrite)
    @debug "Markers correlated with skeleton genes: $(sum(is_correlated_with_skeleton_per_gene))" _group = :mcs_results

    return nothing
end

function fill_vector_of_is_correlated_with_skeleton_per_gene!(;  # NOJET
    min_gene_correlation::AbstractFloat,
    min_gene_correlation_quantile::AbstractFloat,
    genes_correlation_window::Integer,
    indices_of_markers::AbstractVector{<:Integer},
    log_fraction_per_metacell_per_skeleton::AbstractMatrix{<:AbstractFloat},
    log_fraction_per_metacell_per_marker::AbstractMatrix{<:AbstractFloat},
    is_correlated_with_skeleton_per_gene::AbstractVector{Bool},
    correlation_per_skeleton_per_marker::AbstractMatrix{<:AbstractFloat},
    max_correlation_per_marker::AbstractVector{<:AbstractFloat},
    scratch_per_metacell_per_skeleton::AbstractMatrix{<:AbstractFloat},
    scratch_per_metacell_per_marker::AbstractMatrix{<:AbstractFloat},
    scratch_per_skeleton::AbstractVector{<:AbstractFloat},
    scratch_per_marker::AbstractVector{<:AbstractFloat},
    scratch_marker_window::AbstractVector{<:AbstractFloat},
)::Nothing
    @assert length(scratch_marker_window) >= genes_correlation_window

    n_markers = size(log_fraction_per_metacell_per_marker, 2)

    zero_cor_between_matrices_columns(
        log_fraction_per_metacell_per_skeleton,
        log_fraction_per_metacell_per_marker;
        result = correlation_per_skeleton_per_marker,
        left_scratch_matrix = scratch_per_metacell_per_skeleton,
        right_scratch_matrix = scratch_per_metacell_per_marker,
        left_scratch_row = scratch_per_skeleton,
        right_scratch_row = scratch_per_marker,
    )

    for marker_position in 1:n_markers
        max_correlation_per_marker[marker_position] =
            maximum(@view correlation_per_skeleton_per_marker[:, marker_position])
    end
    @assert_vector(max_correlation_per_marker, n_markers)

    quantile_correlation_per_marker = scratch_per_marker
    rolling_quantile!(;
        results = quantile_correlation_per_marker,
        data = max_correlation_per_marker,
        scratch = scratch_marker_window,
        window_size = genes_correlation_window,
        quantile_fraction = min_gene_correlation_quantile,
    )

    for (marker_position, marker_index) in enumerate(indices_of_markers)
        max_correlation = max_correlation_per_marker[marker_position]
        is_correlated_with_skeleton_per_gene[marker_index] =
            ((max_correlation >= quantile_correlation_per_marker[marker_position]) & (max_correlation > 0)) |
            (max_correlation >= min_gene_correlation)
    end

    return nothing
end

# TODO: Move to `TanayLabUtilities`?
function non_allocating_quantile(;
    data::AbstractVector{T},
    scratch::AbstractVector{T},
    quantile_fraction::AbstractFloat,
)::T where {T}
    n_values = length(data)
    @assert length(scratch) == n_values
    scratch .= data

    float_index = quantile_fraction * (n_values - 1)
    k_low = clamp(floor(Int, float_index) + 1, 1, n_values)
    k_high = clamp(k_low + 1, 1, n_values)

    if k_low == k_high
        partialsort!(scratch, k_low)
        return scratch[k_low]
    end

    partialsort!(scratch, k_low)
    low_value = scratch[k_low]
    @views high_values = scratch[(k_low + 1):end]
    high_value = minimum(high_values)

    weight = float_index - (k_low - 1)
    return T(low_value + weight * (high_value - low_value))
end

# TODO: Move to `TanayLabUtilities`?
function rolling_quantile!(;
    results::AbstractVector{T},
    data::AbstractVector{T},
    scratch::AbstractVector{T},
    window_size::Int,
    quantile_fraction::AbstractFloat,
)::Nothing where {T}
    @assert window_size > 1
    @assert length(results) == length(data)
    @assert length(scratch) >= window_size

    n_values = length(data)
    window_left_of_center = div(window_size - 1, 2)
    window_right_of_center = window_size - window_left_of_center - 1

    if n_values <= window_size
        @views scratch_view = scratch[1:n_values]
        fill!(results, non_allocating_quantile(; data, scratch = scratch_view, quantile_fraction))
        return nothing
    end

    @views scratch_view = scratch[1:window_size]

    for index in 1:(n_values - window_size + 1)
        data_view = @view(data[index:(index + window_size - 1)])
        results[window_left_of_center + index] =
            non_allocating_quantile(; data = data_view, scratch = scratch_view, quantile_fraction)
    end

    first_result = results[window_left_of_center + 1]
    results[1:window_left_of_center] .= first_result

    last_result = results[n_values - window_right_of_center]
    results[(end - window_right_of_center + 1):end] .= last_result

    return nothing
end

"""
    function compute_vector_of_marker_rank_per_gene!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_marker_rank_per_gene`](@ref). To compute this:

 1. Non-marker genes are given a rank of `typemax(UInt32)` regardless of their expression level in the metacells (see
    [`vector_of_is_marker_per_gene`](@ref)).
 2. Compute the median of the (log base 2) of the gene expression level for each marker gene across all the metacells.
 3. Compute the per-marker-per-metacell fold factor (absolute difference of the log expression from the median).
 4. Invoke [`rank_variables`](@ref) to rank the marker genes from most to least significant.

Non-marker genes are given a rank of `typemax(UInt32)` regardless of their expression level in the metacells.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_marker_rank_per_gene(CreatedOutput),
    ],
) function compute_vector_of_marker_rank_per_gene!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_genes = axis_length(daf, "gene")

    is_marker_per_gene = get_vector(daf, "gene", "is_marker")
    n_markers = sum(is_marker_per_gene)

    log_fraction_per_metacell_per_marker = daf["@ metacell @ gene [ is_marker ] :: log_linear_fraction"].array
    @assert_matrix(log_fraction_per_metacell_per_marker, n_metacells, n_markers, Columns)

    median_log_fraction_per_marker = daf["@ metacell @ gene [ is_marker ] :: log_linear_fraction >- Median"].array
    @assert_vector(median_log_fraction_per_marker, n_markers)

    abs_fold_per_metacell_per_marker =  # NOJET
        abs.(log_fraction_per_metacell_per_marker .- transpose(median_log_fraction_per_marker))
    abs_fold_per_marker_per_metacell = flipped(abs_fold_per_metacell_per_marker)

    rank_per_marker = rank_variables(abs_fold_per_marker_per_metacell)

    marker_rank_per_gene = fill(typemax(UInt32), n_genes)
    marker_rank_per_gene[is_marker_per_gene] .= rank_per_marker

    set_vector!(daf, "gene", "marker_rank", marker_rank_per_gene; overwrite)
    return nothing
end

"""
    rank_variables(score_per_variable_per_observation::AbstractMatrix{<:AbstractFloat})::AbstractVector{<:Integer}

Given some `score_per_variable_per_observation` matrix, return a vector of the rank of each variable such that the
"most significant" variables are first.

 1. Rank the variables for each observation (1 having the highest score).
 2. For each variable, give it a priority which is a tuple of (1) the minimal rank it has in all observations (2) the maximal
    score it has in observations where it has that rank (negated).
 3. Sort the variables according to this priority.

This heuristic is useful for focusing on the "most significant" variables in a data set. It is used by
[`compute_vector_of_marker_rank_per_gene!`](@ref).

# TODO: Move to `TanayLabUtilities`?
"""
function rank_variables(score_per_variable_per_observation::AbstractMatrix{<:AbstractFloat})::AbstractVector{<:Integer}
    @assert_matrix(score_per_variable_per_observation, Columns)
    n_variables, n_observations = size(score_per_variable_per_observation)

    @views rank_per_variable_per_observation = hcat(
        [
            invperm(sortperm(vec(score_per_variable_per_observation[:, observation_index]); rev = true)) for
            observation_index in 1:n_observations
        ]...,
    )
    @assert_matrix(rank_per_variable_per_observation, n_variables, n_observations, Columns)

    min_rank_per_variable = vec(minimum(rank_per_variable_per_observation; dims = 2))
    @assert_vector(min_rank_per_variable, n_variables)

    # TODO: This can be done in parallel.
    maximal_score_per_variable = Vector{Float32}(undef, n_variables)
    for variable_index in 1:n_variables
        @views maximal_score_per_observation = score_per_variable_per_observation[
            variable_index,
            rank_per_variable_per_observation[variable_index, :] .== min_rank_per_variable[variable_index],
        ]
        maximal_score_per_variable[variable_index] = mean(maximal_score_per_observation)  # NOLINT
    end

    priority_per_variable = collect(zip(min_rank_per_variable, -maximal_score_per_variable))
    rank_per_variable = invperm(sortperm(priority_per_variable))
    @assert_vector(rank_per_variable, n_variables)
    return rank_per_variable
end

"""
    fetch_gmara_transcription_factors!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString = $(DEFAULT.namespace),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = $(DEFAULT.cache_dir),
        timeout::Real = $(DEFAULT.timeout),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Fetch and set [`vector_of_is_transcription_factor_per_gene`](@ref) from `Gmara`.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [vector_of_is_transcription_factor_per_gene(CreatedOutput)],
) function fetch_gmara_vector_of_is_transcription_factor_per_gene!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    version::AbstractString = function_default(set_gmara_genes_mask!, :version),
    cache_dir = function_default(set_gmara_genes_mask!, :cache_dir),
    timeout::Real = function_default(set_gmara_genes_mask!, :timeout),
    overwrite::Bool = false,
)::Nothing
    n_transcription_factors = set_gmara_genes_mask!(
        daf;
        species,
        namespace,
        list = "transcription_factor",
        version,
        cache_dir,
        timeout,
        overwrite,
    )
    @debug "Transcription factors: $(n_transcription_factors)" _group = :mcs_results
    return nothing
end

"""
    fetch_gmara_vector_of_is_regulator_per_gene!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString = $(DEFAULT.namespace),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = $(DEFAULT.cache_dir),
        timeout::Real = $(DEFAULT.timeout),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Fetch and set [`vector_of_is_transcription_factor_per_gene`](@ref) from `Gmara`.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [vector_of_is_regulator_per_gene(CreatedOutput)],
) function fetch_gmara_vector_of_is_regulator_per_gene!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    version::AbstractString = function_default(set_gmara_genes_mask!, :version),
    cache_dir = function_default(set_gmara_genes_mask!, :cache_dir),
    timeout::Real = function_default(set_gmara_genes_mask!, :timeout),
    overwrite::Bool = false,
)::Nothing
    n_regulators =
        set_gmara_genes_mask!(daf; species, namespace, list = "regulator", version, cache_dir, timeout, overwrite)
    @debug "Regulators: $(n_regulators)" _group = :mcs_results
    return nothing
end

"""
    function compute_is_skeleton_per_gene_vector!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_is_skeleton_per_gene`](@ref). We just pick the genes that are regulators, and are not
excluded, lateral or forbidden.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_forbidden_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(CreatedOutput),
    ],
) function compute_vector_of_is_skeleton_per_gene!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_regulator_per_gene = get_vector(daf, "gene", "is_regulator").array
    is_excluded_per_gene = get_vector(daf, "gene", "is_excluded").array
    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    is_forbidden_per_gene = get_vector(daf, "gene", "is_forbidden").array

    is_skeleton_per_gene =
        is_regulator_per_gene .& .!is_excluded_per_gene .& .!is_lateral_per_gene .& .!is_forbidden_per_gene .& is_marker_per_gene

    set_vector!(daf, "gene", "is_skeleton", is_skeleton_per_gene; overwrite)

    name_per_gene = axis_vector(daf, "gene")
    @debug "Skeletons: $(length(is_skeleton_per_gene)) [ $(join(sort(name_per_gene[is_skeleton_per_gene]), ", ")) ]" _group =
        :mcs_results

    @assert any(is_skeleton_per_gene)
    return nothing
end

end  # module
