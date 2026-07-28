"""
Do simple per-gene analysis.
"""
module AnalyzeGenes

export compute_gene_report
export compute_skeleton_report
export compute_vector_of_is_correlated_with_skeleton_per_gene!
export compute_vector_of_is_marker_per_gene!
export compute_vector_of_is_skeleton_per_gene!
export compute_vector_of_marker_rank_per_gene!
export fetch_gmara_vector_of_is_regulator_per_gene!
export fetch_gmara_vector_of_is_transcription_factor_per_gene!
export rank_variables

using DataAxesFormats
using DataFrames
using LinearAlgebra
using LoopVectorization
using TanayLabUtilities
using StatsBase
using Random

using ..AnalyzeMetacells
using ..Contracts
using ..Defaults
using ..Gmara

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.base_block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_projected_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
import Metacells.Contracts.matrix_of_correlation_between_markers_per_gene_per_gene
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
    axes = [gene_axis(RequiredInput)],
    data = [
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
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

    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array
    indices_of_skeletons = findall(is_skeleton_per_gene)

    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)
    n_markers = length(indices_of_markers)

    correlation_per_gene_per_gene = get_matrix(daf, "gene", "gene", "correlation_between_markers").array
    correlation_per_skeleton_per_marker =
        Matrix{Float32}(correlation_per_gene_per_gene[indices_of_skeletons, indices_of_markers])

    is_correlated_with_skeleton_per_gene = zeros(Bool, n_genes)
    max_correlation_per_marker = Vector{Float32}(undef, n_markers)
    scratch_per_marker = Vector{Float32}(undef, n_markers)
    scratch_marker_window = Vector{Float32}(undef, genes_correlation_window)

    mark_is_correlated_with_skeleton_per_gene!(;
        min_gene_correlation,
        min_gene_correlation_quantile,
        genes_correlation_window,
        indices_of_markers,
        correlation_per_skeleton_per_marker,
        is_correlated_with_skeleton_per_gene,
        max_correlation_per_marker,
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
    zero_cor_between_matrices_columns(
        log_fraction_per_metacell_per_skeleton,
        log_fraction_per_metacell_per_marker;
        result = correlation_per_skeleton_per_marker,
        left_scratch_matrix = scratch_per_metacell_per_skeleton,
        right_scratch_matrix = scratch_per_metacell_per_marker,
        left_scratch_row = scratch_per_skeleton,
        right_scratch_row = scratch_per_marker,
    )
    mark_is_correlated_with_skeleton_per_gene!(;
        min_gene_correlation,
        min_gene_correlation_quantile,
        genes_correlation_window,
        indices_of_markers,
        correlation_per_skeleton_per_marker,
        is_correlated_with_skeleton_per_gene,
        max_correlation_per_marker,
        scratch_per_marker,
        scratch_marker_window,
    )
    return nothing
end

# Given the correlation of each marker with each skeleton, mark the markers correlated with the skeleton: those with a
# high enough maximal skeleton correlation, either absolutely (`min_gene_correlation`) or relative to the
# `min_gene_correlation_quantile` of a rolling window of `genes_correlation_window` markers with a similar maximal
# skeleton correlation.
function mark_is_correlated_with_skeleton_per_gene!(;  # NOJET
    min_gene_correlation::AbstractFloat,
    min_gene_correlation_quantile::AbstractFloat,
    genes_correlation_window::Integer,
    indices_of_markers::AbstractVector{<:Integer},
    correlation_per_skeleton_per_marker::AbstractMatrix{<:AbstractFloat},
    is_correlated_with_skeleton_per_gene::AbstractVector{Bool},
    max_correlation_per_marker::AbstractVector{<:AbstractFloat},
    scratch_per_marker::AbstractVector{<:AbstractFloat},
    scratch_marker_window::AbstractVector{<:AbstractFloat},
)::Nothing
    @assert length(scratch_marker_window) >= genes_correlation_window
    n_markers = size(correlation_per_skeleton_per_marker, 2)

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

    log_fraction_per_metacell_per_marker =
        mutable_array(densify(daf["@ metacell @ gene [ is_marker ] :: log_linear_fraction"].array))
    @assert_matrix(log_fraction_per_metacell_per_marker, n_metacells, n_markers, Columns)

    median_log_fraction_per_marker =
        mutable_array(densify(daf["@ metacell @ gene [ is_marker ] :: log_linear_fraction >- Median"].array))
    @assert_vector(median_log_fraction_per_marker, n_markers)

    abs_fold_per_metacell_per_marker =  # NOJET
        Matrix{Float32}(undef, size(log_fraction_per_metacell_per_marker))
    @check_turbo_matrix(abs_fold_per_metacell_per_marker)
    @check_turbo_matrix(log_fraction_per_metacell_per_marker)
    @check_turbo_vector(median_log_fraction_per_marker)
    n_metacells, n_markers = size(log_fraction_per_metacell_per_marker)
    parallel_loop_wo_rng(
        1:n_markers;
        name = "abs_fold_per_metacell_per_marker",
        progress = DebugProgress(n_markers; group = :mcs_loops, desc = "abs_fold_per_metacell_per_marker"),
    ) do marker_index
        @turbo for metacell_index in 1:n_metacells
            abs_fold_per_metacell_per_marker[metacell_index, marker_index] = abs(
                log_fraction_per_metacell_per_marker[metacell_index, marker_index] -
                median_log_fraction_per_marker[marker_index],
            )
        end
        return nothing
    end
    abs_fold_per_marker_per_metacell = flipped(abs_fold_per_metacell_per_marker)

    rank_per_marker = rank_variables(abs_fold_per_marker_per_metacell)  # NOJET

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

    rank_per_variable_per_observation = Matrix{Int}(undef, n_variables, n_observations)
    parallel_loop_wo_rng(
        1:n_observations;
        name = "rank_per_variable_per_observation",
        progress = DebugProgress(n_observations; group = :mcs_loops, desc = "rank_per_variable_per_observation"),
    ) do observation_index
        @views rank_per_variable_per_observation[:, observation_index] .=
            invperm(sortperm(vec(score_per_variable_per_observation[:, observation_index]); rev = true))
        return nothing
    end
    @assert_matrix(rank_per_variable_per_observation, n_variables, n_observations, Columns)

    min_rank_per_variable = Vector{Int}(undef, n_variables)
    maximal_score_per_variable = Vector{Float32}(undef, n_variables)
    parallel_loop_wo_rng(
        1:n_variables;
        name = "min_rank_and_maximal_score_per_variable",
        progress = DebugProgress(n_variables; group = :mcs_loops, desc = "min_rank_and_maximal_score_per_variable"),
    ) do variable_index
        @views min_rank_per_variable[variable_index] = minimum(rank_per_variable_per_observation[variable_index, :])
        @views maximal_score_per_observation = score_per_variable_per_observation[
            variable_index,
            rank_per_variable_per_observation[variable_index, :] .== min_rank_per_variable[variable_index],
        ]
        maximal_score_per_variable[variable_index] = mean(maximal_score_per_observation)  # NOLINT
        return nothing
    end
    @assert_vector(min_rank_per_variable, n_variables)

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

"""
    compute_gene_report(;
        daf::DafReader,
        base_daf::Maybe{DafReader} = nothing,
    )::DataFrame

Return a per-marker-gene report as a `DataFrame`, one row per marker gene, sorted by `marker_rank`. Its columns are:

  - `gene` - the gene name.
  - `marker_rank` - the gene's rank as a marker (`1` being the most significant; see
    [`vector_of_marker_rank_per_gene`](@ref)).
  - `is_lateral`, `is_tf`, `is_regulator`, `is_skeleton` - whether the gene is lateral, a transcription factor, a
    regulator, or a skeleton gene.
  - `top_marker` and `top_marker_cor`, `top_lateral` and `top_lateral_cor`, `top_regulator` and `top_regulator_cor`,
    `top_skeleton` and `top_skeleton_cor` - the most correlated marker / lateral / regulator / skeleton gene (other than
    the gene itself) and that correlation, from [`matrix_of_correlation_between_markers_per_gene_per_gene`](@ref). Since
    only markers are correlated, e.g. `top_lateral` is the most correlated lateral *marker*. An empty set (e.g. no
    regulators) yields an empty name and a zero correlation.

When `base_daf` (the round-0 model) is given, six more columns compare each gene's base-neighborhood correlation in `daf`
against `base_daf`, over the base blocks where the gene is scored (has a non-zero correlation) at round 0, for the
punctuated (`punct_*`) and projected-punctuated (`proj_punct_*`) correlations:

  - `punct_Δ`, `proj_punct_Δ` - the mean change (`Δ`) in the correlation from `base_daf`.
  - `punct_++`, `proj_punct_++` - the fraction of the base blocks where the correlation significantly improved
    (increased by at least `0.05`).
  - `punct_--`, `proj_punct_--` - the fraction of the base blocks where the correlation significantly degraded
    (decreased by at least `0.05`).

These six columns are omitted when `base_daf` is `nothing` (that is, when `daf` is itself the round-0 model).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    name = "daf",
    axes = [gene_axis(RequiredInput), base_block_axis(OptionalInput)],
    data = [
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_marker_rank_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_transcription_factor_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            OptionalInput,
        ),
        matrix_of_correlation_between_base_neighborhood_cells_and_projected_punctuated_metacells_per_gene_per_base_block(
            OptionalInput,
        ),
    ],
) function compute_gene_report(; daf::DafReader, base_daf::Maybe{DafReader} = nothing)::DataFrame
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)

    gene_name_per_gene = axis_vector(daf, "gene")
    gene_name_per_marker = gene_name_per_gene[indices_of_markers]
    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_transcription_factor_per_gene = get_vector(daf, "gene", "is_transcription_factor").array
    is_regulator_per_gene = get_vector(daf, "gene", "is_regulator").array
    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array

    correlation_per_gene_per_gene = get_matrix(daf, "gene", "gene", "correlation_between_markers").array
    correlation_per_marker_per_marker =
        Matrix{Float32}(correlation_per_gene_per_gene[indices_of_markers, indices_of_markers])

    top_marker, top_marker_cor = most_correlated_gene_in_set_per_marker(
        correlation_per_marker_per_marker,
        gene_name_per_marker,
        indices_of_markers,
        is_marker_per_gene,
    )
    top_lateral, top_lateral_cor = most_correlated_gene_in_set_per_marker(
        correlation_per_marker_per_marker,
        gene_name_per_marker,
        indices_of_markers,
        is_lateral_per_gene,
    )
    top_regulator, top_regulator_cor = most_correlated_gene_in_set_per_marker(
        correlation_per_marker_per_marker,
        gene_name_per_marker,
        indices_of_markers,
        is_regulator_per_gene,
    )
    top_skeleton, top_skeleton_cor = most_correlated_gene_in_set_per_marker(
        correlation_per_marker_per_marker,
        gene_name_per_marker,
        indices_of_markers,
        is_skeleton_per_gene,
    )

    data_frame = DataFrame(;
        gene = gene_name_per_marker,
        marker_rank = get_vector(daf, "gene", "marker_rank").array[indices_of_markers],
        is_lateral = is_lateral_per_gene[indices_of_markers],
        is_tf = is_transcription_factor_per_gene[indices_of_markers],
        is_regulator = is_regulator_per_gene[indices_of_markers],
        is_skeleton = is_skeleton_per_gene[indices_of_markers],
        top_marker,
        top_marker_cor,
        top_lateral,
        top_lateral_cor,
        top_regulator,
        top_regulator_cor,
        top_skeleton,
        top_skeleton_cor,
    )

    if base_daf !== nothing
        punct_delta, punct_better, punct_worse = correlation_change_per_marker(
            daf,
            base_daf,
            "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
            indices_of_markers,
        )
        data_frame[!, "punct_Δ"] = punct_delta
        data_frame[!, "punct_++"] = punct_better
        data_frame[!, "punct_--"] = punct_worse
        proj_punct_delta, proj_punct_better, proj_punct_worse = correlation_change_per_marker(
            daf,
            base_daf,
            "correlation_between_base_neighborhood_cells_and_projected_punctuated_metacells",
            indices_of_markers,
        )
        data_frame[!, "proj_punct_Δ"] = proj_punct_delta
        data_frame[!, "proj_punct_++"] = proj_punct_better
        data_frame[!, "proj_punct_--"] = proj_punct_worse
    end

    sort!(data_frame, :marker_rank)
    return data_frame
end

# For each marker gene, the gene in the set (restricted to markers, since only markers are correlated) most correlated
# with it and that correlation, excluding the gene itself. An empty set yields an empty name and a zero correlation.
function most_correlated_gene_in_set_per_marker(
    correlation_per_marker_per_marker::AbstractMatrix{<:Real},
    gene_name_per_marker::AbstractVector{<:AbstractString},
    indices_of_markers::AbstractVector{<:Integer},
    is_in_set_per_gene::AbstractVector{Bool},
)::Tuple{Vector{AbstractString}, Vector{Float64}}
    n_markers = length(indices_of_markers)
    set_marker_positions =
        [position for (position, gene_index) in enumerate(indices_of_markers) if is_in_set_per_gene[gene_index]]
    top_gene_per_marker = Vector{AbstractString}(undef, n_markers)
    top_correlation_per_marker = Vector{Float64}(undef, n_markers)
    parallel_loop_wo_rng(1:n_markers; name = "most_correlated_gene_in_set_per_marker") do marker_position
        best_position = 0
        best_correlation = -Inf
        for set_position in set_marker_positions
            if set_position != marker_position
                correlation = correlation_per_marker_per_marker[marker_position, set_position]
                if correlation > best_correlation
                    best_correlation = correlation
                    best_position = set_position
                end
            end
        end
        if best_position == 0
            top_gene_per_marker[marker_position] = ""
            top_correlation_per_marker[marker_position] = 0.0
        else
            top_gene_per_marker[marker_position] = gene_name_per_marker[best_position]
            top_correlation_per_marker[marker_position] = best_correlation
        end
        return nothing
    end
    return (top_gene_per_marker, top_correlation_per_marker)
end

# For a base-neighborhood correlation matrix (gene X base_block), the per-marker mean change from the round-0 `base_daf`
# and the fraction of the base blocks in which the gene significantly improved (by at least `0.05`) or degraded, over the
# base blocks where the gene is scored (has a non-zero correlation) at round 0.
function correlation_change_per_marker(
    daf::DafReader,
    base_daf::DafReader,
    property_name::AbstractString,
    indices_of_markers::AbstractVector{<:Integer},
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    correlation_per_gene_per_base_block = get_matrix(daf, "gene", "base_block", property_name).array
    base_correlation_per_gene_per_base_block = get_matrix(base_daf, "gene", "base_block", property_name).array
    @assert axis_vector(daf, "base_block") == axis_vector(base_daf, "base_block")
    n_base_blocks = size(base_correlation_per_gene_per_base_block, 2)

    n_markers = length(indices_of_markers)
    delta_per_marker = zeros(Float64, n_markers)
    better_per_marker = zeros(Float64, n_markers)
    worse_per_marker = zeros(Float64, n_markers)
    for (marker_position, gene_index) in enumerate(indices_of_markers)
        @views correlation_per_base_block = correlation_per_gene_per_base_block[gene_index, :]
        @views base_correlation_per_base_block = base_correlation_per_gene_per_base_block[gene_index, :]
        n_relevant = 0
        sum_delta = 0.0
        n_better = 0
        n_worse = 0
        for base_block in 1:n_base_blocks
            base_correlation = base_correlation_per_base_block[base_block]
            if base_correlation != 0
                n_relevant += 1
                delta = correlation_per_base_block[base_block] - base_correlation
                sum_delta += delta
                if delta >= 0.05
                    n_better += 1
                elseif delta <= -0.05
                    n_worse += 1
                end
            end
        end
        if n_relevant > 0
            delta_per_marker[marker_position] = sum_delta / n_relevant
            better_per_marker[marker_position] = n_better / n_relevant
            worse_per_marker[marker_position] = n_worse / n_relevant
        end
    end
    return (delta_per_marker, better_per_marker, worse_per_marker)
end

"""
    compute_skeleton_report(;
        daf::DafReader,
        base_daf::Maybe{DafReader} = nothing,
    )::DataFrame

Return a per-skeleton-gene report as a `DataFrame`, one row per skeleton gene, sorted by the total number of markers most
correlated with it (descending). Each marker gene is assigned to the single skeleton gene it is most correlated with
(from [`matrix_of_correlation_between_markers_per_gene_per_gene`](@ref), excluding itself), and split into lateral
(`lat`) markers and pertinent (`pert`, non-lateral) markers. Apart from the `gene` (skeleton name) column, every
statistic is emitted as a `lat_<name>`, `lat_<name>_f`, `pert_<name>` trio - the lateral value, its fraction of the
lateral-plus-pertinent total, and the pertinent value. The statistics are:

  - `#` - the number of markers most correlated with this skeleton.

When `base_daf` (the round-0 model) is given, four more statistics per each of the punctuated (`punct`) and
projected-punctuated (`proj_punct`) base-neighborhood correlations aggregate the [`compute_gene_report`](@ref) values of
those markers - where each marker's delta correlation (`Δ`) is the mean over its base blocks (with a non-zero round-0
correlation) of `round_N - round_0`, and `++` / `--` are the fractions of its base blocks that significantly improved /
degraded (by at least `0.05`):

  - `<corr>_++` - the mean (`μ`) of the markers' improved fractions.
  - `<corr>_--` - the mean of the markers' degraded fractions.
  - `<corr>_μΔ` - the mean of the markers' delta correlations.
  - `<corr>_ΣΔ` - the sum (`Σ`) of the markers' delta correlations.

These comparison statistics are omitted when `base_daf` is `nothing`.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    name = "daf",
    axes = [gene_axis(RequiredInput), base_block_axis(OptionalInput)],
    data = [
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_skeleton_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
        matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
            OptionalInput,
        ),
        matrix_of_correlation_between_base_neighborhood_cells_and_projected_punctuated_metacells_per_gene_per_base_block(
            OptionalInput,
        ),
    ],
) function compute_skeleton_report(; daf::DafReader, base_daf::Maybe{DafReader} = nothing)::DataFrame
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    indices_of_markers = findall(is_marker_per_gene)
    n_markers = length(indices_of_markers)

    is_skeleton_per_gene = get_vector(daf, "gene", "is_skeleton").array
    indices_of_skeletons = findall(is_skeleton_per_gene)
    n_skeletons = length(indices_of_skeletons)

    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    gene_name_per_gene = axis_vector(daf, "gene")

    correlation_per_gene_per_gene = get_matrix(daf, "gene", "gene", "correlation_between_markers").array
    correlation_per_marker_per_skeleton =
        Matrix{Float32}(correlation_per_gene_per_gene[indices_of_markers, indices_of_skeletons])

    # Assign each marker to the skeleton it is most correlated with (excluding itself when it is a skeleton), and note
    # whether it is lateral.
    top_skeleton_per_marker = Vector{Int}(undef, n_markers)
    is_lateral_per_marker = Vector{Bool}(undef, n_markers)
    for marker_position in 1:n_markers
        marker_gene_index = indices_of_markers[marker_position]
        is_lateral_per_marker[marker_position] = is_lateral_per_gene[marker_gene_index]
        best_skeleton_position = 0
        best_correlation = -Inf
        for skeleton_position in 1:n_skeletons
            if indices_of_skeletons[skeleton_position] == marker_gene_index
                continue
            end
            correlation = correlation_per_marker_per_skeleton[marker_position, skeleton_position]
            if correlation > best_correlation
                best_correlation = correlation
                best_skeleton_position = skeleton_position
            end
        end
        top_skeleton_per_marker[marker_position] = best_skeleton_position
    end

    change_per_marker_per_correlation = nothing
    if base_daf !== nothing
        punct_delta, punct_imp, punct_deg = correlation_change_per_marker(
            daf,
            base_daf,
            "correlation_between_base_neighborhood_cells_and_punctuated_metacells",
            indices_of_markers,
        )
        proj_punct_delta, proj_punct_imp, proj_punct_deg = correlation_change_per_marker(
            daf,
            base_daf,
            "correlation_between_base_neighborhood_cells_and_projected_punctuated_metacells",
            indices_of_markers,
        )
        change_per_marker_per_correlation = [
            ("punct", punct_delta, punct_imp, punct_deg),
            ("proj_punct", proj_punct_delta, proj_punct_imp, proj_punct_deg),
        ]
    end

    lateral_count_per_skeleton = zeros(Int, n_skeletons)
    pertinent_count_per_skeleton = zeros(Int, n_skeletons)
    for marker_position in 1:n_markers
        skeleton = top_skeleton_per_marker[marker_position]
        if is_lateral_per_marker[marker_position]
            lateral_count_per_skeleton[skeleton] += 1
        else
            pertinent_count_per_skeleton[skeleton] += 1
        end
    end
    total_count_per_skeleton = lateral_count_per_skeleton .+ pertinent_count_per_skeleton

    data_frame = DataFrame(; gene = gene_name_per_gene[indices_of_skeletons])
    add_skeleton_trio!(data_frame, "#", lateral_count_per_skeleton, pertinent_count_per_skeleton)

    # Each statistic is a `lat_<name>`, `lat_<name>_f`, `pert_<name>` trio (the lateral value, its fraction of the total,
    # and the pertinent value). Over the grouped markers it aggregates the gene report's per-marker values: the mean
    # (`μ`) of the improved (`++`) and degraded (`--`) base-block fractions and of the delta correlation (`μΔ`), and the
    # sum (`Σ`) of the delta correlations (`ΣΔ`).
    if change_per_marker_per_correlation !== nothing
        for (correlation_name, delta_per_marker, imp_per_marker, deg_per_marker) in change_per_marker_per_correlation
            lateral_sum_imp =
                sum_per_skeleton(top_skeleton_per_marker, is_lateral_per_marker, imp_per_marker, n_skeletons)
            pertinent_sum_imp =
                sum_per_skeleton(top_skeleton_per_marker, .!is_lateral_per_marker, imp_per_marker, n_skeletons)
            lateral_sum_deg =
                sum_per_skeleton(top_skeleton_per_marker, is_lateral_per_marker, deg_per_marker, n_skeletons)
            pertinent_sum_deg =
                sum_per_skeleton(top_skeleton_per_marker, .!is_lateral_per_marker, deg_per_marker, n_skeletons)
            lateral_sum_delta =
                sum_per_skeleton(top_skeleton_per_marker, is_lateral_per_marker, delta_per_marker, n_skeletons)
            pertinent_sum_delta =
                sum_per_skeleton(top_skeleton_per_marker, .!is_lateral_per_marker, delta_per_marker, n_skeletons)

            add_skeleton_trio!(
                data_frame,
                "$(correlation_name)_++",
                mean_per_skeleton(lateral_sum_imp, lateral_count_per_skeleton),
                mean_per_skeleton(pertinent_sum_imp, pertinent_count_per_skeleton),
            )
            add_skeleton_trio!(
                data_frame,
                "$(correlation_name)_--",
                mean_per_skeleton(lateral_sum_deg, lateral_count_per_skeleton),
                mean_per_skeleton(pertinent_sum_deg, pertinent_count_per_skeleton),
            )
            add_skeleton_trio!(
                data_frame,
                "$(correlation_name)_μΔ",
                mean_per_skeleton(lateral_sum_delta, lateral_count_per_skeleton),
                mean_per_skeleton(pertinent_sum_delta, pertinent_count_per_skeleton),
            )
            add_skeleton_trio!(data_frame, "$(correlation_name)_ΣΔ", lateral_sum_delta, pertinent_sum_delta)
        end
    end

    return data_frame[sortperm(total_count_per_skeleton; rev = true), :]
end

# The per-skeleton sum of a per-marker value over the markers of a category (grouped by their most-correlated skeleton).
function sum_per_skeleton(
    top_skeleton_per_marker::AbstractVector{<:Integer},
    is_in_category_per_marker::AbstractVector{Bool},
    value_per_marker::AbstractVector{<:Real},
    n_skeletons::Integer,
)::Vector{Float64}
    n_markers = length(top_skeleton_per_marker)
    sums_per_skeleton = zeros(Float64, n_skeletons)
    for marker_position in 1:n_markers
        if is_in_category_per_marker[marker_position]
            sums_per_skeleton[top_skeleton_per_marker[marker_position]] += value_per_marker[marker_position]
        end
    end
    return sums_per_skeleton
end

# The mean per skeleton (the sum divided by the count, or zero where the count is zero).
function mean_per_skeleton(
    sum_per_skeleton::AbstractVector{<:Real},
    count_per_skeleton::AbstractVector{<:Integer},
)::Vector{Float64}
    n_skeletons = length(count_per_skeleton)
    return [
        count_per_skeleton[skeleton] > 0 ? sum_per_skeleton[skeleton] / count_per_skeleton[skeleton] : 0.0 for
        skeleton in 1:n_skeletons
    ]
end

# Add a `lat_<name>`, `lat_<name>_f` (the lateral value's fraction of the lateral-plus-pertinent total), `pert_<name>`
# trio of columns to the report.
function add_skeleton_trio!(
    data_frame::DataFrame,
    name::AbstractString,
    lateral_per_skeleton::AbstractVector{<:Real},
    pertinent_per_skeleton::AbstractVector{<:Real},
)::Nothing
    n_skeletons = length(lateral_per_skeleton)
    data_frame[!, "lat_$(name)"] = lateral_per_skeleton
    data_frame[!, "lat_$(name)_f"] = [
        if lateral_per_skeleton[skeleton] + pertinent_per_skeleton[skeleton] != 0
            lateral_per_skeleton[skeleton] / (lateral_per_skeleton[skeleton] + pertinent_per_skeleton[skeleton])
        else
            0.0
        end for skeleton in 1:n_skeletons
    ]
    data_frame[!, "pert_$(name)"] = pertinent_per_skeleton
    return nothing
end

end  # module
