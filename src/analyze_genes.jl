"""
Do simple per-gene analysis.
"""
module AnalyzeGenes

export fetch_transcription_factors!
export fetch_regulators!
export identify_genes_correlated_with_skeletons!
export identify_marker_genes!
export identify_skeleton_genes!
export rank_marker_genes!
export rank_markers

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
import Metacells.Contracts.gene_is_excluded_vector
import Metacells.Contracts.gene_is_forbidden_vector
import Metacells.Contracts.gene_is_lateral_vector
import Metacells.Contracts.gene_is_marker_vector
import Metacells.Contracts.gene_is_regulator_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.gene_is_transcription_factor_vector
import Metacells.Contracts.gene_marker_rank_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix

"""
    function identify_marker_genes!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_quantile::Real = $(DEFAULT.min_marker_quantile),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

Marker genes are genes which:

  - Have a maximal gene expression level of at least `min_marker_gene_max_fraction`, and
  - Have a range of expression of at least `min_marker_gene_range_fold` between the maximal and minimal
    (log base 2) of the expression level. For the minimal expression, we take the `min_marker_quantile`
    of the expression levels.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        gene_is_marker_vector(GuaranteedOutput),
    ],
) function identify_marker_genes!(  # UNTESTED
    daf::DafWriter;
    min_marker_gene_max_fraction::AbstractFloat = 1e-4,
    min_marker_gene_range_fold::Real = 2,
    min_marker_quantile::Real = 0,
    overwrite::Bool = false,
)::Nothing
    @assert min_marker_gene_range_fold >= 0
    @assert 0 <= min_marker_gene_max_fraction <= 1
    @assert 0 <= min_marker_quantile < 1

    maximal_fraction_per_gene = daf["/ metacell / gene : linear_fraction %> Max"].array
    if min_marker_quantile == 0
        log_minimal_fraction_per_gene = daf["/ metacell / gene : log_linear_fraction %> Min"].array
    else
        log_minimal_fraction_per_gene =
            daf["/ metacell / gene : log_linear_fraction %> Quantile p $(min_marker_quantile)"].array
    end
    log_maximal_fraction_per_gene = daf["/ metacell / gene : log_linear_fraction %> Max"].array
    range_fold_per_gene = log_maximal_fraction_per_gene .- log_minimal_fraction_per_gene

    is_marker_per_gene =
        (maximal_fraction_per_gene .>= min_marker_gene_max_fraction) .&
        (range_fold_per_gene .>= min_marker_gene_range_fold)

    @debug "Marker genes: $(sum(is_marker_per_gene))"

    set_vector!(daf, "gene", "is_marker", is_marker_per_gene; overwrite)
    return nothing
end

"""
    function identify_genes_correlated_with_skeletons!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_correlation_quantile::AbstractFloat = $(DEFAULT.min_marker_correlation_quantile),
        genes_correlation_window::Integer = $(DEFAULT.genes_correlation_window),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

TODOX

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        gene_is_skeleton_vector(RequiredInput),
        gene_is_correlated_with_skeleton_vector(GuaranteedOutput),
    ],
) function identify_genes_correlated_with_skeletons!(  # UNTESTED
    daf::DafWriter;
    min_gene_max_fraction::AbstractFloat = 1e-4,
    min_gene_range_fold::Real = 0.5,
    min_gene_correlation::AbstractFloat = 0.60,
    min_gene_correlation_quantile::AbstractFloat = 0.90,
    genes_correlation_window::Integer = 100,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= min_gene_max_fraction <= 1
    @assert 0 <= min_gene_range_fold <= 1
    @assert 0 <= min_gene_correlation_quantile < 1
    @assert genes_correlation_window > 0

    n_genes = axis_length(daf, "gene")

    maximal_fraction_per_gene = daf["/ metacell / gene : linear_fraction %> Max"].array
    is_strong_per_gene = maximal_fraction_per_gene .>= min_gene_max_fraction
    index_per_strong_gene = findall(is_strong_per_gene)
    n_strong_genes = length(index_per_strong_gene)

    average_fraction_per_gene = daf["/ metacell / gene : linear_fraction %> Mean"].array
    average_fraction_per_strong_gene = average_fraction_per_gene[index_per_strong_gene]

    position_per_ordered_strong_gene = sortperm(average_fraction_per_strong_gene)
    index_per_ordered_strong_gene = index_per_strong_gene[position_per_ordered_strong_gene]
    average_fraction_per_ordered_strong_gene = average_fraction_per_strong_gene[position_per_ordered_strong_gene]

    log_fraction_per_metacell_per_gene = daf["/ metacell / gene : log_linear_fraction"].array
    @views log_fraction_per_metacell_per_ordered_strong_gene =
        log_fraction_per_metacell_per_gene[:, index_per_ordered_strong_gene]

    log_fraction_per_metacell_per_skeleton = daf["/ metacell / gene & is_skeleton : log_linear_fraction"].array
    n_skeletons = size(log_fraction_per_metacell_per_skeleton, 2)

    correlation_per_skeleton_per_ordered_strong_gene =
        cor(log_fraction_per_metacell_per_skeleton, log_fraction_per_metacell_per_ordered_strong_gene)
    correlation_per_skeleton_per_ordered_strong_gene[isnan.(correlation_per_skeleton_per_ordered_strong_gene)] .= 0.0
    @assert_matrix(correlation_per_skeleton_per_ordered_strong_gene, n_skeletons, n_strong_genes)

    max_correlation_per_ordered_strong_gene = vec(maximum(correlation_per_skeleton_per_ordered_strong_gene; dims = 1))
    @assert_vector(max_correlation_per_ordered_strong_gene, n_strong_genes)

    quantile_correlation_per_ordered_strong_gene = rolling_quantile(
        max_correlation_per_ordered_strong_gene,
        genes_correlation_window,
        min_gene_correlation_quantile,
    )
    is_quantile_correlation_per_ordered_strong_gene =
        max_correlation_per_ordered_strong_gene .>= quantile_correlation_per_ordered_strong_gene
    is_strong_correlation_per_ordered_strong_gene = max_correlation_per_ordered_strong_gene .>= min_gene_correlation

    is_correlated_with_skeleton_per_gene = zeros(Bool, n_genes)
    is_correlated_with_skeleton_per_gene[index_per_ordered_strong_gene[is_quantile_correlation_per_ordered_strong_gene]] .=
        true
    is_correlated_with_skeleton_per_gene[index_per_ordered_strong_gene[is_strong_correlation_per_ordered_strong_gene]] .=
        true

    min_log_fraction_per_metacell_per_gene = daf["/ metacell / gene : log_linear_fraction %> Min"].array
    max_log_fraction_per_metacell_per_gene = daf["/ metacell / gene : log_linear_fraction %> Max"].array
    range_fold_log_fraction_per_metacell_per_gene =
        max_log_fraction_per_metacell_per_gene .- min_log_fraction_per_metacell_per_gene
    is_range_fold_per_gene = range_fold_log_fraction_per_metacell_per_gene .>= min_gene_range_fold
    is_correlated_with_skeleton_per_gene .&= is_range_fold_per_gene

    @debug "Correlated with skeleton genes: $(sum(is_correlated_with_skeleton_per_gene)) out of strong genes: $(n_strong_genes)"

    set_vector!(daf, "gene", "is_correlated_with_skeleton", is_correlated_with_skeleton_per_gene; overwrite)
    return nothing
end

function rolling_quantile(
    data::AbstractVector{T},
    window_size::Int,
    quantile_fraction::AbstractFloat,
)::AbstractVector{T} where {T}
    @assert window_size > 1

    n_values = length(data)
    if n_values <= window_size
        value = quantile(data, quantile_fraction)
        return fill(value, n_values)
    end

    results = Vector{T}(undef, n_values)

    window_left_of_center = div(window_size - 1, 2)
    window_right_of_center = window_size - window_left_of_center - 1

    for index in 1:(n_values - window_size + 1)
        @views window_data = data[index:(index + window_size - 1)]
        results[window_left_of_center + index] = quantile(window_data, quantile_fraction)
    end

    first_median = results[window_left_of_center + 1]
    results[1:window_left_of_center] .= first_median

    last_median = results[n_values - window_right_of_center]
    results[(end - window_right_of_center + 1):end] .= last_median

    @assert !any(isnan.(results)) # TODOX

    return results
end

"""
    function rank_marker_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the relative ranks of marker genes.

 1. Compute the median of the (log base 2) of the gene expression level for each marker gene across all the metacells.
 2. Compute the per-marker-per-metacell fold factor (absolute difference of the log expression from the median).
 3. Rank the markers for each metacell (1 having the largest fold factor relative to the median).
 4. For each marker, give it a priority which is a tuple of (1) the minimal rank it has in all metacells (2) the maximal
    fold it has in metacells where it has that rank (negated).
 5. Sort the markers according to this priority.

Non-marker genes are given a rank of `typemax(UInt32)` regardless of their expression level in the metacells.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        gene_is_marker_vector(RequiredInput),
        gene_marker_rank_vector(GuaranteedOutput),
    ],
) function rank_marker_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_metacells = axis_length(daf, "metacell")
    n_genes = axis_length(daf, "gene")

    is_marker_per_gene = get_vector(daf, "gene", "is_marker")
    n_markers = sum(is_marker_per_gene)

    log_fraction_per_metacell_per_marker = daf["/ metacell / gene & is_marker : log_linear_fraction"].array
    @assert_matrix(log_fraction_per_metacell_per_marker, n_metacells, n_markers, Columns)

    median_log_fraction_per_marker = daf["/ metacell / gene & is_marker : log_linear_fraction %> Median"].array
    @assert_vector(median_log_fraction_per_marker, n_markers)

    abs_fold_per_metacell_per_marker =  # NOJET
        abs.(log_fraction_per_metacell_per_marker .- transpose(median_log_fraction_per_marker))
    abs_fold_per_marker_per_metacell = flipped(abs_fold_per_metacell_per_marker)

    rank_per_marker = rank_markers(abs_fold_per_marker_per_metacell)

    marker_rank_per_gene = fill(typemax(UInt32), n_genes)
    marker_rank_per_gene[is_marker_per_gene] .= rank_per_marker

    set_vector!(daf, "gene", "marker_rank", marker_rank_per_gene; overwrite)
    return nothing
end

"""
TODOX
"""
function rank_markers(score_per_gene_per_something::AbstractMatrix{<:AbstractFloat})::AbstractVector{<:Integer}
    @assert_matrix(score_per_gene_per_something, Columns)
    n_genes, n_somethings = size(score_per_gene_per_something)

    @views rank_per_gene_per_something = hcat(
        [
            invperm(sortperm(vec(score_per_gene_per_something[:, something_index]); rev = true)) for
            something_index in 1:n_somethings
        ]...,
    )
    @assert_matrix(rank_per_gene_per_something, n_genes, n_somethings, Columns)

    min_rank_per_gene = vec(minimum(rank_per_gene_per_something; dims = 2))
    @assert_vector(min_rank_per_gene, n_genes)

    maximal_score_per_gene = Vector{Float32}(undef, n_genes)
    for gene_index in 1:n_genes
        @views maximal_score_per_something = score_per_gene_per_something[
            gene_index,
            rank_per_gene_per_something[gene_index, :] .== min_rank_per_gene[gene_index],
        ]
        maximal_score_per_gene[gene_index] = mean(maximal_score_per_something)
    end

    priority_per_gene = collect(zip(min_rank_per_gene, -maximal_score_per_gene))
    rank_per_gene = invperm(sortperm(priority_per_gene))
    @assert_vector(rank_per_gene, n_genes)
    return rank_per_gene
end

"""
    function identify_skeleton_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the skeleton genes that will be used to predict the rest of the genes. We just pick the genes that are
regulators, and are not excluded, lateral or forbidden.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        gene_is_regulator_vector(RequiredInput),
        gene_is_excluded_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        gene_is_forbidden_vector(RequiredInput),
        gene_is_marker_vector(RequiredInput),
        gene_is_skeleton_vector(GuaranteedOutput),
    ],
) function identify_skeleton_genes!(  # UNTESTED
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

    name_per_gene = axis_vector(daf, "gene")
    @info "Skeleton genes: [ $(join(sort(name_per_gene[is_skeleton_per_gene]), ", ")) ]"

    set_vector!(daf, "gene", "is_skeleton", is_skeleton_per_gene; overwrite)
    @assert any(is_skeleton_per_gene)
    return nothing
end

"""
    fetch_regulators!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString = $(DEFAULT.namespace),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = $(DEFAULT.cache_dir),
        timeout::Real = $(DEFAULT.timeout),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Fetch the names of transcription factor genes from `Gmara` and create a mask based on them. These genes will be assumed
to be the core genes driving gene programs that describe cell behavior.
"""
@logged @computation Contract(axes = [gene_axis(RequiredInput)], data = [gene_is_regulator_vector(GuaranteedOutput)]) function fetch_regulators!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    version::AbstractString = function_default(set_gmara_genes_mask!, :version),
    cache_dir = function_default(set_gmara_genes_mask!, :cache_dir),
    timeout::Real = function_default(set_gmara_genes_mask!, :timeout),
    overwrite::Bool = false,
)::Nothing
    set_gmara_genes_mask!(daf; species, namespace, list = "regulator", version, cache_dir, timeout, overwrite)
    return nothing
end

"""
    fetch_transcription_factors!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString = $(DEFAULT.namespace),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = $(DEFAULT.cache_dir),
        timeout::Real = $(DEFAULT.timeout),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Fetch the names of transcription factor genes from `Gmara` and create a mask based on them. These genes bind with the
DNA, but most aren't meaningful regulator of cell behavior.
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [gene_is_transcription_factor_vector(GuaranteedOutput)],
) function fetch_transcription_factors!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    version::AbstractString = function_default(set_gmara_genes_mask!, :version),
    cache_dir = function_default(set_gmara_genes_mask!, :cache_dir),
    timeout::Real = function_default(set_gmara_genes_mask!, :timeout),
    overwrite::Bool = false,
)::Nothing
    set_gmara_genes_mask!(
        daf;
        species,
        namespace,
        list = "transcription_factor",
        version,
        cache_dir,
        timeout,
        overwrite,
    )
    return nothing
end

end  # module
