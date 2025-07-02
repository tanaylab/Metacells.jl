"""
Do simple per-gene analysis.
"""
module AnalyzeGenes

export fetch_regulators!
export identify_covered_genes!
export identify_marker_genes!
export identify_skeleton_genes!
export rank_marker_genes!

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
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_forbidden_vector
import Metacells.Contracts.gene_is_lateral_vector
import Metacells.Contracts.gene_is_marker_vector
import Metacells.Contracts.gene_is_regulator_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.gene_marker_rank_vector
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.metacell_gene_linear_fraction_matrix
import Metacells.Contracts.metacell_gene_log_linear_fraction_matrix

"""
    function identify_marker_genes!(
        daf::DafWriter;
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

Marker genes are genes which:

  - Have a maximal gene expression level of at least `min_marker_gene_max_fraction`, and
  - Have a range of expression of at least `min_marker_gene_range_fold` between the maximal and minimal
    (log base 2) of the expression level.

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
    min_marker_gene_range_fold::Real = 2,
    min_marker_gene_max_fraction::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert min_marker_gene_range_fold >= 0
    @assert 0 <= min_marker_gene_max_fraction <= 1

    maximal_fraction_per_gene = daf["/ metacell / gene : linear_fraction %> Max"]
    log_minimal_fraction_per_gene = daf["/ metacell / gene : log_linear_fraction %> Min"]
    log_maximal_fraction_per_gene = daf["/ metacell / gene : log_linear_fraction %> Max"]
    range_fold_per_gene = log_maximal_fraction_per_gene .- log_minimal_fraction_per_gene

    is_marker_per_gene =
        (maximal_fraction_per_gene .>= min_marker_gene_max_fraction) .&
        (range_fold_per_gene .>= min_marker_gene_range_fold)

    @debug "Marker genes: $(sum(is_marker_per_gene))"

    set_vector!(daf, "gene", "is_marker", is_marker_per_gene; overwrite)
    return nothing
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

    negative_fold_per_metacell_per_marker =  # NOJET
        .-abs.(log_fraction_per_metacell_per_marker .- transpose(median_log_fraction_per_marker))
    negative_fold_per_marker_per_metacell = transposer(negative_fold_per_metacell_per_marker)

    @views rank_per_marker_per_metacell = hcat(
        [
            invperm(sortperm(vec(negative_fold_per_marker_per_metacell[:, metacell_index]))) for
            metacell_index in 1:n_metacells
        ]...,
    )
    @assert_matrix(rank_per_marker_per_metacell, n_markers, n_metacells, Columns)

    min_rank_per_marker = vec(minimum(rank_per_marker_per_metacell; dims = 2))
    @assert_vector(min_rank_per_marker, n_markers)

    @views most_negative_fold_per_marker = [
        mean(  # NOLINT
            negative_fold_per_metacell_per_marker[
                rank_per_marker_per_metacell[marker_index, :] .== min_rank_per_marker[marker_index],
                marker_index,
            ],
        ) for marker_index in 1:n_markers
    ]

    priority_per_marker = collect(zip(min_rank_per_marker, most_negative_fold_per_marker))
    rank_per_marker = invperm(sortperm(priority_per_marker))

    marker_rank_per_gene = fill(typemax(UInt32), n_genes)
    marker_rank_per_gene[is_marker_per_gene] .= rank_per_marker

    set_vector!(daf, "gene", "marker_rank", marker_rank_per_gene; overwrite)
    return nothing
end

"""
    function identify_covered_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that will be approximated by the local linear programs. Picking them is simple: we cover all marker
genes, that are not lateral.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        gene_is_covered_vector(GuaranteedOutput),
    ],
) function identify_covered_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_covered_per_gene = is_marker_per_gene .& .!is_lateral_per_gene
    set_vector!(daf, "gene", "is_covered", is_covered_per_gene; overwrite)
    return nothing
end

"""
    function identify_skeleton_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the skeleton genes that will be used to predict the rest of the (covered) genes. We just pick the covered genes
that are also regulators, and are not forbidden.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        gene_is_regulator_vector(RequiredInput),
        gene_is_forbidden_vector(RequiredInput),
        gene_is_skeleton_vector(GuaranteedOutput),
    ],
) function identify_skeleton_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_covered_per_gene = get_vector(daf, "gene", "is_covered").array
    is_forbidden_per_gene = get_vector(daf, "gene", "is_forbidden").array
    is_regulator_per_gene = get_vector(daf, "gene", "is_regulator").array

    is_skeleton_per_gene = is_covered_per_gene .& .!is_forbidden_per_gene .& is_regulator_per_gene

    name_per_gene = axis_vector(daf, "gene")
    @info "Skeleton genes: [ $(join(sort(name_per_gene[is_skeleton_per_gene]), ", ")) ]"

    set_vector!(daf, "gene", "is_skeleton", is_skeleton_per_gene; overwrite)
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

end  # module
