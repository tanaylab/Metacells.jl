"""
Identify special genes.
"""
module IdentifyGenes

export compute_genes_divergence!
export rank_marker_genes!
export identify_covered_genes!
export identify_marker_genes!
export identify_skeleton_genes!
export identify_uncorrelated_genes!

using DataAxesFormats
using LinearAlgebra
using Random
using Statistics
using TanayLabUtilities

using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_covered_vector
import Metacells.Contracts.gene_is_lateral_vector
import Metacells.Contracts.gene_is_marker_vector
import Metacells.Contracts.gene_is_skeleton_vector
import Metacells.Contracts.gene_is_uncorrelated_vector
import Metacells.Contracts.gene_marker_rank_vector
import Metacells.Contracts.gene_marker_rank_vector
import Metacells.Contracts.metacell_gene_fraction_matrix
import Metacells.Contracts.metacell_axis

"""
    function identify_marker_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_marker_gene_range_fold::Real = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

 1. Compute the minimal and maximal expression level of each gene.
 2. Compute the fold factor (log2 of maximal over minimal value, using the `gene_fraction_regularization`.
 4. Identify as markers genes whose adjusted fold factor is at least `min_marker_gene_range_fold`, and whose maximal
    expression is at least `min_marker_gene_max_fraction`.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [metacell_gene_fraction_matrix(RequiredInput), gene_is_marker_vector(GuaranteedOutput)],
) function identify_marker_genes!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    min_marker_gene_range_fold::Real = 2,
    min_marker_gene_max_fraction::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_marker_gene_range_fold >= 0
    @assert min_marker_gene_max_fraction >= 0

    log_minimal_fraction_per_gene =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"]
    log_maximal_fraction_per_gene =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"]

    range_fold_per_gene = log_maximal_fraction_per_gene .- log_minimal_fraction_per_gene

    is_marker_per_gene =
        (range_fold_per_gene .>= min_marker_gene_range_fold) .&
        (log_maximal_fraction_per_gene .>= log2.(min_marker_gene_max_fraction + gene_fraction_regularization))

    set_vector!(daf, "gene", "is_marker", is_marker_per_gene; overwrite)
    return nothing
end

"""
    function rank_marker_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute the relative ranks of marker genes.

1. Compute the log base 2 of the expression of the marker genes in the metacells (using the `gene_fraction_regularization`).
2. Compute the median of this for each gene across all the metacells.
3. Compute the per-marker-per-metacell fold factor (absolute difference of the log expression from the median).
4. Rank the markers for each metacell (1 having the largest fold factor relative to the median).
5. For each markers, give it a priority which is a tuple of (1) the minimal rank it has in all metacells (2) the maximal
   fold it has in metacells where it has that rank (negated).
6. Sort the markers according to this priority.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        metacell_gene_fraction_matrix(RequiredInput),
        gene_is_marker_vector(RequiredInput),
        gene_marker_rank_vector(GuaranteedOutput),
    ],
) function rank_marker_genes!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    overwrite::Bool = false,
)::Nothing
    fraction_per_metacell_per_marker = daf["/ metacell / gene & is_marker : fraction"].array
    n_metacells, n_markers = size(fraction_per_metacell_per_marker)
    @assert n_metacells == axis_length(daf, "metacell")

    log_fraction_per_metacell_per_marker = log2.(fraction_per_metacell_per_marker .+ gene_fraction_regularization)

    median_log_fraction_per_marker = vec(median(log_fraction_per_metacell_per_marker; dims = 1))
    @assert_vector(median_log_fraction_per_marker, n_markers)

    negative_fold_per_metacell_per_marker =
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
        mean(
            negative_fold_per_metacell_per_marker[
                rank_per_marker_per_metacell[marker_index, :] .== min_rank_per_marker[marker_index],
                marker_index,
            ],
        ) for marker_index in 1:n_markers
    ]

    priority_per_marker = collect(zip(min_rank_per_marker, most_negative_fold_per_marker))
    rank_per_marker = invperm(sortperm(priority_per_marker))

    n_genes = axis_length(daf, "gene")
    is_marker_per_gene = get_vector(daf, "gene", "is_marker")
    marker_rank_per_gene = fill(typemax(UInt32), n_genes)
    marker_rank_per_gene[is_marker_per_gene] .= rank_per_marker
    set_vector!(daf, "gene", "marker_rank", marker_rank_per_gene; overwrite)
    return nothing
end

"""
    function compute_genes_divergence!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_divergent_gene_range_fold::Real = $(DEFAULT.min_divergent_gene_range_fold),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute a divergence factor for all genes. This factor is typically 0.0 for most genes (and is always below 1.0). Genes
that have a wide range of expression are given a higher divergence factor. When looking for significant differences in
gene expressions, the divergence factor is used to scale down the difference, to compensate for the gene's inherent
wider expression range. Specifically, we take the raw fold factor and multiply it by (1.0 - divergence).

This works as follows:

 1. Compute the minimal and maximal expression level of each gene.
 2. Compute the fold factor (log2 of maximal over minimal value, using the `gene_fraction_regularization`.
 3. The raw factor of the gene is the ratio between the fold factor and the `min_divergent_gene_range_fold`.
    If this is higher than 1, we set the divergence factor accordingly.

!!! note

    Ideally, all code that uses the manual `is_noisy` mask should be modified to use the divergence factor instead.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [metacell_gene_fraction_matrix(RequiredInput), gene_divergence_vector(GuaranteedOutput)],
) function compute_genes_divergence!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    min_divergent_gene_range_fold::Real = 6,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_divergent_gene_range_fold >= 0

    log_minimal_fraction_per_gene =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"].array
    log_maximal_fraction_per_gene =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"].array

    range_fold_per_gene = log_maximal_fraction_per_gene .- log_minimal_fraction_per_gene
    range_fold_per_gene[range_fold_per_gene .< min_divergent_gene_range_fold] .= min_divergent_gene_range_fold
    divergence_per_gene = Vector{Float32}(1.0 .- min.(1.0, min_divergent_gene_range_fold ./ range_fold_per_gene))  # NOJET
    @assert all(divergence_per_gene .>= 0)
    @assert all(divergence_per_gene .< 1.0)

    set_vector!(daf, "gene", "divergence", divergence_per_gene; overwrite)
    return nothing
end

"""
    function identify_uncorrelated_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        correlation_confidence::AbstractFloat = $(DEFAULT.correlation_confidence),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify genes that are correlated with other gene(s). Such genes are good candidates for looking for groups of genes
that act together. If `overwrite`, will overwrite an existing `is_correlated` mask.

 1. Compute the log base 2 of the genes expression in each metacell (using the `gene_fraction_regularization`).
 2. Correlate this between all the pairs of genes.
 3. For each gene, shuffle its values along all metacells, and again correlate this between all the pairs of genes.
 4. Find the maximal absolute correlation for each gene in both cases (that is, strong anti-correlation also counts).
 5. Find the `correlation_confidence` quantile correlation of the shuffled data.
 6. Identify the genes that have at least that level of correlations in the unshuffled data.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [metacell_gene_fraction_matrix(RequiredInput), gene_is_uncorrelated_vector(GuaranteedOutput)],
) function identify_uncorrelated_genes!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    correlation_confidence::AbstractFloat = 0.99,
    overwrite::Bool = false,
    rng::AbstractRNG = default_rng(),
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert 0 <= correlation_confidence <= 1

    n_genes = axis_length(daf, "gene")
    log_fractions_per_gene_in_metacells =
        daf["/ metacell / gene : fraction % Log base 2 eps $(gene_fraction_regularization)"].array

    correlations_between_genes = cor(log_fractions_per_gene_in_metacells)
    correlations_between_genes .= abs.(correlations_between_genes)
    correlations_between_genes[diagind(correlations_between_genes)] .= 0  # NOJET
    correlations_between_genes[isnan.(correlations_between_genes)] .= 0
    max_correlations_per_gene = vec(maximum(correlations_between_genes; dims = 1))
    @assert length(max_correlations_per_gene) == n_genes

    shuffled_log_fraction_per_gene_per_metacell = copy_array(log_fractions_per_gene_in_metacells)
    for gene_index in 1:n_genes
        @views shuggled_log_fractions_of_gene_per_metacell = shuffled_log_fraction_per_gene_per_metacell[:, gene_index]
        shuffle!(rng, shuggled_log_fractions_of_gene_per_metacell)
    end

    shuffled_correlations_between_genes = cor(shuffled_log_fraction_per_gene_per_metacell)
    shuffled_correlations_between_genes .= abs.(shuffled_correlations_between_genes)
    shuffled_correlations_between_genes[diagind(shuffled_correlations_between_genes)] .= 0  # NOJET
    shuffled_correlations_between_genes[isnan.(shuffled_correlations_between_genes)] .= 0
    max_shuffled_correlations_per_gene = vec(maximum(shuffled_correlations_between_genes; dims = 1))

    @debug "mean: $(mean(max_shuffled_correlations_per_gene))"
    @debug "stdev: $(std(max_shuffled_correlations_per_gene))"
    threshold = quantile(max_shuffled_correlations_per_gene, correlation_confidence)
    @debug "threshold: $(threshold)"

    is_uncorrelated_per_gene = max_correlations_per_gene .< threshold

    set_vector!(daf, "gene", "is_uncorrelated", is_uncorrelated_per_gene; overwrite)
    return nothing
end

"""
    function identify_covered_genes!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that will be approximated by the local linear programs. Picking them is simple: we cover all marker
genes that are not lateral and not uncorrelated.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        gene_is_marker_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        gene_is_uncorrelated_vector(RequiredInput),
        gene_is_covered_vector(GuaranteedOutput),
    ],
) function identify_covered_genes!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    is_marker_per_gene = get_vector(daf, "gene", "is_marker").array
    is_lateral_per_gene = get_vector(daf, "gene", "is_lateral").array
    is_uncorrelated_per_gene = get_vector(daf, "gene", "is_uncorrelated").array
    is_covered_per_gene = is_marker_per_gene .& .!is_lateral_per_gene .& .!is_uncorrelated_per_gene
    set_vector!(daf, "gene", "is_covered", is_covered_per_gene; overwrite)
    return nothing
end

"""
    function identify_skeleton_genes!(
        daf::DafWriter;
        max_skeleton_genes::Integer = $(DEFAULT.max_skeleton_genes),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the skeleton genes that will be used to predict the rest of the (covered) genes. We just pick the
`max_skeleton_genes` that have the lowest `marker_rank` out of the covered genes.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput)],
    data = [
        gene_is_covered_vector(RequiredInput),
        gene_marker_rank_vector(RequiredInput),
        gene_is_skeleton_vector(GuaranteedOutput),
    ],
) function identify_skeleton_genes!(  # UNTESTED
    daf::DafWriter;
    max_skeleton_genes::Integer = 200,
    overwrite::Bool = false,
)::Nothing
    @assert max_skeleton_genes > 0

    is_covered_per_gene = get_vector(daf, "gene", "is_covered").array
    n_covered = sum(is_covered_per_gene)

    marker_rank_per_gene = get_vector(daf, "gene", "marker_rank").array

    priority_per_gene = copy_array(marker_rank_per_gene)
    priority_per_gene[.!is_covered_per_gene] .= typemax(eltype(priority_per_gene))

    n_skeletons = min(max_skeleton_genes, n_covered)
    index_per_skeleton_gene = partialsortperm(priority_per_gene, 1:n_skeletons)  # NOJET

    n_genes = axis_length(daf, "gene")
    is_skeleton_of_gene = zeros(Bool, n_genes)
    is_skeleton_of_gene[index_per_skeleton_gene] .= true

    name_per_gene = axis_vector(daf, "gene")
    @info "Skeleton genes: [ $(join(name_per_gene[index_per_skeleton_gene], ", ")) ]"

    set_vector!(daf, "gene", "is_skeleton", is_skeleton_of_gene; overwrite)
    return nothing
end

end  # module
