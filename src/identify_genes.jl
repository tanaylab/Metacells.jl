"""
Identify special genes.
"""
module IdentifyGenes

export compute_genes_divergence!
export identify_marker_genes!

using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using LinearAlgebra
using Random
using Statistics

using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.gene_axis
import Metacells.Contracts.gene_divergence_vector
import Metacells.Contracts.gene_is_marker_vector
import Metacells.Contracts.gene_metacell_fraction_matrix
import Metacells.Contracts.metacell_axis

"""
    function identify_marker_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_marker_gene_range_fold::AbstractFloat = $(DEFAULT.min_marker_gene_range_fold),
        min_marker_gene_max_fraction::AbstractFloat = $(DEFAULT.min_marker_gene_max_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

 1. Compute the minimal and maximal expression level of each gene.
 2. Compute the fold factor (log2 of maximal over minimal value, using the `gene_fraction_regularization`.
 3. Reduce this fold factor using the `divergence` factor of the genes.
 4. Identify as markers genes whose adjusted fold factor is at least `min_marker_gene_range_fold`, and whose maximal
    expression is at least `min_marker_gene_max_fraction`.

$(CONTRACT)
"""
@logged @computation Contract(
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        gene_metacell_fraction_matrix(RequiredInput),
        gene_divergence_vector(RequiredInput),
        gene_is_marker_vector(GuaranteedOutput),
    ],
) function identify_marker_genes!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    min_marker_gene_range_fold::AbstractFloat = 2.0,
    min_marker_gene_max_fraction::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_marker_gene_range_fold >= 0
    @assert min_marker_gene_max_fraction >= 0

    log_minimal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"]
    log_maximal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"]
    divergence_of_genes = daf["/ gene : divergence"]

    range_fold_of_genes = log_maximal_fraction_of_genes .- log_minimal_fraction_of_genes
    range_fold_of_genes .*= (1 .- divergence_of_genes)

    is_marker_of_genes =
        (range_fold_of_genes .>= min_marker_gene_range_fold) .&
        (log_maximal_fraction_of_genes .>= log2.(min_marker_gene_max_fraction + gene_fraction_regularization))

    set_vector!(daf, "gene", "is_marker", is_marker_of_genes; overwrite = overwrite)
    return nothing
end

"""
    function compute_genes_divergence!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_divergent_gene_range_fold::AbstractFloat = $(DEFAULT.min_divergent_gene_range_fold),
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
    data = [gene_metacell_fraction_matrix(RequiredInput), gene_divergence_vector(GuaranteedOutput)],
) function compute_genes_divergence!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION,
    min_divergent_gene_range_fold::AbstractFloat = 6.0,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_divergent_gene_range_fold >= 0

    log_minimal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"].array
    log_maximal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"].array

    range_fold_of_genes = log_maximal_fraction_of_genes .- log_minimal_fraction_of_genes
    range_fold_of_genes[range_fold_of_genes .< min_divergent_gene_range_fold] .= min_divergent_gene_range_fold
    divergence_of_genes = Vector{Float32}(1.0 .- min.(1.0, min_divergent_gene_range_fold ./ range_fold_of_genes))  # NOJET
    @assert all(divergence_of_genes .>= 0)
    @assert all(divergence_of_genes .< 1.0)

    set_vector!(daf, "gene", "divergence", divergence_of_genes; overwrite = overwrite)
    return nothing
end

end  # module
