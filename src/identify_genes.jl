"""
Identify special genes.
"""
module IdentifyGenes

export compute_genes_divergence!
export identify_correlated_genes!
export identify_marker_genes!

using Daf
using Daf.GenericLogging
using LinearAlgebra
using Statistics

"""
    function compute_genes_divergence!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = 1e-5,
        min_divergent_gene_range_fold::AbstractFloat = 6.0,
        overwrite::Bool = false,
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

CONTRACT
"""
@logged @computation Contract(
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (RequiredInput, "The genes to consider."),
    ],
    data = [
        ("metacell", "gene", "fraction") =>
            (RequiredInput, AbstractFloat, "The fraction of the UMIs of each gene in each metacell."),
        ("gene", "divergence") => (GuaranteedOutput, AbstractFloat, "How to scale fold factors for this gene."),
    ],
) function compute_genes_divergence!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 1e-5,
    min_divergent_gene_range_fold::AbstractFloat = 6.0,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    @assert min_divergent_gene_range_fold > 0

    log_minimal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"].array
    log_maximal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"].array

    range_fold_of_genes = log_maximal_fraction_of_genes .- log_minimal_fraction_of_genes
    range_fold_of_genes[range_fold_of_genes .< min_divergent_gene_range_fold] .= min_divergent_gene_range_fold
    divergence_of_genes = Vector{Float32}(1.0 .- min.(1.0, min_divergent_gene_range_fold ./ range_fold_of_genes))
    @assert all(divergence_of_genes .>= 0)
    @assert all(divergence_of_genes .< 1.0)

    set_vector!(daf, "gene", "divergence", divergence_of_genes; overwrite = overwrite)
    return nothing
end

"""
    function identify_marker_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = 1e-5,
        min_gene_range_fold::AbstractFloat = 2.0,
        min_max_marker_gene_fraction::AbstractFloat = 1e-4,
        overwrite::Bool = false,
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

 1. Compute the minimal and maximal expression level of each gene.
 2. Compute the fold factor (log2 of maximal over minimal value, using the `gene_fraction_regularization`.
 3. Reduce this fold factor using the `divergence` factor of the genes.
 4. Identify as markers genes whose adjusted fold factor is at least `min_marker_gene_range_fold`, and whose maximal
    expression is at least `min_max_marker_gene_fraction`.

CONTRACT
"""
@logged @computation Contract(
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (RequiredInput, "The genes to consider."),
    ],
    data = [
        ("metacell", "gene", "fraction") =>
            (RequiredInput, AbstractFloat, "The fraction of the UMIs of each gene in each metacell."),
        ("gene", "divergence") => (RequiredInput, AbstractFloat, "How to scale fold factors for this gene."),
        ("gene", "is_marker") => (GuaranteedOutput, Bool, "A mask of genes that distinguish between metacells."),
    ],
) function identify_marker_genes!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 1e-5,
    min_marker_gene_range_fold::AbstractFloat = 2.0,
    min_max_marker_gene_fraction::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    @assert min_marker_gene_range_fold > 0
    @assert min_max_marker_gene_fraction >= 0

    log_minimal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"]
    log_maximal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"]
    divergence_of_genes = daf["/ gene : divergence"]

    range_fold_of_genes = log_maximal_fraction_of_genes .- log_minimal_fraction_of_genes
    range_fold_of_genes .*= (1 .- divergence_of_genes)

    is_marker_of_genes =
        (range_fold_of_genes .>= min_marker_gene_range_fold) .&
        (log_maximal_fraction_of_genes .>= log2.(min_max_marker_gene_fraction + gene_fraction_regularization))

    set_vector!(daf, "gene", "is_marker", is_marker_of_genes; overwrite = overwrite)
    return nothing
end

"""
    function identify_correlated_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = 1e-5,
        min_gene_correlation::AbstractFloat = 0.5,
        overwrite::Bool = false,
    )::Nothing

Identify genes that are correlated with other gene(s). Such genes are good candidates for looking for groups of genes
that act together. If `overwrite`, will overwrite an existing `is_correlated` mask.

 1. Compute the log base 2 of the genes expression in each metacell (using the `gene_fraction_regularization`).
 2. Correlate this between all the pairs of genes.
 3. Find the maximal absolute correlation for each gene (that is, strong anti-correlation also counts).
 4. Identify the genes which have at least one gene with a correlation of at least `min_gene_correlation`.

CONTRACT
"""
@logged @computation Contract(
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (RequiredInput, "The genes to consider (typically, only marker genes)."),
    ],
    data = [
        ("metacell", "gene", "fraction") =>
            (RequiredInput, AbstractFloat, "The fraction of the UMIs of each gene in each metacell."),
        ("gene", "is_correlated") =>
            (GuaranteedOutput, Bool, "A mask of genes that are correlated with other gene(s)."),
    ],
) function identify_correlated_genes!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 1e-5,
    min_gene_correlation::AbstractFloat = 0.5,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    @assert 0 <= min_gene_correlation <= 1

    log_fractions = daf["/ metacell / gene : fraction % Log base 2 eps $(gene_fraction_regularization)"].array
    genes_correlations = cor(log_fractions)
    genes_correlations .= abs.(genes_correlations)

    genes_correlations[diagind(genes_correlations)] .= 0  # NOJET
    is_correlated_of_genes = vec(any(abs.(genes_correlations) .>= min_gene_correlation; dims = 1))

    set_vector!(daf, "gene", "is_correlated", is_correlated_of_genes; overwrite = overwrite)
    return nothing
end

end  # module
