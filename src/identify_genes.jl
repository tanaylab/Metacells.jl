"""
Identify special genes.
"""
module IdentifyGenes

export identify_correlated_genes!
export identify_marker_genes!

using Daf
using Daf.GenericLogging
using LinearAlgebra
using Statistics

"""
    function identify_marker_genes!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = 1e-5,
        min_gene_range_fold::AbstractFloat = 2.0,
        noisy_gene_fold::AbstractFloat = 1.0,
        min_max_marker_gene_fraction::AbstractFloat = 1e-4,
        overwrite::Bool = false,
    )::Nothing

Identify the genes that distinguish at least one metacell from the rest. Such genes are called "marker" genes as they
(potentially) mark specific cell states. If `overwrite`, will overwrite an existing `is_marker` mask.

 1. Compute the minimal and maximal expression level of each gene.

 2. Select the genes whose fold factor (log2 of maximal over minimal value, using the `gene_fraction_regularization`
    is at least `min_marker_gene_range_fold`. For `is_noisy` genes, we require an additional `noisy_gene_fold`.
 3. Identify the genes whose maximal expression is at least `min_max_marker_gene_fraction`.

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
        ("gene", "is_noisy") =>
            (OptionalInput, Bool, "A mask of noisy genes to be given additional diameter when grouped."),
        ("gene", "is_marker") => (GuaranteedOutput, Bool, "A mask of genes that distinguish between metacells."),
    ],
) function identify_marker_genes!(  # untested
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = 1e-5,
    min_marker_gene_range_fold::AbstractFloat = 2.0,
    noisy_gene_fold::AbstractFloat = 1.0,
    min_max_marker_gene_fraction::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization > 0
    @assert min_marker_gene_range_fold > 0
    @assert noisy_gene_fold >= 0
    @assert min_max_marker_gene_fraction >= 0

    log_minimal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Min % Log base 2 eps $(gene_fraction_regularization)"]
    log_maximal_fraction_of_genes =
        daf["/ metacell / gene : fraction %> Max % Log base 2 eps $(gene_fraction_regularization)"]

    range_fold_of_genes = log_maximal_fraction_of_genes .- log_minimal_fraction_of_genes
    is_noisy_of_genes = daf[q"/ gene : is_noisy || false"].array
    range_fold_of_genes[is_noisy_of_genes] .-= noisy_gene_fold

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
