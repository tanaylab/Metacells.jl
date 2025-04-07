"""
Functions for defining a `Contract` for a metacells `Daf` data set `@computation`. This also serves as a vocabulary
describing the data we keep when doing metacells based analysis, which is a solid basis for understanding what is going
on. As Fred Brooks said: "Show me your flowcharts and conceal your tables, and I shall continue to be mystified. Show me
your tables, and I won’t usually need your flowcharts; they’ll be obvious.".

In the descriptions below, "fold factor" refers to the log base 2 of the ratio between expression levels. For fold
factors between gene RNA expression levels (fractions), we typically use a regularization factor of 1e-5 to avoid
division by zero.
"""
module Contracts

export block_axis
export block_block_is_in_environment_matrix
export block_block_is_in_neighborhood_matrix
export block_gene_fraction_matrix
export block_gene_mean_scaled_log_covered_fraction_matrix
export block_gene_total_UMIs_matrix
export block_linear_RMSE_vector
export block_linear_XRMSE_vector
export block_n_principal_components_vector
export block_principal_component_gene_covered_coefficient_tensor
export block_principal_component_gene_skeleton_coefficient_tensor
export block_principal_component_is_used_matrix
export cell_axis
export cell_gene_UMIs_matrix
export cell_metacell_vector
export cell_new_metacell_vector
export cell_old_metacell_vector
export cell_total_UMIs_vector
export gene_axis
export gene_divergence_vector
export gene_is_covered_vector
export gene_is_forbidden_vector
export gene_is_lateral_vector
export gene_is_marker_vector
export gene_is_skeleton_vector
export gene_is_uncorrelated_vector
export gene_marker_rank_vector
export metacell_axis
export metacell_block_vector
export metacell_covered_UMIs_vector
export metacell_gene_covered_fraction_matrix
export metacell_gene_fraction_matrix
export metacell_gene_scaled_log_covered_fraction_matrix
export metacell_gene_total_UMIs_matrix
export metacell_total_UMIs_vector
export metacell_type_vector
export new_metacell_axis
export new_metacell_block_vector
export old_metacell_axis
export old_metacell_block_vector
export old_metacell_gene_scaled_log_covered_fraction_matrix
export principal_component_axis
export type_axis
export type_color_vector

using DataAxesFormats

"""
    gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced genes. By convention we use gene symbols as the namespace of the genes, but this may be different
depending on the data set.
"""
function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "gene" => (expectation, "Sequenced genes.")
end

"""
    cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced single cells. There's no convention for cell names, as long as they are unique. Typically some
sort of barcode is used, possibly combined with a batch and/or plate and/or experiment identification. In the latter
case it is recommended that `batch` and/or `plate` and/or `experiment` would also be created as explicit axes, to allow
associating metadata with them instead of repeating it for each cell.
"""
function cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "cell" => (expectation, "Sequenced single cells.")
end

"""
    cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique metacell each cell belongs to. All the cells in the same metacell are assumed to have "the same" (relevant)
biological state. This is the empty string if the cell does not belong to any metacell (is excluded or outlier).
"""
function cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "metacell") => (expectation, AbstractString, "The unique metacell each cell belongs to.")
end

"""
    cell_old_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique old metacell each cell belongs to. This just uses the renamed old metacells axis when computing a new set of
metacells for some data based on the old metacells.
"""
function cell_old_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "old_metacell") => (expectation, AbstractString, "The unique old metacell each cell belongs to.")
end

"""
    cell_new_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique new metacell each cell belongs to. This just uses the renamed new metacells axis when computing a new set of
metacells for some data based on the old metacells.
"""
function cell_new_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "new_metacell") => (expectation, AbstractString, "The unique new metacell each cell belongs to.")
end

"""
    metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of metacells, which are minimal-sized groups of cells for robust point estimates. Each metacell is considered
to be a robustly estimated point in the multi-dimensional manifold of cell states. Metacells may be very similar or very
different from each other depending on the data set and the manifold region.
"""
function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "metacell" => (expectation, "Minimal-sized groups of cells for robust point estimates.")
end

"""
    old_metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

Old metacells when computing new ones. This is the same as the regular metacells axis, but renamed to avoid confusion
when computing a new set of metacells for some data based on the old metacells.
"""
function old_metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "old_metacell" => (expectation, "Old metacells when computing new ones.")
end

"""
    new_metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

New metacells when improving old ones. This is the same as the regular metacells axis, but renamed to avoid confusion
when computing a new set of metacells for some data based on the old metacells.
"""
function new_metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "new_metacell" => (expectation, "New metacells when improving new ones.")
end

"""
    block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of blocks, which are distinct groups of metacells with "very close" estimated cell state. The metacells in each
block all have very close estimates of the skeleton gene expressions (maximal fold factor up to some maximal span).
"""
function block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "block" => (expectation, "Distinct groups of metacells with \"very close\" estimated cell state.")
end

"""
    type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of types, which are distinct named biological cell states. Types are convenient labels manually assigned to
large groups of cells, metacells, blocks, neighborhoods, etc. The resolution of the type labels depends on the data set
and the type of analysis. In particular, types are not typically associated with a specific biological cell state, but
rather with a set of related biological cell states (possibly along a gradient of such states).
"""
function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "type" => (expectation, "Distinct named biological cell states.")
end

"""
    gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

For each gene, we scale fold factors of each gene by multiplying with (1 - divergence) of the gene. In particular this
is used when considering the distance between gene expressions. Therefore genes with a non-zero divergence will require
a higher fold factor
"""
function gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "divergence") =>
        (expectation, StorageFloat, "Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.")
end

"""
    gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that distinguish between cell states. These genes have a significant expression level at some cell
state, as well as a significant range of expression across all cell states, so can be used to distinguish between cell
states. Non-marker genes are by definition not useful for such analysis, but marker genes aren't necessarily useful due
to other considerations.
"""
function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_marker") => (expectation, Bool, "A mask of genes that distinguish between cell states.")
end

"""
    gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The relative ranks of the marker genes. The more the gene distinguishes between different cell states, the better
(lower) rank it has. That is, 1 is for the "most" marker gene. Non-marker genes are given an extremely high rank (that
maximal the data type allows).
"""
function gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "marker_rank") => (expectation, StorageUnsigned, "The ralative ranks of the marker genes.")
end

"""
    gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are lateral to the biological behaviors of interest. These genes may satisfy all criteria for being
in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors of
interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological
behaviors of interest. They are therefore masked out during the analysis.
"""
function gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_lateral") =>
        (expectation, Bool, "A mask of genes that are lateral to the biological behaviors of interest.")
end

"""
    gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are excluded from consideration. These genes are completely unrelated to the biological behaviors
of interest. Not only that, they have strong and variable expression levels; enough to have an global impact on the
expression level of the rest of the genes - for example, mitochondrial genes. Such genes make it difficult to estimate
the relative expression level of genes between different cell states. Therefore, such genes aren't even counted in the
total UMIs of each cell.
"""
function gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_lateral") => (expectation, Bool, "A mask of genes that are totally excluded from the analysis.")
end

"""
    gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are forbidden from being used as skeleton factors. When searching for a set of genes that predict
the expression of the rest of the genes, we do not consider the forbidden genes.
"""
function gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_forbidden") =>
        (expectation, Bool, "A mask of genes that are forbidden from being used as skeleton genes.")
end

"""
    gene_is_uncorrelated_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are not correlated with other gene(s). We typically search for groups of genes that act together. Genes
that have no correlation with other genes aren't useful for this sort of analysis, even if they are marker genes.
"""
function gene_is_uncorrelated_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_uncorrelated") =>
        (expectation, Bool, "A mask of genes that are not correlated with other gene(s).")
end

"""
    gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are used to predict the values of the rest of the (covered) genes. The local programs that we use
to approximate the manifold are linear combinations of the skeleton genes. A linear combination of a small set of
(different) such programs is used to predict the rest of the genes.
"""
function gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_skeleton") =>
        (expectation, Bool, "A mask of genes that are used to predict the values of the rest of the (covered) genes.")
end

"""
    cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of cells that are excluded from consideration. This can be due to any number of reasons - doublets, too
low a number of UMIs, to high a percentage of excluded gene UMIs, etc.
"""
function cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("cell", "is_excluded") =>
        (expectation, Bool, "A mask of cells that are totally excluded from the analysis.")
end

"""
    gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are covered by the local linear program. When we approximate the manifold using these local linear
programs, we only model these genes, ignoring the rest. This allows the approximation to be unaffected by "irrelevant"
genes.
"""
function gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_covered") =>
        (expectation, Bool, "A mask of genes that are covered by the local linear program.")
end

"""
    cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the genes in each cell. This doesn't include the UMIs of genes excluded for any reason.
"""
function cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the genes in each cell.")
end

"""
    metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each metacell. This is used to estimate
the robustness of the gene fraction estimates.
"""
function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each metacell.",
    )
end

"""
    metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each metacell. This is used to estimate the robustness of the
gene fraction estimates.
"""
function metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each metacell.")
end

"""
    block_linear_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The root mean squared error of predicting the covered genes using the local linear models. This model computes principal components
from the skeleton genes and then uses these components to predict the values of all covered genes.
"""
function block_linear_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "linear_RMSE") => (
        expectation,
        StorageFloat,
        "The root mean squared error of predicting the covered genes using the linear model.",
    )
end

"""
    block_linear_XRMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The cross-validated root mean squared error of predicting the covered genes using the linear model. This is similar to
the RMSE, but is computed using cross-validation (training a linear model on a subset of the local metacells and
evaluating it on the rest). This will be higher than the RMSE; the difference is an indicator of the "overfitting" of
the linear model, which should be low relative to the total RMSE.
"""
function block_linear_XRMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "linear_XRMSE") => (
        expectation,
        StorageFloat,
        "The cross-validated root mean squared error of predicting the covered genes using the linear model.",
    )
end

"""
    block_n_principal_components_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of used principal components of the local linear model for each block. Since `Daf` can only deal with
fixed-sized axes, we create a "too large" principal components axis. Each block actually uses a smaller (different)
number of principal components that hopefully reflects the actual dimensionality of the manifold at that point. The
coefficients of the rest of the principal components are set to zero.
"""
function block_n_principal_components_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "n_principal_components") => (
        expectation,
        StorageUnsigned,
        "The number of used principal components of the local linear model for each block.",
    )
end

"""
    block_principal_component_is_used_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Whether each principal component is used by each block. This mask is `true` for the first used principal components
of each block and `false` for the rest. It is useful for constructing `Daf` queries that only return coefficients for
the used principal components for a specific block (e.g.,
`gene & is_skeleton / principal_component & is_used ; block = \$(block_name) : \$(block_name)_skeleton_coefficient`).
"""
function block_principal_component_is_used_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "principal_component", "is_used") =>
        (expectation, Bool, "Whether each principal component is used by each block.")
end

"""
    metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each metacell belongs to.
"""
function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "type") => (expectation, AbstractString, "The type each metacell belongs to.")
end

"""
    metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each metacell belongs to. All the metacells in the same block are assumed to be "very close" to each
other, so much so that it is excusable to treat them as "the same" cell state.
"""
function metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "block") => (expectation, AbstractString, "The unique block each metacell belongs to.")
end

"""
    old_metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each old metacell belongs to. This just uses the renamed old metacells axis when computing a new set of
metacells for some data based on the old metacells.
"""
function old_metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("old_metacell", "block") => (expectation, AbstractString, "The unique block each old metacell belongs to.")
end

"""
    new_metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each new metacell belongs to. This just uses the renamed old metacells axis when computing a new set of
metacells for some data based on the old metacells.

!!! note

    The old and new metacells are associated with the same set of blocks, at least initially. At some point you would
    compute brand new blocks for the new metacells, but as these are not based on the old blocks, there's no need to
    access both sets of blocks at the same time so no need for "old blocks" and "new blocks".
"""
function new_metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("new_metacell", "block") => (expectation, AbstractString, "The unique block each new metacell belongs to.")
end

"""
    type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A unique color for each type for graphs.
"""
function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("type", "color") => (expectation, AbstractString, "A unique color for each type for graphs.")
end

"""
    metacell_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated fraction of the UMIs of each gene in each metacell. We assume that each metacell is a sample of the
manifold, representing a real biological state, regardless to its distance to other metacells (subject to cleaning up
batch effects, purging doublets, and compensating for any other technical artifacts).
"""
function metacell_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "fraction") =>
        (expectation, StorageFloat, "The estimated fraction of the UMIs of each gene in each metacell.")
end

"""
    cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of UMIs collected for each gene for each cell. This is the "ground truth" everything else is built on. The
total number of UMIs per cell differs (sometimes wildly) based on the scRNA-seq technology and protocol.
"""
function cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "cell", "UMIs") =>
        (expectation, StorageUnsigned, "The number of UMIs collected for each gene for each cell.")
end

"""
    metacell_gene_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs used to estimate the fraction of each gene in each metacell. This is used to estimate the
robustness of the estimate. When computing fold factors, we require the total number of UMIs (from both compared
estimates) to be some minimum, and possibly adjust the fold factor according to some confidence level (assuming a
multinomial sampling distribution).
"""
function metacell_gene_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of each gene in each metacell.",
    )
end

"""
    metacell_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The relative fraction of covered genes. This ensures the sum of the fractions of the covered genes is one in each
metacell, to ignore the impact of any residual uncovered gene programs.
"""
function metacell_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "covered_fraction") =>
        (expectation, StorageFloat, "The relative fraction of covered genes.")
end

"""
    metacell_gene_scaled_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log of the relative fraction of covered genes, scaled by divergence. This computes the log base 2 of the
`metacell_gene_covered_fraction_matrix`, with some regularization, and then scales the result using the
[`gene_divergence_vector`](@ref) so that the range of expression of "wild" genes is reduced. The end result is what we
actually try to fit with local linear programs.
"""
function metacell_gene_scaled_log_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "scaled_log_covered_fraction") =>
        (expectation, StorageFloat, "The log of the relative fraction of covered genes, scaled by divergence.")
end

"""
    old_metacell_gene_scaled_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log of the relative fraction of covered genes, scaled by divergence. This just uses the renamed old metacells axis
when computing a new set of metacells for some data based on the old metacells.
"""
function old_metacell_gene_scaled_log_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "old_metacell", "scaled_log_covered_fraction") =>
        (expectation, StorageFloat, "The log of the relative fraction of covered genes, scaled by divergence.")
end

"""
    block_gene_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs used to estimate the fraction of each gene in each block. This is used to estimate the
robustness of the estimate. When computing fold factors, we require the total number of UMIs (from both compared
estimates) to be some minimum, and possibly adjust the fold factor according to some confidence level (assuming a
multinomial sampling distribution).
"""
function block_gene_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of each gene in each block.",
    )
end

"""
    block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the mask of nearby blocks in its linear environment. We expect the covered genes in the metacells in the
blocks of the environment to reasonably fit a local linear model (of the log of the expression).
"""
function block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_environment") =>
        (expectation, Bool, "For each block, the mask of nearby blocks in its linear environment.")
end

"""
    block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the mask of nearby blocks in its immediate neighborhood. The neighborhood consists of a small number of
very close blocks, which we use as the basis for evaluating linear approximations (of the log of the expression of the
genes).
"""
function block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_neighborhood") =>
        (expectation, Bool, "For each block, the mask of nearby blocks in its immediate neighborhood.")
end

"""
    block_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated fraction of the UMIs of each gene in each block. This assumes the metacells in each block are samples of
the "same" biological state, based on the fact that all these metacells have very similar expression levels for all
the global predictive transcription factors.
"""
function block_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "fraction") =>
        (expectation, StorageFloat, "The estimated fraction of the UMIs of each gene in each block.")
end

"""
    principal_component_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of principal components of local linear programs.
"""
function principal_component_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "principal_component" => (expectation, "A principal component of a local linear program.")
end

"""
    block_gene_mean_scaled_log_covered_fraction_matrix(
        expectation::ContractExpectation,
    )::Pair{VectorKey, DataSpecification}

The mean in the metacells of the block of the log of the relative fraction of covered genes, scaled by divergence. This is used
by the linear model; specifically, the coefficients are applied to the difference from this mean, and the results are
added back to this mean. This is zero for uncovered genes.
"""
function block_gene_mean_scaled_log_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "gene", "mean_scaled_log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The mean in the metacells of the block of the log of the relative fraction of covered genes, scaled by divergence.",
    )
end

"""
    block_principal_component_gene_skeleton_coefficient_tensor(
        expectation::ContractExpectation
    )::Pair{TensorKey, DataSpecification}

The coefficient of each skeleton gene for computing a principle component. This is zero for non-skeleton genes. The
coefficients can be either positive or negative. Combined with the covered genes coefficients this gives a local linear
model for estimating all the covered genes.
"""
function block_principal_component_gene_skeleton_coefficient_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "principal_component", "gene", "skeleton_coefficient") =>
        (expectation, StorageFloat, "The coefficient of each skeleton gene for computing a principle component.")
end

"""
    block_principal_component_gene_covered_coefficient_tensor(
        expectation::ContractExpectation
    )::Pair{TensorKey, DataSpecification}

The coefficient of each principal component for computing each covered gene. This is zero for uncovered genes. Combined
with the skeleton genes coefficients this gives a local linear model for estimating all the covered genes (*including*
the skeleton genes).
"""
function block_principal_component_gene_covered_coefficient_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "principal_component", "gene", "covered_coefficient") => (
        expectation,
        StorageFloat,
        "The coefficient of each principal component for computing each covered (non-skeleton) gene.",
    )
end

end  # module
