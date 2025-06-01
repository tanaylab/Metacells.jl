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
export block_block_mean_euclidean_skeleton_distance
export block_block_max_skeleton_fold_distance
export block_covered_UMIs_vector
export block_gene_base_covered_fraction_matrix
export block_gene_covered_fraction_matrix
export block_gene_environment_marker_rank_matrix
export block_gene_fraction_matrix
export block_gene_is_anchor_matrix
export block_gene_is_environment_marker_matrix
export block_gene_is_environment_skeleton_matrix
export block_gene_linear_covered_fraction_matrix
export block_gene_linear_fraction_matrix
export block_gene_log_covered_fraction_matrix
export block_gene_log_fraction_matrix
export block_gene_log_linear_covered_fraction_matrix
export block_gene_log_linear_fraction_matrix
export block_gene_log_scaled_linear_covered_fraction_matrix
export block_gene_log_scaled_linear_fraction_matrix
export block_gene_module_matrix
export block_gene_scaled_linear_covered_fraction_matrix
export block_gene_scaled_linear_fraction_matrix
export block_gene_UMIs_matrix
export block_module_anchor_gene_matrix
export block_module_base_covered_fraction_matrix
export block_module_gene_covered_coefficient_tensor
export block_module_is_found_matrix
export block_module_is_used_matrix
export block_module_n_genes_matrix
export block_module_total_UMIs_matrix
export block_modules_block_RMSE_vector
export block_modules_neighborhood_RMSE_vector
export block_modules_neighborhood_XRMSE_vector
export block_n_cells_vector
export block_n_environment_blocks_vector
export block_n_environment_cells_vector
export block_n_environment_metacells_vector
export block_n_found_modules_vector
export block_n_metacells_vector
export block_n_neighborhood_blocks_vector
export block_n_neighborhood_cells_vector
export block_n_neighborhood_metacells_vector
export block_n_used_modules_vector
export block_scaled_covered_UMIs_vector
export block_scaled_total_UMIs_vector
export block_total_UMIs_vector
export block_type_vector
export cell_axis
export cell_covered_UMIs_vector
export cell_excluded_UMIs_vector
export cell_gene_UMIs_matrix
export cell_is_excluded_vector
export cell_linear_covered_metacell_cross_entropy_vector
export cell_linear_covered_metacell_kl_divergence_vector
export cell_linear_metacell_cross_entropy_vector
export cell_linear_metacell_kl_divergence_vector
export cell_metacell_vector
export cell_mitochondrial_UMIs_vector
export cell_ribosomal_UMIs_vector
export cell_total_UMIs_vector
export cell_type_vector
export gene_axis
export gene_divergence_vector
export gene_is_covered_vector
export gene_is_excluded_vector
export gene_is_forbidden_vector
export gene_is_lateral_vector
export gene_is_marker_vector
export gene_is_regulator_vector
export gene_is_skeleton_vector
export gene_is_transcription_factor_vector
export gene_is_uncorrelated_vector
export gene_marker_rank_vector
export metacell_axis
export metacell_block_vector
export metacell_covered_UMIs_vector
export metacell_gene_approximated_covered_fraction_matrix
export metacell_gene_approximated_linear_covered_fraction_matrix
export metacell_gene_approximated_log_covered_fraction_matrix
export metacell_gene_approximated_log_linear_covered_fraction_matrix
export metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix
export metacell_gene_approximated_scaled_linear_covered_fraction_matrix
export metacell_gene_covered_fraction_matrix
export metacell_gene_fraction_matrix
export metacell_gene_geomean_fraction_matrix
export metacell_gene_linear_covered_fraction_matrix
export metacell_gene_linear_fraction_matrix
export metacell_gene_log_covered_fraction_matrix
export metacell_gene_log_fraction_matrix
export metacell_gene_log_geomean_fraction_matrix
export metacell_gene_log_linear_covered_fraction_matrix
export metacell_gene_log_linear_fraction_matrix
export metacell_gene_log_scaled_linear_covered_fraction_matrix
export metacell_gene_log_scaled_linear_fraction_matrix
export metacell_gene_scaled_linear_covered_fraction_matrix
export metacell_gene_scaled_linear_fraction_matrix
export metacell_gene_UMIs_matrix
export metacell_mean_cells_linear_covered_cross_entropy_vector
export metacell_mean_cells_linear_covered_kl_divergence_vector
export metacell_mean_cells_linear_cross_entropy_vector
export metacell_mean_cells_linear_kl_divergence_vector
export metacell_metacell_euclidean_skeleton_distance
export metacell_metacell_max_skeleton_fold_distance
export metacell_module_covered_fraction_matrix
export metacell_module_covered_UMIs_matrix
export metacell_module_fraction_matrix
export metacell_module_linear_covered_fraction_matrix
export metacell_module_linear_fraction_matrix
export metacell_module_log_covered_fraction_matrix
export metacell_module_log_fraction_matrix
export metacell_module_log_linear_covered_fraction_matrix
export metacell_module_log_linear_fraction_matrix
export metacell_module_log_scaled_linear_covered_fraction_matrix
export metacell_module_log_scaled_linear_fraction_matrix
export metacell_module_scaled_covered_UMIs_matrix
export metacell_module_scaled_linear_covered_fraction_matrix
export metacell_module_scaled_linear_fraction_matrix
export metacell_module_scaled_total_UMIs_matrix
export metacell_module_total_UMIs_matrix
export metacell_n_cells_vector
export metacell_scaled_covered_UMIs_vector
export metacell_scaled_total_UMIs_vector
export metacell_total_UMIs_vector
export metacell_type_vector
export module_axis
export type_axis
export type_color_vector

using DataAxesFormats

## Axes

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
    metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of metacells, which are minimal-sized groups of cells for robust point estimates. Each metacell is considered
to be a robustly estimated point in the multi-dimensional manifold of cell states. Metacells may be very similar or very
different from each other depending on the data set and the manifold region.
"""
function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "metacell" => (expectation, "Minimal-sized groups of cells for robust point estimates.")
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
    module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of local (per block) gene modules. There's no relationship between the same-name/index gene module in different
blocks; we just use this axis to be able to create matrices per-block-per-gene-module.
"""
function module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "module" => (expectation, "A local (per block) gene module.")
end

## Vectors

### Gene Vectors

"""
    gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are excluded from consideration. These genes are completely unrelated to the biological behaviors
of interest. Not only that, they have strong and variable expression levels; enough to have an global impact on the
expression level of the rest of the genes - for example, mitochondrial genes. Such genes make it difficult to estimate
the relative expression level of genes between different cell states. Therefore, such genes aren't even counted in the
total UMIs of each cell. Properties with "total" in their name ignore excluded genes, e.g. `total_UMIs` do not
include the UMIs of excluded genes.

This vector is created in a supervised way based on biological and technical considerations.
"""
function gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_excluded") =>
        (expectation, Bool, "A mask of genes that are totally excluded from the analysis.")
end

"""
    gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are lateral to the biological behaviors of interest. These genes may satisfy all criteria for being
in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors of
interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological
behaviors of interest. They are therefore masked out during the analysis.

This vector is created in a supervised way based on biological considerations.
"""
function gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_lateral") =>
        (expectation, Bool, "A mask of genes that are lateral to the biological behaviors of interest.")
end

"""
    gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that distinguish between cell states. These genes have a significant expression level at some cell
state, as well as a significant range of expression across all cell states, so can be used to distinguish between cell
states. Non-marker genes are by definition not useful for such analysis, but marker genes aren't necessarily useful due
to other considerations.

This vector is populated by [`identify_marker_genes!`](@ref Metacells.AnalyzeGenes.identify_marker_genes!).
"""
function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_marker") => (expectation, Bool, "A mask of genes that distinguish between cell states.")
end

"""
    gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The relative ranks of the marker genes. The more the gene distinguishes between different cell states, the better
(lower) rank it has. That is, 1 is for the "most" marker gene. Non-marker genes are given an extremely high rank (that
maximal the data type allows).

This vector is populated by [`rank_marker_genes!`](@ref Metacells.AnalyzeGenes.rank_marker_genes!).
"""
function gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "marker_rank") => (expectation, StorageUnsigned, "The ralative ranks of the marker genes.")
end

"""
    gene_is_uncorrelated_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are not correlated with other gene(s). We typically search for groups of genes that act together. Genes
that have no correlation with other genes aren't useful for this sort of analysis, even if they are marker genes.

This vector is populated by [`identify_uncorrelated_genes!`](@ref Metacells.AnalyzeGenes.identify_uncorrelated_genes!).
"""
function gene_is_uncorrelated_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_uncorrelated") =>
        (expectation, Bool, "A mask of genes that are not correlated with other gene(s).")
end

"""
    gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are covered by the local linear programs. When we approximate the manifold using these local linear
programs, we only model these genes, ignoring the rest. This allows the approximation to be unaffected by "irrelevant"
genes. The covered genes are typically some subset of the marker genes of the local region; in particular, they don't
include lateral genes. Properties that are computed using only the covered genes contain `_covered_` in their name,
e.g., `linear_covered_fraction` vs. `linear_fraction`.

This vector is populated by [`identify_covered_genes!`](@ref Metacells.AnalyzeGenes.identify_covered_genes!).
"""
function gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_covered") =>
        (expectation, Bool, "A mask of genes that are covered by the local linear program.")
end

"""
    gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are used to predict the values of the rest of the (covered) genes. We assume that knowing the
expression level of the skeleton genes is sufficient to reasonably estimate the expression level of the rest of the
(covered) genes. For example, if two metacells have "very close" expression level of all the skeleton genes, we assume
that these metacells are "very similar", without looking at the rest of the genes.

This vector is populated by [`identify_skeleton_genes!`](@ref Metacells.AnalyzeGenes.identify_skeleton_genes!).
"""
function gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_skeleton") =>
        (expectation, Bool, "A mask of genes that are used to predict the values of the rest of the (covered) genes.")
end

"""
    gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are forbidden from being used as skeleton factors. When searching for a set of genes that predict
the expression of the rest of the genes, we do not consider the forbidden genes.

This vector is created in a supervised way based on biological considerations.
"""
function gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_forbidden") =>
        (expectation, Bool, "A mask of genes that are forbidden from being used as skeleton genes.")
end

"""
    gene_is_regulator_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that regulate the expression level of other genes. These genes are expected to be at the core of the gene
programs that describe cell behaviors.

This vector is populated by [`fetch_regulators!`](@ref Metacells.AnalyzeGenes.fetch_regulators!).

!!! note

    These are **not** all the genes that bind to DNA.
"""
function gene_is_regulator_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_regulator") =>
        (expectation, Bool, "A mask of genes that regulate the expression level of other genes.")
end

"""
    gene_is_transcription_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that bind to the DNA.

This vector is populated by [`fetch_transcription_factors!`](@ref Metacells.AnalyzeGenes.fetch_transcription_factors!).
"""
function gene_is_transcription_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_transcription_factor") => (expectation, Bool, "A mask of genes that bind to the DNA.")
end

"""
    gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

For each gene, we scale the fractions of each gene by multiplying with `(1 - divergence)` of the gene. In particular
this is used when considering the distance between gene expressions. When scaling we also renormalize the result so the
fractions will sum to one again (for the chosen subset of the genes). Properties that take this scaling into account
have `scaled` in their name, e.g., `scaled_linear_fraction` vs. `linear_fraction`.

This vector is populated by [`compute_genes_divergence!`](@ref Metacells.AnalyzeGenes.compute_genes_divergence!).
"""
function gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "divergence") =>
        (expectation, StorageFloat, "Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.")
end

### Cell Vectors

#### Cell UMIs Vectors

"""
    cell_excluded_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the excluded genes in each cell.

This vector is created in a supervised way based on biological considerations.
"""
function cell_excluded_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "excluded_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the excluded genes in each cell.")
end

"""
    cell_mitochondrial_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the mitochondrial genes in each cell.

This vector is created in a supervised way, typically based on the naming convention that mitochondrial genes names
match the pattern `^MT-.*\$`.
"""
function cell_mitochondrial_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "mitochondrial_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the mitochondrial genes in each cell.")
end

"""
    cell_ribosomal_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the ribosomal genes in each cell.

This vector is created in a supervised way, typically based on the naming convention that mitochondrial genes names
match the pattern `^RP[LS]-.*\$`.
"""
function cell_ribosomal_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "ribosomal_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the ribosomal genes in each cell.")
end

"""
    cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the genes in each cell. This doesn't include the UMIs of genes excluded for any reason.

This vector is populated by [`compute_cells_total_UMIs!`](@ref Metacells.AnalyzeCells.compute_cells_total_UMIs!).
"""
function cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the genes in each cell.")
end

"""
    cell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the covered genes in each cell.

This vector is populated by [`compute_cells_covered_UMIs!`](@ref Metacells.AnalyzeCells.compute_cells_covered_UMIs!).
"""
function cell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the covered genes in each cell.")
end

#### Cell Metadata Vectors

"""
    cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of cells that are excluded from consideration. This can be due to any number of reasons - doublets, too low a
number of UMIs, to high a percentage of excluded gene UMIs, etc. Excluded cells are not grouped into metacells etc., so
they are ignored in all properties that sum UMIs per group (e.g., `total_UMIs` per metacell).

This vector is created in a supervised way based on biological and technical considerations.
"""
function cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("cell", "is_excluded") =>
        (expectation, Bool, "A mask of cells that are totally excluded from the analysis.")
end

"""
    cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique metacell each cell belongs to. All the cells in the same metacell are assumed to have "the same" (relevant)
biological state. This is the empty string if the cell does not belong to any metacell (is excluded or outlier).

This vector can be populated by any metacells-like algorithm (the original R metacells (1) algorithm, the Python
metacells (2) algorithm, other similar algorithms). It can also be populated by [`sharpen_metacells!`](@ref
Metacells.SharpenMetacells.sharpen_metacells!) which creates a hopefully-improved version of a given grouping into
metacells.
"""
function cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "metacell") => (expectation, AbstractString, "The unique metacell each cell belongs to.")
end

"""
    cell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each cell belongs to. This is typically deduced from the type of the metacell the cell belongs to,
unless we import the data from an annotated cell atlas.

This vector can be populated by [`compute_cells_types_by_metacells!`](@ref
Metacells.AnalyzeMetacells.compute_cells_types_by_metacells!) which assigns to each cell the type of the metacell it
belongs to.
"""
function cell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "type") => (expectation, AbstractString, "The type each cell belongs to.")
end

"""
    cell_linear_metacell_cross_entropy_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The cross entropy of the gene fractions of a cell given their fractions in the rest of the metacell cells.
"""
function cell_linear_metacell_cross_entropy_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "linear_metacell_cross_entropy") => (
        expectation,
        StorageFloat,
        "The cross entropy of the gene fractions of a cell given their fractions in the rest of the metacell cells.",
    )
end

"""
    cell_linear_covered_metacell_cross_entropy_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The cross entropy of the covered gene fractions of a cell given their fractions in the rest of the metacell cells.
"""
function cell_linear_covered_metacell_cross_entropy_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("cell", "linear_covered_metacell_cross_entropy") => (
        expectation,
        StorageFloat,
        "The cross entropy of the covered gene fractions of a cell given their fractions in the rest of the metacell cells.",
    )
end

"""
    cell_linear_metacell_kl_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The KL divergence of the gene fractions of a cell given their fractions in the rest of the metacell cells.
"""
function cell_linear_metacell_kl_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "linear_metacell_kl_divergence") => (
        expectation,
        StorageFloat,
        "The KL divergence of the gene fractions of a cell given their fractions in the rest of the metacell cells.",
    )
end

"""
    cell_linear_covered_metacell_kl_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The KL divergence of the covered gene fractions of a cell given their fractions in the rest of the metacell cells.
"""
function cell_linear_covered_metacell_kl_divergence_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("cell", "linear_covered_metacell_kl_divergence") => (
        expectation,
        StorageFloat,
        "The KL divergence of the covered gene fractions of a cell given their fractions in the rest of the metacell cells.",
    )
end

### Metacell Vectors

#### Metacell UMIs Vectors

"""
    metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each metacell. This is used to estimate
the robustness of the gene fraction estimates.

This vector is populated by [`compute_metacells_total_UMIs!`](@ref Metacells.AnalyzeMetacells.compute_metacells_total_UMIs!).
"""
function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each metacell.",
    )
end

"""
    metacell_scaled_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each metacell, scaled by divergence. This is
rounded to the nearest UMI.

This vector is populated by [`compute_metacells_scaled_total_UMIs!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_scaled_total_UMIs!).
"""
function metacell_scaled_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "scaled_total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each metacell, scaled by divergence.",
    )
end

"""
    metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each metacell.

This vector is populated by [`compute_metacells_covered_UMIs!`](@ref Metacells.AnalyzeMetacells.compute_metacells_covered_UMIs!).
"""
function metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each metacell.")
end

"""
    metacell_scaled_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each metacell, scaled by divergence. This is rounded to the nearest UMI.

This vector is populated by [`compute_metacells_scaled_covered_UMIs!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_scaled_covered_UMIs!).
"""
function metacell_scaled_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "scaled_covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each metacell, scaled by divergence.")
end

#### Metacell Counts Vectors

"""
    metacell_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of cells in each metacell.

This vector is populated by [`compute_metacells_n_cells!`](@ref Metacells.AnalyzeMetacells.compute_metacells_n_cells!).
"""
function metacell_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "n_cells") => (expectation, StorageUnsigned, "The number of cells in each metacell.")
end

#### Metacell Metadata Vectors

"""
    metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each metacell belongs to. All the metacells in the same block are assumed to be "very close" to each
other, so much so that it is excusable to treat them as "the same" cell state.

This vector is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "block") => (expectation, AbstractString, "The unique block each metacell belongs to.")
end

"""
    metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each metacell belongs to. This can be assigned in many ways, from manual annotation to using projection methods
on an annotated atlas. Often there are multiple type annotations due to different methods. By convention these alternate
type vectors are named `type.something` (e.g., `type.projected`, `type.manual`) and the simple `type` is used for "the"
type of each metacell.

This vector is created in a supervised way based on biological considerations.
"""
function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "type") => (expectation, AbstractString, "The type each metacell belongs to.")
end

"""
    metacell_mean_cells_linear_metacell_cross_entropy_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The mean cross entropy of the gene fractions of each metacell cell, given their fractions in the rest of the metacell
cells.
"""
function metacell_mean_cells_linear_cross_entropy_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "mean_cells_linear_cross_entropy") => (
        expectation,
        StorageFloat,
        "The mean cross entropy of the gene fractions of each metacell cell, given their fractions in the rest of the metacell cells.",
    )
end

"""
    metacell_mean_cells_linear_covered_metacell_cross_entropy_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The mean cross entropy of the covered gene fractions of each metacell cell, given their fractions in the rest of the
metacell cells.
"""
function metacell_mean_cells_linear_covered_cross_entropy_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "mean_cells_linear_covered_cross_entropy") => (
        expectation,
        StorageFloat,
        "The mean cross entropy of the covered gene fractions of each metacell cell, given their fractions in the rest of the metacell cells.",
    )
end

"""
    metacell_mean_cells_linear_metacell_kl_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The mean KL divergence of the gene fractions of each metacell cell, given their fractions in the rest of the metacell
cells.
"""
function metacell_mean_cells_linear_kl_divergence_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "mean_cells_linear_kl_divergence") => (
        expectation,
        StorageFloat,
        "The mean KL divergence of the gene fractions of each metacell cell, given their fractions in the rest of the metacell cells.",
    )
end

"""
    metacell_mean_cells_linear_covered_metacell_kl_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The mean KL divergence of the covered gene fractions of each metacell cell, given their fractions in the rest of the
metacell cells.
"""
function metacell_mean_cells_linear_covered_kl_divergence_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "mean_cells_linear_covered_kl_divergence") => (
        expectation,
        StorageFloat,
        "The mean KL divergence of the covered gene fractions of each metacell cell, given their fractions in the rest of the metacell cells.",
    )
end

### Block Vectors

#### Block UMIs Vectors

"""
    block_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each block. This is used to estimate
the robustness of the gene fraction estimates.

This vector is populated by [`compute_blocks_total_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_total_UMIs!).
"""
function block_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each block.",
    )
end

"""
    block_scaled_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each block, scaled by divergence. This is
rounded to the nearest UMI.

This vector is populated by [`compute_blocks_scaled_total_UMIs!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_scaled_total_UMIs!).
"""
function block_scaled_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "scaled_total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each block, scaled by divergence.",
    )
end

"""
    block_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each block.

This vector is populated by [`compute_blocks_covered_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_covered_UMIs!).
"""
function block_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each block.")
end

"""
    block_scaled_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each block, scaled by divergence.

This vector is populated by [`compute_blocks_scaled_covered_UMIs!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_scaled_covered_UMIs!).
"""
function block_scaled_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "scaled_covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each block, scaled by divergence.")
end

#### Block Counts Vectors

"""
    block_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of cells in the metacells in each block.

This vector is populated by [`compute_blocks_n_cells!`](@ref Metacells.AnalyzeBlocks.compute_blocks_n_cells!).
"""
function block_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_cells") => (expectation, StorageUnsigned, "The number of cells in the metacells in each block.")
end

"""
    block_n_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of metacells in each block.

This vector is populated by [`compute_blocks_n_metacells!`](@ref Metacells.AnalyzeBlocks.compute_blocks_n_metacells!).
"""
function block_n_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_metacells") => (expectation, StorageUnsigned, "The number of metacells in each block.")
end

"""
    block_n_neighborhood_blocks_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of blocks in the neighborhood centered at each block.

This vector is populated by [`compute_blocks_n_neighborhood_blocks!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_n_neighborhood_blocks!).
"""
function block_n_neighborhood_blocks_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_blocks") =>
        (expectation, StorageUnsigned, "The total number of blocks in the neighborhood centered at each block.")
end

"""
    block_n_environment_blocks_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of blocks in the environment centered at each block.

This vector is populated by [`compute_blocks_n_environment_blocks!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_n_environment_blocks!).
"""
function block_n_environment_blocks_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_environment_blocks") =>
        (expectation, StorageUnsigned, "The total number of blocks in the environment centered at each block.")
end

"""
    block_n_neighborhood_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of metacells in the blocks of the neighborhood centered at each block.

This vector is populated by [`compute_blocks_n_neighborhood_metacells!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_n_neighborhood_metacells!).
"""
function block_n_neighborhood_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_metacells") => (
        expectation,
        StorageUnsigned,
        "The total number of metacells in the blocks of the neighborhood centered at each block.",
    )
end

"""
    block_n_environment_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of metacells in the blocks of the environment centered at each block.

This vector is populated by [`compute_blocks_n_environment_metacells!`](@ref Metacells.AnalyzeBlocks.compute_blocks_n_environment_metacells!).
"""
function block_n_environment_metacells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_environment_metacells") => (
        expectation,
        StorageUnsigned,
        "The total number of metacells in the blocks of the environment centered at each block.",
    )
end

"""
    block_n_neighborhood_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of cells in the metacells in the blocks of the neighborhood centered at each block.

This vector is populated by [`compute_blocks_n_neighborhood_cells!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_n_neighborhood_cells!).
"""
function block_n_neighborhood_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_cells") => (
        expectation,
        StorageUnsigned,
        "The total number of cells in the metacells in the blocks of the neighborhood centered at each block.",
    )
end

"""
    block_n_environment_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of cells in the metacells in the blocks of the environment centered at each block.

This vector is populated by [`compute_blocks_n_environment_cells!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_n_environment_cells!).
"""
function block_n_environment_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_environment_cells") => (
        expectation,
        StorageUnsigned,
        "The total number of cells in the metacells in the blocks of the environment centered at each block.",
    )
end

#### Block Metadata Vectors

"""
    block_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each block belongs to. This is typically taken to be the most frequent metacell type in the block.

This vector is populated by [`compute_blocks_types!`](@ref Metacells.AnalyzeBlocks.compute_blocks_types!) based on the
metacell types in each block.
"""
function block_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "type") => (expectation, AbstractString, "The type each block belongs to.")
end

#### Gene Modules Analysis

"""
    block_n_found_modules_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of found gene modules for each block. We identify a different number of gene modules for each block. Due to
`Daf` limitations, we set the size  of the `modules` axis to the maximal number of gene modules in a block. The
coefficients for the extra gene modules for each block are set to zero.

This vector is populated by [`compute_blocks_n_found_modules!`](@ref
Metacells.AnalyzeModules.compute_blocks_n_found_modules!).
"""
function block_n_found_modules_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "n_found_modules") => (expectation, StorageUnsigned, "The number of found modules for each block.")
end

"""
    block_n_used_modules_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of used gene modules for each block. We disqualify some of the found gene modules, if they contain a high
fraction of lateral genes.

This vector is populated by [`compute_blocks_n_used_modules!`](@ref
Metacells.AnalyzeModules.compute_blocks_n_used_modules!).
"""
function block_n_used_modules_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "n_used_modules") => (expectation, StorageUnsigned, "The number of used modules for each block.")
end

"""
    block_modules_block_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block,
evaluated at the block metacells.

This vector is populated by
[`compute_blocks_modules_block_RMSE!`](@ref Metacells.AnalyzeApproximation.compute_blocks_modules_block_RMSE!).
"""
function block_modules_block_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "modules_block_RMSE") => (
        expectation,
        StorageFloat,
        "The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block, evaluated at the block metacells.",
    )
end

"""
    block_modules_neighborhood_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block,
evaluated at the neighborhood metacells.

This vector is populated by
[`compute_blocks_modules_neighborhood_RMSE!`](@ref Metacells.AnalyzeApproximation.compute_blocks_modules_neighborhood_RMSE!).
"""
function block_modules_neighborhood_RMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "modules_neighborhood_RMSE") => (
        expectation,
        StorageFloat,
        "The root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block, evaluated at the neighborhood metacells.",
    )
end

"""
    block_modules_neighborhood_XRMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The cross-validated root mean squared error of predicting the covered genes using the gene modules in the vicinity of
the block, evaluated at the neighborhood metacells. This is similar to the `modules_neighborhood_RMSE`, but is computed
using cross-validation (training the model on a subset of the local metacells, and evaluating it on the rest). This will
be higher than the RMSE; the difference is an indicator of the overfitting of the linear model, which should be low
relative to the total RMSE. The cross-validated value is the one used to determine the environment size.

This vector is populated by [`compute_approximation!`](@ref Metacells.ComputeApproximation.compute_approximation!).
"""
function block_modules_neighborhood_XRMSE_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "modules_neighborhood_XRMSE") => (
        expectation,
        StorageFloat,
        "The cross-validated root mean squared error of predicting the covered genes using the gene modules in the vicinity of the block, evaluated at the neighborhood metacells.",
    )
end

### Type Vectors

"""
    type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A unique color for each type for graphs.

This vector is created in a supervised way based on biological considerations and conventions (e.g., red blood cells are
often given some red color). It is also possible to use `Chameleon` to automatically assign colors to cell types based
on gene expression levels.
"""
function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("type", "color") => (expectation, AbstractString, "A unique color for each type for graphs.")
end

## Matrices

### Cell Matrices

"""
    cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of UMIs collected for each gene for each cell. This is the "ground truth" everything else is built on. The
total number of UMIs is different (sometimes wildly) in each cell, based on the scRNA-seq technology and protocol.

This data is obtained from scRNA-seq experiments.
"""
function cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "cell", "UMIs") =>
        (expectation, StorageUnsigned, "The number of UMIs collected for each gene for each cell.")
end

### Metacell Matrices

#### Metacell UMIs

"""
    metacell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs used to estimate the fraction of each gene in each metacell. This can be used to estimate the
robustness of the fraction.

This matrix is populated by [`compute_metacells_genes_UMIs!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_UMIs!).
"""
function metacell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of each gene in each metacell.",
    )
end

#### Metacell Geomean Fractions

"""
    metacell_gene_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

An estimated geomean fraction of the UMIs of each gene in each metacell. We use geomean in an attempt to combat the
disproportionate effect of a few cells with very high gene expression ("bursty" genes), and then normalizes the
fractions to sum to one. While effective, this has the unfortunate effect of inflating the value of weak genes.

This matrix is populated by [`compute_metacells_genes_geomean_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_geomean_fractions!).
"""
function metacell_gene_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "geomean_fraction") =>
        (expectation, StorageFloat, "The estimated geomean fraction of the UMIs of each gene in each metacell.")
end

"""
    metacell_gene_log_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated geomean fraction of the UMIs of each gene in each metacell. This adds some gene fraction
regularization to deal with zero fractions. Using the log makes it easier to visualize, and the difference between log
values (the "fold factor", log base 2 of the ratio between the expression levels) is a good measure of difference
between gene expression levels.

This matrix is populated by [`compute_metacells_genes_log_geomean_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_geomean_fractions!).
"""
function metacell_gene_log_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_geomean_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated geomean fraction of the UMIs of each gene in each metacell.",
    )
end

#### Metacell Linear Fractions

"""
    metacell_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

An estimated linear fraction of the UMIs of each gene in each metacell, out of the total UMIs. This is the "best"
estimate assuming multinomial sampling noise. However, this is sensitive to a few cells with very high expression levels
("bursty" genes), as these impact the denominator.

This matrix is populated by [`compute_metacells_genes_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_linear_fractions!).
"""
function metacell_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each gene in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_gene_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of each gene in each metacell, out of the total UMIs, scaled by divergence.

This matrix is populated by [`compute_metacells_genes_scaled_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_scaled_linear_fractions!).
"""
function metacell_gene_scaled_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each gene in each metacell, scaled by divergence.",
    )
end

"""
    metacell_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of each gene in each metacell, out of the total UMIs. This
adds some gene fraction regularization to deal with zero fractions. Using the log makes it easier to visualize, and the
difference between log values (the "fold factor", log base 2 of the ratio between the expression levels) is a good
measure of difference between gene expression levels.

This matrix is populated by [`compute_metacells_genes_log_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_linear_fractions!).
"""
function metacell_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of each gene in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_gene_log_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the scaled estimated linear fraction of the UMIs of each gene in each metacell. This adds some
gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_genes_log_scaled_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_scaled_linear_fractions!).
"""
function metacell_gene_log_scaled_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the scaled estimated linear fraction of the UMIs of each gene in each metacell.",
    )
end

#### Metacell Linear Covered Fractions

"""
    metacell_gene_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

An estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs. By
considering only the covered genes this avoid the impact of highly-expressed lateral genes (e.g., cell cycle).

This matrix is populated by [`compute_metacells_genes_linear_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_linear_covered_fractions!).
"""
function metacell_gene_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs.",
    )
end

"""
    metacell_gene_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs, scaled
by divergence. By considering only the covered genes this avoid the impact of highly-expressed lateral genes (e.g., cell
cycle).

This matrix is populated by [`compute_metacells_genes_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_scaled_linear_covered_fractions!).
"""
function metacell_gene_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each gene in each metacell, scaled by divergence.",
    )
end

"""
    metacell_gene_log_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total
covered UMIs. By considering only the covered genes this avoid the impact of highly-expressed lateral genes (e.g., cell
cycle).

This matrix is populated by [`compute_metacells_genes_log_linear_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_linear_covered_fractions!).
"""
function metacell_gene_log_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each metacell.",
    )
end

"""
    metacell_gene_log_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total
covered UMIs, scaled by divergence. By considering only the covered genes this avoid the impact of highly-expressed
lateral genes (e.g., cell cycle).

This matrix is populated by [`compute_metacells_genes_log_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_scaled_linear_covered_fractions!).
"""
function metacell_gene_log_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs, scaled by divergence.",
    )
end

#### Metacell Virtual Fractions

"""
    metacell_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each gene in each metacell. This is a "virtual" matrix, that is, generic code
uses it and expects an `adapter` to map it to one of the concrete `metacell_gene_*_fraction_matrix` matrices.
"""
function metacell_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "fraction") =>
        (expectation, StorageFloat, "Some estimated fraction of the UMIs of each gene in each metacell.")
end

"""
    metacell_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each covered gene in each metacell. This is a "virtual" matrix, that is, generic
code uses it and expects an `adapter` to map it to one of the concrete `metacell_gene_*_covered_fraction_matrix`
matrices.
"""
function metacell_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "covered_fraction") =>
        (expectation, StorageFloat, "Some estimated fraction of the UMIs of each covered gene in each metacell.")
end

"""
    metacell_gene_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the UMIs of each gene in each metacell. This is a "virtual" matrix,
that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`metacell_gene_log_*_fraction_matrix` matrices.
"""
function metacell_gene_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated fraction of the UMIs of each gene in each metacell.",
    )
end

"""
    metacell_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the UMIs of each covered gene in each metacell. This is a "virtual"
matrix, that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`metacell_gene_log_*_covered_fraction_matrix` matrices.
"""
function metacell_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated fraction of the UMIs of each covered gene in each metacell.",
    )
end

#### Metacell Distances

"""
    metacell_metacell_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The maximal fold factor between skeleton genes between the metacells. This is a symmetric distances matrix with zeros at
the diagonal. The exact semantics depend on the estimated fractions used to compute the distances.

This matrix is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function metacell_metacell_max_skeleton_fold_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "max_skeleton_fold_distance") =>
        (expectation, StorageFloat, "The maximal fold factor between skeleton genes between the metacells.")
end

"""
    metacell_metacell_euclidean_skeleton_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

TODOX

This matrix is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function metacell_metacell_euclidean_skeleton_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "euclidean_skeleton_distance") => (expectation, StorageFloat, "TODOX.")
end

#### Gene Modules Analysis

##### Gene Modules UMIs

"""
    metacell_module_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the genes of each module in each metacell.

This matrix is populated by [`compute_metacells_modules_total_UMIs!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_total_UMIs!).
"""
function metacell_module_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total UMIs of the genes of each module in each metacell.")
end

"""
    metacell_module_scaled_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the genes of each module in each metacell, scaled by divergence.

This matrix is populated by [`compute_metacells_modules_scaled_total_UMIs!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_scaled_total_UMIs!).
"""
function metacell_module_scaled_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "scaled_total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total UMIs of the genes of each module in each metacell, scaled by divergence.",
    )
end

"""
    metacell_module_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the covered genes of each module in each metacell.

This matrix is populated by [`compute_metacells_modules_covered_UMIs!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_covered_UMIs!).
"""
function metacell_module_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total UMIs of the covered genes of each module in each metacell.")
end

"""
    metacell_module_scaled_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the covered genes of each module in each metacell, scaled by divergence.

This matrix is populated by [`compute_metacells_modules_scaled_covered_UMIs!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_scaled_covered_UMIs!).
"""
function metacell_module_scaled_covered_UMIs_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "scaled_covered_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total UMIs of the covered genes of each module in each metacell, scaled by divergence.",
    )
end

##### Gene Modules Linear Fractions

"""
    metacell_module_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs.

This matrix is populated by [`compute_metacells_modules_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_linear_fractions!).
"""
function metacell_module_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_module_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs, scaled by divergence.

This matrix is populated by [`compute_metacells_modules_scaled_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_scaled_linear_fractions!).
"""
function metacell_module_scaled_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs, scaled by divergence.",
    )
end

"""
    metacell_module_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the
total UMIs. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_modules_log_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_log_linear_fractions!).
"""
function metacell_module_log_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_module_log_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the
total UMIs, scaled by divergence. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_modules_log_scaled_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_log_scaled_linear_fractions!).
"""
function metacell_module_log_scaled_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "log_scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of the genes of each module in each metacell, out of the total UMIs, scaled by divergence.",
    )
end

##### Gene Modules Linear Covered Fractions

"""
    metacell_module_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs.

This matrix is populated by [`compute_metacells_modules_linear_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_linear_covered_fractions!).
"""
function metacell_module_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs.",
    )
end

"""
    metacell_module_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered
UMIs, scaled by divergence.

This matrix is populated by [`compute_metacells_modules_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_scaled_linear_covered_fractions!).
"""
function metacell_module_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs, scaled by divergence.",
    )
end

"""
    metacell_module_log_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of
the total covered UMIs. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_modules_log_linear_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_log_linear_covered_fractions!).
"""
function metacell_module_log_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "log_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs.",
    )
end

"""
    metacell_module_log_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of
the total covered UMIs, scaled by divergence. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_modules_log_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_metacells_modules_log_scaled_linear_covered_fractions!).
"""
function metacell_module_log_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "module", "log_scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of the covered genes of each module in each metacell, out of the total covered UMIs, scaled by divergence.",
    )
end

##### Gene Modules Approximated Covered Fractions

"""
    metacell_gene_approximated_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model.

This matrix is populated by [`compute_metacells_genes_approximated_linear_covered_fractions!`](@ref
Metacells.AnalyzeApproximation.compute_metacells_genes_approximated_linear_covered_fractions!).
"""
function metacell_gene_approximated_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "gene", "approximated_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model.",
    )
end

"""
    metacell_gene_approximated_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of each covered gene UMIs out of the total covered UMIs, scaled by divergence, as approximated by the linear model.

This matrix is populated by [`compute_metacells_genes_approximated_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeApproximation.compute_metacells_genes_approximated_scaled_linear_covered_fractions!).
"""
function metacell_gene_approximated_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "gene", "approximated_scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of each covered gene UMIs out of the total covered UMIs, scaled by divergence, as approximated by the linear model.",
    )
end

"""
    metacell_gene_approximated_log_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the
linear model. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_genes_approximated_log_linear_covered_fractions!`](@ref
Metacells.AnalyzeApproximation.compute_metacells_genes_approximated_log_linear_covered_fractions!).
"""
function metacell_gene_approximated_log_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "gene", "approximated_log_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model.",
    )
end

"""
    metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs, scaled by divergence, as
approximated by the linear model. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_genes_approximated_log_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeApproximation.compute_metacells_genes_approximated_log_scaled_linear_covered_fractions!).
"""
function metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "gene", "approximated_log_scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of each covered gene UMIs out of the total covered UMIs, scaled by divergence, as approximated by the linear model.",
    )
end

#### Gene Modules Virtual Fractions

"""
    metacell_module_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each gene in each module. This is a "virtual" matrix, that is, generic code
uses it and expects an `adapter` to map it to one of the concrete `metacell_module_*_fraction_matrix` matrices.
"""
function metacell_module_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("module", "metacell", "fraction") =>
        (expectation, StorageFloat, "Some estimated fraction of the UMIs of each module in each metacell.")
end

"""
    metacell_module_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each covered gene in each module. This is a "virtual" matrix, that is, generic
code uses it and expects an `adapter` to map it to one of the concrete `metacell_module_*_covered_fraction_matrix`
matrices.
"""
function metacell_module_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("module", "metacell", "covered_fraction") =>
        (expectation, StorageFloat, "Some estimated covered fraction of the UMIs of each module in each metacell.")
end

"""
    metacell_module_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the UMIs of each module in each metacell. This is a "virtual" matrix,
that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`metacell_module_log_*_fraction_matrix` matrices.
"""
function metacell_module_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("module", "metacell", "log_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated linear fraction of the UMIs of each module in each metacell.",
    )
end

"""
    metacell_module_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the covered UMIs of each module in each metacell. This is a "virtual" matrix,
that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`metacell_module_log_*_covered_fraction_matrix` matrices.
"""
function metacell_module_log_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("module", "metacell", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated linear fraction of the covered UMIs of each module in each metacell.",
    )
end

"""
    metacell_gene_approximated_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model. This is a
"virtual" matrix, that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`metacell_gene_approximated_*_fraction_matrix` matrices.
"""
function metacell_gene_approximated_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "approximated_covered_fraction") => (
        expectation,
        StorageFloat,
        "The fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model.",
    )
end

"""
    metacell_gene_approximated_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear
model. This is a "virtual" matrix, that is, generic code uses it and expects an `adapter` to map it to one of the
concrete `metacell_gene_approximated_log_*_fraction_matrix` matrices.
"""
function metacell_gene_approximated_log_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("module", "metacell", "approximated_log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the fraction of each covered gene UMIs out of the total covered UMIs as approximated by the linear model.",
    )
end

### Block Matrices

#### Block Distances

"""
    block_block_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The maximal fold factor between skeleton genes between the metacells of the blocks. This is a symmetric distances matrix. The
diagonal contains the maximal fold factor between skeleton genes between metacells in a single block, which is only zero
for single-metacell blocks.

This matrix is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function block_block_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "max_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The maximal fold factor between skeleton genes between the metacells of the blocks.",
    )
end

"""
    block_block_mean_euclidean_skeleton_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

TODOX

This matrix is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function block_block_mean_euclidean_skeleton_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "mean_euclidean_skeleton_distance") => (expectation, StorageFloat, "TODOX.")
end

#### Block Vicinities

"""
    block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the (column) mask of nearby blocks in its immediate neighborhood. The neighborhood consists of a small
number of very close blocks, which we use as the basis for evaluating local approximations of the manifold.

This matrix is populated by [`compute_blocks_is_in_neighborhood!`](@ref
Metacells.ComputeBlocks.compute_blocks_is_in_neighborhood!).
"""
function block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_neighborhood") =>
        (expectation, Bool, "For each block, the (column) mask of nearby blocks in its immediate neighborhood.")
end

"""
    block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the (column) mask of nearby blocks in its linear environment. The environment consists of a larger
region than the neighborhood, which we use as the basis for computing local approximations of the manifold.

This matrix is populated by [`compute_blocks_is_in_environment!`](@ref
Metacells.ComputeBlocks.compute_blocks_is_in_environment!).
"""
function block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_environment") =>
        (expectation, Bool, "For each block, the (column) mask of nearby blocks in its expanded environment.")
end

#### Block UMIs

"""
    block_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs used to estimate the fraction of each gene in each block. This can used to estimate the
robustness of the fraction.

This matrix is populated by [`compute_blocks_genes_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_genes_UMIs!).
"""
function block_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of each gene in each block.",
    )
end

#### Block Linear Fractions

"""
    block_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

An estimated linear fraction of the UMIs of each gene out of the total UMIs in each block.

This matrix is populated by [`compute_blocks_genes_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_linear_fractions!).
"""
function block_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each gene out of the total UMIs in each block.",
    )
end

"""
    block_gene_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of each gene out of the total UMIs in each block, scaled by divergence.

This matrix is populated by [`compute_blocks_genes_scaled_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_scaled_linear_fractions!).
"""
function block_gene_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each gene out of the total UMIs in each block, scaled by divergence.",
    )
end

"""
    block_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of each gene out of the total UMIs in each block.

This matrix is populated by [`compute_blocks_genes_log_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_linear_fractions!).
"""
function block_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of each gene out of the total UMIs in each block.",
    )
end

"""
    block_gene_log_scaled_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the scaled estimated linear fraction of the UMIs of each gene out of the total UMIs in each block,
scaled by divergence.

This matrix is populated by [`compute_blocks_genes_log_scaled_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_scaled_linear_fractions!).
"""
function block_gene_log_scaled_linear_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_scaled_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the scaled estimated linear fraction of the UMIs of each gene out of the total UMIs in each block, scaled by divergence.",
    )
end

#### Block Linear Covered Fractions

"""
    block_gene_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

An estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.

This matrix is populated by [`compute_blocks_genes_linear_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_linear_covered_fractions!).
"""
function block_gene_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.",
    )
end

"""
    block_gene_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block, scaled by divergence.

This matrix is populated by [`compute_blocks_genes_log_linear_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_scaled_linear_covered_fractions!).
"""
function block_gene_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block, scaled by divergence.",
    )
end

"""
    block_gene_log_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.

This matrix is populated by [`compute_blocks_genes_log_linear_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_linear_covered_fractions!).
"""
function block_gene_log_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.",
    )
end

"""
    block_gene_log_scaled_linear_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the scaled estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.

This matrix is populated by [`compute_blocks_genes_log_scaled_linear_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_scaled_linear_covered_fractions!).
"""
function block_gene_log_scaled_linear_covered_fraction_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_scaled_linear_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the scaled estimated linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.",
    )
end

#### Block Virtual Fractions

"""
    block_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each gene in each block. This is a "virtual" matrix, that is, generic code
uses it and expects an `adapter` to map it to one of the concrete `block_gene_*_fraction_matrix` matrices.
"""
function block_gene_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "fraction") =>
        (expectation, StorageFloat, "Some estimated fraction of the UMIs of each gene in each block.")
end

"""
    block_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Some estimated fraction of the UMIs of each covered gene in each block. This is a "virtual" matrix, that is, generic
code uses it and expects an `adapter` to map it to one of the concrete `block_gene_*_covered_fraction_matrix` matrices.
"""
function block_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "covered_fraction") =>
        (expectation, StorageFloat, "Some estimated fraction of the UMIs of each covered gene in each block.")
end

"""
    block_gene_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the UMIs of each gene in each block. This is a "virtual" matrix,
that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`block_gene_log_*_fraction_matrix` matrices.
"""
function block_gene_log_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated linear fraction of the UMIs of each gene in each block.",
    )
end

"""
    block_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of some estimated linear fraction of the UMIs of each covered gene in each block. This is a "virtual"
matrix, that is, generic code uses it and expects an `adapter` to map it to one of the concrete
`block_gene_log_*_covered_fraction_matrix` matrices.
"""
function block_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of some estimated linear fraction of the UMIs of each covered gene in each block.",
    )
end

#### Block Gene Masks

"""
    block_gene_is_environment_marker_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

A mask of genes that distinguish between cell states in each environment.

This matrix is populated by [`compute_blocks_genes_is_environment_markers!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_is_environment_markers!).
"""
function block_gene_is_environment_marker_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "is_environment_marker") =>
        (expectation, Bool, "A mask of genes that distinguish between cell states in each environment.")
end

"""
    block_gene_environment_marker_rank_matrix(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The relative ranks of the marker genes in each environment.

This vector is populated by [`compute_blocks_is_in_environment!`](@ref Metacells.ComputeBlocks.compute_blocks_is_in_environment!).
"""
function block_gene_environment_marker_rank_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "gene", "environment_marker_rank") =>
        (expectation, StorageUnsigned, "The ralative ranks of the marker genes in each environment.")
end

#### Gene Modules Analysis

##### Gene Modules Counts

"""
    block_module_n_genes_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of genes in each gene module in each block. This is zero for the extra (not-found) gene modules.

This matrix is populated by [`compute_blocks_modules_n_genes!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_n_genes!).
"""
function block_module_n_genes_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "n_genes") =>
        (expectation, StorageUnsigned, "The number of genes in each gene module in each block.")
end

"""
    block_module_n_covered_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of covered genes in each gene module in each block. This is zero for the extra (not-found) gene modules.

This matrix is populated by [`compute_blocks_modules_n_covered!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_n_covered!).
"""
function block_module_n_covered_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "n_covered") =>
        (expectation, StorageUnsigned, "The number of covered genes in each gene module in each block.")
end

##### Gene Modules UMIs

"""
    block_module_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the genes in each gene module in each block.

This matrix is populated by [`compute_blocks_modules_total_UMIs!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_total_UMIs!).
"""
function block_module_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total UMIs of the genes in each gene module in each block.")
end

"""
    block_module_scaled_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the genes in each gene module in each block, scaled by divergence.

This matrix is populated by [`compute_blocks_modules_scaled_covered_UMIs!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_scaled_covered_UMIs!).
"""
function block_module_scaled_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "scaled_total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total UMIs of the genes in each gene module in each bloc, scaled by divergence.",
    )
end

"""
    block_module_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the covered genes in each gene module in each block.

This matrix is populated by [`compute_blocks_modules_covered_UMIs!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_covered_UMIs!).
"""
function block_module_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total UMIs of the covered genes in each gene module in each block.")
end

"""
    block_module_scaled_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total UMIs of the covered genes in each gene module in each block, scaled by divergence.

This matrix is populated by [`compute_blocks_modules_scaled_covered_UMIs!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_scaled_covered_UMIs!).
"""
function block_module_scaled_covered_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "scaled_covered_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total UMIs of the covered genes in each gene module in each block, scaled by divergence.",
    )
end

##### Gene Modules Masks

"""
    block_gene_module_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The gene module each gene belongs to in the environment of each block. This is the empty string for genes that do not
belong to any module.

This matrix is populated by [`compute_blocks_modules!`](@ref
Metacells.ComputeModules.compute_blocks_modules!).
"""
function block_gene_module_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "module") =>
        (expectation, AbstractString, "The gene module each gene belongs to in the environment of each block.")
end

"""
    block_module_anchor_gene_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The anchor gene for each gene module for each block. This is the transcription factor that all the rest of the genes
in the module are associated with.

This matrix is populated by [`compute_blocks_modules!`](@ref
Metacells.ComputeModules.compute_blocks_modules!).
"""
function block_module_anchor_gene_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "gene.anchor") =>
        (expectation, AbstractString, "The anchor gene for each gene module for each block.")
end

"""
    block_gene_is_anchor_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Whether a gene is the anchor for one of the modules for each block.
"""
function block_gene_is_anchor_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "is_anchor") =>
        (expectation, Bool, "Whether a gene is the anchor for one of the modules for each block.")
end

"""
    block_module_is_found_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Whether each gene module is found in each block. This mask is `true` for the first `n_found_modules` of each block
and `false` for the rest. Due to `Daf` limitations, we create a `module` axis with the maximal number of gene
modules in each block, so only the first `n_found_modules` of these actually contain genes.

This matrix is populated by [`compute_blocks_modules!`](@ref
Metacells.ComputeModules.compute_blocks_modules!).
"""
function block_module_is_found_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "is_found") => (expectation, Bool, "Whether each gene module is found in each block.")
end

##### Gene Modules Model

"""
    block_gene_base_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The covered base fractions used by the gene module analysis. The models actually predicts the difference of the genes
from these levels, based on the difference in the fraction of the gene modules from their base levels,
using the `covered_coefficient` matrix.

This matrix is populated by [`compute_blocks_is_in_environment!`](@ref
Metacells.ComputeBlocks.compute_blocks_is_in_environment!).
"""
function block_gene_base_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "gene", "base_covered_fraction") =>
        (expectation, StorageFloat, "The base covered fractions used by the modules analysis.")
end

"""
    block_module_base_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The base covered fractions used by the gene modules analysis. The model computes the difference of each gene module
fraction in each metacell (or cell) from these levels, then uses them to predict the difference from their
base levels, using the `covered_coefficient` matrix.

This matrix is populated by [`compute_metacells_genes_approximated_linear_covered_fractions!`](@ref
Metacells.AnalyzeApproximation.compute_metacells_genes_approximated_linear_covered_fractions!).

!!! note

    The model only considers the fractions of the covered genes in each gene module. Lateral genes in the modules are
    ignored.
"""
function block_module_base_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "module", "base_covered_fraction") =>
        (expectation, StorageFloat, "The base covered fractions used by the gene modules analysis.")
end

"""
    block_module_gene_covered_coefficient_tensor(
        expectation::ContractExpectation
    )::Pair{TensorKey, DataSpecification}

The coefficient of each gene module for computing each covered gene. This is zero for unused modules and uncovered
genes. Combined with the `module_fraction` and the `covered_coefficient` matrices, this gives a local linear model for
estimating all the covered genes (including the modules genes).

This matrix is populated by [`compute_approximation!`](@ref Metacells.ComputeApproximation.compute_approximation!).
"""
function block_module_gene_covered_coefficient_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "module", "gene", "covered_coefficient") =>
        (expectation, StorageFloat, "The coefficient of each gene module for computing each covered gene.")
end

end  # module
