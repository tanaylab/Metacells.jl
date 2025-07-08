"""
Functions for defining a `Contract` for a metacells `Daf` data set `@computation`. This also serves as a vocabulary
describing the data we keep when doing metacells based analysis, which is a solid basis for understanding what is going
on. As Fred Brooks said: "Show me your flowcharts and conceal your tables, and I shall continue to be mystified. Show me
your tables, and I won’t usually need your flowcharts; they’ll be obvious.".

We often look at the log base 2 of the gene expression level. We use a regularization factor of 1e-5 to avoid a log of
zero. This means the minimal log value is -16.6 and the maximal value for a well-behaved gene (~0.1% of the total UMIs
of the cell) is -10. We typically look at the range -16 "basically no expression" to -10 "respectably high expression".

In the descriptions below, "fold factor" refers to the log base 2 of the ratio between expression levels. That is, we
look at an interesting fold factor of ~6 between "basically no expression" and "respectably high expression". A fold
factor (difference in logs) of ~3 (that is, ~8x expression level) is typically taken to be "significant" for a single
gene between cells by the metacell algorithm. However lower fold factors may be significant, especially as the number of
UMIs used to estimate the expression levels are larger (e.g., between expression levels of metacells). As a quick and
useful heuristic, we don't consider fold factors between metacell genes to be meaningful unless the sum of the gene UMIs
in both compared metacells is too low (below 40).

Genes that are higher than -10 (that is, more than 0.1%) are "bursty", e.g. Hemoglobins and Cytokines. Sometimes such
genes take more than 10% of the total UMIs of a cell! TODO: Special handling for such genes, as they effect the
denominator for linear fraction estimations and overwhelm least-square approximations.
"""
module Contracts

export block_axis
export block_block_is_in_environment_matrix
export block_block_is_in_neighborhood_matrix
export block_block_max_skeleton_fold_distance
export block_block_mean_euclidean_skeleton_distance
export block_covered_UMIs_vector
export block_gene_covered_fraction_matrix
export block_gene_is_correlated_with_metacells
export block_gene_is_environment_marker_matrix
export block_gene_linear_fraction_matrix
export block_gene_log_covered_fraction_matrix
export block_gene_log_linear_fraction_matrix
export block_gene_module_matrix
export block_gene_neighborhood_correlation_matrix
export block_gene_UMIs_matrix
export block_mean_covered_neighborhood_correlation_vector
export block_mean_modules_neighborhood_correlation_vector
export block_mean_neighborhood_correlation_vector
export block_metacell_module_environment_covered_fraction_tensor
export block_metacell_module_environment_linear_fraction_tensor
export block_metacell_module_environment_log_covered_fraction_tensor
export block_metacell_module_environment_log_linear_fraction_tensor
export block_metacell_module_environment_total_UMIs_tensor
export block_module_n_genes_matrix
export block_module_n_skeletons_matrix
export block_n_cells_vector
export block_n_environment_blocks_vector
export block_n_environment_cells_vector
export block_n_environment_metacells_vector
export block_n_metacells_vector
export block_n_neighborhood_blocks_vector
export block_n_neighborhood_cells_vector
export block_n_neighborhood_metacells_vector
export block_total_UMIs_vector
export block_type_vector
export cell_axis
export cell_covered_UMIs_vector
export cell_excluded_UMIs_vector
export cell_gene_UMIs_matrix
export cell_is_excluded_vector
export cell_metacell_vector
export cell_mitochondrial_UMIs_vector
export cell_ribosomal_UMIs_vector
export cell_total_UMIs_vector
export cell_type_vector
export gene_axis
export gene_is_covered_vector
export gene_is_excluded_vector
export gene_is_forbidden_vector
export gene_is_lateral_vector
export gene_is_marker_vector
export gene_is_mitochondrial_vector
export gene_is_regulator_vector
export gene_is_ribosomal_vector
export gene_is_skeleton_vector
export gene_marker_rank_vector
export metacell_axis
export metacell_block_vector
export metacell_covered_UMIs_vector
export metacell_gene_covered_fraction_matrix
export metacell_gene_geomean_fraction_matrix
export metacell_gene_linear_fraction_matrix
export metacell_gene_log_covered_fraction_matrix
export metacell_gene_log_geomean_fraction_matrix
export metacell_gene_log_linear_fraction_matrix
export metacell_gene_UMIs_matrix
export metacell_metacell_euclidean_skeleton_distance
export metacell_metacell_max_skeleton_fold_distance
export metacell_n_cells_vector
export metacell_total_UMIs_vector
export metacell_type_vector
export type_axis
export type_color_vector

using DataAxesFormats

## Genes

"""
    gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced genes. By convention we use gene symbols as the namespace of the genes, but this may be different
depending on the data set.
"""
function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "gene" => (expectation, "Sequenced genes.")
end

### Genes Masks

"""
    gene_is_mitochondrial_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of mitochondrial genes. These genes are typically excluded from the analysis as they are (mostly) unrelated to
the biological behaviors of interest, and also have such a large and variable expression level that including them would
skew the denominator when estimating linear gene expression level fractions.

This vector is created in a supervised way based on biological and technical considerations.
"""
function gene_is_mitochondrial_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_mitochondrial") => (expectation, Bool, "A mask of mitochondrial genes.")
end

"""
    gene_is_ribosomal_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of mitochondrial genes. These genes are typically excluded from the analysis, as they aren't always unrelated to
the biological behaviors of interest, they have such a large and variable expression level that including them would
skew the denominator when estimating linear gene expression level fractions.

This vector is created in a supervised way based on biological and technical considerations.
"""
function gene_is_ribosomal_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_mitochondrial") => (expectation, Bool, "A mask of ribosomal genes.")
end

"""
    gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are excluded from consideration. These genes are either unrelated to the biological behaviors of
interest and/or have such a large and variable expression level that including them would skew the denominator when
estimating linear genes expression level fractions.

This vector is created in a supervised way based on biological and technical considerations.
"""
function gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_excluded") => (expectation, Bool, "A mask of genes that are excluded from consideration.")
end

"""
    gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are lateral to the biological behaviors of interest. These genes may satisfy all criteria for being
in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors of
interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological
behaviors of interest. They are therefore masked out during some phases of the analysis. However, unlike excluded genes,
we still track their expression level and take it into account (e.g., when removing outlier cells from metacells).

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
to other considerations (e.g., they may be lateral genes).

This vector is populated by [`identify_marker_genes!`](@ref Metacells.AnalyzeGenes.identify_marker_genes!).
"""
function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_marker") => (expectation, Bool, "A mask of genes that distinguish between cell states.")
end

"""
    gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The relative ranks of the marker genes. The more the gene distinguishes between different cell states, the better
(lower) rank it has. That is, 1 is the "most" marker gene. Non-marker genes are given an extremely high rank (that
maximal the data type allows).

This vector is populated by [`rank_marker_genes!`](@ref Metacells.AnalyzeGenes.rank_marker_genes!).
"""
function gene_marker_rank_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "marker_rank") => (expectation, StorageUnsigned, "The ralative ranks of the marker genes.")
end

"""
    gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are covered by the model(s). When we construct some model(s) for approximating the manifold
describing the cell behaviors, we focus only on these genes, ignoring the rest. This allows the model(s) to be more
sensitive to variations in the covered genes, without the distraction of the noise from the rest. The covered genes are
typically some subset of the marker genes; in particular, they don't include lateral genes. Properties that are computed
using only the covered genes contain `_covered_` in their name, e.g., `covered_fraction` vs. `linear_fraction`.

This vector is populated by [`identify_covered_genes!`](@ref Metacells.AnalyzeGenes.identify_covered_genes!).
"""
function gene_is_covered_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_covered") => (expectation, Bool, "A mask of genes that are covered by the model(s).")
end

"""
    gene_is_regulator_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that regulate the expression level of other genes. These genes are expected to be at the core of the
gene programs that describe the cell behaviors manifold.

This vector is populated by [`fetch_regulators!`](@ref Metacells.AnalyzeGenes.fetch_regulators!).

!!! note

    These are **not** all the genes that bind to DNA (aka "transcription factors").
"""
function gene_is_regulator_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_regulator") =>
        (expectation, Bool, "A mask of genes that regulate the expression level of other genes.")
end

"""
    gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are used to predict the values of the rest of the (covered) genes. We assume that knowing the
expression level of the skeleton genes is sufficient to reasonably estimate the expression level of the rest of the
covered genes. For example, if two metacells have "very close" expression level of all the skeleton genes, we assume
that these metacells are "very similar", without looking at the rest of the covered genes.

This vector is populated by [`identify_skeleton_genes!`](@ref Metacells.AnalyzeGenes.identify_skeleton_genes!).
"""
function gene_is_skeleton_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_skeleton") =>
        (expectation, Bool, "A mask of genes that are used to predict the values of the rest of the (covered) genes.")
end

"""
    gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are forbidden from being used as skeleton genes.

This vector is created in a supervised way based on biological and technical considerations.
"""
function gene_is_forbidden_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_forbidden") =>
        (expectation, Bool, "A mask of genes that are forbidden from being used as skeleton genes.")
end

## Cells

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

### Cells UMIs

"""
    cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of UMIs collected for each gene for each cell. This is the measurement everything else is built on. This is
some sampling of the ground truth of the number of RNA molecules of each gene in each cell. This suffers from randomness
(only a small subset of the molecules are actually counted), and biases (some molecules are favored over others,
depending on the specific technology and protocol used). More fundamentally, the total number of RNA molecules varies
(sometimes wildly) between different cells.

The UMIs count is therefore very sparse and very noisy. The challenge of scRNA-seq analysis is to extract meaningful
patterns (cooperating gene programs) that describe the cell behaviors based on this data.

To analyze this data, we are forced to normalize the UMIs counts to a fraction of the total. This is mostly meaningless
in the single-cell level because we have so few UMIs per each cell (x,000 or x0,000 UMIs out of the x00,000 RNA
molecules) and low-expressing (but still important) genes may have one or at most few UMIs.

Therefore fractions are only computed for collection of cells. Metacells are defined as the smallest such collections
that can be used to give meaningful fraction estimates (therefore maximizing sensitivity for capturing all cell
behaviors). This is in contrast with classical methods that look at very large collections (cell "types").

Using fractions is sensitive to the choice of the denominator; a gene (or gene program) with very high expression would
artificially reduce the fraction of the rest of the genes. We therefore exclude known troublesome genes from the
denominator (e.g., mitochondrial and ribosomal genes). When constructing parametric models, we recompute the fractions
to just the covered genes, to reduce the impact on the denominator of lateral genes (e.g., cell cycle). And when
comparing fractions between different data sets, one must *always* renormalize the fractions to the common genes set.

Ideally, we could somehow create an estimate of the concentration of each RNA molecule as this is in absolute units so
is directly comparable across analysis methods and data sets. However, in addition to being dependent on the unknown
total number of molecules in each cell, this also depends on the unknown and cell volume.

This data is obtained from scRNA-seq experiments.
"""
function cell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "cell", "UMIs") =>
        (expectation, StorageUnsigned, "The number of UMIs collected for each gene for each cell.")
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
    cell_excluded_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the excluded genes in each cell.

This vector is created in a supervised way based on biological and technical considerations. It typically includes
mitochondrial and ribosomal genes, as well as other genes which are highly correlated with other excluded genes.
"""
function cell_excluded_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "excluded_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the excluded genes in each cell.")
end

"""
    cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each cell.

This vector is populated by [`compute_cells_total_UMIs!`](@ref Metacells.AnalyzeCells.compute_cells_total_UMIs!).
"""
function cell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each cell.")
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

### Cells masks

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

## Metacells

"""
    metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of metacells, which are minimal-sized groups of cells for robust fraction estimates. Each metacell is
considered to be a robustly estimated point in the multi-dimensional manifold of cell states. Metacells may be very
similar or very different from each other depending on how well a region of the manifold was sampled; that is, highly
sampled cell states would have a large number of higher quality metacells that are more similar (dare one say, even
"identical") to each other, while sparsely sampled cell states would have a smaller number of lower quality metacells
that are more different from one another. At the threshold of the sensitivity of the algorithm, very rare cell states
would be discarded as outliers, or not even appear in the data in the first place.
"""
function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "metacell" => (expectation, "Minimal-sized groups of cells for robust point estimates.")
end

"""
    cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique metacell each cell belongs to. All the cells in the same metacell are assumed to have "the same" (relevant)
biological state. This is the empty string if the cell does not belong to any metacell (is excluded or outlier).

This vector can be populated by any metacells-like algorithm (the original R metacells (1) algorithm, the Python
metacells (2) algorithm, other similar algorithms).
"""
function cell_metacell_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "metacell") => (expectation, AbstractString, "The unique metacell each cell belongs to.")
end

"""
    metacell_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of cells in each metacell.

This vector is populated by [`compute_metacells_n_cells!`](@ref Metacells.AnalyzeMetacells.compute_metacells_n_cells!).
"""
function metacell_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "n_cells") => (expectation, StorageUnsigned, "The number of cells in each metacell.")
end

"""
    metacell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs of each gene in the cells of each metacell.

This matrix is populated by [`compute_metacells_genes_UMIs!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_UMIs!).
"""
function metacell_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of each gene in the cells of each metacell.")
end

"""
    metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each metacell. This is used as the denominator for
`linear_fraction`.

This vector is populated by [`compute_metacells_total_UMIs!`](@ref Metacells.AnalyzeMetacells.compute_metacells_total_UMIs!).
"""
function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each metacell.")
end

"""
    metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of the covered genes in each metacell. This is used as the denominator for `covered_fraction`.

This vector is populated by [`compute_metacells_covered_UMIs!`](@ref Metacells.AnalyzeMetacells.compute_metacells_covered_UMIs!).
"""
function metacell_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the covered genes in each metacell.")
end

### Metacells Fractions

"""
    metacell_gene_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The geomean fraction of the UMIs of each gene in each metacell. We use geomean in an attempt to combat the
disproportionate effect of a few cells with very high gene expression ("bursty" genes), and then normalizes the
fractions to sum to one. While effective, this has the unfortunate effect of inflating the value of weak genes.

TODO: Get rid of this.

This matrix is populated by [`compute_metacells_genes_geomean_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_geomean_fractions!).
"""
function metacell_gene_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "geomean_fraction") =>
        (expectation, StorageFloat, "The geomean fraction of the UMIs of each gene in each metacell.")
end

"""
    metacell_gene_log_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the geomean fraction of the UMIs of each gene in each metacell. This adds some gene fraction
regularization to deal with zero fractions.

TODO: Get rid of this.

This matrix is populated by [`compute_metacells_genes_log_geomean_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_geomean_fractions!).
"""
function metacell_gene_log_geomean_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_geomean_fraction") =>
        (expectation, StorageFloat, "The log base 2 of the geomean fraction of the UMIs of each gene in each metacell.")
end

"""
    metacell_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs. This is the "best"
estimate assuming unbiased multinomial sampling. However, this is sensitive to a few cells with very high expression
levels ("bursty" genes), as these impact the denominator.

This matrix is populated by [`compute_metacells_genes_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_linear_fractions!).
"""
function metacell_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total
UMIs. This adds some gene fraction regularization to deal with zero fractions. Using the log makes it easier to
visualize, and the difference between log values (the "fold factor", log base 2 of the ratio between the expression
levels) is a good measure of difference between gene expression levels.

This matrix is populated by [`compute_metacells_genes_log_linear_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_linear_fractions!).
"""
function metacell_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs.",
    )
end

"""
    metacell_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs. By considering
only the covered genes this avoid the impact of highly-expressed lateral genes (e.g., cell cycle).

This matrix is populated by [`compute_metacells_genes_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_covered_fractions!).
"""
function metacell_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "covered_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs.",
    )
end

"""
    metacell_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each covered gene in each metacell, out of the total
covered UMIs. This adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_metacells_genes_log_covered_fractions!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_genes_log_covered_fractions!).
"""
function metacell_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "gene", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each covered gene in each metacell, out of the total covered UMIs.",
    )
end

### Metacells Distances

"""
    metacell_metacell_euclidean_skeleton_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The Euclidean distance between the log of the covered fraction of the skeleton genes between the metacells.

This matrix is computed by [`compute_metacells_euclidean_distances!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_euclidean_distances!).

TODO: Get rid of this?
"""
function metacell_metacell_euclidean_skeleton_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "euclidean_skeleton_distance") => (
        expectation,
        StorageFloat,
        "The Euclidean distance between the log of the covered fraction of the skeleton genes between the metacells.",
    )
end

"""
    metacell_metacell_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The maximal significant fold factor between the covered fraction of skeleton genes between the metacells. This uses
heuristics to require the fold factor be based on a sufficient number of UMIs to be robust.

This matrix may be populated by [`compute_metacells_max_skeleton_fold_distances!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_max_skeleton_fold_distances!).

TODO: Get rid of this?
"""
function metacell_metacell_max_skeleton_fold_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "max_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The maximal significant fold factor between the covered fraction of skeleton genes between the metacells.",
    )
end

## Blocks

"""
    block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of blocks, which are distinct groups of metacells with "very close" cell state. The metacells in each
block all have "very close" estimates of the covered fractions of the skeleton genes. Using blocks instead of metacells
allows us to more uniformly sample the overall manifold. Highly sampled cell states would have blocks with larger number
of higher-quality metacells, while sparsely sampled cell states would have blocks with fewer lower-quality metacells.
"""
function block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "block" => (expectation, "Distinct groups of metacells with \"very close\" cell state.")
end

"""
    metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each metacell belongs to.

This vector is populated by [`compute_blocks!`](@ref Metacells.ComputeBlocks.compute_blocks!).
"""
function metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "block") => (expectation, AbstractString, "The unique block each metacell belongs to.")
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
    block_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The number of cells in the metacells in each block.

This vector is populated by [`compute_blocks_n_cells!`](@ref Metacells.AnalyzeBlocks.compute_blocks_n_cells!).
"""
function block_n_cells_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_cells") => (expectation, StorageUnsigned, "The number of cells in the metacells in each block.")
end

### Blocks UMIs

"""
    block_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs of each gene in each block.

This matrix is populated by [`compute_blocks_genes_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_genes_UMIs!).
"""
function block_gene_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of each gene in each block.")
end

"""
    block_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each block. This is used as the denominator for
`linear_fraction`.

This vector is populated by [`compute_blocks_total_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_total_UMIs!).
"""
function block_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each block.")
end

"""
    block_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the covered genes in each block. This is used as the denominator for `covered_fraction`.

This vector is populated by [`compute_blocks_covered_UMIs!`](@ref Metacells.AnalyzeBlocks.compute_blocks_covered_UMIs!).
"""
function block_covered_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "covered_UMIs") =>
        (expectation, StorageUnsigned, "The total number of the all covered genes in each block.")
end

### Blocks Fractions

"""
    block_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.

This matrix is populated by [`compute_blocks_genes_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_linear_fractions!).
"""
function block_gene_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.",
    )
end

"""
    block_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block. This
adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_blocks_genes_log_linear_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_linear_fractions!).
"""
function block_gene_log_linear_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.",
    )
end

"""
    block_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block.

This matrix is populated by [`compute_blocks_genes_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_covered_fractions!).
"""
function block_gene_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "covered_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.",
    )
end

"""
    block_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each covered gene out of the total covered UMIs in each block. This
adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_blocks_genes_log_covered_fractions!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_log_covered_fractions!).
"""
function block_gene_log_covered_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each covered gene out of the total UMIs in each block.",
    )
end

### Blocks Distances

"""
    block_block_mean_euclidean_skeleton_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The mean Euclidean skeleton genes covered fractions distance between the metacells of the blocks.

This matrix may be computed by [`compute_blocks_mean_euclidean_distances!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_mean_euclidean_distances!).

TODO: Get rid of this?
"""
function block_block_mean_euclidean_skeleton_distance(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "mean_euclidean_skeleton_distance") => (
        expectation,
        StorageFloat,
        "The mean Euclidean skeleton genes covered fractions distance between the metacells of the blocks.",
    )
end

"""
    block_block_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The maximal significant skeleton genes covered fractions fold factor between metacells of the blocks.

This matrix may be populated by [`compute_blocks_max_skeleton_fold_distances!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_max_skeleton_fold_distances!).

TODO: Get rid of this?
"""
function block_block_max_skeleton_fold_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "max_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The maximal significant skeleton genes covered fractions fold factor between metacells of the blocks.",
    )
end

### Blocks Neighborhoods

"""
    block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the (column) mask of nearby blocks in its immediate neighborhood. The neighborhood consists of a small
number of very close blocks.

This matrix is populated by [`compute_blocks_is_in_neighborhood!`](@ref
Metacells.ComputeBlocks.compute_blocks_is_in_neighborhood!).
"""
function block_block_is_in_neighborhood_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_neighborhood") =>
        (expectation, Bool, "For each block, the (column) mask of nearby blocks in its immediate neighborhood.")
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

### Blocks Environments

"""
    block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

For each block, the (column) mask of nearby blocks in its expanded environment. The environment consists of a larger
region than the neighborhood, which is still close to the center block.

This matrix is populated by [`compute_blocks_is_in_environment!`](@ref
Metacells.ComputeBlocks.compute_blocks_is_in_environment!).
"""
function block_block_is_in_environment_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_environment") =>
        (expectation, Bool, "For each block, the (column) mask of nearby blocks in its expanded environment.")
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

"""
    block_gene_is_environment_marker_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

A mask of covered genes that distinguish between cell states in each environment.

This matrix is populated by [`compute_blocks_genes_is_environment_markers!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_is_environment_markers!).
"""
function block_gene_is_environment_marker_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "is_environment_marker") =>
        (expectation, Bool, "A mask of covered genes that distinguish between cell states in each environment.")
end

### Blocks Correlations

"""
    block_gene_neighborhood_correlation_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The correlation between cells and metacells gene expression levels in each block's neighborhood. This is zero for excluded genes.

This matrix is populated by [`compute_blocks_genes_neighborhood_correlations!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_genes_neighborhood_correlations!).
"""
function block_gene_neighborhood_correlation_matrix(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "neighborhood_correlation") => (
        expectation,
        StorageFloat,
        "The correlation between cells and metacells gene expression levels in each block's neighborhood.",
    )
end

"""
    block_mean_neighborhood_correlation_vector(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The mean correlation between cells and metacells gene expression levels in each block's neighborhood.

This matrix is populated by [`compute_blocks_mean_neighborhood_correlations!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_mean_neighborhood_correlations!).
"""
function block_mean_neighborhood_correlation_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "mean_neighborhood_correlation") => (
        expectation,
        StorageFloat,
        "The mean correlation between cells and metacells gene expression levels in each block's neighborhood.",
    )
end

"""
    block_mean_neighborhood_covered_correlation_vector(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The mean correlation between cells and metacells covered gene expression levels in each block's neighborhood.

This matrix is populated by [`compute_blocks_mean_covered_neighborhood_correlations!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_mean_covered_neighborhood_correlations!).
"""
function block_mean_covered_neighborhood_correlation_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "mean_covered_neighborhood_correlation") => (
        expectation,
        StorageFloat,
        "The mean correlation between cells and metacells covered gene expression levels in each block's neighborhood.",
    )
end

"""
    block_mean_neighborhood_modules_correlation_vector(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The mean correlation between cells and metacells modules gene expression levels in each block's neighborhood.

This matrix is populated by [`compute_blocks_mean_modules_neighborhood_correlations!`](@ref
Metacells.AnalyzeBlocks.compute_blocks_mean_modules_neighborhood_correlations!).
"""
function block_mean_modules_neighborhood_correlation_vector(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "mean_modules_neighborhood_correlation") => (
        expectation,
        StorageFloat,
        "The mean correlation between cells and metacells modules gene expression levels in each block's neighborhood.",
    )
end

"""
    block_gene_is_correlated_with_metacells(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

Whether each gene is strongly correlated between cells and metacells in each neighborhood.
"""
function block_gene_is_correlated_with_metacells(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "is_correlated_with_metacells") => (
        expectation,
        Bool,
        "Whether each gene is strongly correlated between cells and metacells in each neighborhood.",
    )
end

## Type Annotations

"""
    type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of types, which are distinct named biological cell states. Types are convenient labels manually assigned to
large groups of cells, metacells, blocks, etc. Types are not associated with an exact biological cell state, but rather
with a set of related biological cell states, possibly along a gradient of such states. The resolution of the type
labels therefore depends on the data set and the type of analysis.
"""
function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "type" => (expectation, "Distinct named biological cell states.")
end

"""
    type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A unique color for each type, for visualizations.

This vector is created in a supervised way based on biological considerations and conventions (e.g., red blood cells are
often given some red color). It is also possible to use `Chameleon` to automatically assign colors to cell types based
on some gene expression levels.
"""
function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("type", "color") => (expectation, AbstractString, "A unique color for each type for graphs.")
end

"""
    metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each metacell belongs to. This can be assigned in many ways, from manual annotation to using projection methods
on an annotated atlas. Often there are multiple type annotations due to different methods. By convention these alternate
type vectors are named `type.something` (e.g., `type.projected`, `type.manual`) and the simple `type` is used for "the"
type of each metacell.

This vector is created in a supervised way based on biological considerations. Alternatively, if one has a cell type
vector from somewhere, it can be populated by [`compute_metacells_types_by_cells!`](@ref
Metacells.AnalyzeMetacells.compute_metacells_types_by_cells!) which assigns to each metacell the most frequent type of
the cells it contains.
"""
function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "type") => (expectation, AbstractString, "The type each metacell belongs to.")
end

"""
    cell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each cell belongs to. By convention `type` is used for "the" type of the cell, which is compatible with the
type of the metacell it belongs to, and other `type.something` are used for types computed in other ways, such as
importing data from a type-annotated single-cell atlas, or by running some classifier directly on the cells.

This vector can be populated by [`compute_cells_types_by_metacells!`](@ref
Metacells.AnalyzeMetacells.compute_cells_types_by_metacells!) which assigns to each cell the type of the metacell it
belongs to.
"""
function cell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "type") => (expectation, AbstractString, "The type each cell belongs to.")
end

"""
    block_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each block belongs to. This is typically taken to be the most frequent metacell type in the block.

This vector is populated by [`compute_blocks_types!`](@ref Metacells.AnalyzeBlocks.compute_blocks_types!) based on the
metacell types in each block.
"""
function block_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "type") => (expectation, AbstractString, "The type each block belongs to.")
end

## Gene Modules

"""
    module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of modules, which are groups of genes that together can predict the rest of the covered genes. Each module is
based on a single anchor regulator gene, and may contain additional regulator and regular genes. The set and composition
of modules varies between blocks across the manifold. The names of the modules are the names of their anchors (for
interpretability) followed by a `.MOD` suffix (to clarify this is the name of a module and not a gene).
"""
function module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "module" => (expectation, "Groups of genes that together predict the rest of the covered genes.")
end

"""
    module_anchor_gene_vector(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The anchor gene each module is based on. Even though the "same" module in different blocks is based on the same anchor
gene, the rest of the composition of the module may vary wildly, even between neighboring blocks.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function module_anchor_gene_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("module", "gene.anchor") => (expectation, AbstractString, "The anchor gene each module is based on.")
end

"""
    block_module_is_found_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

A mask of the modules that were found for each block. Due to `Daf` limitations, the modules axis must cover all the
modules of all the blocks. This mask specifies which of the modules were actually found for each of the blocks.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function block_module_is_found_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "is_found") =>
        (expectation, Bool, "A mask of the modules that were found for each block.")
end

"""
    block_gene_module_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The module each gene belongs to in each block. Most genes do not belong to any module and therefore
have an empty string as the value.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function block_gene_module_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "module") =>
        (expectation, AbstractString, "The module each gene belongs to in each block.")
end

"""
    block_module_is_strong_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

A mask of the strong modules that have enough UMIs in enough cells. Since cells have much less UMIs than metacells,
not all modules have a high enough number of UMIs in enough cells to be useful.

This matrix is populated by [`compute_blocks_modules_is_strong!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_is_strong!).
"""
function block_module_is_strong_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "is_strong") =>
        (expectation, Bool, "A mask of the modules that have enough UMIs in enough cells.")
end

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
    block_module_n_skeletons_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The number of skeleton genes in each gene module in each block. This is zero for the extra (not-found) gene modules.

This matrix is populated by [`compute_blocks_modules_n_skeletons!`](@ref
Metacells.AnalyzeModules.compute_blocks_modules_n_skeletons!).
"""
function block_module_n_skeletons_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "n_skeletons") =>
        (expectation, StorageUnsigned, "The number of skeleton genes in each gene module in each block.")
end

"""
    block_metacell_module_environment_total_UMIs_tensor(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total UMIs of each environment module in each environment metacell. This is zero for metacells that are not part of
the block's environment, or modules that aren't found in that environment.

This vector is populated by [`compute_block_metacell_module_environment_total_UMIs!`](@ref
Metacells.AnalyzeModules.compute_block_metacell_module_environment_total_UMIs!).
"""
function block_metacell_module_environment_total_UMIs_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total UMIs of each environment module in each environment metacell.")
end

"""
    block_metacell_module_environment_linear_fraction_tensor(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total UMIs. This is zero
for metacells that are not part of the block's environment, or modules that aren't found in that environment.

This vector is populated by [`compute_block_metacell_module_environment_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_block_metacell_module_environment_linear_fractions!).
"""
function block_metacell_module_environment_linear_fraction_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total UMIs.",
    )
end

"""
    block_metacell_module_environment_log_linear_fraction_tensor(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of
the total UMIs. This adds some gene fraction regularization to deal with zero fractions. This is zero for metacells that
are not part of the block's environment, or modules that aren't found in that environment.

This vector is populated by [`compute_block_metacell_module_environment_log_linear_fractions!`](@ref
Metacells.AnalyzeModules.compute_block_metacell_module_environment_log_linear_fractions!).
"""
function block_metacell_module_environment_log_linear_fraction_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of the total UMIs.",
    )
end

"""
    block_metacell_module_environment_covered_fraction_tensor(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total covered
UMIs. This is zero for metacells that are not part of the block's environment, or modules that aren't found in that
environment.

This vector is populated by [`compute_block_metacell_module_environment_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_block_metacell_module_environment_covered_fractions!).
"""
function block_metacell_module_environment_covered_fraction_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "covered_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the total UMIs of each environment module in each environment metacell, out of the total covered UMIs.",
    )
end

"""
    block_metacell_module_environment_log_covered_fraction_tensor(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of
the total covered UMIs. This adds some gene fraction regularization to deal with zero fractions. This is zero for
metacells that are not part of the block's environment, or modules that aren't found in that environment.

This vector is populated by [`compute_block_metacell_module_environment_log_covered_fractions!`](@ref
Metacells.AnalyzeModules.compute_block_metacell_module_environment_log_covered_fractions!).
"""
function block_metacell_module_environment_log_covered_fraction_tensor(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "log_covered_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the total UMIs of each environment module in each environment metacell, out of the total covered UMIs.",
    )
end

end  # module
