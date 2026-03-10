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

Naming convention is `vector_of_something_per_axis` for vectors (e.g., `vector_of_metacell_per_cell`) and
`matrix_of_something_per_axis_per_axis` for matrices (e.g., `vector_of_UMIs_per_gene_per_cell`). For square symmetric
matrices, we use `matrix_of_something_between_axis` (e.g., `matrix_of_euclidean_skeleton_fold_distance_between_metacells`).
We also allow for `tensor_of_something_per_axis_per_axis_per_axis` - these are simply a set of matrices, one for each
entry of the first axis, whose names are `entry_something`.
"""
module Contracts

export base_block_axis
export block_axis
export cell_axis
export gene_axis
export matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block
export matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
export matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block
export matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block
export matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block
export matrix_of_euclidean_skeleton_fold_distance_between_metacells
export matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block
export matrix_of_is_found_per_module_per_block
export matrix_of_is_in_neighborhood_per_block_per_block
export matrix_of_is_neighborhood_distinct_per_gene_per_block
export matrix_of_is_neighborhood_marker_per_gene_per_block
export matrix_of_linear_fraction_per_gene_per_block
export matrix_of_linear_fraction_per_gene_per_metacell
export matrix_of_log_linear_fraction_per_gene_per_block
export matrix_of_log_linear_fraction_per_gene_per_metacell
export matrix_of_max_skeleton_fold_distance_between_metacells
export matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks
export matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
export matrix_of_module_per_gene_per_block
export matrix_of_module_status_per_gene_per_block
export matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
export matrix_of_n_genes_per_module_per_block
export matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
export matrix_of_UMIs_per_gene_per_block
export matrix_of_UMIs_per_gene_per_cell
export metacell_axis
export module_axis
export projected_block_axis
export tensor_of_linear_fraction_per_block_per_module_per_metacell
export type_axis
export vector_of_anchor_per_module
export vector_of_base_block_per_metacell
export vector_of_block_closest_by_pertinent_markers_per_cell
export vector_of_block_per_metacell
export vector_of_color_per_type
export vector_of_correlation_between_cells_and_projected_metacells_per_gene
export vector_of_correlation_between_cells_and_punctuated_metacells_per_gene
export vector_of_excluded_UMIs_per_cell
export vector_of_is_correlated_with_skeleton_per_gene
export vector_of_is_excluded_per_cell
export vector_of_is_excluded_per_gene
export vector_of_is_forbidden_per_gene
export vector_of_is_lateral_per_gene
export vector_of_is_marker_per_gene
export vector_of_is_mitochondrial_per_gene
export vector_of_is_regulator_per_gene
export vector_of_is_ribosomal_per_gene
export vector_of_is_skeleton_per_gene
export vector_of_is_transcription_factor_per_gene
export vector_of_marker_rank_per_gene
export vector_of_mean_euclidean_modules_cells_distance_per_metacell
export vector_of_metacell_per_cell
export vector_of_mitochondrial_UMIs_per_cell
export vector_of_n_cells_per_block
export vector_of_n_cells_per_metacell
export vector_of_n_metacells_per_block
export vector_of_n_modules_per_block
export vector_of_n_neighborhood_blocks_per_block
export vector_of_n_neighborhood_cells_per_block
export vector_of_n_neighborhood_metacells_per_block
export vector_of_ribosomal_UMIs_per_cell
export vector_of_std_euclidean_modules_cells_distance_per_metacell
export vector_of_total_neighborhood_UMIs_per_block
export vector_of_total_UMIs_per_block
export vector_of_total_UMIs_per_cell
export vector_of_total_UMIs_per_metacell
export vector_of_type_per_block
export vector_of_type_per_cell
export vector_of_type_per_metacell

using DataAxesFormats

## Genes

"""
    gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced genes. By convention we use gene symbols as the namespace of the genes, but this may be different
depending on the data set.

This axis is typically created when importing data.
"""
function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "gene" => (expectation, "Sequenced genes.")
end

### Genes Masks

"""
    vector_of_is_mitochondrial_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of mitochondrial genes. These genes are typically excluded from the analysis as they are (mostly) unrelated to
the biological behaviors of interest, and also have such a large and variable expression level that including them would
skew the denominator when estimating linear gene expression level fractions.

This vector is created in a supervised way based on biological and technical considerations. TODO: Add this list in Gmara.
"""
function vector_of_is_mitochondrial_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_mitochondrial") => (expectation, Bool, "A mask of mitochondrial genes.")
end

"""
    vector_of_is_ribosomal_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of mitochondrial genes. These genes are typically excluded from the analysis, as they aren't always unrelated to
the biological behaviors of interest, they have such a large and variable expression level that including them would
skew the denominator when estimating linear gene expression level fractions.

This vector is created in a supervised way based on biological and technical considerations. TODO: Add this list in Gmara.
"""
function vector_of_is_ribosomal_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_mitochondrial") => (expectation, Bool, "A mask of ribosomal genes.")
end

"""
    vector_of_is_excluded_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that are excluded from consideration. These genes are either unrelated to the biological behaviors of
interest and/or have such a large and variable expression level that including them would skew the denominator when
estimating linear genes expression level fractions.

This vector is created in a supervised way based on biological and technical considerations.
"""
function vector_of_is_excluded_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_excluded") => (expectation, Bool, "A mask of genes that are excluded from consideration.")
end

"""
    vector_of_is_lateral_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that are lateral to the biological behaviors of interest. These genes may satisfy all criteria for being
in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors of
interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological
behaviors of interest. They are therefore masked out during some phases of the analysis. However, unlike excluded genes,
we still track their expression level and take it into account (e.g., when removing outlier cells from metacells).

This vector is created in a supervised way based on biological considerations. TODO: Add lists in Gmara (cell cycle,
stress, etc.) so this can be initialized by them.
"""
function vector_of_is_lateral_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_lateral") =>
        (expectation, Bool, "A mask of genes that are lateral to the biological behaviors of interest.")
end

"""
    vector_of_is_marker_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that distinguish between cell states. These genes have a significant expression level at some cell
state, as well as a significant range of expression across all cell states, so can be used to distinguish between cell
states. Non-marker genes are by definition not useful for such analysis, but marker genes aren't necessarily useful due
to other considerations (e.g., they may be lateral genes).

This vector is populated by [`compute_vector_of_is_marker_per_gene!`](@ref
Metacells.AnalyzeGenes.compute_vector_of_is_marker_per_gene!).
"""
function vector_of_is_marker_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_marker") => (expectation, Bool, "A mask of genes that distinguish between cell states.")
end

"""
    vector_of_marker_rank_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The relative ranks of the marker genes. The more the gene distinguishes between different cell states, the better
(lower) rank it has. That is, 1 is the "most" marker gene. Non-marker genes are given an extremely high rank (that
maximal the data type allows).

This vector is populated by [`compute_vector_of_marker_rank_per_gene!`](@ref
Metacells.AnalyzeGenes.compute_vector_of_marker_rank_per_gene!).
"""
function vector_of_marker_rank_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "marker_rank") => (expectation, StorageUnsigned, "The ralative ranks of the marker genes.")
end

"""
    vector_of_is_transcription_factor_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that bind to DNA from `Gmara`. Most such genes aren't meaningful regulator of cell behavior.

This vector is populated by [`fetch_gmara_vector_of_is_transcription_factor_per_gene!`](@ref
Metacells.AnalyzeGenes.fetch_gmara_vector_of_is_transcription_factor_per_gene!).

!!! note

    Regardless of the name, these are all the genes that bind to DNA, regardless of their role in actual transcription.
    See [`vector_of_is_regulator_per_gene`](@ref).
"""
function vector_of_is_transcription_factor_per_gene(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_transcription_factor") =>
        (expectation, Bool, "A mask of genes that regulate the expression level of other genes.")
end

"""
    vector_of_is_regulator_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that regulate the expression level of other genes. These genes are expected to be at the core of the
gene programs that describe the cell behaviors manifold.

This vector is populated by [`fetch_gmara_vector_of_is_regulator_per_gene!`](@ref
Metacells.AnalyzeGenes.fetch_gmara_vector_of_is_regulator_per_gene!).

!!! note

    These are **not** all the genes that bind to DNA. See [`vector_of_is_transcription_factor_per_gene`](@ref).
"""
function vector_of_is_regulator_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_regulator") =>
        (expectation, Bool, "A mask of genes that regulate the expression level of other genes.")
end

"""
    vector_of_is_skeleton_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that are used to predict the values of the rest of the genes. We assume that knowing the expression
level of the skeleton genes is sufficient to reasonably estimate the expression level of the rest of the rest of the
(non-excluded) genes. Specifically, if two metacells have "very close" expression level of all the skeleton genes, we
assume that these metacells are "very similar" (will be in the same block), without looking at the rest of the genes.

This vector is populated by [`compute_vector_of_is_skeleton_per_gene!`](@ref
Metacells.AnalyzeGenes.compute_vector_of_is_skeleton_per_gene!).
"""
function vector_of_is_skeleton_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_skeleton") =>
        (expectation, Bool, "A mask of genes that are used to predict the values of the rest of the genes.")
end

"""
    vector_of_is_forbidden_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of genes that are forbidden from being used as skeleton genes.

This vector is created in a supervised way based on biological and technical considerations. If not set, no gene is forbidden.
"""
function vector_of_is_forbidden_per_gene(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_forbidden") =>
        (expectation, Bool, "A mask of genes that are forbidden from being used as skeleton genes.")
end

"""
    vector_of_is_correlated_with_skeleton_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of (marker) genes that have significant correlation with any of the skeleton genes. What is interesting are
marker genes which do not have such correlation; these are potential candidates for looking for missing skeleton genes.

This vector is populated by [`compute_vector_of_is_correlated_with_skeleton_per_gene!`](@ref
Metacells.AnalyzeGenes.compute_vector_of_is_correlated_with_skeleton_per_gene!).
"""
function vector_of_is_correlated_with_skeleton_per_gene(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_correlated_with_skeleton") => (
        expectation,
        Bool,
        "A mask of (marker) genes that have significant correlation with any of the skeleton genes.",
    )
end

## Cells

"""
    cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced single cells. There's no convention for cell names, as long as they are unique. Typically some
sort of barcode is used, possibly combined with a batch and/or plate and/or experiment identification. In the latter
case it is recommended that `batch` and/or `plate` and/or `experiment` would also be created as explicit axes, to allow
associating metadata with them instead of repeating it for each cell.

This axis is typically created when importing data.
"""
function cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "cell" => (expectation, "Sequenced single cells.")
end

### Cells UMIs

"""
    matrix_of_UMIs_per_gene_per_cell(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

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

Therefore fractions are only (well, mostly) computed for collection of cells. Metacells are defined as the smallest such
collections that can be used to give meaningful fraction estimates (therefore maximizing sensitivity for capturing all
cell behaviors). This is in contrast with classical methods that look at very large collections (cell "types").

Using fractions is sensitive to the choice of the denominator; a gene (or gene program) with very high expression would
artificially reduce the fraction of the rest of the genes. We therefore exclude known troublesome genes from the
denominator (e.g., mitochondrial and ribosomal genes). And when comparing fractions between different data sets, one
must *always* renormalize the fractions to the common genes set.

Ideally, we could somehow create an estimate of the concentration of each RNA molecule as this is in absolute units so
is directly comparable across analysis methods and data sets. However, in addition to being dependent on the unknown
total number of molecules in each cell, this also depends on the unknown and cell volume.

This data is obtained from scRNA-seq experiments.
"""
function matrix_of_UMIs_per_gene_per_cell(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "cell", "UMIs") =>
        (expectation, StorageUnsigned, "The number of UMIs collected for each gene for each cell.")
end

"""
    vector_of_mitochondrial_UMIs_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the mitochondrial genes in each cell.

TODO: Implement `compute_vector_of_mitochondrial_UMIs_per_cell` based on `vector_of_is_mitochondrial_per_gene`.
"""
function vector_of_mitochondrial_UMIs_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "mitochondrial_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the mitochondrial genes in each cell.")
end

"""
    vector_of_ribosomal_UMIs_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the ribosomal genes in each cell.

TODO: Implement `compute_vector_of_mitochondrial_UMIs_per_cell` based on `vector_of_is_ribosomal_per_gene`.
"""
function vector_of_ribosomal_UMIs_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "ribosomal_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the ribosomal genes in each cell.")
end

"""
    vector_of_excluded_UMIs_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the excluded genes in each cell.

TODO: Implement `compute_vector_of_excluded_UMIs_per_cell` based on `vector_of_is_excluded_per_gene`.
"""
function vector_of_excluded_UMIs_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "excluded_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the excluded genes in each cell.")
end

"""
    vector_of_total_UMIs_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each cell.

This vector is populated by [`compute_vector_of_total_UMIs_per_cell!`](@ref
Metacells.AnalyzeCells.compute_vector_of_total_UMIs_per_cell!).
"""
function vector_of_total_UMIs_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each cell.")
end

### Cells masks

"""
    vector_of_is_excluded_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A mask of cells that are excluded from consideration. This can be due to any number of reasons - doublets, too low a
number of UMIs, to high a percentage of excluded gene UMIs, etc. Excluded cells are not grouped into metacells etc., so
they are ignored in all properties that sum UMIs per group (e.g., `total_UMIs` per metacell).

This vector is created in a supervised way based on biological and technical considerations.
"""
function vector_of_is_excluded_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
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

This axis is typically created when importing data prepared by the Python metacells package, or by
[`sharpen_metacells!`](@ref Metacells.SharpenMetacells.sharpen_metacells!).
"""
function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "metacell" => (expectation, "Minimal-sized groups of cells for robust point estimates.")
end

"""
    vector_of_metacell_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The unique metacell each cell belongs to. All the cells in the same metacell are assumed to have "the same" (relevant)
biological state. This is the empty string if the cell does not belong to any metacell (is excluded or outlier).

This vector can be populated by any metacells-like algorithm (the original R metacells (1) algorithm, the Python
metacells (2) algorithm, other similar algorithms).
"""
function vector_of_metacell_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "metacell") => (expectation, AbstractString, "The unique metacell each cell belongs to.")
end

"""
    vector_of_n_cells_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The number of cells in each metacell.

This vector is populated by [`compute_vector_of_n_cells_per_metacell!`](@ref
Metacells.AnalyzeMetacells.compute_vector_of_n_cells_per_metacell!).
"""
function vector_of_n_cells_per_metacell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "n_cells") => (expectation, StorageUnsigned, "The number of cells in each metacell.")
end

"""
    matrix_of_UMIs_per_gene_per_metacell(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The total number of UMIs of each gene in the cells of each metacell.

This matrix is populated by [`compute_matrix_of_UMIs_per_gene_per_metacell!`](@ref
Metacells.AnalyzeMetacells.compute_matrix_of_UMIs_per_gene_per_metacell!).
"""
function matrix_of_UMIs_per_gene_per_metacell(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of each gene in the cells of each metacell.")
end

"""
    metacell_total_UMIs_vector(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each metacell. This is used as the denominator for
`linear_fraction`.

This vector is populated by [`compute_vector_of_total_UMIs_per_metacell!`](@ref
Metacells.AnalyzeMetacells.compute_vector_of_total_UMIs_per_metacell!).
"""
function vector_of_total_UMIs_per_metacell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each metacell.")
end

### Metacells Fractions

"""
    matrix_of_linear_fraction_per_gene_per_metacell(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs. This is the "best"
estimate assuming unbiased multinomial sampling. However, this is sensitive to a few cells with very high expression
levels ("bursty" genes), as these impact the denominator.

This matrix is populated by [`compute_matrix_of_linear_fraction_per_gene_per_metacell!`](@ref
Metacells.AnalyzeMetacells.compute_matrix_of_linear_fraction_per_gene_per_metacell!).
"""
function matrix_of_linear_fraction_per_gene_per_metacell(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs.",
    )
end

"""
    matrix_of_log_linear_fraction_per_gene_per_metacell(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs.
This adds some gene fraction regularization to deal with zero fractions. Using the log makes it easier to visualize, and
the difference between log values (the "fold factor", log base 2 of the ratio between the expression levels) is a good
measure of difference between gene expression levels.

This matrix is populated by [`compute_matrix_of_log_linear_fraction_per_gene_per_metacell!`](@ref
Metacells.AnalyzeMetacells.compute_matrix_of_log_linear_fraction_per_gene_per_metacell!).
"""
function matrix_of_log_linear_fraction_per_gene_per_metacell(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each non-excluded gene in each metacell, out of the total UMIs.",
    )
end

### Metacells Distances

"""
    matrix_of_euclidean_skeleton_fold_distance_between_metacells(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The Euclidean distance between the log of the fraction of the skeleton genes between the metacells. This is a symmetric
matrix with zeros in the diagonal.

This matrix is computed by [`compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!`](@ref
Metacells.AnalyzeMetacells.compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!).
"""
function matrix_of_euclidean_skeleton_fold_distance_between_metacells(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "euclidean_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The Euclidean distance between the log of the fraction of the skeleton genes between the metacells.",
    )
end

"""
    matrix_of_max_skeleton_fold_distance_between_metacells(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The maximal significant fold factor between the fraction of skeleton genes between the metacells. This uses heuristics
to require the fold factor be based on a sufficient number of UMIs to be robust. This is a symmetric matrix with
zeros in the diagonal.

This matrix may be populated by [`compute_matrix_of_max_skeleton_fold_distance_between_metacells!`](@ref
Metacells.AnalyzeMetacells.compute_matrix_of_max_skeleton_fold_distance_between_metacells!).
"""
function matrix_of_max_skeleton_fold_distance_between_metacells(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("metacell", "metacell", "max_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The maximal significant fold factor between the fraction of skeleton genes between the metacells.",
    )
end

### Metacells Correlations

"""
    vector_of_correlation_between_cells_and_punctuated_metacells_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The correlation between cells and their metacells (minus the correlated cell) of each gene expression levels. This is
zero for excluded genes. This reflects the quality of the grouping of cells into metacells, and verifying that all
relevant genes are properly described by the metacell model.

This vector may be populated by [`compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!`](@ref
Metacells.AnalyzeMetacells.compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!).
"""
function vector_of_correlation_between_cells_and_punctuated_metacells_per_gene(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("gene", "correlation_between_cells_and_punctuated_metacells") => (
        expectation,
        StorageFloat,
        "The correlation between cells and their metacells (minus the correlated cell) gene expression levels.",
    )
end

## Blocks

"""
    block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of blocks, which are distinct groups of metacells with "very close" cell state. The metacells in each block all
have "very close" estimates of the fractions of the skeleton genes. Using blocks instead of metacells allows us to more
uniformly sample the overall manifold, where each block corresponds to a (coarse) cell state. Highly sampled cell states
would have blocks with larger number of higher-quality metacells, while sparsely sampled cell states would have blocks
with fewer lower-quality metacells.

This axis is typically created by [`compute_metacells_blocks!`](@ref Metacells.ComputeBlocks.compute_metacells_blocks!).
"""
function block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "block" => (expectation, "Distinct groups of metacells with \"very close\" cell state.")
end

"""
    vector_of_block_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The unique block each metacell belongs to.

This vector is populated by [`compute_metacells_blocks!`](@ref
Metacells.ComputeBlocks.compute_metacells_blocks!).
"""
function vector_of_block_per_metacell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "block") => (expectation, AbstractString, "The unique block each metacell belongs to.")
end

"""
    vector_of_base_block_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The unique base block each sharpened metacell belongs to. This refers to the blocks of the base metacells,
not to the new blocks computed for the sharpened metacells (if any).

This vector is populated by [`sharpen_metacells!`](@ref Metacells.SharpenMetacells.sharpen_metacells!).
"""
function vector_of_base_block_per_metacell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("metacell", "base_block") =>
        (expectation, AbstractString, "The unique base block each sharpened metacell belongs to.")
end

"""
    vector_of_n_metacells_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The number of metacells in each block.

This vector is populated by [`compute_vector_of_n_metacells_per_block!`](@ref Metacells.AnalyzeBlocks.compute_vector_of_n_metacells_per_block!).
"""
function vector_of_n_metacells_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_metacells") => (expectation, StorageUnsigned, "The number of metacells in each block.")
end

"""
    vector_of_n_cells_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The number of cells in each block.

This vector is populated by [`compute_vector_of_n_cells_per_block!`](@ref Metacells.AnalyzeBlocks.compute_vector_of_n_cells_per_block!).
"""
function vector_of_n_cells_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_cells") => (expectation, StorageUnsigned, "The number of cells in each block.")
end

### Blocks UMIs

"""
    matrix_of_UMIs_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The total number of UMIs of each gene in each block.

This matrix is populated by [`compute_matrix_of_UMIs_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_UMIs_per_gene_per_block!).
"""
function matrix_of_UMIs_per_gene_per_block(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of each gene in each block.")
end

"""
    vector_of_total_UMIs_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of UMIs of all the non-excluded genes in each block. This is used as the denominator for
`linear_fraction`.

This vector is populated by [`compute_vector_of_total_UMIs_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_total_UMIs_per_block!).
"""
function vector_of_total_UMIs_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "total_UMIs") =>
        (expectation, StorageUnsigned, "The total number of UMIs of all the non-excluded genes in each block.")
end

### Blocks Fractions

"""
    matrix_of_linear_fraction_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.

This matrix is populated by [`compute_matrix_of_linear_fraction_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_linear_fraction_per_gene_per_block!).
"""
function matrix_of_linear_fraction_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.",
    )
end

"""
    matrix_of_log_linear_fraction_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The log base 2 of the linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block. This
adds some gene fraction regularization to deal with zero fractions.

This matrix is populated by [`compute_matrix_of_log_linear_fraction_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_log_linear_fraction_per_gene_per_block!).
"""
function matrix_of_log_linear_fraction_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "log_linear_fraction") => (
        expectation,
        StorageFloat,
        "The log base 2 of the linear fraction of the UMIs of each non-excluded gene out of the total UMIs in each block.",
    )
end

### Blocks Distances

"""
    matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The mean Euclidean skeleton genes fractions distance between the metacells of the blocks. This is a symmetric matrix,
and the diagonal is 0.

This matrix may be computed by [`compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!).
"""
function matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "mean_euclidean_skeleton_fold_distance") => (
        expectation,
        StorageFloat,
        "The mean Euclidean skeleton genes fractions distance between the metacells of the blocks.",
    )
end

### Blocks Confusion

"""
    vector_of_block_closest_by_pertinent_markers_per_cell(
        expectation::ContractExpectation,
    )::Pair{VectorKey, DataSpecification}

The block each cell is closest to considering only the non-lateral marker genes. Distance is in fold factors (difference
in log base 2 of the gene expression level). This adds some gene fraction regularization to deal with zero fractions.

This matrix may be computed by [`compute_vector_of_block_closest_by_pertinent_markers_per_cell!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_block_closest_by_pertinent_markers_per_cell!).
"""
function vector_of_block_closest_by_pertinent_markers_per_cell(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("cell", "block.closest_by_pertinent_markers") => (
        expectation,
        AbstractString,
        "The block each cell is closest to considering only the non-lateral marker genes.",
    )
end

"""
    matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The number of cells which are included in each (column) block and are closest to each (row) block.

This matrix may be computed by [`compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!).
"""
function matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "confusion_by_closest_by_pertinent_markers") => (
        expectation,
        StorageUnsigned,
        "The number of cells which are included in each column (block) and are closest to each (row) block.",
    )
end

### Blocks Neighborhoods

"""
    matrix_of_is_in_neighborhood_per_block_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

For each (column) block, whether each (row) block is in its neighborhood. This is **not** a symmetric matrix.
Neighborhoods consist of a small number of "close" blocks. Each such small (neighborhood) region of the manifold allows
us to investigate the local gradients of the coarse (block) cell states.

This matrix is populated by [`compute_matrix_of_is_in_neighborhood_per_block_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_is_in_neighborhood_per_block_per_block!).
"""
function matrix_of_is_in_neighborhood_per_block_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "block", "is_in_neighborhood") =>
        (expectation, Bool, "For each (column) block, whether each (row) block is in its immediate neighborhood.")
end

"""
    vector_of_n_neighborhood_blocks_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of (rows) blocks in each (column) block's neighborhood.

This vector is populated by [`compute_vector_of_n_neighborhood_blocks_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_blocks_per_block!).
"""
function vector_of_n_neighborhood_blocks_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_blocks") =>
        (expectation, StorageUnsigned, "The total number of (rows) blocks in each (column) block's neighborhood.")
end

"""
    vector_of_n_neighborhood_metacells_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of metacells in each block's neighborhood.

This vector is populated by [`compute_vector_of_n_neighborhood_metacells_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_metacells_per_block!).
"""
function vector_of_n_neighborhood_metacells_per_block(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_metacells") =>
        (expectation, StorageUnsigned, "The total number of metacells in the blocks in each block's neighborhood.")
end

"""
    vector_of_n_neighborhood_cells_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of cells in the metacells in the each block's neighborhood.

This vector is populated by [`compute_vector_of_n_neighborhood_cells_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_cells_per_block!).
"""
function vector_of_n_neighborhood_cells_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("block", "n_neighborhood_cells") =>
        (expectation, StorageUnsigned, "The total number of cells in the metacells in each block's neighborhood.")
end

"""
    vector_of_total_neighborhood_UMIs_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The total number of non-excluded gene UMIs in the cells in the metacells in each block's neighborhood.

This vector is populated by [`compute_vector_of_total_neighborhood_UMIs_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_total_neighborhood_UMIs_per_block!).
"""
function vector_of_total_neighborhood_UMIs_per_block(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("block", "total_neighborhood_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of non-excluded gene UMIs in the cells in the metacells in each block's neighborhood.",
    )
end

"""
    matrix_of_is_neighborhood_marker_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

A mask of genes that distinguish between cell states in each block's neighborhood. Such genes do not necessarily distinguish the
neighborhood from the rest of the population - see [`matrix_of_is_neighborhood_distinct_per_gene_per_block`](@ref).

This matrix is populated by [`compute_matrix_of_is_neighborhood_marker_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_is_neighborhood_marker_per_gene_per_block!).
"""
function matrix_of_is_neighborhood_marker_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "is_neighborhood_marker") =>
        (expectation, Bool, "A mask of genes that distinguish between cell states in each block's neighborhood.")
end

"""
    matrix_of_is_neighborhood_distinct_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

A mask of genes that distinguish between cell states of each block's neighborhood and the rest of the manifold.

This matrix is populated by [`compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!).
"""
function matrix_of_is_neighborhood_distinct_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "is_neighborhood_distinct") => (
        expectation,
        Bool,
        "A mask of genes that distinguish between cell states of each block's neighborhood and the rest of the manifold.",
    )
end

### Neighborhoods Correlations

"""
    matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

Whether each gene has strong correlation with the skeleton genes in each block's neighborhood. If relevant genes are not
correlated with any skeleton gene, the list of skeleton (or, more likely, regulator) genes may be incomplete and need to
be fixed.

This matrix is populated by [`compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!).
"""
function matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "is_correlated_with_skeleton_in_neighborhood") => (
        expectation,
        Bool,
        "Whether each gene has strong correlation with the skeleton genes in each block's neighborhood.",
    )
end

"""
    matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The correlation between cells and their metacells (minus the correlated cell) of each gene's expression levels in each
block's neighborhood. This is zero for excluded genes. This reflects the overall quality of the neighborhood model.

This matrix is populated by
[`compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!).
"""
function matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "block", "correlation_between_neighborhood_cells_and_punctuated_metacells") => (
        expectation,
        StorageFloat,
        "The correlation between cells and their metacells (minus the correlated cell) gene expression levels in each block's neighborhood.",
    )
end

## Type Annotations

"""
    type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of types, which are distinct named biological cell states. Types are convenient labels manually assigned to
large groups of cells, metacells, blocks, etc. Types are not associated with an exact biological cell state, but rather
with a set of related biological cell states, possibly along a gradient of such states. The resolution of the type
labels therefore depends on the data set and the type of analysis.

This axis is typically created manually, or when importing data.

!!! note

    There are often multiple simultaneous type assignments for the same entities during the analysis process. For
    example, we may have per-cell types imported from single-cell analysis of a data set and a different set of types
    from metacell analysis of the same data. By convention, "the" type property is called `type` and alternative
    properties are called `type.something` (e.g., `type.imported`, `type.by_metacells`). All the functions that deal
    with type look at and/or set "the" `type` property; specify an `adapter` to redirect them to an alternative type
    property.
"""
function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "type" => (expectation, "Distinct named biological cell states.")
end

"""
    vector_of_color_per_type(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

A unique color for each type for visualizations.

This vector is created in a supervised way based on biological considerations and conventions (e.g., red blood cells are
often given some red color). It is also possible to use `Chameleon` to automatically assign colors to cell types based
on some gene expression levels.
"""
function vector_of_color_per_type(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("type", "color") => (expectation, AbstractString, "A unique color for each type for visualizations.")
end

"""
    vector_of_type_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The type each cell belongs to.

When importing an analyzed single-cell data set, it typically comes with per-cell type annotation. When performing
metacell based analysis for a data set, then type annotation is done at the metacell level. In this case, the per-cell
type vector can be populated by [`compute_vector_of_type_per_cell_by_metacells!`](@ref
Metacells.AnalyzeMetacells.compute_vector_of_type_per_cell_by_metacells!). If the type annotation is done at the block
level, the per-cell type vector can be populated by [`compute_vector_of_type_per_cell_by_blocks!`](@ref
"""
function vector_of_type_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "type") => (expectation, AbstractString, "The type each cell belongs to.")
end

"""
    vector_of_type_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The type each metacell belongs to.

Type annotation for metacells is typically a supervised process based on biological and technical considerations. When
importing an analyzed single-cell data set, it typically comes with per-cell type annotation. This allows populating (or
merely initializing) the per-metacell type vector by [`compute_vector_of_type_per_metacell_by_cells!`](@ref
Metacells.AnalyzeMetacells.compute_vector_of_type_per_metacell_by_cells!). If manually annotating the types of blocks
instead, it is also possible to populate the per-metacell type vector based on the blocks by
[`compute_vector_of_type_per_metacell_by_blocks!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_type_per_metacell_by_blocks!).
"""
function vector_of_type_per_metacell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "type") => (expectation, AbstractString, "The type each metacell belongs to.")
end

"""
    vector_of_type_per_block(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The type each block belongs to.

Type annotation for blocks can be done in a supervised process (instead of annotating the individual metacells).
Alternatively, is can be populated by [`compute_vector_of_type_per_block_by_metacells!`](@ref
Metacells.AnalyzeBlocks.compute_vector_of_type_per_block_by_metacells!) based on the metacell types in each block, or by
[`compute_vector_of_type_per_block_by_cells!`](@ref Metacells.AnalyzeBlocks.compute_vector_of_type_per_block_by_cells!)
based on the cell types in each block,
"""
function vector_of_type_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "type") => (expectation, AbstractString, "The type each block belongs to.")
end

## Gene Modules

"""
    module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

Groups of non-lateral genes that together predict the rest of the genes at a specific region of the manifold. Each
module is based on a single anchor skeleton gene, and may contain additional skeleton, regulator and/or marker genes.
The set and composition of modules varies between blocks across the manifold. For convenient, the names of the modules
are the names of their anchors followed by a `.MOD` suffix (to clarify this is the name of a module and not of a gene).

This axis is typically created by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function module_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "module" => (
        expectation,
        "Groups of non-lateral genes that together predict the rest of the genes at a specific region of the manifold.",
    )
end

"""
    vector_of_anchor_per_module(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The skeleton anchor gene each module is based on. Even though the "same" module in different blocks is based on the same
anchor gene, the rest of the composition of the module may vary wildly, even between neighboring blocks.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function vector_of_anchor_per_module(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("module", "gene.anchor") => (expectation, AbstractString, "The anchor gene each module is based on.")
end

"""
    matrix_of_is_found_per_module_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

A mask of the modules that were found for each block. Due to `Daf` limitations, the modules axis must cover all the
modules of all the blocks. This mask specifies which of the modules were actually found for each of the blocks.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function matrix_of_is_found_per_module_per_block(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "is_found") =>
        (expectation, Bool, "A mask of the modules that were found for each block.")
end

"""
    matrix_of_module_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The module each gene belongs to in each block. Most genes do not belong to any module and therefore
have an empty string as the value.

This matrix is populated by [`compute_blocks_modules!`](@ref Metacells.ComputeModules.compute_blocks_modules!).
"""
function matrix_of_module_per_gene_per_block(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "module") =>
        (expectation, AbstractString, "The module each gene belongs to in each block.")
end

"""
    matrix_of_module_status_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

What happened to each gene when computing modules in each block. This is a textual description that is helpful
in debugging.

This matrix is (optionally) populated by [`compute_blocks_modules!`](@ref
Metacells.ComputeModules.compute_blocks_modules!).
"""
function matrix_of_module_status_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "gene", "module_status") =>
        (expectation, AbstractString, "What happened to each gene when computing modules in each block.")
end

"""
    vector_of_n_modules_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The number of found modules in each block.

This vector is populated by [`compute_vector_of_n_modules_per_block!`](@ref
Metacells.AnalyzeModules.compute_vector_of_n_modules_per_block!).
"""
function vector_of_n_modules_per_block(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("block", "n_modules") => (expectation, StorageUnsigned, "The number of found modules in each block.")
end

"""
    matrix_of_n_genes_per_module_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The number of genes in each module in each block. This is zero for the gene modules which were not found in the block.

This matrix is populated by [`compute_matrix_of_n_genes_per_module_per_block!`](@ref
Metacells.AnalyzeModules.compute_matrix_of_n_genes_per_module_per_block!).
"""
function matrix_of_n_genes_per_module_per_block(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "module", "n_genes") =>
        (expectation, StorageUnsigned, "The number of genes in each module in each block.")
end

"""
    tensor_of_linear_fraction_per_block_per_module_per_metacell(
        expectation::ContractExpectation
    )::Pair{TensorKey, DataSpecification}

The linear fraction of the UMIs of each block's found module genes out of the total UMIs in each metacell. We actually
only compute this for metacells in the neighborhood of the block, as each block's modules are not expected to apply
outside this neighborhood. The rest of the entries are zero. Still, this is a tensor with a lot of data, and even with
sparse representation, it has non-trivial storage overhead. It is needed for projecting query cells onto an atlas.

`compute_linear_fraction_of_modules_in_metacells!`

This matrix is populated by [`compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!`](@ref
Metacells.AnalyzeModules.compute_tensor_of_linear_fraction_per_block_per_module_per_metacell!).
"""
function tensor_of_linear_fraction_per_block_per_module_per_metacell(
    expectation::ContractExpectation,
)::Pair{TensorKey, DataSpecification}
    return ("block", "metacell", "module", "linear_fraction") => (
        expectation,
        StorageFloat,
        "The linear fraction of the UMIs of each block's found module genes out of the total UMIs in each metacell.",
    )
end

"""
    matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The mean of the linear fraction of each module in the cells of the neighborhood of each block.

This matrix is populated by [`compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!`](@ref
Metacells.AnalyzeModules.compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!).
"""
function matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "module", "mean_linear_fraction_in_neighborhood_cells") => (
        expectation,
        StorageFloat,
        "The mean of the linear fraction of each module in the cells of the neighborhood of each block.",
    )
end

"""
    matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The standard deviation of the linear fraction of each module in the cells of the neighborhood of each block.

This matrix is populated by [`compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!`](@ref
Metacells.AnalyzeModules.compute_stats_of_linear_fraction_in_neighborhood_cells_per_module_per_block!).
"""
function matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("block", "module", "std_linear_fraction_in_neighborhood_cells") => (
        expectation,
        StorageFloat,
        "The standard deviation of the linear fraction of each module in the cells of the neighborhood of each block.",
    )
end

"""
    vector_of_mean_euclidean_modules_cells_distance_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

In each metacell, the mean of the euclidean distances of between its cells and the metacell's linear fraction of the
modules of the metacell's block's neighborhood.

This matrix is populated by [`compute_stats_of_euclidean_modules_cells_distance_per_metacell!`](@ref
Metacells.AnalyzeModules.compute_stats_of_euclidean_modules_cells_distance_per_metacell!).
"""
function vector_of_mean_euclidean_modules_cells_distance_per_metacell(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "mean_euclidean_modules_cells_distance") => (
        expectation,
        StorageFloat,
        "In each metacell, the mean of the euclidean distances of between its cells and the metacell's linear fraction of the modules of the metacell's block's neighborhood.",
    )
end

"""
    vector_of_std_euclidean_modules_cells_distance_per_metacell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

In each metacell, the standard deviation of the euclidean distance of between its cells and the metacell's linear
fraction of the modules of the metacell's block's neighborhood.

This matrix is populated by [`compute_stats_of_euclidean_modules_cells_distance_per_metacell!`](@ref
Metacells.AnalyzeModules.compute_stats_of_euclidean_modules_cells_distance_per_metacell!).
"""
function vector_of_std_euclidean_modules_cells_distance_per_metacell(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("metacell", "std_euclidean_modules_cells_distance") => (
        expectation,
        StorageFloat,
        "In each metacell, the standard deviation of the euclidean distances of between its cells and the metacell's linear fraction of the modules of the metacell's block's neighborhood.",
    )
end

## Projection

"""
    vector_of_projected_metacell_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The projected atlas metacell for each cell. This is our "best match" of the cell. It does **not** mean that it is
necessarily a *good* match. If you project blood cells on an atlas containing only neuron metacells, there will still be
a "best match" ("least bad match" would be a better term).

This vector is populated by [`compute_cells_projection`](@ref Metacells.ProjectCells.compute_cells_projection!).
"""
function vector_of_projected_metacell_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "projected_metacell") =>
        (expectation, AbstractString, "The projected atlas cell for each metacell.")
end

"""
    vector_of_projected_block_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The atlas block of the projected metacell for each query cell.

This vector is populated by [`compute_cells_projection`](@ref Metacells.ProjectCells.compute_cells_projection!).
"""
function vector_of_projected_block_per_cell(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("cell", "projected_block") =>
        (expectation, AbstractString, "The atlas block of the projected metacell for each query cell.")
end

"""
    vector_of_projected_modules_z_score_per_cell(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The z score of the distance between each query cell and its projected atlas metacell modules. If this is low, then it is
likely that the cell is actually a good match for the metacell.

This vector is populated by [`compute_cells_projection`](@ref Metacells.ProjectCells.compute_cells_projection!).
"""
function vector_of_projected_modules_z_score_per_cell(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("cell", "projected_metacell_modules_z_score") => (
        expectation,
        StorageFloat,
        "The z score of the distance between each query cell and its projected atlas metacell modules.",
    )
end

"""
    vector_of_correlation_between_cells_and_projected_metacells_per_gene(
        expectation::ContractExpectation
    )::Pair{VectorKey, DataSpecification}

The correlation between query cells and their projected atlas metacells of each gene expression levels. This is zero for
excluded genes. This reflects the overall quality of the projection of the query cells to the atlas, and verifying that
all relevant genes are properly described by the atlas metacell model.

This vector may be populated by [`compute_vector_of_correlation_between_cells_and_projected_metacells!`](@ref
Metacells.ProjectCells.compute_vector_of_correlation_between_cells_and_projected_metacells!).
"""
function vector_of_correlation_between_cells_and_projected_metacells_per_gene(
    expectation::ContractExpectation,
)::Pair{VectorKey, DataSpecification}
    return ("gene", "correlation_between_cells_and_projected_metacells") => (
        expectation,
        StorageFloat,
        "The correlation between query cells and their projected atlas metacells gene expression levels.",
    )
end

"""
    matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The correlation between query cells and their projected atlas metacells of each gene's expression levels in each atlas
block's neighborhood. This is zero for excluded genes. This reflects the quality of the projection of the query cells to
the atlas in each region of the manifold.

This matrix is populated by
[`compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!`](@ref
Metacells.ProjectCells.compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!).
"""
function matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "projected_block", "correlation_between_neighborhood_cells_and_projected_metacells") => (
        expectation,
        StorageFloat,
        "The correlation between query cells and their projected atlas metacells of each gene's expression levels in each block's neighborhood.",
    )
end

"""
    projected_block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

A copy of the atlas [`block_axis`](@ref), copied into the query. This is needed to store per-atlas-block quality control
data.

This axis is typically created by
[`compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!`](@ref
Metacells.ProjectCells.compute_matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block!).
"""
function projected_block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "projected_block" => (expectation, "A copy of the atlas [`block_axis`](@ref), copied into the query.")
end

## Quality Control

"""
    matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(
        expectation::ContractExpectation
    )::Pair{MatrixKey, DataSpecification}

The most correlated neighborhood marker gene of each (base) neighborhood marker gene per block. This is
intended for use for computing [`matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block`](@ref).

This is sparse, we actually only compute it for the non-lateral neighborhood markers.

This matrix may be populated by [`compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!).
"""
function matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "block", "most_correlated_gene_in_neighborhood") => (
        expectation,
        AbstractString,
        "The most correlated neighborhood marker gene of each (base) neighborhood marker gene per block.",
    )
end

"""
    base_block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

A copy of the base [`block_axis`](@ref), copied into the alternative. This is needed to store per-base-block quality
control data. Such data is more directly comparable between the alternative metacells and the base ones.

This axis is typically created by
[`compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!).
"""
function base_block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "base_block" => (expectation, "A copy of the base [`block_axis`](@ref), copied into the alternative.")
end

"""
    matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block(
        expectation::ContractExpectation,
    )::Pair{MatrixKey, DataSpecification}

The correlation of each gene with the most correlated other gene between the cells of a base neighborhood and their
metacells of each base neighborhood.

We use this as a poor man's cross validation; we use the set of cells in the base neighborhoods. We measure the
correlation between the expression of these genes in the base and alternative metacells. The absolute correlation values
are of less interest than the difference between different genes, and the difference between the alternative metacells.
Ideally "relevant" genes would have high correlations while "irrelevant" genes would not, and "improved" metacells would
have higher correlation than the base metacells.

These matrices can be populated by
[`compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!).
"""
function matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "base_block", "correlation_with_most_between_base_neighborhood_cells_and_metacells") => (
        expectation,
        StorageFloat,
        "The correlation of each gene with the most correlated other gene between the cells of a base neighborhood and their metacells of each base neighborhood.",
    )
end

"""
    matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(

The correlation between cells and their metacells (minus the correlated cell) of each gene's expression levels in each
base block's neighborhood. This is zero for excluded genes. This allows comparing the correlation between alternative
and base metacells in different locations of the manifold.

This matrix is populated by
[`compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!`](@ref
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!).
"""
function matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block(
    expectation::ContractExpectation,
)::Pair{MatrixKey, DataSpecification}
    return ("gene", "base_block", "correlation_between_base_neighborhood_cells_and_punctuated_metacells") => (
        expectation,
        StorageFloat,
        "The correlation between cells and their metacells (minus the correlated cell) of each gene's expression levels in each base block's neighborhood.",
    )
end

end  # module
