"""
Functions for defining a `Contract` for a metacells `Daf` data set `@computation`. This also serves as a quick
vocabulary of the contents of a metacells `Daf` data set.

In the descriptions below, "fold factor" refers to the log base 2 of the ratio between expression levels. For gene RNA
expression levels, we typically use a regularization factor of 1e-5; fold factors of 1 (2x) is typically considered to
be irrelevant, fold factors of 2 (4x) are considered to be potentially significant, fold factor of 3 (8x) are considered
to point to a definite significant difference, etc. When used as thresholds, fold factors are adjusted using confidence
levels (based on the number of UMIs used for the estimates) and the [`gene_divergence_vector`](@ref).
"""
module Contracts

export block_axis
export block_block_distance_matrix
export cell_axis
export gene_axis
export gene_divergence_vector
export gene_is_global_predictive_factor_vector
export gene_is_lateral_vector
export gene_is_marker_vector
export gene_is_transcription_factor_vector
export gene_metacell_fraction_matrix
export gene_metacell_total_UMIs_matrix
export metacell_axis
export metacell_block_vector
export metacell_total_UMIs_vector
export metacell_type_vector
export type_axis
export type_color_vector

using Daf

"""
    function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced genes. By convention we use TODO as the namespace of the genes, but this may be different
depending on the data set.
"""
function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "gene" => (expectation, "Sequenced genes.")
end

"""
    function cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of sequenced single cells. There's no convention for cell names, as long as they are unique. Typically some
sort of barcode is used, possibly combined with a batch and/or plate and/or experiment identification. In the latter
case it is recommended that `batch` and/or `plate` and/or `experiment` would also be created as explicit axes, to allow
associating metadata with them instead of repeating it for each cell.
"""
function cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "cell" => (expectation, "Sequenced single cells.")
end

"""
    function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of metacells, which are minimal-sized groups of cells for robust point estimates. That is, each metacell is
considered to be a robustly estimated point in the multi-dimensional manifold of cell states. Metacells may be very
similar or very different from each other depending on the data set and the manifold region.
"""
function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}
    return "metacell" => (expectation, "Minimal-sized groups of cells for robust point estimates.")
end

"""
    function block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of blocks, which are distinct groups of metacells with "very close" estimated cell state. That is, for some
chosen set of genes, the metacells in each block all have very close estimates of gene expressions (maximal fold
factor up to some maximal span). This maximal span is small and identical for all the blocks.
"""
function block_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "block" => (expectation, "Distinct groups of metacells with \"very close\" estimated cell state.")
end

"""
    function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}

The axis of types, which are distinct named biological cell states. That is, types are convenient labels manually
assigned to large groups of cells, metacells, blocks, neighborhoods, etc. The resolution of the type labels depends on
the data set and the type of analysis. In particular, types are not typically associated with a specific biological cell
state, but rather with a set of related biological cell states (possibly along a gradient of such states).
"""
function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}  # untested
    return "type" => (expectation, "Distinct named biological cell states.")
end

"""
    function gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

For each gene, we scale fold factors of each gene by multiplying with (1 - divergence) of the gene. In particular this
is used when considering the distance between gene expressions. Therefore genes with a non-zero divergence will require
a higher fold factor
"""
function gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "divergence") =>
        (expectation, StorageFloat, "Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.")
end

"""
    function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that distinguish between cell states. That is, these genes have a significant expression level at some
cell state, as well as a significant range of expression across all cell states, so can be used to distinguish between
cell states. Non-marker genes are by definition not useful for such analysis, but marker genes aren't necessarily useful
due to other considerations.
"""
function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}
    return ("gene", "is_marker") => (expectation, Bool, "A mask of genes that distinguish between cell states.")
end

"""
    function gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that are lateral to the biological behaviors of interest. That is, these genes may satisfy all criteria for
being in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors
of interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological
behaviors of interest. They are therefore masked out during the analysis.
"""
function gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_lateral") =>
        (expectation, Bool, "A mask of genes that are lateral to the biological behaviors of interest.")
end

"""
    function gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

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
    function gene_is_transcription_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of genes that bind to the DNA. That is, these genes are good candidates for controlling the expression of other
genes.
"""
function gene_is_transcription_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_transcription_factor") => (expectation, Bool, "A mask of genes that bind to the DNA.")
end

"""
    function gene_is_global_predictive_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of globally predictive transcription factors. That is, knowing the values of all these genes allows us to predict the
values of the rest of the genes (but not as well as when using the locally predictive transcription factors).
"""
function gene_is_global_predictive_factor_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("gene", "is_global_predictive_factor") =>
        (expectation, Bool, "A mask of globally predictive transcription factors.")
end

"""
    function cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A mask of cells that are excluded from consideration. This can be due to any number of reasons - doublets, too
low a number of UMIs, to high a percentage of excluded gene UMIs, etc.
"""
function cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("cell", "is_excluded") =>
        (expectation, Bool, "A mask of cells that are totally excluded from the analysis.")
end

"""
    function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The total number of UMIs used to estimate the fraction of all the genes in each metacell. This is used to estimate
the robustness of the estimates.
"""
function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of all the genes in each metacell.",
    )
end

"""
    function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The type each metacell belongs to.
"""
function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "type") => (expectation, AbstractString, "The type each metacell belongs to.")
end

"""
    function metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

The unique block each metacell belongs to.
"""
function metacell_block_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("metacell", "block") => (expectation, AbstractString, "The unique block each metacell belongs to.")
end

"""
    function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}

A unique color for each type for graphs.
"""
function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}  # untested
    return ("type", "color") => (expectation, AbstractString, "A unique color for each type for graphs.")
end

"""
    function gene_metacell_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The estimated fraction of the UMIs of each gene in each metacell. We assume that each metacell is a sample of the
manifold, representing a real biological state, regardless to its distance to other metacells (subject to cleaning up
batch effects, purging doublets, and compensating for any other technical artifacts).
"""
function gene_metacell_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}
    return ("gene", "metacell", "fraction") =>
        (expectation, StorageFloat, "The estimated fraction of the UMIs of each gene in each metacell.")
end

"""
    function gene_metacell_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The total number of UMIs used to estimate the fraction of each gene in each metacell. This is used to estimate the
robustness of the estimate. When computing fold factors, we require the total number of UMIs (from both compared
estimates) to be some minimum, and possibly adjust the fold factor according to some confidence level (assuming a
multinomial sampling distribution).
"""
function gene_metacell_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("gene", "metacell", "total_UMIs") => (
        expectation,
        StorageUnsigned,
        "The total number of UMIs used to estimate the fraction of each gene in each metacell.",
    )
end

"""
    function block_block_distance_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}

The distance (fold factor) between the most different metacell genes between the blocks. This is the fold factor between
the most different gene expression in a pair of metacells, one in each block. This considers only the global predictive
genes.
"""
function block_block_distance_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}  # untested
    return ("block", "block", "distance") => (
        expectation,
        StorageFloat,
        "The distance (fold factor) between the most different metacell genes between the blocks.",
    )
end

end  # module
