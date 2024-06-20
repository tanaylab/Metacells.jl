var documenterSearchIndex = {"docs":
[{"location":"identify_genes.html#Identify-Genes","page":"Identify Genes","title":"Identify Genes","text":"","category":"section"},{"location":"identify_genes.html","page":"Identify Genes","title":"Identify Genes","text":"Metacells.IdentifyGenes\nMetacells.IdentifyGenes.compute_genes_divergence!\nMetacells.IdentifyGenes.identify_marker_genes!\nMetacells.IdentifyGenes.identify_correlated_genes!","category":"page"},{"location":"identify_genes.html#Metacells.IdentifyGenes","page":"Identify Genes","title":"Metacells.IdentifyGenes","text":"Identify special genes.\n\n\n\n\n\n","category":"module"},{"location":"identify_genes.html#Metacells.IdentifyGenes.compute_genes_divergence!","page":"Identify Genes","title":"Metacells.IdentifyGenes.compute_genes_divergence!","text":"function compute_genes_divergence!(\n    daf::DafWriter;\n    gene_fraction_regularization::AbstractFloat = 1.0e-5,\n    min_divergent_gene_range_fold::AbstractFloat = 6.0,\n    overwrite::Bool = false,\n)::Nothing\n\nCompute a divergence factor for all genes. This factor is typically 0.0 for most genes (and is always below 1.0). Genes that have a wide range of expression are given a higher divergence factor. When looking for significant differences in gene expressions, the divergence factor is used to scale down the difference, to compensate for the gene's inherent wider expression range. Specifically, we take the raw fold factor and multiply it by (1.0 - divergence).\n\nThis works as follows:\n\nCompute the minimal and maximal expression level of each gene.\nCompute the fold factor (log2 of maximal over minimal value, using the gene_fraction_regularization.\nThe raw factor of the gene is the ratio between the fold factor and the min_divergent_gene_range_fold. If this is higher than 1, we set the divergence factor accordingly.\n\nnote: Note\nIdeally, all code that uses the manual is_noisy mask should be modified to use the divergence factor instead.\n\nInputs\n\nAxes\n\ngene (required): Sequenced genes.\n\nmetacell (required): Minimal-sized groups of cells for robust point estimates.\n\nMatrices\n\ngene, metacell @ fraction::AbstractFloat (required): The estimated fraction of the UMIs of each gene in each metacell.\n\nOutputs\n\nVectors\n\ngene @ divergence::AbstractFloat (guaranteed): Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.\n\n\n\n\n\n","category":"function"},{"location":"identify_genes.html#Metacells.IdentifyGenes.identify_marker_genes!","page":"Identify Genes","title":"Metacells.IdentifyGenes.identify_marker_genes!","text":"function identify_marker_genes!(\n    daf::DafWriter;\n    gene_fraction_regularization::AbstractFloat = 1.0e-5,\n    min_marker_gene_range_fold::AbstractFloat = 2.0,\n    min_max_marker_gene_fraction::AbstractFloat = 0.0001,\n    overwrite::Bool = false,\n)::Nothing\n\nIdentify the genes that distinguish at least one metacell from the rest. Such genes are called \"marker\" genes as they (potentially) mark specific cell states. If overwrite, will overwrite an existing is_marker mask.\n\nCompute the minimal and maximal expression level of each gene.\nCompute the fold factor (log2 of maximal over minimal value, using the gene_fraction_regularization.\nReduce this fold factor using the divergence factor of the genes.\nIdentify as markers genes whose adjusted fold factor is at least min_marker_gene_range_fold, and whose maximal expression is at least min_max_marker_gene_fraction.\n\nInputs\n\nAxes\n\ngene (required): Sequenced genes.\n\nmetacell (required): Minimal-sized groups of cells for robust point estimates.\n\nVectors\n\ngene @ divergence::AbstractFloat (required): Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.\n\nMatrices\n\ngene, metacell @ fraction::AbstractFloat (required): The estimated fraction of the UMIs of each gene in each metacell.\n\nOutputs\n\nVectors\n\ngene @ is_marker::Bool (guaranteed): A mask of genes that distinguish between cell states.\n\n\n\n\n\n","category":"function"},{"location":"identify_genes.html#Metacells.IdentifyGenes.identify_correlated_genes!","page":"Identify Genes","title":"Metacells.IdentifyGenes.identify_correlated_genes!","text":"function identify_correlated_genes!(\n    daf::DafWriter;\n    gene_fraction_regularization::AbstractFloat = 1.0e-5,\n    correlation_confidence::AbstractFloat = 0.99,\n    overwrite::Bool = false,\n)::Nothing\n\nIdentify genes that are correlated with other gene(s). Such genes are good candidates for looking for groups of genes that act together. If overwrite, will overwrite an existing is_correlated mask.\n\nCompute the log base 2 of the genes expression in each metacell (using the gene_fraction_regularization).\nCorrelate this between all the pairs of genes.\nFor each gene, shuffle its values along all metacells, and again correlate this between all the pairs of genes.\nFind the maximal absolute correlation for each gene in both cases (that is, strong anti-correlation also counts).\nFind the correlation_confidence quantile correlation of the shuffled data.\nIdentify the genes that have at least that level of correlations in the unshuffled data.\n\nInputs\n\nAxes\n\ngene (required): Sequenced genes.\n\nmetacell (required): Minimal-sized groups of cells for robust point estimates.\n\nMatrices\n\ngene, metacell @ fraction::AbstractFloat (required): The estimated fraction of the UMIs of each gene in each metacell.\n\nOutputs\n\nVectors\n\ngene @ is_correlated::Bool (guaranteed): A mask of genes that are correlated with other gene(s).\n\n\n\n\n\n","category":"function"},{"location":"identify_genes.html#Index","page":"Identify Genes","title":"Index","text":"","category":"section"},{"location":"identify_genes.html","page":"Identify Genes","title":"Identify Genes","text":"Pages = [\"identify_genes.md\"]","category":"page"},{"location":"boxes.html#Boxes","page":"Boxes","title":"Boxes","text":"","category":"section"},{"location":"boxes.html","page":"Boxes","title":"Boxes","text":"Metacells.Boxes\nMetacells.Boxes.compute_boxes!","category":"page"},{"location":"boxes.html#Metacells.Boxes","page":"Boxes","title":"Metacells.Boxes","text":"Given a set of raw metacells, partition them into boxes such that all metacells in the same box are within some (fold factor) radius of each other. The centroids of these boxes can serve as a representation of the cell state manifold which is less sensitive to oversampling of common cell states. Group these boxes in overlapping neighborhoods of \"similar\" boxes for further analysis.\n\n\n\n\n\n","category":"module"},{"location":"boxes.html#Metacells.Boxes.compute_boxes!","page":"Boxes","title":"Metacells.Boxes.compute_boxes!","text":"function compute_boxes!(\n    daf::DafWriter;\n    min_significant_gene_UMIs::Integer = 40,\n    gene_fraction_regularization::AbstractFloat = 1.0e-5,\n    fold_confidence::AbstractFloat = 0.9,\n    max_box_span::AbstractFloat = 2.0,\n    max_neighborhood_span::AbstractFloat = 0.01,\n    correlation_confidence::AbstractFloat = 0.9,\n    max_deviant_genes_fraction::AbstractFloat = 0.01,\n    overwrite::Bool = false,\n)::Nothing\n\nPartition raw metacells into distinct boxes, and boxes into overlapping neighborhoods.\n\nInitial boxes and neighborhoods are computed in a first round, and then refined in a series of followup rounds.\nIn each round, we compute a distance between each two metacells. This is based on the fold factor between the expression level of each (relevant) gene in the metacells. The fold factor is the absolute value of the difference in the log (base 2) of the fraction of the gene in the metacells. This log is computed with the gene_fraction_regularization (by default, 1e-5). Since the fraction of the gene is a random variable, we decrease the high fraction and increase the low fraction by a factor based on the fold_confidence of the test (by default, 0.9), assuming a multinomial distribution. In addition, if the sum of the total UMIs of the gene in both metacells is less than min_significant_gene_UMIs (by default, 40), we ignore this fold factor as insignificant. Finally, genes, we reduce the fold factor using the gene's divergence factor. In the first round, we simply count the number of genes whose fold factor is more than max_box_span (for computing boxes) and max_box_span + max_neighborhood_span (for computing neighborhoods). In the followup rounds, we use the maximal gene fold, for genes that are correlated in the vicinity of the metacells (see below). Finally, in the followup rounds, we only compare a base metacell with other metacells which were in the same neighborhood in the previous round.\nWe use hierarchical clustering to partition the metacells to distinct boxes, such that the maximal distance between any metacells in the box is bounded. In the first round, this bound is the max_deviant_genes_fraction out of the total number of genes. In the followup rounds, this is the max_box_span.\nFor each box, we compute a main neighborhood of other boxes such that the maximal distance between any metacells in the neighborhood is bounded. In the first round, this bound is again the maximal number of deviant genes (this time, using the increased fold distance computed above). In the followup rounds, this is the max_box_span plus the max_neighborhood_span. These neighborhoods may overlap. The main neighborhoods of different boxes may even be identical.\nFor each box, we compute the set of genes which are correlated (using the correlation_confidence) with some other gene(s) in its main neighborhood.\nIf the new sets of correlated genes only differ up to max_convergence_fraction, consider this to have converged. Otherwise, repeat from step 2.\n\nIf overwrite is set, the results will replace any previously computed boxes and neighborhoods.\n\nInputs\n\nAxes\n\ngene (required): Sequenced genes.\n\nmetacell (required): Minimal-sized groups of cells for robust point estimates.\n\nVectors\n\ngene @ divergence::AbstractFloat (required): Scale fold factors of each gene by multiplying with (1 - divergence) of the gene.\n\nmetacell @ total_UMIs::Unsigned (required): The total number of UMIs used to estimate the fraction of all the genes in each metacell.\n\nMatrices\n\ngene, metacell @ fraction::AbstractFloat (required): The estimated fraction of the UMIs of each gene in each metacell.\n\ngene, metacell @ total_UMIs::Unsigned (required): The total number of UMIs used to estimate the fraction of each gene in each metacell.\n\nOutputs\n\nAxes\n\nbox (guaranteed): Distinct groups of metacells with \"very close\" estimated cell state.\n\nneighborhood (guaranteed): Overlapping groups of boxes with \"close\" estimated cell states.\n\nVectors\n\nmetacell @ box::AbstractString (guaranteed): The unique box each metacell belongs to.\n\nbox @ neighborhood.main::AbstractString (guaranteed): The unique main neighborhood of each box.\n\nneighborhood @ span::AbstractFloat (guaranteed): The span (fold factor) used to compute the neighborhood.\n\nMatrices\n\ngene, neighborhood @ is_correlated::Bool (guaranteed): Which genes are correlated in each neighborhood.\n\nbox, box @ distance::AbstractFloat (guaranteed): The distance (fold factor) between the most different metacell genes between the boxes.\n\nbox, neighborhood @ is_member::Bool (guaranteed): A mask of the member boxes of each neighborhood.\n\n\n\n\n\n","category":"function"},{"location":"boxes.html#Index","page":"Boxes","title":"Index","text":"","category":"section"},{"location":"boxes.html","page":"Boxes","title":"Boxes","text":"Pages = [\"boxes.md\"]","category":"page"},{"location":"anndata_format.html#AnnData-Format","page":"AnnData Format","title":"AnnData Format","text":"","category":"section"},{"location":"anndata_format.html","page":"AnnData Format","title":"AnnData Format","text":"Metacells.AnnDataFormat\nMetacells.AnnDataFormat.import_h5ads!\nMetacells.AnnDataFormat.CopyAnnData","category":"page"},{"location":"anndata_format.html#Metacells.AnnDataFormat","page":"AnnData Format","title":"Metacells.AnnDataFormat","text":"Import and export metacells data from/to h5ad files. This allows moving data between the old Python/C++ based AnnData world and the new Julia based Daf world.\n\n\n\n\n\n","category":"module"},{"location":"anndata_format.html#Metacells.AnnDataFormat.import_h5ads!","page":"AnnData Format","title":"Metacells.AnnDataFormat.import_h5ads!","text":"function import_h5ads!(\n    destination::DafWriter;\n    raw_cells_h5ad::Maybe{AbstractString} = nothing,\n    clean_cells_h5ad::AbstractString,\n    metacells_h5ad::AbstractString,\n    copy_clean_data::Maybe{CopyAnnData} = nothing,\n    type_property::Maybe{AbstractString} = nothing,\n    rename_type::Maybe{AbstractString} = \"type\",\n    type_colors_csv::Maybe{AbstractString} = nothing,\n    type_properties::Maybe{AbstractSet{<:AbstractString}} = nothing,\n    properties_defaults::Maybe{Dict} = nothing,\n)::Nothing\n\nImport an AnnData based metacells dataset into a Daf destination data set. Ideally, the input must include clean_cells_h5ad and the metacells_h5ad computed for them, and optionally also the raw_cells_h5ad including the excluded cells and genes.\n\nIf type annotations were assigned to the metacells, then the name of the type_property should be specified. This can be further enhanced by specifying a type_colors_csv file mapping type names to colors. This should be a comma or tab separated file containing at least two columns, one named \"color\" and one with the same name as the type_property. For consistency, by default the type_property is renamed to the value of rename_type (by default, \"type\"). You can disable this by setting rename_type to nothing. We also call reconstruct_axis! to build the type axis; you can therefore specify an empty_type name, which will be converted to the empty string, to match the Daf convention of \"no value\" for string data, and specify an explicit set of type_properties (by default, any per-metacell property that has the same value for all metacells of each type will be converted to a type property) and properties_defaults.\n\nThis will mostly just read all the specified h5ad files and copy the data into the destination, with the following changes to match the Daf capabilities and conventions:\n\nThe X matrix of the cells is renamed to UMIs, and the X matrix of the metacells is renamed to fraction.\nMatrices and vectors of counts (UMIs, zeros) or module indices are converted to an unsigned type.\nThe __name__ scalar is not copied.\nThe excluded_gene and excluded_cell masks are not copied. Instead, if raw_cells_h5ad is specified, an is_excluded mask is created for both cells and genes, marking these that exist only in the raw_cells_h5ad and not in clean_cells_h5ad and metacells_h5ad.\nThe full_gene_index is not copied.\nThe properly_sampled_gene mask is renamed to the per-gene is_properly_sampled mask.\nThe bursty_lonely_gene mask is renamed to the per-gene is_bursty_lonely mask.\nThe lateral_gene mask is renamed to the per-gene is_lateral mask.\nThe noisy_gene mask is renamed to the per-gene is_noisy mask.\nThe rare_gene mask is renamed to the per-gene is_rare mask.\nThe rare_gene_module has 1 added to it (that is, \"no module\" is 0 in Daf) and is renamed to rare_module.\nThe lateral_genes_module has 1 added to it (that is, \"no module\" is 0 in Daf) and is renamed to lateral_module.\nThe marker_gene mask is renamed to the per-gene is_marker mask.\nThe selected_gene mask is renamed to the per-gene is_selected mask.\nThe ignored_gene mask is renamed to the per-gene is_ignored mask.\nThe ignored_gene_of_<type> masks are converted to an is_ignored mask per-gene-per-type.\nThe projected_noisy_gene mask is renamed to the per-gene is_projected_noisy mask.\nThe atlas_gene, atlas_lateral_gene, atlas_noisy_gene, atlas_marker_gene masks are renamed to the is_atlas, is_atlas_lateral, is_atlas_noisy and is_atlas_marker per-gene masks.\nThe essential_gene_of_<type> masks are converted to an is_essential mask per-gene-per-type.\nThe atlas_essential_gene_of_<type> masks are converted to an is_atlas_essential mask per-gene-per-type.\nThe fitted_gene_of_<type> masks are converted to an is_fitted mask per-gene-per-type.\nThe fitted mask per-gene-per-metacell is renamed to is_fitted.\nThe misfit mask per-gene-per-metacell is renamed to is_misfit.\nThe essential mask per-gene-per-metacell is renamed to is_essential.\nThe full_cell_index is not copied.\nThe properly_sampled_cell mask is renamed to the per-cell is_properly_sampled mask.\nThe rare_cell mask is renamed to the per-cell is_rare mask.\nThe cells_rare_gene_module has 1 added to it (that is, \"no module\" is 0 in Daf) and is renamed to rare_gene_module.\nThe per-cell dissolve mask is renamed to is_dissolved.\nThe per-cell metacell integer annotation is not copied, and the metacell_name string annotation is renamed to metacell.\nThe per-cell most_similar integer annotation is not copied, and the most_similar_name string annotation is renamed to metacell.most_similar.\nThe rare_metacell mask is renamed to the per-metacell is_rare mask.\nThe per-metacell metacells_level is renamed to level.\nThe per-metacell similar mask is renamed to is_similar.\n\nnote: Note\nThere is much duplication of data between the three h5ad files (in particular, per-gene data). Data in raw_cells_h5ad will override data in clean_cells_h5ad, which will override data in metacells_h5ad.\n\nData that exists only in clean_cells_h5ad poses a question when being copied into the full data set, which includes the full raw set of cells and genes. If copy_clean_data is nothing (the default), this is simply an error. Otherwise, data that is listed in copy_clean_data is copied using the specified name and the default value is applied to the raw-only genes or cells.\n\nnote: Note\nIt is common to call reconstruct_axis! on the result (e.g., if the cells were collected from a set of batches).\n\n\n\n\n\n","category":"function"},{"location":"anndata_format.html#Metacells.AnnDataFormat.CopyAnnData","page":"AnnData Format","title":"Metacells.AnnDataFormat.CopyAnnData","text":"Specify how to copy data from AnnData to Daf. The key is simply a vector or matrix name (ignoring axes), and the value is either nothing to ignore the data, or a tuple with the name of the destination Daf property and an optional value to use for missing entries (raw-only cells and/or genes).\n\n\n\n\n\n","category":"type"},{"location":"anndata_format.html#Index","page":"AnnData Format","title":"Index","text":"","category":"section"},{"location":"anndata_format.html","page":"AnnData Format","title":"AnnData Format","text":"Pages = [\"anndata_format.md\"]","category":"page"},{"location":"index.html#Metacells","page":"Metacells","title":"Metacells","text":"","category":"section"},{"location":"index.html","page":"Metacells","title":"Metacells","text":"Metacells.Metacells","category":"page"},{"location":"index.html#Metacells.Metacells","page":"Metacells","title":"Metacells.Metacells","text":"The Metacells.jl package provides computational services for the metacells package, using Daf to hold the data. In the future, we'll ideally migrate all of the metacellspackage computations to this package, converting the Python package to a thin wrapper, and provide a similar thin R wrapper to provide metacell analysis from R as well. For now,Metacells.jlonly provides a subset of the features of the Pythonmetacellspackage, which requires users to convert data fromAnnData(for the old features) toDaf (to the new features).\n\n\n\n\n\n","category":"module"},{"location":"index.html#Index","page":"Metacells","title":"Index","text":"","category":"section"},{"location":"index.html","page":"Metacells","title":"Metacells","text":"","category":"page"},{"location":"contracts.html#Contracts","page":"Contracts","title":"Contracts","text":"","category":"section"},{"location":"contracts.html","page":"Contracts","title":"Contracts","text":"Metacells.Contracts","category":"page"},{"location":"contracts.html#Metacells.Contracts","page":"Contracts","title":"Metacells.Contracts","text":"Functions for defining a Contract for a metacells Daf data set @computation. This also serves as a quick vocabulary of the contents of a metacells Daf data set.\n\nIn the descriptions below, \"fold factor\" refers to the log base 2 of the ratio between expression levels. For gene RNA expression levels, we typically use a regularization factor of 1e-5; fold factors of 1 (2x) is typically considered to be irrelevant, fold factors of 2 (4x) are considered to be potentially significant, fold factor of 3 (8x) are considered to point to a definite significant difference, etc. When used as thresholds, fold factors are adjusted using confidence levels (based on the number of UMIs used for the estimates) and the gene_divergence_vector.\n\n\n\n\n\n","category":"module"},{"location":"contracts.html#Axes","page":"Contracts","title":"Axes","text":"","category":"section"},{"location":"contracts.html","page":"Contracts","title":"Contracts","text":"Metacells.Contracts.gene_axis\nMetacells.Contracts.cell_axis\nMetacells.Contracts.metacell_axis\nMetacells.Contracts.box_axis\nMetacells.Contracts.neighborhood_axis\nMetacells.Contracts.type_axis","category":"page"},{"location":"contracts.html#Metacells.Contracts.gene_axis","page":"Contracts","title":"Metacells.Contracts.gene_axis","text":"function gene_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of sequenced genes. By convention we use TODO as the namespace of the genes, but this may be different depending on the data set.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.cell_axis","page":"Contracts","title":"Metacells.Contracts.cell_axis","text":"function cell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of sequenced single cells. There's no convention for cell names, as long as they are unique. Typically some sort of barcode is used, possibly combined with a batch and/or plate and/or experiment identification. In the latter case it is recommended that batch and/or plate and/or experiment would also be created as explicit axes, to allow associating metadata with them instead of repeating it for each cell.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.metacell_axis","page":"Contracts","title":"Metacells.Contracts.metacell_axis","text":"function metacell_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of metacells, which are minimal-sized groups of cells for robust point estimates. That is, each metacell is considered to be a robustly estimated point in the multi-dimensional manifold of cell states. Metacells may be very similar or very different from each other depending on the data set and the manifold region.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_axis","page":"Contracts","title":"Metacells.Contracts.box_axis","text":"function box_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of boxes, which are distinct groups of metacells with \"very close\" estimated cell state. That is, for some chosen set of genes, the metacells in each box all have very close estimates of gene expressions (maximal fold factor up to some maximal span). This maximal span is small and identical for all the boxes.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.neighborhood_axis","page":"Contracts","title":"Metacells.Contracts.neighborhood_axis","text":"function neighborhood_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of neighborhoods, which are overlapping groups of boxes with \"close\" estimated cell states. That is, for some chosen set of genes, the metacells in all the boxes of the neighborhood will all have close estimates of gene expression (maximal fold factor up to some maximal span). This maximal span is moderate and different for each neighborhood to ensure neighborhood sizes are not too big.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.type_axis","page":"Contracts","title":"Metacells.Contracts.type_axis","text":"function type_axis(expectation::ContractExpectation)::Pair{AxisKey, AxisSpecification}\n\nThe axis of types, which are distinct named biological cell states. That is, types are convenient labels manually assigned to large groups of cells, metacells, boxes, neighborhoods, etc. The resolution of the type labels depends on the data set and the type of analysis. In particular, types are not typically associated with a specific biological cell state, but rather with a set of related biological cell states (possibly along a gradient of such states).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Vectors","page":"Contracts","title":"Vectors","text":"","category":"section"},{"location":"contracts.html","page":"Contracts","title":"Contracts","text":"Metacells.Contracts.gene_is_excluded_vector\nMetacells.Contracts.gene_is_lateral_vector\nMetacells.Contracts.gene_is_marker_vector\nMetacells.Contracts.gene_is_correlated_vector\nMetacells.Contracts.gene_divergence_vector\nMetacells.Contracts.cell_is_excluded_vector\nMetacells.Contracts.metacell_total_UMIs_vector\nMetacells.Contracts.metacell_box_vector\nMetacells.Contracts.metacell_type_vector\nMetacells.Contracts.box_main_neighborhood_vector\nMetacells.Contracts.box_total_UMIs_vector\nMetacells.Contracts.box_type_vector\nMetacells.Contracts.neighborhood_span_vector\nMetacells.Contracts.type_color_vector","category":"page"},{"location":"contracts.html#Metacells.Contracts.gene_is_excluded_vector","page":"Contracts","title":"Metacells.Contracts.gene_is_excluded_vector","text":"function gene_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA mask of genes that are excluded from consideration. These genes are completely unrelated to the biological behaviors of interest. Not only that, they have strong and variable expression levels; enough to have an global impact on the expression level of the rest of the genes - for example, mitochondrial genes. Such genes make it difficult to estimate the relative expression level of genes between different cell states. Therefore, such genes aren't even counted in the total UMIs of each cell.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_is_lateral_vector","page":"Contracts","title":"Metacells.Contracts.gene_is_lateral_vector","text":"function gene_is_lateral_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA mask of genes that are lateral to the biological behaviors of interest. That is, these genes may satisfy all criteria for being in a group of cooperating genes, but the biological behavior they participate in isn't relevant to the behaviors of interest - for example, genes related to cell cycle or stress. Such genes make it harder to focus on the biological behaviors of interest. They are therefore masked out during the analysis.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_is_marker_vector","page":"Contracts","title":"Metacells.Contracts.gene_is_marker_vector","text":"function gene_is_marker_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA mask of genes that distinguish between cell states. That is, these genes have a significant expression level at some cell state, as well as a significant range of expression across all cell states, so can be used to distinguish between cell states. Non-marker genes are by definition not useful for such analysis, but marker genes aren't necessarily useful due to other considerations.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_is_correlated_vector","page":"Contracts","title":"Metacells.Contracts.gene_is_correlated_vector","text":"function gene_is_correlated_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA mask of genes that are correlated with other gene(s). We typically search for groups of genes that act together. Genes that have no correlation with other genes aren't useful for this sort of analysis, even if they are marker genes.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_divergence_vector","page":"Contracts","title":"Metacells.Contracts.gene_divergence_vector","text":"function gene_divergence_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nFor each gene, we scale fold factors of each gene by multiplying with (1 - divergence) of the gene. In particular this is used when considering the distance between gene expressions. Therefore genes with a non-zero divergence will require a higher fold factor\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.cell_is_excluded_vector","page":"Contracts","title":"Metacells.Contracts.cell_is_excluded_vector","text":"function cell_is_excluded_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA mask of cells that are excluded from consideration. This can be due to any number of reasons - doublets, too low a number of UMIs, to high a percentage of excluded gene UMIs, etc.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.metacell_total_UMIs_vector","page":"Contracts","title":"Metacells.Contracts.metacell_total_UMIs_vector","text":"function metacell_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe total number of UMIs used to estimate the fraction of all the genes in each metacell. This is used to estimate the robustness of the estimates.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.metacell_box_vector","page":"Contracts","title":"Metacells.Contracts.metacell_box_vector","text":"function metacell_box_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe unique box each metacell belongs to.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.metacell_type_vector","page":"Contracts","title":"Metacells.Contracts.metacell_type_vector","text":"function metacell_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe type each metacell belongs to.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_main_neighborhood_vector","page":"Contracts","title":"Metacells.Contracts.box_main_neighborhood_vector","text":"function box_main_neighborhood_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe unique main neighborhood of each box. Ideally, each box is the center of its own main neighborhood, and also belongs to overlapping neighborhood of some other nearby boxes. However, where the manifold is sparsely sampled, a few nearby boxes may share the same main neighborhood. If the samples are sufficiently sparse, the main neighborhood may include only just the single box (which itself may include just a single metacell).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_total_UMIs_vector","page":"Contracts","title":"Metacells.Contracts.box_total_UMIs_vector","text":"function box_total_UMIs_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe unique main neighborhood of each box. Ideally, each box is the center of its own main neighborhood, and also belongs to overlapping neighborhood of some other nearby boxes. However, where the manifold is sparsely sampled, a few nearby boxes may share the same main neighborhood. If the samples are sufficiently sparse, the main neighborhood may include only just the single box (which itself may include just a single metacell).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_type_vector","page":"Contracts","title":"Metacells.Contracts.box_type_vector","text":"function box_type_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe type each box belongs to.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.neighborhood_span_vector","page":"Contracts","title":"Metacells.Contracts.neighborhood_span_vector","text":"function neighborhood_span_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nThe span (fold factor) used to compute the neighborhood. This is different for each neighborhood (but is never too big), to ensure that the sizes of the neighborhoods are not too large even for densely sampled regions of the manifold (a simple fixed span will not do, due to the curse of multi-dimensionality).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.type_color_vector","page":"Contracts","title":"Metacells.Contracts.type_color_vector","text":"function type_color_vector(expectation::ContractExpectation)::Pair{VectorKey, DataSpecification}\n\nA unique color for each type for graphs.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Matrices","page":"Contracts","title":"Matrices","text":"","category":"section"},{"location":"contracts.html","page":"Contracts","title":"Contracts","text":"Metacells.Contracts.gene_metacell_fraction_matrix\nMetacells.Contracts.gene_metacell_total_UMIs_matrix\nMetacells.Contracts.gene_neighborhood_is_correlated_matrix\nMetacells.Contracts.box_box_distance\nMetacells.Contracts.gene_box_fraction_matrix\nMetacells.Contracts.gene_box_total_UMIs_matrix\nMetacells.Contracts.box_neighborhood_is_member_matrix","category":"page"},{"location":"contracts.html#Metacells.Contracts.gene_metacell_fraction_matrix","page":"Contracts","title":"Metacells.Contracts.gene_metacell_fraction_matrix","text":"function gene_metacell_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nThe estimated fraction of the UMIs of each gene in each metacell. We assume that each metacell is a sample of the manifold, representing a real biological state, regardless to its distance to other metacells (subject to cleaning up batch effects, purging doublets, and compensating for any other technical artifacts).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_metacell_total_UMIs_matrix","page":"Contracts","title":"Metacells.Contracts.gene_metacell_total_UMIs_matrix","text":"function gene_metacell_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nThe total number of UMIs used to estimate the fraction of each gene in each metacell. This is used to estimate the robustness of the estimate. When computing fold factors, we require the total number of UMIs (from both compared estimates) to be some minimum, and possibly adjust the fold factor according to some confidence level (assuming a multinomial sampling distribution).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_neighborhood_is_correlated_matrix","page":"Contracts","title":"Metacells.Contracts.gene_neighborhood_is_correlated_matrix","text":"function gene_neighborhood_is_correlated_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nWhich genes are correlated in each neighborhood. This is a smaller set than the set of globally correlated genes, since some genes may be a strong marker for some cell type (and thus have a strong correlation with other genes specific to this type), but may lack correlation with any other genes when considering only cell states of this type. As such, they hamper analysis of groups of cooperating genes within this specific cell type. However, instead of relying on manual type annotations, we make use of the computed neighborhoods to obtain the set of correlated genes in each local region of the manifold.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_box_distance","page":"Contracts","title":"Metacells.Contracts.box_box_distance","text":"function box_box_distance(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nThe distance (fold factor) between the most different metacell genes between the boxes. This is the fold factor between the most different gene expression in a pair of metacells, one in each box. This considers only the chosen genes (marker genes that are also correlated in the main neighborhood of either of the boxes).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_box_fraction_matrix","page":"Contracts","title":"Metacells.Contracts.gene_box_fraction_matrix","text":"function gene_box_fraction_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nThe estimated fraction of the UMIs of each gene in each box. Each box is a sample of the manifold, representing a real biological state, which is different from the state of any other box.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.gene_box_total_UMIs_matrix","page":"Contracts","title":"Metacells.Contracts.gene_box_total_UMIs_matrix","text":"function gene_box_total_UMIs_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nThe total number of UMIs used to estimate the fraction of each gene in each box. This is used to estimate the robustness of the estimate. When computing fold factors, we require the total number of UMIs (from both compared estimates) to be some minimum, and possibly adjust the fold factor according to some confidence level (assuming a multinomial sampling distribution).\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Metacells.Contracts.box_neighborhood_is_member_matrix","page":"Contracts","title":"Metacells.Contracts.box_neighborhood_is_member_matrix","text":"function box_neighborhood_is_member_matrix(expectation::ContractExpectation)::Pair{MatrixKey, DataSpecification}\n\nA mask of the member boxes of each neighborhood. This is needed since each box may belong to multiple neighborhoods, and each neighborhood contains multiple boxes.\n\n\n\n\n\n","category":"function"},{"location":"contracts.html#Index","page":"Contracts","title":"Index","text":"","category":"section"},{"location":"contracts.html","page":"Contracts","title":"Contracts","text":"Pages = [\"contracts.md\"]","category":"page"}]
}
