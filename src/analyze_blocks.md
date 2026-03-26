# Analyze Blocks

```@docs
Metacells.AnalyzeBlocks
```

## Neighborhoods

```@docs
Metacells.AnalyzeBlocks.compute_matrix_of_is_in_neighborhood_per_block_per_block!
```

## Counts

```@docs
Metacells.AnalyzeBlocks.compute_vector_of_n_cells_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_cells_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_n_metacells_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_metacells_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_n_neighborhood_blocks_per_block!
```

## UMIs

```@docs
Metacells.AnalyzeBlocks.compute_matrix_of_UMIs_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_total_UMIs_per_block!
Metacells.AnalyzeBlocks.compute_vector_of_total_neighborhood_UMIs_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_linear_fraction_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_log_linear_fraction_per_gene_per_block!
```

## Genes

```@docs
Metacells.AnalyzeBlocks.compute_matrix_of_is_neighborhood_marker_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_is_neighborhood_distinct_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!
```

## Types

```@docs
Metacells.AnalyzeBlocks.compute_vector_of_type_per_block_by_metacells!
Metacells.AnalyzeBlocks.compute_vector_of_type_per_metacell_by_blocks!
Metacells.AnalyzeBlocks.compute_vector_of_type_per_block_by_cells!
Metacells.AnalyzeBlocks.compute_vector_of_type_per_cell_by_blocks!
```

## Distances

```@docs
Metacells.AnalyzeBlocks.compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!
Metacells.AnalyzeBlocks.compute_vector_of_block_closest_by_pertinent_markers_per_cell!
Metacells.AnalyzeBlocks.compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!
```

## Correlations

```@docs
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block!
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!
Metacells.AnalyzeBlocks.compute_matrix_of_correlation_with_most_between_base_neighborhood_cells_and_metacells_per_gene_per_base_block!
```

## UMAP

```@docs
Metacells.AnalyzeBlocks.compute_blocks_2d_umap_by_metacells!
Metacells.AnalyzeBlocks.compute_blocks_3d_umap_by_metacells!
```

## Index

```@index
Pages = ["analyze_blocks.md"]
```
