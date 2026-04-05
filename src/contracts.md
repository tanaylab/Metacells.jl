# Contracts

```@docs
Metacells.Contracts
```

## Genes

```@docs
Metacells.Contracts.gene_axis
Metacells.Contracts.vector_of_is_mitochondrial_per_gene
Metacells.Contracts.vector_of_is_ribosomal_per_gene
Metacells.Contracts.vector_of_is_excluded_per_gene
Metacells.Contracts.vector_of_is_lateral_per_gene
Metacells.Contracts.vector_of_is_marker_per_gene
Metacells.Contracts.vector_of_marker_rank_per_gene
Metacells.Contracts.vector_of_is_transcription_factor_per_gene
Metacells.Contracts.vector_of_is_regulator_per_gene
Metacells.Contracts.vector_of_is_skeleton_per_gene
Metacells.Contracts.vector_of_is_forbidden_per_gene
Metacells.Contracts.vector_of_is_correlated_with_skeleton_per_gene
```

## Cells

```@docs
Metacells.Contracts.cell_axis
```

### Cells UMIs

```@docs
Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
Metacells.Contracts.vector_of_mitochondrial_UMIs_per_cell
Metacells.Contracts.vector_of_ribosomal_UMIs_per_cell
Metacells.Contracts.vector_of_excluded_UMIs_per_cell
Metacells.Contracts.vector_of_total_UMIs_per_cell
```

### Cells masks

```@docs
Metacells.Contracts.vector_of_is_excluded_per_cell
```

## Metacells

```@docs
Metacells.Contracts.metacell_axis
Metacells.Contracts.vector_of_metacell_per_cell
Metacells.Contracts.vector_of_n_cells_per_metacell
Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
Metacells.Contracts.vector_of_total_UMIs_per_metacell
```

### Metacells Fractions

```@docs
Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
```

### Metacells Distances

```@docs
Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
```

### Metacells Correlations

```@docs
Metacells.Contracts.vector_of_correlation_between_cells_and_punctuated_metacells_per_gene
```

### Metacells UMAP

```@docs
Metacells.Contracts.vector_of_umap_x_per_metacell
Metacells.Contracts.vector_of_umap_y_per_metacell
Metacells.Contracts.vector_of_umap_u_per_metacell
Metacells.Contracts.vector_of_umap_v_per_metacell
Metacells.Contracts.vector_of_umap_w_per_metacell
```

## Blocks

```@docs
Metacells.Contracts.block_axis
Metacells.Contracts.vector_of_block_per_metacell
Metacells.Contracts.vector_of_base_block_per_metacell
Metacells.Contracts.vector_of_n_metacells_per_block
Metacells.Contracts.vector_of_n_cells_per_block
```

### Blocks UMIs

```@docs
Metacells.Contracts.matrix_of_UMIs_per_gene_per_block
Metacells.Contracts.vector_of_total_UMIs_per_block
```

### Blocks Fractions

```@docs
Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_block
Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_block
```

### Blocks Distances

```@docs
Metacells.Contracts.matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks
```

### Blocks Confusion

```@docs
Metacells.Contracts.vector_of_block_closest_by_pertinent_markers_per_cell
Metacells.Contracts.matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block
```

### Blocks Neighborhoods

```@docs
Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
Metacells.Contracts.vector_of_n_neighborhood_blocks_per_block
Metacells.Contracts.vector_of_n_neighborhood_metacells_per_block
Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
Metacells.Contracts.matrix_of_is_neighborhood_marker_per_gene_per_block
Metacells.Contracts.matrix_of_is_neighborhood_distinct_per_gene_per_block
```

### Neighborhoods Correlations

```@docs
Metacells.Contracts.matrix_of_is_correlated_with_skeleton_in_neighborhood_per_gene_per_block
Metacells.Contracts.matrix_of_correlation_between_neighborhood_cells_and_punctuated_metacells_per_gene_per_block
```

### Blocks UMAP

```@docs
Metacells.Contracts.vector_of_umap_x_per_block
Metacells.Contracts.vector_of_umap_y_per_block
Metacells.Contracts.vector_of_umap_u_per_block
Metacells.Contracts.vector_of_umap_v_per_block
Metacells.Contracts.vector_of_umap_w_per_block
```

## Type Annotations

```@docs
Metacells.Contracts.type_axis
Metacells.Contracts.vector_of_color_per_type
Metacells.Contracts.vector_of_type_per_cell
Metacells.Contracts.vector_of_type_per_metacell
Metacells.Contracts.vector_of_type_per_block
```

## Gene Modules

```@docs
Metacells.Contracts.module_axis
Metacells.Contracts.vector_of_anchor_per_module
Metacells.Contracts.matrix_of_is_found_per_module_per_block
Metacells.Contracts.matrix_of_module_per_gene_per_block
Metacells.Contracts.matrix_of_module_status_per_gene_per_block
Metacells.Contracts.vector_of_n_modules_per_block
Metacells.Contracts.matrix_of_n_genes_per_module_per_block
Metacells.Contracts.tensor_of_linear_fraction_per_block_per_module_per_metacell
Metacells.Contracts.matrix_of_mean_linear_fraction_in_neighborhood_cells_per_module_per_block
Metacells.Contracts.matrix_of_std_linear_fraction_in_neighborhood_cells_per_module_per_block
Metacells.Contracts.vector_of_mean_euclidean_modules_cells_distance_per_metacell
Metacells.Contracts.vector_of_std_euclidean_modules_cells_distance_per_metacell
```

## Projection

```@docs
Metacells.Contracts.vector_of_projected_metacell_per_cell
Metacells.Contracts.vector_of_projected_block_per_cell
Metacells.Contracts.vector_of_projected_modules_z_score_per_cell
Metacells.Contracts.vector_of_correlation_between_cells_and_projected_metacells_per_gene
Metacells.Contracts.matrix_of_correlation_between_neighborhood_cells_and_projected_metacells_per_gene_per_projected_block
Metacells.Contracts.projected_block_axis
```

## Quality Control

```@docs
Metacells.Contracts.matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
Metacells.Contracts.matrix_of_correlation_with_most_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
Metacells.Contracts.matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block
Metacells.Contracts.base_block_axis
```

## Index

```@index
Pages = ["contracts.md"]
```
