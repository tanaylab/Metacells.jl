# Contracts

```@docs
Metacells.Contracts
```

## Genes

```@docs
Metacells.Contracts.gene_axis
```

### Genes Masks

```@docs
Metacells.Contracts.gene_is_mitochondrial_vector
Metacells.Contracts.gene_is_ribosomal_vector
Metacells.Contracts.gene_is_excluded_vector
Metacells.Contracts.gene_is_lateral_vector
Metacells.Contracts.gene_is_marker_vector
Metacells.Contracts.gene_marker_rank_vector
Metacells.Contracts.gene_is_regulator_vector
Metacells.Contracts.gene_is_skeleton_vector
Metacells.Contracts.gene_is_transcription_factor_vector
Metacells.Contracts.gene_is_forbidden_vector
Metacells.Contracts.gene_is_correlated_with_skeleton_vector
```

## Cells

```@docs
Metacells.Contracts.cell_axis
```

### Cells UMIs

```@docs
Metacells.Contracts.cell_gene_UMIs_matrix
Metacells.Contracts.cell_mitochondrial_UMIs_vector
Metacells.Contracts.cell_ribosomal_UMIs_vector
Metacells.Contracts.cell_excluded_UMIs_vector
Metacells.Contracts.cell_total_UMIs_vector
```

### Cells masks

```@docs
Metacells.Contracts.cell_is_excluded_vector
```

## Metacells

```@docs
Metacells.Contracts.metacell_axis
Metacells.Contracts.cell_metacell_vector
Metacells.Contracts.metacell_n_cells_vector
```

### Metacells UMIs

```@docs
Metacells.Contracts.metacell_gene_UMIs_matrix
Metacells.Contracts.metacell_total_UMIs_vector
```

### Metacells Fractions

```@docs
Metacells.Contracts.metacell_gene_linear_fraction_matrix
Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
Metacells.Contracts.metacell_gene_geomean_fraction_matrix
Metacells.Contracts.metacell_gene_log_geomean_fraction_matrix
```

### Metacells Distances

```@docs
Metacells.Contracts.metacell_metacell_euclidean_skeleton_distance
Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
```

### Metacells Correlations

```@docs
Metacells.Contracts.cell_genes_correlation_with_metacells_vector
Metacells.Contracts.cell_markers_correlation_with_metacells_vector
Metacells.Contracts.cell_non_lateral_markers_correlation_with_metacells_vecto
Metacells.Contracts.metacell_mean_cells_genes_correlation_with_metacells_vector
Metacells.Contracts.metacell_mean_cells_markers_correlation_with_metacells_vector
Metacells.Contracts.metacell_mean_cells_non_lateral_markers_correlation_with_metacells_vecto
Metacells.Contracts.gene_correlation_between_cells_and_metacells_vector
Metacells.Contracts.metacell_gene_most_significant_correlated_gene_matrix
Metacells.Contracts.metacell_gene_most_significant_correlation_matrix
Metacells.Contracts.block_gene_most_correlated_gene_matrix
```

### Metacells Sharpening

```@docs
Metacells.Contracts.metacell_sharpening_rounds_vector
Metacells.Contracts.metacell_original_block_vector
Metacells.Contracts.metacell_radius
```

## Blocks

```@docs
Metacells.Contracts.block_axis
Metacells.Contracts.metacell_block_vector
Metacells.Contracts.block_n_metacells_vector
Metacells.Contracts.block_n_cells_vector
```

### Blocks UMIs

```@docs
Metacells.Contracts.block_gene_UMIs_matrix
Metacells.Contracts.block_total_UMIs_vector
```

### Blocks Fractions

```@docs
Metacells.Contracts.block_gene_linear_fraction_matrix
Metacells.Contracts.block_gene_log_linear_fraction_matrix
```

### Blocks Distances

```@docs
Metacells.Contracts.block_block_mean_euclidean_skeleton_distance
Metacells.Contracts.block_block_max_skeleton_fold_distance
```

### Blocks Neighborhoods

```@docs
Metacells.Contracts.block_block_is_in_neighborhood_matrix
Metacells.Contracts.block_n_neighborhood_blocks_vector
Metacells.Contracts.block_n_neighborhood_metacells_vector
Metacells.Contracts.block_n_neighborhood_cells_vector
Metacells.Contracts.block_gene_is_neighborhood_marker_matrix
Metacells.Contracts.block_gene_is_neighborhood_distinct_matrix
```

### Blocks Environments

```@docs
Metacells.Contracts.block_block_is_in_environment_matrix
Metacells.Contracts.block_n_environment_blocks_vector
Metacells.Contracts.block_n_environment_metacells_vector
Metacells.Contracts.block_n_environment_cells_vector
Metacells.Contracts.block_gene_is_environment_marker_matrix
```

### Blocks Correlations

```@docs
Metacells.Contracts.block_gene_correlation_between_neighborhood_cells_and_metacells_matrix
Metacells.Contracts.block_mean_cells_genes_correlation_with_metacells_vector
Metacells.Contracts.block_mean_cells_markers_correlation_with_metacells_vector
Metacells.Contracts.block_mean_cells_non_lateral_markers_correlation_with_metacells_vector
Metacells.Contracts.block_gene_is_correlated_with_skeleton_in_neighborhood_matrix
```

## Type Annotations

```@docs
Metacells.Contracts.type_axis
Metacells.Contracts.type_color_vector
Metacells.Contracts.metacell_type_vector
Metacells.Contracts.cell_type_vector
Metacells.Contracts.block_type_vector
```

## Modules

```@docs
Metacells.Contracts.module_axis
Metacells.Contracts.module_anchor_gene_vector
Metacells.Contracts.block_gene_module_matrix
Metacells.Contracts.block_module_is_found_matrix
Metacells.Contracts.block_module_is_strong_matrix
Metacells.Contracts.block_module_is_marker_in_environment_matrix
Metacells.Contracts.block_module_n_genes_matrix
Metacells.Contracts.block_module_n_skeletons_matrix
Metacells.Contracts.block_metacell_module_environment_total_UMIs_tensor
Metacells.Contracts.block_metacell_module_environment_linear_fraction_tensor
Metacells.Contracts.block_metacell_module_environment_log_linear_fraction_tensor
Metacells.Contracts.metacell_module_variance_over_mean_matrix
Metacells.Contracts.block_module_neighborhood_mean_linear_fraction_matrix
Metacells.Contracts.block_module_neighborhood_std_linear_fraction_matrix
```

## Chosens

```@docs
Metacells.Contracts.chosen_axis
Metacells.Contracts.chosen_block_vector
Metacells.Contracts.chosen_module_vector
Metacells.Contracts.block_module_chosen_matrix
Metacells.Contracts.gene_chosen_is_member_matrix
Metacells.Contracts.chosen_n_genes_vector
Metacells.Contracts.metacell_chosen_total_UMIs_matrix
Metacells.Contracts.metacell_chosen_linear_fraction_matrix
Metacells.Contracts.metacell_chosen_log_linear_fraction_matrix
Metacells.Contracts.metacell_chosen_variance_over_mean_matrix
Metacells.Contracts.block_chosen_variance_over_mean_matrix
```

## Explanations

```@docs
Metacells.Contracts.block_gene_most_correlated_gene_in_environment_matrix
Metacells.Contracts.block_gene_most_correlation_in_environment_matrix
Metacells.Contracts.block_gene_most_correlated_lateral_gene_in_environment_matrix
Metacells.Contracts.block_gene_most_correlation_with_lateral_in_environment_matrix
Metacells.Contracts.block_gene_most_correlated_regulator_gene_in_environment_matrix
Metacells.Contracts.block_gene_most_correlation_with_regulator_in_environment_matrix
Metacells.Contracts.block_gene_unexplained_correlation_in_environment_matrix
Metacells.Contracts.block_gene_is_unexplained_in_environment_matrix
Metacells.Contracts.block_gene_gene_correlation_in_environment_tensor
Metacells.Contracts.block_n_unexplained_genes_in_environment_vector
```

## Fitting

```@docs
Metacells.Contracts.block_gene_most_correlated_gene_in_neighborhood_cells_matrix
Metacells.Contracts.block_gene_gene_is_most_correlated_in_neighborhood_cells_tensor
Metacells.Contracts.block_gene_gene_most_correlation_in_neighborhood_cells_tensor
```

## Projection

```@docs
Metacells.Contracts.block_cell_correlation_matrix
```

## Index

```@index
Pages = ["contracts.md"]
```
