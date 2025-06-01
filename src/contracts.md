# Contracts

```@docs
Metacells.Contracts
```

## Axes

```@docs
Metacells.Contracts.gene_axis
Metacells.Contracts.cell_axis
Metacells.Contracts.metacell_axis
Metacells.Contracts.block_axis
Metacells.Contracts.type_axis
Metacells.Contracts.component_axis
Metacells.Contracts.module_axis
```

## Vectors

### Gene Vectors

```@docs
Metacells.Contracts.gene_is_excluded_vector
Metacells.Contracts.gene_is_lateral_vector
Metacells.Contracts.gene_is_marker_vector
Metacells.Contracts.gene_marker_rank_vector
Metacells.Contracts.gene_is_uncorrelated_vector
Metacells.Contracts.gene_is_covered_vector
Metacells.Contracts.gene_is_skeleton_vector
Metacells.Contracts.gene_is_forbidden_vector
Metacells.Contracts.gene_is_regulator_vector
Metacells.Contracts.gene_is_transcription_factor_vector
Metacells.Contracts.gene_divergence_vector
```

### Cell Vectors

#### Cell UMIs Vectors

```@docs
Metacells.Contracts.cell_excluded_UMIs_vector
Metacells.Contracts.cell_total_UMIs_vector
Metacells.Contracts.cell_ribosomal_UMIs_vector
Metacells.Contracts.cell_mitochondrial_UMIs_vector
Metacells.Contracts.cell_covered_UMIs_vector
```

#### Cell Metadata Vectors

```@docs
Metacells.Contracts.cell_is_excluded_vector
Metacells.Contracts.cell_metacell_vector
Metacells.Contracts.cell_type_vector
Metacells.Contracts.cell_linear_metacell_cross_entropy_vector
Metacells.Contracts.cell_linear_covered_metacell_cross_entropy_vector
Metacells.Contracts.cell_linear_metacell_kl_divergence_vector
Metacells.Contracts.cell_linear_covered_metacell_kl_divergence_vector
```

### Metacell Vectors

#### Metacell UMIs Vectors

```@docs
Metacells.Contracts.metacell_total_UMIs_vector
Metacells.Contracts.metacell_scaled_total_UMIs_vector
Metacells.Contracts.metacell_covered_UMIs_vector
Metacells.Contracts.metacell_scaled_covered_UMIs_vector
```

#### Metacell Counts Vectors

```@docs
Metacells.Contracts.metacell_n_cells_vector
```

#### Metacell Metadata Vectors

```@docs
Metacells.Contracts.metacell_block_vector
Metacells.Contracts.metacell_type_vector
Metacells.Contracts.metacell_mean_cells_linear_cross_entropy_vector
Metacells.Contracts.metacell_mean_cells_linear_covered_cross_entropy_vector
Metacells.Contracts.metacell_mean_cells_linear_kl_divergence_vector
Metacells.Contracts.metacell_mean_cells_linear_covered_kl_divergence_vector
```

### Block Vectors

#### Block UMIs Vectors

```@docs
Metacells.Contracts.block_total_UMIs_vector
Metacells.Contracts.block_scaled_total_UMIs_vector
Metacells.Contracts.block_covered_UMIs_vector
Metacells.Contracts.block_scaled_covered_UMIs_vector
```

#### Block Counts Vectors

```@docs
Metacells.Contracts.block_n_cells_vector
Metacells.Contracts.block_n_metacells_vector
Metacells.Contracts.block_n_neighborhood_blocks_vector
Metacells.Contracts.block_n_environment_blocks_vector
Metacells.Contracts.block_n_neighborhood_metacells_vector
Metacells.Contracts.block_n_environment_metacells_vector
Metacells.Contracts.block_n_neighborhood_cells_vector
Metacells.Contracts.block_n_environment_cells_vector
```

#### Block Metadata Vectors

```@docs
Metacells.Contracts.block_type_vector
```

#### Gene Components Analysis

```@docs
Metacells.Contracts.block_n_found_components_vector
Metacells.Contracts.block_components_neighborhood_RMSE_vector
Metacells.Contracts.block_components_neighborhood_XRMSE_vector
```

#### Gene Modules Analysis

```@docs
Metacells.Contracts.block_n_found_modules_vector
Metacells.Contracts.block_n_used_modules_vector
Metacells.Contracts.block_modules_neighborhood_RMSE_vector
Metacells.Contracts.block_modules_neighborhood_XRMSE_vector
Metacells.Contracts.block_modules_block_RMSE_vector
```

### Type Vectors

```@docs
Metacells.Contracts.type_color_vector
```

## Matrices

### Cell Matrices

```@docs
Metacells.Contracts.cell_gene_UMIs_matrix
```

### Metacell Matrices

#### Metacell UMIs

```@docs
Metacells.Contracts.metacell_gene_UMIs_matrix
```

#### Metacell Geomean Fractions

```@docs
Metacells.Contracts.metacell_gene_geomean_fraction_matrix
Metacells.Contracts.metacell_gene_log_geomean_fraction_matrix
```

#### Metacell Linear Fractions

```@docs
Metacells.Contracts.metacell_gene_linear_fraction_matrix
Metacells.Contracts.metacell_gene_scaled_linear_fraction_matrix
Metacells.Contracts.metacell_gene_log_linear_fraction_matrix
Metacells.Contracts.metacell_gene_log_scaled_linear_fraction_matrix
```

#### Metacell Linear Covered Fractions

```@docs
Metacells.Contracts.metacell_gene_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_scaled_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_log_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_log_scaled_linear_covered_fraction_matrix
```

#### Metacell Virtual Fractions

```@docs
Metacells.Contracts.metacell_gene_fraction_matrix
Metacells.Contracts.metacell_gene_log_fraction_matrix
Metacells.Contracts.metacell_gene_covered_fraction_matrix
Metacells.Contracts.metacell_gene_log_covered_fraction_matrix
```

#### Metacell Distances

```@docs
Metacells.Contracts.metacell_metacell_max_skeleton_fold_distance
```

#### Gene Modules Analysis

##### Gene Modules UMIs

```@docs
Metacells.Contracts.metacell_module_total_UMIs_matrix
Metacells.Contracts.metacell_module_scaled_total_UMIs_matrix
Metacells.Contracts.metacell_module_covered_UMIs_matrix
Metacells.Contracts.metacell_module_scaled_covered_UMIs_matrix
```

##### Gene Modules Linear Fractions

```@docs
Metacells.Contracts.metacell_module_linear_fraction_matrix
Metacells.Contracts.metacell_module_scaled_linear_fraction_matrix
Metacells.Contracts.metacell_module_log_linear_fraction_matrix
Metacells.Contracts.metacell_module_log_scaled_linear_fraction_matrix
```

##### Gene Modules Linear Covered Fractions

```@docs
Metacells.Contracts.metacell_module_linear_covered_fraction_matrix
Metacells.Contracts.metacell_module_scaled_linear_covered_fraction_matrix
Metacells.Contracts.metacell_module_log_linear_covered_fraction_matrix
Metacells.Contracts.metacell_module_log_scaled_linear_covered_fraction_matrix
```

##### Gene Modules Approximated Covered Fractions

```@docs
Metacells.Contracts.metacell_gene_approximated_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_approximated_scaled_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_approximated_log_linear_covered_fraction_matrix
Metacells.Contracts.metacell_gene_approximated_log_scaled_linear_covered_fraction_matrix
```

##### Gene Modules Virtual Fractions

```@docs
Metacells.Contracts.metacell_module_fraction_matrix
Metacells.Contracts.metacell_module_log_fraction_matrix
Metacells.Contracts.metacell_module_covered_fraction_matrix
Metacells.Contracts.metacell_module_log_covered_fraction_matrix
Metacells.Contracts.metacell_gene_approximated_covered_fraction_matrix
Metacells.Contracts.metacell_gene_approximated_log_covered_fraction_matrix
```

### Block Matrices

#### Block Distances

```@docs
Metacells.Contracts.block_block_max_skeleton_fold_distance
```

#### Block Vicinities

```@docs
Metacells.Contracts.block_block_is_in_neighborhood_matrix
Metacells.Contracts.block_block_is_in_environment_matrix
```

#### Block UMIs

```@docs
Metacells.Contracts.block_gene_UMIs_matrix
```

#### Block Linear Fractions

```@docs
Metacells.Contracts.block_gene_linear_fraction_matrix
Metacells.Contracts.block_gene_scaled_linear_fraction_matrix
Metacells.Contracts.block_gene_log_linear_fraction_matrix
Metacells.Contracts.block_gene_log_scaled_linear_fraction_matrix
```

#### Block Linear Covered Fractions

```@docs
Metacells.Contracts.block_gene_linear_covered_fraction_matrix
Metacells.Contracts.block_gene_scaled_linear_covered_fraction_matrix
Metacells.Contracts.block_gene_log_linear_covered_fraction_matrix
Metacells.Contracts.block_gene_log_scaled_linear_covered_fraction_matrix
```

#### Block Virtual Fractions

```@docs
Metacells.Contracts.block_gene_fraction_matrix
Metacells.Contracts.block_gene_log_fraction_matrix
Metacells.Contracts.block_gene_covered_fraction_matrix
Metacells.Contracts.block_gene_log_covered_fraction_matrix
```

#### Block Gene Masks

```@docs
Metacells.Contracts.block_gene_is_environment_marker_matrix
Metacells.Contracts.block_gene_environment_marker_rank_matrix
Metacells.Contracts.block_gene_is_environment_skeleton_matrix
```

#### Gene Modules Analysis

##### Gene Modules Counts

```@docs
Metacells.Contracts.block_module_n_genes_matrix
Metacells.Contracts.block_module_n_covered_matrix
```

##### Gene Modules UMIs

```@docs
Metacells.Contracts.block_module_total_UMIs_matrix
Metacells.Contracts.block_module_scaled_total_UMIs_matrix
Metacells.Contracts.block_module_covered_UMIs_matrix
Metacells.Contracts.block_module_scaled_covered_UMIs_matrix
```

##### Gene Modules Masks

```@docs
Metacells.Contracts.block_gene_module_index_matrix
Metacells.Contracts.block_module_is_found_matrix
Metacells.Contracts.block_module_is_used_matrix
```

##### Gene Modules Metadata

```@docs
Metacells.Contracts.block_module_min_gene_correlation_matrix
```

##### Gene Components Model

```@docs
Metacells.Contracts.block_gene_component_index_matrix
Metacells.Contracts.block_component_is_found_matrix
Metacells.Contracts.block_component_n_genes_matrix
Metacells.Contracts.block_component_base_covered_fraction_matrix
Metacells.Contracts.block_component_gene_covered_coefficient_tensor
```

##### Gene Modules Model

```@docs
Metacells.Contracts.block_gene_base_covered_fraction_matrix
Metacells.Contracts.block_module_base_covered_fraction_matrix
Metacells.Contracts.block_module_gene_covered_coefficient_tensor
```

## Index

```@index
Pages = ["contracts.md"]
```
