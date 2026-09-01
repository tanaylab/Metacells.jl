"""
Run whole stages of the metacells pipeline.

The functions here run multiple computation functions in sequence. They allow performing the computational pipeline of
sharpening metacells with few high-level calls instead of working through each low-level step on its own. None of the
functions here do anything other than invoking the lower level computations. As such they are purely convenience
functions. One is free to call the lower level steps directly in case the pipeline needs to be tweaked for some special
case.
"""
module Pipeline

export analyze_metacells!
export import_base_metacells!
export prepare_metacells!

using DataAxesFormats
using Random
using TanayLabUtilities

using ..AnalyzeBlocks
using ..AnalyzeCells
using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..AnalyzeModules
using ..ComputeBlocks
using ..ComputeModules
using ..Contracts
using ..ProjectCells

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_cells_dispersion_per_metacell_per_module
import Metacells.Contracts.matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block
import Metacells.Contracts.matrix_of_is_environment_distinct_per_gene_per_block
import Metacells.Contracts.matrix_of_is_environment_marker_per_gene_per_block
import Metacells.Contracts.matrix_of_is_found_per_module_per_block
import Metacells.Contracts.matrix_of_is_in_environment_per_metacell_per_block
import Metacells.Contracts.matrix_of_is_in_neighborhood_per_block_per_block
import Metacells.Contracts.matrix_of_is_neighborhood_marker_per_gene_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_block
import Metacells.Contracts.matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks
import Metacells.Contracts.matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block
import Metacells.Contracts.matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_module_per_gene_per_block
import Metacells.Contracts.matrix_of_module_status_per_gene_per_block
import Metacells.Contracts.matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_most_correlated_quantile_per_gene_in_neighborhood_per_gene_per_block
import Metacells.Contracts.matrix_of_n_genes_per_module_per_block
import Metacells.Contracts.matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_block
import Metacells.Contracts.module_axis
import Metacells.Contracts.vector_of_anchor_per_module
import Metacells.Contracts.vector_of_block_closest_by_pertinent_markers_per_cell
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_n_cells_per_block
import Metacells.Contracts.vector_of_n_environment_cells_per_block
import Metacells.Contracts.vector_of_n_environment_metacells_per_block
import Metacells.Contracts.vector_of_n_metacells_per_block
import Metacells.Contracts.vector_of_n_modules_per_block
import Metacells.Contracts.vector_of_n_neighborhood_blocks_per_block
import Metacells.Contracts.vector_of_n_neighborhood_cells_per_block
import Metacells.Contracts.vector_of_n_neighborhood_metacells_per_block
import Metacells.Contracts.vector_of_total_environment_UMIs_per_block
import Metacells.Contracts.vector_of_total_neighborhood_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_block
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_type_per_block
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_is_correlated_with_skeleton_per_gene
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_forbidden_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_regulator_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_marker_rank_per_gene
import Metacells.Contracts.vector_of_is_base_outlier_per_cell
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell
import Metacells.Contracts.vector_of_type_per_cell
import Metacells.Contracts.vector_of_type_per_metacell
import Metacells.Contracts.vector_of_umap_x_per_metacell
import Metacells.Contracts.vector_of_umap_y_per_metacell

"""
    prepare_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Aggregate the cells of each metacell, given which cells belong to which metacell. This computes everything about the
metacells which follows from the cells alone, per-metacell properties, and basic gene properties derived from metacells
alone. Gene properties that depend on gene masks (other than exclusion) are not computed here.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_cell(OptionalInput),
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_type_per_metacell(OptionalOutput),
        matrix_of_UMIs_per_gene_per_metacell(CreatedOutput),
        vector_of_total_UMIs_per_metacell(CreatedOutput),
        vector_of_n_cells_per_metacell(CreatedOutput),
        matrix_of_linear_fraction_per_gene_per_metacell(CreatedOutput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(CreatedOutput),
        vector_of_is_marker_per_gene(CreatedOutput),
        vector_of_marker_rank_per_gene(CreatedOutput),
    ],
) function prepare_metacells!(daf::DafWriter; overwrite::Bool = false)::Nothing
    # The types of the metacells come from the types of their cells, so without the one there is not the other.
    if has_vector(daf, "cell", "type")
        compute_vector_of_type_per_metacell_by_cells!(daf; overwrite)
    end
    compute_matrix_of_UMIs_per_gene_per_metacell!(daf; overwrite)
    compute_vector_of_total_UMIs_per_metacell!(daf; overwrite)
    compute_vector_of_n_cells_per_metacell!(daf; overwrite)
    compute_matrix_of_linear_fraction_per_gene_per_metacell!(daf; overwrite)
    compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(daf; overwrite)
    compute_vector_of_is_marker_per_gene!(daf; overwrite)
    compute_vector_of_marker_rank_per_gene!(daf; overwrite)
    return nothing
end

"""
    analyze_metacells!(
        daf::DafWriter;
        prefix::AbstractString = $(DEFAULT.prefix),
        prev_daf::Maybe{DafReader} = nothing,
        module_status::Bool = $(DEFAULT.module_status),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute more advanced properties based on the metacells, taking into account gene masks. Specifically this depends on
the lateral and regulator gene masks, and the forbidden gene masks. These are used to determine the set of skeleton
genes which are then used to drive the rest of the analysis, starting with grouping metacells into blocks and ending
with local gene modules.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [
        gene_axis(RequiredInput),
        cell_axis(RequiredInput),
        metacell_axis(RequiredInput),
        block_axis(GuaranteedOutput),
        module_axis(GuaranteedOutput),
    ],
    data = [
        # What this is given: the cells, the metacells they were aggregated into, and the gene masks.
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_n_cells_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_forbidden_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),

        # Types are optional throughout, so what is computed from them is optional as well.
        vector_of_type_per_metacell(OptionalInput),
        vector_of_type_per_block(OptionalOutput),

        # The skeleton genes, and the geometry of the metacells which follows from them.
        vector_of_is_skeleton_per_gene(CreatedOutput),
        matrix_of_max_skeleton_fold_distance_between_metacells(CreatedOutput),
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(CreatedOutput),
        vector_of_umap_x_per_metacell(CreatedOutput),
        vector_of_umap_y_per_metacell(CreatedOutput),

        # The blocks, and what each is made of.
        vector_of_block_per_metacell(CreatedOutput),
        vector_of_block_closest_by_pertinent_markers_per_cell(CreatedOutput),
        vector_of_n_metacells_per_block(CreatedOutput),
        vector_of_n_cells_per_block(CreatedOutput),
        vector_of_total_UMIs_per_block(CreatedOutput),
        matrix_of_UMIs_per_gene_per_block(CreatedOutput),
        matrix_of_linear_fraction_per_gene_per_block(CreatedOutput),
        matrix_of_log_linear_fraction_per_gene_per_block(CreatedOutput),
        matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block(CreatedOutput),
        matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks(CreatedOutput),
        matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block(CreatedOutput),

        # The neighborhood of each block.
        matrix_of_is_in_neighborhood_per_block_per_block(CreatedOutput),
        vector_of_n_neighborhood_blocks_per_block(CreatedOutput),
        vector_of_n_neighborhood_metacells_per_block(CreatedOutput),
        vector_of_n_neighborhood_cells_per_block(CreatedOutput),
        vector_of_total_neighborhood_UMIs_per_block(CreatedOutput),
        matrix_of_is_neighborhood_marker_per_gene_per_block(CreatedOutput),
        matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block(CreatedOutput),
        matrix_of_most_correlated_quantile_per_gene_in_neighborhood_per_gene_per_block(CreatedOutput),

        # The environment of each block, which the gene modules are estimated over.
        matrix_of_is_in_environment_per_metacell_per_block(CreatedOutput),
        vector_of_n_environment_metacells_per_block(CreatedOutput),
        vector_of_n_environment_cells_per_block(CreatedOutput),
        vector_of_total_environment_UMIs_per_block(CreatedOutput),
        matrix_of_is_environment_marker_per_gene_per_block(CreatedOutput),
        matrix_of_is_environment_distinct_per_gene_per_block(CreatedOutput),
        matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block(CreatedOutput),

        # The gene modules of each block.
        vector_of_anchor_per_module(CreatedOutput),
        vector_of_n_modules_per_block(CreatedOutput),
        matrix_of_module_per_gene_per_block(CreatedOutput),
        matrix_of_is_found_per_module_per_block(CreatedOutput),
        matrix_of_module_status_per_gene_per_block(OptionalOutput),
        matrix_of_n_genes_per_module_per_block(CreatedOutput),
        matrix_of_mean_linear_fraction_in_environment_cells_per_module_per_block(CreatedOutput),
        matrix_of_std_linear_fraction_in_environment_cells_per_module_per_block(CreatedOutput),
        matrix_of_cells_dispersion_per_metacell_per_module(CreatedOutput),
    ],
) function analyze_metacells!(
    daf::DafWriter;
    prefix::AbstractString = "B",
    prev_daf::Maybe{DafReader} = nothing,
    module_status::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    # Which genes predict the rest, and the geometry of the metacells which follows from them.
    compute_vector_of_is_skeleton_per_gene!(daf; overwrite)
    compute_matrix_of_max_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_metacells_2d_umap!(daf; prev_daf, rng, overwrite)

    # The blocks - regions of the manifold the metacells fall into - and what each is made of.
    compute_metacells_blocks!(daf; prefix)
    compute_matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block!(daf; overwrite)
    compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!(daf; overwrite)
    compute_vector_of_n_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_cells_per_block!(daf; overwrite)
    compute_matrix_of_UMIs_per_gene_per_block!(daf; overwrite)
    compute_vector_of_total_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_linear_fraction_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_log_linear_fraction_per_gene_per_block!(daf; overwrite)

    # A block has a type only when the metacells have one, which is optional data.
    if has_vector(daf, "metacell", "type")
        compute_vector_of_type_per_block_by_metacells!(daf; overwrite)
    end

    compute_vector_of_block_closest_by_pertinent_markers_per_cell!(daf; overwrite)
    compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!(daf; overwrite)

    # The neighborhood of a block is the blocks close enough to it to be describing the same local behavior.
    compute_matrix_of_is_in_neighborhood_per_block_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_blocks_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_cells_per_block!(daf; overwrite)
    compute_vector_of_total_neighborhood_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_is_neighborhood_marker_per_gene_per_block!(daf; overwrite)

    # The environment extends the neighborhood with metacells close enough to the block, which gives the gene modules
    # below more metacells to be estimated from.
    compute_matrix_of_is_in_environment_per_metacell_per_block!(daf; overwrite)
    compute_vector_of_n_environment_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_environment_cells_per_block!(daf; overwrite)
    compute_vector_of_total_environment_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_is_environment_marker_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_is_environment_distinct_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block!(daf; overwrite)

    # Always recomputed rather than inherited: this describes the metacells of this round, and a repository chained onto
    # an earlier one would otherwise keep that round's answer.
    compute_matrix_of_most_correlated_gene_in_neighborhood_per_gene_per_block!(daf; overwrite = true)

    # The gene modules of each block - the local programs the sharpening clusters cells by.
    compute_blocks_modules!(daf; module_status, rng, overwrite)
    compute_vector_of_n_modules_per_block!(daf; overwrite)
    compute_matrix_of_n_genes_per_module_per_block!(daf; overwrite)
    compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!(daf; overwrite)
    compute_matrix_of_cells_dispersion_per_metacell_per_module!(daf; overwrite)

    return nothing
end

"""
    function import_base_metacells!(;
        cells_daf::DafWriter,
        metacells_daf::DafWriter,
        metacell_per_cell::AbstractVector{<:AbstractString},
        empty_metacells::Maybe{EmptyImplicit} = nothing,
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Bring in the base metacells the sharpening pipeline starts with. The cells are expected to have been imported already,
but their assignment to metacells is supplied explicitly. It is filtered through `empty_metacells` to allow for values
like `Outliers` to be safely converted to the empty string before being applied.

# Cells

$(CONTRACT1)

# Metacells

$(CONTRACT2)
"""
@logged :mcs_ops @computation Contract(;
    name = "cells_daf",
    axes = [cell_axis(RequiredInput)],
    data = [vector_of_is_base_outlier_per_cell(GuaranteedOutput)],
) Contract(;
    name = "metacells_daf",
    axes = [cell_axis(RequiredInput)],
    data = [vector_of_metacell_per_cell(GuaranteedOutput)],
) function import_base_metacells!(;  # UNTESTED
    cells_daf::DafWriter,
    metacells_daf::DafWriter,
    metacell_per_cell::AbstractVector{<:AbstractString},
    empty_metacells::Maybe{EmptyImplicit} = nothing,
    overwrite::Bool = false,
)::Nothing
    set_vector!(metacells_daf, "cell", "metacell", metacell_per_cell; overwrite)
    if empty_metacells !== nothing
        unify_empty_vector_values!(metacells_daf; axis = "cell", property = "metacell", empty_values = empty_metacells)
    end
    compute_vector_of_is_base_outlier_per_cell!(; cells_daf, metacells_daf, overwrite)
    return nothing
end

end  # module
