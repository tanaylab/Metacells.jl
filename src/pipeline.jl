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
export qc_metacells!

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
import Metacells.Contracts.cell_axis
import Metacells.Contracts.vector_of_metacell_per_cell

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
@logged :mcs_ops @computation (
    optional_contract(function_contract(compute_vector_of_type_per_metacell_by_cells!)) |>
    function_contract(compute_matrix_of_UMIs_per_gene_per_metacell!) |>
    function_contract(compute_vector_of_total_UMIs_per_metacell!) |>
    function_contract(compute_vector_of_n_cells_per_metacell!) |>
    function_contract(compute_matrix_of_linear_fraction_per_gene_per_metacell!) |>
    function_contract(compute_matrix_of_log_linear_fraction_per_gene_per_metacell!) |>
    function_contract(compute_vector_of_is_marker_per_gene!) |>
    function_contract(compute_vector_of_marker_rank_per_gene!) |>
    function_contract(compute_matrix_of_correlation_between_markers_per_gene_per_gene!)
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
    compute_matrix_of_correlation_between_markers_per_gene_per_gene!(daf; overwrite)
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

# Metacells

$(CONTRACT1)

# Previous Metacells

Only if `prev_daf` is given.

$(CONTRACT2)
"""
@logged :mcs_ops @computation (
    function_contract(compute_vector_of_is_skeleton_per_gene!) |>
    function_contract(compute_matrix_of_max_skeleton_fold_distance_between_metacells!) |>
    function_contract(compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!) |>
    function_contract(compute_metacells_2d_umap!) |>
    function_contract(compute_metacells_blocks!) |>
    function_contract(compute_matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block!) |>
    function_contract(compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!) |>
    function_contract(compute_vector_of_n_metacells_per_block!) |>
    function_contract(compute_vector_of_n_cells_per_block!) |>
    function_contract(compute_matrix_of_UMIs_per_gene_per_block!) |>
    function_contract(compute_vector_of_total_UMIs_per_block!) |>
    function_contract(compute_matrix_of_linear_fraction_per_gene_per_block!) |>
    function_contract(compute_matrix_of_log_linear_fraction_per_gene_per_block!) |>
    optional_contract(function_contract(compute_vector_of_type_per_block_by_metacells!)) |>
    function_contract(compute_vector_of_block_closest_by_pertinent_markers_per_cell!) |>
    function_contract(compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!) |>
    function_contract(compute_matrix_of_is_in_neighborhood_per_block_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_blocks_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_metacells_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_cells_per_block!) |>
    function_contract(compute_vector_of_total_neighborhood_UMIs_per_block!) |>
    function_contract(compute_matrix_of_is_neighborhood_marker_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_in_environment_per_metacell_per_block!) |>
    function_contract(compute_vector_of_n_environment_metacells_per_block!) |>
    function_contract(compute_vector_of_n_environment_cells_per_block!) |>
    function_contract(compute_vector_of_total_environment_UMIs_per_block!) |>
    function_contract(compute_matrix_of_is_environment_marker_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_environment_distinct_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block!) |>
    function_contract(compute_blocks_modules!) |>
    function_contract(compute_vector_of_n_modules_per_block!) |>
    function_contract(compute_matrix_of_n_genes_per_module_per_block!) |>
    function_contract(compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!) |>
    function_contract(compute_matrix_of_cells_dispersion_per_metacell_per_module!)
) function_contract(compute_metacells_2d_umap!, 2) function analyze_metacells!(
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
    compute_metacells_blocks!(daf; prefix, overwrite)
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
@logged :mcs_ops @computation function_contract(compute_vector_of_is_base_outlier_per_cell!, 1) renamed_contract(
    # The metacell of each cell is set directly rather than computed, so it has no contract to combine.
    Contract(;
        name = "metacells_daf",
        axes = [cell_axis(RequiredInput)],
        data = [vector_of_metacell_per_cell(GuaranteedOutput)],
    ) |> function_contract(compute_vector_of_is_base_outlier_per_cell!, 2),
    "metacells_daf",
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

"""
    function qc_metacells!(;
        daf::DafWriter,
        base_daf::DafReader,
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Say how well the metacells describe the cells they were aggregated from, judged in the locations of the manifold a base
repository laid out.

Each cell is correlated against its own metacell minus itself, over the cells of each base block's neighborhood, and
the result is averaged over the environment marker genes of that base block. A higher number means the metacells of
that location are a better account of the cells there. Only metacells scored against the same `base_daf` can be
compared to each other, since it is the `base_daf` which decides both the locations and the genes each of them is read
by. A repository is free to be its own base, which is how the metacells the sharpening started from are scored.

# Daf

$(CONTRACT1)

# Base

$(CONTRACT2)
"""
@logged :mcs_ops @computation renamed_contract(
    function_contract(
        compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!,
        1,
    ) |> function_contract(
        compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!,
        1,
    ),
    "daf",
) renamed_contract(
    function_contract(
        compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!,
        2,
    ) |> function_contract(
        compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!,
        2,
    ),
    "base_daf",
) function qc_metacells!(; daf::DafWriter, base_daf::DafReader, overwrite::Bool = false)::Nothing
    compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!(;
        other_daf = daf,
        base_daf,
        overwrite,
    )
    compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!(;
        other_daf = daf,
        base_daf,
        overwrite,
    )
    return nothing
end

end  # module
