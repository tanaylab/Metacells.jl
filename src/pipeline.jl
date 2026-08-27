"""
Run whole stages of the metacells pipeline.

The functions here compute nothing themselves. Each runs a sequence of the computations of the rest of the package, in
the order their inputs and outputs demand, which is otherwise something every caller would have to know and repeat.
"""
module Pipeline

export prepare_markers!
export prepare_metacells!
export prepare_skeletons!

using DataAxesFormats
using Random
using TanayLabUtilities

using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_correlation_between_markers_per_gene_per_gene
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
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
metacells which follows from the cells alone - their UMIs, their sizes, and the fraction of the UMIs of each gene in
each of them.

Types are optional. If the cells have a type, then so does each metacell, by the types of its cells; if they do not,
then neither do the metacells, and everything else here is computed just the same.

Nothing here depends on the gene masks, so a repository this was run on can be shared by several analyses which use
different masks. Do run it again whenever the assignment of cells to metacells changes, that is, after each round of
[`sharpen_metacells!`](@ref Metacells.SharpenMetacells.sharpen_metacells!).

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
    return nothing
end

"""
    prepare_markers!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Find the marker genes - the genes which distinguish between the metacells - rank them, and correlate each of them with
every other.

Like [`prepare_metacells!`](@ref), and for the same reason, nothing here depends on the gene masks. It does depend on
the metacells, so run it again whenever they change.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        vector_of_is_marker_per_gene(CreatedOutput),
        vector_of_marker_rank_per_gene(CreatedOutput),
        matrix_of_correlation_between_markers_per_gene_per_gene(CreatedOutput),
    ],
) function prepare_markers!(daf::DafWriter; overwrite::Bool = false)::Nothing
    compute_vector_of_is_marker_per_gene!(daf; overwrite)
    compute_vector_of_marker_rank_per_gene!(daf; overwrite)
    compute_matrix_of_correlation_between_markers_per_gene_per_gene!(daf; overwrite)
    return nothing
end

"""
    prepare_skeletons!(
        daf::DafWriter;
        prev_daf::Maybe{DafReader} = nothing,
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Choose the skeleton genes - the genes whose expression is taken to predict the rest - and compute what follows from
them: which marker genes correlate with them, how far the metacells are from each other, and a 2D layout of the
metacells based on that distance.

This is the part which the gene masks decide. Two analyses of the same metacells which differ only in their
`is_lateral`, `is_forbidden` or `is_regulator` masks differ from here onwards, and only from here onwards, so the
results of [`prepare_metacells!`](@ref) and [`prepare_markers!`](@ref) can be shared between them.

Give `prev_daf` to seed the layout from an earlier one, so that a metacell of the same cells is placed close to where
it was, and successive rounds of [`sharpen_metacells!`](@ref Metacells.SharpenMetacells.sharpen_metacells!) can be
compared by eye.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput), cell_axis(OptionalInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_is_forbidden_per_gene(RequiredInput),
        vector_of_is_regulator_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        matrix_of_correlation_between_markers_per_gene_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_metacell_per_cell(OptionalInput),
        vector_of_is_skeleton_per_gene(CreatedOutput),
        vector_of_is_correlated_with_skeleton_per_gene(CreatedOutput),
        matrix_of_max_skeleton_fold_distance_between_metacells(CreatedOutput),
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(CreatedOutput),
        vector_of_umap_x_per_metacell(CreatedOutput),
        vector_of_umap_y_per_metacell(CreatedOutput),
    ],
) function prepare_skeletons!(
    daf::DafWriter;
    prev_daf::Maybe{DafReader} = nothing,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    compute_vector_of_is_skeleton_per_gene!(daf; overwrite)
    compute_vector_of_is_correlated_with_skeleton_per_gene!(daf; overwrite)
    compute_matrix_of_max_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_metacells_2d_umap!(daf; prev_daf, rng, overwrite)
    return nothing
end

end  # module
