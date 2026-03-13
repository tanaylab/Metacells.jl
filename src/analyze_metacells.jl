"""
Do simple metacells analysis.
"""
module AnalyzeMetacells

export compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!
export compute_matrix_of_linear_fraction_per_gene_per_metacell!
export compute_matrix_of_log_linear_fraction_per_gene_per_metacell!
export compute_matrix_of_max_skeleton_fold_distance_between_metacells!
export compute_matrix_of_UMIs_per_gene_per_metacell!
export compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!
export compute_vector_of_n_cells_per_metacell!
export compute_vector_of_total_UMIs_per_metacell!
export compute_vector_of_type_per_cell_by_metacells!
export compute_vector_of_type_per_metacell_by_cells!

using Base.Threads
using DataAxesFormats
using Distances
using Distributions
using MultipleTesting
using TanayLabUtilities
using Random
using StatsBase

using ..Defaults
using ..Contracts

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_euclidean_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_log_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_max_skeleton_fold_distance_between_metacells
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_correlation_between_cells_and_punctuated_metacells_per_gene
import Metacells.Contracts.vector_of_is_excluded_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_is_lateral_per_gene
import Metacells.Contracts.vector_of_is_marker_per_gene
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_metacell_per_cell
import Metacells.Contracts.vector_of_n_cells_per_metacell
import Metacells.Contracts.vector_of_total_UMIs_per_cell
import Metacells.Contracts.vector_of_total_UMIs_per_metacell
import Metacells.Contracts.vector_of_type_per_cell
import Metacells.Contracts.vector_of_type_per_metacell

"""
    function compute_vector_of_n_cells_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_n_cells_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [vector_of_metacell_per_cell(RequiredInput), vector_of_n_cells_per_metacell(CreatedOutput)],
) function compute_vector_of_n_cells_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_cells_per_metacell = daf["@ cell / metacell =@ >> Count"].array
    set_vector!(daf, "metacell", "n_cells", n_cells_per_metacell; overwrite)
    @debug "Mean cells in metacell: $(mean(n_cells_per_metacell))" _group = :mcs_results
    return nothing
end

"""
    function compute_matrix_of_UMIs_per_gene_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_UMIs_per_gene_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_UMIs_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    UMIs_per_gene_per_metacell = daf["@ gene @ cell :: UMIs |/ metacell =@ >| Sum"].array
    set_matrix!(daf, "gene", "metacell", "UMIs", bestify(UMIs_per_gene_per_metacell); overwrite)  # NOJET
    return nothing
end

"""
    function compute_vector_of_total_UMIs_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_total_UMIs_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_total_UMIs_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_metacell = daf["@ metacell @ gene [ ! is_excluded ] :: UMIs >| Sum"].array
    set_vector!(daf, "metacell", "total_UMIs", total_UMIs_per_metacell; overwrite)
    @debug "Mean UMIs in metacell: $(mean(total_UMIs_per_metacell))" _group = :mcs_results
    return nothing
end

"""
    function compute_matrix_of_linear_fraction_per_gene_per_metacell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_linear_fraction_per_gene_per_metacell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_linear_fraction_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")

    indices_of_included_genes = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    UMIs_per_metacell_per_included_gene = daf["@ metacell @ gene [ ! is_excluded ] :: UMIs"].array
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    linear_fraction_per_metacell_per_gene = zeros(Float32, n_metacells, n_genes)
    @. linear_fraction_per_metacell_per_gene[:, indices_of_included_genes] =
        UMIs_per_metacell_per_included_gene / total_UMIs_per_metacell

    set_matrix!(  # NOJET
        daf,
        "metacell",
        "gene",
        "linear_fraction",
        bestify(linear_fraction_per_metacell_per_gene);
        overwrite,
    )  # NOJET

    return nothing
end

"""
    function compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_log_linear_fraction_per_gene_per_metacell`](@ref). This adds the
`gene_fraction_regularization` to deal with zero fractions.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [gene_axis(RequiredInput), metacell_axis(RequiredInput)],
    data = [
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(CreatedOutput),
    ],
) function compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "linear_fraction").array
    empty_dense_matrix!(
        daf,
        "metacell",
        "gene",
        "log_linear_fraction",
        Float32;
        overwrite,
    ) do log_fraction_per_metacell_per_gene
        @. log_fraction_per_metacell_per_gene = log2(fraction_per_metacell_per_gene + gene_fraction_regularization)
        return nothing
    end
    return nothing
end

"""
    function compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_euclidean_skeleton_fold_distance_between_metacells`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_log_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_euclidean_skeleton_fold_distance_between_metacells(CreatedOutput),
    ],
) function compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    log_linear_fraction_per_metacell_per_skeleton =
        daf["@ metacell @ gene [ is_skeleton ] :: log_linear_fraction"].array
    log_linear_fraction_per_skeleton_per_metacell = flipped(log_linear_fraction_per_metacell_per_skeleton)
    n_metacells = axis_length(daf, "metacell")
    distances_between_metacells = parallel_pairwise(  # NOJET
        Euclidean(),
        log_linear_fraction_per_skeleton_per_metacell;
        dims = 2,
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "euclidean_skeleton_fold_distance_between_metacells",
        ),
    )
    set_matrix!(daf, "metacell", "metacell", "euclidean_skeleton_fold_distance", distances_between_metacells; overwrite)
    return nothing
end

"""
    function compute_matrix_of_max_skeleton_fold_distance_between_metacells!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = $(DEFAULT.gene_fraction_regularization),
        min_significant_gene_UMIs::Integer = $(DEFAULT.min_significant_gene_UMIs),
        fold_confidence::AbstractFloat = $(DEFAULT.fold_confidence),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`matrix_of_max_skeleton_fold_distance_between_metacells`](@ref). For computing this fold factor, we add
the `gene_fraction_regularization` to deal with zero fractions. We ignore "insignificant" genes whose total UMIs in the
compared metacells isn't at least `min_significant_gene_UMIs`. We also we reduce the distance using the
`fold_confidence` based on the number of UMIs used to estimate the expression in the metacells. That is, when we have a
high fold factor, we are pretty certain the metacells truly have a significant difference in the expression of the
skeleton gene.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        matrix_of_max_skeleton_fold_distance_between_metacells(CreatedOutput),
    ],
) function compute_matrix_of_max_skeleton_fold_distance_between_metacells!(
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    min_significant_gene_UMIs::Integer = 40,
    fold_confidence::AbstractFloat = 0.9,
    overwrite::Bool = false,
)::Nothing
    @assert gene_fraction_regularization >= 0
    @assert min_significant_gene_UMIs >= 0
    @assert 0 < fold_confidence < 1

    linear_fraction_per_metacell_per_skeleton = daf["@ metacell @ gene [ is_skeleton ] :: linear_fraction"].array
    linear_fraction_per_skeleton_per_metacell = flipped(linear_fraction_per_metacell_per_skeleton)
    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array

    confidence_stds = quantile(Normal(), fold_confidence)

    confidence_linear_fractions_per_skeleton_per_metacells =  # NOJET
        confidence_stds .* sqrt.(linear_fraction_per_skeleton_per_metacell .* transpose(total_UMIs_per_metacell)) ./
        transpose(total_UMIs_per_metacell)

    low_log_linear_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            max.(
                linear_fraction_per_skeleton_per_metacell .- confidence_linear_fractions_per_skeleton_per_metacells,
                0.0,
            ) .+ gene_fraction_regularization,
        )

    high_log_linear_fraction_per_skeleton_per_metacell = # NOJET
        log2.(
            linear_fraction_per_skeleton_per_metacell .+ confidence_linear_fractions_per_skeleton_per_metacells .+
            gene_fraction_regularization,
        )

    UMIs_per_metacell_per_skeleton = daf["@ metacell @ gene [ is_skeleton ] :: UMIs"].array
    UMIs_per_skeleton_per_metacell = flip(UMIs_per_metacell_per_skeleton)

    n_skeletons, n_metacells = size(UMIs_per_skeleton_per_metacell)

    distances_between_metacells = Matrix{Float32}(undef, n_metacells, n_metacells)
    distances_between_metacells[1, 1] = 0.0

    parallel_loop_wo_rng(
        reverse(2:(n_metacells));
        progress = DebugProgress(
            n_metacells - 1;
            group = :mcs_loops,
            desc = "max_skeleton_fold_distance_between_metacells",
        ),
    ) do base_metacell_index
        distances_between_metacells[base_metacell_index, base_metacell_index] = 0.0

        @views base_metacell_UMIs_per_skeleton = vec(UMIs_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_low_log_linear_fraction_per_skeleton =
            vec(low_log_linear_fraction_per_skeleton_per_metacell[:, base_metacell_index])
        @views base_metacell_high_log_linear_fraction_per_skeleton =
            vec(high_log_linear_fraction_per_skeleton_per_metacell[:, base_metacell_index])

        n_other_metacells = base_metacell_index - 1
        index_per_other_metacell = 1:(base_metacell_index - 1)
        @views UMIs_per_skeleton_per_other_metacells = UMIs_per_skeleton_per_metacell[:, index_per_other_metacell]
        @views low_log_linear_fraction_per_skeleton_per_other_metacells =
            low_log_linear_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]
        @views high_log_linear_fraction_per_skeleton_per_other_metacells =
            high_log_linear_fraction_per_skeleton_per_metacell[:, index_per_other_metacell]

        significant_fold_per_skeleton_per_other_metacell = confident_gene_distance.(
            min_significant_gene_UMIs,
            base_metacell_UMIs_per_skeleton,
            base_metacell_low_log_linear_fraction_per_skeleton,
            base_metacell_high_log_linear_fraction_per_skeleton,
            UMIs_per_skeleton_per_other_metacells,
            low_log_linear_fraction_per_skeleton_per_other_metacells,
            high_log_linear_fraction_per_skeleton_per_other_metacells,
        )
        @assert_matrix(significant_fold_per_skeleton_per_other_metacell, n_skeletons, n_other_metacells, Columns)

        distances_between_base_and_other_metacells =
            vec(maximum(significant_fold_per_skeleton_per_other_metacell; dims = 1))
        @assert_vector(distances_between_base_and_other_metacells, n_other_metacells)

        distances_between_metacells[index_per_other_metacell, base_metacell_index] .=
            distances_between_base_and_other_metacells
        distances_between_metacells[base_metacell_index, index_per_other_metacell] .=
            distances_between_base_and_other_metacells
        return nothing
    end

    set_matrix!(daf, "metacell", "metacell", "max_skeleton_fold_distance", distances_between_metacells; overwrite)
    return nothing
end

@inline function confident_gene_distance(  # UNTESTED
    min_significant_gene_UMIs::Integer,
    base_metacell_total_UMIs_of_gene::Integer,
    base_metacell_low_log_linear_fraction_of_gene::AbstractFloat,
    base_metacell_high_log_linear_fraction_of_gene::AbstractFloat,
    other_metacell_total_UMIs_of_gene::Integer,
    other_metacell_low_log_linear_fraction_of_gene::AbstractFloat,
    other_metacell_high_log_linear_fraction_of_gene::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = base_metacell_total_UMIs_of_gene + other_metacell_total_UMIs_of_gene
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = base_metacell_high_log_linear_fraction_of_gene < other_metacell_high_log_linear_fraction_of_gene

    highest_low_log_linear_fraction_of_gene =
        is_base_low * base_metacell_high_log_linear_fraction_of_gene +
        !is_base_low * other_metacell_high_log_linear_fraction_of_gene

    lowest_high_log_linear_fraction_of_gene =
        is_base_low * other_metacell_low_log_linear_fraction_of_gene +
        !is_base_low * base_metacell_low_log_linear_fraction_of_gene

    confident_gap = lowest_high_log_linear_fraction_of_gene - highest_low_log_linear_fraction_of_gene

    return (is_significant * max(confident_gap, 0.0))
end

"""
    function compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!(
        daf::DafWriter;
        gene_fraction_regularization::AbstractFloat = (DEFAULT.gene_fraction_regularization),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_correlation_between_cells_and_punctuated_metacells_per_gene`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        vector_of_is_marker_per_gene(RequiredInput),
        vector_of_is_lateral_per_gene(RequiredInput),
        vector_of_metacell_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_correlation_between_cells_and_punctuated_metacells_per_gene(CreatedOutput),
    ],
) function compute_vector_of_correlation_between_cells_and_punctuated_metacells_per_gene!(  # UNTESTED
    daf::DafWriter;
    gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_CELLS,
    overwrite::Bool = false,
)::Nothing
    n_genes = axis_length(daf, "gene")

    index_per_included_gene = get_query(daf, "@ gene [ ! is_excluded ] : index").array
    n_included_genes = length(index_per_included_gene)

    indices_of_grouped_cells = get_query(daf, "@ cell [ metacell ] : index").array
    n_grouped_cells = length(indices_of_grouped_cells)

    UMIs_per_cell_per_gene = get_matrix(daf, "cell", "gene", "UMIs").array
    UMIs_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "UMIs").array

    total_UMIs_per_metacell = get_vector(daf, "metacell", "total_UMIs").array
    total_UMIs_per_cell = get_vector(daf, "cell", "total_UMIs").array
    total_UMIs_per_grouped_cell = total_UMIs_per_cell[indices_of_grouped_cells]

    metacell_index_per_grouped_cell = daf["@ cell : metacell ?? : index"].array

    total_punctuated_metacell_UMIs_per_grouped_cell =
        total_UMIs_per_metacell[metacell_index_per_grouped_cell] .- total_UMIs_per_grouped_cell

    correlation_between_cells_and_punctuated_metacells_per_gene = Vector{Float32}(undef, n_genes)

    gene_cell_log_fraction_per_grouped_cell_per_thread = [Vector{Float32}(undef, n_grouped_cells) for _ in 1:nthreads()]
    gene_punctuated_metacell_log_fraction_per_grouped_cell_per_thread =
        [Vector{Float32}(undef, n_grouped_cells) for _ in 1:nthreads()]

    parallel_loop_wo_rng(
        1:n_included_genes;
        policy = :static,
        progress = DebugProgress(
            n_included_genes;
            group = :mcs_loops,
            desc = "correlation_between_cells_and_punctuated_metacells_per_gene",
        ),
        progress_chunk = 100,
    ) do included_gene_position
        gene_index = index_per_included_gene[included_gene_position]

        gene_cell_log_fraction_per_grouped_cell = gene_cell_log_fraction_per_grouped_cell_per_thread[threadid()]
        gene_punctuated_metacell_log_fraction_per_grouped_cell =
            gene_punctuated_metacell_log_fraction_per_grouped_cell_per_thread[threadid()]

        for (grouped_cell_position, cell_index) in enumerate(indices_of_grouped_cells)
            gene_cell_UMIs = UMIs_per_cell_per_gene[cell_index, gene_index]
            gene_cell_log_fraction_per_grouped_cell[grouped_cell_position] =
                log2(gene_cell_UMIs / total_UMIs_per_grouped_cell[grouped_cell_position] + gene_fraction_regularization)
            metacell_index = metacell_index_per_grouped_cell[grouped_cell_position]
            gene_metacell_UMIs = UMIs_per_metacell_per_gene[metacell_index, gene_index]
            gene_punctuated_metacell_log_fraction_per_grouped_cell[grouped_cell_position] = log2(
                (
                    (gene_metacell_UMIs .- gene_cell_UMIs) /
                    total_punctuated_metacell_UMIs_per_grouped_cell[grouped_cell_position]
                ) .+ gene_fraction_regularization,
            )
        end

        correlation_between_cells_and_punctuated_metacells_per_gene[gene_index] = zero_cor_between_vectors(
            gene_cell_log_fraction_per_grouped_cell,
            gene_punctuated_metacell_log_fraction_per_grouped_cell,
        )

        return nothing
    end

    set_vector!(
        daf,
        "gene",
        "correlation_between_cells_and_punctuated_metacells",
        correlation_between_cells_and_punctuated_metacells_per_gene;
        overwrite,
    )

    is_pertinent_marker_per_gene =
        get_vector(daf, "gene", "is_marker").array .& .!get_vector(daf, "gene", "is_lateral").array
    @debug (
        "Mean correlation of pertinent marker genes between cells and their punctuated metacells: " *
        "$(mean(correlation_between_cells_and_punctuated_metacells_per_gene[is_pertinent_marker_per_gene]))"
    ) _group = :mcs_results

    return nothing
end

"""
    function compute_vector_of_type_per_cell_by_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_cell`](@ref) by using the type of the metacell each cell belongs to. Cells that do
not belong to any metacells (including excluded cells) are given an empty type. This method is used when metacells are
annotated with types (which is typically a supervised process).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_metacell(RequiredInput),
        vector_of_type_per_cell(CreatedOutput),
    ],
) function compute_vector_of_type_per_cell_by_metacells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_cell = daf["@ cell : metacell ?? '' : type"].array
    set_vector!(daf, "cell", "type", type_per_cell; overwrite)
    return nothing
end

"""
    function compute_vector_of_type_per_metacell_by_cells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_type_per_metacell`](@ref) by using the type of the cells grouped into each metacell. The
most frequent ("mode") of the type is used for each metacell. This method is used when cells are annotated with types
(typically when importing single-cell data with type annotations).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), cell_axis(RequiredInput)],
    data = [
        vector_of_metacell_per_cell(RequiredInput),
        vector_of_type_per_cell(RequiredInput),
        vector_of_type_per_metacell(CreatedOutput),
    ],
) function compute_vector_of_type_per_metacell_by_cells!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    type_per_metacell = daf["@ cell : type / metacell =@ >> Mode"].array
    set_vector!(daf, "metacell", "type", type_per_metacell; overwrite)
    return nothing
end

end  # module
