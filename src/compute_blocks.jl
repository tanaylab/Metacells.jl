"""
Group "very close" metacells into disjoint blocks.
"""
module ComputeBlocks

export compute_metacells_blocks!

using Base.Threads
using Clustering
using DataAxesFormats
using Distributions
using LinearAlgebra
using LoopVectorization
using MultivariateStats
using Random
using StatsBase
using TanayLabUtilities

using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..Contracts
using ..Defaults

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.block_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_linear_fraction_per_gene_per_metacell
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_metacell
import Metacells.Contracts.metacell_axis
import Metacells.Contracts.vector_of_block_per_metacell
import Metacells.Contracts.vector_of_is_skeleton_per_gene
import Metacells.Contracts.vector_of_total_UMIs_per_metacell

# The low and high log of the confidence interval of the linear fraction of each skeleton gene in each metacell.
@kwdef struct ConfidenceIntervals
    low_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
    high_log_linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat}
end

"""
    function compute_metacells_blocks!(
        daf::DafWriter;
        n_blocks::Integer = $(DEFAULT.n_blocks),
        max_block_span::Real = $(DEFAULT.max_block_span),
        max_fold_gene_fraction_regularization::AbstractFloat = $(DEFAULT.max_fold_gene_fraction_regularization),
        max_fold_min_significant_gene_UMIs::Integer = $(DEFAULT.max_fold_min_significant_gene_UMIs),
        max_fold_confidence::AbstractFloat = $(DEFAULT.max_fold_confidence),
        max_fold_min_confidence_UMIs::Real = $(DEFAULT.max_fold_min_confidence_UMIs),
        prefix::AbstractString = $(DEFAULT.prefix),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set the [`block_axis`](@ref) and [`vector_of_block_per_metacell`](@ref). We group the metacells into blocks
where the metacells in each one are "similar" in all the skeleton genes. That is, each block is an approximation of a
single cell state, assuming the skeleton genes adequately predict the rest of the genes.

If `n_blocks` is not specified, we decide on the number of blocks by using hierarchical clustering (`hclust` with
`:complete` linkage) based on the maximal skeleton fold distance between the metacells, such that the distance between
the metacells in each cluster is at most `max_block_span`. For computing this fold factor, we add the
`max_fold_gene_fraction_regularization` to deal with zero fractions. We ignore "insignificant" genes whose total UMIs in
the compared metacells isn't at least `max_fold_min_significant_gene_UMIs`. We also reduce the distance using the
`max_fold_confidence` based on the number of UMIs used to estimate the expression in the metacells (but at least
`max_fold_min_confidence_UMIs`). That is, when we have a high fold factor, we are pretty certain the metacells truly
have a significant difference in the expression of the skeleton gene. This distances matrix is only computed for this
purpose and is not stored.

The actual blocks are computed using K-means on the log of the skeleton genes fractions of the metacells, using the
[`GENE_FRACTION_REGULARIZATION_FOR_CELLS`](@ref) regularization. If `n_blocks` is not specified, K-means starts from the
centroids of the hierarchical clustering blocks. Otherwise, it starts from `:kmpp` seeds chosen using the `rng`.
Therefore, specifying `n_blocks` gives different blocks than the automatic method, even if the number of blocks is the
same. Blocks which end up empty are dropped.

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [metacell_axis(RequiredInput), gene_axis(RequiredInput), block_axis(CreatedOutput)],
    data = [
        vector_of_is_skeleton_per_gene(RequiredInput),
        matrix_of_linear_fraction_per_gene_per_metacell(RequiredInput),
        matrix_of_UMIs_per_gene_per_metacell(RequiredInput),
        vector_of_total_UMIs_per_metacell(RequiredInput),
        vector_of_block_per_metacell(CreatedOutput),
    ],
) function compute_metacells_blocks!(  # UNTESTED
    daf::DafWriter;
    n_blocks::Maybe{Integer} = nothing,
    max_block_span::Real = function_default(compute_vector_of_is_marker_per_gene!, :min_marker_gene_range_fold),
    max_fold_gene_fraction_regularization::AbstractFloat = GENE_FRACTION_REGULARIZATION_FOR_METACELLS,
    max_fold_min_significant_gene_UMIs::Integer = MIN_SIGNIFICANT_GENE_UMIS,
    max_fold_confidence::AbstractFloat = 0.9,
    max_fold_min_confidence_UMIs::Real = 3,
    prefix::AbstractString = "B",
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    @assert max_fold_gene_fraction_regularization >= 0
    @assert max_fold_min_significant_gene_UMIs >= 0
    @assert 0 < max_fold_confidence < 1

    linear_fraction_per_metacell_per_skeleton = plain_array(daf["@ metacell @ gene [ is_skeleton ] :: linear_fraction"])
    log_linear_fraction_per_skeleton_per_metacell = Matrix{Float64}(flipped(linear_fraction_per_metacell_per_skeleton))
    log_linear_fraction_per_skeleton_per_metacell .=
        log2.(log_linear_fraction_per_skeleton_per_metacell .+ GENE_FRACTION_REGULARIZATION_FOR_CELLS)

    if n_blocks === nothing
        distances_between_metacells = compute_max_skeleton_fold_distance_between_metacells(
            daf;
            gene_fraction_regularization = max_fold_gene_fraction_regularization,
            min_significant_gene_UMIs = max_fold_min_significant_gene_UMIs,
            fold_confidence = max_fold_confidence,
            min_confidence_UMIs = max_fold_min_confidence_UMIs,
        )
        clusters = hclust(distances_between_metacells; linkage = :complete)
        seed_block_index_per_metacell = cutree(clusters; h = max_block_span)
        mean_log_linear_fraction_per_skeleton_per_block = reduce(
            hcat,
            [
                vec(mean(view(log_linear_fraction_per_skeleton_per_metacell, :, metacell_indices); dims = 2)) for  # NOLINT
                metacell_indices in collect_group_members(seed_block_index_per_metacell)
            ],
        )
    else
        seed_metacell_index_per_block = initseeds(:kmpp, log_linear_fraction_per_skeleton_per_metacell, n_blocks; rng)
        mean_log_linear_fraction_per_skeleton_per_block =
            log_linear_fraction_per_skeleton_per_metacell[:, seed_metacell_index_per_block]
    end

    block_index_per_metacell = assignments(
        kmeans!(log_linear_fraction_per_skeleton_per_metacell, mean_log_linear_fraction_per_skeleton_per_block; rng),
    )
    n_blocks = compact_groups!(block_index_per_metacell)
    metacell_indices_per_block = collect_group_members(block_index_per_metacell)
    name_per_block = group_names(axis_vector(daf, "metacell"), metacell_indices_per_block; prefix)

    add_axis!(daf, "block", name_per_block; overwrite)
    set_vector!(daf, "metacell", "block", name_per_block[block_index_per_metacell]; overwrite)

    @debug "Blocks: $(n_blocks)" _group = :mcs_results

    return nothing
end

# The maximal confident fold distance of any significant skeleton gene between each pair of metacells. See
# `compute_metacells_blocks!` for the meaning of the parameters.
function compute_max_skeleton_fold_distance_between_metacells(  # UNTESTED
    daf::DafReader;
    gene_fraction_regularization::AbstractFloat,
    min_significant_gene_UMIs::Integer,
    fold_confidence::AbstractFloat,
    min_confidence_UMIs::Real,
)::Matrix{Float32}
    linear_fraction_per_metacell_per_skeleton =  # NOJET
        plain_array(daf["@ metacell @ gene [ is_skeleton ] :: linear_fraction"])
    linear_fraction_per_skeleton_per_metacell = Matrix(flip(plain_array(linear_fraction_per_metacell_per_skeleton)))
    total_UMIs_per_metacell = plain_array(densify(get_vector(daf, "metacell", "total_UMIs")))
    confidence_intervals = compute_confidence_intervals(;
        linear_fraction_per_skeleton_per_metacell,
        total_UMIs_per_metacell,
        gene_fraction_regularization,
        fold_confidence,
        min_confidence_UMIs,
    )
    low_log_linear_fraction_per_skeleton_per_metacell =
        confidence_intervals.low_log_linear_fraction_per_skeleton_per_metacell
    high_log_linear_fraction_per_skeleton_per_metacell =
        confidence_intervals.high_log_linear_fraction_per_skeleton_per_metacell

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

    return distances_between_metacells
end

# The low and high log of the confidence interval of the linear fraction of each skeleton gene in each metacell. See
# `compute_metacells_blocks!` for the meaning of the parameters.
function compute_confidence_intervals(;  # UNTESTED
    linear_fraction_per_skeleton_per_metacell::AbstractMatrix{<:AbstractFloat},
    total_UMIs_per_metacell::AbstractVector{<:Real},
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
    min_confidence_UMIs::Real,
)::ConfidenceIntervals
    confidence_stds = quantile(Normal(), fold_confidence)

    confidence_linear_fractions_per_skeleton_per_metacells =
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell))
    @check_turbo_matrix(confidence_linear_fractions_per_skeleton_per_metacells)
    @check_turbo_matrix(linear_fraction_per_skeleton_per_metacell)
    @check_turbo_vector(total_UMIs_per_metacell)
    n_skeletons, n_metacells = size(linear_fraction_per_skeleton_per_metacell)
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "confidence_linear_fractions_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "confidence_linear_fractions_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index] =
                confidence_stds * max(
                    sqrt(
                        linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] *
                        total_UMIs_per_metacell[metacell_index],
                    ),
                    min_confidence_UMIs,
                ) / total_UMIs_per_metacell[metacell_index]
        end
        return nothing
    end

    low_log_linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell))
    @check_turbo_matrix(low_log_linear_fraction_per_skeleton_per_metacell)
    @check_turbo_matrix(linear_fraction_per_skeleton_per_metacell)
    @check_turbo_matrix(confidence_linear_fractions_per_skeleton_per_metacells)
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "low_log_linear_fraction_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "low_log_linear_fraction_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            low_log_linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] = log2(
                max(
                    linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] -
                    confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index],
                    0.0,
                ) + gene_fraction_regularization,
            )
        end
        return nothing
    end

    high_log_linear_fraction_per_skeleton_per_metacell =
        Matrix{Float32}(undef, size(linear_fraction_per_skeleton_per_metacell))
    @check_turbo_matrix(high_log_linear_fraction_per_skeleton_per_metacell)
    @check_turbo_matrix(linear_fraction_per_skeleton_per_metacell)
    @check_turbo_matrix(confidence_linear_fractions_per_skeleton_per_metacells)
    parallel_loop_wo_rng(
        1:n_metacells;
        name = "high_log_linear_fraction_per_skeleton_per_metacell",
        progress = DebugProgress(
            n_metacells;
            group = :mcs_loops,
            desc = "high_log_linear_fraction_per_skeleton_per_metacell",
        ),
    ) do metacell_index
        @turbo for skeleton_index in 1:n_skeletons
            high_log_linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] = log2(
                linear_fraction_per_skeleton_per_metacell[skeleton_index, metacell_index] +
                confidence_linear_fractions_per_skeleton_per_metacells[skeleton_index, metacell_index] +
                gene_fraction_regularization,
            )
        end
        return nothing
    end

    return ConfidenceIntervals(;
        low_log_linear_fraction_per_skeleton_per_metacell,
        high_log_linear_fraction_per_skeleton_per_metacell,
    )
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

end
