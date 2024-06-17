"""
Given a set of raw metacells, partition them into boxes such that all metacells in the same box are within some
(fold factor) radius of each other. The centroids of these boxes can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states. Group these boxes in overlapping neighborhoods
of "similar" boxes for further analysis.
"""
module Boxes

export compute_boxes!

using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using Distributions
using Random
using SparseArrays
using Statistics

"""
    function compute_boxes!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = 40,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        fold_confidence::AbstractFloat = 0.9,
        max_box_span::AbstractFloat = 2.0,
        max_neighborhood_span::AbstractFloat = 2.0,
        correlation_confidence::AbstractFloat = 0.99,
        max_deviant_genes_fraction::AbstractFloat = 0.01,
        overwrite::Bool = false,
    )::Nothing

Partition raw metacells into distinct boxes, and boxes into overlapping neighborhoods.

 1. Initial boxes and neighborhoods are computed in a first round, and then refined in a series of followup rounds.
 2. In each round, we compute a distance between each two metacells. This is based on the fold factor between the
    expression level of each (relevant) gene in the metacells. The fold factor is the absolute value of the difference
    in the log (base 2) of the fraction of the gene in the metacells. This log is computed with the
    `gene_fraction_regularization` (by default, `1e-5`). Since the fraction of the gene is a random variable, we
    decrease the high fraction and increase the low fraction by a factor based on the `fold_confidence` of the test (by
    default, 0.9), assuming a multinomial distribution. In addition, if the sum of the total UMIs of the gene in both
    metacells is less than `min_significant_gene_UMIs` (by default, `40`), we ignore this fold factor as insignificant.
    Finally, genes, we reduce the fold factor using the gene's divergence factor. In the first round, we simply count
    the number of genes whose fold factor is more than `max_box_span` (for computing boxes) and `max_box_span` +
    `max_neighborhood_span` (for computing neighborhoods). In the followup rounds, we use the maximal gene fold, for
    genes that are correlated in the vicinity of the metacells (see below). Finally, in the followup rounds, we only
    compare a base metacell with other metacells which were in the same neighborhood in the previous round.
 3. We use hierarchical clustering to partition the metacells to distinct boxes, such that the maximal distance between
    any metacells in the box is bounded. In the first round, this bound is the `max_deviant_genes_fraction` out of the
    total number of genes. In the followup rounds, this is the `max_box_span`.
 4. For each box, we compute a main neighborhood of other boxes such that the maximal distance between any metacells
    in the neighborhood is bounded. In the first round, this bound is again the maximal number of deviant genes (this
    time, using the increased fold distance computed above). In the followup rounds, this is the `max_box_span`
    plus the `max_neighborhood_span`. These neighborhoods may overlap. The main neighborhoods of different boxes
    may even be identical.
 5. For each box, we compute the set of genes which are correlated (using the `correlation_confidence`) with some other
    gene(s) in its main neighborhood.
 6. If the new sets of correlated genes only differ up to `max_convergence_fraction`, consider this to have converged.
    Otherwise, repeat from step 2.

If `overwrite` is set, the results will replace any previously computed boxes and neighborhoods.

CONTRACT
"""
@logged @computation Contract(;
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (
            RequiredInput,
            "The genes to consider (typically, only non-lateral globally marker and correlated genes).",
        ),
        "box" => (GuaranteedOutput, "A partition of the metacells into distinct boxes."),
        "neighborhood" => (GuaranteedOutput, "A grouping of boxes into overlapping neighborhoods."),
    ],
    data = [
        ("metacell", "gene", "fraction") =>
            (RequiredInput, AbstractFloat, "The fraction of the UMIs of each gene in each metacell."),
        ("metacell", "total_UMIs") => (
            RequiredInput,
            Unsigned,
            "The total number of UMIs used to estimate the fraction of all the genes in each metacell.",
        ),
        ("metacell", "gene", "total_UMIs") => (
            RequiredInput,
            Unsigned,
            "The total number of UMIs used to estimate the fraction of each gene in each metacell.",
        ),
        ("gene", "divergence") => (OptionalInput, AbstractFloat, "How to scale fold factors for this gene."),
        ("metacell", "box") => (GuaranteedOutput, AbstractString, "The unique box each metacell belongs to."),
        ("gene", "neighborhood", "is_correlated") =>
            (GuaranteedOutput, Bool, "Which genes are correlated in each neighborhood."),
        ("box", "neighborhood.main") => (GuaranteedOutput, AbstractString, "The main neighborhood of each box."),
        ("neighborhood", "span") => (GuaranteedOutput, AbstractFloat, "The span used to compute the neighborhood."),
        ("box", "neighborhood", "is_member") =>
            (GuaranteedOutput, Bool, "Membership matrix for boxes and neighborhoods."),
        ("box", "box", "distance") =>
            (GuaranteedOutput, Float32, "For each box, the distances of the boxes in its main neighborhood."),
    ],
) function compute_boxes!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    fold_confidence::AbstractFloat = 0.9,
    max_box_span::AbstractFloat = 2.0,
    max_neighborhood_span::AbstractFloat = 2.0,
    target_boxes_in_neighborhood::Integer = 20,
    correlation_confidence::AbstractFloat = 0.9,
    max_deviant_genes_fraction::AbstractFloat = 0.01,
    max_convergence_fraction::AbstractFloat = 0.0015,
    overwrite::Bool = false,
    rng::Maybe{AbstractRNG} = nothing,
)::Nothing
    @assert min_significant_gene_UMIs >= 0
    @assert gene_fraction_regularization >= 0
    @assert 0 <= fold_confidence <= 1
    @assert max_box_span > 0
    @assert max_neighborhood_span >= 0
    @assert target_boxes_in_neighborhood > 0
    @assert 0 <= correlation_confidence <= 1
    @assert 0 <= max_deviant_genes_fraction <= 1

    total_UMIs_of_genes_in_metacells, fraction_of_genes_in_metacells, total_UMIs_of_metacells, divergence_of_genes =
        read_data(daf)

    log_decreased_fraction_of_genes_in_metacells, log_increased_fraction_of_genes_in_metacells =
        compute_confidence_log_fraction_of_genes_in_metacells(;
            gene_fraction_regularization = gene_fraction_regularization,
            fold_confidence = fold_confidence,
            fraction_of_genes_in_metacells = fraction_of_genes_in_metacells,
            total_UMIs_of_metacells = total_UMIs_of_metacells,
        )

    distance_for_boxes_between_metacells, distance_for_neighborhoods_between_metacells =
        compute_initial_metacells_distances(;
            total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
            log_decreased_fraction_of_genes_in_metacells = log_decreased_fraction_of_genes_in_metacells,
            log_increased_fraction_of_genes_in_metacells = log_increased_fraction_of_genes_in_metacells,
            divergence_of_genes = divergence_of_genes,
            min_significant_gene_UMIs = min_significant_gene_UMIs,
            max_box_span = max_box_span,
            max_neighborhood_span = max_neighborhood_span,
        )

    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")
    is_correlated_of_genes_in_metacells = fill(true, n_genes, n_metacells)
    previous_n_correlated_genes_in_metacells = n_genes * n_metacells
    max_deviant_genes = max_deviant_genes_fraction * n_genes
    @assert max_deviant_genes > 0
    boxes_threshold = max_deviant_genes
    neighborhood_threshold = max_deviant_genes

    boxes_of_metacells = nothing
    names_of_boxes = nothing
    distances_between_boxes = nothing
    main_neighborhoods_of_boxes = nothing
    names_of_neighborhoods = nothing
    boxes_of_neighborhoods = nothing
    spans_of_neighborhoods = nothing
    is_correlated_of_genes_in_neighborhoods = nothing

    round_index = 0
    metacells_of_boxes = nothing
    while true
        round_index += 1

        boxes_of_metacells = cluster_metacells_to_boxes(;
            distances_between_metacells = distance_for_boxes_between_metacells,
            threshold = boxes_threshold,
        )
        metacells_of_boxes = collect_group_members(boxes_of_metacells)
        names_of_boxes = group_names(daf, "metacell", metacells_of_boxes; prefix = "B")
        @debug "round: $(round_index) n_boxes: $(length(names_of_boxes))"

        distances_between_boxes = compute_box_distances(;
            distances_between_metacells = distance_for_neighborhoods_between_metacells,
            metacells_of_boxes = metacells_of_boxes,
        )

        spans_of_neighborhoods, boxes_of_neighborhoods, main_neighborhoods_of_boxes = collect_boxes_of_neighborhoods(;
            distances_between_boxes = distances_between_boxes,
            threshold = neighborhood_threshold,
            target_boxes_in_neighborhood = target_boxes_in_neighborhood,
            is_first = round_index == 1,
        )

        names_of_neighborhoods = group_names(daf, "metacell", boxes_of_neighborhoods; prefix = "N")
        @debug "round: $(round_index) n_neighborhoods: $(length(names_of_neighborhoods))"
        @assert length(names_of_neighborhoods) <= length(names_of_boxes)

        intermediate_daf = chain_writer([daf, MemoryDaf(; name = "intermediate")])
        add_axis!(intermediate_daf, "box", names_of_boxes)
        add_axis!(intermediate_daf, "neighborhood", names_of_neighborhoods)
        set_vector!(intermediate_daf, "metacell", "box", names_of_boxes[boxes_of_metacells]; overwrite = overwrite)
        set_vector!(intermediate_daf, "box", "neighborhood.main", names_of_neighborhoods[main_neighborhoods_of_boxes])

        is_correlated_of_genes_in_neighborhoods = identify_neighborhoods_correlated_genes(
            intermediate_daf;
            gene_fraction_regularization = gene_fraction_regularization,
            correlation_confidence = correlation_confidence,
            is_correlated_of_genes_in_metacells = is_correlated_of_genes_in_metacells,
            rng = rng,
        )
        next_n_correlated_genes_in_metacells = sum(is_correlated_of_genes_in_metacells)
        delta =
            (next_n_correlated_genes_in_metacells - previous_n_correlated_genes_in_metacells) /
            previous_n_correlated_genes_in_metacells

        @debug "round: $(round_index) is_correlated_genes_in_metacells: $(depict(is_correlated_of_genes_in_metacells)) (delta: $(delta))"
        if abs(delta) <= max_convergence_fraction
            break
        end
        previous_n_correlated_genes_in_metacells = next_n_correlated_genes_in_metacells

        distance_for_boxes_between_metacells =
            distance_for_neighborhoods_between_metacells = compute_followup_metacells_distances(;
                total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
                log_decreased_fraction_of_genes_in_metacells = log_decreased_fraction_of_genes_in_metacells,
                log_increased_fraction_of_genes_in_metacells = log_increased_fraction_of_genes_in_metacells,
                metacells_of_boxes = metacells_of_boxes,
                main_neighborhoods_of_boxes = main_neighborhoods_of_boxes,
                is_correlated_of_genes_in_neighborhoods = is_correlated_of_genes_in_neighborhoods,
                divergence_of_genes = divergence_of_genes,
                min_significant_gene_UMIs = min_significant_gene_UMIs,
            )

        boxes_threshold = max_box_span
        neighborhood_threshold = max_box_span + max_neighborhood_span
    end

    is_member_of_boxes_in_neighborhoods = compute_membership_matrix(boxes_of_neighborhoods, length(names_of_boxes))

    write_data(
        daf;
        names_of_boxes = names_of_boxes,
        names_of_neighborhoods = names_of_neighborhoods,
        boxes_of_metacells = boxes_of_metacells,
        main_neighborhoods_of_boxes = main_neighborhoods_of_boxes,
        is_member_of_boxes_in_neighborhoods = is_member_of_boxes_in_neighborhoods,
        is_correlated_of_genes_in_neighborhoods = is_correlated_of_genes_in_neighborhoods,
        spans_of_neighborhoods = spans_of_neighborhoods,
        distances_between_boxes = distances_between_boxes,
        overwrite = overwrite,
    )

    return nothing
end

@logged function read_data(  # untested
    daf::DafReader,
)::Tuple{
    AbstractMatrix{<:Unsigned},
    AbstractMatrix{<:AbstractFloat},
    AbstractVector{<:Unsigned},
    AbstractVector{<:AbstractFloat},
}
    total_UMIs_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    fraction_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    total_UMIs_of_metacells = get_vector(daf, "metacell", "total_UMIs").array
    divergence_of_genes = get_vector(daf, "gene", "divergence").array

    return total_UMIs_of_genes_in_metacells,
    fraction_of_genes_in_metacells,
    total_UMIs_of_metacells,
    divergence_of_genes
end

@logged function write_data(  # untested
    daf::DafWriter;
    names_of_boxes::AbstractVector{<:AbstractString},
    names_of_neighborhoods::AbstractVector{<:AbstractString},
    boxes_of_metacells::Vector{UInt32},
    main_neighborhoods_of_boxes::Vector{UInt32},
    is_member_of_boxes_in_neighborhoods::AbstractMatrix{Bool},
    is_correlated_of_genes_in_neighborhoods::AbstractMatrix{Bool},
    spans_of_neighborhoods::Vector{Float32},
    distances_between_boxes::Matrix{<:AbstractFloat},
    overwrite::Bool,
)::Nothing
    check_efficient_action("compute_boxes", Columns, "distances_between_boxes", distances_between_boxes)

    if overwrite
        delete_axis!(daf, "box"; must_exist = false)
        delete_axis!(daf, "neighborhood"; must_exist = false)
    end
    add_axis!(daf, "box", names_of_boxes)
    add_axis!(daf, "neighborhood", names_of_neighborhoods)

    set_vector!(daf, "metacell", "box", names_of_boxes[boxes_of_metacells]; overwrite = overwrite)
    set_vector!(daf, "box", "neighborhood.main", names_of_neighborhoods[main_neighborhoods_of_boxes])
    set_matrix!(daf, "box", "neighborhood", "is_member", SparseMatrixCSC(is_member_of_boxes_in_neighborhoods))
    set_matrix!(daf, "gene", "neighborhood", "is_correlated", SparseMatrixCSC(is_correlated_of_genes_in_neighborhoods))
    set_matrix!(daf, "box", "box", "distance", distances_between_boxes)
    set_vector!(daf, "neighborhood", "span", spans_of_neighborhoods)

    return nothing
end

@logged function compute_confidence_log_fraction_of_genes_in_metacells(;  # untested
    gene_fraction_regularization::AbstractFloat,
    fold_confidence::AbstractFloat,
    fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_metacells::AbstractVector{<:Unsigned},
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action("compute_boxes", Columns, "fraction_of_genes_in_metacells", fraction_of_genes_in_metacells)

    confidence_stdevs = quantile(Normal(), fold_confidence)

    confidence_fraction_of_genes_in_metacells =
        confidence_stdevs .* sqrt.(transpose(total_UMIs_of_metacells) .* fraction_of_genes_in_metacells) ./
        transpose(total_UMIs_of_metacells)

    log_decreased_fraction_of_genes_in_metacells =
        log2.(
            max.(fraction_of_genes_in_metacells .- confidence_fraction_of_genes_in_metacells, 0.0) .+
            gene_fraction_regularization
        )

    log_increased_fraction_of_genes_in_metacells =
        log2.(
            fraction_of_genes_in_metacells .+ confidence_fraction_of_genes_in_metacells .+ gene_fraction_regularization
        )

    return log_decreased_fraction_of_genes_in_metacells, log_increased_fraction_of_genes_in_metacells
end

@logged function compute_initial_metacells_distances(;  # untested
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    log_decreased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    log_increased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    divergence_of_genes::AbstractVector{<:AbstractFloat},
    min_significant_gene_UMIs::Integer,
    max_box_span::AbstractFloat,
    max_neighborhood_span::AbstractFloat,
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action(
        "compute_boxes",
        Columns,
        "total_UMIs_of_genes_in_metacells",
        total_UMIs_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_boxes",
        Columns,
        "log_decreased_fraction_of_genes_in_metacells",
        log_decreased_fraction_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_boxes",
        Columns,
        "log_increased_fraction_of_genes_in_metacells",
        log_increased_fraction_of_genes_in_metacells,
    )
    @assert size(log_decreased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)
    @assert size(log_increased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)

    n_genes = size(total_UMIs_of_genes_in_metacells, 1)
    n_metacells = size(total_UMIs_of_genes_in_metacells, 2)

    distance_for_boxes_between_metacells = zeros(Float32, n_metacells, n_metacells)
    distance_for_neighborhoods_between_metacells = zeros(Float32, n_metacells, n_metacells)

    @threads for base_metacell in reverse(2:n_metacells)
        @views total_UMIs_of_genes_in_base_metacell = total_UMIs_of_genes_in_metacells[:, base_metacell]
        @views log_decreased_fraction_of_genes_in_base_metacell =
            log_decreased_fraction_of_genes_in_metacells[:, base_metacell]
        @views log_increased_fraction_of_genes_in_base_metacell =
            log_increased_fraction_of_genes_in_metacells[:, base_metacell]

        @views total_UMIs_of_genes_in_other_metacells = total_UMIs_of_genes_in_metacells[:, 1:(base_metacell - 1)]
        @views log_decreased_fraction_of_genes_in_other_metacells =
            log_decreased_fraction_of_genes_in_metacells[:, 1:(base_metacell - 1)]
        @views log_increased_fraction_of_genes_in_other_metacells =
            log_increased_fraction_of_genes_in_metacells[:, 1:(base_metacell - 1)]

        significant_fold_factors_of_genes_in_other_metacells =
            gene_distance.(
                min_significant_gene_UMIs,
                total_UMIs_of_genes_in_base_metacell,
                log_decreased_fraction_of_genes_in_base_metacell,
                log_increased_fraction_of_genes_in_base_metacell,
                total_UMIs_of_genes_in_other_metacells,
                log_decreased_fraction_of_genes_in_other_metacells,
                log_increased_fraction_of_genes_in_other_metacells,
                divergence_of_genes,
            )
        @assert size(significant_fold_factors_of_genes_in_other_metacells) == (n_genes, base_metacell - 1)

        is_deviant_box_of_gene_in_other_metacells = significant_fold_factors_of_genes_in_other_metacells .> max_box_span
        is_deviant_neighborhood_of_genes_in_other_metacells =
            significant_fold_factors_of_genes_in_other_metacells .> max_box_span + max_neighborhood_span

        n_deviant_box_genes_of_other_metacells = vec(sum(is_deviant_box_of_gene_in_other_metacells; dims = 1))
        n_deviant_neighborhood_genes_of_other_metacells =
            vec(sum(is_deviant_neighborhood_of_genes_in_other_metacells; dims = 1))

        @assert length(n_deviant_box_genes_of_other_metacells) == base_metacell - 1
        @assert length(n_deviant_neighborhood_genes_of_other_metacells) == base_metacell - 1

        distance_for_boxes_between_metacells[1:(base_metacell - 1), base_metacell] .=
            n_deviant_box_genes_of_other_metacells
        distance_for_neighborhoods_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_neighborhood_genes_of_other_metacells

        distance_for_boxes_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_box_genes_of_other_metacells
        distance_for_neighborhoods_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_neighborhood_genes_of_other_metacells
    end

    return distance_for_boxes_between_metacells, distance_for_neighborhoods_between_metacells
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Integer,
    total_UMIs_of_gene_in_base_metacell::Integer,
    log_decreased_fraction_of_gene_in_base_metacell::AbstractFloat,
    log_increased_fraction_of_gene_in_base_metacell::AbstractFloat,
    total_UMIs_of_genes_in_other_metacell::Integer,
    log_decreased_fraction_of_gene_in_other_metacell::AbstractFloat,
    log_increased_fraction_of_gene_in_other_metacell::AbstractFloat,
    divergence::AbstractFloat,
)::AbstractFloat
    total_UMIs_of_gene = total_UMIs_of_gene_in_base_metacell + total_UMIs_of_genes_in_other_metacell
    is_significant = total_UMIs_of_gene >= min_significant_gene_UMIs

    is_base_low = log_increased_fraction_of_gene_in_base_metacell < log_increased_fraction_of_gene_in_other_metacell

    log_increased_low_fraction_of_gene =
        is_base_low * log_increased_fraction_of_gene_in_base_metacell +
        !is_base_low * log_increased_fraction_of_gene_in_other_metacell

    log_decreased_high_fraction_of_gene =
        is_base_low * log_decreased_fraction_of_gene_in_other_metacell +
        !is_base_low * log_decreased_fraction_of_gene_in_base_metacell

    return (
        is_significant *
        max.((log_decreased_high_fraction_of_gene - log_increased_low_fraction_of_gene) * (1.0 - divergence), 0.0)
    )
end

@logged function cluster_metacells_to_boxes(;  # untested
    distances_between_metacells::Matrix{Float32},
    threshold::AbstractFloat,
)::Vector{UInt32}
    check_efficient_action("compute_boxes", Columns, "distances_between_metacells", distances_between_metacells)
    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    boxes_of_metacells = Vector{UInt32}(cutree(clusters; h = threshold))
    return boxes_of_metacells
end

@logged function compute_box_distances(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacells_of_boxes::Vector{Vector{UInt32}},
)::Matrix{Float32}
    check_efficient_action("compute_boxes", Columns, "distances_between_metacells", distances_between_metacells)

    n_metacells = size(distances_between_metacells, 1)
    @assert size(distances_between_metacells, 2) == n_metacells

    n_boxes = length(metacells_of_boxes)
    distances_between_boxes = Matrix{Float32}(undef, n_boxes, n_boxes)

    @threads for base_box in reverse(1:n_boxes)
        metacells_of_base_box = metacells_of_boxes[base_box]

        @views distance_of_metacells_from_base_box_metacells = distances_between_metacells[:, metacells_of_base_box]

        distance_of_metacells_from_base_box = maximum(distance_of_metacells_from_base_box_metacells; dims = 2)
        @assert length(distance_of_metacells_from_base_box) == n_metacells

        for other_box in 1:base_box
            metacells_of_other_box = metacells_of_boxes[other_box]

            @views distance_of_other_box_metacells_from_base_box =
                distance_of_metacells_from_base_box[metacells_of_other_box]

            distance_between_other_and_base_boxes = maximum(distance_of_other_box_metacells_from_base_box)

            distances_between_boxes[base_box, other_box] = distance_between_other_and_base_boxes
            distances_between_boxes[other_box, base_box] = distance_between_other_and_base_boxes
        end
    end

    return distances_between_boxes
end

@logged function collect_boxes_of_neighborhoods(;  # untested
    distances_between_boxes::Matrix{Float32},
    threshold::AbstractFloat,
    target_boxes_in_neighborhood::Integer,
    is_first::Bool,
)::Tuple{Vector{Float32}, Vector{Vector{UInt32}}, Vector{UInt32}}
    check_efficient_action("compute_boxes", Columns, "distances_between_boxes", distances_between_boxes)

    n_boxes = size(distances_between_boxes, 1)
    @assert size(distances_between_boxes, 2) == n_boxes

    neighborhood_spans_of_boxes = Vector{Float32}(undef, n_boxes)
    boxes_of_neighborhoods = Vector{Vector{UInt32}}(undef, n_boxes)

    boxes_quantile = min(target_boxes_in_neighborhood / n_boxes, 1.0)
    box_threshold = threshold
    for box_index in 1:n_boxes
        @views distances_from_box = distances_between_boxes[:, box_index]
        if !is_first
            target_threshold = quantile(distances_from_box, boxes_quantile)
            box_threshold = min(target_threshold, threshold)
        end
        boxes_of_neighborhoods[box_index] = findall(distances_from_box .<= box_threshold)
        neighborhood_spans_of_boxes[box_index] = box_threshold
    end

    main_neighborhoods_of_boxes = Vector{UInt32}(undef, n_boxes)
    spans_of_neighborhoods = Vector{Float32}(undef, n_boxes)

    if n_boxes <= 1
        fill!(main_neighborhoods_of_boxes, 1)
        unique_boxes_of_neighborhoods = boxes_of_neighborhoods
        spans_of_neighborhoods = neighborhood_spans_of_boxes
    else
        seen_boxes_of_neighborhoods = Dict(boxes_of_neighborhoods[1] => 1)
        unique_boxes_of_neighborhoods = similar(boxes_of_neighborhoods)
        unique_boxes_of_neighborhoods[1] = boxes_of_neighborhoods[1]
        spans_of_neighborhoods = similar(neighborhood_spans_of_boxes)
        main_neighborhoods_of_boxes[1] = 1
        spans_of_neighborhoods[1] = neighborhood_spans_of_boxes[1]
        for box_index in 2:n_boxes
            boxes_of_neighborhood = boxes_of_neighborhoods[box_index]
            neighborhood_index = get(seen_boxes_of_neighborhoods, boxes_of_neighborhood, nothing)
            if neighborhood_index === nothing
                neighborhood_index = length(seen_boxes_of_neighborhoods) + 1
                seen_boxes_of_neighborhoods[boxes_of_neighborhood] = neighborhood_index
                unique_boxes_of_neighborhoods[neighborhood_index] = boxes_of_neighborhood
                spans_of_neighborhoods[neighborhood_index] = neighborhood_spans_of_boxes[box_index]
            end
            main_neighborhoods_of_boxes[box_index] = neighborhood_index
        end
        resize!(unique_boxes_of_neighborhoods, length(seen_boxes_of_neighborhoods))
        resize!(spans_of_neighborhoods, length(seen_boxes_of_neighborhoods))
    end

    return spans_of_neighborhoods, unique_boxes_of_neighborhoods, main_neighborhoods_of_boxes
end

@logged function identify_neighborhoods_correlated_genes(  # untested
    daf::DafReader;
    gene_fraction_regularization::AbstractFloat,
    correlation_confidence::AbstractFloat,
    is_correlated_of_genes_in_metacells::Matrix{Bool},
    rng::Maybe{AbstractRNG},
)::Matrix{Bool}
    check_efficient_action(
        "compute_boxes",
        Columns,
        "is_correlated_of_genes_in_metacells",
        is_correlated_of_genes_in_metacells,
    )
    n_genes = axis_length(daf, "gene")
    n_neighborhoods = axis_length(daf, "neighborhood")
    names_of_neighborhoods = axis_array(daf, "neighborhood")

    is_correlated_of_genes_in_neighborhoods = zeros(Bool, n_genes, n_neighborhoods)
    indices_of_metacells = axis_dict(daf, "metacell")

    @threads for neighborhood_index in 1:n_neighborhoods
        neighborhood_name = names_of_neighborhoods[neighborhood_index]

        neighborhood_daf = chain_writer(
            [
                viewer(
                    daf;
                    name = "intermediate.$(neighborhood_name)",
                    axes = [
                        "metacell" => "/ metacell & box => neighborhood.main = $(neighborhood_name)",
                        "gene" => "=",
                    ],
                    data = [("metacell", "gene", "fraction") => "="],
                ),
                MemoryDaf(; name = "intermediate.$(neighborhood_name).identified"),
            ];
            name = "intermediate.$(neighborhood_name).chain",
        )
        names_of_metacells_in_neighborhood = axis_array(neighborhood_daf, "metacell")
        indices_of_metacells_in_neighborhood =
            [indices_of_metacells[metacell_name] for metacell_name in names_of_metacells_in_neighborhood]

        if length(indices_of_metacells_in_neighborhood) > 1
            identify_correlated_genes!(
                neighborhood_daf;
                gene_fraction_regularization = gene_fraction_regularization,
                correlation_confidence = correlation_confidence,
                rng = rng,
            )
            is_correlated_of_genes_in_neighborhood = neighborhood_daf["/ gene : is_correlated"].array
        else
            @assert length(indices_of_metacells_in_neighborhood) == 1
            metacell_index = indices_of_metacells_in_neighborhood[1]
            is_correlated_of_genes_in_neighborhood = vec(is_correlated_of_genes_in_metacells[:, metacell_index])
        end

        is_correlated_of_genes_in_neighborhoods[:, neighborhood_index] .= is_correlated_of_genes_in_neighborhood
        is_correlated_of_genes_in_metacells[:, indices_of_metacells_in_neighborhood] .=
            is_correlated_of_genes_in_neighborhood
    end

    return is_correlated_of_genes_in_neighborhoods
end

@logged function compute_followup_metacells_distances(;  # untested
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    log_decreased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    log_increased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    metacells_of_boxes::Vector{Vector{UInt32}},
    main_neighborhoods_of_boxes::Vector{UInt32},
    is_correlated_of_genes_in_neighborhoods::Matrix{Bool},
    divergence_of_genes::AbstractVector{<:AbstractFloat},
    min_significant_gene_UMIs::Integer,
)::Matrix{Float32}
    check_efficient_action(
        "compute_boxes",
        Columns,
        "total_UMIs_of_genes_in_metacells",
        total_UMIs_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_boxes",
        Columns,
        "log_decreased_fraction_of_genes_in_metacells",
        log_decreased_fraction_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_boxes",
        Columns,
        "log_increased_fraction_of_genes_in_metacells",
        log_increased_fraction_of_genes_in_metacells,
    )
    @assert size(log_decreased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)
    @assert size(log_increased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)

    n_metacells = size(total_UMIs_of_genes_in_metacells, 2)
    distances_between_metacells = fill(Inf32, n_metacells, n_metacells)

    n_boxes = length(metacells_of_boxes)
    did_compute_distance_between_boxes = zeros(Bool, n_boxes, n_boxes)

    @threads for base_box_index in 1:n_boxes
        metacell_indices_of_base_box = metacells_of_boxes[base_box_index]
        main_neighborhood_of_base_box = main_neighborhoods_of_boxes[base_box_index]
        @views is_correlated_of_genes_in_base_box =
            is_correlated_of_genes_in_neighborhoods[:, main_neighborhood_of_base_box]

        for other_box_index in 1:base_box_index
            if did_compute_distance_between_boxes[base_box_index, other_box_index]
                continue
            end
            did_compute_distance_between_boxes[base_box_index, other_box_index] = true
            did_compute_distance_between_boxes[other_box_index, base_box_index] = true

            metacell_indices_of_other_box = metacells_of_boxes[other_box_index]
            n_metacells_of_other_box = length(metacell_indices_of_other_box)

            main_neighborhood_of_other_box = main_neighborhoods_of_boxes[other_box_index]
            @views is_correlated_of_genes_in_other_box =
                is_correlated_of_genes_in_neighborhoods[:, main_neighborhood_of_other_box]
            is_correlated_of_genes_for_distance =
                is_correlated_of_genes_in_base_box .| is_correlated_of_genes_in_other_box
            if !any(is_correlated_of_genes_for_distance)
                distances_between_metacells[metacell_indices_of_base_box, metacell_indices_of_other_box] .= Inf
                distances_between_metacells[metacell_indices_of_other_box, metacell_indices_of_base_box] .= Inf
                continue
            end
            divergence_of_genes_for_distance = divergence_of_genes[is_correlated_of_genes_for_distance]

            for base_metacell_index in metacell_indices_of_base_box
                total_UMIs_of_genes_in_base_metacell =
                    total_UMIs_of_genes_in_metacells[is_correlated_of_genes_for_distance, base_metacell_index]
                log_decreased_fraction_of_genes_in_base_metacell = log_decreased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    base_metacell_index,
                ]
                log_increased_fraction_of_genes_in_base_metacell = log_increased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    base_metacell_index,
                ]

                total_UMIs_of_genes_in_other_metacells =
                    total_UMIs_of_genes_in_metacells[is_correlated_of_genes_for_distance, metacell_indices_of_other_box]
                log_decreased_fraction_of_genes_in_other_metacells = log_decreased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    metacell_indices_of_other_box,
                ]
                log_increased_fraction_of_genes_in_other_metacells = log_increased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    metacell_indices_of_other_box,
                ]

                significant_fold_factors_of_genes_in_other_metacells =
                    gene_distance.(
                        min_significant_gene_UMIs,
                        total_UMIs_of_genes_in_base_metacell,
                        log_decreased_fraction_of_genes_in_base_metacell,
                        log_increased_fraction_of_genes_in_base_metacell,
                        total_UMIs_of_genes_in_other_metacells,
                        log_decreased_fraction_of_genes_in_other_metacells,
                        log_increased_fraction_of_genes_in_other_metacells,
                        divergence_of_genes_for_distance,
                    )
                @assert size(significant_fold_factors_of_genes_in_other_metacells, 2) == n_metacells_of_other_box

                maximal_fold_factors_of_genes_in_other_metacells =
                    vec(maximum(significant_fold_factors_of_genes_in_other_metacells; dims = 1))
                @assert length(maximal_fold_factors_of_genes_in_other_metacells) == n_metacells_of_other_box

                distances_between_metacells[base_metacell_index, metacell_indices_of_other_box] .=
                    maximal_fold_factors_of_genes_in_other_metacells
                distances_between_metacells[metacell_indices_of_other_box, base_metacell_index] .=
                    maximal_fold_factors_of_genes_in_other_metacells
            end
        end
    end
    return distances_between_metacells
end

@logged function compute_membership_matrix(  # untested
    boxes_of_neighborhoods::Vector{Vector{UInt32}},
    n_boxes::Integer,
)::AbstractMatrix{Bool}
    n_neighborhoods = length(boxes_of_neighborhoods)
    is_member_of_boxes_in_neighborhoods = zeros(Bool, n_boxes, n_neighborhoods)
    for (neighborhood, boxes_of_neighborhood) in enumerate(boxes_of_neighborhoods)
        is_member_of_boxes_in_neighborhoods[boxes_of_neighborhood, neighborhood] .= true  # NOJET
    end
    return is_member_of_boxes_in_neighborhoods
end

end  # module
