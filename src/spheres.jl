"""
Given a set of raw metacells, partition them into spheres such that all metacells in the same sphere are within some
(fold factor) radius of each other. The centroids of these spheres can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states. Group these spheres in overlapping neighborhoods
of "similar" spheres for further analysis.
"""
module Spheres

export compute_improved_spheres!
export compute_spheres!
export improve_spheres!

using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using Distributions
using Statistics

SPHERE_CONTRACT = Contract(;
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (RequiredInput, "The genes to consider (typically, only marker genes)."),
        "sphere" => (GuaranteedOutput, "A partition of the metacells into distinct spheres."),
        "neighborhood" => (GuaranteedOutput, "A grouping of spheres into overlapping neighborhoods."),
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
        ("metacell", "gene", "is_marker") =>
            (OptionalInput, Bool, "The genes that are distinguish metacells in the vicinity of each metacell."),
        ("metacell", "gene", "is_lonely") =>
            (OptionalInput, Bool, "The genes that are uncorrelated in the vicinity of each metacell."),
        ("gene", "is_noisy") =>
            (OptionalInput, Bool, "A mask of noisy genes to be given additional diameter when grouped."),
        ("metacell", "sphere") => (GuaranteedOutput, AbstractString, "The unique sphere each metacell belongs to."),
        ("sphere", "neighborhood.main") =>
            (GuaranteedOutput, AbstractString, "The main neighborhood of each sphere."),
        ("sphere", "neighborhood", "is_member") =>
            (GuaranteedOutput, Bool, "Membership matrix for spheres and neighborhoods."),
        ("sphere", "sphere", "distance") =>
            (GuaranteedOutput, Float32, "For each sphere, the distances of the spheres in its main neighborhood."),
    ],
)

IMPROVED_SPHERE_CONTRACT = Contract(;
    axes = SPHERE_CONTRACT.axes,
    data = vcat(
        SPHERE_CONTRACT.data,
        [
            ("metacell", "gene", "is_marker") => (
                GuaranteedOutput,
                Bool,
                "A mask of genes that distinguish between metacells in the vicinity of the metacell.",
            ),
            ("metacell", "gene", "is_lonely") => (
                GuaranteedOutput,
                Bool,
                "A mask of genes that are uncorrelated with other genes in the vicinity of the metacell.",
            ),
        ],
    ),
)

"""
    function compute_spheres!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = 40,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        confidence::AbstractFloat = 0.9,
        max_sphere_diameter::AbstractFloat = 2.0,
        max_neighborhood_diameter::AbstractFloat = 2.0,
        noisy_gene_fold::AbstractFloat = 1.0,
        max_deviant_genes::Integer = 0,
        overwrite::Bool = false,
    )::Nothing

Partition raw metacells into distinct spheres, and spheres into overlapping neighborhoods.

!!! note

    This is a lower-level function; you probably want to call [`compute_improved_spheres!`](@ref) instead.

 1. Compute a distance between each two metacells, defined as the maximal fold factor of any of the genes. This fold
    factor is the absolute value of the difference in the log base 2 of the fraction of the gene in both metacells,
    using the `gene_fraction_regularization` (by default, `1e-5`). Since the fraction of the gene is a random variable,
    we decrease the high fraction and increase the low fraction by a factor based on the `confidence` of the test (by
    default, 0.9), assuming a multinomial distribution. In addition, if the sum of the total UMIs of the gene in both
    metacells is less than `min_significant_gene_UMIs` (by default, `40`), we ignore this fold factor as insignificant.
    We also ignore the fold factors for genes that aren't `is_marker` in the vicinity of at least one of the metacells,
    or that are `is_lonely` in the vicinity of both. For `is_noisy_gene_fold` genes, we reduce the distance by
    `noisy_gene_fold` to allow for their bursty nature.
 2. If `max_deviant_genes` is not zero, we redefine the distance between metacells to be the number of genes with a maximal
    fold distance above the `max_sphere_diameter` (for computing spheres) and `max_neighborhood_diameter` (for computing
    neighborhoods). Otherwise, we use the raw maximal fold distance.
 3. Compute hierarchical clustering of the metacells using these distances and the `complete` linkage; that is, each
    cluster can be seen as a sphere of some diameter where all the metacells are within that sphere.
 4. Use this hierarchical clustering to partition the metacells into as-large as possible spheres. If `max_deviant_genes`
    is zero (the default), then this is `cutree` using the `max_sphere_diameter` (by default, `2.0`, or 4x); otherwise,
    it is a `cutree` using the `max_deviant_genes`.
 5. Compute for each sphere the set of other spheres in its main neighborhood. If `max_deviant_genes` is zero, this
    means the maximal distance between any pair of metacells is at most an additional `max_neighborhood_diameter` (by
    default, `2.0`, for a total neighborhood diameter of 16x). Otherwise, it means the maximal number of genes whose
    distance is more than this threshold in the metacells pair is at most `max_deviant_genes`. Neighborhoods can
    overlap; also, if the main neighborhoods defined by two spheres end up being identical, we unify them, so we may end
    up with less neighborhoods than spheres.
 6. The distances to the other spheres in the main neighborhood of each sphere form a natural graph representation of
    the structure of the manifold of the distinct biological cell states captured by the metacells.

If `overwrite` is set, the results will replace any previously computed spheres and neighborhoods.

CONTRACT
"""
@logged @computation SPHERE_CONTRACT function compute_spheres!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    confidence::AbstractFloat = 0.9,
    max_sphere_diameter::AbstractFloat = 2.0,
    max_neighborhood_diameter::AbstractFloat = 2.0,
    noisy_gene_fold::AbstractFloat = 1.0,
    max_deviant_genes::Integer = 0,
    overwrite::Bool = false,
)::Nothing
    @assert min_significant_gene_UMIs >= 0
    @assert gene_fraction_regularization > 0
    @assert confidence > 0.0
    @assert confidence < 1.0
    @assert max_sphere_diameter > 0
    @assert max_neighborhood_diameter >= 0.0
    @assert max_deviant_genes >= 0

    genes_metacells_total_UMIs,
    genes_metacells_fraction,
    genes_metacells_is_marker,
    genes_metacells_is_lonely,
    metacells_total_UMIs,
    genes_is_noisy = read_data(daf)

    genes_metacells_log_decreased_fraction, genes_metacells_log_increased_fraction =
        prepare_metacells_data(gene_fraction_regularization, confidence, genes_metacells_fraction, metacells_total_UMIs)

    metacells_distances_for_spheres, metacells_distances_for_neighborhoods = compute_metacells_distances(  # NOJET
        min_significant_gene_UMIs,
        genes_metacells_total_UMIs,
        genes_metacells_is_marker,
        genes_metacells_is_lonely,
        genes_metacells_log_decreased_fraction,
        genes_metacells_log_increased_fraction,
        max_sphere_diameter,
        max_neighborhood_diameter,
        noisy_gene_fold,
        genes_is_noisy,
        max_deviant_genes,
    )

    spheres_of_metacells =
        compute_spheres_of_metacells(metacells_distances_for_spheres, max_sphere_diameter, max_deviant_genes)
    metacells_of_spheres = collect_group_members(spheres_of_metacells)
    sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")

    spheres_distances = compute_sphere_distances(
        metacells_of_spheres,
        metacells_distances_for_neighborhoods,
        max_deviant_genes == 0 ? max_sphere_diameter : nothing,
    )
    spheres_of_neighborhoods, main_neighborhoods_of_spheres = collect_spheres_of_neighborhoods(
        spheres_distances,
        max_sphere_diameter,
        max_neighborhood_diameter,
        max_deviant_genes,
    )
    neighborhood_names = group_names(daf, "metacell", spheres_of_neighborhoods; prefix = "N")
    spheres_neighborhoods_is_member = compute_membership_matrix(spheres_of_neighborhoods, length(sphere_names))
    @assert length(neighborhood_names) <= length(sphere_names)

    write_data(
        daf,
        sphere_names,
        neighborhood_names,
        spheres_of_metacells,
        main_neighborhoods_of_spheres,
        spheres_neighborhoods_is_member,
        spheres_of_neighborhoods,
        spheres_distances,
        overwrite,
    )

    if max_deviant_genes > 0
        return nothing
    end

    return nothing
end

@logged function read_data(  # untested
    daf::DafReader,
)::Tuple{
    AbstractMatrix{<:Unsigned},
    AbstractMatrix{<:AbstractFloat},
    AbstractMatrix{Bool},
    AbstractMatrix{Bool},
    AbstractVector{<:Unsigned},
    AbstractVector{Bool},
}
    genes_metacells_total_UMIs = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    genes_metacells_fraction = get_matrix(daf, "gene", "metacell", "fraction").array
    genes_metacells_is_marker = get_matrix(daf, "gene", "metacell", "is_marker"; default = true)
    genes_metacells_is_lonely = get_matrix(daf, "gene", "metacell", "is_lonely"; default = false)
    metacells_total_UMIs = get_vector(daf, "metacell", "total_UMIs").array
    genes_is_noisy = get_vector(daf, "gene", "is_noisy"; default = false).array
    return (
        genes_metacells_total_UMIs,
        genes_metacells_fraction,
        genes_metacells_is_marker,
        genes_metacells_is_lonely,
        metacells_total_UMIs,
        genes_is_noisy,
    )
end

@logged function prepare_metacells_data(  # untested
    gene_fraction_regularization::AbstractFloat,
    confidence::AbstractFloat,
    genes_metacells_fraction::AbstractMatrix{<:AbstractFloat},
    metacells_total_UMIs::AbstractVector{<:Unsigned},
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action("compute_spheres", Columns, "genes_metacells_fraction", genes_metacells_fraction)

    confidence_stdevs = quantile(Normal(), confidence)

    genes_metacells_confidence =
        confidence_stdevs .* sqrt.(transpose(metacells_total_UMIs) .* genes_metacells_fraction) ./
        transpose(metacells_total_UMIs)

    genes_metacells_log_decreased_fraction =
        log2.(max.(genes_metacells_fraction .- genes_metacells_confidence, 0.0) .+ gene_fraction_regularization)

    genes_metacells_log_increased_fraction =
        log2.(genes_metacells_fraction .+ genes_metacells_confidence .+ gene_fraction_regularization)

    return genes_metacells_log_decreased_fraction, genes_metacells_log_increased_fraction
end

@logged function compute_metacells_distances(  # untested
    min_significant_gene_UMIs::Integer,
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
    genes_metacells_is_marker::AbstractMatrix{Bool},
    genes_metacells_is_lonely::AbstractMatrix{Bool},
    genes_metacells_log_decreased_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_log_increased_fraction::AbstractMatrix{<:AbstractFloat},
    max_sphere_diameter::AbstractFloat,
    max_neighborhood_diameter::AbstractFloat,
    noisy_gene_fold::AbstractFloat,
    genes_is_noisy::AbstractVector{Bool},
    max_deviant_genes::Integer,
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action("compute_spheres", Columns, "genes_metacells_total_UMIs", genes_metacells_total_UMIs)
    check_efficient_action(
        "compute_spheres",
        Columns,
        "genes_metacells_log_decreased_fraction",
        genes_metacells_log_decreased_fraction,
    )
    check_efficient_action(
        "compute_spheres",
        Columns,
        "genes_metacells_log_increased_fraction",
        genes_metacells_log_increased_fraction,
    )
    @assert size(genes_metacells_log_decreased_fraction) == size(genes_metacells_total_UMIs)
    @assert size(genes_metacells_log_increased_fraction) == size(genes_metacells_total_UMIs)

    n_genes = size(genes_metacells_total_UMIs, 1)
    n_metacells = size(genes_metacells_total_UMIs, 2)
    metacells_distances_for_spheres = zeros(Float32, n_metacells, n_metacells)
    if max_deviant_genes > 0
        metacells_distances_for_neighborhoods = zeros(Float32, n_metacells, n_metacells)
    else
        metacells_distances_for_neighborhoods = metacells_distances_for_spheres
    end

    @threads for base_metacell in reverse(2:n_metacells)
        @views genes_total_UMIs_of_base_metacell = genes_metacells_total_UMIs[:, base_metacell]
        @views genes_total_UMIs_of_other_metacells = genes_metacells_total_UMIs[:, 1:(base_metacell - 1)]
        @views genes_log_decreased_fraction_of_base_metacell = genes_metacells_log_decreased_fraction[:, base_metacell]
        @views genes_log_increased_fraction_of_base_metacell = genes_metacells_log_increased_fraction[:, base_metacell]
        @views genes_log_decreased_fraction_of_other_metacells =
            genes_metacells_log_decreased_fraction[:, 1:(base_metacell - 1)]
        @views genes_log_increased_fraction_of_other_metacells =
            genes_metacells_log_increased_fraction[:, 1:(base_metacell - 1)]

        @views base_metacell_is_marker_of_genes = genes_metacells_is_marker[:, base_metacell]
        @views base_metacell_is_lonely_of_genes = genes_metacells_is_lonely[:, base_metacell]
        @views other_metacells_is_marker_of_genes = genes_metacells_is_marker[:, 1:(base_metacell - 1)]
        @views other_metacells_is_lonely_of_genes = genes_metacells_is_lonely[:, 1:(base_metacell - 1)]

        significant_fold_factors_of_other_metacells =
            gene_distance.(
                min_significant_gene_UMIs,
                genes_total_UMIs_of_base_metacell,
                genes_total_UMIs_of_other_metacells,
                base_metacell_is_marker_of_genes,
                base_metacell_is_lonely_of_genes,
                other_metacells_is_marker_of_genes,
                other_metacells_is_lonely_of_genes,
                genes_log_decreased_fraction_of_base_metacell,
                genes_log_increased_fraction_of_base_metacell,
                genes_log_decreased_fraction_of_other_metacells,
                genes_log_increased_fraction_of_other_metacells,
                noisy_gene_fold,
                genes_is_noisy,
            )
        @assert size(significant_fold_factors_of_other_metacells) == (n_genes, base_metacell - 1)

        if max_deviant_genes > 0
            far_sphere_genes_of_other_metacells = significant_fold_factors_of_other_metacells .> max_sphere_diameter
            far_neighborhood_genes_of_other_metacells =
                significant_fold_factors_of_other_metacells .> max_sphere_diameter + max_neighborhood_diameter

            count_far_sphere_genes_of_other_metacells = vec(sum(far_sphere_genes_of_other_metacells; dims = 1))
            count_far_neighborhood_genes_of_other_metacells =
                vec(sum(far_neighborhood_genes_of_other_metacells; dims = 1))

            @assert length(count_far_sphere_genes_of_other_metacells) == base_metacell - 1
            @assert length(count_far_neighborhood_genes_of_other_metacells) == base_metacell - 1

            metacells_distances_for_spheres[1:(base_metacell - 1), base_metacell] .=
                count_far_sphere_genes_of_other_metacells
            metacells_distances_for_spheres[base_metacell, 1:(base_metacell - 1)] .=
                count_far_sphere_genes_of_other_metacells

            metacells_distances_for_neighborhoods[1:(base_metacell - 1), base_metacell] .=
                count_far_neighborhood_genes_of_other_metacells
            metacells_distances_for_neighborhoods[base_metacell, 1:(base_metacell - 1)] .=
                count_far_neighborhood_genes_of_other_metacells

        else
            maximal_fold_factors_of_other_metacells =
                vec(maximum(significant_fold_factors_of_other_metacells; dims = 1))
            @assert length(maximal_fold_factors_of_other_metacells) == base_metacell - 1
            metacells_distances_for_spheres[1:(base_metacell - 1), base_metacell] .=
                maximal_fold_factors_of_other_metacells
            metacells_distances_for_spheres[base_metacell, 1:(base_metacell - 1)] .=
                maximal_fold_factors_of_other_metacells
        end
    end

    return metacells_distances_for_spheres, metacells_distances_for_neighborhoods
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Integer,
    gene_total_UMIs_of_base_metacell::Unsigned,
    gene_total_UMIs_of_other_metacell::Unsigned,
    base_metacell_gene_is_marker::Bool,
    base_metacell_gene_is_lonely::Bool,
    other_metacells_gene_is_marker::Bool,
    other_metacells_gene_is_lonely::Bool,
    gene_log_decreased_fraction_of_base_metacell::AbstractFloat,
    gene_log_increased_fraction_of_base_metacell::AbstractFloat,
    gene_log_decreased_fraction_of_other_metacell::AbstractFloat,
    gene_log_increased_fraction_of_other_metacell::AbstractFloat,
    noisy_gene_fold::AbstractFloat,
    gene_is_noisy::Bool,
)::AbstractFloat
    gene_total_UMIs = gene_total_UMIs_of_base_metacell + gene_total_UMIs_of_other_metacell
    is_significant =
        (gene_total_UMIs >= min_significant_gene_UMIs) &
        (base_metacell_gene_is_marker | other_metacells_gene_is_marker) &
        (!base_metacell_gene_is_lonely | !other_metacells_gene_is_lonely)
    is_base_low = gene_log_increased_fraction_of_base_metacell < gene_log_increased_fraction_of_other_metacell
    gene_log_low_fraction =
        is_base_low * gene_log_increased_fraction_of_base_metacell +
        !is_base_low * gene_log_increased_fraction_of_other_metacell
    gene_log_high_fraction =
        !is_base_low * gene_log_decreased_fraction_of_base_metacell +
        is_base_low * gene_log_decreased_fraction_of_other_metacell
    return is_significant * max.(gene_log_high_fraction - gene_log_low_fraction - gene_is_noisy * noisy_gene_fold, 0.0)
end

@logged function compute_sphere_distances(  # untested
    metacells_of_spheres::Vector{Vector{UInt32}},
    metacells_distances::Matrix{Float32},
    max_sphere_diameter::Maybe{AbstractFloat},
)::Matrix{Float32}
    check_efficient_action("compute_spheres", Columns, "metacells_distances", metacells_distances)

    n_spheres = length(metacells_of_spheres)
    spheres_distances = Matrix{eltype(metacells_distances)}(undef, n_spheres, n_spheres)

    @threads for base_sphere in reverse(1:n_spheres)
        metacells_of_base_sphere = metacells_of_spheres[base_sphere]
        @views base_sphere_distances_matrix = metacells_distances[:, metacells_of_base_sphere]
        base_sphere_distances_vector = maximum(base_sphere_distances_matrix; dims = 2)
        for other_sphere in 1:base_sphere
            metacells_of_other_sphere = metacells_of_spheres[other_sphere]
            @views other_sphere_distances_vector = base_sphere_distances_vector[metacells_of_other_sphere]
            distance = maximum(other_sphere_distances_vector)
            if max_sphere_diameter !== nothing
                if other_sphere == base_sphere
                    @assert distance <= max_sphere_diameter
                else
                    @assert distance > max_sphere_diameter
                end
            end
            spheres_distances[base_sphere, other_sphere] = distance
            spheres_distances[other_sphere, base_sphere] = distance
        end
    end

    return spheres_distances
end

@logged function compute_spheres_of_metacells(  # untested
    metacells_distances::Matrix{Float32},
    max_sphere_diameter::AbstractFloat,
    max_deviant_genes::Integer,
)::Vector{UInt32}
    check_efficient_action("compute_spheres", Columns, "metacells_distances", metacells_distances)
    clusters = hclust(metacells_distances; linkage = :complete)  # NOJET
    threshold = max_deviant_genes == 0 ? max_sphere_diameter : max_deviant_genes + 0.5
    spheres_of_metacells = Vector{UInt32}(cutree(clusters; h = threshold))
    return spheres_of_metacells
end

@logged function collect_spheres_of_neighborhoods(  # untested
    spheres_distances::Matrix{Float32},
    max_sphere_diameter::AbstractFloat,
    max_neighborhood_diameter::AbstractFloat,
    max_deviant_genes::Integer,
)::Tuple{Vector{Vector{<:Integer}}, Vector{UInt32}}
    n_spheres = size(spheres_distances, 1)
    threshold = max_deviant_genes == 0 ? max_sphere_diameter + max_neighborhood_diameter : max_deviant_genes
    spheres_of_neighborhoods = [findall(spheres_distances[:, sphere] .<= threshold) for sphere in 1:n_spheres]
    main_neighborhoods_of_spheres = Vector{UInt32}(undef, n_spheres)

    if n_spheres <= 1
        fill!(main_neighborhoods_of_spheres, 1)
        unique_spheres_of_neighborhoods = spheres_of_neighborhoods
    else
        seen_spheres_of_neighborhoods = Dict(spheres_of_neighborhoods[1] => 1)
        unique_spheres_of_neighborhoods = similar(spheres_of_neighborhoods)
        unique_spheres_of_neighborhoods[1] = spheres_of_neighborhoods[1]
        main_neighborhoods_of_spheres[1] = 1
        for sphere_index in 2:n_spheres
            spheres_of_neighborhood = spheres_of_neighborhoods[sphere_index]
            neighborhood_index = get(seen_spheres_of_neighborhoods, spheres_of_neighborhood, nothing)
            if neighborhood_index === nothing
                neighborhood_index = length(seen_spheres_of_neighborhoods) + 1
                seen_spheres_of_neighborhoods[spheres_of_neighborhood] = neighborhood_index
                unique_spheres_of_neighborhoods[neighborhood_index] = spheres_of_neighborhood
            end
            main_neighborhoods_of_spheres[sphere_index] = neighborhood_index
        end
        resize!(unique_spheres_of_neighborhoods, length(seen_spheres_of_neighborhoods))
    end

    return unique_spheres_of_neighborhoods, main_neighborhoods_of_spheres
end

@logged function compute_membership_matrix(  # untested
    spheres_of_neighborhoods::Vector{Vector{<:Integer}},
    n_spheres::Integer,
)::AbstractMatrix{Bool}
    n_neighborhoods = length(spheres_of_neighborhoods)
    spheres_neighborhoods_is_member = zeros(Bool, n_spheres, n_neighborhoods)
    for (neighborhood, spheres_of_neighborhood) in enumerate(spheres_of_neighborhoods)
        spheres_neighborhoods_is_member[spheres_of_neighborhood, neighborhood] .= true  # NOJET
    end
    return spheres_neighborhoods_is_member
end

@logged function write_data(  # untested
    daf::DafWriter,
    sphere_names::AbstractStringVector,
    neighborhood_names::AbstractStringVector,
    spheres_of_metacells::Vector{UInt32},
    main_neighborhoods_of_spheres::Vector{UInt32},
    spheres_neighborhoods_is_member::AbstractMatrix{Bool},
    spheres_of_neighborhoods::Vector{Vector{<:Integer}},
    spheres_distances::Matrix{<:AbstractFloat},
    overwrite::Bool,
)::Nothing
    check_efficient_action("compute_spheres", Columns, "spheres_distances", spheres_distances)

    if overwrite
        delete_axis!(daf, "sphere"; must_exist = false)
        delete_axis!(daf, "neighborhood"; must_exist = false)
    end
    add_axis!(daf, "sphere", sphere_names)
    add_axis!(daf, "neighborhood", neighborhood_names)

    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells]; overwrite = overwrite)
    set_vector!(daf, "sphere", "neighborhood.main", neighborhood_names[main_neighborhoods_of_spheres])
    set_matrix!(daf, "sphere", "neighborhood", "is_member", spheres_neighborhoods_is_member)

    n_spheres = length(sphere_names)
    nnz = sum([length(spheres_of_neighborhoods[main_neighborhoods_of_spheres[sphere]]) - 1 for sphere in 1:n_spheres])
    empty_sparse_matrix!(daf, "sphere", "sphere", "distance", Float32, nnz) do colptr, rowval, nzval
        colptr[1] = 1
        next_value_index = 1
        for source_sphere in 1:n_spheres
            previous_destination_sphere = 0
            did_see_source_sphere = false
            for destination_sphere in spheres_of_neighborhoods[main_neighborhoods_of_spheres[source_sphere]]
                @assert destination_sphere > previous_destination_sphere
                previous_destination_sphere = destination_sphere
                if destination_sphere != source_sphere
                    rowval[next_value_index] = destination_sphere
                    nzval[next_value_index] = spheres_distances[destination_sphere, source_sphere]
                    next_value_index += 1
                else
                    did_see_source_sphere = true
                end
            end
            @assert did_see_source_sphere
            colptr[source_sphere + 1] = next_value_index
        end
        @assert next_value_index == nnz + 1
        return nothing
    end

    return nothing
end

"""
    function improve_spheres!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = 40,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        confidence::AbstractFloat = 0.9,
        max_sphere_diameter::AbstractFloat = 2.0,
        max_neighborhood_diameter::AbstractFloat = 2.0,
        noisy_gene_fold::AbstractFloat = 1.0,
        min_gene_correlation::AbstractFloat = 0.5,
    )::Nothing

Improve the partition raw metacells into distinct spheres, and spheres into overlapping neighborhoods. The improved
partition ensures no gene is not "too different" between the metacells in each sphere and spheres in each neighborhood;
however, this only applies to genes that have at least `min_gene_correlation` with at least one other gene in the
neighborhood.

 1. Invoke [`identify_marker_genes!`](@ref) and [`identify_lonely_genes!`](@ref) for each neighborhood and store the
    results in `is_marker` and `is_lonely` of the genes for the metacells in sphere for which this is the main
    neighborhood.
 2. Recompute the spheres with `max_deviant_genes = 0` and `overwrite = true`. When computing the distance between each
    two metacells, this will only consider genes that are markers in either metacell and are not lonely in both. The
    result should be "improved" spheres as they will not be affected by genes that aren't part of some multi-gene
    behavior relevant to each sphere's neighborhood.

CONTRACT
"""
@logged @computation IMPROVED_SPHERE_CONTRACT function improve_spheres!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    confidence::AbstractFloat = 0.9,
    max_sphere_diameter::AbstractFloat = 2.0,
    max_neighborhood_diameter::AbstractFloat = 2.0,
    noisy_gene_fold::AbstractFloat = 1.0,
    min_marker_gene_range_fold::AbstractFloat = 2.0,
    min_max_marker_gene_fraction::AbstractFloat = 1e-4,
    max_lonely_gene_correlation::AbstractFloat = 0.5,
    max_lonely_gene_correlated::Integer = 0,
)::Nothing
    identify_neighborhoods_genes!(
        daf;
        gene_fraction_regularization = gene_fraction_regularization,
        noisy_gene_fold = noisy_gene_fold,
        min_marker_gene_range_fold = min_marker_gene_range_fold,
        min_max_marker_gene_fraction = min_max_marker_gene_fraction,
        max_lonely_gene_correlation = max_lonely_gene_correlation,
        max_lonely_gene_correlated = max_lonely_gene_correlated,
    )

    compute_spheres!(
        daf;
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        gene_fraction_regularization = gene_fraction_regularization,
        confidence = confidence,
        max_sphere_diameter = max_sphere_diameter,
        max_neighborhood_diameter = max_neighborhood_diameter,
        noisy_gene_fold = noisy_gene_fold,
        max_deviant_genes = 0,
        overwrite = true,
    )
    return nothing
end

@logged function identify_neighborhoods_genes!(  # untested
    daf::DafReader;
    gene_fraction_regularization::AbstractFloat,
    noisy_gene_fold::AbstractFloat,
    min_marker_gene_range_fold::AbstractFloat,
    min_max_marker_gene_fraction::AbstractFloat,
    max_lonely_gene_correlation::AbstractFloat,
    max_lonely_gene_correlated::Integer,
)::Nothing
    n_genes = axis_length(daf, "gene")
    n_neighborhoods = axis_length(daf, "neighborhood")
    n_metacells = axis_length(daf, "metacell")
    genes_metacells_is_marker = zeros(Bool, n_genes, n_metacells)
    genes_metacells_is_lonely = zeros(Bool, n_genes, n_metacells)
    metacell_index_of_names = daf["/ metacell : name"].dicts[1]

    neighborhood_names = axis_array(daf, "neighborhood")
    @threads for neighborhood_index in 1:n_neighborhoods
        neighborhood_name = neighborhood_names[neighborhood_index]
        neighborhood_daf = chain_writer(
            [
                viewer(
                    daf;
                    name = "view.$(neighborhood_name)",
                    axes = [
                        "metacell" => "/ metacell & sphere => neighborhood.main = $(neighborhood_name)",
                        "gene" => "=",
                    ],
                    data = [("metacell", "gene", "fraction") => "="],
                ),
                MemoryDaf(; name = "identified.$(neighborhood_name)"),
            ];
            name = "chain.$(neighborhood_name)",
        )

        identify_marker_genes!(
            neighborhood_daf;
            gene_fraction_regularization = gene_fraction_regularization,
            min_marker_gene_range_fold = min_marker_gene_range_fold,
            noisy_gene_fold = noisy_gene_fold,
            min_max_marker_gene_fraction = min_max_marker_gene_fraction,
        )
        identify_lonely_genes!(
            neighborhood_daf;
            gene_fraction_regularization = gene_fraction_regularization,
            max_lonely_gene_correlation = max_lonely_gene_correlation,
            max_lonely_gene_correlated = max_lonely_gene_correlated,
        )

        is_marker_of_neighborhood_genes = neighborhood_daf["/ gene : is_marker"]
        is_lonely_of_neighborhood_genes = neighborhood_daf["/ gene : is_lonely"]

        metacell_names_of_neighborhood = axis_array(neighborhood_daf, "metacell")
        metacell_indices_of_neighborhood =
            [metacell_index_of_names[metacell_name] for metacell_name in metacell_names_of_neighborhood]

        genes_metacells_is_marker[:, metacell_indices_of_neighborhood] .= is_marker_of_neighborhood_genes
        genes_metacells_is_lonely[:, metacell_indices_of_neighborhood] .= is_lonely_of_neighborhood_genes
    end

    set_matrix!(daf, "gene", "metacell", "is_marker", genes_metacells_is_marker; overwrite = true)
    set_matrix!(daf, "gene", "metacell", "is_lonely", genes_metacells_is_lonely; overwrite = true)

    return nothing
end

"""
    function compute_improved_spheres!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = 40,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        confidence::AbstractFloat = 0.9,
        max_sphere_diameter::AbstractFloat = 2.0,
        max_neighborhood_diameter::AbstractFloat = 2.0,
        noisy_gene_fold::AbstractFloat = 1.0,
        min_gene_correlation::AbstractFloat = 0.5,
    )::Nothing

Partition raw metacells into distinct spheres, and spheres into overlapping neighborhoods, such that each sphere is
homogeneous in the genes that participate in multi-gene behaviors in its main neighborhood.

This calls [`compute_spheres!`](@ref) and then [`improve_spheres!](@ref). In the 1st call to [`compute_spheres!`](@ref) call, we allow up to `max_deviant_genes_fraction` of the genes to be deviant; in [`improve_spheres!`](@ref), we forbid
all genes from being deviant, as we expect to ignore all "irrelevant" genes as non-markers or lonely.

CONTRACT
"""
@logged @computation IMPROVED_SPHERE_CONTRACT function compute_improved_spheres!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    confidence::AbstractFloat = 0.9,
    max_sphere_diameter::AbstractFloat = 2.0,
    max_neighborhood_diameter::AbstractFloat = 2.0,
    noisy_gene_fold::AbstractFloat = 1.0,
    min_marker_gene_range_fold::AbstractFloat = 2.0,
    min_max_marker_gene_fraction::AbstractFloat = 1e-4,
    max_lonely_gene_correlation::AbstractFloat = 0.5,
    max_lonely_gene_correlated::Integer = 0,
    max_deviant_genes_fraction::AbstractFloat = 0.01,
)::Nothing
    @assert 0.0 <= max_deviant_genes_fraction <= 1.0

    compute_spheres!(
        daf;
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        gene_fraction_regularization = gene_fraction_regularization,
        confidence = confidence,
        max_sphere_diameter = max_sphere_diameter,
        max_neighborhood_diameter = max_neighborhood_diameter,
        noisy_gene_fold = noisy_gene_fold,
        max_deviant_genes = Int(round(axis_length(daf, "gene") * max_deviant_genes_fraction)),
        overwrite = false,
    )
    improve_spheres!(
        daf;
        min_significant_gene_UMIs = min_significant_gene_UMIs,
        gene_fraction_regularization = gene_fraction_regularization,
        confidence = confidence,
        max_sphere_diameter = max_sphere_diameter,
        max_neighborhood_diameter = max_neighborhood_diameter,
        noisy_gene_fold = noisy_gene_fold,
        min_marker_gene_range_fold = min_marker_gene_range_fold,
        min_max_marker_gene_fraction = min_max_marker_gene_fraction,
        max_lonely_gene_correlation = max_lonely_gene_correlation,
        max_lonely_gene_correlated = max_lonely_gene_correlated,
    )
    return nothing
end

end  # module
