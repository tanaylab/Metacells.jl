"""
Given a set of raw metacells, partition them into spheres such that all metacells in the same sphere are within some
(fold factor) radius of each other. The centroids of these spheres can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states.
"""
module Spheres

export compute_spheres!

using Base.Iterators
using Base.Threads
using Clustering
using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using Distributions
using Statistics

"""
    function compute_spheres(
        daf::DafWriter;
        min_significant_gene_UMIs::Unsigned = Unsigned(40),
        gene_fraction_regularization::AbstractFloat = 1e-5,
        confidence::AbstractFloat = 0.9,
        max_sphere_diameter::AbstractFloat = 2.0,
        max_neighborhood_diameter::AbstractFloat = 2.0,
        noisy_diameter::AbstractFloat = 1.0,
    )::Nothing

Partition raw metacells into distinct spheres, and spheres into overlapping neighborhoods.

 1. Compute a distance between each two metacells, defined as the maximal fold factor of any of the genes. This fold
    factor is the absolute value of log base 2 of the fraction of the gene in both metacells, using the
    `gene_fraction_regularization` (by default, `1e-5`). Since the fraction of the gene is a random variable, we
    decrease the high fraction and increase the low fraction by a factor based on the `confidence` of the test (by
    default, 0.9), assuming a multinomial distribution. In addition, if the sum of the total UMIs of the gene in both
    metacells is less than `min_significant_gene_UMIs` (by default, `40`), we ignore this fold factor as insignificant.
    For `is_noisy` genes, we reduce the distance by `noisy_diameter` to allow for their bursty nature.
 2. Compute hierarchical clustering of the metacells using these distances and the `complete` linkage; that is, each
    cluster can be seen as a sphere of some diameter where all the metacells are within that sphere.
 3. Use this hierarchical clustering to partition the metacells into as-small as possible spheres, where each sphere
    diameter is at most `max_sphere_diameter` (by default, `2.0`, or 4x).
 4. Compute for each sphere the set of other spheres with a maximal distance of an additional
    `max_neighborhood_diameter` (by default, `2.0`, for a total neighborhood diameter of 16x). For each sphere, the main
    neighborhood is the one defined by that sphere; however, each sphere may belong to other neighborhoods as well.
 5. The distances to the other spheres in the main neighborhood of each sphere form a natural graph representation of
    the structure of the manifold of the distinct biological cell states captured by the metacells.

CONTRACT
"""
@logged @computation Contract(
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
        ("gene", "is_noisy") =>
            (OptionalInput, Bool, "A mask of noisy genes to be given additional diameter when grouped."),
        ("metacell", "sphere") => (GuaranteedOutput, AbstractString, "The unique sphere each metacell belongs to."),
        ("sphere", "main_neighborhood") => (GuaranteedOutput, AbstractString, "The main neighborhood of each sphere."),
        ("sphere", "neighborhood", "is_member") =>
            (GuaranteedOutput, Bool, "Membership matrix for spheres and neighborhoods."),
        ("sphere", "sphere", "distance") =>
            (GuaranteedOutput, Float32, "For each sphere, the distances of the spheres in its main neighborhood."),
    ],
) function compute_spheres!(
    daf::DafWriter;
    min_significant_gene_UMIs::Unsigned = Unsigned(40),
    gene_fraction_regularization::AbstractFloat = 1e-5,
    confidence::AbstractFloat = 0.9,
    max_sphere_diameter::AbstractFloat = 2.0,
    max_neighborhood_diameter::AbstractFloat = 2.0,
    noisy_diameter::AbstractFloat = 1.0,
    metacell_distances::Maybe{Matrix{Float32}},
)::Matrix{Float32}
    @assert gene_fraction_regularization > 0
    @assert confidence > 0.0
    @assert confidence < 1.0
    @assert max_sphere_diameter > 0
    @assert max_neighborhood_diameter >= 0.0

    genes_metacells_total_UMIs, genes_metacells_fraction, metacells_total_UMIs, genes_is_noisy = read_data(daf)

    genes_metacells_log_decreased_fraction, genes_metacells_log_increased_fraction =
        prepare_metacells_data(gene_fraction_regularization, confidence, genes_metacells_fraction, metacells_total_UMIs)

    if metacell_distances === nothing
        metacell_distances = compute_metacell_distances(  # NOJET
            min_significant_gene_UMIs,
            genes_metacells_total_UMIs,
            genes_metacells_log_decreased_fraction,
            genes_metacells_log_increased_fraction,
            noisy_diameter,
            genes_is_noisy,
        )
    end

    spheres_of_metacells = compute_spheres_of_metacells(metacell_distances, max_sphere_diameter)
    metacells_of_spheres = collect_group_members(spheres_of_metacells)
    sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")

    spheres_distances = compute_sphere_distances(metacells_of_spheres, metacell_distances, max_sphere_diameter)
    spheres_of_neighborhoods, main_neighborhoods_of_spheres =
        collect_spheres_of_neighborhoods(spheres_distances, max_sphere_diameter, max_neighborhood_diameter)
    neighborhood_names = group_names(daf, "metacell", spheres_of_neighborhoods; prefix = "N")
    spheres_neighborhoods_is_member = compute_membership_matrix(spheres_of_neighborhoods, length(sphere_names))

    write_metacells_data(
        daf,
        sphere_names,
        neighborhood_names,
        spheres_of_metacells,
        main_neighborhoods_of_spheres,
        spheres_neighborhoods_is_member,
        spheres_of_neighborhoods,
        spheres_distances,
    )

    return metacell_distances
end

@logged function read_data(  # untested
    daf::DafReader,
)::Tuple{AbstractMatrix{<:Unsigned}, AbstractMatrix{<:AbstractFloat}, AbstractVector{<:Unsigned}, AbstractVector{Bool}}
    genes_metacells_total_UMIs = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    genes_metacells_fraction = get_matrix(daf, "gene", "metacell", "fraction").array
    metacells_total_UMIs = get_vector(daf, "metacell", "total_UMIs").array
    genes_is_noisy = get_vector(daf, "gene", "is_noisy"; default = false).array
    return genes_metacells_total_UMIs, genes_metacells_fraction, metacells_total_UMIs, genes_is_noisy
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

@logged function compute_metacell_distances(  # untested
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
    genes_metacells_log_decreased_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_log_increased_fraction::AbstractMatrix{<:AbstractFloat},
    noisy_diameter::AbstractFloat,
    genes_is_noisy::AbstractVector{Bool},
)::Matrix{Float32}
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

    n_metacells = size(genes_metacells_total_UMIs, 2)
    metacell_distances = zeros(Float32, n_metacells, n_metacells)

    @threads for base_metacell in reverse(2:n_metacells)
        @views genes_total_UMIs_of_base_metacell = genes_metacells_total_UMIs[:, base_metacell]
        @views genes_total_UMIs_of_other_metacell = genes_metacells_total_UMIs[:, 1:(base_metacell - 1)]
        @views genes_log_decreased_fraction_of_base_metacell = genes_metacells_log_decreased_fraction[:, base_metacell]
        @views genes_log_increased_fraction_of_base_metacell = genes_metacells_log_increased_fraction[:, base_metacell]
        @views genes_log_decreased_fraction_of_other_metacells =
            genes_metacells_log_decreased_fraction[:, 1:(base_metacell - 1)]
        @views genes_log_increased_fraction_of_other_metacells =
            genes_metacells_log_increased_fraction[:, 1:(base_metacell - 1)]

        significant_fold_factors_of_other_metacells =
            gene_distance.(
                min_significant_gene_UMIs,
                genes_total_UMIs_of_base_metacell,
                genes_total_UMIs_of_other_metacell,
                genes_log_decreased_fraction_of_base_metacell,
                genes_log_increased_fraction_of_base_metacell,
                genes_log_decreased_fraction_of_other_metacells,
                genes_log_increased_fraction_of_other_metacells,
                noisy_diameter,
                genes_is_noisy,
            )

        metacell_distances[1:(base_metacell - 1), base_metacell] .=
            transpose(maximum(significant_fold_factors_of_other_metacells; dims = 1))
    end

    metacell_distances .+= transpose(metacell_distances)
    return metacell_distances
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Unsigned,
    gene_total_UMIs_of_base_metacell::Unsigned,
    gene_total_UMIs_of_other_metacell::Unsigned,
    gene_log_decreased_fraction_of_base_metacell::AbstractFloat,
    gene_log_increased_fraction_of_base_metacell::AbstractFloat,
    gene_log_decreased_fraction_of_other_metacell::AbstractFloat,
    gene_log_increased_fraction_of_other_metacell::AbstractFloat,
    noisy_diameter::AbstractFloat,
    gene_is_noisy::Bool,
)::AbstractFloat
    gene_total_UMIs = gene_total_UMIs_of_base_metacell + gene_total_UMIs_of_other_metacell
    is_significant = gene_total_UMIs >= min_significant_gene_UMIs
    is_base_low = gene_log_increased_fraction_of_base_metacell < gene_log_increased_fraction_of_other_metacell
    gene_log_low_fraction =
        is_base_low * gene_log_increased_fraction_of_base_metacell +
        !is_base_low * gene_log_increased_fraction_of_other_metacell
    gene_log_high_fraction =
        !is_base_low * gene_log_decreased_fraction_of_base_metacell +
        is_base_low * gene_log_decreased_fraction_of_other_metacell
    return is_significant * max.(gene_log_high_fraction - gene_log_low_fraction - gene_is_noisy * noisy_diameter, 0.0)
end

@logged function compute_sphere_distances(  # untested
    metacells_of_spheres::Vector{Vector{I}},
    metacells_distances::AbstractMatrix{F},
    max_sphere_diameter::AbstractFloat,
)::Matrix{F} where {I <: Integer, F <: AbstractFloat}
    check_efficient_action("compute_spheres", Columns, "metacells_distances", metacells_distances)

    n_spheres = length(metacells_of_spheres)
    spheres_distances = zeros(eltype(metacells_distances), n_spheres, n_spheres)

    @threads for base_sphere in reverse(1:n_spheres)
        metacells_of_base_sphere = metacells_of_spheres[base_sphere]
        @views base_sphere_distances_matrix = metacells_distances[:, metacells_of_base_sphere]
        base_sphere_distances_vector = maximum(base_sphere_distances_matrix; dims = 2)
        for other_sphere in 1:base_sphere
            metacells_of_other_sphere = metacells_of_spheres[other_sphere]
            @views other_sphere_distances_vector = base_sphere_distances_vector[metacells_of_other_sphere]
            distance = maximum(other_sphere_distances_vector)
            if base_sphere == other_sphere
                @assert distance <= max_sphere_diameter
            else
                @assert distance > max_sphere_diameter
            end
            spheres_distances[base_sphere, other_sphere] = distance
            spheres_distances[other_sphere, base_sphere] = distance
        end
    end

    return spheres_distances
end

@logged function compute_spheres_of_metacells(  # untested
    metacell_distances::Matrix{Float32},
    max_sphere_diameter::AbstractFloat,
)::Vector{UInt32}
    check_efficient_action("compute_spheres", Columns, "metacell_distances", metacell_distances)
    clusters = hclust(metacell_distances; linkage = :complete)  # NOJET
    spheres_of_metacells = Vector{UInt32}(cutree(clusters; h = max_sphere_diameter))
    return spheres_of_metacells
end

@logged function collect_spheres_of_neighborhoods(  # untested
    spheres_distances::AbstractMatrix{F},
    max_sphere_diameter::AbstractFloat,
    max_neighborhood_diameter::AbstractFloat,
)::Tuple{Vector{Vector{<:Integer}}, Vector{UInt32}} where {F <: AbstractFloat}
    n_spheres = size(spheres_distances, 1)
    spheres_of_neighborhoods = [
        findall(spheres_distances[:, sphere] .<= max_sphere_diameter + max_neighborhood_diameter) for
        sphere in 1:n_spheres
    ]
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

@logged function write_metacells_data(  # untested
    daf::DafWriter,
    sphere_names::AbstractStringVector,
    neighborhood_names::AbstractStringVector,
    spheres_of_metacells::Vector{UInt32},
    main_neighborhoods_of_spheres::Vector{UInt32},
    spheres_neighborhoods_is_member::AbstractMatrix{Bool},
    spheres_of_neighborhoods::Vector{Vector{<:Integer}},
    spheres_distances::Matrix{<:AbstractFloat},
)::Nothing
    check_efficient_action("compute_spheres", Columns, "spheres_distances", spheres_distances)

    add_axis!(daf, "sphere", sphere_names)
    add_axis!(daf, "neighborhood", neighborhood_names)
    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells])
    set_vector!(daf, "sphere", "main_neighborhood", neighborhood_names[main_neighborhoods_of_spheres])
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

end  # module
