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
using Statistics

"""
    function compute_spheres(
        daf::DafWriter;
        max_sphere_fold_factor::AbstractFloat = 3.0,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        min_significant_gene_UMIs:: Unsigned = 40,
    )::Nothing

Partition raw metacells into distinct spheres, using a variant of union-find:

 1. Compute a distance between each two metacells, defined as the maximal fold factor of any of the genes. This fold
    factor is the absolute value of log base 2 of the fraction of the gene in both metacells, using the
    `gene_fraction_regularization` (by default, `1e-5`). If the sum of the total UMIs of the gene in both metacells is
    less than `min_significant_gene_UMIs` (by default, `40`), we ignore this fold factor as insignificant.
 2. Compute hierarchical clustering of the metacells using these distances and the `complete` linkage; that is, each
    cluster can be seen as a sphere of some diameter where all the metacells are within that sphere.
 3. Use this hierarchical clustering to partition the metacells into as-small as possible spheres, where each sphere
    diameter is at most `max_sphere_fold_factor` (by default, `3.0`, or 8x).

CONTRACT
"""
@logged @computation Contract(
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (RequiredInput, "The genes to consider (typically, only marker genes)."),
        "sphere" => (GuaranteedOutput, "A partition of the metacells into distinct spheres."),
    ],
    data = [
        ("metacell", "gene", "fraction") =>
            (RequiredInput, AbstractFloat, "The fraction of the UMIs of each gene in each metacell."),
        ("metacell", "gene", "total_UMIs") => (
            RequiredInput,
            Unsigned,
            "The total number of UMIs used to estimate the fraction of each gene in each metacell.",
        ),
        ("metacell", "sphere") => (GuaranteedOutput, AbstractString, "The unique sphere each metacell belongs to."),
        ("metacell", "metacell", "distance") => (
            GuaranteedOutput,
            AbstractFloat,
            "The distance (fold factor between the most different gene) between metacells - upper triangular, zero diagonal.",
        ),
        ("sphere", "sphere", "distance") => (
            GuaranteedOutput,
            AbstractFloat,
            """
            The distance (fold factor between the most different gene) between metacells of spheres - upper triangular.
            The diagonal holds the fold factor between the most different gene between two metacells inside each sphere.
            """,
        ),
    ],
) function compute_spheres!(
    daf::DafWriter;
    max_sphere_fold_factor::AbstractFloat = 3.0,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    min_significant_gene_UMIs::Unsigned = Unsigned(40),
)::Nothing
    @assert max_sphere_fold_factor > 0
    @assert gene_fraction_regularization > 0

    @debug "read..."

    genes_metacells_log_fraction =
        daf["/ gene / metacell : fraction % Log base 2 eps $(gene_fraction_regularization)"].array
    genes_metacells_total_UMIs = get_matrix(daf, "gene", "metacell", "total_UMIs").array

    @debug "compute..."

    metacells_distances =
        compute_metacell_distances(min_significant_gene_UMIs, genes_metacells_log_fraction, genes_metacells_total_UMIs)
    clusters = hclust(metacells_distances; linkage = :complete)  # NOJET
    spheres_of_metacells = Vector{UInt32}(cutree(clusters; h = max_sphere_fold_factor))
    metacells_of_spheres = collect_group_members(spheres_of_metacells)
    sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")
    spheres_distances = compute_sphere_distances(metacells_of_spheres, metacells_distances, max_sphere_fold_factor)

    @debug "write..."

    add_axis!(daf, "sphere", sphere_names)
    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells])
    set_matrix!(daf, "metacell", "metacell", "distance", metacells_distances)
    set_matrix!(daf, "sphere", "sphere", "distance", spheres_distances)

    return nothing
end

@logged function compute_metacell_distances(  # untested
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Matrix{Float32}
    check_efficient_action("compute_spheres", Columns, "genes_metacells_log_fraction", genes_metacells_log_fraction)
    check_efficient_action("compute_spheres", Columns, "genes_metacells_total_UMIs", genes_metacells_total_UMIs)
    @assert size(genes_metacells_log_fraction) == size(genes_metacells_total_UMIs)

    n_metacells = size(genes_metacells_log_fraction, 2)
    distances = zeros(Float32, n_metacells, n_metacells)

    @threads for base_metacell in reverse(2:n_metacells)
        @views genes_log_fraction_of_base_metacell = genes_metacells_log_fraction[:, base_metacell]
        @views genes_log_fraction_of_other_metacells = genes_metacells_log_fraction[:, 1:(base_metacell - 1)]
        @views genes_total_UMIs_of_base_metacell = genes_metacells_total_UMIs[:, base_metacell]
        @views genes_total_UMIs_of_other_metacells = genes_metacells_total_UMIs[:, 1:(base_metacell - 1)]

        fold_factors_of_other_metacells =
            abs.(genes_log_fraction_of_other_metacells .- genes_log_fraction_of_base_metacell) .*
            (genes_total_UMIs_of_other_metacells .+ genes_total_UMIs_of_base_metacell .>= min_significant_gene_UMIs)
        distances_of_other_metacells = maximum(fold_factors_of_other_metacells; dims = 1)  # NOJET
        distances[base_metacell, base_metacell] = 0.0
        distances[1:(base_metacell - 1), base_metacell] = distances_of_other_metacells
    end
    distances .+= transpose(distances)
    return distances
end

@logged function compute_sphere_distances(  # untested
    metacells_of_spheres::Vector{Vector{I}},
    metacells_distances::AbstractMatrix{F},
    max_sphere_fold_factor::AbstractFloat,
)::Matrix{F} where {I <: Integer, F <: AbstractFloat}
    check_efficient_action("compute_spheres", Columns, "metacells_distances", metacells_distances)
    n_spheres = length(metacells_of_spheres)  # NOJET
    spheres_distances = zeros(eltype(metacells_distances), n_spheres, n_spheres)
    @threads for base_sphere in reverse(1:n_spheres)
        metacells_of_base_sphere = metacells_of_spheres[base_sphere]
        @views base_sphere_distances_matrix = metacells_distances[:, metacells_of_base_sphere]
        base_sphere_distances_vector = maximum(base_sphere_distances_matrix; dims = 2)
        for other_sphere in 1:base_sphere
            metacells_of_other_sphere = metacells_of_spheres[other_sphere]
            @views other_sphere_distances_vector = base_sphere_distances_vector[metacells_of_other_sphere]
            distance = maximum(other_sphere_distances_vector)
            if other_sphere == base_sphere
                @assert distance <= max_sphere_fold_factor
            else
                @assert distance > max_sphere_fold_factor
            end
            spheres_distances[base_sphere, other_sphere] = maximum(other_sphere_distances_vector)
        end
    end
    return spheres_distances
end

end  # module
