"""
Given a set of raw metacells, partition them into spheres such that all metacells in the same sphere are within some
(fold factor) radius of each other. The centroids of these spheres can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states.
"""
module Spheres

export compute_spheres!

using Base.Iterators
using Base.Threads
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

 2. Initially place each metacell in its own sphere.
 3. Pass on all metacell pairs, ordered by increasing distance. If the metacells belong to different spheres, merge
    both spheres into a single one if the maximal distance between the metacells of the combined sphere is no more than
    `max_sphere_fold_factor` (by default, `3.0`, or 8x).
 4. Compute the centroid (geomean) of each gene's fraction of each sphere, using the same
    `gene_fraction_regularization` as above.

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
        ("sphere", "gene", "fraction") =>
            (GuaranteedOutput, AbstractFloat, "The centroid fraction of each gene in each sphere."),
        ("sphere", "gene", "total_UMIs") => (
            GuaranteedOutput,
            AbstractFloat,
            "The total number of UMIs used to estimate the fraction of each gene in each sphere.",
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

    genes_metacells_log_fraction = daf["/ gene / metacell : fraction % Log base 2 eps $(gene_fraction_regularization)"]
    genes_metacells_total_UMIs = get_matrix(daf, "gene", "metacell", "total_UMIs")

    check_efficient_action("compute_spheres", Columns, "genes_metacells_log_fraction", genes_metacells_log_fraction)
    check_efficient_action("compute_spheres", Columns, "genes_metacells_total_UMIs", genes_metacells_total_UMIs)

    sorted_distance_pairs = compute_sorted_distance_pairs(
        daf,
        max_sphere_fold_factor,
        min_significant_gene_UMIs,
        genes_metacells_log_fraction,
        genes_metacells_total_UMIs,
    )

    spheres_of_metacells = compute_sphere_of_metacells(
        sorted_distance_pairs,
        max_sphere_fold_factor,
        min_significant_gene_UMIs,
        genes_metacells_log_fraction,
        genes_metacells_total_UMIs,
    )

    metacells_of_spheres = collect_group_members(spheres_of_metacells)

    sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")
    add_axis!(daf, "sphere", sphere_names)
    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells])

    genes_spheres_fraction, genes_spheres_total_UMIs = compute_spheres_data(
        gene_fraction_regularization,
        metacells_of_spheres,
        genes_metacells_log_fraction,
        genes_metacells_total_UMIs,
    )
    set_matrix!(daf, "gene", "sphere", "fraction", genes_spheres_fraction)
    set_matrix!(daf, "gene", "sphere", "total_UMIs", genes_spheres_total_UMIs)

    return nothing
end

function compute_sorted_distance_pairs(  # untested
    daf::DafReader,
    max_sphere_fold_factor::AbstractFloat,
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Vector{Tuple{Float32, UInt32, UInt32}}
    distance_pairs = Vector{Tuple{Float32, UInt32, UInt32}}()

    n_metacells = axis_length(daf, "metacell")
    for first_metacell in 1:n_metacells
        for second_metacell in 1:(first_metacell - 1)
            distance = compute_distance(
                first_metacell,
                second_metacell,
                min_significant_gene_UMIs,
                genes_metacells_log_fraction,
                genes_metacells_total_UMIs,
            )
            if distance <= max_sphere_fold_factor
                push!(distance_pairs, (distance, first_metacell, second_metacell))
            end
        end
    end

    sort!(distance_pairs)
    return distance_pairs
end

function compute_sphere_of_metacells(  # untested
    sorted_distance_pairs::Vector{Tuple{Float32, UInt32, UInt32}},
    max_sphere_fold_factor::AbstractFloat,
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Vector{UInt32}
    n_metacells = size(genes_metacells_log_fraction, 2)
    spheres_of_metacells = collect(UInt32, 1:n_metacells)
    distinct_spheres_set = Set{Tuple{UInt32, UInt32}}()
    metacells_of_spheres = Vector{Vector{UInt32}}()
    for metacell in 1:n_metacells
        push!(metacells_of_spheres, UInt32[metacell])
    end

    for (_, first_metacell, second_metacell) in sorted_distance_pairs
        first_sphere = spheres_of_metacells[first_metacell]
        second_sphere = spheres_of_metacells[second_metacell]
        if first_sphere != second_sphere
            low_sphere = min(first_sphere, second_sphere)
            high_sphere = max(first_sphere, second_sphere)
            spheres_key = (low_sphere, high_sphere)
            if spheres_key in distinct_spheres_set
                continue
            end
            if !can_merge(
                low_sphere,
                high_sphere,
                metacells_of_spheres,
                max_sphere_fold_factor,
                min_significant_gene_UMIs,
                genes_metacells_log_fraction,
                genes_metacells_total_UMIs,
            )
                push!(distinct_spheres_set, spheres_key)
            end
            merge_spheres(low_sphere, high_sphere, spheres_of_metacells, metacells_of_spheres)
        end
    end

    compact_groups!(spheres_of_metacells)
    return spheres_of_metacells
end

function can_merge(  # untested
    low_sphere::Integer,
    high_sphere::Integer,
    metacells_of_spheres::Vector{Vector{UInt32}},
    max_sphere_fold_factor::AbstractFloat,
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Bool
    low_metacells = metacells_of_spheres[low_sphere]
    high_metacells = metacells_of_spheres[high_sphere]
    for (low_metacell, high_metacell) in product(low_metacells, high_metacells)
        if compute_distance(
            low_metacell,
            high_metacell,
            min_significant_gene_UMIs,
            genes_metacells_log_fraction,
            genes_metacells_total_UMIs,
        ) > max_sphere_fold_factor
            return false
        end
    end

    return true
end

function compute_distance(  # untested
    first_metacell::Integer,
    second_metacell::Integer,
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Float32
    @views genes_log_fraction_of_first_metacell = genes_metacells_log_fraction[:, first_metacell]  # NOJET
    @views genes_log_fraction_of_second_metacell = genes_metacells_log_fraction[:, second_metacell]
    fold_factors = abs.(genes_log_fraction_of_first_metacell .- genes_log_fraction_of_second_metacell)

    @views genes_total_UMIs_of_first_metacell = genes_metacells_total_UMIs[:, first_metacell]  # NOJET
    @views genes_total_UMIs_of_second_metacell = genes_metacells_total_UMIs[:, second_metacell]
    significant_fold_factors_mask =
        (genes_total_UMIs_of_first_metacell .+ genes_total_UMIs_of_second_metacell) .>= min_significant_gene_UMIs

    @views distance = max(fold_factors[significant_fold_factors_mask])
    return distance
end

function merge_spheres(  # untested
    low_sphere::Integer,
    high_sphere::Integer,
    spheres_of_metacells::Vector{UInt32},
    metacells_of_spheres::Vector{Vector{UInt32}},
)::Nothing
    low_metacells = metacells_of_spheres[low_sphere]
    high_metacells = metacells_of_spheres[high_sphere]
    spheres_of_metacells[high_metacells] .= low_sphere
    append!(low_metacells, high_metacells)
    empty!(high_metacells)
    return nothing
end

function compute_spheres_data(  # untested
    gene_fraction_regularization::AbstractFloat,
    metacells_of_spheres::AbstractVector{<:AbstractVector{<:Integer}},
    genes_metacells_log_fraction::AbstractMatrix{F},
    genes_metacells_total_UMIs::AbstractMatrix{T},
)::Tuple{Matrix{F}, Matrix{T}} where {F <: AbstractFloat, T <: Unsigned}
    n_genes = size(genes_metacells_log_fraction, 1)
    n_spheres = length(metacells_of_spheres)

    genes_spheres_fraction = Matrix{F}(undef, (n_genes, n_spheres))
    genes_spheres_total_UMIs = Matrix{T}(undef, (n_genes, n_spheres))

    check_efficient_action("compute_spheres_data", Columns, "genes_spheres_fraction", genes_spheres_fraction)
    check_efficient_action("compute_spheres_data", Columns, "genes_spheres_total_UMIs", genes_spheres_total_UMIs)

    @threads for (sphere, metacells_of_sphere) in enumerate(metacells_of_spheres)
        @views genes_sphere_metacells_total_UMIs = genes_metacells_total_UMIs[:, metacells_of_sphere]
        @views sphere_total_UMIs_of_genes = genes_spheres_total_UMIs[:, sphere]
        sphere_total_UMIs_of_genes .= sum(genes_sphere_metacells_total_UMIs; dims = 2)

        @views genes_sphere_metacells_log_fraction = genes_metacells_log_fraction[:, metacells_of_sphere]
        @views sphere_fraction_of_genes = genes_spheres_fraction[:, sphere]
        sphere_fraction_of_genes .=  # NOJET
            2 .^ mean(genes_sphere_metacells_log_fraction; dims = 2) - gene_fraction_regularization
    end

    return (genes_spheres_fraction, genes_spheres_total_UMIs)
end

end  # module
