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

    @debug "group..."

    sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")
    add_axis!(daf, "sphere", sphere_names)
    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells])

    return nothing
end

@logged function compute_sorted_distance_pairs(  # untested
    daf::DafReader,
    max_sphere_fold_factor::AbstractFloat,
    min_significant_gene_UMIs::Unsigned,
    genes_metacells_log_fraction::AbstractMatrix{<:AbstractFloat},
    genes_metacells_total_UMIs::AbstractMatrix{<:Unsigned},
)::Vector{Tuple{Float32, UInt32, UInt32}}
    n_metacells = axis_length(daf, "metacell")
    n_pairs = div(n_metacells * (n_metacells - 1), 2)

    distance_pairs = Vector{Tuple{Float32, UInt32, UInt32}}(undef, n_pairs)
    next_pair_index = Atomic{Int64}(1)

    @debug "collect close distance pairs..."
    @threads for base_metacell in reverse(2:n_metacells)
        @views genes_log_fraction_of_base_metacell = genes_metacells_log_fraction[:, base_metacell]
        @views genes_log_fraction_of_other_metacells = genes_metacells_log_fraction[:, 1:(base_metacell - 1)]
        @views genes_total_UMIs_of_base_metacell = genes_metacells_total_UMIs[:, base_metacell]
        @views genes_total_UMIs_of_other_metacells = genes_metacells_total_UMIs[:, 1:(base_metacell - 1)]

        fold_factors_of_other_metacells =
            abs.(genes_log_fraction_of_other_metacells .- genes_log_fraction_of_base_metacell) .*
            (genes_total_UMIs_of_other_metacells .+ genes_total_UMIs_of_base_metacell .>= min_significant_gene_UMIs)
        distances_of_other_metacells = maximum(fold_factors_of_other_metacells; dims = 1)  # NOJET
        is_close_of_other_metacells = distances_of_other_metacells .<= max_sphere_fold_factor
        n_close = sum(is_close_of_other_metacells)
        pair_index = atomic_add!(next_pair_index, n_close)
        for (other_metacell, is_close) in enumerate(is_close_of_other_metacells)
            if is_close
                distance = distances_of_other_metacells[other_metacell]
                distance_pairs[pair_index] = (distance, base_metacell, other_metacell)
                pair_index += 1
            end
        end
    end

    n_close = next_pair_index[] - 1
    resize!(distance_pairs, n_close)
    @debug "collected $(n_close) close distance pairs, sort..."
    sort!(distance_pairs)
    return distance_pairs
end

@logged function compute_sphere_of_metacells(  # untested
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
    @views genes_total_UMIs_of_first_metacell = genes_metacells_total_UMIs[:, first_metacell]  # NOJET
    @views genes_total_UMIs_of_second_metacell = genes_metacells_total_UMIs[:, second_metacell]

    return maximum(
        abs.(
            genes_log_fraction_of_first_metacell .- genes_log_fraction_of_second_metacell  #
        )[genes_total_UMIs_of_first_metacell .+ genes_total_UMIs_of_second_metacell .>= min_significant_gene_UMIs],
    )
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

end  # module
