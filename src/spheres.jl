"""
Given a set of raw metacells, partition them into spheres such that all metacells in the same sphere are within some
(fold factor) radius of each other. The centroids of these spheres can serve as a representation of the cell state
manifold which is less sensitive to oversampling of common cell states. Group these spheres in overlapping neighborhoods
of "similar" spheres for further analysis.
"""
module Spheres

export compute_spheres!

using ..IdentifyGenes

using Base.Iterators
using Base.Threads
using Clustering
using Daf
using Daf.GenericLogging
using Daf.GenericTypes
using Distributions
using SparseArrays
using Statistics

"""
    function compute_spheres!(
        daf::DafWriter;
        min_significant_gene_UMIs::Integer = 40,
        gene_fraction_regularization::AbstractFloat = 1e-5,
        confidence::AbstractFloat = 0.9,
        max_sphere_diameter::AbstractFloat = 2.0,
        max_neighborhood_diameter::AbstractFloat = 2.0,
        noisy_gene_fold::AbstractFloat = 1.0,
        min_gene_correlation::AbstractFloat = 0.5,
        max_deviant_genes_fraction::AbstractFloat = 0.01,
        overwrite::Bool = false,
    )::Nothing

Partition raw metacells into distinct spheres, and spheres into overlapping neighborhoods.

 1. Initial spheres and neighborhoods are computed in a first round, and then refined in a series of followup rounds.
 2. In each round, we compute a distance between each two metacells. This is based on the fold factor between the
    expression level of each (relevant) gene in the metacells. The fold factor is the absolute value of the difference in
    the log (base 2) of the fraction of the gene in the metacells. This log is computed with the
    `gene_fraction_regularization` (by default, `1e-5`). Since the fraction of the gene is a random variable, we decrease
    the high fraction and increase the low fraction by a factor based on the `confidence` of the test (by default, 0.9),
    assuming a multinomial distribution. In addition, if the sum of the total UMIs of the gene in both metacells is less
    than `min_significant_gene_UMIs` (by default, `40`), we ignore this fold factor as insignificant. Finally, for noisy
    genes, we reduce the fold factor by `noisy_gene_fold`. In the first round, we simply count the number of genes whose
    fold factor is more than `max_sphere_diameter` (for computing spheres) and `max_sphere_diameter` +
    `max_neighborhood_diameter` (for computing neighborhoods). In the followup rounds, we use the maximal gene fold, for
    genes that are correlated in the vicinity of the metacells (see below).
 3. We use hierarchical clustering to partition the metacells to distinct spheres, such that the maximal distance between
    any metacells in the sphere is bounded. In the first round, this bound is the `max_deviant_genes_fraction` out of the
    total number of genes. In the followup rounds, this is the `max_sphere_diameter`.
 4. For each sphere, we compute a main neighborhood of other spheres such that the maximal distance between any metacells
    in the neighborhood is bounded. In the first round, this bound is again the maximal number of deviant genes (this
    time, using the increased fold distance computed above). In the followup rounds, this is the `max_sphere_diameter`
    plus the `max_neighborhood_diameter`. These neighborhoods may overlap. The main neighborhoods of different spheres
    may even be identical.
 5. For each sphere, we compute the set of genes which have at least the `min_gene_correlation` with some other gene(s)
    in its main neighborhood. We restrict the correlated set of genes of each metacell to be the intersection of this set
    with the set from its sphere in the previous round.
 6. If the new sets of correlated genes are identical to the previous round, we are done. Otherwise we repeat the round,
    using the more restricted sets of correlated genes.

If `overwrite` is set, the results will replace any previously computed spheres and neighborhoods.

CONTRACT
"""
@logged @computation Contract(;
    axes = [
        "metacell" => (RequiredInput, "The metacells to group into neighborhoods."),
        "gene" => (
            RequiredInput,
            "The genes to consider (typically, only non-lateral globally marker and correlated genes).",
        ),
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
        ("gene", "sphere", "is_correlated") =>
            (GuaranteedOutput, Bool, "Which genes are correlated in the vicinity of each sphere."),
        ("sphere", "neighborhood.main") => (GuaranteedOutput, AbstractString, "The main neighborhood of each sphere."),
        ("sphere", "neighborhood", "is_member") =>
            (GuaranteedOutput, Bool, "Membership matrix for spheres and neighborhoods."),
        ("sphere", "sphere", "distance") =>
            (GuaranteedOutput, Float32, "For each sphere, the distances of the spheres in its main neighborhood."),
    ],
) function compute_spheres!(  # untested
    daf::DafWriter;
    min_significant_gene_UMIs::Integer = 40,
    gene_fraction_regularization::AbstractFloat = 1e-5,
    confidence::AbstractFloat = 0.9,
    max_sphere_diameter::AbstractFloat = 2.0,
    max_neighborhood_diameter::AbstractFloat = 2.0,
    noisy_gene_fold::AbstractFloat = 1.0,
    min_gene_correlation::AbstractFloat = 0.5,
    max_deviant_genes_fraction::AbstractFloat = 0.01,
    overwrite::Bool = false,
)::Nothing
    @assert min_significant_gene_UMIs >= 0
    @assert gene_fraction_regularization > 0
    @assert 0 < confidence < 1
    @assert max_sphere_diameter > 0
    @assert max_neighborhood_diameter >= 0
    @assert noisy_gene_fold >= 0.0
    @assert 0 < min_gene_correlation < 1
    @assert 0 < max_deviant_genes_fraction < 1

    total_UMIs_of_genes_in_metacells, fraction_of_genes_in_metacells, total_UMIs_of_metacells, is_noisy_of_genes =
        read_data(daf)

    log_decreased_fraction_of_genes_in_metacells, log_increased_fraction_of_genes_in_metacells =
        compute_confidence_log_fraction_of_genes_in_metacells(;
            gene_fraction_regularization = gene_fraction_regularization,
            confidence = confidence,
            fraction_of_genes_in_metacells = fraction_of_genes_in_metacells,
            total_UMIs_of_metacells = total_UMIs_of_metacells,
        )

    distance_for_spheres_between_metacells, distance_for_neighborhoods_between_metacells =
        compute_initial_metacells_distances(;
            total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
            log_decreased_fraction_of_genes_in_metacells = log_decreased_fraction_of_genes_in_metacells,
            log_increased_fraction_of_genes_in_metacells = log_increased_fraction_of_genes_in_metacells,
            is_noisy_of_genes = is_noisy_of_genes,
            min_significant_gene_UMIs = min_significant_gene_UMIs,
            max_sphere_diameter = max_sphere_diameter,
            max_neighborhood_diameter = max_neighborhood_diameter,
            noisy_gene_fold = noisy_gene_fold,
        )

    n_genes = axis_length(daf, "gene")
    n_metacells = axis_length(daf, "metacell")
    is_correlated_of_genes_in_metacells = fill(true, n_genes, n_metacells)
    previous_n_correlated_genes_in_metacells = n_genes * n_metacells
    max_deviant_genes = max_deviant_genes_fraction * n_genes
    @assert max_deviant_genes > 0
    spheres_threshold = max_deviant_genes
    neighborhood_threshold = max_deviant_genes

    spheres_of_metacells = nothing
    sphere_names = nothing
    distances_between_spheres = nothing
    main_neighborhoods_of_spheres = nothing
    neighborhood_names = nothing
    spheres_of_neighborhoods = nothing
    is_correlated_of_genes_in_spheres = nothing

    round_index = 0
    while true
        round_index += 1

        spheres_of_metacells = cluster_metacells_to_spheres(;
            distances_between_metacells = distance_for_spheres_between_metacells,
            threshold = spheres_threshold,
        )
        metacells_of_spheres = collect_group_members(spheres_of_metacells)
        sphere_names = group_names(daf, "metacell", metacells_of_spheres; prefix = "S")
        @debug "round: $(round_index) n_spheres: $(length(sphere_names))"

        distances_between_spheres = compute_sphere_distances(;
            distances_between_metacells = distance_for_neighborhoods_between_metacells,
            metacells_of_spheres = metacells_of_spheres,
        )

        spheres_of_neighborhoods, main_neighborhoods_of_spheres = collect_spheres_of_neighborhoods(;
            distances_between_spheres = distances_between_spheres,
            threshold = neighborhood_threshold,
        )

        neighborhood_names = group_names(daf, "metacell", spheres_of_neighborhoods; prefix = "N")
        @debug "round: $(round_index) n_neighborhoods: $(length(neighborhood_names))"
        @assert length(neighborhood_names) <= length(sphere_names)

        intermediate_daf = chain_writer([daf, MemoryDaf(; name = "intermediate")])
        add_axis!(intermediate_daf, "sphere", sphere_names)
        add_axis!(intermediate_daf, "neighborhood", neighborhood_names)
        set_vector!(intermediate_daf, "metacell", "sphere", sphere_names[spheres_of_metacells]; overwrite = overwrite)
        set_vector!(intermediate_daf, "sphere", "neighborhood.main", neighborhood_names[main_neighborhoods_of_spheres])

        is_correlated_of_genes_in_spheres = identify_spheres_correlated_genes(
            intermediate_daf;
            gene_fraction_regularization = gene_fraction_regularization,
            min_gene_correlation = min_gene_correlation,
            is_correlated_of_genes_in_metacells = is_correlated_of_genes_in_metacells,
            metacells_of_spheres = metacells_of_spheres,
        )
        next_n_correlated_genes_in_metacells = sum(is_correlated_of_genes_in_metacells)
        @debug "round: $(round_index) is_correlated_genes_in_metacells: $(depict(is_correlated_of_genes_in_metacells))"
        if next_n_correlated_genes_in_metacells == previous_n_correlated_genes_in_metacells
            break
        end
        previous_n_correlated_genes_in_metacells = next_n_correlated_genes_in_metacells

        distance_for_spheres_between_metacells =
            distance_for_neighborhoods_between_metacells = compute_followup_metacells_distances(;
                total_UMIs_of_genes_in_metacells = total_UMIs_of_genes_in_metacells,
                log_decreased_fraction_of_genes_in_metacells = log_decreased_fraction_of_genes_in_metacells,
                log_increased_fraction_of_genes_in_metacells = log_increased_fraction_of_genes_in_metacells,
                metacells_of_spheres = metacells_of_spheres,
                spheres_of_neighborhoods = spheres_of_neighborhoods,
                main_neighborhoods_of_spheres = main_neighborhoods_of_spheres,
                is_correlated_of_genes_in_spheres = is_correlated_of_genes_in_spheres,
                is_noisy_of_genes = is_noisy_of_genes,
                min_significant_gene_UMIs = min_significant_gene_UMIs,
                noisy_gene_fold = noisy_gene_fold,
            )

        spheres_threshold = max_sphere_diameter
        neighborhood_threshold = max_sphere_diameter + max_neighborhood_diameter
    end

    is_member_of_spheres_in_neighborhoods = compute_membership_matrix(spheres_of_neighborhoods, length(sphere_names))

    write_data(
        daf;
        sphere_names = sphere_names,
        neighborhood_names = neighborhood_names,
        spheres_of_metacells = spheres_of_metacells,
        main_neighborhoods_of_spheres = main_neighborhoods_of_spheres,
        is_member_of_spheres_in_neighborhoods = is_member_of_spheres_in_neighborhoods,
        is_correlated_of_genes_in_spheres = is_correlated_of_genes_in_spheres,
        distances_between_spheres = distances_between_spheres,
        overwrite = overwrite,
    )

    return nothing
end

@logged function read_data(  # untested
    daf::DafReader,
)::Tuple{AbstractMatrix{<:Unsigned}, AbstractMatrix{<:AbstractFloat}, AbstractVector{<:Unsigned}, AbstractVector{Bool}}
    total_UMIs_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "total_UMIs").array
    fraction_of_genes_in_metacells = get_matrix(daf, "gene", "metacell", "fraction").array
    total_UMIs_of_metacells = get_vector(daf, "metacell", "total_UMIs").array
    is_noisy_of_genes = get_vector(daf, "gene", "is_noisy"; default = false).array

    return total_UMIs_of_genes_in_metacells, fraction_of_genes_in_metacells, total_UMIs_of_metacells, is_noisy_of_genes
end

@logged function write_data(  # untested
    daf::DafWriter;
    sphere_names::AbstractVector{<:AbstractString},
    neighborhood_names::AbstractVector{<:AbstractString},
    spheres_of_metacells::Vector{UInt32},
    main_neighborhoods_of_spheres::Vector{UInt32},
    is_member_of_spheres_in_neighborhoods::AbstractMatrix{Bool},
    is_correlated_of_genes_in_spheres::AbstractMatrix{Bool},
    distances_between_spheres::Matrix{<:AbstractFloat},
    overwrite::Bool,
)::Nothing
    check_efficient_action("compute_spheres", Columns, "distances_between_spheres", distances_between_spheres)

    if overwrite
        delete_axis!(daf, "sphere"; must_exist = false)
        delete_axis!(daf, "neighborhood"; must_exist = false)
    end
    add_axis!(daf, "sphere", sphere_names)
    add_axis!(daf, "neighborhood", neighborhood_names)

    set_vector!(daf, "metacell", "sphere", sphere_names[spheres_of_metacells]; overwrite = overwrite)
    set_vector!(daf, "sphere", "neighborhood.main", neighborhood_names[main_neighborhoods_of_spheres])
    set_matrix!(daf, "sphere", "neighborhood", "is_member", SparseMatrixCSC(is_member_of_spheres_in_neighborhoods))
    set_matrix!(daf, "gene", "sphere", "is_correlated", SparseMatrixCSC(is_correlated_of_genes_in_spheres))
    set_matrix!(daf, "sphere", "sphere", "distance", distances_between_spheres)

    return nothing
end

@logged function compute_confidence_log_fraction_of_genes_in_metacells(;  # untested
    gene_fraction_regularization::AbstractFloat,
    confidence::AbstractFloat,
    fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    total_UMIs_of_metacells::AbstractVector{<:Unsigned},
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action("compute_spheres", Columns, "fraction_of_genes_in_metacells", fraction_of_genes_in_metacells)

    confidence_stdevs = quantile(Normal(), confidence)

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
    is_noisy_of_genes::AbstractVector{Bool},
    min_significant_gene_UMIs::Integer,
    max_sphere_diameter::AbstractFloat,
    max_neighborhood_diameter::AbstractFloat,
    noisy_gene_fold::AbstractFloat,
)::Tuple{Matrix{Float32}, Matrix{Float32}}
    check_efficient_action(
        "compute_spheres",
        Columns,
        "total_UMIs_of_genes_in_metacells",
        total_UMIs_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_spheres",
        Columns,
        "log_decreased_fraction_of_genes_in_metacells",
        log_decreased_fraction_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_spheres",
        Columns,
        "log_increased_fraction_of_genes_in_metacells",
        log_increased_fraction_of_genes_in_metacells,
    )
    @assert size(log_decreased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)
    @assert size(log_increased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)

    n_genes = size(total_UMIs_of_genes_in_metacells, 1)
    n_metacells = size(total_UMIs_of_genes_in_metacells, 2)

    distance_for_spheres_between_metacells = zeros(Float32, n_metacells, n_metacells)
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
                noisy_gene_fold,
                is_noisy_of_genes,
            )
        @assert size(significant_fold_factors_of_genes_in_other_metacells) == (n_genes, base_metacell - 1)

        is_deviant_sphere_of_gene_in_other_metacells =
            significant_fold_factors_of_genes_in_other_metacells .> max_sphere_diameter
        is_deviant_neighborhood_of_genes_in_other_metacells =
            significant_fold_factors_of_genes_in_other_metacells .> max_sphere_diameter + max_neighborhood_diameter

        n_deviant_sphere_genes_of_other_metacells = vec(sum(is_deviant_sphere_of_gene_in_other_metacells; dims = 1))
        n_deviant_neighborhood_genes_of_other_metacells =
            vec(sum(is_deviant_neighborhood_of_genes_in_other_metacells; dims = 1))

        @assert length(n_deviant_sphere_genes_of_other_metacells) == base_metacell - 1
        @assert length(n_deviant_neighborhood_genes_of_other_metacells) == base_metacell - 1

        distance_for_spheres_between_metacells[1:(base_metacell - 1), base_metacell] .=
            n_deviant_sphere_genes_of_other_metacells
        distance_for_neighborhoods_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_neighborhood_genes_of_other_metacells

        distance_for_spheres_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_sphere_genes_of_other_metacells
        distance_for_neighborhoods_between_metacells[base_metacell, 1:(base_metacell - 1)] .=
            n_deviant_neighborhood_genes_of_other_metacells
    end

    return distance_for_spheres_between_metacells, distance_for_neighborhoods_between_metacells
end

@inline function gene_distance(  # untested
    min_significant_gene_UMIs::Integer,
    total_UMIs_of_gene_in_base_metacell::Integer,
    log_decreased_fraction_of_gene_in_base_metacell::AbstractFloat,
    log_increased_fraction_of_gene_in_base_metacell::AbstractFloat,
    total_UMIs_of_genes_in_other_metacell::Integer,
    log_decreased_fraction_of_gene_in_other_metacell::AbstractFloat,
    log_increased_fraction_of_gene_in_other_metacell::AbstractFloat,
    noisy_gene_fold::AbstractFloat,
    gene_is_noisy::Bool,
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
        max.(
            log_decreased_high_fraction_of_gene - log_increased_low_fraction_of_gene - gene_is_noisy * noisy_gene_fold,
            0.0,
        )
    )
end

@logged function cluster_metacells_to_spheres(;  # untested
    distances_between_metacells::Matrix{Float32},
    threshold::AbstractFloat,
)::Vector{UInt32}
    check_efficient_action("compute_spheres", Columns, "distances_between_metacells", distances_between_metacells)
    clusters = hclust(distances_between_metacells; linkage = :complete)  # NOJET
    spheres_of_metacells = Vector{UInt32}(cutree(clusters; h = threshold))
    return spheres_of_metacells
end

@logged function compute_sphere_distances(;  # untested
    distances_between_metacells::Matrix{Float32},
    metacells_of_spheres::Vector{Vector{UInt32}},
)::Matrix{Float32}
    check_efficient_action("compute_spheres", Columns, "distances_between_metacells", distances_between_metacells)

    n_metacells = size(distances_between_metacells, 1)
    @assert size(distances_between_metacells, 2) == n_metacells

    n_spheres = length(metacells_of_spheres)
    distances_between_spheres = Matrix{Float32}(undef, n_spheres, n_spheres)

    @threads for base_sphere in reverse(1:n_spheres)
        metacells_of_base_sphere = metacells_of_spheres[base_sphere]

        @views distance_of_metacells_from_base_sphere_metacells =
            distances_between_metacells[:, metacells_of_base_sphere]

        distance_of_metacells_from_base_sphere = maximum(distance_of_metacells_from_base_sphere_metacells; dims = 2)
        @assert length(distance_of_metacells_from_base_sphere) == n_metacells

        for other_sphere in 1:base_sphere
            metacells_of_other_sphere = metacells_of_spheres[other_sphere]

            @views distance_of_other_sphere_metacells_from_base_sphere =
                distance_of_metacells_from_base_sphere[metacells_of_other_sphere]

            distance_between_other_and_base_spheres = maximum(distance_of_other_sphere_metacells_from_base_sphere)

            distances_between_spheres[base_sphere, other_sphere] = distance_between_other_and_base_spheres
            distances_between_spheres[other_sphere, base_sphere] = distance_between_other_and_base_spheres
        end
    end

    return distances_between_spheres
end

@logged function collect_spheres_of_neighborhoods(;  # untested
    distances_between_spheres::Matrix{Float32},
    threshold::AbstractFloat,
)::Tuple{Vector{Vector{<:Integer}}, Vector{UInt32}}
    check_efficient_action("compute_spheres", Columns, "distances_between_spheres", distances_between_spheres)

    n_spheres = size(distances_between_spheres, 1)
    @assert size(distances_between_spheres, 2) == n_spheres

    @views spheres_of_neighborhoods =  # NOJET
        [findall(distances_between_spheres[:, sphere] .<= threshold) for sphere in 1:n_spheres]
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

@logged function identify_spheres_correlated_genes(  # untested
    daf::DafReader;
    gene_fraction_regularization::AbstractFloat,
    min_gene_correlation::AbstractFloat,
    is_correlated_of_genes_in_metacells::Matrix{Bool},
    metacells_of_spheres::Vector{Vector{UInt32}},
)::Matrix{Bool}
    check_efficient_action(
        "compute_spheres",
        Columns,
        "is_correlated_of_genes_in_metacells",
        is_correlated_of_genes_in_metacells,
    )
    n_genes = axis_length(daf, "gene")
    n_spheres = axis_length(daf, "sphere")
    sphere_names = axis_array(daf, "sphere")

    is_correlated_of_genes_in_spheres = zeros(Bool, n_genes, n_spheres)
    main_neighborhoods_of_spheres = get_vector(daf, "sphere", "neighborhood.main")

    @threads for sphere_index in 1:n_spheres
        sphere_name = sphere_names[sphere_index]
        neighborhood_name = main_neighborhoods_of_spheres[sphere_index]
        neighborhood_daf = chain_writer(
            [
                viewer(
                    daf;
                    name = "intermediate.$(sphere_name)",
                    axes = [
                        "metacell" => "/ metacell & sphere => neighborhood.main = $(neighborhood_name)",
                        "gene" => "=",
                    ],
                    data = [("metacell", "gene", "fraction") => "="],
                ),
                MemoryDaf(; name = "intermediate.$(sphere_name).identified"),
            ];
            name = "intermediate.$(neighborhood_name).chain",
        )

        metacell_indices_of_sphere = metacells_of_spheres[sphere_index]
        @views is_correlated_of_genes_in_metacells_of_sphere =
            is_correlated_of_genes_in_metacells[:, metacell_indices_of_sphere]
        is_correlated_of_genes_in_sphere = vec(maximum(is_correlated_of_genes_in_metacells_of_sphere; dims = 2))
        @assert length(is_correlated_of_genes_in_sphere) == n_genes

        n_metacells_of_neighborhood = axis_length(neighborhood_daf, "metacell")
        if n_metacells_of_neighborhood > 1
            identify_correlated_genes!(
                neighborhood_daf;
                gene_fraction_regularization = gene_fraction_regularization,
                min_gene_correlation = min_gene_correlation,
            )
            is_correlated_of_genes_in_neighborhood = neighborhood_daf["/ gene : is_correlated || true"].array
            is_correlated_of_genes_in_sphere .&= is_correlated_of_genes_in_neighborhood
        end

        is_correlated_of_genes_in_spheres[:, sphere_index] .= is_correlated_of_genes_in_sphere
        is_correlated_of_genes_in_metacells[:, metacell_indices_of_sphere] .= is_correlated_of_genes_in_sphere
    end

    return is_correlated_of_genes_in_spheres
end

@logged function compute_followup_metacells_distances(;  # untested
    total_UMIs_of_genes_in_metacells::AbstractMatrix{<:Unsigned},
    log_decreased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    log_increased_fraction_of_genes_in_metacells::AbstractMatrix{<:AbstractFloat},
    metacells_of_spheres::Vector{Vector{UInt32}},
    spheres_of_neighborhoods::Vector{Vector{<:Integer}},
    main_neighborhoods_of_spheres::Vector{UInt32},
    is_correlated_of_genes_in_spheres::Matrix{Bool},
    is_noisy_of_genes::AbstractVector{Bool},
    min_significant_gene_UMIs::Integer,
    noisy_gene_fold::AbstractFloat,
)::Matrix{Float32}
    check_efficient_action(
        "compute_spheres",
        Columns,
        "total_UMIs_of_genes_in_metacells",
        total_UMIs_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_spheres",
        Columns,
        "log_decreased_fraction_of_genes_in_metacells",
        log_decreased_fraction_of_genes_in_metacells,
    )
    check_efficient_action(
        "compute_spheres",
        Columns,
        "log_increased_fraction_of_genes_in_metacells",
        log_increased_fraction_of_genes_in_metacells,
    )
    @assert size(log_decreased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)
    @assert size(log_increased_fraction_of_genes_in_metacells) == size(total_UMIs_of_genes_in_metacells)

    n_metacells = size(total_UMIs_of_genes_in_metacells, 2)
    distances_between_metacells = fill(Inf32, n_metacells, n_metacells)

    n_spheres = length(metacells_of_spheres)
    did_compute_distance_between_spheres = zeros(Bool, n_spheres, n_spheres)

    @threads for base_sphere_index in 1:n_spheres
        metacell_indices_of_base_sphere = metacells_of_spheres[base_sphere_index]

        @views is_correlated_of_genes_in_base_sphere = is_correlated_of_genes_in_spheres[:, base_sphere_index]

        base_neighborhood_index = main_neighborhoods_of_spheres[base_sphere_index]
        other_sphere_indices = spheres_of_neighborhoods[base_neighborhood_index]
        for other_sphere_index in other_sphere_indices
            if did_compute_distance_between_spheres[base_sphere_index, other_sphere_index]
                continue
            end
            did_compute_distance_between_spheres[base_sphere_index, other_sphere_index] = true
            did_compute_distance_between_spheres[other_sphere_index, base_sphere_index] = true

            metacell_indices_of_other_sphere = metacells_of_spheres[other_sphere_index]
            n_metacells_of_other_sphere = length(metacell_indices_of_other_sphere)

            @views is_correlated_of_genes_in_other_sphere = is_correlated_of_genes_in_spheres[:, other_sphere_index]
            is_correlated_of_genes_for_distance =
                is_correlated_of_genes_in_base_sphere .| is_correlated_of_genes_in_other_sphere
            if !any(is_correlated_of_genes_for_distance)
                distances_between_metacells[metacell_indices_of_base_sphere, metacell_indices_of_other_sphere] .= 0
                distances_between_metacells[metacell_indices_of_other_sphere, metacell_indices_of_base_sphere] .= 0
                continue
            end
            is_noisy_of_genes_for_distance = is_noisy_of_genes[is_correlated_of_genes_for_distance]

            for base_metacell_index in metacell_indices_of_base_sphere
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

                total_UMIs_of_genes_in_other_metacells = total_UMIs_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    metacell_indices_of_other_sphere,
                ]
                log_decreased_fraction_of_genes_in_other_metacells = log_decreased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    metacell_indices_of_other_sphere,
                ]
                log_increased_fraction_of_genes_in_other_metacells = log_increased_fraction_of_genes_in_metacells[
                    is_correlated_of_genes_for_distance,
                    metacell_indices_of_other_sphere,
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
                        noisy_gene_fold,
                        is_noisy_of_genes_for_distance,
                    )
                @assert size(significant_fold_factors_of_genes_in_other_metacells, 2) == n_metacells_of_other_sphere

                maximal_fold_factors_of_genes_in_other_metacells =
                    vec(maximum(significant_fold_factors_of_genes_in_other_metacells; dims = 1))
                @assert length(maximal_fold_factors_of_genes_in_other_metacells) == n_metacells_of_other_sphere

                distances_between_metacells[base_metacell_index, metacell_indices_of_other_sphere] .=
                    maximal_fold_factors_of_genes_in_other_metacells
                distances_between_metacells[metacell_indices_of_other_sphere, base_metacell_index] .=
                    maximal_fold_factors_of_genes_in_other_metacells
            end
        end
    end
    return distances_between_metacells
end

@logged function compute_membership_matrix(  # untested
    spheres_of_neighborhoods::Vector{Vector{<:Integer}},
    n_spheres::Integer,
)::AbstractMatrix{Bool}
    n_neighborhoods = length(spheres_of_neighborhoods)
    is_member_of_spheres_in_neighborhoods = zeros(Bool, n_spheres, n_neighborhoods)
    for (neighborhood, spheres_of_neighborhood) in enumerate(spheres_of_neighborhoods)
        is_member_of_spheres_in_neighborhoods[spheres_of_neighborhood, neighborhood] .= true  # NOJET
    end
    return is_member_of_spheres_in_neighborhoods
end

end  # module
