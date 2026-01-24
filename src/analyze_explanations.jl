"""
Check which genes are well-explained by the regulators and lateral genes.
"""
module AnalyzeExplanations

export compute_blocks_genes_explanations!

using Base.Threads
using LinearAlgebra
using TanayLabUtilities
using StatsBase

using ..Contracts

using DataAxesFormats
struct GeneGeneStorage
    correlation_per_gene_per_gene::Matrix{Float32}
end

function TanayLabUtilities.reset_reusable_storage!(storage::GeneGeneStorage)::Nothing
    fill!(storage.correlation_per_gene_per_gene, 0)
    return nothing
end

"""
    compute_blocks_genes_explanations!(
        daf::DafWriter;
        min_explainer_gene_correlation::AbstractFloat = $(DEFAULT.min_explainer_gene_correlation),
        min_unexplained_gene_correlation::AbstractFloat = $(DEFAULT.min_unexplained_gene_correlation),
        min_unexplained_correlation_fraction::AbstractFloat = $(DEFAULT.min_unexplained_correlation_fraction),
        min_unexplained_gene_fraction_in_environment::AbstractFloat = $(DEFAULT.min_unexplained_gene_fraction_in_environment),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute for each block for each gene how well the gene is explained by the regulators and lateral genes in the environment.

TODOX
"""
@logged @computation Contract(
    axes = [block_axis(RequiredInput), metacell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        gene_is_regulator_vector(RequiredInput),
        gene_is_lateral_vector(RequiredInput),
        block_gene_is_environment_marker_matrix(RequiredInput),
        metacell_block_vector(RequiredInput),
        metacell_gene_linear_fraction_matrix(RequiredInput),
        metacell_gene_log_linear_fraction_matrix(RequiredInput),
        block_block_is_in_environment_matrix(RequiredInput),
        block_gene_most_correlated_gene_in_environment_matrix(GuaranteedOutput),
        block_gene_most_correlation_in_environment_matrix(GuaranteedOutput),
        block_gene_most_correlated_lateral_gene_in_environment_matrix(GuaranteedOutput),
        block_gene_most_correlation_with_lateral_in_environment_matrix(GuaranteedOutput),
        block_gene_most_correlated_regulator_gene_in_environment_matrix(GuaranteedOutput),
        block_gene_most_correlation_with_regulator_in_environment_matrix(GuaranteedOutput),
        block_gene_unexplained_correlation_in_environment_matrix(GuaranteedOutput),
        block_gene_is_unexplained_in_environment_matrix(GuaranteedOutput),
        # TODOX block_gene_gene_correlation_in_environment_tensor(GuaranteedOutput),
        block_n_unexplained_genes_in_environment_vector(GuaranteedOutput),
    ],
) function compute_blocks_genes_explanations!(
    daf::DafWriter;
    min_explainer_gene_correlation::AbstractFloat = 0.6,
    min_unexplained_gene_correlation::AbstractFloat = 0.4,
    min_unexplained_correlation_fraction::AbstractFloat = 0.2,
    min_unexplained_gene_fraction_in_environment::AbstractFloat = 1e-4,
    overwrite::Bool = false,
)::Nothing
    @assert 0 <= min_explainer_gene_correlation <= 1
    @assert 0 <= min_unexplained_gene_correlation <= 1
    @assert 0 <= min_unexplained_correlation_fraction <= 1

    n_genes = axis_length(daf, "gene")
    name_per_gene = axis_vector(daf, "gene")

    n_blocks = axis_length(daf, "block")
    name_per_block = axis_vector(daf, "block")

    most_correlation_per_gene_per_block = zeros(Float32, n_genes, n_blocks)
    most_correlation_with_lateral_per_gene_per_block = zeros(Float32, n_genes, n_blocks)
    most_correlation_with_regulator_per_gene_per_block = zeros(Float32, n_genes, n_blocks)

    most_correlated_gene_per_gene_per_block = Matrix{AbstractString}(undef, n_genes, n_blocks)
    most_correlated_lateral_gene_per_gene_per_block = Matrix{AbstractString}(undef, n_genes, n_blocks)
    most_correlated_regulator_gene_per_gene_per_block = Matrix{AbstractString}(undef, n_genes, n_blocks)

    fill!(most_correlated_gene_per_gene_per_block, "")
    fill!(most_correlated_lateral_gene_per_gene_per_block, "")
    fill!(most_correlated_regulator_gene_per_gene_per_block, "")

    unexplained_correlation_per_gene_per_block = Matrix{Float32}(undef, n_genes, n_blocks)
    is_unexplained_per_gene_per_block = zeros(Bool, n_genes, n_blocks)

    n_unexplained_genes_per_block = Vector{UInt32}(undef, n_blocks)

    log_linear_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "log_linear_fraction")
    linear_fraction_per_metacell_per_gene = get_matrix(daf, "metacell", "gene", "linear_fraction")

    reusable_storage = ReusableStorage() do
        return GeneGeneStorage(zeros(Float32, n_genes, n_genes))
    end

    parallel_loop_wo_rng(1:n_blocks) do block_index
        block_name = name_per_block[block_index]

        indices_of_environment_metacells = daf["/ metacell & block => is_in_environment ;= $(block_name) : index"].array
        indices_of_environment_markers = daf["/ gene & is_environment_marker ; block = $(block_name) : index"].array
        n_environment_markers = length(indices_of_environment_markers)

        is_lateral_per_environment_marker =
            daf["/ gene & is_environment_marker ; block = $(block_name) : is_lateral"].array
        is_regulator_per_environment_marker =
            daf["/ gene & is_environment_marker ; block = $(block_name) : is_regulator"].array .&
            .! is_lateral_per_environment_marker

        linear_fraction_per_environment_metacell_per_environment_marker =
            linear_fraction_per_metacell_per_gene[indices_of_environment_metacells, indices_of_environment_markers]
        maximal_linear_fraction_per_environment_marker =
            vec(maximum(linear_fraction_per_environment_metacell_per_environment_marker; dims = 1))
        @assert_vector(maximal_linear_fraction_per_environment_marker, n_environment_markers)
        is_significant_per_environment_marker =
            maximal_linear_fraction_per_environment_marker .>= min_unexplained_gene_fraction_in_environment

        indices_of_significant_genes = indices_of_environment_markers[is_significant_per_environment_marker]
        n_significant_genes = length(indices_of_significant_genes)
        @assert length(n_significant_genes) > 0
        is_lateral_per_significant_gene = is_lateral_per_environment_marker[is_significant_per_environment_marker]
        is_regulator_per_significant_gene = is_regulator_per_environment_marker[is_significant_per_environment_marker]

        log_linear_fraction_per_environment_metacell_per_significant_gene =
            log_linear_fraction_per_metacell_per_gene[indices_of_environment_metacells, indices_of_significant_genes]

        correlation_between_significant_genes = flame_timed("cor") do
            return cor(log_linear_fraction_per_environment_metacell_per_significant_gene)  # NOLINT
        end
        @assert_matrix(correlation_between_significant_genes, n_significant_genes, n_significant_genes)
        correlation_between_significant_genes[diagind(correlation_between_significant_genes)] .= 0

        most_correlated_significant_gene_position_per_significant_gene =
            argmax.(eachcol(correlation_between_significant_genes))
        most_correlated_gene_per_gene_per_block[indices_of_significant_genes, block_index] .=
            name_per_gene[indices_of_significant_genes[most_correlated_significant_gene_position_per_significant_gene]]

        most_correlation_per_significant_gene = vec(maximum(correlation_between_significant_genes; dims = 1))
        @assert_vector(most_correlation_per_significant_gene, n_significant_genes)
        most_correlation_per_gene_per_block[indices_of_significant_genes, block_index] .=
            most_correlation_per_significant_gene

        indices_of_lateral_significant_genes = indices_of_significant_genes[is_lateral_per_significant_gene]
        if length(indices_of_lateral_significant_genes) == 0
            most_correlation_with_lateral_per_significant_gene = fill(0.0, n_significant_genes)
        else
            @views correlation_between_lateral_and_significant_genes =
                correlation_between_significant_genes[is_lateral_per_significant_gene, :]
            most_correlated_lateral_significant_gene_index_per_significant_gene =
                argmax.(eachcol(correlation_between_lateral_and_significant_genes))
            indices_of_most_correlated_lateral_genes_per_significant_gene =
                indices_of_lateral_significant_genes[most_correlated_lateral_significant_gene_index_per_significant_gene]
            most_correlated_lateral_gene_per_gene_per_block[indices_of_significant_genes, block_index] .=
                name_per_gene[indices_of_most_correlated_lateral_genes_per_significant_gene]
            most_correlation_with_lateral_per_significant_gene =
                vec(maximum(correlation_between_lateral_and_significant_genes; dims = 1))
            most_correlation_with_lateral_per_gene_per_block[indices_of_significant_genes, block_index] .=
                most_correlation_with_lateral_per_significant_gene
        end

        indices_of_regulator_significant_genes = indices_of_significant_genes[is_regulator_per_significant_gene]
        if length(indices_of_regulator_significant_genes) == 0
            most_correlation_with_regulator_per_significant_gene = fill(0.0, n_significant_genes)
        else
            @views correlation_between_regulator_and_significant_genes =
                correlation_between_significant_genes[is_regulator_per_significant_gene, :]
            most_correlated_regulator_significant_gene_index_per_significant_gene =
                argmax.(eachcol(correlation_between_regulator_and_significant_genes))
            indices_of_most_correlated_regulator_genes_per_significant_gene =
                indices_of_regulator_significant_genes[most_correlated_regulator_significant_gene_index_per_significant_gene]
            most_correlated_regulator_gene_per_gene_per_block[indices_of_significant_genes, block_index] .=
                name_per_gene[indices_of_most_correlated_regulator_genes_per_significant_gene]
            most_correlation_with_regulator_per_significant_gene =
                vec(maximum(correlation_between_regulator_and_significant_genes; dims = 1))
            most_correlation_with_regulator_per_gene_per_block[indices_of_significant_genes, block_index] .=
                most_correlation_with_regulator_per_significant_gene
        end

        explained_correlation_per_significant_gene = max.(
            most_correlation_with_lateral_per_significant_gene,
            most_correlation_with_regulator_per_significant_gene,
        )
        unexplained_correlation_per_significant_gene =
            most_correlation_per_significant_gene .- explained_correlation_per_significant_gene
        unexplained_correlation_per_gene_per_block[indices_of_significant_genes, block_index] .=
            unexplained_correlation_per_significant_gene

        is_unexplained_per_significant_gene =
            (most_correlation_per_significant_gene .>= min_unexplained_gene_correlation) .&
            (explained_correlation_per_significant_gene .< min_explainer_gene_correlation) .& (
                unexplained_correlation_per_significant_gene .>=
                most_correlation_per_significant_gene .* min_unexplained_correlation_fraction
            ) .& (.!is_regulator_per_significant_gene) .& (.!is_lateral_per_significant_gene)
        is_unexplained_per_gene_per_block[indices_of_significant_genes, block_index] .=
            is_unexplained_per_significant_gene

        indices_of_unexplained_genes = indices_of_significant_genes[is_unexplained_per_significant_gene]
        n_unexplained_genes = length(indices_of_unexplained_genes)

        correlation_between_significant_genes[diagind(correlation_between_significant_genes)] .= 1

        if false  # TODOX
            with_reusable(reusable_storage) do storage
                storage.correlation_per_gene_per_gene[indices_of_significant_genes, indices_of_significant_genes] .=
                    correlation_between_significant_genes
                return set_matrix!(
                    daf,
                    "gene",
                    "gene",
                    "$(block_name)_correlation_in_environment",
                    sparsify(storage.correlation_per_gene_per_gene);
                    overwrite,
                )
            end
        end

        n_unexplained_genes_per_block[block_index] = n_unexplained_genes

        @debug (
            "- Block: $(block_name)" *
            " Markers: $(n_environment_markers) " *
            " Significant: $(n_significant_genes) " *
            " Unexplained: $(n_unexplained_genes)"
        )
        return nothing
    end

    set_matrix!(
        daf,
        "gene",
        "block",
        "gene.most_correlated_in_environment",
        most_correlated_gene_per_gene_per_block;
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "most_correlation_in_environment",
        bestify(most_correlation_per_gene_per_block);
        overwrite,
    )

    set_matrix!(
        daf,
        "gene",
        "block",
        "gene.most_correlated_lateral_in_environment",
        most_correlated_lateral_gene_per_gene_per_block;
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "most_correlation_with_lateral_in_environment",
        bestify(most_correlation_with_lateral_per_gene_per_block);
        overwrite,
    )

    set_matrix!(
        daf,
        "gene",
        "block",
        "gene.most_correlated_regulator_in_environment",
        most_correlated_regulator_gene_per_gene_per_block;
        overwrite,
    )
    set_matrix!(
        daf,
        "gene",
        "block",
        "most_correlation_with_regulator_in_environment",
        bestify(most_correlation_with_regulator_per_gene_per_block);
        overwrite,
    )

    set_matrix!(
        daf,
        "gene",
        "block",
        "unexplained_correlation_in_environment",
        unexplained_correlation_per_gene_per_block;
        overwrite,
    )
    set_matrix!(daf, "gene", "block", "is_unexplained_in_environment", is_unexplained_per_gene_per_block; overwrite)

    set_vector!(daf, "block", "n_unexplained_genes_in_environment", n_unexplained_genes_per_block; overwrite)

    @debug "Mean Unexplained Genes per Block: $(mean(n_unexplained_genes_per_block))"  # NOLINT

    return nothing
end

end
