"""
Run whole stages of the metacells pipeline.

The functions here run multiple computation functions in sequence. They allow performing the computational pipeline of
sharpening metacells with few high-level calls instead of working through each low-level step on its own. None of the
functions here do anything other than invoking the lower level computations. As such they are purely convenience
functions. One is free to call the lower level steps directly in case the pipeline needs to be tweaked for some special
case.
"""
module Pipeline

export analyze_metacells!
export import_base_metacells!
export prepare_metacells!
export qc_metacells!
export run!
export sharpen_round!
export SharpeningRound
export SharpeningRounds
export sharpening_rounds

using DataAxesFormats
using Random
using TanayLabUtilities

using ..AnalyzeBlocks
using ..AnalyzeCells
using ..AnalyzeGenes
using ..AnalyzeMetacells
using ..AnalyzeModules
using ..ComputeBlocks
using ..ComputeModules
using ..Contracts
using ..ProjectCells
using ..SharpenMetacells

import Random.default_rng

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.vector_of_metacell_per_cell

"""
    prepare_metacells!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Aggregate the cells of each metacell, given which cells belong to which metacell. This computes everything about the
metacells which follows from the cells alone, per-metacell properties, and basic gene properties derived from metacells
alone. Gene properties that depend on gene masks (other than exclusion) are not computed here.

$(CONTRACT)
"""
@logged :mcs_ops @computation (
    optional_contract(function_contract(compute_vector_of_type_per_metacell_by_cells!)) |>
    function_contract(compute_matrix_of_UMIs_per_gene_per_metacell!) |>
    function_contract(compute_vector_of_total_UMIs_per_metacell!) |>
    function_contract(compute_vector_of_n_cells_per_metacell!) |>
    function_contract(compute_matrix_of_linear_fraction_per_gene_per_metacell!) |>
    function_contract(compute_matrix_of_log_linear_fraction_per_gene_per_metacell!) |>
    function_contract(compute_vector_of_is_marker_per_gene!) |>
    function_contract(compute_vector_of_marker_rank_per_gene!) |>
    function_contract(compute_matrix_of_correlation_between_markers_per_gene_per_gene!)
) function prepare_metacells!(daf::DafWriter; overwrite::Bool = false)::Nothing  # UNTESTED
    # The types of the metacells come from the types of their cells, so without the one there is not the other.
    if has_vector(daf, "cell", "type")
        compute_vector_of_type_per_metacell_by_cells!(daf; overwrite)
    end
    compute_matrix_of_UMIs_per_gene_per_metacell!(daf; overwrite)
    compute_vector_of_total_UMIs_per_metacell!(daf; overwrite)
    compute_vector_of_n_cells_per_metacell!(daf; overwrite)
    compute_matrix_of_linear_fraction_per_gene_per_metacell!(daf; overwrite)
    compute_matrix_of_log_linear_fraction_per_gene_per_metacell!(daf; overwrite)
    compute_vector_of_is_marker_per_gene!(daf; overwrite)
    compute_vector_of_marker_rank_per_gene!(daf; overwrite)
    compute_matrix_of_correlation_between_markers_per_gene_per_gene!(daf; overwrite)
    return nothing
end

"""
    analyze_metacells!(
        daf::DafWriter;
        prefix::AbstractString = $(DEFAULT.prefix),
        prev_daf::Maybe{DafReader} = nothing,
        module_status::Bool = $(DEFAULT.module_status),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute more advanced properties based on the metacells, taking into account gene masks. Specifically this depends on
the lateral and regulator gene masks, and the forbidden gene masks. These are used to determine the set of skeleton
genes which are then used to drive the rest of the analysis, starting with grouping metacells into blocks and ending
with local gene modules.

# Metacells

$(CONTRACT1)

# Previous Metacells

Only if `prev_daf` is given.

$(CONTRACT2)
"""
@logged :mcs_ops @computation (
    function_contract(compute_vector_of_is_skeleton_per_gene!) |>
    function_contract(compute_matrix_of_max_skeleton_fold_distance_between_metacells!) |>
    function_contract(compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!) |>
    function_contract(compute_metacells_2d_umap!) |>
    function_contract(compute_metacells_blocks!) |>
    function_contract(compute_matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block!) |>
    function_contract(compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!) |>
    function_contract(compute_vector_of_n_metacells_per_block!) |>
    function_contract(compute_vector_of_n_cells_per_block!) |>
    function_contract(compute_matrix_of_UMIs_per_gene_per_block!) |>
    function_contract(compute_vector_of_total_UMIs_per_block!) |>
    function_contract(compute_matrix_of_linear_fraction_per_gene_per_block!) |>
    function_contract(compute_matrix_of_log_linear_fraction_per_gene_per_block!) |>
    optional_contract(function_contract(compute_vector_of_type_per_block_by_metacells!)) |>
    function_contract(compute_vector_of_block_closest_by_pertinent_markers_per_cell!) |>
    function_contract(compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!) |>
    function_contract(compute_matrix_of_is_in_neighborhood_per_block_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_blocks_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_metacells_per_block!) |>
    function_contract(compute_vector_of_n_neighborhood_cells_per_block!) |>
    function_contract(compute_vector_of_total_neighborhood_UMIs_per_block!) |>
    function_contract(compute_matrix_of_is_neighborhood_marker_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_in_environment_per_metacell_per_block!) |>
    function_contract(compute_vector_of_n_environment_metacells_per_block!) |>
    function_contract(compute_vector_of_n_environment_cells_per_block!) |>
    function_contract(compute_vector_of_total_environment_UMIs_per_block!) |>
    function_contract(compute_matrix_of_is_environment_marker_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_environment_distinct_per_gene_per_block!) |>
    function_contract(compute_matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block!) |>
    function_contract(compute_blocks_modules!) |>
    function_contract(compute_vector_of_n_modules_per_block!) |>
    function_contract(compute_matrix_of_n_genes_per_module_per_block!) |>
    function_contract(compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!) |>
    function_contract(compute_matrix_of_cells_dispersion_per_metacell_per_module!)
) function_contract(compute_metacells_2d_umap!, 2) function analyze_metacells!(  # UNTESTED
    daf::DafWriter;
    prefix::AbstractString = "B",
    prev_daf::Maybe{DafReader} = nothing,
    module_status::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    # Which genes predict the rest, and the geometry of the metacells which follows from them.
    compute_vector_of_is_skeleton_per_gene!(daf; overwrite)
    compute_matrix_of_max_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_matrix_of_euclidean_skeleton_fold_distance_between_metacells!(daf; overwrite)
    compute_metacells_2d_umap!(daf; prev_daf, rng, overwrite)

    # The blocks - regions of the manifold the metacells fall into - and what each is made of.
    compute_metacells_blocks!(daf; prefix, overwrite)
    compute_matrix_of_mean_euclidean_skeleton_fold_distance_per_metacell_per_block!(daf; overwrite)
    compute_matrix_of_mean_euclidean_skeleton_fold_distance_between_blocks!(daf; overwrite)
    compute_vector_of_n_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_cells_per_block!(daf; overwrite)
    compute_matrix_of_UMIs_per_gene_per_block!(daf; overwrite)
    compute_vector_of_total_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_linear_fraction_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_log_linear_fraction_per_gene_per_block!(daf; overwrite)

    # A block has a type only when the metacells have one, which is optional data.
    if has_vector(daf, "metacell", "type")
        compute_vector_of_type_per_block_by_metacells!(daf; overwrite)
    end

    compute_vector_of_block_closest_by_pertinent_markers_per_cell!(daf; overwrite)
    compute_matrix_of_confusion_by_closest_by_pertinent_markers_per_block_per_block!(daf; overwrite)

    # The neighborhood of a block is the blocks close enough to it to be describing the same local behavior.
    compute_matrix_of_is_in_neighborhood_per_block_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_blocks_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_neighborhood_cells_per_block!(daf; overwrite)
    compute_vector_of_total_neighborhood_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_is_neighborhood_marker_per_gene_per_block!(daf; overwrite)

    # The environment extends the neighborhood with metacells close enough to the block, which gives the gene modules
    # below more metacells to be estimated from.
    compute_matrix_of_is_in_environment_per_metacell_per_block!(daf; overwrite)
    compute_vector_of_n_environment_metacells_per_block!(daf; overwrite)
    compute_vector_of_n_environment_cells_per_block!(daf; overwrite)
    compute_vector_of_total_environment_UMIs_per_block!(daf; overwrite)
    compute_matrix_of_is_environment_marker_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_is_environment_distinct_per_gene_per_block!(daf; overwrite)
    compute_matrix_of_is_correlated_with_skeleton_in_environment_per_gene_per_block!(daf; overwrite)

    # The gene modules of each block - the local programs the sharpening clusters cells by.
    compute_blocks_modules!(daf; module_status, rng, overwrite)
    compute_vector_of_n_modules_per_block!(daf; overwrite)
    compute_matrix_of_n_genes_per_module_per_block!(daf; overwrite)
    compute_stats_of_linear_fraction_in_environment_cells_per_module_per_block!(daf; overwrite)
    compute_matrix_of_cells_dispersion_per_metacell_per_module!(daf; overwrite)

    return nothing
end

"""
    function import_base_metacells!(;
        cells_daf::DafWriter,
        metacells_daf::DafWriter,
        metacell_per_cell::AbstractVector{<:AbstractString},
        empty_metacells::Maybe{EmptyImplicit} = nothing,
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Bring in the base metacells the sharpening pipeline starts with. The cells are expected to have been imported already,
but their assignment to metacells is supplied explicitly. It is filtered through `empty_metacells` to allow for values
like `Outliers` to be safely converted to the empty string before being applied.

# Cells

$(CONTRACT1)

# Metacells

$(CONTRACT2)
"""
@logged :mcs_ops @computation function_contract(compute_vector_of_is_base_outlier_per_cell!, 1) renamed_contract(
    # The metacell of each cell is set directly rather than computed, so it has no contract to combine.
    Contract(;
        name = "metacells_daf",
        axes = [cell_axis(RequiredInput)],
        data = [vector_of_metacell_per_cell(GuaranteedOutput)],
    ) |> function_contract(compute_vector_of_is_base_outlier_per_cell!, 2),
    "metacells_daf",
) function import_base_metacells!(;  # UNTESTED
    cells_daf::DafWriter,
    metacells_daf::DafWriter,
    metacell_per_cell::AbstractVector{<:AbstractString},
    empty_metacells::Maybe{EmptyImplicit} = nothing,
    overwrite::Bool = false,
)::Nothing
    set_vector!(metacells_daf, "cell", "metacell", metacell_per_cell; overwrite)
    if empty_metacells !== nothing
        unify_empty_vector_values!(metacells_daf; axis = "cell", property = "metacell", empty_values = empty_metacells)
    end
    compute_vector_of_is_base_outlier_per_cell!(; cells_daf, metacells_daf, overwrite)
    return nothing
end

"""
    function qc_metacells!(;
        daf::DafWriter,
        score_daf::DafReader,
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Say how well the metacells describe the cells they were aggregated from, judged in the locations of the manifold the
`score_daf` repository laid out.

Each cell is correlated against its own metacell minus itself, over the cells of each base block's neighborhood, and
the result is averaged over the environment marker genes of that base block. The base blocks are the blocks of the
`score_daf`. A higher number means the metacells of that location are a better account of the cells there. Only
metacells scored against the same `score_daf` can be compared to each other, since it is the `score_daf` which decides
both the locations and the genes each of them is read by. A repository is free to be its own `score_daf`, which is how
the metacells the sharpening started from are scored.

# Daf

$(CONTRACT1)

# Score

$(CONTRACT2)
"""
@logged :mcs_ops @computation renamed_contract(
    function_contract(
        compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!,
        1,
    ) |> function_contract(
        compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!,
        1,
    ),
    "daf",
) renamed_contract(
    function_contract(
        compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!,
        2,
    ) |> function_contract(
        compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!,
        2,
    ),
    "score_daf",
) function qc_metacells!(; daf::DafWriter, score_daf::DafReader, overwrite::Bool = false)::Nothing  # UNTESTED
    compute_matrix_of_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_gene_per_base_block!(;
        other_daf = daf,
        base_daf = score_daf,
        overwrite,
    )
    compute_vector_of_mean_correlation_between_base_neighborhood_cells_and_punctuated_metacells_per_base_block!(;
        other_daf = daf,
        base_daf = score_daf,
        overwrite,
    )
    return nothing
end

SHARPEN_ROUND_PREV_CONTRACT =
    renamed_contract(function_contract(sharpen_metacells!, 2) |> function_contract(analyze_metacells!, 2), "prev_daf")

"""
    function sharpen_round!(;
        sharp_daf::DafWriter,
        prev_daf::DafReader,
        score_daf::DafReader,
        sharpening_round::Integer,
        metacells_prefix::AbstractString = $(DEFAULT.metacells_prefix),
        blocks_prefix::AbstractString = $(DEFAULT.blocks_prefix),
        module_status::Bool = $(DEFAULT.module_status),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Run one round of sharpening. This regroups the cells into new metacells in `sharp_daf`, using what `prev_daf` says about
the manifold ([`sharpen_metacells!`](@ref)). It then aggregates the new metacells ([`prepare_metacells!`](@ref)) and
analyzes them ([`analyze_metacells!`](@ref)). Finally, it scores them against `score_daf` ([`qc_metacells!`](@ref)).

The new metacells are named using the `metacells_prefix`, and their blocks using the `blocks_prefix`. The prefixes
[`sharpening_rounds`](@ref) gives each round are one way to choose them.

# Sharpened Metacells

$(CONTRACT1)

# Previous Metacells

$(CONTRACT2)

# Score

$(CONTRACT3)
"""
@logged :mcs_ops @computation renamed_contract(
    function_contract(sharpen_metacells!, 1) |>
    function_contract(prepare_metacells!) |>
    function_contract(analyze_metacells!, 1) |>
    function_contract(qc_metacells!, 1),
    "sharp_daf",
) SHARPEN_ROUND_PREV_CONTRACT function_contract(qc_metacells!, 2) function sharpen_round!(;  # UNTESTED
    sharp_daf::DafWriter,
    prev_daf::DafReader,
    score_daf::DafReader,
    sharpening_round::Integer,
    metacells_prefix::AbstractString = "M",
    blocks_prefix::AbstractString = "B",
    module_status::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::Nothing
    sharpen_metacells!(; sharp_daf, prev_daf, prefix = metacells_prefix, sharpening_round, rng, overwrite)
    prepare_metacells!(sharp_daf; overwrite)
    analyze_metacells!(sharp_daf; prefix = blocks_prefix, prev_daf, module_status, rng, overwrite)
    qc_metacells!(; daf = sharp_daf, score_daf, overwrite)
    return nothing
end

# The number of letters a round prefix is spelled with, starting at its base letter. Eleven letters give rounds 0 to 10
# a single letter each, and keep the letters of the metacells base `M` apart from those of the blocks base `B`.
const PREFIX_LETTERS = 11

# The letter a base prefix is. It must leave room for all the `PREFIX_LETTERS` letters before `Z`.
function base_prefix_letter(base_prefix::AbstractString)::Char
    @assert length(base_prefix) == 1 "the base prefix: $(base_prefix) is not a single letter"
    base_letter = base_prefix[1]
    @assert 'A' <= base_letter && base_letter + (PREFIX_LETTERS - 1) <= 'Z' "invalid base prefix: $(base_prefix)"
    return base_letter
end

# The prefix of the names of a round. This is the `sharpening_round`-th string (counting from 0) in shortlex order over
# the `PREFIX_LETTERS` letters starting at the base letter. Shortlex order lists the shorter strings first, and
# alphabetically within the same length. Round 0 is the base letter itself, and rounds 1 to 10 are the letters after it.
function round_prefix(base_prefix::AbstractString, sharpening_round::Integer)::String
    @assert sharpening_round >= 0
    base_letter = base_prefix_letter(base_prefix)

    n_letters = 1
    remaining = sharpening_round
    while remaining >= PREFIX_LETTERS^n_letters
        remaining -= PREFIX_LETTERS^n_letters
        n_letters += 1
    end

    letters = Vector{Char}(undef, n_letters)
    for position in n_letters:-1:1
        letters[position] = base_letter + remaining % PREFIX_LETTERS
        remaining ÷= PREFIX_LETTERS
    end
    return String(letters)
end

"""
    mutable struct SharpeningRounds ... end

The sharpening rounds, as returned by [`sharpening_rounds`](@ref). Iterating on it gives a [`SharpeningRound`](@ref)
for each round in turn. Each must be [`run!`](@ref) before the next one is taken.
"""
mutable struct SharpeningRounds
    base_daf::DafReader
    score_daf::DafReader
    directory::AbstractString
    metacells_prefix::AbstractString
    blocks_prefix::AbstractString
    module_status::Bool
    rng::AbstractRNG
    overwrite::Bool
    previous_daf::DafReader
    n_run_rounds::Int
end

"""
    struct SharpeningRound
        index::Int
        previous_daf::DafReader
        base_daf::DafReader
        score_daf::DafReader
        directory::AbstractString
        metacells_prefix::AbstractString
        blocks_prefix::AbstractString
        module_status::Bool
        rng::AbstractRNG
        overwrite::Bool
    end

One round of the [`sharpening_rounds`](@ref). The `index` is the number of the round, starting at 1. The
`previous_daf` is the metacells of the round before this one (for the first round, the `initial_daf` the rounds started
from). The `metacells_prefix` and the `blocks_prefix` are this round's prefixes. The rest are as given to
[`sharpening_rounds`](@ref).

Running the round using [`run!`](@ref) uses all of these. Any of the parameters can be overridden when calling it.
"""
struct SharpeningRound
    index::Int
    previous_daf::DafReader
    base_daf::DafReader
    score_daf::DafReader
    directory::AbstractString
    metacells_prefix::AbstractString
    blocks_prefix::AbstractString
    module_status::Bool
    rng::AbstractRNG
    overwrite::Bool
    rounds::SharpeningRounds
end

"""
    sharpening_rounds(;
        initial_daf::DafReader,
        base_daf::DafReader,
        score_daf::DafReader,
        directory::AbstractString,
        metacells_prefix::AbstractString = $(DEFAULT.metacells_prefix),
        blocks_prefix::AbstractString = $(DEFAULT.blocks_prefix),
        module_status::Bool = $(DEFAULT.module_status),
        rng::AbstractRNG = default_rng(),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::SharpeningRounds

Iterate on rounds of sharpening, starting from the `initial_daf` metacells. Each round is a [`SharpeningRound`](@ref).
Running it using [`run!`](@ref) creates a new repository in the `directory`, chained on the `base_daf`. It then fills
it using [`sharpen_round!`](@ref), which scores it against the `score_daf`. The iteration never ends by itself. The
caller decides when to stop:

```julia
for sharpening_round in sharpening_rounds(; initial_daf, base_daf, score_daf, directory = "dafs")
    sharpening_round.index > 2 && break
    metacells_daf = run!(sharpening_round, "metacells.R\$(sharpening_round.index)")
end
```

The `metacells_prefix` and the `blocks_prefix` are the base letters of the prefixes of each round. The `initial_daf`
metacells are round 0, and are expected to be named using the base letters. Each round's prefix is the next string, in
shortlex order, over the 11 letters starting at the base letter. That is, for `B`, the rounds are `C`, `D`, ... `L`,
then `BB`, `BC`, ... `LL`, then `BBB`, and so on. The letters of the two base prefixes must not overlap, so a block
prefix can never be mistaken for a metacells prefix.
"""
@documented function sharpening_rounds(;
    initial_daf::DafReader,
    base_daf::DafReader,
    score_daf::DafReader,
    directory::AbstractString,
    metacells_prefix::AbstractString = "M",
    blocks_prefix::AbstractString = "B",
    module_status::Bool = false,
    rng::AbstractRNG = default_rng(),
    overwrite::Bool = false,
)::SharpeningRounds
    distance = abs(base_prefix_letter(metacells_prefix) - base_prefix_letter(blocks_prefix))
    @assert distance >= PREFIX_LETTERS (
        "overlapping metacells prefix: $(metacells_prefix) and blocks prefix: $(blocks_prefix)"
    )
    return SharpeningRounds(
        base_daf,
        score_daf,
        directory,
        metacells_prefix,
        blocks_prefix,
        module_status,
        rng,
        overwrite,
        initial_daf,
        0,
    )
end

Base.IteratorSize(::Type{SharpeningRounds}) = Base.IsInfinite()

Base.eltype(::Type{SharpeningRounds}) = SharpeningRound

function Base.iterate(rounds::SharpeningRounds, index::Int = 1)::Tuple{SharpeningRound, Int}
    @assert rounds.n_run_rounds == index - 1 "the sharpening round: $(index - 1) was not run"
    sharpening_round = SharpeningRound(
        index,
        rounds.previous_daf,
        rounds.base_daf,
        rounds.score_daf,
        rounds.directory,
        round_prefix(rounds.metacells_prefix, index),
        round_prefix(rounds.blocks_prefix, index),
        rounds.module_status,
        rounds.rng,
        rounds.overwrite,
        rounds,
    )
    return (sharpening_round, index + 1)
end

"""
    run!(
        sharpening_round::SharpeningRound,
        name::AbstractString;
        metacells_prefix::AbstractString = sharpening_round.metacells_prefix,
        blocks_prefix::AbstractString = sharpening_round.blocks_prefix,
        module_status::Bool = sharpening_round.module_status,
        rng::AbstractRNG = sharpening_round.rng,
        overwrite::Bool = sharpening_round.overwrite,
    )::DafWriter

Run one of the [`sharpening_rounds`](@ref). This creates a repository with the `name` in the round's `directory`,
chained on the round's `base_daf`. It fills it using [`sharpen_round!`](@ref), and returns it. The next round starts
from it. Any of the round's parameters can be overridden for this round only.
"""
function run!(  # UNTESTED
    sharpening_round::SharpeningRound,
    name::AbstractString;
    metacells_prefix::AbstractString = sharpening_round.metacells_prefix,
    blocks_prefix::AbstractString = sharpening_round.blocks_prefix,
    module_status::Bool = sharpening_round.module_status,
    rng::AbstractRNG = sharpening_round.rng,
    overwrite::Bool = sharpening_round.overwrite,
)::DafWriter
    rounds = sharpening_round.rounds
    @assert rounds.n_run_rounds == sharpening_round.index - 1 "the sharpening round: $(sharpening_round.index) was run"

    sharp_daf = complete_chain!(;
        base_daf = sharpening_round.base_daf,
        new_daf = FilesDaf(joinpath(sharpening_round.directory, name), "w"; name),
        name,
    )
    sharpen_round!(;
        sharp_daf,
        prev_daf = sharpening_round.previous_daf,
        score_daf = sharpening_round.score_daf,
        sharpening_round = sharpening_round.index,
        metacells_prefix,
        blocks_prefix,
        module_status,
        rng,
        overwrite,
    )

    rounds.previous_daf = sharp_daf
    rounds.n_run_rounds = sharpening_round.index
    return sharp_daf
end

end  # module
