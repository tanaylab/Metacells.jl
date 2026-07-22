"""
Do simple per-cell analysis.
"""
module AnalyzeCells

export cell_group_membership_from_block_neighborhoods
export compute_vector_of_total_UMIs_per_cell!
export gather_gene_UMIs_per_region_cell!
export own_block_punctuated_correlation_per_gene_per_block!
export sparse_cell_reference_correlation_per_gene_per_group!
export sum_sparse_UMIs_per_gene_group_per_cell!

using Base.Threads
using DataAxesFormats
using SparseArrays
using StatsBase
using TanayLabUtilities

import Base.Threads.maxthreadid

using ..Contracts

# Needed because of JET:
import Metacells.Contracts.cell_axis
import Metacells.Contracts.gene_axis
import Metacells.Contracts.matrix_of_UMIs_per_gene_per_cell
import Metacells.Contracts.vector_of_is_excluded_per_gene
import Metacells.Contracts.vector_of_total_UMIs_per_cell

"""
    function compute_vector_of_total_UMIs_per_cell!(
        daf::DafWriter;
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Compute and set [`vector_of_total_UMIs_per_cell`](@ref).

$(CONTRACT)
"""
@logged :mcs_ops @computation Contract(;
    axes = [cell_axis(RequiredInput), gene_axis(RequiredInput)],
    data = [
        vector_of_is_excluded_per_gene(RequiredInput),
        matrix_of_UMIs_per_gene_per_cell(RequiredInput),
        vector_of_total_UMIs_per_cell(CreatedOutput),
    ],
) function compute_vector_of_total_UMIs_per_cell!(  # UNTESTED
    daf::DafWriter;
    overwrite::Bool = false,
)::Nothing
    total_UMIs_per_cell = daf["@ cell @ gene [ ! is_excluded ] :: UMIs >| Sum"].array
    set_vector!(daf, "cell", "total_UMIs", total_UMIs_per_cell; overwrite)
    @debug "Mean (included) UMIs per cell: $(mean(total_UMIs_per_cell))" _group = :mcs_results  # NOLINT
    return nothing
end

# Gather each region cell's UMIs of gene `gene_index` into the first `n_region_cells` entries of `UMIs_per_region_cell`
# (zeroed, then filled). `region_position_per_cell[cell_index]` is the cell's 1-based position in the region, or 0 if
# the cell is not in the region. `sparse_UMIs_per_cell_per_gene` is the CSC cell-UMIs matrix (rows = cells, columns =
# genes); walking the gene's column once replaces per-cell random `UMIs_per_cell_per_gene[cell_index, gene_index]`
# access, each of which is a binary search into the sparse column. Returns whether any region cell has nonzero UMIs of
# the gene (the all-zero early-exit check the callers need).
function gather_gene_UMIs_per_region_cell!(
    UMIs_per_region_cell::AbstractVector{<:Real},
    n_region_cells::Integer,
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC,
    gene_index::Integer,
    region_position_per_cell::AbstractVector{<:Integer},
)::Bool
    @views fill!(UMIs_per_region_cell[1:n_region_cells], 0)
    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)
    any_nonzero = false
    @inbounds for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
        region_position = region_position_per_cell[row_index_per_stored[stored_index]]
        if region_position > 0
            stored_UMIs = UMIs_per_stored[stored_index]
            UMIs_per_region_cell[region_position] = stored_UMIs
            any_nonzero |= stored_UMIs > 0
        end
    end
    return any_nonzero
end

# Build the (group X cell) boolean membership matrix for a block-neighborhood grouping: cell `cell` in block
# `block_index_per_cell[cell]` is a member of group `group` iff `is_in_neighborhood_per_block_per_group[block, group]`,
# for cells with `is_participating_per_cell`. Columns are cells, so a cell's groups are a single CSC column walk in
# `sparse_cell_reference_correlation_per_gene_per_group!`. Non-participating cells (no metacell / outside all
# neighborhoods) are members of no group.
function cell_group_membership_from_block_neighborhoods(
    block_index_per_cell::AbstractVector{<:Integer},
    is_in_neighborhood_per_block_per_group::AbstractMatrix{Bool},
    is_participating_per_cell::Union{AbstractVector{Bool}, BitVector},
)::SparseMatrixCSC{Bool, Int32}
    n_cells = length(block_index_per_cell)
    n_blocks, n_groups = size(is_in_neighborhood_per_block_per_group)
    @assert length(is_participating_per_cell) == n_cells

    groups_per_block = [Int32[] for _ in 1:n_blocks]
    for block in 1:n_blocks, group in 1:n_groups
        if is_in_neighborhood_per_block_per_group[block, group]
            push!(groups_per_block[block], group)
        end
    end

    group_index_per_membership = Int32[]
    cell_index_per_membership = Int32[]
    for cell in 1:n_cells
        if !is_participating_per_cell[cell]
            continue
        end
        block = block_index_per_cell[cell]
        for group in groups_per_block[block]
            push!(group_index_per_membership, group)
            push!(cell_index_per_membership, cell)
        end
    end

    return sparse(
        group_index_per_membership,
        cell_index_per_membership,
        trues(length(group_index_per_membership)),
        n_groups,
        n_cells,
    )
end

# Correlate, per gene, each cell's log expression fraction against its reference's (metacell's) log expression fraction,
# over arbitrary cell groups, filling `correlation_per_gene_per_group` (rows indexed by the cell gene, prefilled with
# zeros). Each cell belongs to the groups given by `is_member_per_group_per_cell` (a group X cell boolean matrix, one
# CSC column per cell) and maps to a reference via `reference_index_per_cell`. This is a sparse-aware, gene-parallel
# equivalent of gathering each group's dense per-cell vectors and calling `zero_cor_between_vectors`: it walks each
# gene's cell column once, routes the nonzero ("active") cells to their groups, and folds the constant zero-UMI
# background (every such cell has log fraction `log2(regularization)`) into the Pearson sufficient statistics in closed
# form, using each group's cell count and its reference distribution. The result matches the dense correlation up to
# floating point (the sums are reassociated). Groups with fewer than `min_group_cells` cells are left as zero.
function sparse_cell_reference_correlation_per_gene_per_group!(
    correlation_per_gene_per_group::AbstractMatrix{Float32},
    is_member_per_group_per_cell::SparseMatrixCSC{Bool},
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC,
    total_UMIs_per_cell::AbstractVector{<:Integer},
    reference_index_per_cell::AbstractVector{<:Integer},
    UMIs_per_reference_per_gene::AbstractMatrix,
    total_UMIs_per_reference::AbstractVector{<:Integer},
    cell_gene_index_per_correlated_gene::AbstractVector{<:Integer},
    reference_gene_index_per_correlated_gene::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat;
    min_group_cells::Integer = 2,
    name::AbstractString = "sparse_cell_reference_correlation_per_gene_per_group",
)::Nothing
    gene_fraction_regularization = Float32(gene_fraction_regularization)
    n_groups = size(is_member_per_group_per_cell, 1)
    n_cells = size(is_member_per_group_per_cell, 2)
    n_references = length(total_UMIs_per_reference)
    n_correlated_genes = length(cell_gene_index_per_correlated_gene)
    @assert length(reference_gene_index_per_correlated_gene) == n_correlated_genes

    group_index_per_membership = rowvals(is_member_per_group_per_cell)

    # Per-group cell count and reference distribution: `count_per_reference_per_group` (references X groups) sums, for
    # each group, how many of its cells map to each reference; a group's Sy, Syy over the constant-per-reference y are
    # then a single CSC column walk. Duplicate (reference, group) pairs are summed by `sparse`.
    n_cells_per_group = zeros(Int, n_groups)
    reference_index_per_group_membership = Int32[]
    group_index_per_group_membership = Int32[]
    for cell in 1:n_cells
        reference_index = reference_index_per_cell[cell]
        for membership in nzrange(is_member_per_group_per_cell, cell)
            group_index = group_index_per_membership[membership]
            n_cells_per_group[group_index] += 1
            push!(reference_index_per_group_membership, reference_index)
            push!(group_index_per_group_membership, group_index)
        end
    end
    count_per_reference_per_group = sparse(
        reference_index_per_group_membership,
        group_index_per_group_membership,
        ones(Int32, length(reference_index_per_group_membership)),
        n_references,
        n_groups,
    )
    reference_index_per_count = rowvals(count_per_reference_per_group)
    count_per_count = nonzeros(count_per_reference_per_group)

    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)

    background_cell_log_fraction = Float64(log2(gene_fraction_regularization))

    active_count_per_group_per_thread = [zeros(Int, n_groups) for _ in 1:maxthreadid()]
    sum_cell_per_group_per_thread = [zeros(Float64, n_groups) for _ in 1:maxthreadid()]
    sum_cell_squared_per_group_per_thread = [zeros(Float64, n_groups) for _ in 1:maxthreadid()]
    sum_cell_reference_per_group_per_thread = [zeros(Float64, n_groups) for _ in 1:maxthreadid()]
    sum_reference_per_group_per_thread = [zeros(Float64, n_groups) for _ in 1:maxthreadid()]
    touched_groups_per_thread = [Int[] for _ in 1:maxthreadid()]
    reference_log_fraction_per_reference_per_thread = [zeros(Float64, n_references) for _ in 1:maxthreadid()]

    progress = DebugProgress(n_correlated_genes; group = :mcs_loops, desc = name)
    parallel_loop_wo_rng(1:n_correlated_genes; name, progress) do correlated_gene_index
        cell_gene_index = cell_gene_index_per_correlated_gene[correlated_gene_index]
        reference_gene_index = reference_gene_index_per_correlated_gene[correlated_gene_index]

        active_count_per_group = active_count_per_group_per_thread[threadid()]
        sum_cell_per_group = sum_cell_per_group_per_thread[threadid()]
        sum_cell_squared_per_group = sum_cell_squared_per_group_per_thread[threadid()]
        sum_cell_reference_per_group = sum_cell_reference_per_group_per_thread[threadid()]
        sum_reference_per_group = sum_reference_per_group_per_thread[threadid()]
        touched_groups = touched_groups_per_thread[threadid()]
        reference_log_fraction_per_reference = reference_log_fraction_per_reference_per_thread[threadid()]
        empty!(touched_groups)

        # The reference (metacell) log fraction of this gene is constant per reference.
        @inbounds for reference_index in 1:n_references
            reference_log_fraction_per_reference[reference_index] = log2(
                UMIs_per_reference_per_gene[reference_index, reference_gene_index] /
                total_UMIs_per_reference[reference_index] + gene_fraction_regularization,
            )
        end

        # Walk this gene's cell column once; route each active cell to its groups.
        @inbounds for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, cell_gene_index)
            cell_UMIs = UMIs_per_stored[stored_index]
            if cell_UMIs <= 0
                continue
            end
            cell_index = row_index_per_stored[stored_index]
            reference_index = reference_index_per_cell[cell_index]
            if reference_index < 1  # A non-participating cell belongs to no group; skip (and avoid a zero reference).
                continue
            end
            cell_log_fraction = log2(cell_UMIs / total_UMIs_per_cell[cell_index] + gene_fraction_regularization)
            reference_log_fraction = reference_log_fraction_per_reference[reference_index]
            for membership in nzrange(is_member_per_group_per_cell, cell_index)
                group_index = group_index_per_membership[membership]
                if active_count_per_group[group_index] == 0
                    push!(touched_groups, group_index)
                end
                active_count_per_group[group_index] += 1
                sum_cell_per_group[group_index] += cell_log_fraction
                sum_cell_squared_per_group[group_index] += cell_log_fraction * cell_log_fraction
                sum_cell_reference_per_group[group_index] += cell_log_fraction * reference_log_fraction
                sum_reference_per_group[group_index] += reference_log_fraction
            end
        end

        # Finalize each touched group's Pearson correlation from the active sums plus the constant background, then
        # reset its scratch for the next gene on this thread.
        @inbounds for group_index in touched_groups
            n_group_cells = n_cells_per_group[group_index]
            active_count = active_count_per_group[group_index]
            if n_group_cells >= min_group_cells
                sum_reference = 0.0
                sum_reference_squared = 0.0
                for count_index in nzrange(count_per_reference_per_group, group_index)
                    reference_index = reference_index_per_count[count_index]
                    count = count_per_count[count_index]
                    reference_log_fraction = reference_log_fraction_per_reference[reference_index]
                    sum_reference += count * reference_log_fraction
                    sum_reference_squared += count * reference_log_fraction * reference_log_fraction
                end
                n_background_cells = n_group_cells - active_count
                sum_cell = background_cell_log_fraction * n_background_cells + sum_cell_per_group[group_index]
                sum_cell_squared =
                    background_cell_log_fraction * background_cell_log_fraction * n_background_cells +
                    sum_cell_squared_per_group[group_index]
                sum_cell_reference =
                    sum_cell_reference_per_group[group_index] +
                    background_cell_log_fraction * (sum_reference - sum_reference_per_group[group_index])
                cell_variance = n_group_cells * sum_cell_squared - sum_cell * sum_cell
                reference_variance = n_group_cells * sum_reference_squared - sum_reference * sum_reference
                # var >= 0 mathematically (equality iff constant); a relative guard rejects the floating point residual
                # of constant data (matching a two-pass zero correlation).
                if cell_variance > 1e-9 * n_group_cells * sum_cell_squared &&
                   reference_variance > 1e-9 * n_group_cells * sum_reference_squared
                    correlation_per_gene_per_group[cell_gene_index, group_index] =
                        (n_group_cells * sum_cell_reference - sum_cell * sum_reference) /
                        sqrt(cell_variance * reference_variance)
                end
            end
            active_count_per_group[group_index] = 0
            sum_cell_per_group[group_index] = 0.0
            sum_cell_squared_per_group[group_index] = 0.0
            sum_cell_reference_per_group[group_index] = 0.0
            sum_reference_per_group[group_index] = 0.0
        end
        return nothing
    end

    return nothing
end

# Correlate, per gene, each cell's log expression fraction against its own metacell's *punctuated* (leave-one-out) log
# expression fraction, over each target block's neighborhood, filling `correlation_per_gene_per_block` (rows indexed by
# the gene, prefilled with zeros). Punctuation - `(metacell_UMIs - cell_UMIs) / (metacell_total - cell_total)` - makes
# the metacell fraction vary per cell, so the constant-per-reference background collapse of
# `sparse_cell_reference_correlation_per_gene_per_group!` does not apply. Instead this exploits that each cell's `x` and
# `y` are per-(cell, gene) - independent of the target block - and each target block's neighborhood is a union of whole
# own blocks (a cell is in block `b`'s neighborhood iff its own block is in `b`'s neighborhood). So it computes each
# participating cell's `(x, y)` once per gene, accumulates the Pearson sufficient statistics per own block, then combines
# each target block's neighborhood own blocks - exact, and far cheaper than gathering each block's dense vectors. A cell
# participates iff `own_block_index_per_cell[cell] > 0`. Two numeric guards keep the one-pass sufficient statistics
# stable: the log fractions are shifted by `log2(regularization)` (so zero-UMI background cells contribute exactly 0),
# and a relative-variance guard rejects the floating point residual of constant data (matching a two-pass zero
# correlation). Target blocks with fewer than `min_block_cells` neighborhood cells are left as zero. When
# `is_punctuated_per_cell` is given, a cell's metacell fraction is punctuated only where its flag is set (elsewhere it
# uses the plain `metacell_UMIs / metacell_total`), so `metacell_index_per_cell` need not be the cell's own metacell.
function own_block_punctuated_correlation_per_gene_per_block!(
    correlation_per_gene_per_block::AbstractMatrix{Float32},
    own_block_index_per_cell::AbstractVector{<:Integer},
    is_in_neighborhood_per_own_block_per_target_block::Union{AbstractMatrix{Bool}, BitMatrix},
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC,
    total_UMIs_per_cell::AbstractVector{<:Integer},
    metacell_index_per_cell::AbstractVector{<:Integer},
    UMIs_per_metacell_per_gene::AbstractMatrix,
    total_UMIs_per_metacell::AbstractVector{<:Integer},
    gene_index_per_correlated_gene::AbstractVector{<:Integer},
    gene_fraction_regularization::AbstractFloat;
    is_punctuated_per_cell::Maybe{Union{AbstractVector{Bool}, BitVector}} = nothing,
    min_block_cells::Integer = 2,
    name::AbstractString = "own_block_punctuated_correlation_per_gene_per_block",
)::Nothing
    gene_fraction_regularization = Float32(gene_fraction_regularization)
    n_cells = length(own_block_index_per_cell)
    n_own_blocks, n_target_blocks = size(is_in_neighborhood_per_own_block_per_target_block)
    n_correlated_genes = length(gene_index_per_correlated_gene)
    background_log_fraction = Float64(log2(gene_fraction_regularization))

    own_blocks_per_target_block =
        [findall(@view is_in_neighborhood_per_own_block_per_target_block[:, target]) for target in 1:n_target_blocks]
    participating_cell_indices = findall(>(0), own_block_index_per_cell)

    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)

    gene_UMIs_per_cell_per_thread = [zeros(eltype(UMIs_per_stored), n_cells) for _ in 1:maxthreadid()]
    sum_cell_per_own_block_per_thread = [zeros(Float64, n_own_blocks) for _ in 1:maxthreadid()]
    sum_cell_squared_per_own_block_per_thread = [zeros(Float64, n_own_blocks) for _ in 1:maxthreadid()]
    sum_metacell_per_own_block_per_thread = [zeros(Float64, n_own_blocks) for _ in 1:maxthreadid()]
    sum_metacell_squared_per_own_block_per_thread = [zeros(Float64, n_own_blocks) for _ in 1:maxthreadid()]
    sum_cell_metacell_per_own_block_per_thread = [zeros(Float64, n_own_blocks) for _ in 1:maxthreadid()]
    n_cells_per_own_block_per_thread = [zeros(Int, n_own_blocks) for _ in 1:maxthreadid()]

    progress = DebugProgress(n_correlated_genes; group = :mcs_loops, desc = name)
    parallel_loop_wo_rng(1:n_correlated_genes; name, progress) do correlated_gene_index
        gene_index = gene_index_per_correlated_gene[correlated_gene_index]
        gene_UMIs_per_cell = gene_UMIs_per_cell_per_thread[threadid()]
        sum_cell = sum_cell_per_own_block_per_thread[threadid()]
        sum_cell_squared = sum_cell_squared_per_own_block_per_thread[threadid()]
        sum_metacell = sum_metacell_per_own_block_per_thread[threadid()]
        sum_metacell_squared = sum_metacell_squared_per_own_block_per_thread[threadid()]
        sum_cell_metacell = sum_cell_metacell_per_own_block_per_thread[threadid()]
        n_cells_per_own_block = n_cells_per_own_block_per_thread[threadid()]
        fill!(sum_cell, 0.0)
        fill!(sum_cell_squared, 0.0)
        fill!(sum_metacell, 0.0)
        fill!(sum_metacell_squared, 0.0)
        fill!(sum_cell_metacell, 0.0)
        fill!(n_cells_per_own_block, 0)

        @inbounds for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
            gene_UMIs_per_cell[row_index_per_stored[stored_index]] = UMIs_per_stored[stored_index]
        end

        # Accumulate each participating cell's Pearson sufficient statistics into its own block.
        @inbounds for cell_index in participating_cell_indices
            metacell_index = metacell_index_per_cell[cell_index]
            if metacell_index < 1
                continue
            end
            own_block = own_block_index_per_cell[cell_index]
            cell_total_UMIs = total_UMIs_per_cell[cell_index]
            cell_UMIs = gene_UMIs_per_cell[cell_index]
            metacell_UMIs = UMIs_per_metacell_per_gene[metacell_index, gene_index]
            cell_log_fraction =
                log2(cell_UMIs / cell_total_UMIs + gene_fraction_regularization) - background_log_fraction
            if is_punctuated_per_cell === nothing || is_punctuated_per_cell[cell_index]
                metacell_log_fraction =
                    log2(
                        (metacell_UMIs - cell_UMIs) / (total_UMIs_per_metacell[metacell_index] - cell_total_UMIs) +
                        gene_fraction_regularization,
                    ) - background_log_fraction
            else
                metacell_log_fraction =
                    log2(metacell_UMIs / total_UMIs_per_metacell[metacell_index] + gene_fraction_regularization) -
                    background_log_fraction
            end
            sum_cell[own_block] += cell_log_fraction
            sum_cell_squared[own_block] += cell_log_fraction * cell_log_fraction
            sum_metacell[own_block] += metacell_log_fraction
            sum_metacell_squared[own_block] += metacell_log_fraction * metacell_log_fraction
            sum_cell_metacell[own_block] += cell_log_fraction * metacell_log_fraction
            n_cells_per_own_block[own_block] += 1
        end

        # Combine each target block's neighborhood own blocks and finalize the Pearson correlation.
        @inbounds for target_block in 1:n_target_blocks
            combined_cell = 0.0
            combined_cell_squared = 0.0
            combined_metacell = 0.0
            combined_metacell_squared = 0.0
            combined_cell_metacell = 0.0
            combined_n_cells = 0
            for own_block in own_blocks_per_target_block[target_block]
                combined_cell += sum_cell[own_block]
                combined_cell_squared += sum_cell_squared[own_block]
                combined_metacell += sum_metacell[own_block]
                combined_metacell_squared += sum_metacell_squared[own_block]
                combined_cell_metacell += sum_cell_metacell[own_block]
                combined_n_cells += n_cells_per_own_block[own_block]
            end
            if combined_n_cells >= min_block_cells
                cell_variance = combined_n_cells * combined_cell_squared - combined_cell * combined_cell
                metacell_variance = combined_n_cells * combined_metacell_squared - combined_metacell * combined_metacell
                if cell_variance > 1e-9 * combined_n_cells * combined_cell_squared &&
                   metacell_variance > 1e-9 * combined_n_cells * combined_metacell_squared
                    correlation_per_gene_per_block[gene_index, target_block] =
                        (combined_n_cells * combined_cell_metacell - combined_cell * combined_metacell) /
                        sqrt(cell_variance * metacell_variance)
                end
            end
        end

        @inbounds for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
            gene_UMIs_per_cell[row_index_per_stored[stored_index]] = 0
        end
        return nothing
    end

    return nothing
end

# Sum a sparse cell X gene UMIs matrix over gene groups, per target cell, accumulating into
# `summed_UMIs_per_group_per_cell` (rows indexed by the group, columns by the target cell, prefilled with zeros).
# `sparse_UMIs_per_cell_per_gene` is a `SparseMatrixCSC` whose columns are the genes; each group's gene columns are
# walked once - touching only stored non-zeros - and scattered into their cells, rather than gathering a dense per-group
# submatrix. The target cell of output column `position` is `cell_index_per_position[position]`; a gene contributes to a
# column only for cells in this set. `cell_position_per_cell` is caller-provided scratch of length `size(sparse, 1)`, all
# zero on entry and left all zero on exit, so it can be a reused per-thread buffer.
function sum_sparse_UMIs_per_gene_group_per_cell!(
    summed_UMIs_per_group_per_cell::AbstractMatrix{<:AbstractFloat},
    sparse_UMIs_per_cell_per_gene::SparseMatrixCSC,
    cell_index_per_position::AbstractVector{<:Integer},
    gene_indices_per_group::AbstractVector{<:AbstractVector{<:Integer}},
    cell_position_per_cell::AbstractVector{<:Integer},
)::Nothing
    n_positions = length(cell_index_per_position)
    row_index_per_stored = rowvals(sparse_UMIs_per_cell_per_gene)
    UMIs_per_stored = nonzeros(sparse_UMIs_per_cell_per_gene)

    @inbounds for position in 1:n_positions
        cell_position_per_cell[cell_index_per_position[position]] = position
    end

    @inbounds for (group_index, gene_indices_of_group) in enumerate(gene_indices_per_group)
        for gene_index in gene_indices_of_group
            for stored_index in nzrange(sparse_UMIs_per_cell_per_gene, gene_index)
                position = cell_position_per_cell[row_index_per_stored[stored_index]]
                if position > 0
                    summed_UMIs_per_group_per_cell[group_index, position] += UMIs_per_stored[stored_index]
                end
            end
        end
    end

    @inbounds for position in 1:n_positions
        cell_position_per_cell[cell_index_per_position[position]] = 0
    end

    return nothing
end

end  # module
