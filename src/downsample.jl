"""
Downsampling of data.
"""
module Downsample

export downsample
export downsamples

using Base.Threads
using DataAxesFormats
using Random
using Statistics
using TanayLabUtilities

import Random.default_rng

"""
    downsample(
        vector::AbstractVector{<:Integer},
        samples::Integer;
        rng::AbstractRNG = default_rng(),
        output::Maybe{AbstractVector} = nothing,
    )::AbstractVector

    downsample(
        matrix::AbstractMatrix{<:Integer},
        samples::Integer;
        dims::Integer,
        rng::AbstractRNG = default_rng(),
        output::Maybe{AbstractMatrix} = nothing,
    )::AbstractMatrix

Given a `vector` of integer non-negative data values, return a new vector such that the sum of entries in it is
`samples`. Think of the original vector as containing a number of marbles in each entry. We pick `samples` marbles from
this vector; each time we pick a marble we take it out of the original vector and move it to the result.

If the sum of the entries of a vector is less than `samples`, it is copied to the output. If `output` is not specified,
it is allocated automatically using the same element type as the input.

When downsampling a `matrix`, then `dims` must be specified to be `1`/`Rows` (downsample each row vector) or
`2`/`Columns` (downsample each column vector).

We need to downsample before correlating single-cell gene UMI count profiles. If we don't, then less sparse cells (with
more total UMIs) will appear to be more correlated with everything, simply due to chance. Yes, downsampling does
increase the noise level of the correlation results, but that's a price we must pay for being able to meaningfully
compare these correlations.

Correlating metacell gene fraction profiles doesn't require downsampling as these fractions are dense with robust
estimations for all the genes, so the correlation between metacells isn't sensitive to the total number of UMIs in each
metacell.
"""
function downsample(
    vector::AbstractVector{<:Integer},
    samples::Integer;
    rng::AbstractRNG = default_rng(),
    output::Maybe{AbstractVector} = nothing,
)::AbstractVector
    n_values = length(vector)
    if output === nothing
        output = Vector{eltype(vector)}(undef, n_values)
    else
        @assert_vector(output, n_values)
    end

    if n_values > 0
        @assert minimum(vector) >= 0 "Downsampling a vector with negative values"
    end

    if n_values == 1
        output[1] = min(samples, vector[1])

    elseif n_values > 1
        tree = initialize_tree(vector)

        if tree[end] <= samples
            output .= vector

        else
            output .= 0
            for _ in 1:samples
                output[random_sample(tree, rand(rng, 1:tree[end]))] += 1
            end
        end
    end

    return output
end

function downsample(
    matrix::AbstractMatrix{<:Integer},
    samples::Integer;
    dims::Integer,
    rng::AbstractRNG = default_rng(),
    output::Maybe{AbstractMatrix} = nothing,
)::AbstractMatrix
    @assert 1 <= dims <= 2
    n_rows, n_columns = size(matrix)
    @assert_matrix(matrix, n_rows, n_columns, dims)

    if output === nothing
        if dims == Rows
            output = transpose(Matrix{eltype(matrix)}(undef, n_columns, n_rows))
        elseif dims == Columns
            output = Matrix{eltype(matrix)}(undef, n_rows, n_columns)
        else
            @assert false
        end
    end
    @assert_matrix(output, n_rows, n_columns, dims)

    if dims == Rows
        parallel_loop_with_rng(n_rows; rng) do row_index, rng
            @views row_vector = matrix[row_index, :]
            @views output_vector = output[row_index, :]
            downsample(row_vector, samples; rng, output = output_vector)
        end
    elseif dims == Columns
        parallel_loop_with_rng(n_columns; rng) do column_index, rng
            @views column_vector = matrix[:, column_index]
            @views output_vector = output[:, column_index]
            downsample(column_vector, samples; rng, output = output_vector)
        end
    else
        @assert false
    end

    return output
end

function initialize_tree(input::AbstractVector{T})::AbstractVector{T} where {T <: Integer}
    n_values = length(input)
    @assert n_values > 1

    n_values_in_level = ceil_power_of_two(n_values)
    tree_size = 2 * n_values_in_level - 1

    tree = Vector{T}(undef, tree_size)

    tree[1:n_values] .= input
    tree[(n_values + 1):end] .= 0

    tree_of_level = tree

    while (n_values_in_level > 1)
        @assert iseven(n_values_in_level)

        @views input_of_level = tree_of_level[1:n_values_in_level]
        @views tree_of_level = tree_of_level[(n_values_in_level + 1):end]
        n_values_in_level = div(n_values_in_level, 2)

        @assert length(tree_of_level) >= n_values_in_level

        for index_in_level in 1:n_values_in_level
            left_value = input_of_level[index_in_level * 2 - 1]
            right_value = input_of_level[index_in_level * 2]
            tree_of_level[index_in_level] = left_value + right_value
        end
    end

    @assert length(tree_of_level) == 1

    return tree
end

function tree_length(size::Integer)::Integer
    return 2 * ceil_power_of_two(size) - 1
end

function ceil_power_of_two(size::Integer)::Integer
    return 2^Int(ceil(log2(size)))
end

function random_sample(tree::AbstractVector{<:Integer}, random::Integer)::Integer
    size_of_level = 1
    base_of_level = length(tree)

    index_in_level = 1
    index_in_tree = base_of_level + index_in_level - 1

    while true
        @assert tree[index_in_tree] > 0
        tree[index_in_tree] -= 1

        size_of_level *= 2
        base_of_level -= size_of_level

        if base_of_level <= 0
            return index_in_level
        end

        index_in_level = index_in_level * 2 - 1
        index_in_tree = base_of_level + index_in_level - 1
        right_random = random - tree[index_in_tree]

        if right_random > 0
            index_in_level += 1
            index_in_tree += 1
            random = right_random
        end
    end
end

"""
    downsamples(
        samples_vector::AbstractVector{<:Integer};
        min_downsamples::Integer = _TODOX(DEFAULT.min_downsamples),
        min_downsamples_quantile::AbstractFloat = _TODOX(DEFAULT.min_downsamples_quantile),
        max_downsamples_quantile::AbstractFloat = _TODOX(DEFAULT.max_downsamples_quantile),
    )::Integer

When downsampling multiple vectors (the amount of data in each available in `samples_vector`), we need to pick a
"reasonable" number of samples to downsample to. We have conflicting requirements, so this is a compromise. First, we
want most vectors to have at least the target number of samples, so we start with the `min_downsamples_quantile` of the
`samples_vector`. Second, we also want to have at least `min_downsamples` to ensure we don't throw away too much data
even if most vectors are sparse, so we increase the target to this value. Finally, we don't want a target which is too
big for too many vectors, so so we reduce the result to the `max_downsamples_quantile` of the `samples_vector`.

The defaults were chosen to
"""
function downsamples(
    samples_vector::AbstractVector{<:Integer};
    min_downsamples::Integer = 750,
    min_downsamples_quantile::AbstractFloat = 0.05,
    max_downsamples_quantile::AbstractFloat = 0.5,
)::Integer
    @assert 0 <= min_downsamples_quantile <= max_downsamples_quantile <= 1
    return Int(
        round(
            min(
                max(min_downsamples, quantile(samples_vector, min_downsamples_quantile)),
                quantile(samples_vector, max_downsamples_quantile),
            ),
        ),
    )
end

end
