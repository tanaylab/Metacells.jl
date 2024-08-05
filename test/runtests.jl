using Test

using Daf
using Metacells
using NestedTests
using Random
using SparseArrays
using StatsBase

Random.seed!(123456)

test_prefixes(ARGS)
abort_on_first_failure(true)

include("gmara.jl")
