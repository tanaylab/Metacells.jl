using Test

using DataAxesFormats
using Metacells
using NestedTests
using Random
using SparseArrays
using StatsBase
using TanayLabUtilities

Random.seed!(123456)

test_prefixes(ARGS)
abort_on_first_failure(true)

include("anndata_format.jl")
include("gmara.jl")
include("pipeline.jl")
