using Test

using Daf
using Metacells
using NestedTests

test_prefixes(ARGS)
abort_on_first_failure(true)

include("gmara.jl")
