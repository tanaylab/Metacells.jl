using Test

using Metacells
using NestedTests

test_prefixes(ARGS)
abort_on_first_failure(true)

nested_test("nothing") do
    @test 1 == 1
end
