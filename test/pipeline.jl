import Metacells.Pipeline.round_prefix

nested_test("pipeline") do
    nested_test("round_prefix") do
        nested_test("single") do
            @test round_prefix("B", 0) == "B"
            @test round_prefix("B", 1) == "C"
            @test round_prefix("B", 10) == "L"
            @test round_prefix("M", 1) == "N"
            @test round_prefix("M", 10) == "W"
        end

        nested_test("double") do
            @test round_prefix("B", 11) == "BB"
            @test round_prefix("B", 12) == "BC"
            @test round_prefix("B", 22) == "CB"
            @test round_prefix("B", 131) == "LL"
        end

        nested_test("triple") do
            @test round_prefix("B", 132) == "BBB"
        end

        nested_test("!letter") do
            @test_throws "the base prefix: BB is not a single letter" round_prefix("BB", 1)
        end

        nested_test("!room") do
            @test_throws "invalid base prefix: Q" round_prefix("Q", 1)
        end
    end

    nested_test("sharpening_rounds") do
        daf = MemoryDaf(; name = "memory!")

        nested_test("!overlap") do
            @test_throws "overlapping metacells prefix: F and blocks prefix: B" sharpening_rounds(;
                initial_daf = daf,
                base_daf = daf,
                score_daf = daf,
                directory = "unused",
                metacells_prefix = "F",
            )
        end

        nested_test("prefixes") do
            rounds = sharpening_rounds(; initial_daf = daf, base_daf = daf, score_daf = daf, directory = "unused")
            sharpening_round, _ = iterate(rounds)
            @test sharpening_round.index == 1
            @test sharpening_round.previous_daf === daf
            @test sharpening_round.metacells_prefix == "N"
            @test sharpening_round.blocks_prefix == "C"
        end

        nested_test("!run") do
            rounds = sharpening_rounds(; initial_daf = daf, base_daf = daf, score_daf = daf, directory = "unused")
            _, index = iterate(rounds)
            @test_throws "the sharpening round: 1 was not run" iterate(rounds, index)
        end
    end
end
