nested_test("downsample") do
    Random.seed!(123456)

    nested_test("empty") do
        values = Int32[]
        @test isempty(downsample(values, 1))
    end

    nested_test("single") do
        values = [10]
        nested_test("low") do
            @test downsample(values, 1) == [1]
        end

        nested_test("same") do
            @test downsample(values, 10) == [10]
        end

        nested_test("high") do
            @test downsample(values, 20) == [10]
        end
    end

    nested_test("many") do
        values = [0, 1, 2, 3, 4, 5]
        downsampled = downsample(values, 10)
        @test downsampled == [0, 1, 2, 1, 3, 3]
    end

    nested_test("columns") do
        values = Matrix(transpose([0 1 2 3 4 5; 0 1 2 3 4 5]))
        downsampled = downsample(values, 10; dims = Columns)
        @test downsampled == transpose([0 1 1 2 3 3; 0 0 2 2 3 3])
    end

    nested_test("rows") do
        values = transpose(transposer([0 1 2 3 4 5; 0 1 2 3 4 5]))
        downsampled = downsample(values, 10; dims = Rows)
        @test downsampled == [0 1 1 2 3 3; 0 0 2 2 3 3]
    end

    nested_test("samples") do
        nested_test("low") do
            samples = Vector{Int32}(undef, 20)
            samples .= 1000
            samples[1:2] .= 800
            @test downsamples(samples) == 800
        end

        nested_test("middle") do
            samples = Vector{Int32}(undef, 20)
            samples .= 1000
            samples[1:2] .= 100
            @test downsamples(samples) == 750
        end

        nested_test("high") do
            samples = Vector{Int32}(undef, 20)
            samples .= 500
            samples[1:2] .= 80
            @test downsamples(samples) == 500
        end
    end
end
