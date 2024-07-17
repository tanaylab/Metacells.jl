nested_test("programs") do
    n_metacells = 5
    n_genes = 10
    n_factors = 3
    coefficients_of_factors_of_genes = Matrix(sprand(n_factors, n_genes, 0.2))
    values_of_factors_of_metacells = rand(n_metacells, n_factors)
    values_of_genes_of_metacells = values_of_factors_of_metacells * coefficients_of_factors_of_genes
    @assert size(values_of_genes_of_metacells) == (n_metacells, n_genes)

    nested_test("exact") do
        solution = Metacells.Programs.solve(
            Metacells.Programs.Problem(values_of_factors_of_metacells, values_of_genes_of_metacells, nothing);
            regularization = 0.0,
            min_abs_coefficient = 1e-3,
        )
        mean_abs_error = mean(abs.(solution.coefficients_of_factors_of_genes .- coefficients_of_factors_of_genes))
        @test mean_abs_error < 0.2
    end

    nested_test("warmstart") do
        start_of_factors_of_genes = coefficients_of_factors_of_genes
        solution = Metacells.Programs.solve(
            Metacells.Programs.Problem(
                values_of_factors_of_metacells,
                values_of_genes_of_metacells,
                start_of_factors_of_genes,
            );
            regularization = 0.0,
            min_abs_coefficient = 1e-3,
        )
        mean_abs_error = mean(abs.(solution.coefficients_of_factors_of_genes .- coefficients_of_factors_of_genes))
        @test mean_abs_error < 0.2
    end

    nested_test("approximate") do
        values_of_genes_of_metacells .+= rand(n_metacells, n_genes) ./ 10
        solution = Metacells.Programs.solve(
            Metacells.Programs.Problem(values_of_factors_of_metacells, values_of_genes_of_metacells, nothing);
            regularization = 0.01,
            min_abs_coefficient = 1e-3,
        )
        mean_abs_error = mean(abs.(solution.coefficients_of_factors_of_genes .- coefficients_of_factors_of_genes))
        @test mean_abs_error < 0.2
    end
end
