nested_test("gmara") do
    nested_test("normalize_gene_name") do
        nested_test("UCSC") do
            @test normalize_gene_name("UC12345.1"; namespace = "UCSC") == "uc12345.1"
        end

        nested_test("Symbol") do
            @test normalize_gene_name("Abc.1"; namespace = "Symbol") == "ABC"
        end
    end

    nested_test("namespace") do
        mktempdir() do cache_dir
            symbols = gmara_genes(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            @test "SOX4" in symbols
            @test "MT-LSP" in symbols
            @test !("AC008770" in symbols)
            @assert symbols === gmara_genes(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            empty_gmara_cache!()
            @assert symbols !== gmara_genes(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            @assert symbols == gmara_genes(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
        end
    end

    nested_test("list") do
        mktempdir() do cache_dir
            transcription_factors = gmara_genes(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            @test "SOX4" in transcription_factors
            @test !("MT-LSP" in transcription_factors)
            @test !("AC008770" in transcription_factors)
            @assert transcription_factors == gmara_genes(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            empty_gmara_cache!()
            @assert transcription_factors !== gmara_genes(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            @assert transcription_factors == gmara_genes(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
        end
    end

    nested_test("mask") do
        mktempdir() do cache_dir
            daf = MemoryDaf()
            add_axis!(daf, "gene", ["SOX4", "MT-LSP", "AC008770"])
            set_gmara_genes_mask!(daf; species = "human", namespace = "Symbol", list = "transcription_factor")
            @test get_vector(daf, "gene", "is_transcription_factor").array == [true, false, false]
        end
    end
end
