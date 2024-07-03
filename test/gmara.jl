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
            symbols = gmara_namespace(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            @assert symbols === gmara_namespace(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            empty_gmara_cache!()
            @assert symbols !== gmara_namespace(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
            @assert symbols == gmara_namespace(; species = "human", namespace = "Symbol", cache_dir = cache_dir)
        end
    end

    nested_test("list") do
        mktempdir() do cache_dir
            transcription_factors = gmara_list(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factors",
                cache_dir = cache_dir,
            )
            @assert transcription_factors == gmara_list(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factors",
                cache_dir = cache_dir,
            )
            empty_gmara_cache!()
            @assert transcription_factors !== gmara_list(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factors",
                cache_dir = cache_dir,
            )
            @assert transcription_factors == gmara_list(;
                species = "human",
                namespace = "Symbol",
                list = "transcription_factors",
                cache_dir = cache_dir,
            )
        end
    end
end
