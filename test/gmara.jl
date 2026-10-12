nested_test("gmara") do
    nested_test("normalize_gene_name") do
        nested_test("UCSC") do
            @test normalize_gene_name("uc001abc.1"; namespace = "UCSC") == "uc001abc.1"
        end

        nested_test("GeneSymbol") do
            @test normalize_gene_name("AC092171.3"; namespace = "GeneSymbol") == "AC092171.3"
        end

        nested_test("unknown") do
            @test_throws "unknown Gmara namespace: EnsembleGene" normalize_gene_name(
                "ENSG00000124766";
                namespace = "EnsembleGene",
            )
        end

        nested_test("EnsemblGene") do
            nested_test("version") do
                @test normalize_gene_name("ENSG00000124766.5"; namespace = "EnsemblGene") == "ENSG00000124766"
            end

            nested_test("junk") do
                @test normalize_gene_name("ENSG00000124766;"; namespace = "EnsemblGene") == "ENSG00000124766"
            end

            nested_test("invalid") do
                @test normalize_gene_name("ERCC-00002"; namespace = "EnsemblGene") == "ERCC-00002"
            end

            nested_test("validate") do
                @test_throws "invalid EnsemblGene gene name: ERCC-00002" normalize_gene_name(
                    "ERCC-00002";
                    namespace = "EnsemblGene",
                    validate = true,
                )
            end
        end
    end

    nested_test("namespace") do
        mktempdir() do cache_dir
            symbols = gmara_genes(; species = "human", namespace = "GeneSymbol", cache_dir = cache_dir)
            @test "SOX4" in symbols
            @test "MT-LSP" in symbols
            @test !("AC008770" in symbols)
            @assert symbols === gmara_genes(; species = "human", namespace = "GeneSymbol", cache_dir = cache_dir)
            empty_gmara_cache!()
            @assert symbols !== gmara_genes(; species = "human", namespace = "GeneSymbol", cache_dir = cache_dir)
            @assert symbols == gmara_genes(; species = "human", namespace = "GeneSymbol", cache_dir = cache_dir)
        end
    end

    nested_test("list") do
        mktempdir() do cache_dir
            transcription_factors = gmara_genes(;
                species = "human",
                namespace = "GeneSymbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            @test "SOX4" in transcription_factors
            @test !("MT-LSP" in transcription_factors)
            @test !("AC008770" in transcription_factors)
            @assert transcription_factors == gmara_genes(;
                species = "human",
                namespace = "GeneSymbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            empty_gmara_cache!()
            @assert transcription_factors !== gmara_genes(;
                species = "human",
                namespace = "GeneSymbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
            @assert transcription_factors == gmara_genes(;
                species = "human",
                namespace = "GeneSymbol",
                list = "transcription_factor",
                cache_dir = cache_dir,
            )
        end
    end

    nested_test("mask") do
        mktempdir() do cache_dir
            daf = MemoryDaf()
            add_axis!(daf, "gene", ["SOX4", "MT-LSP", "AC008770"])
            set_gmara_genes_mask!(daf; species = "human", namespace = "GeneSymbol", list = "transcription_factor")
            @test get_vector(daf, "gene", "is_transcription_factor").array == [true, false, false]
        end
    end
end
