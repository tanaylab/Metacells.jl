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

            nested_test("generic") do
                @test normalize_gene_name(Test.GenericString("ENSG00000124766.5"); namespace = "EnsemblGene") ==
                      "ENSG00000124766"
            end

            nested_test("junk") do
                @test normalize_gene_name("ENSG00000124766;"; namespace = "EnsemblGene") == "ENSG00000124766"
            end

            nested_test("par") do
                @test normalize_gene_name("ENSG00000182378.15_PAR_Y;"; namespace = "EnsemblGene") == "ENSG00000182378"
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
            set_vector!(daf, "gene", "ensembl", ["ENSG00000124766.5", "", "ERCC-00002"])

            nested_test("name") do
                set_gmara_genes_mask!(daf; species = "human", namespace = "GeneSymbol", list = "transcription_factor")
                @test get_vector(daf, "gene", "is_transcription_factor").array == [true, false, false]
            end

            nested_test("ensembl") do
                set_gmara_genes_mask!(
                    daf;
                    species = "human",
                    namespace = "EnsemblGene",
                    list = "transcription_factor",
                    mask_property = "is_ensembl_transcription_factor",
                )
                @test get_vector(daf, "gene", "is_ensembl_transcription_factor").array == [true, false, false]
            end

            nested_test("gmara_genes_mask") do
                nested_test("symbols") do
                    names = ["RPS9", "SOX4", "MT-CO1"]
                    @test gmara_genes_mask(; names, species = "human", list = "ribosomal") == [true, false, false]
                    @test gmara_genes_mask(; names, species = "human", list = "mitochondrial") == [false, false, true]
                end

                nested_test("ensembl") do
                    names = ["ENSG00000170889.12", "ENSG00000124766"]
                    @test gmara_genes_mask(; names, species = "human", namespace = "EnsemblGene", list = "ribosomal") ==
                          [true, false]
                end

                nested_test("none") do
                    @test gmara_genes_mask(; names = ["SOX4"], species = "human", list = "ribosomal") == [false]
                end
            end

            nested_test("!name_property") do
                @test_throws "no name_property specified for the namespace: HGNC" set_gmara_genes_mask!(
                    daf;
                    species = "human",
                    namespace = "HGNC",
                    list = "transcription_factor",
                )
            end

            nested_test("!mask_property") do
                message = "no mask_property specified for the names of the namespace: GeneSymbol"
                @test_throws message set_gmara_genes_mask!(daf; species = "human", namespace = "GeneSymbol")
            end
        end
    end

    nested_test("ensembl") do
        mktempdir() do cache_dir
            gmara_ensembl_genes = fetch_gmara_ensembl_genes(; species = "human", cache_dir = cache_dir)

            nested_test("fetch") do
                @test gmara_ensembl_genes.ensembls_of_name["SOX4"] == Set(["ENSG00000124766"])
                @test gmara_ensembl_genes.final_ensembls_of_ensembl["ENSG00000124766"] == Set(["ENSG00000124766"])
                @test gmara_ensembl_genes.final_ensembls_of_ensembl["ENSG00000263255"] == Set(["ENSG00000256229"])
                @test fetch_gmara_ensembl_genes(; species = "human", cache_dir = cache_dir).ensembls_of_name ===
                      gmara_ensembl_genes.ensembls_of_name
            end

            nested_test("collect") do
                @test collect_ensembl_genes(;
                    names = ["SOX4", "MTRNR2L6", "NO-SUCH-GENE"],
                    gmara_ensembl_genes = gmara_ensembl_genes,
                ) == ["ENSG00000124766", "", ""]
            end

            nested_test("consistency") do
                consistencies = verify_consistent_ensembl_genes(;
                    names = ["SOX4", "SOX4", "PLCXD1", "SOX4", "RPS9", "", "NO-SUCH-GENE", "SOX4"],
                    ensembl_per_name = [
                        "ENSG00000124766",
                        "ENSG00000124766.5",
                        "ENSG00000182378.15_PAR_Y",
                        "ENSG00000263255.2",
                        "ENSG00000124766",
                        "",
                        "ERCC-00002",
                        "ENSG99999999999",
                    ],
                    gmara_ensembl_genes = gmara_ensembl_genes,
                )
                properties = propertynames(consistencies[1])
                @test properties[1:7] == (
                    :name,
                    :ensembl,
                    :normalized_ensembl,
                    :name_ensembls,
                    :final_ensembls,
                    :is_consistent,
                    :ensembl_is_normalized,
                )
                true_properties = [
                    [property for property in properties[6:end] if getproperty(consistency, property)] for
                    consistency in consistencies
                ]
                @test true_properties == [
                    [:is_consistent, :ensembl_is_normalized],
                    [:is_consistent],
                    [:is_consistent],
                    [:ensembl_is_retired, :ensembl_is_different],
                    [:ensembl_is_normalized, :ensembl_is_different],
                    [:name_is_empty, :ensembl_is_empty],
                    [:name_is_unknown, :ensembl_is_invalid],
                    [:ensembl_is_normalized, :ensembl_is_unknown],
                ]
                @test consistencies[4].final_ensembls == Set(["ENSG00000256229"])
            end

            nested_test("generic") do
                consistencies = verify_consistent_ensembl_genes(;
                    names = [Test.GenericString("SOX4")],
                    ensembl_per_name = [Test.GenericString("ENSG00000124766.5")],
                    gmara_ensembl_genes = gmara_ensembl_genes,
                )
                @test consistencies[1].is_consistent
            end

            nested_test("unique") do
                ensembl_per_name = [
                    "ENSG00000182378.15",
                    "ENSG00000124766",
                    "",
                    "ENSG00000182378.15_PAR_Y",
                    "ERCC-00002",
                    "",
                    "ENSG00000124766",
                ]
                @test indices_of_duplicates(ensembl_per_name; empty_value = "") == Dict("ENSG00000124766" => [2, 7])
            end
        end
    end
end
