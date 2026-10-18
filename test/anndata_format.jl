# Write a small cells `h5ad` file, whose properties have the `AnnData` names the import renames. If the
# `metacell_name_per_cell` is `nothing`, the file has no `metacell_name`. Similarly for `most_similar_name_per_cell`.
function write_cells_h5ad(
    path::AbstractString;
    metacell_name_per_cell::Maybe{AbstractVector{<:AbstractString}} = ["M1", "M1", "M2", ""],
    most_similar_name_per_cell::Maybe{AbstractVector{<:AbstractString}} = ["", "", "", "M2"],
)::Nothing
    daf = MemoryDaf(; name = "cells!")
    add_axis!(daf, "obs", ["C1", "C2", "C3", "C4"])
    add_axis!(daf, "var", ["G1", "G2", "G3"])
    set_matrix!(daf, "obs", "var", "X", Int32[1 0 2; 0 3 0; 4 0 0; 0 0 5]; relayout = true)

    set_scalar!(daf, "__name__", "cells")
    set_scalar!(daf, "batch_size", 4)

    set_vector!(daf, "var", "correction_factor", [1.0, 2.0, 0.5])
    set_vector!(daf, "var", "full_gene_index", [10, 20, 30])
    set_vector!(daf, "var", "fitted", [true, false, true])
    set_vector!(daf, "var", "significant_inner_folds_count", [1, 0, 2])
    set_vector!(daf, "var", "lateral_gene", [false, true, false])
    set_vector!(daf, "var", "noisy_genes", [true, false, false])
    set_vector!(daf, "var", "gene_ids", ["ENSG00000124766.5", "ENSG00000170889", ""])
    set_vector!(daf, "var", "gene_foo_module", [-1, 0, 2])
    set_vector!(daf, "var", "genes_bar_module", [0, -1, 1])
    set_vector!(daf, "var", "baz_gene_module", [1, 0, -1])
    set_vector!(daf, "var", "qux_genes_module", [2, 2, 0])
    set_vector!(daf, "var", "score", [0.5, 1.5, 2.5])

    set_vector!(daf, "obs", "dissolve", [false, true, false, false])
    set_vector!(daf, "obs", "full_cell_index", [1, 2, 3, 4])
    set_vector!(daf, "obs", "metacell", [0, 0, 1, -1])
    if metacell_name_per_cell !== nothing
        set_vector!(daf, "obs", "metacell_name", metacell_name_per_cell)
    end
    if most_similar_name_per_cell !== nothing
        set_vector!(daf, "obs", "most_similar_name", most_similar_name_per_cell)
    end
    set_vector!(daf, "obs", "total_umis", [3, 3, 4, 5])
    set_vector!(daf, "obs", "doublet_cell", [false, false, true, false])
    set_vector!(daf, "obs", "batch", ["B1", "B1", "B2", "B2"])

    daf_as_anndata(daf; h5ad = path)
    return nothing
end

# Write a small metacells `h5ad` file, for the cells of the cells `h5ad` file.
function write_metacells_h5ad(path::AbstractString)::Nothing
    daf = MemoryDaf(; name = "metacells!")
    add_axis!(daf, "obs", ["M1", "M2"])
    add_axis!(daf, "var", ["G1", "G2", "G3"])
    set_matrix!(daf, "obs", "var", "X", [0.1 0.2 0.7; 0.5 0.25 0.25]; relayout = true)
    set_matrix!(daf, "obs", "var", "total_umis", [5 3 2; 4 0 0]; relayout = true)
    set_matrix!(daf, "obs", "var", "zeros", [0 2 1; 3 0 0]; relayout = true)
    set_matrix!(daf, "obs", "var", "custom", [1.5 0.0 2.5; 0.0 3.5 0.0]; relayout = true)
    set_matrix!(daf, "obs", "obs", "obs_outgoing_weights", [0.0 1.0; 0.5 0.0])

    set_scalar!(daf, "__name__", "metacells")
    set_scalar!(daf, "target_metacell_size", 96)

    set_vector!(daf, "var", "fitted", [true, true, false])

    set_vector!(daf, "obs", "metacells_level", [0, 1])
    set_vector!(daf, "obs", "similar", [true, false])
    set_vector!(daf, "obs", "grouped", [2, 1])
    set_vector!(daf, "obs", "type", ["T1", "T2"])
    set_vector!(daf, "obs", "u", [1.0, 2.0])
    set_vector!(daf, "obs", "v", [3.0, 4.0])

    daf_as_anndata(daf; h5ad = path)
    return nothing
end

nested_test("anndata_format") do
    mktempdir() do directory
        cells_h5ad = "$(directory)/cells.h5ad"
        write_cells_h5ad(cells_h5ad)

        nested_test("import_cells_h5ad") do
            daf = MemoryDaf(; name = "cells!")

            nested_test("default") do
                @test import_cells_h5ad!(daf; cells_h5ad) == ["M1", "M1", "M2", ""]

                @test axis_vector(daf, "cell") == ["C1", "C2", "C3", "C4"]
                @test axis_vector(daf, "gene") == ["G1", "G2", "G3"]

                @test scalars_set(daf) == Set(["batch_size"])
                @test get_scalar(daf, "batch_size") == 4

                umis = get_matrix(daf, "cell", "gene", "UMIs").array
                @test eltype(umis) == UInt32
                @test umis == [1 0 2; 0 3 0; 4 0 0; 0 0 5]

                @test vectors_set(daf, "gene") == Set([
                    "correction_factor",
                    "is_fitted",
                    "significant_inner_folds_count",
                    "is_lateral",
                    "is_noisy",
                    "gene_ids",
                    "foo_module",
                    "bar_module",
                    "baz_module",
                    "qux_module",
                    "score",
                ])
                @test eltype(get_vector(daf, "gene", "correction_factor").array) == Float32
                @test get_vector(daf, "gene", "correction_factor").array == [1.0, 2.0, 0.5]
                @test get_vector(daf, "gene", "is_fitted").array == [true, false, true]
                @test eltype(get_vector(daf, "gene", "significant_inner_folds_count").array) == UInt32
                @test get_vector(daf, "gene", "significant_inner_folds_count").array == [1, 0, 2]
                @test get_vector(daf, "gene", "is_lateral").array == [false, true, false]
                @test get_vector(daf, "gene", "is_noisy").array == [true, false, false]
                @test get_vector(daf, "gene", "gene_ids").array == ["ENSG00000124766.5", "ENSG00000170889", ""]
                @test eltype(get_vector(daf, "gene", "foo_module").array) == UInt32
                @test get_vector(daf, "gene", "foo_module").array == [0, 1, 3]
                @test get_vector(daf, "gene", "bar_module").array == [1, 0, 2]
                @test get_vector(daf, "gene", "baz_module").array == [2, 1, 0]
                @test get_vector(daf, "gene", "qux_module").array == [3, 3, 1]
                @test get_vector(daf, "gene", "score").array == [0.5, 1.5, 2.5]

                @test vectors_set(daf, "cell") == Set(["is_dissolved", "total_UMIs", "is_doublet", "batch"])
                @test get_vector(daf, "cell", "is_dissolved").array == [false, true, false, false]
                @test eltype(get_vector(daf, "cell", "total_UMIs").array) == UInt32
                @test get_vector(daf, "cell", "total_UMIs").array == [3, 3, 4, 5]
                @test get_vector(daf, "cell", "is_doublet").array == [false, false, true, false]
                @test get_vector(daf, "cell", "batch").array == ["B1", "B1", "B2", "B2"]
            end

            nested_test("copy_data") do
                copy_data = CopyAnnData(
                    "batch_size" => ("n_batches", nothing),
                    ("gene", "score") => nothing,
                    ("gene", "fitted") => ("=", nothing),
                    ("cell", "batch") => ("plate", nothing),
                    ("gene", "cell", "X") => ("=", nothing),
                )
                import_cells_h5ad!(daf; cells_h5ad, copy_data)
                @test get_matrix(daf, "cell", "gene", "X").array == [1 0 2; 0 3 0; 4 0 0; 0 0 5]
                @test !has_matrix(daf, "cell", "gene", "UMIs")
                @test scalars_set(daf) == Set(["n_batches"])
                @test !has_vector(daf, "gene", "score")
                @test get_vector(daf, "gene", "fitted").array == [true, false, true]
                @test !has_vector(daf, "gene", "is_fitted")
                @test get_vector(daf, "cell", "plate").array == ["B1", "B1", "B2", "B2"]
                @test !has_vector(daf, "cell", "batch")
            end

            nested_test("twice") do
                import_cells_h5ad!(daf; cells_h5ad)
                set_vector!(daf, "gene", "score", [9.0, 9.0, 9.0]; overwrite = true)
                set_matrix!(daf, "cell", "gene", "UMIs", UInt32[9 9 9; 9 9 9; 9 9 9; 9 9 9]; overwrite = true)
                import_cells_h5ad!(daf; cells_h5ad)
                @test get_vector(daf, "gene", "score").array == [9.0, 9.0, 9.0]
                @test get_matrix(daf, "cell", "gene", "UMIs").array == [9 9 9; 9 9 9; 9 9 9; 9 9 9]
            end

            nested_test("ensembl") do
                import_cells_h5ad!(daf; cells_h5ad, ensembl_property = "gene_ids")
                @test get_vector(daf, "gene", "ensembl").array == ["ENSG00000124766.5", "ENSG00000170889", ""]
                @test !has_vector(daf, "gene", "gene_ids")
            end

            nested_test("!ensembl") do
                expected_message = "both the ensembl_property and the copy_data specify the per-gene: gene_ids"
                @test_throws expected_message import_cells_h5ad!(
                    daf;
                    cells_h5ad,
                    copy_data = CopyAnnData(("gene", "gene_ids") => nothing),
                    ensembl_property = "gene_ids",
                )
            end
        end

        metacells_h5ad = "$(directory)/metacells.h5ad"
        write_metacells_h5ad(metacells_h5ad)

        nested_test("import_metacells_h5ad") do
            daf = MemoryDaf(; name = "metacells!")

            nested_test("default") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                import_metacells_h5ad!(daf; cells_h5ad, metacells_h5ad)

                @test axis_vector(daf, "metacell") == ["M1", "M2"]
                @test axis_vector(daf, "gene") == ["G1", "G2", "G3"]

                @test scalars_set(daf) == Set(["target_metacell_size"])
                @test get_scalar(daf, "target_metacell_size") == 96

                @test vectors_set(daf, "cell") == Set(["metacell", "most_similar.metacell", "is_base_outlier"])
                @test get_vector(daf, "cell", "metacell").array == ["M1", "M1", "M2", ""]
                @test get_vector(daf, "cell", "most_similar.metacell").array == ["", "", "", "M2"]
                @test get_vector(daf, "cell", "is_base_outlier").array == [false, false, false, true]

                @test vectors_set(daf, "gene") == Set(["is_fitted"])
                @test get_vector(daf, "gene", "is_fitted").array == [true, true, false]

                @test vectors_set(daf, "metacell") == Set(["level", "is_similar", "n_cell", "type", "umap_u", "umap_v"])
                @test eltype(get_vector(daf, "metacell", "level").array) == UInt32
                @test get_vector(daf, "metacell", "level").array == [0, 1]
                @test get_vector(daf, "metacell", "is_similar").array == [true, false]
                @test eltype(get_vector(daf, "metacell", "n_cell").array) == UInt32
                @test get_vector(daf, "metacell", "n_cell").array == [2, 1]
                @test get_vector(daf, "metacell", "type").array == ["T1", "T2"]
                @test get_vector(daf, "metacell", "umap_u").array == [1.0, 2.0]
                @test get_vector(daf, "metacell", "umap_v").array == [3.0, 4.0]

                @test matrices_set(daf, "metacell", "gene") == Set(["geomean_fraction", "UMIs", "zeros", "custom"])
                geomean_fraction = get_matrix(daf, "metacell", "gene", "geomean_fraction").array
                @test eltype(geomean_fraction) == Float32
                @test geomean_fraction == Float32[0.1 0.2 0.7; 0.5 0.25 0.25]
                umis = get_matrix(daf, "metacell", "gene", "UMIs").array
                @test eltype(umis) == UInt32
                @test umis == [5 3 2; 4 0 0]
                zeros_count = get_matrix(daf, "metacell", "gene", "zeros").array
                @test eltype(zeros_count) == UInt32
                @test zeros_count == [0 2 1; 3 0 0]
                @test get_matrix(daf, "metacell", "gene", "custom").array == [1.5 0.0 2.5; 0.0 3.5 0.0]

                @test matrices_set(daf, "metacell", "metacell") == Set(["outgoing_weights"])
                outgoing_weights = get_matrix(daf, "metacell", "metacell", "outgoing_weights").array
                @test eltype(outgoing_weights) == Float32
                @test outgoing_weights == [0.0 1.0; 0.5 0.0]
            end

            nested_test("more_cells") do
                add_axis!(daf, "cell", ["C0", "C1", "C2", "C3", "C4"])
                import_metacells_h5ad!(daf; cells_h5ad, metacells_h5ad)
                @test get_vector(daf, "cell", "metacell").array == ["", "M1", "M1", "M2", ""]
                @test get_vector(daf, "cell", "most_similar.metacell").array == ["", "", "", "", "M2"]
                @test get_vector(daf, "cell", "is_base_outlier").array == [true, false, false, false, true]
            end

            nested_test("!most_similar_name") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                plain_cells_h5ad = "$(directory)/plain_cells.h5ad"
                write_cells_h5ad(plain_cells_h5ad; most_similar_name_per_cell = nothing)
                import_metacells_h5ad!(daf; cells_h5ad = plain_cells_h5ad, metacells_h5ad)
                @test vectors_set(daf, "cell") == Set(["metacell", "is_base_outlier"])
                @test get_vector(daf, "cell", "metacell").array == ["M1", "M1", "M2", ""]
            end

            nested_test("copy_data") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                import_metacells_h5ad!(
                    daf;
                    cells_h5ad,
                    metacells_h5ad,
                    copy_data = CopyAnnData(("gene", "metacell", "zeros") => nothing),
                )
                @test matrices_set(daf, "metacell", "gene") == Set(["geomean_fraction", "UMIs", "custom"])
            end

            nested_test("twice") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                import_metacells_h5ad!(daf; cells_h5ad, metacells_h5ad)
                set_vector!(daf, "cell", "metacell", ["M2", "M2", "M1", ""]; overwrite = true)
                import_metacells_h5ad!(daf; cells_h5ad, metacells_h5ad)
                @test get_vector(daf, "cell", "metacell").array == ["M2", "M2", "M1", ""]
            end

            nested_test("!metacell_name") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                bad_cells_h5ad = "$(directory)/bad_cells.h5ad"
                write_cells_h5ad(bad_cells_h5ad; metacell_name_per_cell = nothing)
                @test_throws "missing per-obs annotation: metacell_name" import_metacells_h5ad!(
                    daf;
                    cells_h5ad = bad_cells_h5ad,
                    metacells_h5ad,
                )
            end

            nested_test("!metacell") do
                add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
                bad_cells_h5ad = "$(directory)/bad_cells.h5ad"
                write_cells_h5ad(bad_cells_h5ad; metacell_name_per_cell = ["M1", "M1", "M3", ""])
                @test_throws "invalid metacell_name: M3\nfor the cell: C3" import_metacells_h5ad!(
                    daf;
                    cells_h5ad = bad_cells_h5ad,
                    metacells_h5ad,
                )
            end
        end

        nested_test("import_gene_masks_per_type") do
            daf = MemoryDaf(; name = "metacells!")
            add_axis!(daf, "gene", ["G1", "G2", "G3"])
            add_axis!(daf, "type", ["T1", "T2", "T3"])
            set_vector!(daf, "gene", "marker_gene_of_T1", [true, false, false])
            set_vector!(daf, "gene", "marker_gene_of_T2", [false, true, true])
            import_gene_masks_per_type!(daf)
            @test matrices_set(daf, "gene", "type") == Set(["is_marker"])
            @test get_matrix(daf, "gene", "type", "is_marker").array ==
                  [true false false; false true false; false true false]
        end

        nested_test("import_type_colors_csv") do
            daf = MemoryDaf(; name = "cells!")
            add_axis!(daf, "cell", ["C1", "C2", "C3", "C4"])
            set_vector!(daf, "cell", "type", ["T1", "", "T2", "T1"])
            type_colors_csv = "$(directory)/type_colors.csv"

            nested_test("default") do
                write(type_colors_csv, "type,color\nT1,red\nT2,blue\nT3,green\n")
                import_type_colors_csv!(daf; type_colors_csv)
                @test axis_vector(daf, "type") == ["T1", "T2", "T3"]
                @test get_vector(daf, "type", "color").array == ["red", "blue", "green"]
            end

            nested_test("!color") do
                write(type_colors_csv, "type,colour\nT1,red\nT2,blue\n")
                @test_throws "expected two columns, the types and then \"color\": type, colour" import_type_colors_csv!(
                    daf;
                    type_colors_csv,
                )
                @test !has_axis(daf, "type")
            end

            nested_test("!unique") do
                write(type_colors_csv, "type,color\nT2,red\nT1,blue\nT2,green\nT1,pink\n")
                @test_throws "repeated type(s): T1, T2\n" import_type_colors_csv!(daf; type_colors_csv)
                @test !has_axis(daf, "type")
            end

            nested_test("!listed") do
                write(type_colors_csv, "type,color\nT1,red\nT3,green\n")
                @test_throws "missing from the type colors csv: $(type_colors_csv)\nare: T2" import_type_colors_csv!(
                    daf;
                    type_colors_csv,
                )
                @test !has_axis(daf, "type")
            end
        end
    end
end
