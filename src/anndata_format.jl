"""
Import and export `metacells` data from/to `h5ad` files. This allows moving data between the old Python/C++ based
`AnnData` world and the new Julia based `Daf` world.
"""
module AnnDataFormat

export import_h5ads!

using CSV
using Daf
using Daf.Generic
using DataFrames
using SparseArrays

CopyData = Dict{AbstractString, Maybe{Tuple{AbstractString, Maybe{StorageScalar}}}}

GENE_VECTORS_DATA = CopyData([
    "excluded_gene" => nothing,
    "full_gene_index" => nothing,
    "bursty_lonely_gene" => ("is_bursty_lonely", false),
    "properly_sampled_gene" => ("is_properly_sampled", false),
    "lateral_gene" => ("is_lateral", false),
    "noisy_gene" => ("is_noisy", false),
    "rare_gene" => ("is_rare", false),
    "marker_gene" => ("is_marker", false),
    "selected_gene" => ("is_selected", false),
    "ignored_gene" => ("is_ignored", false),
    "projected_noisy_gene" => ("is_projected_noisy", false),
    "atlas_gene" => ("is_atlas", false),
    "atlas_lateral_gene" => ("is_atlas_lateral", false),
    "atlas_noisy_gene" => ("is_atlas_noisy", false),
    "atlas_marker_gene" => ("is_atlas_marker", false),
    "rare_gene_module" => ("rare_module", 0),
    "lateral_gene_module" => ("lateral_module", 0),
])

CELL_VECTORS_DATA = CopyData([
    "full_gene_index" => nothing,
    "excluded_cell" => nothing,
    "properly_sampled_cell" => ("is_properly_sampled", false),
    "rare_cell" => ("is_rare", false),
    "dissolve" => ("is_dissolved", false),
    "metacell" => nothing,
    "metacell_name" => ("metacell", nothing),
    "most_similar" => nothing,
    "most_similar_name" => ("metacell.most_similar", nothing),
    "cells_rare_gene_module" => ("rare_gene_module", 0),
])

METACELL_VECTORS_DATA = CopyData([
    "rare_metacell" => ("is_rare", nothing),
    "metacells_level" => ("level", nothing),
    "similar" => ("is_similar", nothing),
])

NO_MATRICES_DATA = CopyData()

METACELLS_MATRICES_DATA = CopyData([
    "fitted" => ("is_fitted", nothing),
    "misfit" => ("is_misfit", nothing),
    "essential" => ("is_essential", nothing),
])

"""
    function import_h5ads!(;
        destination::DafWriter,
        raw_cells_h5ad::Maybe{AbstractString} = nothing,
        clean_cells_h5ad::AbstractString = nothing,
        metacells_h5ad::AbstractString,
        type_property::Maybe{AbstractString} = nothing,
        rename_type::Maybe{AbstractString} = "type",
        type_colors_csv::Maybe{AbstractString} = nothing,
        type_properties::Maybe{AbstractStringSet} = nothing,
        properties_defaults::Maybe{Dict} = nothing,
    )::Nothing

Import an `AnnData` based metacells dataset into a `Daf` `destination` data set. Ideally, the input must include
`clean_cells_h5ad` and the `metacells_h5ad` computed for them, and optionally also the `raw_cells_h5ad` including the
excluded cells and genes.

If type annotations were assigned to the metacells, then the name of the `type_property` should be specified. This can
be further enhanced by specifying a `type_colors_csv` file mapping type names to colors. This should be a comma or tab
separated file containing at least two columns, one named "color" and one with the same name as the `type_property`. For
consistency, by default the `type_property` is renamed to the value of `rename_type` (by default, "type"). You can
disable this by setting `rename_type` to `nothing`. We also call `reconstruct_axis!` to build the type axis; you can
therefore specify an `empty_type` name, which will be converted to the empty string, to match the `Daf` convention of
"no value" for string data, and specify an explicit set of `type_properties` (by default, any per-metacell property that
has the same value for all metacells of each type will be converted to a type property) and `properties_defaults`.

This will mostly just read all the specified `h5ad` files and copy the data into the `destination`, with the following
changes to match the ``Daf`` capabilities and conventions:

  - The `X` matrix of the cells is renamed to `UMIs`, and the `X` matrix of the metacells is renamed to `fraction`.
  - The `excluded_gene` and `excluded_cell` masks are not copied. Instead, if `raw_cells_h5ad` is specified, an
    `is_excluded` mask is created for both cells and genes, marking these that exist only in the `raw_cells_h5ad` and
    not in `clean_cells_h5ad` and `metacells_h5ad`.
  - The `full_gene_index` is not copied.
  - The `properly_sampled_gene` mask is renamed to the per-gene `is_properly_sampled` mask.
  - The `bursty_lonely_gene` mask is renamed to the per-gene `is_bursty_lonely` mask.
  - The `lateral_gene` mask is renamed to the per-gene `is_lateral` mask.
  - The `noisy_gene` mask is renamed to the per-gene `is_noisy` mask.
  - The `rare_gene` mask is renamed to the per-gene `is_rare` mask.
  - The `rare_gene_module` has 1 added to it (that is, "no module" is 0 in `Daf`) and is renamed to `rare_module`.
  - The `lateral_gene_module` has 1 added to it (that is, "no module" is 0 in `Daf`) and is renamed to `lateral_module`.
  - The `marker_gene` mask is renamed to the per-gene `is_marker` mask.
  - The `selected_gene` mask is renamed to the per-gene `is_selected` mask.
  - The `ignored_gene` mask is renamed to the per-gene `is_ignored` mask.
  - The `ignored_gene_of_<type>` masks are converted to an `is_ignored` mask per-gene-per-type.
  - The `projected_noisy_gene` mask is renamed to the per-gene `is_projected_noisy` mask.
  - The `atlas_gene`, `atlas_lateral_gene`, `atlas_noisy_gene`, `atlas_marker_gene` masks are renamed to
    the `is_atlas`, `is_atlas_lateral`, `is_atlas_noisy` and `is_atlas_marker` per-gene masks.
  - The `essential_gene_of_<type>` masks are converted to an `is_essential` mask per-gene-per-type.
  - The `atlas_essential_gene_of_<type>` masks are converted to an `is_atlas_essential` mask per-gene-per-type.
  - The `fitted_gene_of_<type>` masks are converted to an `is_fitted` mask per-gene-per-type.
  - The `fitted` mask per-gene-per-metacell is renamed to `is_fitted`.
  - The `misfit` mask per-gene-per-metacell is renamed to `is_misfit`.
  - The `essential` mask per-gene-per-metacell is renamed to `is_essential`.
  - The `full_cell_index` is not copied.
  - The `properly_sampled_cell` mask is renamed to the per-cell `is_properly_sampled` mask.
  - The `rare_cell` mask is renamed to the per-cell `is_rare` mask.
  - The `cells_rare_gene_module` has 1 added to it (that is, "no module" is 0 in `Daf`) and is renamed to `rare_gene_module`.
  - The per-cell `dissolve` mask is renamed to `is_dissolved`.
  - The per-cell `metacell` integer annotation is not copied, and the `metacell_name` string annotation is renamed to
    `metacell`.
  - The per-cell `most_similar` integer annotation is not copied, and the `most_similar_name` string annotation is
    renamed to `metacell.most_similar`.
  - The `rare_metacell` mask is renamed to the per-metacell `is_rare` mask.
  - The per-metacell `metacells_level` is renamed to `level`.
  - The per-metacell `similar` mask is renamed to `is_similar`.

!!! note

    There is much duplication of data between the three `h5ad` files (in particular, per-gene data). Data in
    `raw_cells_h5ad` will override data in `clean_cells_h5ad`, which will override data in `metacells_h5ad`.

!!! note

    It is common to call `reconstruct_axis!` on the result (e.g., if the cells were collected from a set of batches).
"""
function import_h5ads!(;
    destination::DafWriter,
    raw_cells_h5ad::Maybe{AbstractString} = nothing,
    clean_cells_h5ad::AbstractString = nothing,
    metacells_h5ad::AbstractString,
    type_property::Maybe{AbstractString} = nothing,
    rename_type::Maybe{AbstractString} = "type",
    empty_type::Maybe{AbstractString} = nothing,
    type_colors_csv::Maybe{AbstractString} = nothing,
    type_properties::Maybe{AbstractStringSet} = nothing,
    properties_defaults::Maybe{Dict} = nothing,
)::Nothing
    if raw_cells_h5ad != nothing
        copy_raw_cells(destination, raw_cells_h5ad)
    end

    copy_clean_cells(destination, clean_cells_h5ad; copy_axes = raw_cells_h5ad == nothing)

    return copy_metacells(
        destination,
        metacells_h5ad,
        type_property,
        rename_type,
        empty_type,
        type_colors_csv,
        type_properties,
        properties_defaults,
    )
end

function copy_raw_cells(destination::DafWriter, raw_cells_h5ad::AbstractString)::Nothing
    raw_cells_daf = anndata_as_daf(raw_cells_h5ad; obs_is = "cell", var_is = "gene", X_is = "UMIs")

    copy_axis!(; destination = destination, source = raw_cells_daf, axis = "cell")
    copy_axis!(; destination = destination, source = raw_cells_daf, axis = "gene")

    copy_scalars_data(destination, raw_cells_daf)

    copy_vectors(destination, raw_cells_daf, "gene", GENE_VECTORS_DATA)
    copy_vectors(destination, raw_cells_daf, "cell", CELL_VECTORS_DATA)

    copy_matrices(destination, raw_cells_daf, "cell", "gene", NO_MATRICES_DATA)
    copy_matrices(destination, raw_cells_daf, "gene", "cell", NO_MATRICES_DATA)

    return nothing
end

function copy_clean_cells(destination::DafWriter, clean_cells_h5ad::AbstractString; copy_axes::Bool)::Nothing
    clean_cells_daf = anndata_as_daf(clean_cells_h5ad; obs_is = "cell", var_is = "gene", X_is = "UMIs")

    if copy_axes
        copy_axis!(; destination = destination, source = clean_cells_daf, axis = "cell")
        copy_axis!(; destination = destination, source = clean_cells_daf, axis = "gene")
    end

    copy_scalars_data(destination, clean_cells_daf)

    for axis in ("cell", "gene")
        copy_vector!(;  # NOJET
            destination = destination,
            source = clean_cells_daf,
            axis = axis,
            name = "is_excluded",
            default = false,
            empty = true,
        )
    end

    copy_vectors(destination, clean_cells_daf, "gene", GENE_VECTORS_DATA)
    copy_vectors(destination, clean_cells_daf, "cell", CELL_VECTORS_DATA)

    copy_matrices(destination, clean_cells_daf, "cell", "gene", NO_MATRICES_DATA)
    copy_matrices(destination, clean_cells_daf, "gene", "cell", NO_MATRICES_DATA)

    return nothing
end

function copy_metacells(
    destination::DafWriter,
    metacells_h5ad::AbstractString,
    type_property::Maybe{AbstractString},
    rename_type::Maybe{AbstractString},
    empty_type::Maybe{AbstractString},
    type_colors_csv::Maybe{AbstractString},
    type_properties::Maybe{AbstractStringSet},
    properties_defaults::Maybe{Dict},
)::Nothing
    metacells_daf = anndata_as_daf(metacells_h5ad; obs_is = "metacell", var_is = "gene", X_is = "fraction")

    copy_axis!(; destination = destination, source = metacells_daf, axis = "metacell")

    copy_scalars_data(destination, metacells_daf)

    copy_vectors(destination, metacells_daf, "gene", GENE_VECTORS_DATA)
    copy_vectors(destination, metacells_daf, "metacell", METACELL_VECTORS_DATA, type_property, rename_type)

    copy_matrices(destination, metacells_daf, "metacell", "gene", METACELLS_MATRICES_DATA)
    copy_matrices(destination, metacells_daf, "gene", "metacell", METACELLS_MATRICES_DATA)

    if type_property != nothing
        copy_metacell_types(
            destination,
            metacells_daf,
            type_property,
            rename_type,
            empty_type,
            type_colors_csv,
            type_properties,
            properties_defaults,
        )
    end

    return nothing
end

function copy_metacell_types(
    destination::DafWriter,
    metacells_daf::DafReader,
    type_property::AbstractString,
    rename_type::Maybe{AbstractString},
    empty_type::Maybe{AbstractString},
    type_colors_csv::Maybe{AbstractString},
    type_properties::Maybe{AbstractStringSet},
    properties_defaults::Maybe{Dict},
)::Nothing
    if rename_type == nothing
        rename_type = type_property
    end

    if type_colors_csv != nothing
        data_frame = CSV.read(type_colors_csv, DataFrame)  # NOJET
        names = data_frame[:, type_property]
        colors = data_frame[:, "color"]
        add_axis!(destination, rename_type, names)
        set_vector!(destination, rename_type, "color", colors)
    end

    reconstruct_axis!(
        destination;
        existing_axis = "metacell",
        implicit_axis = rename_type,
        empty_implicit = empty_type,
        implicit_properties = type_properties,
        properties_defaults = properties_defaults,
    )

    type_names = get_axis(destination, rename_type)
    for prefix in ("ignored", "essential", "atlas_essential", "fitted")
        copy_mask_matrix(destination, metacells_daf, rename_type, type_names, prefix)
    end

    return nothing
end

function copy_scalars_data(destination::DafWriter, source::DafReader)::Nothing
    for scalar_name in scalar_names(source)
        if !has_scalar(destination, scalar_name)
            copy_scalar!(; destination = destination, source = source, name = scalar_name)
        end
    end
    return nothing
end

function copy_vectors(
    destination::DafWriter,
    source::DafWriter,
    axis::AbstractString,
    copy_data::CopyData,
    type_property::Maybe{AbstractString} = nothing,
    rename_type::Maybe{AbstractString} = nothing,
)::Nothing
    for vector_name in vector_names(source, axis)
        if !contains(vector_name, "_gene_of_")
            if vector_name == type_property && rename_type != nothing
                data = (rename_type, nothing)
            else
                data = get(copy_data, vector_name, (vector_name, nothing))
            end
            if data != nothing
                rename, empty = data
                if !has_vector(destination, axis, vector_name)
                    if empty == 0 && !(empty isa Bool)
                        vector = get_vector(source, axis, vector_name)
                        vector .+= 1
                        set_vector!(source, axis, vector_name, SparseVector(vector); overwrite = true)
                    end
                    copy_vector!(;  # NOJET
                        destination = destination,
                        source = source,
                        axis = axis,
                        name = vector_name,
                        rename = rename,
                        empty = empty,
                    )
                end
            end
        end
    end
    return nothing
end

function copy_matrices(
    destination::DafWriter,
    source::DafReader,
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    copy_data::CopyData,
)::Nothing
    for matrix_name in matrix_names(source, rows_axis, columns_axis; relayout = false)
        if !has_matrix(destination, rows_axis, columns_axis, matrix_name; relayout = true)
            data = get(copy_data, matrix_name, nothing)
            if data != nothing
                rename, empty = data
                copy_matrix!(;  # NOJET
                    destination = destination,
                    source = source,
                    rows_axis = rows_axis,
                    columns_axis = columns_axis,
                    name = matrix_name,
                    rename = rename,
                    empty = empty,
                    relayout = true,
                )
            end
        end
    end
    return nothing
end

function copy_mask_matrix(
    destination::DafWriter,
    source::DafReader,
    rename_type::AbstractString,
    type_names::AbstractStringVector,
    prefix::AbstractString,
)::Nothing
    mask_vectors = Vector{SparseVector{Bool}}()
    any_exist = false
    for type_name in type_names
        mask_name = "$(prefix)_gene_of_type_$(type_name)"
        if has_vector(source, "gene", mask_name)
            any_exist = true
        end
        mask_vector = get_vector(source, "gene", mask_name; default = false)
        @assert eltype(mask_vector) == Bool
        if !(mask_name isa SparseVector{Bool})
            mask_vector = SparseVector{Bool}(mask_vector)
        end
        push!(mask_vectors, mask_vector)
    end
    mask_matrix::SparseMatrixCSC{Bool} = hcat(mask_vectors...)  # NOJET
    mask_name = "is_$(prefix)"
    set_matrix!(source, "type", "gene", mask_name, mask_matrix)
    return copy_matrix!(;
        destination = destination,
        source = source,
        row_axis = rename_type,
        columns_axis = "gene",
        name = mask_name,
        empty = false,
        relayout = true,
    )
end

end  # module
