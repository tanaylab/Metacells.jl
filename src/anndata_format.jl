"""
Import and export cells and metacells data from/to `h5ad` files. This allows moving data between the old Python/C++
based [AnnData](https://github.com/scverse/anndata) world and the brave new Julia based
[Daf](https://github.com/tanaylab/DataAxesFormats.jl) world.

The expected flow is as follows:

  - Create a `Daf` repository for the raw input cells and import the raw cells `h5ad` into it using
    [`import_cells_h5ad!`](@ref).

  - Alternatively, import just the clean cells `h5ad` into it, but that would be leaving out some of the data so is not
    recommended.
  - Create another `Daf` repository for the metacells, chain it to the cells repository, and import the metacells into
    it using [`import_metacells_h5ad!`](@ref). Give this the cells-with-metacells data `h5ad` - only the assignment of
    cells to metacells will be imported from it.
  - If you have any per-cell or per-gene computed data in the cells-with-metacells data `h5ad` (unlikely, as computed
    data typically goes into the metacells `h5ad`), import it into the chained (metacells) repository using
    [`import_cells_h5ad!`](@ref).
  - Create a type axis in the chained metacells `Daf` repository using [`reconstruct_type_axis!`](@ref).
"""
module AnnDataFormat

export import_cells_h5ad!
export import_metacells_h5ad!
export reconstruct_type_axis!
export CopyAnnData

using CSV
using DataAxesFormats
using DataFrames
using Muon
using SparseArrays
using TanayLabUtilities

"""
Specify how to copy data from `AnnData` to `Daf`. The key is simply a vector or matrix name (ignoring axes), and the
value is either `nothing` to ignore the data, or a tuple with the name of the destination `Daf` property and an optional
value to use for missing entries (raw-only cells and/or genes).
"""
CopyAnnData = Dict{AbstractString, Maybe{Tuple{AbstractString, Maybe{StorageScalar}}}}

GENE_VECTORS_DATA = CopyAnnData([
    "correction_factor" => ("correction_factor", Float32(0.0)),
    "full_gene_index" => nothing,
    "fitted" => ("is_fitted", false),
    "significant_inner_folds_count" => ("significant_inner_folds_count", UInt32(0)),
    "*" => ("=", nothing),
])

CELL_VECTORS_DATA = CopyAnnData([
    "dissolve" => ("is_dissolved", false),
    "full_cell_index" => nothing,
    "metacell" => nothing,
    "metacell_name" => nothing,
    "metacell_level" => nothing,
    "most_similar" => nothing,
    "most_similar_name" => nothing,
    "*" => ("=", nothing),
])

CELL_MATRICES_DATA = CopyAnnData(["X" => ("UMIs", UInt32(0))])

METACELL_VECTORS_DATA = CopyAnnData([
    "metacells_level" => ("level", UInt32(0)),
    "similar" => ("is_similar", false),
    "grouped" => ("n_cell", UInt32(0)),
    "type" => ("type", nothing),
    "*" => ("=", nothing),
])

METACELL_MATRICES_DATA = CopyAnnData([
    "X" => ("geomean_fraction", Float32(0.0)),
    "corrected_fraction" => ("corrected_geomean_fraction", Float32(0.0)),
    "essential" => ("is_essential", false),
    "fitted" => ("is_fitted", false),
    "inner_fold" => ("inner_fold", Float32(0.0)),
    "inner_stdev_log" => ("inner_std_log", Float32(0.0)),
    "misfit" => ("is_misfit", false),
    "projected_fold" => ("projected_fold", Float32(0.0)),
    "projected_fraction" => ("projected_geomean_fraction", Float32(0.0)),
    "total_umis" => ("UMIs", UInt32(0)),
    "zeros" => ("zeros", UInt32(0)),
    "*" => ("=", nothing),
])

METACELL_SQUARE_DATA = CopyAnnData(["obs_outgoing_weights" => ("outgoing_weights", Float32(0)), "*" => ("=", nothing)])

"""
    function import_cells_h5ad!(
        daf::DafWriter;
        cells_h5ad::AbstractString,
        copy_data::Maybe{CopyAnnData} = $(DEFAULT.copy_data),
        bestify::Bool = $(DEFAULT.bestify),
        min_sparse_saving_fraction::AbstractFloat = $(DEFAULT.min_sparse_saving_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
        insist::Bool = $(DEFAULT.insist),
    )::Nothing

Import an `AnnData` based cells dataset into a destination `daf` data set. Ideally you'd copy the full (raw) cells
into an empty `Daf` repository. Then, you'd treat this repository as read-only, and copy the metacells data using
[`import_metacells_h5ad!`](@ref) into a separate `Daf` repository chained with the read-only cells repository. This
allows separate alternative metacells computations to share the read-only cells data.

You can copy an `h5ad` file containing just the clean cells on top of the raw cells data, to capture any data computed
during or after computing the metacells. You may need to specify the `copy_data` to specify defaults for values of
properties that exist only for clean cells and/or genes. Or, you can skip copying the raw data altogether, copying just
the clean data into the base cells repository, though this is less recommended as you are needlessly discarding data
that may prove to be useful later.

The `bestify`, `min_sparse_saving_fraction`, `overwrite`, and `insist` have their usual meaning from `Daf`'s
copying functions.

When copying, we apply the following general rules:

  - A `something_gene` per-gene property and/or `something_cell` per-cell property are renamed to `is_something`,
    (and given a default of `false`), because `Daf` (unlike `AnnData`) has no problem with properties with the same name
    for different axes.
  - Similarly `{gene,cell}[s]_something_module` and `something_{gene,cell}[s]_module` properties are renamed to
    `something_module`. We add 1 to the value and store the results in a `UInt32`; that is, in `Daf`, module indices are
    1-based, and 0 is "no module".
  - Any `something_umis` is renamed to `something_UMIs`, given a default of 0, and stored as a `UInt32`.

And we make the following special exceptions:

Scalars:

  - We do not copy the `__name__` scalar.
  - All other scalars are copied as-is.

Per-cell-per-gene:

  - The `X` matrix is renamed to `UMIs`, and stored as a `UInt32`.
  - No other per-cell-per-gene matrix is copied by default.

Per-gene:

  - The per-gene `correction_factor` is given the default value 0.
  - The `fitted` per-gene vector is renamed to `is_fitted` and given the default `false`.
  - The `significant_inner_folds_count` per-gene property is stored as a `UInt32` and given the default 0.
  - The `full_gene_index` property is not copied. Either you import the full (raw) data or you don't.
  - All other vectors are copied as-is.

Per-cell:

  - The `full_cell_index` property is not copied. Either you import the full (raw) data or you don't.
  - The `metacell`, `metacell_name`, `metacell_level`, `most_similar` and `most_similar_name` per-cell properties are
    not copied. To import these, use [`import_metacells_h5ad!`](@ref).
  - All other vectors are copied as-is.

!!! note

    It is common to manually call `reconstruct_axis!` on the result to create additional axes (e.g., if the cells were
    collected from a set of batches and some properties are actually per-batch).
"""
@documented @logged function import_cells_h5ad!(
    daf::DafWriter;
    cells_h5ad::AbstractString,
    copy_data::Maybe{CopyAnnData} = nothing,
    bestify::Bool = true,
    min_sparse_saving_fraction::AbstractFloat = function_default(copy_matrix!, :min_sparse_saving_fraction),
    overwrite::Bool = false,
    insist::Bool = false,
)::Nothing
    cells_daf = anndata_as_daf(cells_h5ad; name = "cells", obs_is = "cell", var_is = "gene", X_is = "X")  # NOJET

    copy_axis!(; destination = daf, source = cells_daf, axis = "cell", overwrite, insist)
    copy_axis!(; destination = daf, source = cells_daf, axis = "gene", overwrite, insist)

    import_scalars_data(daf, cells_daf; overwrite, insist)

    import_vectors_data(
        daf,
        cells_daf,
        "gene",
        GENE_VECTORS_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )
    import_vectors_data(
        daf,
        cells_daf,
        "cell",
        CELL_VECTORS_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )

    import_matrices_data(
        daf,
        cells_daf,
        "cell",
        "gene",
        CELL_MATRICES_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )
    import_matrices_data(
        daf,
        cells_daf,
        "gene",
        "cell",
        CELL_MATRICES_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )

    return nothing
end

"""
    function import_metacells_h5ad!(
        daf::DafWriter;
        cells_h5ad::AbstractString,
        metacells_h5ad::AbstractString,
        copy_data::Maybe{CopyAnnData} = $(DEFAULT.copy_data),
        bestify::Bool = $(DEFAULT.bestify),
        min_sparse_saving_fraction::AbstractFloat = $(DEFAULT.min_sparse_saving_fraction),
        overwrite::Bool = $(DEFAULT.overwrite),
        insist::Bool = $(DEFAULT.insist),
    )::Nothing

Import an `AnnData` based metacells dataset into a destination `daf` data set. It is expected that you have first
imported the per-cell data. Ideally you'd create a new empty repository for the metacells data and chained it on top of
the per-cell repository, which you'd keep read-only to allow sharing it when you (inevitably) compute different
metacells for it.

This behaves similarly to [`import_cells_h5ad!`](@ref), specifically the generic rules (except that we copy per-metacell
properties and not per-cell properties so the rules are adjusted accordingly), and we make the following special exceptions:

Per-metacell-per-gene:

  - The `X` matrix is renamed to `fraction` and always stored as `Float32`.
  - The `corrected_fraction` matrix is always stored as `Float32`.
  - The `essential` matrix is renamed to `is_essential`.
  - The `essential`, `fitted`, and `misfit` matrices are renamed to `is_essential`, `is_fitted` and `is_misfit`, respectively.
  - The `inner_fold`, `inner_stdev_log` (renamed to `inner_std_log`), `projected_fold`, `projected_fraction` matrices
    are always stored as `Float32`.
  - The `total_umis` matrix is renamed to `UMIs` and always stored as `UInt32`.
  - The `zeros` matrix is always stored as `UInt32`.
  - All other matrices are copied as-is.

Scalars and Per-gene:

Same as in [`import_cells_h5ad!`](@ref)

Per-cell:

  - The only properties we copy per cell are `metacell_name` (renamed to `metacell` with a default of the empty string),
    and similarly `most_similar_name` (renamed to `most_similar.metacell`, same default). That's the only reason we have
    a `cells_h5ad` parameter. You should therefore pass here the cells-with-metacells and not the clean cells `h5ad`.

Per-metacell:

  - The `metacells_level` property is renamed to `level`.
  - The `similar` property is renamed to `is_similar`.
  - The `type` property is copied. If "the" type property of the metacells is different, use `copy_data` to rename
    it to `type` to match the `Daf` naming convention.
  - All other vectors are copied as-is.

Per-metacell-per-metacell:

  - The `obs_outgoing_weights` matrix is renamed to `outgoing_weights` and always stored as `Float32`.
  - All other matrices are copied as-is.

!!! note

    It is common to manually call `reconstruct_type!` on the result to create a type axis.
"""
@documented @logged function import_metacells_h5ad!(
    daf::DafWriter;
    cells_h5ad::AbstractString,
    metacells_h5ad::AbstractString,
    copy_data::Maybe{CopyAnnData} = nothing,
    bestify::Bool = true,
    min_sparse_saving_fraction::AbstractFloat = function_default(copy_matrix!, :min_sparse_saving_fraction),
    overwrite::Bool = false,
    insist::Bool = false,
)::Nothing
    metacells_daf = anndata_as_daf(metacells_h5ad; name = "metacells", obs_is = "metacell", var_is = "gene", X_is = "X")

    copy_axis!(; destination = daf, source = metacells_daf, axis = "metacell", overwrite, insist)
    copy_axis!(; destination = daf, source = metacells_daf, axis = "gene", overwrite, insist)

    copy_metacells_of_cells(daf, cells_h5ad; overwrite, insist)

    import_scalars_data(daf, metacells_daf; overwrite, insist)

    import_vectors_data(
        daf,
        metacells_daf,
        "gene",
        GENE_VECTORS_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )
    import_vectors_data(
        daf,
        metacells_daf,
        "metacell",
        METACELL_VECTORS_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )

    import_matrices_data(
        daf,
        metacells_daf,
        "metacell",
        "gene",
        METACELL_MATRICES_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )
    import_matrices_data(
        daf,
        metacells_daf,
        "gene",
        "metacell",
        METACELL_MATRICES_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )

    import_matrices_data(
        daf,
        metacells_daf,
        "metacell",
        "metacell",
        METACELL_SQUARE_DATA;
        copy_data,
        bestify,
        min_sparse_saving_fraction,
        overwrite,
        insist,
    )

    return nothing
end

"""
    reconstruct_type!(
        daf::DafWriter,
        base_axis::AbstractString = $(DEFAULT.base_axis),
        type_property::AbstractString = $(DEFAULT.type_property),
        type_axis::AbstractString = $(DEFAULT.type_axis),
        empty_type::Maybe{AbstractString} = $(DEFAULT.empty_type),
        type_colors_csv::Maybe{AbstractString} = $(DEFAULT.type_colors_csv),
        implicit_properties::Maybe{AbstractSet{<:AbstractString}} = $(DEFAULT.implicit_properties),
        skipped_properties::Maybe{AbstractSet{<:AbstractString}} = $(DEFAULT.skipped_properties),
        properties_defaults::Maybe{Dict} = $(DEFAULT.properties_defaults),
        overwrite::Bool = $(DEFAULT.overwrite),
    )::Nothing

Create a type axis after importing data containing type annotations. By default this will look for a type per metacell,
but if you have type annotation per cell (which is **not** simply the type of the metacell they belong to), you can also
use this for the cells.

By default this assumes that you have imported "the" type annotation to a property called "type", and that you would
like the new type axis to be called "type" as well. If you want to import secondary types (or per-cell types), change
these via the `type_property` and `type_axis` parameters.

If the type is equal to `empty_type` it is replaced with the empty string to match the `Daf` conventions for "no value"
for string properties. Any set of per-gene properties named `something_gene_of_type` is converted to a per-gene-per-type
matrix called `is_something` (with a default of `false`).

Otherwise, this is mostly just a wrapper for `reconstruct_axis!`. It can be further enhanced by specifying a
`type_colors_csv` file mapping type names to colors. This should be a comma or tab separated file containing at least
two columns, one named "color" and one with the same name as the `type_property`. If this CSV file contains types that
aren't actually used in the data, you will have to specify a default value for any other per-type property in
`properties_defaults`.

!!! note

    Most metacells data has type annotations and colors associated with types, so it is highly recommended you invoke
    this to capture these into the `Daf` repository. This will enable all types (:-) of downstream processing, coloring
    graphs, etc.
"""
@documented @logged function reconstruct_type_axis!(
    daf::DafWriter;
    base_axis::AbstractString = "metacell",
    type_property::AbstractString = "type",
    type_axis::AbstractString = "type",
    empty_type::Maybe{AbstractString} = nothing,
    type_colors_csv::Maybe{AbstractString} = nothing,
    implicit_properties::Maybe{AbstractSet{<:AbstractString}} = nothing,
    skipped_properties::Maybe{AbstractSet{<:AbstractString}} = Set(["rare_gene_module", "is_rare"]),
    properties_defaults::Maybe{Dict} = nothing,
)::Nothing
    if type_colors_csv !== nothing
        data_frame = CSV.read(type_colors_csv, DataFrame)  # NOJET
        names = data_frame[:, type_property]
        colors = data_frame[:, "color"]
        add_axis!(daf, type_axis, names)
        @info "set $(type_axis) vector: color"
        set_vector!(daf, type_axis, "color", colors)
        delete_vector!(daf, "metacell", "color"; must_exist = false)
        delete_vector!(daf, "cell", "color"; must_exist = false)
    end

    reconstruct_axis!(
        daf;
        existing_axis = base_axis,
        implicit_axis = type_axis,
        empty_implicit = empty_type,
        implicit_properties,
        skipped_properties,
        properties_defaults,
    )

    type_names = axis_vector(daf, type_axis)

    prefixes = Set{AbstractString}()
    for property in vectors_set(daf, "gene")
        parts = split(property, "_gene_of_")
        if length(parts) == 2
            push!(prefixes, parts[1])
        end
    end

    for prefix in prefixes
        import_mask_matrix(daf, type_axis, type_names, prefix)
    end

    return nothing
end

function import_mask_matrix(
    daf::DafWriter,
    type_axis::AbstractString,
    type_names::AbstractVector{<:AbstractString},
    prefix::AbstractString,
)::Nothing
    mask_vectors = Vector{SparseVector{Bool}}()
    for type_name in type_names
        mask_name = "$(prefix)_gene_of_$(type_name)"
        mask_vector = get_vector(daf, "gene", mask_name; default = false)
        @assert eltype(mask_vector) == Bool
        if !(mask_name isa SparseVector{Bool})
            mask_vector = SparseVector{Bool}(mask_vector)
        end
        push!(mask_vectors, mask_vector)
    end

    mask_matrix::SparseMatrixCSC{Bool} = hcat(mask_vectors...)  # NOJET
    @info "reconstruct gene-$(type_axis) matrix: is_$(prefix)"
    set_matrix!(daf, "gene", type_axis, "is_$(prefix)", bestify(mask_matrix); relayout = true, overwrite = true)  # NOJET

    # TODOX for type_name in type_names
    # TODOX mask_name = "$(prefix)_gene_of_$(type_name)"
    # TODOX delete_vector!(daf, "gene", mask_name; must_exist = false)
    # TODOX end

    return nothing
end

function import_scalars_data(daf::DafWriter, source::DafReader; overwrite::Bool, insist::Bool)::Nothing
    for scalar_name in scalars_set(source)
        if scalar_name == "__name__"
            @info "skip scalar: $(scalar_name)"
        else
            @info "copy scalar: $(scalar_name)"
            copy_scalar!(; destination = daf, source, name = scalar_name, overwrite, insist)
        end
    end
    return nothing
end

function import_vectors_data(
    daf::DafWriter,
    source::DafReader,
    axis::AbstractString,
    base_copy_data::CopyAnnData;
    copy_data::Maybe{CopyAnnData},
    bestify::Bool,
    min_sparse_saving_fraction::AbstractFloat,
    overwrite::Bool,
    insist::Bool,
)::Nothing
    for vector_name in vectors_set(source, axis)
        data = nothing

        if copy_data !== nothing
            data = get(copy_data, vector_name, nothing)
        end

        if data === nothing && (copy_data === nothing || !haskey(copy_data, vector_name))
            data = get(base_copy_data, vector_name, nothing)

            if data === nothing
                if endswith(vector_name, "_$(axis)")
                    data = ("is_$(vector_name[1:end - length(axis) - 1])", false)
                elseif endswith(vector_name, "_$(axis)s")
                    data = ("is_$(vector_name[1:end - length(axis) - 2])", false)
                elseif endswith(vector_name, "_umis")
                    data = ("$(vector_name[1:end - 5])_UMIs", UInt32(0))
                elseif endswith(vector_name, "_module")
                    if startswith(vector_name, "$(axis)_")
                        data = (vector_name[(length(axis) + 2):end], UInt32(0))
                    elseif startswith(vector_name, "$(axis)s_")
                        data = (vector_name[(length(axis) + 3):end], UInt32(0))
                    elseif endswith(vector_name, "_$(axis)_module")
                        data = ("$(vector_name[1:end - length(axis) - 8])_module", UInt32(0))
                    elseif endswith(vector_name, "_$(axis)s_module")
                        data = ("$(vector_name[1:end - length(axis) - 9])_module", UInt32(0))
                    end
                end

                if data === nothing && haskey(base_copy_data, "*") && !haskey(base_copy_data, vector_name)
                    data = (vector_name, nothing)
                end
            end
        end

        if data === nothing
            @info "skip $(axis) vector: $(vector_name)"
            continue
        end

        rename, empty = data

        if !overwrite && !insist && has_vector(daf, axis, rename)
            @info "skip existing $(axis) vector: $(vector_name)"
            continue
        end

        vector = get_vector(source, axis, vector_name).array
        if endswith(vector_name, "_module")
            vector = vector .+ 1
            set_vector!(source, axis, vector_name, vector; overwrite = true)
        end

        if empty !== nothing
            eltype = typeof(empty)
        else
            eltype = nothing
        end

        if rename == vector_name
            @info "copy $(axis) vector: $(vector_name)"
        else
            @info "copy $(axis) vector: $(vector_name) to: $(rename)"
        end

        copy_vector!(;  # NOJET
            destination = daf,
            source,
            axis,
            name = vector_name,
            rename,
            eltype,
            empty,
            bestify,
            min_sparse_saving_fraction,
            overwrite,
            insist,
        )
    end

    return nothing
end

function import_matrices_data(
    daf::DafWriter,
    source::DafReader,
    rows_axis::AbstractString,
    columns_axis::AbstractString,
    base_copy_data::CopyAnnData;
    copy_data::Maybe{CopyAnnData},
    bestify::Bool,
    min_sparse_saving_fraction::AbstractFloat,
    overwrite::Bool,
    insist::Bool,
)::Nothing
    for matrix_name in matrices_set(source, rows_axis, columns_axis; relayout = false)
        data = nothing

        if copy_data !== nothing
            data = get(copy_data, matrix_name, nothing)
        end

        if data === nothing && (copy_data === nothing || !haskey(copy_data, matrix_name))
            data = get(base_copy_data, matrix_name, nothing)

            if data === nothing && haskey(base_copy_data, "*") && !haskey(base_copy_data, matrix_name)
                data = (matrix_name, nothing)
            end
        end

        if data === nothing
            @info "skip $(rows_axis)-$(columns_axis) matrix: $(matrix_name)"
            continue
        end

        rename, empty = data

        relayout = !has_matrix(source, columns_axis, rows_axis, matrix_name; relayout = false)

        if !overwrite && !insist && has_matrix(daf, rows_axis, columns_axis, rename; relayout)
            @info "skip existing $(rows_axis)-$(columns_axis) matrix: $(matrix_name) ($(relayout ? "" : "!")relayout)"
            continue
        end

        if empty !== nothing
            eltype = typeof(empty)
        else
            eltype = nothing
        end

        if rename == matrix_name
            @info "copy $(rows_axis)-$(columns_axis) matrix: $(matrix_name) ($(relayout ? "" : "!")relayout)"
        else
            @info "copy $(rows_axis)-$(columns_axis) matrix: $(matrix_name) to: $(rename) ($(relayout ? "" : "!")relayout)"
        end

        copy_matrix!(;  # NOJET
            destination = daf,
            source,
            rows_axis,
            columns_axis,
            name = matrix_name,
            rename,
            eltype,
            empty,
            bestify,
            min_sparse_saving_fraction,
            overwrite,
            insist,
            relayout,
        )
    end

    return nothing
end

function copy_metacells_of_cells(daf::DafWriter, cells_h5ad::AbstractString; overwrite::Bool, insist::Bool)::Nothing
    @debug "readh5ad $(cells_h5ad) {"
    cells_adata = readh5ad(cells_h5ad)
    @debug "readh5ad $(cells_h5ad) }"

    if axis_length(daf, "cell") == size(cells_adata, 1)
        cell_indices = nothing
    else
        cell_indices = axis_indices(daf, "cell", cells_h5ad.obs_names)
    end

    index_per_metacell = axis_dict(daf, "metacell")
    name_per_cell = cells_adata.obs_names
    is_outlier_per_cell = cells_adata.obs[!, "metacell"] .< 0

    for (vector, rename, be_outlier) in
        (("metacell_name", "metacell", false), ("most_similar_name", "most_similar.metacell", true))
        if !(vector in names(cells_adata.obs))
            if vector == "metacell_name"
                error("missing per-obs annotation: metacell_name\nin cells h5ad: $(cells_h5ad)")
            end
            continue
        end

        data = cells_adata.obs[!, vector]
        data[is_outlier_per_cell .!= be_outlier] .= ""

        for (cell_name, metacell_name) in zip(name_per_cell, data)
            if metacell_name != "" && !haskey(index_per_metacell, metacell_name)
                error(
                    "invalid metacell_name: $(metacell_name)\n" *
                    "for the cell: $(cell_name)\n" *
                    "in the cells h5ad: $(cells_h5ad)",
                )
            end
        end

        if !overwrite && !insist && has_vector(daf, "cell", rename)
            @info "skip cell vector: $(vector)"
            continue
        end

        if cell_indices !== nothing
            full_data = fill("", axis_length(daf, "cell"))
            full_data[cell_indices] .= data
            data = full_data
        end

        @info "copy cell vector: $(vector) to: $(rename)"
        set_vector!(daf, "cell", rename, data; overwrite)
    end
end

end  # module
