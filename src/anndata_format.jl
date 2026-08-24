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

  - Any per-gene masks which name a type, which a metacells `h5ad` may hold, become a per-gene-per-type matrix
    using [`import_gene_masks_per_type!`](@ref).
"""
module AnnDataFormat

export CopyAnnData
export import_cells_h5ad!
export import_gene_masks_per_type!
export import_metacells_h5ad!

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
    "u" => ("umap_u", nothing),
    "v" => ("umap_v", nothing),
    "w" => ("umap_w", nothing),
    "x" => ("umap_x", nothing),
    "y" => ("umap_y", nothing),
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
        type_colors_csv::Maybe{AbstractString} = $(DEFAULT.type_colors_csv),
        empty_type::Maybe{EmptyImplicit} = $(DEFAULT.empty_type),
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

If a `type_colors_csv` is given, a `type` axis is created from it. The file must have exactly two columns: the types,
under whatever name, and then a column named `color`. The types are taken in the order they are listed in, which is
typically a meaningful one, and their colors are stored as a per-type `color` property. The file is the authority on
which types exist: it may list types no cell has, but a cell whose type it does not list is an error, in one or the
other of them.

This needs a per-cell `type` property to have been imported, which usually means saying so in the `copy_data`, the
column holding it rarely being called `type`. Data often spells "this cell has no type" as a value of its own -
`Outliers`, `Doublet`, and the like - and which values those are depends on the data set; list them as the
`empty_type` and they become the empty string, which is how `Daf` spells "no value". They are the only values exempt
from having to appear in the CSV.

!!! note

    It is common to manually call `reconstruct_axis!` on the result to create additional axes (e.g., if the cells were
    collected from a set of batches and some properties are actually per-batch).
"""
@logged :mcs_ops @documented function import_cells_h5ad!(
    daf::DafWriter;
    cells_h5ad::AbstractString,
    copy_data::Maybe{CopyAnnData} = nothing,
    type_colors_csv::Maybe{AbstractString} = nothing,
    empty_type::Maybe{EmptyImplicit} = nothing,
    bestify::Bool = true,
    min_sparse_saving_fraction::AbstractFloat = function_default(copy_matrix!, :min_sparse_saving_fraction),
    overwrite::Bool = false,
    insist::Bool = false,
)::Nothing
    if empty_type !== nothing && type_colors_csv === nothing
        error("specified an empty_type without a type_colors_csv, so there are no types for it to be empty of")
    end

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

    if type_colors_csv !== nothing
        import_type_colors_csv(daf, type_colors_csv, empty_type, overwrite)
    end

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
  - The `u`, `v`, `w`, `x` and `y` UMAP coordinates are renamed to `umap_u`, `umap_v`, `umap_w`, `umap_x` and `umap_y`.
  - The `type` property is copied. If "the" type property of the metacells is different, use `copy_data` to rename
    it to `type` to match the `Daf` naming convention.
  - All other vectors are copied as-is.

Per-metacell-per-metacell:

  - The `obs_outgoing_weights` matrix is renamed to `outgoing_weights` and always stored as `Float32`.
  - All other matrices are copied as-is.

!!! note

    It is common to manually call `reconstruct_type!` on the result to create a type axis.
"""
@logged :mcs_ops @documented function import_metacells_h5ad!(
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

    copy_metacells_of_cells(daf, cells_h5ad; bestify, overwrite, insist)

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
    import_gene_masks_per_type!(
        daf::DafWriter;
        type_axis::AbstractString = $(DEFAULT.type_axis),
    )::Nothing

Convert per-gene masks which name a type in their own name into a per-gene-per-type matrix. `AnnData` has only the two
axes, so a mask of genes per type has to be spelled as one property per type: `something_gene_of_Bcell`,
`something_gene_of_Tcell`, and so on. Each such set becomes a single `is_something` matrix of genes by types (with a
default of `false`), which is what it was all along.

The `type_axis` must exist, since it says which types there are; a type with no mask property simply has no genes in
it. In practice these masks come in a metacells `h5ad` rather than a cells one.
"""
@logged :mcs_ops @documented function import_gene_masks_per_type!(
    daf::DafWriter;
    type_axis::AbstractString = "type",
)::Nothing
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
        if !(mask_vector isa SparseVector{Bool})
            mask_vector = SparseVector{Bool}(mask_vector)
        end
        push!(mask_vectors, mask_vector)
    end

    mask_matrix::SparseMatrixCSC{Bool} = hcat(mask_vectors...)  # NOJET
    @debug "reconstruct gene-$(type_axis) matrix: is_$(prefix)" _group = :mcs_details
    set_matrix!(daf, "gene", type_axis, "is_$(prefix)", bestify(mask_matrix); relayout = true, overwrite = true)  # NOJET

    return nothing
end

# Create the type axis out of the CSV, which is the authority on which types there are, what they are called and in
# what order they are listed. The data only says which cells have which type, and is checked against the CSV rather
# than allowed to add to it: a type in the data and not in the CSV is a mistake in one of them, and saying which is
# more useful than silently accepting either.
function import_type_colors_csv(
    daf::DafWriter,
    type_colors_csv::AbstractString,
    empty_type::Maybe{EmptyImplicit},
    overwrite::Bool,
)::Nothing
    data_frame = CSV.read(type_colors_csv, DataFrame)  # NOJET
    column_names = names(data_frame)
    if length(column_names) != 2 || column_names[2] != "color"
        error(chomp("""
            expected two columns, the types and then "color": $(join(column_names, ", "))
            in the type colors csv: $(type_colors_csv)
            """))
    end

    type_names = string.(data_frame[:, 1])
    colors = string.(data_frame[:, 2])

    seen_type_names = Set{AbstractString}()
    repeated_type_names = AbstractString[]
    for type_name in type_names
        if type_name in seen_type_names
            push!(repeated_type_names, type_name)
        else
            push!(seen_type_names, type_name)
        end
    end
    if !isempty(repeated_type_names)
        error(chomp("""
            repeated type(s): $(join(unique(repeated_type_names), ", "))
            in the type colors csv: $(type_colors_csv)
            """))
    end

    # Whatever the data spells "this cell has no type" as - `Outliers`, `Doublet`, and so on - becomes the empty string
    # `Daf` uses for "no value". Which values those are is a property of the data set, so it has to be said rather than
    # guessed at, and they are the only values exempt from being listed in the CSV.
    empty_types = Reconstruction.set_of_empty_implicit(empty_type)
    type_per_cell =
        [type_name in empty_types ? "" : string(type_name) for type_name in get_vector(daf, "cell", "type").array]

    unlisted_type_names =
        [type_name for type_name in unique(type_per_cell) if type_name != "" && !(type_name in seen_type_names)]
    if !isempty(unlisted_type_names)
        error(chomp("""
            type(s) of cells: $(join(unlisted_type_names, ", "))
            missing from the type colors csv: $(type_colors_csv)
            """))
    end

    # Nothing is written until everything above has been verified, so a rejected CSV leaves the data as it was.
    add_axis!(daf, "type", type_names; overwrite)
    set_vector!(daf, "type", "color", colors; overwrite)
    set_vector!(daf, "cell", "type", type_per_cell; overwrite = true)

    return nothing
end

function import_scalars_data(daf::DafWriter, source::DafReader; overwrite::Bool, insist::Bool)::Nothing
    for scalar_name in scalars_set(source)
        if scalar_name == "__name__"
            @debug "skip scalar: $(scalar_name)" _group = :mcs_details
        else
            @debug "copy scalar: $(scalar_name)" _group = :mcs_details
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
            @debug "skip $(axis) vector: $(vector_name)" _group = :mcs_details
            continue
        end

        rename, empty = data

        if !overwrite && !insist && has_vector(daf, axis, rename)
            @debug "skip existing $(axis) vector: $(vector_name)" _group = :mcs_details
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
            @debug "copy $(axis) vector: $(vector_name)" _group = :mcs_details
        else
            @debug "copy $(axis) vector: $(vector_name) to: $(rename)" _group = :mcs_details
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
            @debug "skip $(rows_axis)-$(columns_axis) matrix: $(matrix_name)" _group = :mcs_details
            continue
        end

        rename, empty = data

        relayout = !has_matrix(source, columns_axis, rows_axis, matrix_name; relayout = false)

        if !overwrite && !insist && has_matrix(daf, rows_axis, columns_axis, rename; relayout)
            @debug "skip existing $(rows_axis)-$(columns_axis) matrix: $(matrix_name) ($(relayout ? "" : "!")relayout)" _group =
                :mcs_details
            continue
        end

        if empty !== nothing
            eltype = typeof(empty)
        else
            eltype = nothing
        end

        if rename == matrix_name
            @debug "copy $(rows_axis)-$(columns_axis) matrix: $(matrix_name) ($(relayout ? "" : "!")relayout)" _group =
                :mcs_details
        else
            @debug "copy $(rows_axis)-$(columns_axis) matrix: $(matrix_name) to: $(rename) ($(relayout ? "" : "!")relayout)" _group =
                :mcs_details
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

function copy_metacells_of_cells(
    daf::DafWriter,
    cells_h5ad::AbstractString;
    bestify::Bool,
    overwrite::Bool,
    insist::Bool,
)::Nothing
    cells_adata = flame_timed("readh5ad") do
        return readh5ad(cells_h5ad)  # NOJET
    end

    if axis_length(daf, "cell") == size(cells_adata, 1)
        cell_indices = nothing
    else
        cell_indices = axis_indices(daf, "cell", cells_h5ad.obs_names)
    end

    index_per_metacell = axis_dict(daf, "metacell")
    name_per_cell = cells_adata.obs_names
    is_base_outlier_per_cell = cells_adata.obs[!, "metacell"] .< 0

    for (vector, rename, be_outlier) in
        (("metacell_name", "metacell", false), ("most_similar_name", "most_similar.metacell", true))
        if !(vector in names(cells_adata.obs))
            if vector == "metacell_name"
                error("missing per-obs annotation: metacell_name\nin cells h5ad: $(cells_h5ad)")
            end
            continue
        end

        data = cells_adata.obs[!, vector]
        data[is_base_outlier_per_cell .!= be_outlier] .= ""

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
            @debug "skip cell vector: $(vector)" _group = :mcs_details
            continue
        end

        if cell_indices !== nothing
            full_data = fill("", axis_length(daf, "cell"))
            full_data[cell_indices] .= data
            data = full_data
        end

        @debug "copy cell vector: $(vector) to: $(rename)" _group = :mcs_details
        set_vector!(daf, "cell", rename, data; overwrite)
    end

    if overwrite || insist || !has_vector(daf, "cell", "is_base_outlier")
        if cell_indices !== nothing
            full_is_base_outlier_per_cell = fill(true, axis_length(daf, "cell"))
            full_is_base_outlier_per_cell[cell_indices] .= is_base_outlier_per_cell
            is_base_outlier_per_cell = full_is_base_outlier_per_cell
        end
        if bestify
            is_base_outlier_per_cell = TanayLabUtilities.bestify(is_base_outlier_per_cell)
        end
        @debug "set cell vector: is_base_outlier" _group = :mcs_details
        set_vector!(daf, "cell", "is_base_outlier", is_base_outlier_per_cell; overwrite)
    end
end

end  # module
