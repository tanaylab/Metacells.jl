"""
Access gene names lists from [Gmara](https://github.com/tanaylab/Gmara).

!!! note

    All the functions here are thread-safe and can also be invoked from multiple parallel processes using the same
    [`CACHE_DIR`](@ref). This directory can even be shared between multiple users as long as they have read-write
    permissions to the shared directory. This should even work on NFS mounted volumes shared between multiple servers.
"""
module Gmara

export gmara_genes
export empty_gmara_cache!
export normalize_gene_name
export set_gmara_genes_mask!

using ConcurrentUtils
using DataAxesFormats
using DataAxesFormats.GenericTypes
using GZip
using HTTP
using Serialization

"""
    normalize_gene_name(name::AbstractString; namespace::AbstractString)::AbstractString

Normalize the a gene name in some namespace. In most namespaces, this means removing the `.[0-9]` version suffix from
the name, and converting the name to upper case. To lookup a name in a list or a namespace, you need to normalize the
query gene name accordingly. The UCSC namespace is an exception in that it is all-lower-case and the `.[0-9]` suffix
seems to be an inherent part of the identifier.
"""
function normalize_gene_name(name::AbstractString; namespace::AbstractString)::AbstractString
    if namespace == "UCSC"
        return name
    end
    parts = split(name, ".")
    if length(parts) == 2
        name = parts[1]
    end
    return name
end

"""
The default (`\$HOME/.cache/gmara`) location of the cache of downloaded Gmara data files.

You can override this by setting the `METACELLS_GMARA_CACHE_DIR` environment variable, or by passing
an explicit `cache_dir` parameter to the functions.

The top-level under this is the version indicator, where `main` is always the latest and greatest version. Under each
version we store the files in the same path as in github, with a `.gz` suffix for the compressed raw data, `.jl_set.gz`
for serialized Julia set objects, and `.lock` for temporary lock files for coordinating between parallel processes.
"""
CACHE_DIR::AbstractString = get(ENV, "METACELLS_GMARA_CACHE") do
    return ENV["HOME"] * "/.cache/gmara"
end

"""
The default timeout in seconds (10) for waiting for a lock file in the Gmara cache. If not positive, will wait forever.
If a process crashes very badly then a lock file may be left behind and may need to be removed by hand to allow access
for the data.

You can override this by setting the `METACELLS_GMARA_TIMEOUT` environment variable, or by passing an explicit `timeout`
parameter to the functions.
"""
TIMEOUT::Real = parse(Int, get(ENV, "METACELLS_GMARA_TIMEOUT", "10"))

CACHE_IN_MEMORY = Dict{AbstractString, Any}()

CACHE_LOCK = ReadWriteLock()

"""
    empty_gmara_cache!()::Nothing

All requests are cached in-memory. This makes repeated requests cheap. This consumes some (modest amount of) memory;
also, if the data in the server has been updated (which rarely happens), you will keep getting the old result. This
function releases all the memory and forces all subsequent requests to query the server. In the common case the server
tells us our disk cache data is up to date, we don't re-download it).
"""
function empty_gmara_cache!()::Nothing
    lock(CACHE_LOCK) do
        return empty!(CACHE_IN_MEMORY)
    end
    return nothing
end

"""
    gmara_genes(;
        species::AbstractString,
        namespace::AbstractString,
        [list::Maybe{AbstractString} = nothing,
        version::AbstractString = "main"
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT],
    )::AbstractSet{<:AbstractString}

Return the set of names of a `version` of a `list` in a `namespace` of genes of some `species`. This returns all the
names that are (probably) in the list; it a name isn't in the result, it is almost certain it does not belong in the
list. As usual in Gmara, this includes everything that may be used as name, e.g. for Ensembl it includes genes,
transcripts and proteins; for Symbol it includes approved and alises; etc. If the `list` is `nothing`, this just returns
the set of known gene names in the `namespace`.
"""
function gmara_genes(;
    species::AbstractString,
    namespace::AbstractString,
    list::Maybe{AbstractString} = nothing,
    version::AbstractString = "main",
    cache_dir = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::AbstractSet{<:AbstractString}
    if list === nothing
        path = "genes/$(species)/namespaces/names/$(namespace).tsv"
    else
        path = "genes/$(species)/lists/$(list)/names/$(namespace).tsv"
    end
    return get_set_file(path; version = version, cache_dir = cache_dir, timeout = timeout)
end

function get_set_file(
    path::AbstractString;
    version::AbstractString = "main",
    cache_dir::AbstractString = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::AbstractSet{<:AbstractString}
    cache_key = version * "/" * path
    set = lock_read(CACHE_LOCK) do
        return get(CACHE_IN_MEMORY, cache_key, nothing)
    end
    if set === nothing
        set = lock(CACHE_LOCK) do
            set = get(CACHE_IN_MEMORY, cache_key, nothing)
            if set === nothing
                set = load_set_file(path; version = version, cache_dir = cache_dir, timeout = timeout)
                CACHE_IN_MEMORY[cache_key] = set
            end
            return set
        end
    end
    return set
end

function load_set_file(
    path::AbstractString;
    version::AbstractString = "main",
    cache_dir::AbstractString = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::AbstractSet{<:AbstractString}
    cache_version_path = cache_dir * "/" * version * "/" * path
    with_lock_file(cache_version_path; timeout = timeout) do
        etag_path = cache_version_path * ".etag"
        if ispath(etag_path)
            etag = read(etag_path, String)
        else
            etag = nothing
        end

        if (
            ensure_is_downloaded(path; version = version, etag = etag, cache_dir = cache_dir) ||
            !ispath(cache_version_path * ".jl_set.gz")
        )
            set = write_set_to_file(path; version = version, cache_dir = cache_dir)
        else
            set = read_set_from_file(path; version = version, cache_dir = cache_dir)
        end

        return set
    end
end

function ensure_is_downloaded(
    path::AbstractString;
    version::AbstractString = "main",
    etag::Maybe{AbstractString},
    cache_dir::AbstractString = CACHE_DIR,
)::Bool
    url = github_url(path; version = version)

    headers = ["Accept-Encoding" => "gzip"]
    if etag !== nothing
        push!(headers, "If-None-Match" => etag)
    end

    response =  # NOJET
        HTTP.request("GET", url, headers; canonicalize_headers = true, status_exception = false, decompress = false)

    if response.status == 304
        return false
    end
    if response.status != 200
        error("Error GET $(url): HTTP status $(response.status)")  # untested
    end

    response_headers = Dict(response.headers)
    cache_data_path = cache_dir * "/" * version * "/" * path * ".gz"
    # TODO: Seems that LFS files are never served zipped.
    # Strange, you would think it would be a priotity for them.
    if get(response_headers, "Content-Encoding", nothing) == "gzip"
        open(cache_data_path, "w") do file  # untested
            return write(file, response.body)
        end
    else
        GZip.open(cache_data_path, "w") do file
            return write(file, response.body)
        end
    end

    cache_etag_path = cache_dir * "/" * version * "/" * path * ".etag"
    etag = get(response_headers, "Etag", nothing)
    if etag !== nothing
        open(cache_etag_path, "w") do file
            print(file, etag)
            return nothing
        end
    elseif ispath(cache_etag_path)  # untested
        rm(cache_etag_path)  # untested
    end

    return true
end

function github_url(path::AbstractString; version::AbstractString = "main")::AbstractString
    # PRE-LFS: return "https://raw.githubusercontent.com/tanaylab/Gmara/$(version)/$(path)"
    return "https://media.githubusercontent.com/media/tanaylab/Gmara/$(version)/$(path)"
end

function write_set_to_file(
    path::AbstractString;
    version = "main",
    cache_dir::AbstractString = CACHE_DIR,
)::AbstractSet{<:AbstractString}
    set = GZip.open(cache_dir * "/" * version * "/" * path * ".gz") do file  # NOJET
        text = read(file, String)
        lines = split(text, "\n")
        last_line = pop!(lines)
        @assert last_line == ""
        return Set([split(line, "\t")[1] for line in lines])
    end

    GZip.open(cache_dir * "/" * version * "/" * path * ".jl_set.gz", "w") do file
        return serialize(file, set)
    end

    return set
end

function read_set_from_file(
    path::AbstractString;
    version::AbstractString = "main",
    cache_dir::AbstractString = CACHE_DIR,
)::AbstractSet{<:AbstractString}
    return GZip.open(cache_dir * "/" * version * "/" * path * ".jl_set.gz") do file
        return deserialize(file)
    end
end

function with_lock_file(action::Function, path::AbstractString; timeout::Real = TIMEOUT)::Any
    file = lock_file(path; timeout = timeout)
    try
        return action()
    finally
        unlock_file(file, path)
    end
end

function lock_file(path::AbstractString; timeout::Real = TIMEOUT)::Base.Filesystem.File
    mkpath(dirname(path))  # NOJET
    start_time = time()
    lock_path = path * ".lock"
    while true
        try
            return Base.Filesystem.open(lock_path, Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL)
        catch exception
            if exception isa Base.IOError  # untested
                if timeout > 0  # untested
                    elapsed_time = time() - start_time  # untested
                    @assert elapsed_time < timeout "Waiting for more than $(timeout) seconds for $(lock_path)"  # untested
                end
                sleep(1)  # untested
            else
                rethrow()  # untested
            end
        end
    end
end

function unlock_file(file::Base.Filesystem.File, path::AbstractString)::Nothing
    close(file)
    return rm(path * ".lock")
end

"""
    set_gmara_genes_mask!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString,
        [list::Maybe{AbstractString} = nothing,
        gene_names::AbstractString = "name",
        property::Maybe{AbstractString} = nothing,
        overwrite::Bool = false,
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT],
    )::Nothing

Set a gene property mask in `daf` based on some `version` of a Gmara `list` of some `namespace` for some `species`. We
match the `gene_names` (by default, just the unique names in the gene axis) with the list names and set the result mask
as a per-gene `property` (by default, ``is_``_list_). If `list` is `nothing`, this just marks the gene names that exist
in the namespace. If `overwrite`, this will overwrite an existing property of the.array same name.
"""
function set_gmara_genes_mask!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString,
    list::Maybe{AbstractString} = nothing,
    version::AbstractString = "main",
    gene_names::AbstractString = "name",
    property::Maybe{AbstractString} = nothing,
    overwrite::Bool = false,
    cache_dir = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::Nothing
    if property === nothing
        property = "is_$(list)"
    end
    gene_names = get_vector(daf, "gene", gene_names).array
    gmara_names = gmara_genes(;
        species = species,
        namespace = namespace,
        list = list,
        version = version,
        cache_dir = cache_dir,
        timeout = timeout,
    )
    gene_mask = [normalize_gene_name(gene_name; namespace = namespace) in gmara_names for gene_name in gene_names]
    set_vector!(daf, "gene", property, gene_mask; overwrite = overwrite)
    return nothing
end

end  # module
