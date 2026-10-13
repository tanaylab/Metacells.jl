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
using GZip
using HTTP
using Serialization
using TanayLabUtilities

# The namespaces whose names have no pattern.
const FREE_FORM_NAMESPACES = ("GeneSymbol",)

# The pattern of the names in each of the other namespaces.
const NAMESPACE_PATTERNS = Dict(
    "EnsemblGene" => raw"(ENS[A-Z]*G|ASMPATCHG)[0-9]{11}",
    "EnsemblTranscript" => raw"(ENS[A-Z]*T|ASMPATCHT)[0-9]{11}",
    "EnsemblProtein" => raw"(ENS[A-Z]*P|ASMPATCHP)[0-9]{11}",
    "HGNC" => raw"HGNC:[0-9]+",
    "MGI" => raw"MGI:[0-9]+",
    "RefSeq" => raw"(NM|NR|NP|XM|XR|XP|YP)_[0-9]+",
    "UCSC" => raw"uc[0-9]{3}[a-z]{3}\.[0-9]+",
)

# The namespaces whose names may carry a `.<version>` suffix, which is not part of the identity of the gene.
const VERSIONED_NAMESPACES = ("EnsemblGene", "EnsemblTranscript", "EnsemblProtein", "RefSeq")

# The namespaces whose names may carry a `_PAR_Y` suffix (after the version). GENCODE uses it for the copy of a gene in
# the pseudoautosomal region of chromosome Y, while Ensembl has a single gene for both copies.
const PAR_NAMESPACES = ("EnsemblGene", "EnsemblTranscript", "EnsemblProtein")

# A name in a namespace with a pattern is the pattern, then a `.<version>` suffix in the `VERSIONED_NAMESPACES`, then a
# `_PAR_Y` suffix in the `PAR_NAMESPACES`, then junk which contains no letters or digits (such as a trailing `;`). Only
# the pattern is the actual name.
const NAMESPACE_REGEXES = Dict(
    namespace => Regex(
        "^(" *
        pattern *
        ")" *
        (namespace in VERSIONED_NAMESPACES ? raw"(\.[0-9]+)?" : "") *
        (namespace in PAR_NAMESPACES ? "(_PAR_Y)?" : "") *
        raw"[^A-Za-z0-9]*$",
    ) for (namespace, pattern) in NAMESPACE_PATTERNS
)

"""
    normalize_gene_name(
        name::AbstractString;
        namespace::AbstractString,
        validate::Bool = false
    )::AbstractString

Normalize a gene name in some namespace. Each namespace other than `GeneSymbol` has a pattern for its names (for
example, `ENSG00000124766` for `EnsemblGene`). In such a namespace, this removes the `.[0-9]` version suffix (in the
Ensembl and RefSeq namespaces), the `_PAR_Y` suffix (in the Ensembl namespaces, which GENCODE uses for the copy of a
gene in the pseudoautosomal region of chromosome Y), and any trailing junk which contains no letters or digits (such as
a trailing `;`). If the name does not match the pattern, then if `validate` this fails, otherwise this returns the name
unchanged. To lookup a name in a list or a namespace, you need to normalize the query gene name accordingly. In other
namespaces, the `.[0-9]` suffix is part of the name. For example, `AC092171.3` and `AC092171.5` are different genes. An
unknown namespace is an error.

The case of a name is never changed, as Gmara names are in mixed case, differing between namespaces; for example,
`SOX4` for the human gene and `Sox4` for the mouse one.
"""
function normalize_gene_name(name::AbstractString; namespace::AbstractString, validate::Bool = false)::AbstractString
    if namespace in FREE_FORM_NAMESPACES
        return name
    end
    regex = get(NAMESPACE_REGEXES, namespace, nothing)
    if regex === nothing
        error("unknown Gmara namespace: $(namespace)")
    end
    name_match = match(regex, name)
    if name_match === nothing
        @assert !validate "invalid $(namespace) gene name: $(name)"
        return name
    end
    return name_match.captures[1]
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
        [namespace::AbstractString = $(DEFAULT.namespace),
        list::Maybe{AbstractString} = $(DEFAULT.list),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT,
    )::AbstractSet{<:AbstractString}

Return the set of names of a `version` of a `list` in a `namespace` of genes of some `species`. This returns all the
names that are (probably) in the list; it a name isn't in the result, it is almost certain it does not belong in the
list. As usual in Gmara, this includes everything that may be used as name, e.g. for Ensembl it includes genes,
transcripts and proteins; for Symbol it includes approved and alises; etc. If the `list` is `nothing`, this just returns
the set of known gene names in the `namespace`.
"""
@logged :mcs_ops @documented function gmara_genes(;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    list::Maybe{AbstractString} = nothing,
    version::AbstractString = "main",
    cache_dir::AbstractString = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::AbstractSet{<:AbstractString}
    if list === nothing
        path = "genes/$(species)/namespaces/names/$(namespace).tsv"
    else
        path = "genes/$(species)/lists/$(list)/names/$(namespace).tsv"
    end
    return get_set_file(path; version, cache_dir, timeout)
end

function get_set_file(
    path::AbstractString;
    version::AbstractString,
    cache_dir::AbstractString,
    timeout::Real,
)::AbstractSet{<:AbstractString}
    cache_key = version * "/" * path
    set = lock_read(CACHE_LOCK) do
        return get(CACHE_IN_MEMORY, cache_key, nothing)
    end
    if set === nothing
        set = lock(CACHE_LOCK) do
            set = get(CACHE_IN_MEMORY, cache_key, nothing)
            if set === nothing
                set = load_set_file(path; version, cache_dir, timeout)
                CACHE_IN_MEMORY[cache_key] = set
            end
            return set
        end
    end
    return set
end

function load_set_file(
    path::AbstractString;
    version::AbstractString,
    cache_dir::AbstractString,
    timeout::Real,
)::AbstractSet{<:AbstractString}
    cache_version_path = cache_dir * "/" * version * "/" * path
    with_lock_file(cache_version_path; timeout) do
        etag_path = cache_version_path * ".etag"
        if ispath(etag_path)
            etag = read(etag_path, String)
        else
            etag = nothing
        end

        if ensure_is_downloaded(path; version, etag, cache_dir) || !ispath(cache_version_path * ".jl_set.gz")
            set = write_set_to_file(path; version, cache_dir)
        else
            set = read_set_from_file(path; version, cache_dir)
        end

        return set
    end
end

function ensure_is_downloaded(
    path::AbstractString;
    version::AbstractString,
    etag::Maybe{AbstractString},
    cache_dir::AbstractString,
)::Bool
    url = github_url(path; version)

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
        error("Error GET $(url): HTTP status $(response.status)")  # UNTESTED
    end

    response_headers = Dict(response.headers)
    cache_data_path = cache_dir * "/" * version * "/" * path * ".gz"
    # The cache holds the data gzipped, so a body the server has already gzipped is written as it came - `decompress =
    # false` above is what keeps it that way - and any other body is compressed on the way in.
    if get(response_headers, "Content-Encoding", nothing) == "gzip"
        open(cache_data_path, "w") do file
            return write(file, response.body)
        end
    else
        GZip.open(cache_data_path, "w") do file  # UNTESTED
            return write(file, response.body)  # UNTESTED
        end
    end

    cache_etag_path = cache_dir * "/" * version * "/" * path * ".etag"
    etag = get(response_headers, "Etag", nothing)
    if etag !== nothing
        open(cache_etag_path, "w") do file
            print(file, etag)
            return nothing
        end
    elseif ispath(cache_etag_path)  # UNTESTED
        rm(cache_etag_path)  # UNTESTED
    end

    return true
end

function github_url(path::AbstractString; version::AbstractString = "main")::AbstractString
    # If Gmara switches back to using LFS, this may need to change to:
    # `https://media.githubusercontent.com/media/tanaylab/Gmara/$(version)/$(path)`
    return "https://raw.githubusercontent.com/tanaylab/Gmara/$(version)/$(path)"
end

function write_set_to_file(
    path::AbstractString;
    version::AbstractString,
    cache_dir::AbstractString,
)::AbstractSet{<:AbstractString}
    set = GZip.open(cache_dir * "/" * version * "/" * path * ".gz") do file
        text = read(file, String)
        lines = split(text, "\n")
        last_line = pop!(lines)
        @assert last_line == ""
        return Set([split(line, "\t")[1] for line in lines])
    end

    GZip.open(cache_dir * "/" * version * "/" * path * ".jl_set.gz", "w") do file
        return serialize(file, set)  # NOJET
    end

    return set
end

function read_set_from_file(
    path::AbstractString;
    version::AbstractString,
    cache_dir::AbstractString,
)::AbstractSet{<:AbstractString}
    return GZip.open(cache_dir * "/" * version * "/" * path * ".jl_set.gz") do file
        return deserialize(file)
    end
end

function with_lock_file(action::Function, path::AbstractString; timeout::Real)::Any
    file = lock_file(path; timeout)
    try
        return action()
    finally
        unlock_file(file, path)
    end
end

function lock_file(path::AbstractString; timeout::Real)::Base.Filesystem.File
    mkpath(dirname(path))  # NOJET
    start_time = time()
    lock_path = path * ".lock"
    while true
        try
            return Base.Filesystem.open(lock_path, Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL)
        catch exception
            if exception isa Base.IOError  # UNTESTED
                if timeout > 0  # UNTESTED
                    elapsed_time = time() - start_time  # UNTESTED
                    @assert elapsed_time < timeout "Waiting for more than $(timeout) seconds for $(lock_path)"  # UNTESTED
                end
                sleep(1)  # UNTESTED
            else
                rethrow()  # UNTESTED
            end
        end
    end
end

function unlock_file(file::Base.Filesystem.File, path::AbstractString)::Nothing  # UNTESTED
    close(file)
    return rm(path * ".lock")
end

"""
    set_gmara_genes_mask!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString,
        [list::Maybe{AbstractString} = $(DEFAULT.list),
        gene_name_property::AbstractString = $(DEFAULT.gene_name_property),
        property::Maybe{AbstractString} = $(DEFAULT.property),
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT,
        overwrite::Bool = $(DEFAULT.overwrite)]
    )::Integer

Set a gene property mask in `daf` based on some `version` of a Gmara `list` of some `namespace` for some `species`. We
match the `gene_name_property` (by default, just the unique `name` in the gene axis) with the list names and set the
result mask as a per-gene `property` (by default, ``is_``_list_). If `list` is `nothing`, this just marks the gene names
that exist in the namespace. If `overwrite`, this will overwrite an existing property of the.array same name.

Returns the number of genes that are set in the mask.
"""
@logged :mcs_ops @documented function set_gmara_genes_mask!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = function_default(gmara_genes, :namespace),
    list::Maybe{AbstractString} = function_default(gmara_genes, :list),
    version::AbstractString = function_default(gmara_genes, :version),
    gene_name_property::AbstractString = "name",
    property::Maybe{AbstractString} = nothing,
    cache_dir = CACHE_DIR,
    timeout::Real = TIMEOUT,
    overwrite::Bool = false,
)::Integer
    if property === nothing
        property = "is_$(list)"
    end
    gene_names = get_vector(daf, "gene", gene_name_property).array
    gmara_names = gmara_genes(; species, namespace, list, version, cache_dir, timeout)
    is_in_list_per_gene = [normalize_gene_name(gene_name; namespace) in gmara_names for gene_name in gene_names]
    if !any(is_in_list_per_gene)
        error("empty $(property) genes mask loaded from Gmara")  # UNTESTED
    end
    set_vector!(daf, "gene", property, is_in_list_per_gene; overwrite)
    return sum(is_in_list_per_gene)
end

end  # module
