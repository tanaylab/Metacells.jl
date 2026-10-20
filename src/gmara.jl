"""
Access gene names lists from [Gmara](https://github.com/tanaylab/Gmara).

!!! note

    All the functions here are thread-safe and can also be invoked from multiple parallel processes using the same
    [`CACHE_DIR`](@ref). This directory can even be shared between multiple users as long as they have read-write
    permissions to the shared directory. This should even work on NFS mounted volumes shared between multiple servers.
"""
module Gmara

export EnsemblConsistency
export GmaraEnsemblGenes
export collect_ensembl_genes
export empty_gmara_cache!
export fetch_gmara_ensembl_genes
export gmara_genes
export gmara_genes_mask
export normalize_gene_name
export set_gmara_genes_mask!
export verify_consistent_ensembl_genes

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
    name_match = match(regex, String(name))
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
for serialized Julia set objects, `.jl_dict.gz` for serialized Julia dictionary objects, and `.lock` for temporary lock
files for coordinating between parallel processes.
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

const CACHE_IN_MEMORY = Dict{AbstractString, Any}()

const CACHE_LOCK = ReadWriteLock()

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
    return get_parsed_file(path; kind = "set", version, cache_dir, timeout)
end

"""
    struct GmaraEnsemblGenes
        ensembls_of_name::AbstractDict{<:AbstractString, <:AbstractSet{<:AbstractString}}
        final_ensembls_of_ensembl::AbstractDict{<:AbstractString, <:AbstractSet{<:AbstractString}}
    end

The mappings of gene names to EnsemblGenes of some version of Gmara for some species, as returned by
[`fetch_gmara_ensembl_genes`](@ref). These are mappings rather than per-axis-entry vectors, so their names use `of`
rather than `per`.

A final EnsemblGene is the end of a chain of replacements of retired identifiers. It is a current gene, or a retired
identifier which has no replacement. It is the best identifier there is for the gene.

  - `ensembls_of_name` maps each name of some namespace to its final EnsemblGenes. A name which Gmara knows, but which
    has no EnsemblGene, maps to an empty set.
  - `final_ensembls_of_ensembl` maps each EnsemblGene identifier Gmara knows to its final EnsemblGenes. A final
    identifier maps to itself. A retired identifier which was replaced maps to the final EnsemblGenes it was replaced by
    (several, if it was split).
"""
struct GmaraEnsemblGenes
    ensembls_of_name::AbstractDict{<:AbstractString, <:AbstractSet{<:AbstractString}}
    final_ensembls_of_ensembl::AbstractDict{<:AbstractString, <:AbstractSet{<:AbstractString}}
end

"""
    fetch_gmara_ensembl_genes(;
        species::AbstractString,
        [namespace::AbstractString = $(DEFAULT.namespace),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT],
    )::GmaraEnsemblGenes

Fetch the mappings of gene names to EnsemblGenes of some `version` of Gmara for some `species`. This returns both the
EnsemblGenes of each name of the `namespace`, and the final EnsemblGenes of each EnsemblGene identifier, so both always
come from the same species and version. See [`GmaraEnsemblGenes`](@ref).
"""
@logged :mcs_ops @documented function fetch_gmara_ensembl_genes(;
    species::AbstractString,
    namespace::AbstractString = "GeneSymbol",
    version::AbstractString = "main",
    cache_dir::AbstractString = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::GmaraEnsemblGenes
    ensembls_of_name = get_parsed_file(
        "genes/$(species)/namespaces/names/$(namespace).tsv";
        kind = "dict",
        version,
        cache_dir,
        timeout,
    )
    final_ensembls_of_ensembl =
        get_parsed_file("genes/$(species)/namespaces/names/EnsemblGene.tsv"; kind = "dict", version, cache_dir, timeout)
    return GmaraEnsemblGenes(ensembls_of_name, final_ensembls_of_ensembl)
end

"""
    collect_ensembl_genes(;
        names::AbstractVector{<:AbstractString},
        gmara_ensembl_genes::GmaraEnsemblGenes,
    )::AbstractVector{<:AbstractString}

Return the EnsemblGene of each of the `names`, using the `ensembls_of_name` of the `gmara_ensembl_genes` (see
[`fetch_gmara_ensembl_genes`](@ref)). This is the one EnsemblGene the name maps to, or empty if it maps to none or to
several. The names are looked up as given, so in a namespace other than `GeneSymbol` they should be normalized first
(see [`normalize_gene_name`](@ref)).

This can be used to set [`vector_of_ensembl_per_gene`](@ref Metacells.Contracts.vector_of_ensembl_per_gene) when the
data has no Ensembl identifiers.
"""
@logged :mcs_ops @documented function collect_ensembl_genes(;
    names::AbstractVector{<:AbstractString},
    gmara_ensembl_genes::GmaraEnsemblGenes,
)::AbstractVector{<:AbstractString}
    ensembls_of_name = gmara_ensembl_genes.ensembls_of_name
    return [single_ensembl_of_name(ensembls_of_name, name) for name in names]
end

function single_ensembl_of_name(
    ensembls_of_name::AbstractDict{<:AbstractString, <:AbstractSet{<:AbstractString}},
    name::AbstractString,
)::AbstractString
    ensembls = get(ensembls_of_name, name, nothing)
    if ensembls === nothing || length(ensembls) != 1
        return ""
    end
    return first(ensembls)
end

"""
    struct EnsemblConsistency
        name::AbstractString
        ensembl::AbstractString
        normalized_ensembl::Maybe{AbstractString}
        name_ensembls::AbstractSet{<:AbstractString}
        final_ensembls::AbstractSet{<:AbstractString}
    end

The consistency of the Ensembl identifier of a gene name with the mappings of Gmara, as returned by
[`verify_consistent_ensembl_genes`](@ref).

  - `name` and `ensembl` are as given.
  - `normalized_ensembl` is the `ensembl` without its suffixes (see [`normalize_gene_name`](@ref)), or `nothing` if
    it is not in the form of an Ensembl identifier.
  - `name_ensembls` are the final EnsemblGenes Gmara maps the `name` to (see [`GmaraEnsemblGenes`](@ref)).
  - `final_ensembls` are the final EnsemblGenes Gmara maps the `normalized_ensembl` to.

It also has computed Boolean properties. The first two are true when all is well:

  - `is_consistent`: none of the inconsistency properties listed below is true, so the `ensembl` identifies the gene
    of the `name`. The identifier may still have suffixes (such as a version).
  - `ensembl_is_normalized`: the `ensembl` is in the form of an Ensembl identifier, and has no suffixes.

The rest say what the inconsistency is. Several of them may be true at once (e.g. for a retired identifier which was
replaced by a different gene):

  - `name_is_empty`: the `name` is empty.
  - `name_is_unknown`: the `name` is not empty, and Gmara maps it to no EnsemblGene.
  - `ensembl_is_empty`: the `ensembl` is empty.
  - `ensembl_is_invalid`: the `ensembl` is not empty, and is not in the form of an Ensembl identifier.
  - `ensembl_is_unknown`: Gmara does not know the `normalized_ensembl`.
  - `ensembl_is_retired`: Gmara maps the `normalized_ensembl` to other final EnsemblGenes.
  - `ensembl_is_different`: the `name_ensembls` and the `final_ensembls` are not empty, and have no EnsemblGene in
    common.
"""
struct EnsemblConsistency
    name::AbstractString
    ensembl::AbstractString
    normalized_ensembl::Maybe{AbstractString}
    name_ensembls::AbstractSet{<:AbstractString}
    final_ensembls::AbstractSet{<:AbstractString}
end

# The computed properties of an `EnsemblConsistency` which say what the inconsistency is.
const ENSEMBL_INCONSISTENCY_PROPERTIES = (
    :name_is_empty,
    :name_is_unknown,
    :ensembl_is_empty,
    :ensembl_is_invalid,
    :ensembl_is_unknown,
    :ensembl_is_retired,
    :ensembl_is_different,
)

function Base.getproperty(consistency::EnsemblConsistency, property::Symbol)::Any
    name = getfield(consistency, :name)
    ensembl = getfield(consistency, :ensembl)
    normalized_ensembl = getfield(consistency, :normalized_ensembl)
    name_ensembls = getfield(consistency, :name_ensembls)
    final_ensembls = getfield(consistency, :final_ensembls)
    if property == :is_consistent
        return !any(getproperty(consistency, inconsistency) for inconsistency in ENSEMBL_INCONSISTENCY_PROPERTIES)
    elseif property == :ensembl_is_normalized
        return normalized_ensembl !== nothing && normalized_ensembl == ensembl
    elseif property == :name_is_empty
        return name == ""
    elseif property == :name_is_unknown
        return name != "" && isempty(name_ensembls)
    elseif property == :ensembl_is_empty
        return ensembl == ""
    elseif property == :ensembl_is_invalid
        return ensembl != "" && normalized_ensembl === nothing
    elseif property == :ensembl_is_unknown
        return normalized_ensembl !== nothing && isempty(final_ensembls)
    elseif property == :ensembl_is_retired
        return !isempty(final_ensembls) && final_ensembls != Set([normalized_ensembl])
    elseif property == :ensembl_is_different
        return !isempty(name_ensembls) && !isempty(final_ensembls) && isempty(intersect(name_ensembls, final_ensembls))
    else
        return getfield(consistency, property)
    end
end

function Base.propertynames(::EnsemblConsistency, private::Bool = false)::Tuple
    return (
        fieldnames(EnsemblConsistency)...,
        :is_consistent,
        :ensembl_is_normalized,
        ENSEMBL_INCONSISTENCY_PROPERTIES...,
    )
end

# Return the `ensembl` without its suffixes, or `nothing` if it is not in the form of an Ensembl identifier.
function normalized_ensembl_gene(ensembl::AbstractString)::Maybe{AbstractString}
    ensembl_match = match(NAMESPACE_REGEXES["EnsemblGene"], String(ensembl))
    if ensembl_match === nothing
        return nothing
    end
    return ensembl_match.captures[1]
end

"""
    verify_consistent_ensembl_genes(;
        names::AbstractVector{<:AbstractString},
        ensembl_per_name::AbstractVector{<:AbstractString},
        gmara_ensembl_genes::GmaraEnsemblGenes,
    )::AbstractVector{EnsemblConsistency}

Verify that the Ensembl identifier of each of the `names` is consistent with the `gmara_ensembl_genes` (see
[`fetch_gmara_ensembl_genes`](@ref)). Returns the [`EnsemblConsistency`](@ref) of each name, in the same order. This
only lists the inconsistencies; it is up to the caller to decide what to do about them.

The names are looked up as given, so in a namespace other than `GeneSymbol` they should be normalized first (see
[`normalize_gene_name`](@ref)).

This does not verify that the Ensembl identifiers are unique. To list the identifiers given to more than one name (and
the indices of these names), use `indices_of_duplicates` from `TanayLabUtilities`, ignoring the empty identifiers:
`indices_of_duplicates(ensembl_per_name; empty_value = "")`.
"""
@logged :mcs_ops @documented function verify_consistent_ensembl_genes(;
    names::AbstractVector{<:AbstractString},
    ensembl_per_name::AbstractVector{<:AbstractString},
    gmara_ensembl_genes::GmaraEnsemblGenes,
)::AbstractVector{EnsemblConsistency}
    @assert length(names) == length(ensembl_per_name)
    return [ensembl_consistency(name, ensembl, gmara_ensembl_genes) for (name, ensembl) in zip(names, ensembl_per_name)]
end

const NO_ENSEMBLS = Set{String}()

function ensembl_consistency(
    name::AbstractString,
    ensembl::AbstractString,
    gmara_ensembl_genes::GmaraEnsemblGenes,
)::EnsemblConsistency
    normalized_ensembl = normalized_ensembl_gene(ensembl)
    name_ensembls = get(gmara_ensembl_genes.ensembls_of_name, name, NO_ENSEMBLS)
    if normalized_ensembl === nothing
        final_ensembls = NO_ENSEMBLS
    else
        final_ensembls = get(gmara_ensembl_genes.final_ensembls_of_ensembl, normalized_ensembl, NO_ENSEMBLS)
    end
    return EnsemblConsistency(name, ensembl, normalized_ensembl, name_ensembls, final_ensembls)
end

# Return the parsed content of a Gmara names file, cached in memory and on disk. The `kind` says how to parse it: `set`
# for the set of names, or `dict` for the set of EnsemblGenes of each name.
function get_parsed_file(
    path::AbstractString;
    kind::AbstractString,
    version::AbstractString,
    cache_dir::AbstractString,
    timeout::Real,
)::Any
    cache_key = version * "/" * path * "." * kind
    cached = lock_read(CACHE_LOCK) do
        return get(CACHE_IN_MEMORY, cache_key, nothing)
    end
    if cached !== nothing
        return cached
    end
    return lock(CACHE_LOCK) do
        parsed = get(CACHE_IN_MEMORY, cache_key, nothing)
        if parsed === nothing
            parsed = load_parsed_file(path; kind, version, cache_dir, timeout)
            CACHE_IN_MEMORY[cache_key] = parsed
        end
        return parsed
    end
end

function load_parsed_file(
    path::AbstractString;
    kind::AbstractString,
    version::AbstractString,
    cache_dir::AbstractString,
    timeout::Real,
)::Any
    cache_version_path = cache_dir * "/" * version * "/" * path
    with_lock_file(cache_version_path; timeout) do
        etag_path = cache_version_path * ".etag"
        if ispath(etag_path)
            etag = read(etag_path, String)
        else
            etag = nothing
        end

        # The parsed file may be older than the downloaded data, if this was downloaded when parsing another kind.
        parsed_path = cache_version_path * ".jl_$(kind).gz"
        if ensure_is_downloaded(path; version, etag, cache_dir) ||
           !ispath(parsed_path) ||
           mtime(parsed_path) < mtime(cache_version_path * ".gz")
            parsed = write_parsed_to_file(path; kind, version, cache_dir)
        else
            parsed = read_parsed_from_file(path; kind, version, cache_dir)
        end

        return parsed
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

function write_parsed_to_file(
    path::AbstractString;
    kind::AbstractString,
    version::AbstractString,
    cache_dir::AbstractString,
)::Any
    parsed = GZip.open(cache_dir * "/" * version * "/" * path * ".gz") do file
        text = read(file, String)
        lines = split(text, "\n")
        last_line = pop!(lines)
        @assert last_line == ""
        header = popfirst!(lines)
        @assert startswith(header, "name\tensembl_gene\t") "unexpected header of: $(path)"
        if kind == "set"
            return Set{String}([split(line, "\t")[1] for line in lines])
        else
            @assert kind == "dict"
            return parse_ensembls_of_name(lines)
        end
    end

    GZip.open(cache_dir * "/" * version * "/" * path * ".jl_$(kind).gz", "w") do file
        return serialize(file, parsed)  # NOJET
    end

    return parsed
end

# Return the set of EnsemblGenes of each name, given the lines of a names file. A name whose only EnsemblGenes are
# placeholders (`ENS!...`, for names with no EnsemblGene) maps to an empty set.
function parse_ensembls_of_name(lines::AbstractVector{<:AbstractString})::Dict{String, Set{String}}
    ensembls_of_name = Dict{String, Set{String}}()
    for line in lines
        fields = split(line, "\t")
        ensembls = get!(Set{String}, ensembls_of_name, fields[1])
        if !startswith(fields[2], "ENS!")
            push!(ensembls, fields[2])
        end
    end
    return ensembls_of_name
end

function read_parsed_from_file(
    path::AbstractString;
    kind::AbstractString,
    version::AbstractString,
    cache_dir::AbstractString,
)::Any
    return GZip.open(cache_dir * "/" * version * "/" * path * ".jl_$(kind).gz") do file
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
    gmara_genes_mask(;
        names::AbstractVector{<:AbstractString},
        species::AbstractString,
        [namespace::AbstractString = $(DEFAULT.namespace),
        list::Maybe{AbstractString} = $(DEFAULT.list),
        version::AbstractString = $(DEFAULT.version),
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT],
    )::AbstractVector{Bool}

Return a mask with, for each of the `names` in the `namespace`, whether it is in some `version` of a Gmara `list` for
some `species`. If `list` is `nothing`, this marks the names that exist in the namespace. The names are normalized
before they are matched (see [`normalize_gene_name`](@ref)). For example, pass `axis_vector(daf, "gene")` as the
`names` of the `GeneSymbol` namespace, or the [`vector_of_ensembl_per_gene`](@ref
Metacells.Contracts.vector_of_ensembl_per_gene) as the `names` of the `EnsemblGene` namespace.

This is useful for comparing a mask with a Gmara list, without changing the data. To set the mask as a per-gene
property, use [`set_gmara_genes_mask!`](@ref).
"""
@logged :mcs_ops @documented function gmara_genes_mask(;
    names::AbstractVector{<:AbstractString},
    species::AbstractString,
    namespace::AbstractString = function_default(gmara_genes, :namespace),
    list::Maybe{AbstractString} = function_default(gmara_genes, :list),
    version::AbstractString = function_default(gmara_genes, :version),
    cache_dir::AbstractString = CACHE_DIR,
    timeout::Real = TIMEOUT,
)::AbstractVector{Bool}
    gmara_names = gmara_genes(; species, namespace, list, version, cache_dir, timeout)
    return [normalize_gene_name(name; namespace) in gmara_names for name in names]
end

"""
    set_gmara_genes_mask!(
        daf::DafWriter;
        species::AbstractString,
        namespace::AbstractString,
        [list::Maybe{AbstractString} = $(DEFAULT.list),
        name_property::Maybe{AbstractString} = $(DEFAULT.name_property),
        mask_property::Maybe{AbstractString} = $(DEFAULT.mask_property),
        cache_dir = CACHE_DIR,
        timeout::Real = TIMEOUT,
        overwrite::Bool = $(DEFAULT.overwrite)]
    )::Integer

Set a gene property mask in `daf` based on some `version` of a Gmara `list` of some `namespace` for some `species`. We
match the names of the genes in the `namespace` with the list names, and set the result as a per-gene `mask_property`
(by default, ``is_``_list_). If `list` is `nothing`, this just marks the gene names that exist in the namespace, and the
`mask_property` must be specified. If `overwrite`, this will overwrite an existing property of the same name.

The names of the genes in the `namespace` are the per-gene `name_property`. By default, this is `name` (the gene axis
entries) for the `GeneSymbol` namespace, and `ensembl` (see [`vector_of_ensembl_per_gene`](@ref
Metacells.Contracts.vector_of_ensembl_per_gene)) for the `EnsemblGene` namespace. It must be specified for any other
namespace. The names are normalized before they are matched (see [`normalize_gene_name`](@ref)).

Returns the number of genes that are set in the mask.
"""
@logged :mcs_ops @documented function set_gmara_genes_mask!(
    daf::DafWriter;
    species::AbstractString,
    namespace::AbstractString = function_default(gmara_genes, :namespace),
    list::Maybe{AbstractString} = function_default(gmara_genes, :list),
    version::AbstractString = function_default(gmara_genes, :version),
    name_property::Maybe{AbstractString} = nothing,
    mask_property::Maybe{AbstractString} = nothing,
    cache_dir = CACHE_DIR,
    timeout::Real = TIMEOUT,
    overwrite::Bool = false,
)::Integer
    name_property = effective_name_property(namespace, name_property)
    if mask_property === nothing
        if list === nothing
            error("no mask_property specified for the names of the namespace: $(namespace)")
        end
        mask_property = "is_$(list)"
    end
    names = plain_array(get_vector(daf, "gene", name_property))
    is_in_list_per_gene = gmara_genes_mask(; names, species, namespace, list, version, cache_dir, timeout)
    if !any(is_in_list_per_gene)
        error("empty $(mask_property) genes mask loaded from Gmara")  # UNTESTED
    end
    set_vector!(daf, "gene", mask_property, is_in_list_per_gene; overwrite)
    return sum(is_in_list_per_gene)
end

# Return the per-gene property holding the names of the genes in the `namespace`. If the `name_property` is specified,
# this is it. Otherwise, it is `name` for the `GeneSymbol` namespace and `ensembl` for the `EnsemblGene` namespace, and
# it is an error for any other namespace.
function effective_name_property(namespace::AbstractString, name_property::Maybe{AbstractString})::AbstractString
    if name_property !== nothing
        return name_property  # UNTESTED
    elseif namespace == "GeneSymbol"
        return "name"
    elseif namespace == "EnsemblGene"
        return "ensembl"
    else
        error("no name_property specified for the namespace: $(namespace)")
    end
end

end  # module
