"""
The `Metacells.jl` package provides computational services for the [metacells](https://github.com/tanaylab/metacells)
package, using [Daf](https://github.com/tanaylab/DataAxesFormats.jl) to hold the data. In the future, we'll ideally
migrate all of the `metacells`package computations to this package, converting the Python package to a thin wrapper, and
provide a similar thin R wrapper to provide metacell analysis from R as well. For now,`Metacells.jl`only provides a
subset of the features of the Python`metacells`package, which requires users to convert data from`AnnData`(for the old
features) to `Daf` (to the new features).
"""
module Metacells

using Reexport

include("defaults.jl")
@reexport using .Defaults

include("gmara.jl")
@reexport using .Gmara

include("anndata_format.jl")
@reexport using .AnnDataFormat

include("contracts.jl")
@reexport using .Contracts

include("identify_genes.jl")
@reexport using .IdentifyGenes

include("programs.jl")
@reexport using .Programs

end # module
