# Metacells - Single-cell RNA Sequencing Analysis

The `Metacells.jl` package provides computational services for the [metacells](https://github.com/tanaylab/metacells)
package, using [Daf](https://github.com/tanaylab/Daf.jl) to hold the data. In the future, we'll ideally migrate all of
the `metacells`package computations to this package, converting the Python package to a thin wrapper, and provide a
similar thin R wrapper to provide metacell analysis from R as well. For now,`Metacells.jl` only provides a subset of the
features of the Python`metacells`package, which requires users to convert data from`AnnData`(for the old features)
to `Daf` (to the new features).

See the [v0.1.0 documentation](https://tanaylab.github.io/Metacells.jl/v0.1.0) for details.

## Motivation

The Python `metacells` package uses C++ for some of the underlying computations. In addition, higher level computations
are performed directly in Python, with parallelism using multi-processing (based on `fork` to avoid GIL issues and use
shared memory for at least the inputs). Finally, the current Python `metacells` is based on the restricted `AnnData`
format.

This combination of technologies is inefficient (e.g., the use of parallelism is severely restricted, and using GPUs
would require a major overhaul of the code) and awkward (the representation of per-cell, per-metacell, and per-type data
requires a combination of multiple `AnnData` objects and CSV files). This forced tools like
[MCView](https://github.com/tanaylab/MCView) to implement their own storage mechanism.

To solve these problems, we will be (incrementally) migrating to an implementation where all the computations would be
implemented in Julia using `Daf`. This would allow us to make the best use of multiple CPUs (and, hopefully, GPUs), and
provide a coherent representation of the data in a single data set, which would allow `MCView` and similar tools to work
directly from this data set.

For the duration of the transition, we'll maintain the existing `metacells` implementation to minimize the impact on
existing users. New features will only be available in the new technology stack, and old features will be migrated (one
at a time) to it, until the new system will be complete and the old system will become deprecated.

This process will take a while (months - probably all of 2024), so "Keep calm and carry on" - just use the current
Python `metacells` package API for now. That said, you would probably want to convert the data from `AnnData` to `Daf`
at some point of your processing pipeline, in order to use the new features, and/or use the (updated) `MCView` tool.

## Installation

Just `Pkg.add("Metacells")`, like installing any other Julia package. However, until this package is complete, there's
little point in directly installing it. You should probably install the Python package instead (which will automatically
install the Julia package as a dependency).

To install the Python `metacells` [package](https://github.com/tanaylab/metacells), just `pip install metacells`, like
installing any other Python package.

TODO: To install the R wrappers...

## License (MIT)

Copyright © 2024 Weizmann Institute of Science

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
