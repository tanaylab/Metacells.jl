#!/bin/bash
set -e -o pipefail
python3 deps/modules_to_dot.py | dot -Tsvg > src/assets/modules.svg
julia --project=deps/document_env -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
JULIA_DEBUG="" julia --project=deps/document_env --color=no deps/document.jl
sed -i 's:<img src="assets/modules.svg":<embed src="assets/modules.svg":;s: on <span class="colophon-date" title="[^"]*">[^<]*</span>::;s:<:\n<:g' docs/v0.1.0/*html
rm -rf docs/*/*.{cov,jl}
