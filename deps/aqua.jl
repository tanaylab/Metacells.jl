push!(LOAD_PATH, ".")

using Aqua
using Metacells
Aqua.test_ambiguities([Metacells])
Aqua.test_all(Metacells; ambiguities = false, unbound_args = false, deps_compat = false, persistent_tasks = false)
