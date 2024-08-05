using Pkg

for pkg in (
    "Aqua",
    "Coverage",
    "Documenter",
    "JET",
    "JuliaFormatter",
    "LanguageServer",
    "Logging",
    "LoggingExtras",
    "StaticLint",
    "SymbolServer",
)
    println("Adding $(pkg):")
    Pkg.add(pkg)
end
