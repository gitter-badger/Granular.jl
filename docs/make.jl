using Documenter, SeaIce

makedocs(
    modules = [SeaIce],
    clean = false,
    format = :html,
    sitename = "SeaIce.jl",
    authors = "Anders Damsgaard",
    pages = Any[ # Compat: `Any` for 0.4 compat
        "Home" => "index.md",
        "Manual" => Any[
            "man/installation.md",
            "man/simple_example.md",
        ],
        "Library" => Any[
            "Public API" => "lib/public.md",
            hide("Internals" => "lib/internals.md", Any[
              "lib/internals.md",
             ])
        ]
    ],
)

deploydocs(
    repo = "github.com/anders-dc/SeaIce.jl.git",
    julia = "0.6",
    deps = nothing,
    make = nothing,
    target = "build",
)
