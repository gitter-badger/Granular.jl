using Documenter, SeaIce

makedocs(
    modules = [SeaIce],
    clean   = false,
    format   = :html,
    sitename = "SeaIce.jl",
    authors = "Anders Damsgaard",
    pages    = Any[ # Compat: `Any` for 0.4 compat
      "Home" => "index.md",
      "Manual" => Any[
        "installation.md"
      ]
    ]
)

deploydocs(
    repo   = "github.com/anders-dc/SeaIce.jl.git",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
    target = "build"
)
