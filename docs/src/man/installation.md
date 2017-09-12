# Installation
SeaIce.jl can be installed directly from the Julia shell by:

```julia-repl
julia> Pkg.clone("git://github.com/anders-dc/SeaIce.jl.git")
```

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/SeaIce` and install its requirements.  The package [JLD](https://github.com/JuliaIO/JLD.jl) 
is used for model restarts and is recommended but not required, and is thus not 
automatically installed.

You can run the package tests, which are contained in the `test/` directory, with
the following command:

```julia-repl
julia> Pkg.test("SeaIce")
```

The package can be updated from this repository using:

```julia-repl
julia> Pkg.update("SeaIce")
```