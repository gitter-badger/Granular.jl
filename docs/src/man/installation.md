# Installation
Granular.jl can be installed directly from the Julia shell by:

```julia-repl
julia> Pkg.add("Granular")
```

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/Granular` and install its requirements.  The 
package [JLD](https://github.com/JuliaIO/JLD.jl) is used for model restarts and 
is recommended but not required, and is thus not automatically installed.

You can run the package tests, which are contained in the `test/` directory, with
the following command:

```julia-repl
julia> Pkg.test("Granular")
```
