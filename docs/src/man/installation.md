# Installation

## Stable installation (recommended)
The latest stable release of Granular.jl can be installed directly from the 
Julia shell by:

```julia-repl
julia> Pkg.add("Granular")
```

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/Granular` and install its requirements.  The 
package [JLD](https://github.com/JuliaIO/JLD.jl) is used for model restarts and 
is recommended but not required, and is thus not automatically installed.

If desired, the current developmental version of the [Granular.jl Github 
repository](https://github.com/anders-dc/Granular.jl) can be installed with the 
command:

```julia-repl
julia> Pkg.clone("git://github.com/anders-dc/Granular.jl")
```

*Please note:* The developmental version is considered unstable and should only 
be used over the stable version if there is a compelling reason to do so.

## Package tests
The Granular.jl package contains many tests that verify that the functionality 
works as intended.  The extent of test coverage of the source code is monitored 
and published with [CodeCov](https://codecov.io/gh/anders-dc/Granular.jl).

The package tests are during development continuously run with 
[Travis-CI](https://travis-ci.org/anders-dc/Granular.jl) for Mac (latest stable 
release) and Linux (Ubuntu stable (trusty)), and 
[AppVeyor](https://ci.appveyor.com/project/anders-dc/seaice-jl) for Windows.

The test scripts are contained in the `test/` directory, can be run locally 
with the following command:

```julia-repl
julia> Pkg.test("Granular")
```

In case any of these tests fail, please open a [Github 
Issue](https://github.com/anders-dc/Granular.jl/issues) so it can be 
investigated and diagnosed further.
