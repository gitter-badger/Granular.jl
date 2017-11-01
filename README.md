# Granular

![Granular.jl 
logo](https://github.com/anders-dc/Granular.jl/raw/master/docs/src/assets/logo.gif)

A [Julia](https://julialang.org) package for granular mechanics.

| Documentation | Build Status (Linux/Mac) | Build Status (Win) | Test Coverage |
|:-------------:|:------------------------:|:------------------:|:-------------:|
|[![Granular.jl Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://anders-dc.github.io/Granular.jl/latest) | [![Build Status](https://travis-ci.org/anders-dc/Granular.jl.svg?branch=master)](https://travis-ci.org/anders-dc/Granular.jl) | [![Build Status](https://ci.appveyor.com/api/projects/status/github/anders-dc/Granular.jl?svg=true)](https://ci.appveyor.com/project/anders-dc/seaice-jl/) | [![codecov.io](http://codecov.io/github/anders-dc/Granular.jl/coverage.svg?branch=master)](http://codecov.io/github/anders-dc/Granular.jl?branch=master) |

## Installation
[Granular.jl](https://github.com/anders-dc/Granular.jl) is not yet registered 
in the [official Julia package repository](https://pkg.julialang.org), but can 
be installed directly from the Julia shell by:

```julia
julia> Pkg.clone("git://github.com/anders-dc/Granular.jl")
```

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/Granular`, and install the packages specified as 
[requirements](REQUIRE). You can run the package tests, which are contained in
the [test/ directory](test/), with the following command:

```julia
julia> Pkg.test("Granular")
```

The package can be updated from this repository using:

```julia
julia> Pkg.update("Granular")
```

For more information on installation and usage, please refer to the [documentation](https://anders-dc.github.io/Granular.jl/latest).

## Author
[Anders Damsgaard](https://adamsgaard.dk), Geophysical Fluid Dynamics Laboratory, Princeton University.
