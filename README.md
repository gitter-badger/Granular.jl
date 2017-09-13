# SeaIce

![SeaIce.jl 
logo](https://github.com/anders-dc/SeaIce.jl/raw/master/docs/src/assets/logo.gif)

A [Julia](https://julialang.org) package for sea-ice granular mechanics.

| Documentation | Build Status (Linux/Mac) | Build Status (Win) | Test Coverage |
|:-------------:|:------------------------:|:------------------:|:-------------:|
|[![Seaice.jl Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://anders-dc.github.io/SeaIce.jl/latest) | [![Build Status](https://travis-ci.org/anders-dc/SeaIce.jl.svg?branch=master)](https://travis-ci.org/anders-dc/SeaIce.jl) | [![Build Status](https://ci.appveyor.com/api/projects/status/github/anders-dc/SeaIce.jl?svg=true)](https://ci.appveyor.com/project/anders-dc/seaice-jl/) | [![codecov.io](http://codecov.io/github/anders-dc/SeaIce.jl/coverage.svg?branch=master)](http://codecov.io/github/anders-dc/SeaIce.jl?branch=master) |

## Installation
[SeaIce.jl](https://github.com/anders-dc/SeaIce.jl) is in heavy development and 
not yet ready for third-party use.  If you know better install it directly from 
the Julia shell by:

```julia
julia> Pkg.clone("git://github.com/anders-dc/SeaIce.jl.git")
```

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/SeaIce`, and install the packages specified as 
[requirements](REQUIRE). You can run the package tests, which are contained in
the [test/ directory](test/), with the following command:

```julia
julia> Pkg.test("SeaIce")
```

The package can be updated from this repository using:

```julia
julia> Pkg.update("SeaIce")
```

For more information on installation and usage, please refer to the [documentation](https://anders-dc.github.io/SeaIce.jl/latest).

## Author
[Anders Damsgaard](mailto:anders.damsgaard@noaa.gov),
[www.adamsgaard.dk](https://adamsgaard.dk)
