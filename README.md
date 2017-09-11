# SeaIce
A [Julia](https://julialang.org) package for sea-ice thermodynamics and granular 
mechanics.

| Documentation | Build Status (Linux/Mac) | Build Status (Win) | Test Coverage |
|:-------------:|:------------------------:|:------------------:|:-------------:|
|[![Seaice.jl Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://anders-dc.github.io/SeaIce.jl/latest) | [![Build Status](https://travis-ci.org/anders-dc/SeaIce.jl.svg?branch=master)](https://travis-ci.org/anders-dc/SeaIce.jl) | [![Build Status](https://ci.appveyor.com/api/projects/status/github/anders-dc/SeaIce.jl?svg=true)](https://ci.appveyor.com/project/anders-dc/seaice-jl/) | [![codecov.io](http://codecov.io/github/anders-dc/SeaIce.jl/coverage.svg?branch=master)](http://codecov.io/github/anders-dc/SeaIce.jl?branch=master) |

## Installation
[SeaIce.jl](https://github.com/anders-dc/SeaIce.jl) is in heavy development and 
not yet ready for third-party use.  If you know better install it directly from 
the Julia shell by:

    Pkg.clone("git://github.com/anders-dc/SeaIce.jl.git")

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/SeaIce`, and install the packages specified as 
[requirements](REQUIRE).  The package [JLD](https://github.com/JuliaIO/JLD.jl) 
is used for model restarts and is recommended but not required, and thus is not 
automatically installed.

Import the package contents into the current Julia session with:

    import SeaIce

This will import all functions and data types in the `SeaIce` namespace.  You 
can run the package tests, which are contained in the [test/ directory](test/), 
with the following command:

    Pkg.test("SeaIce")

The package can be updated from this repository using:

    Pkg.update("SeaIce")

## Author
[Anders Damsgaard](mailto:anders.damsgaard@noaa.gov),
[www.adamsgaard.dk](https://adamsgaard.dk)
