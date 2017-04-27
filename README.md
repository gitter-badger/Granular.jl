# SeaIce

[![Build Status](https://travis-ci.org/anders-dc/SeaIce.jl.svg?branch=master)](https://travis-ci.org/anders-dc/SeaIce.jl) [![codecov.io](http://codecov.io/github/anders-dc/SeaIce.jl/coverage.svg?branch=master)](http://codecov.io/github/anders-dc/SeaIce.jl?branch=master)

A [Julia](https://julialang.org) package for sea-ice thermodynamics and granular 
mechanics.

## Installation
[SeaIce.jl](https://github.com/anders-dc/SeaIce.jl) is in heavy development and 
not yet ready for third-party use.  If you know better install it directly from 
the Julia shell by:

    Pkg.clone("git://github.com/anders-dc/SeaIce.jl.git")

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/SeaIce`, and install the packages specified as 
[requirements](REQUIRE).  Import the package contents into the current Julia 
session with:

    import SeaIce

This will import all functions and data types in the `SeaIce` namespace.  You 
can run the package tests, which are contained in the [test/ directory](test/), 
with the following command:

    Pkg.test("SeaIce")

The package can be updated from this repository using:

    Pkg.update("SeaIce")

## Documentation
All functions and types are documented via docstrings.  The documentation can be 
displayed in the Julia shell by typing `?` followed by the function or type 
name.

## Author
[Anders Damsgaard](mailto:anders.damsgaard@noaa.gov),
[www.adamsgaard.dk](https://adamsgaard.dk)
