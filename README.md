# SeaIce

[![Build Status](https://travis-ci.org/anders-dc/SeaIce.jl.svg?branch=master)](https://travis-ci.org/anders-dc/SeaIce.jl) [![Coverage Status](https://coveralls.io/repos/anders-dc/SeaIce.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/anders-dc/SeaIce.jl?branch=master) [![codecov.io](http://codecov.io/github/anders-dc/SeaIce.jl/coverage.svg?branch=master)](http://codecov.io/github/anders-dc/SeaIce.jl?branch=master)

Toy model for sea-ice thermodynamics and granular mechanics.

## Installation
Clone the [SeaIce.jl](https://github.com/anders-dc/SeaIce.jl) repository and add 
the source folder to your julia path, for example by:

    push!(LOAD_PATH, "/home/user/src/SeaIce.jl/src/")

If this statement is added to `~/.juliarc.jl`, it will become persistent between
julia sessions. Note that the `~` symbol for the home folder does not seem to
work (julia v. 0.4.1) in the `.juliarc.jl` file.

Import package contents with:
    import SeaIce

## Author
[Anders Damsgaard](mailto:andersd@riseup.net) 
[www.adamsgaard.dk](https://adamsgaard.dk)
