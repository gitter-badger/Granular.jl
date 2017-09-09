# Installation
SeaIce.jl can be installed directly from the Julia shell by:

    Pkg.clone("git://github.com/anders-dc/SeaIce.jl.git")

This will install the contents of this repository in the folder 
`~/.julia/v$(JULIA_VERSION)/SeaIce`, and install the packages specified as 
[requirements](REQUIRE).  The package [JLD](https://github.com/JuliaIO/JLD.jl) 
is used for model restarts and is recommended but not required, and thus is not 
automatically installed.

Import the package contents into the current Julia session or script with:

    import SeaIce

This will import all functions and data types in the `SeaIce` namespace.  You 
can run the package tests, which are contained in the [test/ directory](test/), 
with the following command:

    Pkg.test("SeaIce")

The package can be updated from this repository using:

    Pkg.update("SeaIce")
