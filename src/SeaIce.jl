#!/usr/bin/env julia

##############################
# Sea-ice dynamics simulator #
##############################

#=
Add SeaIce to your julia path, e.g. by:
    push!(LOAD_PATH, "/home/user/src/seaice/")

If this statement is added to `~/.juliarc.jl`, it will become persistent between
julia sessions. Note that the `~` symbol for the home folder does not seem to
work (julia v. 0.4.1) in the `.juliarc.jl` file.

Import package contents with:
    import SeaIce

If this file is changed, reimport it with:
    reload("SeaIce")
=#

module SeaIce

include("datatypes.jl")
include("icefloe.jl")
include("simulation.jl")
include("arrays.jl")
include("grid.jl")
include("packing.jl")
include("contact_search.jl")
include("interaction.jl")
include("temporal.jl")
include("temporal_integration.jl")
include("io.jl")

end # module end
