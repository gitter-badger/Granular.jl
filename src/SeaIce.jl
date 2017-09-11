#!/usr/bin/env julia

"""
# SeaIce.jl
Offline sea-ice dynamics simulator module.
"""
module SeaIce

include("datatypes.jl")
include("icefloe.jl")
include("simulation.jl")
include("grid.jl")
include("packing.jl")
include("contact_search.jl")
include("interaction.jl")
include("temporal.jl")
include("temporal_integration.jl")
include("io.jl")
include("ocean.jl")
include("atmosphere.jl")
include("util.jl")

end # module end
