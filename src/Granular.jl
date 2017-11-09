#!/usr/bin/env julia
if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
else
    using Test
end

"""
# Granular.jl
Offline granular dynamics simulator module.
"""
module Granular

include("datatypes.jl")
include("grain.jl")
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
