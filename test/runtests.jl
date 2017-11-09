using Compat
if VERSION < v"0.7.0-DEV.2004"
    using Base.Test
else
    using Test
end
import Granular

include("grain.jl")
include("packing.jl")
include("util.jl")
include("temporal.jl")
include("contact-search-and-geometry.jl")
include("collision-2floes-normal.jl")
include("collision-5floes-normal.jl")
include("collision-2floes-oblique.jl")
include("cohesion.jl")
include("netcdf.jl")
include("vtk.jl")
include("jld.jl")
include("grid.jl")
include("grid-boundaries.jl")
include("ocean.jl")
include("atmosphere.jl")
include("memory-management.jl")
