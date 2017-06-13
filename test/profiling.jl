#!/usr/bin/env julia
import Plots

info("#### $(basename(@__FILE__)) ####")

verbose=false

info("Testing performance with many interacting ice floes")

function timeSingleStepInDenseSimulation(nx::Int; verbose::Bool=true,
                                         profile::Bool=false,
                                         grid_sorting::Bool=true)

    sim = SeaIce.createSimulation()
    #nx, ny = 25, 25
    #nx, ny = 250, 250
    ny = nx
    dx, dy = 40., 40.
    sim.ocean = SeaIce.createRegularOceanGrid([nx, ny, 2], [nx*dx, ny*dy, 10.])
    if !grid_sorting
        sim.ocean.input_file = false  # fallback to all-to-all contact search
    end
    r = min(dx, dy)/2.

    # add ice floes in regular packing
    for iy=1:ny
        for ix=1:nx
            x = r + (ix - 1)*dx
            y = r + (iy - 1)*dy
            fixed = false
            if ix == 1 || iy == 1 || ix == nx || iy == ny
                fixed = true
            end
            SeaIce.addIceFloeCylindrical(sim, [x, y], r*1.1, 1.,
                                         fixed=fixed, verbose=false)
        end
    end
    info("number of ice floes: $(length(sim.ice_floes))")

    SeaIce.setTotalTime!(sim, 1.0)
    SeaIce.setTimeStep!(sim)
    SeaIce.run!(sim, single_step=true, verbose=true)
    SeaIce.run!(sim, single_step=true, verbose=true)
    SeaIce.run!(sim, single_step=true, verbose=true)
    if profile
        @profile SeaIce.run!(sim, single_step=true, verbose=true)
        if verbose
            Profile.print()
            #@time SeaIce.run!(sim, single_step=true, verbose=true)
        end
    end
    tic()
    SeaIce.run!(sim, single_step=true, verbose=true)
    t_elapsed = toc()

    #SeaIce.writeVTK(sim)

    @test sim.ice_floes[1].n_contacts == 0
    @test sim.ice_floes[2].n_contacts == 1
    @test sim.ice_floes[3].n_contacts == 1
    @test sim.ice_floes[nx].n_contacts == 0
    @test sim.ice_floes[nx + 1].n_contacts == 1
    @test sim.ice_floes[nx + 2].n_contacts == 4
    return t_elapsed
end

#nx = Int[4 8 16 32 64 128]
nx = Int[4 8 16 32 64]
t_elapsed = zeros(length(nx))
for i=1:length(nx)
    info("nx = $(nx[i])")
    t_elapsed[i] = timeSingleStepInDenseSimulation(nx[i])
end

#Plots.gr()
Plots.pyplot()
Plots.title!("Performance analysis of dense granular system")
Plots.plot(nx.*nx, t_elapsed,
           xscale=:log10,
           yscale=:log10)
Plots.xaxis!("Number of ice floes")
Plots.yaxis!("Elapsed time")
Plots.savefig("profiling.pdf")
