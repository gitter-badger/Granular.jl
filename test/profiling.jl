#!/usr/bin/env julia
import Plots
import SeaIce
using Base.Test

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
    if profile
        @profile SeaIce.run!(sim, single_step=true, verbose=true)
        if verbose
            Profile.print()
        end
    end
    n_runs = 4
    t_elapsed = 1e12
    for i=1:n_runs
        tic()
        SeaIce.run!(sim, single_step=true, verbose=true)
        t = toc()
        if t < t_elapsed
            t_elapsed = t
        end
    end

    #SeaIce.writeVTK(sim)

    @test sim.ice_floes[1].n_contacts == 0
    @test sim.ice_floes[2].n_contacts == 1
    @test sim.ice_floes[3].n_contacts == 1
    @test sim.ice_floes[nx].n_contacts == 0
    @test sim.ice_floes[nx + 1].n_contacts == 1
    @test sim.ice_floes[nx + 2].n_contacts == 4
    return t_elapsed
end

#nx = Int[4 8 10 12 16 19 24 28 32 36 42 50 64]
nx = round(logspace(1, 2, 16))
elements = zeros(length(nx))
t_elapsed = zeros(length(nx))
t_elapsed_all_to_all = zeros(length(nx))
t_elapsed_cell_sorting = zeros(length(nx))
for i=1:length(nx)
    info("nx = $(nx[i])")
    t_elapsed_all_to_all[i] =
        timeSingleStepInDenseSimulation(Int(nx[i]), grid_sorting=false)
    t_elapsed_cell_sorting[i] =
        timeSingleStepInDenseSimulation(Int(nx[i]), grid_sorting=true)
    elements[i] = nx[i]*nx[i]
end

#Plots.gr()
Plots.pyplot()
Plots.scatter(elements, t_elapsed_all_to_all,
              xscale=:log10,
              yscale=:log10,
              label="All to all")
Plots.scatter!(elements, t_elapsed_cell_sorting,
               xscale=:log10,
               yscale=:log10,
               label="Cell-based spatial decomposition")
Plots.title!("Dense granular system " * "(host: $(gethostname())")
Plots.xaxis!("Number of ice floes")
Plots.yaxis!("Wall time per time step [s]")
Plots.savefig("profiling.pdf")
