#!/usr/bin/env julia
import Plots
import SeaIce
import CurveFit
using Base.Test

info("#### $(basename(@__FILE__)) ####")

verbose=false

info("Testing performance with many interacting ice floes")

function timeSingleStepInDenseSimulation(nx::Int; verbose::Bool=true,
                                         profile::Bool=false,
                                         grid_sorting::Bool=true,
                                         include_atmosphere::Bool=false)

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
    if include_atmosphere
        sim.atmosphere = SeaIce.createRegularAtmosphereGrid([nx, ny, 2],
                                                            [nx*dx, ny*dy, 10.])
    end

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
    if grid_sorting
        info("using cell-based spatial decomposition (ocean)")
        if include_atmosphere
            info("using cell-based spatial decomposition (ocean + atmosphere)")
        end
    else
        info("using all-to-all contact search")
    end

    SeaIce.setTotalTime!(sim, 1.0)
    SeaIce.setTimeStep!(sim)
    SeaIce.run!(sim, single_step=true, verbose=true)
    if profile
        @profile SeaIce.run!(sim, single_step=true, verbose=true)
        if verbose
            Profile.print()
        end
        SeaIce.run!(sim, single_step=true, verbose=true)
    end
    n_runs = 4
    t_elapsed = 1e12
    for i=1:n_runs
        tic()
        @time SeaIce.run!(sim, single_step=true, verbose=true)
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

#nx = Int[4 8 16 32 64 96]
nx = round(logspace(1, 2, 16))
elements = zeros(length(nx))
t_elapsed = zeros(length(nx))
t_elapsed_all_to_all = zeros(length(nx))
t_elapsed_cell_sorting = zeros(length(nx))
t_elapsed_cell_sorting2 = zeros(length(nx))
for i=1:length(nx)
    info("nx = $(nx[i])")
    t_elapsed_all_to_all[i] =
        timeSingleStepInDenseSimulation(Int(nx[i]), grid_sorting=false)
    t_elapsed_cell_sorting[i] =
        timeSingleStepInDenseSimulation(Int(nx[i]), grid_sorting=true)
    t_elapsed_cell_sorting2[i] =
        timeSingleStepInDenseSimulation(Int(nx[i]), grid_sorting=true, 
                                        include_atmosphere=true)
    elements[i] = nx[i]*nx[i]
end

#Plots.gr()
Plots.pyplot()
Plots.scatter(elements, t_elapsed_all_to_all,
              xscale=:log10,
              yscale=:log10,
              label="All to all")
fit_all_to_all = CurveFit.curve_fit(CurveFit.PowerFit,
                                    elements, t_elapsed_all_to_all)
label_all_to_all = @sprintf "%1.3g n^%3.2f" fit_all_to_all.coefs[1] fit_all_to_all.coefs[2]
Plots.plot!(elements, fit_all_to_all(elements),
            xscale=:log10,
            yscale=:log10,
            label=label_all_to_all)

Plots.scatter!(elements, t_elapsed_cell_sorting,
               xscale=:log10,
               yscale=:log10,
               label="Cell-based spatial decomposition (ocean only)")
fit_cell_sorting = CurveFit.curve_fit(CurveFit.PowerFit,
                                    elements, t_elapsed_cell_sorting)
label_cell_sorting = @sprintf "%1.3g n^%3.2f" fit_cell_sorting.coefs[1] fit_cell_sorting.coefs[2]
Plots.plot!(elements, fit_cell_sorting(elements),
            xscale=:log10,
            yscale=:log10,
            label=label_cell_sorting)

Plots.scatter!(elements, t_elapsed_cell_sorting2,
               xscale=:log10,
               yscale=:log10,
               label="Cell-based spatial decomposition (ocean + atmosphere)")
fit_cell_sorting2 = CurveFit.curve_fit(CurveFit.PowerFit,
                                       elements, t_elapsed_cell_sorting2)
label_cell_sorting2 = @sprintf "%1.3g n^%3.2f" fit_cell_sorting2.coefs[1] fit_cell_sorting2.coefs[2]
Plots.plot!(elements, fit_cell_sorting2(elements),
            xscale=:log10,
            yscale=:log10,
            label=label_cell_sorting2)

Plots.title!("Dense granular system " * "(host: $(gethostname()))")
Plots.xaxis!("Number of ice floes")
Plots.yaxis!("Wall time per time step [s]")
Plots.savefig("profiling.pdf")
