#!/usr/bin/env julia
import SeaIce
using Base.Test

info("#### $(basename(@__FILE__)) ####")

verbose=false

info("Testing memory footprint of SeaIce types")

sim = SeaIce.createSimulation()
empty_sim_size = 96
empty_sim_size_recursive = 472

@test sizeof(sim) == empty_sim_size
@test Base.summarysize(sim) == empty_sim_size_recursive

size_per_icefloe = 352
size_per_icefloe_recursive = 1136

info("Testing memory usage when adding ice floes")
for i=1:100
    SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1., 1., verbose=false)

    @test sizeof(sim) == empty_sim_size

    @test sizeof(sim.ice_floes) == sizeof(Int)*i
    @test sizeof(sim.ice_floes[:]) == sizeof(Int)*i
    @test Base.summarysize(sim.ice_floes) == size_per_icefloe_recursive*i + 
        sizeof(Int)*i

    @test Base.summarysize(sim) == empty_sim_size_recursive + sizeof(Int)*i + 
        size_per_icefloe_recursive*i

    @test Base.summarysize(sim.ice_floes[i]) == size_per_icefloe_recursive

    for j=1:i
        @test sizeof(sim.ice_floes[j]) == size_per_icefloe
        @test Base.summarysize(sim.ice_floes[j]) == size_per_icefloe_recursive
    end

end

info("Checking memory footprint when overwriting simulation object")
sim = SeaIce.createSimulation()
@test sizeof(sim) == empty_sim_size
@test Base.summarysize(sim) == empty_sim_size_recursive

info("Check memory usage when stepping time for empty simulation object")
sim = SeaIce.createSimulation()
sim.time_step = 1.0
for i=1:10
    SeaIce.run!(sim, single_step=true, verbose=false)
    @test sizeof(sim) == empty_sim_size
    @test Base.summarysize(sim) == empty_sim_size_recursive
end

info("Check memory when stepping time with single ice floe")
sim = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1., 1., verbose=false)
sim.time_step = 1.0
for i=1:10
    SeaIce.run!(sim, single_step=true, verbose=false)
    @test sizeof(sim) == empty_sim_size
    @test Base.summarysize(sim) == empty_sim_size_recursive + 
        sizeof(Int)*length(sim.ice_floes) + 
        size_per_icefloe_recursive*length(sim.ice_floes)
end

info("Check memory when stepping time with two separate ice floes")
sim = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 3., 1., verbose=false)
sim.time_step = 1.0
for i=1:10
    SeaIce.run!(sim, single_step=true, verbose=false)
    @test sizeof(sim) == empty_sim_size
    @test Base.summarysize(sim) == empty_sim_size_recursive + 
        sizeof(Int)*length(sim.ice_floes) + 
        size_per_icefloe_recursive*length(sim.ice_floes)
end

info("Check memory when stepping time with two interacting ice floes (all to all)")
sim = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1.9, 1., verbose=false)
sim.time_step = 1.0
for i=1:10
    SeaIce.run!(sim, single_step=true, verbose=false)
    @test sizeof(sim) == empty_sim_size
    @test Base.summarysize(sim) == empty_sim_size_recursive + 
        sizeof(Int)*length(sim.ice_floes) + 
        size_per_icefloe_recursive*length(sim.ice_floes)
end

info("Check memory when stepping time with two interacting ice floes (cell sorting)")
sim = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [1., 1.], 1.9, 1., verbose=false)
nx, ny, nz = 5, 5, 1
sim.ocean = SeaIce.createRegularOceanGrid([nx, ny, nz], [10., 10., 10.])
sim.time_step = 1e-6
SeaIce.run!(sim, single_step=true, verbose=false)
original_size_recursive = Base.summarysize(sim)
for i=1:10
    SeaIce.run!(sim, single_step=true, verbose=false)
    @test Base.summarysize(sim) == original_size_recursive
end

info("Checking if memory is freed after ended collision (all to all)")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1
SeaIce.setTotalTime!(sim, 10.0)
SeaIce.setTimeStep!(sim, epsilon=0.07)
original_size_recursive = Base.summarysize(sim)
SeaIce.run!(sim, verbose=false)
@test Base.summarysize(sim) == original_size_recursive

info("Checking if memory is freed after ended collision (cell sorting)")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.ocean = SeaIce.createRegularOceanGrid([2, 2, 2], [40., 40., 10.])
sim.ice_floes[1].lin_vel[1] = 0.1
SeaIce.setTotalTime!(sim, 10.0)
SeaIce.setTimeStep!(sim, epsilon=0.07)
SeaIce.run!(sim, single_step=true, verbose=false)
original_sim_size_recursive = Base.summarysize(sim)
original_icefloes_size_recursive = Base.summarysize(sim.ice_floes)
original_ocean_size_recursive = Base.summarysize(sim.ocean)
original_atmosphere_size_recursive = Base.summarysize(sim.atmosphere)
SeaIce.run!(sim, verbose=false)
@test Base.summarysize(sim.ice_floes) == original_icefloes_size_recursive
@test Base.summarysize(sim.ocean) == original_ocean_size_recursive
@test Base.summarysize(sim.atmosphere) == original_atmosphere_size_recursive
@test Base.summarysize(sim) == original_sim_size_recursive
