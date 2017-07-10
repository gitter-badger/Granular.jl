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

