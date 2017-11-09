#!/usr/bin/env julia

verbose = true

info("#### $(basename(@__FILE__)) ####")

info("Testing regular packing generation (power law GSD)")
sim = Granular.createSimulation()
Granular.regularPacking!(sim, [2, 2], 1., 1., size_distribution="powerlaw")
@test 4 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius ≈ 1.
end

sim = Granular.createSimulation()
Granular.regularPacking!(sim, [10, 10], 1., 10., size_distribution="powerlaw")
@test 100 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius >= 1.
    @test grain.contact_radius <= 10.
end

info("Testing regular packing generation (uniform GSD)")
sim = Granular.createSimulation()
Granular.regularPacking!(sim, [2, 2], 1., 1., size_distribution="uniform")
@test 4 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius ≈ 1.
end

sim = Granular.createSimulation()
Granular.regularPacking!(sim, [10, 10], 1., 10., size_distribution="uniform")
@test 100 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius >= 1.
    @test grain.contact_radius <= 10.
end
