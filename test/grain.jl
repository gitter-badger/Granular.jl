#!/usr/bin/env julia

# Check the basic icefloe functionality

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.printGrainInfo(sim.grains[1])


Test.@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1, .1],
                                                          10., 1.)
Test.@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_vel=[.2,.2,.2])
Test.@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_acc=[.2,.2,.2])
Test.@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          0., 1.)
Test.@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., density=-2.)
Test.@test_throws ErrorException Granular.disableGrain!(sim, 0)

sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.compareGrains(sim.grains[1], sim.grains[2])

if typeof(Pkg.installed("PyPlot")) == VersionNumber
    Granular.plotGrainSizeDistribution(sim)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, skip_fixed=false)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, log_y=false)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, size_type="areal")
    rm("test-grain-size-distribution.png")
    Test.@test_throws ErrorException Granular.plotGrainSizeDistribution(sim, size_type="asdf")
else
    Test.@test_throws ErrorException Granular.plotGrainSizeDistribution(sim)
end

sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.setBodyForce!(sim.grains[1], [1., 2.])
Test.@test sim.grains[1].external_body_force ≈ [1., 2.]
Granular.addBodyForce!(sim.grains[1], [1., 2.])
Test.@test sim.grains[1].external_body_force ≈ [2., 4.]
