#!/usr/bin/env julia

# Check the basic icefloe functionality

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.printGrainInfo(sim.grains[1])

info("Testing grain value checks ")
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1, .1],
                                                          10., 1.)
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_vel=[.2,.2,.2])
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_acc=[.2,.2,.2])
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          0., 1.)
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., density=-2.)
@test_throws ErrorException Granular.disableGrain!(sim, 0)

info("Testing grain comparison ")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.compareGrains(sim.grains[1], sim.grains[2])

if typeof(Pkg.installed("PyPlot")) == VersionNumber
    info("Testing GSD plotting ")
    Granular.plotGrainSizeDistribution(sim)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, skip_fixed=false)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, log_y=false)
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, size_type="areal")
    rm("test-grain-size-distribution.png")
    @test_throws ErrorException Granular.plotGrainSizeDistribution(sim, size_type="asdf")
else
    @test_throws ErrorException Granular.plotGrainSizeDistribution(sim)
end

info("Testing external body force routines")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.setBodyForce!(sim.grains[1], [1., 2.])
@test sim.grains[1].external_body_force ≈ [1., 2.]
Granular.addBodyForce!(sim.grains[1], [1., 2.])
@test sim.grains[1].external_body_force ≈ [2., 4.]

info("Testing zeroKinematics!()")
sim.grains[1].force .= ones(2)
sim.grains[1].lin_acc .= ones(2)
sim.grains[1].lin_vel .= ones(2)
sim.grains[1].torque = 1.
sim.grains[1].ang_acc = 1.
sim.grains[1].ang_vel = 1.
Granular.zeroKinematics!(sim)
@test Granular.totalGrainKineticTranslationalEnergy(sim) ≈ 0.
@test Granular.totalGrainKineticRotationalEnergy(sim) ≈ 0.
