#!/usr/bin/env julia

# Check the basic dynamic wall functionality

info("#### $(basename(@__FILE__)) ####")

sim = Granular.createSimulation(id="test")
info("Testing grain value checks ")
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim,
                                                                [.1, .1, .1],
                                                                1.)
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim,
                                                                [1., 1.],
                                                                1.)


info("Check that wall mass equals total grain mass and max. thickness")
sim = Granular.createSimulation(id="test")
@test length(sim.walls) == 0
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 2., verbose=false)
sim.grains[1].mass = 1.0
Granular.addWallLinearFrictionless!(sim, [1., 0.], 1., verbose=true)
@test length(sim.walls) == 1
@test sim.walls.mass ≈ 1.0
@test sim.walls.mass ≈ 2.0

