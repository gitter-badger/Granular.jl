#!/usr/bin/env julia

# Check the basic dynamic wall functionality

info("#### $(basename(@__FILE__)) ####")

info("Testing argument value checks")
sim = Granular.createSimulation()
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 2., verbose=false)
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim,
                                                                [.1, .1, .1],
                                                                1.)
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim,
                                                                [1., 1.],
                                                                1.)
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim,
                                                                [.1, .1, .1],
                                                                1.)
sim = Granular.createSimulation()
@test_throws ErrorException Granular.addWallLinearFrictionless!(sim, [1., 0.],
                                                                1.)


info("Check that wall mass equals total grain mass and max. thickness")
sim = Granular.createSimulation()
@test length(sim.walls) == 0
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 2., verbose=false)
sim.grains[1].mass = 1.0
Granular.addWallLinearFrictionless!(sim, [1., 0.], 1., verbose=true)
@test length(sim.walls) == 1
@test sim.walls[1].mass ≈ 1.0
@test sim.walls[1].thickness ≈ 2.0

info("Test wall surface area and defined normal stress")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 2., verbose=false)
sim.grains[1].mass = 1.0
Granular.addWallLinearFrictionless!(sim, [1., 0.], 1., verbose=false)
Granular.addWallLinearFrictionless!(sim, [0., 1.], 1., verbose=false)
@test length(sim.walls) == 2
@test sim.walls[1].mass ≈ 1.0
@test sim.walls[1].thickness ≈ 2.0
@test sim.walls[2].mass ≈ 1.0
@test sim.walls[2].thickness ≈ 2.0
@test Granular.getWallSurfaceArea(sim, 1) ≈ 20.0*2.0
@test Granular.getWallSurfaceArea(sim, 2) ≈ 10.0*2.0

info("Test wall-grain interaction")

info("...wall present but no contact")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [1., 0.], -1.01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force ≈ 0.
@test sim.grains[1].force[1] ≈ 0.
@test sim.grains[1].force[2] ≈ 0.

info("...wall present but no contact")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [-1., 0.], +2.01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force ≈ 0.
@test sim.grains[1].force[1] ≈ 0.
@test sim.grains[1].force[2] ≈ 0.

info("...wall at -x")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [1., 0.], -1. + .01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force > 0.
@test sim.grains[1].force[1] > 0.
@test sim.grains[1].force[2] ≈ 0.

info("...wall at +x")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [-1., 0.], +1. - .01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force > 0.
@test sim.grains[1].force[1] < 0.
@test sim.grains[1].force[2] ≈ 0.

info("...wall at -y")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [0., 1.], -1. + .01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force > 0.
@test sim.grains[1].force[1] ≈ 0.
@test sim.grains[1].force[2] > 0.

info("...wall at +y")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [10., 20., 1.0])
Granular.addGrainCylindrical!(sim, [ 0., 0.], 1., 2., verbose=false)
Granular.addWallLinearFrictionless!(sim, [0., -1.], +1. - .01, verbose=false)
Granular.interactWalls!(sim)
@test sim.walls[1].force > 0.
@test sim.grains[1].force[1] ≈ 0.
@test sim.grains[1].force[2] < 0.
