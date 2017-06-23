info("Testing temporal functionality")

sim = SeaIce.createSimulation()
@test_throws ErrorException SeaIce.setTimeStep!(sim)

SeaIce.setOutputFileInterval!(sim, 1e-9)
@test_throws ErrorException SeaIce.setTotalTime!(sim, 0.)
@test_throws ErrorException SeaIce.setCurrentTime!(sim, 0.)
SeaIce.setCurrentTime!(sim, 1.)
@test sim.time ≈ 1.0
@test_throws ErrorException SeaIce.incrementCurrentTime!(sim, 0.)
SeaIce.setOutputFileInterval!(sim, 0.)
SeaIce.disableOutputFiles!(sim)
@test_throws ErrorException SeaIce.checkTimeParameters(sim)
SeaIce.addIceFloeCylindrical!(sim, [.1,.1], 2., 2.)
sim.ice_floes[1].mass = 0.
@test_throws ErrorException SeaIce.setTimeStep!(sim)

sim = SeaIce.createSimulation()
sim2 = SeaIce.createSimulation()
SeaIce.compareSimulations(sim, sim2)
