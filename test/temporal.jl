info("Testing temporal functionality")

sim = Granular.createSimulation()
Test.@test_throws ErrorException Granular.setTimeStep!(sim)

Granular.setOutputFileInterval!(sim, 1e-9)
Test.@test_throws ErrorException Granular.setTotalTime!(sim, 0.)
Test.@test_throws ErrorException Granular.setCurrentTime!(sim, 0.)
Granular.setCurrentTime!(sim, 1.)
Test.@test sim.time ≈ 1.0
Test.@test_throws ErrorException Granular.incrementCurrentTime!(sim, 0.)
Granular.setOutputFileInterval!(sim, 0.)
Granular.disableOutputFiles!(sim)
Test.@test_throws ErrorException Granular.checkTimeParameters(sim)
Granular.addGrainCylindrical!(sim, [.1,.1], 2., 2.)
sim.grains[1].mass = 0.
Test.@test_throws ErrorException Granular.setTimeStep!(sim)

sim = Granular.createSimulation()
sim2 = Granular.createSimulation()
Granular.compareSimulations(sim, sim2)
