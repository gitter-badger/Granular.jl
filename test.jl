#!/usr/bin/env julia

push!(LOAD_PATH, "./src/")
import SeaIce

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1.)
SeaIce.addIceFloeCylindrical(sim, [20., 0.], 10., 1.)
sim.ice_floes[1].lin_vel[1] = 1.

SeaIce.setTimeStep!(sim)
SeaIce.setOutputFileInterval!(sim, 0.1)
#SeaIce.setTotalTime!(sim, 10.0)
SeaIce.setTotalTime!(sim, 1.0)
SeaIce.run!(sim)
