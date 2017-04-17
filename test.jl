#!/usr/bin/env julia

push!(LOAD_PATH, "./src/")
import SeaIce

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0., 0.], 10., 1.)
SeaIce.addIceFloeCylindrical(sim, [20., 0., 0.], 10., 1.)

findTimeStep(sim)
setOutputFileInterval(sim, 0.1)
setTotalTime(sim, 1.0)
run(sim)
