#!/usr/bin/env julia

push!(LOAD_PATH, "./src/")
import SeaIce

SeaIce.id("test")
SeaIce.addIceFloeCylindrical([0., 0., 0.], 0.5)
SeaIce.addIceFloeCylindrical([1., 0., 0.], 0.5)

SeaIce.findTimeStep()
SeaIce.setOutputFileInterval(0.1)
SeaIce.setTotalTime(1.0)
SeaIce.run()
