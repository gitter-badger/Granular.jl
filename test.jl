#!/usr/bin/env julia

push!(LOAD_PATH, "./seaice/")
import SeaIce

SeaIce.id("test")
SeaIce.addCylindricalIceFloe([0., 0., 0.], 0.5)
SeaIce.addCylindricalIceFloe([1., 0., 0.], 0.5)

SeaIce.findTimeStep()
SeaIce.setOutputFileInterval(0.1)
SeaIce.setTotalTime(1.0)
SeaIce.run()
