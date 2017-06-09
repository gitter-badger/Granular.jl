#!/usr/bin/env julia

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to JLD file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)
sim.ocean = SeaIce.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
SeaIce.findContacts!(sim, method="all to all")
SeaIce.writeVTK(sim, verbose=false)

SeaIce.writeSimulation(sim)

sim2 = SeaIce.readSimulation("./test/test.1.jld")
SeaIce.compareSimulations(sim, sim2)
