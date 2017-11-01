#!/usr/bin/env julia

info("#### $(basename(@__FILE__)) ####")

info("Determining if JLD is installed")
if typeof(Pkg.installed("JLD")) == VersionNumber
    info("JLD found, proceeding with JLD-specific tests")

    info("Writing simple simulation to JLD file")
    sim = Granular.createSimulation(id="test")
    Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
    Granular.addGrainCylindrical!(sim, [18., 0.], 10., 1., verbose=false)
    sim.ocean = Granular.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
    Granular.findContacts!(sim, method="all to all")
    Granular.writeVTK(sim, verbose=false)

    Granular.writeSimulation(sim)

    sim2 = Granular.readSimulation("./test/test.1.jld")
    Granular.compareSimulations(sim, sim2)
end
