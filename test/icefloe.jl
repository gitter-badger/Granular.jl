#!/usr/bin/env julia

# Check the basic icefloe functionality

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.printIceFloeInfo(sim.ice_floes[1])


@test_throws ErrorException SeaIce.addIceFloeCylindrical!(sim, [.1, .1, .1],
                                                          10., 1.)
@test_throws ErrorException SeaIce.addIceFloeCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_vel=[.2,.2,.2])
@test_throws ErrorException SeaIce.addIceFloeCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_acc=[.2,.2,.2])
@test_throws ErrorException SeaIce.addIceFloeCylindrical!(sim, [.1, .1],
                                                          0., 1.)
@test_throws ErrorException SeaIce.addIceFloeCylindrical!(sim, [.1, .1],
                                                          10., 1., density=-2.)
@test_throws ErrorException SeaIce.disableIceFloe!(sim, 0)

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.compareIceFloes(sim.ice_floes[1], sim.ice_floes[2])

if typeof(Pkg.installed("PyPlot")) == VersionNumber
    SeaIce.plotIceFloeSizeDistribution(sim)
    rm("test-ice-floe-size-distribution.png")
    SeaIce.plotIceFloeSizeDistribution(sim, skip_fixed=false)
    rm("test-ice-floe-size-distribution.png")
    SeaIce.plotIceFloeSizeDistribution(sim, log_y=false)
    rm("test-ice-floe-size-distribution.png")
    SeaIce.plotIceFloeSizeDistribution(sim, size_type="areal")
    rm("test-ice-floe-size-distribution.png")
    @test_throws ErrorException SeaIce.plotIceFloeSizeDistribution(sim, size_type="asdf")
else
    @test_throws ErrorException SeaIce.plotIceFloeSizeDistribution(sim)
end
