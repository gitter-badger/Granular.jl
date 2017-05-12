#!/usr/bin/env julia

# Check the basic icefloe functionality

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.printIceFloeInfo(sim.ice_floes[1])


@test_throws ErrorException SeaIce.addIceFloeCylindrical(sim, [.1, .1, .1], 10., 1.)
@test_throws ErrorException SeaIce.addIceFloeCylindrical(sim, [.1, .1], 10., 1., 
    lin_vel=[.2,.2,.2])
@test_throws ErrorException SeaIce.addIceFloeCylindrical(sim, [.1, .1], 10., 1., 
    lin_acc=[.2,.2,.2])
@test_throws ErrorException SeaIce.addIceFloeCylindrical(sim, [.1, .1], 0., 1.)
@test_throws ErrorException SeaIce.addIceFloeCylindrical(sim, [.1, .1], 10., 1., 
density=-2.)
