#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction
import Base.Test
import SeaIce

info("#### $(basename(@__FILE__)) ####")

info("Testing interIceFloePositionVector(...) and findOverlap(...)")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)

position_ij = SeaIce.interIceFloePositionVector(sim, 1, 2)
overlap_ij = SeaIce.findOverlap(sim, 1, 2, position_ij)

Base.Test.@test_approx_eq [18., 0.] position_ij
Base.Test.@test_approx_eq -2. overlap_ij


info("Testing findContactsAllToAll(...)")
sim_copy = deepcopy(sim)
SeaIce.findContactsAllToAll!(sim)

Base.Test.@test 1 == length(sim.overlaps)
Base.Test.@test 1 == length(sim.contact_pairs)
Base.Test.@test_approx_eq [1, 2] sim.contact_pairs[1]
Base.Test.@test_approx_eq [-2., 0.] sim.overlaps[1]
