#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Testing interIceFloePositionVector(...) and findOverlap(...)")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)

position_ij = SeaIce.interIceFloePositionVector(sim, 1, 2)
overlap_ij = SeaIce.findOverlap(sim, 1, 2, position_ij)

@test_approx_eq [18., 0.] position_ij
@test_approx_eq -2. overlap_ij


info("Testing findContactsAllToAll(...)")
sim_copy = deepcopy(sim)
SeaIce.findContactsAllToAll!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]


info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)

sim.ice_floes[1].fixed = true
@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

@test_throws ErrorException SeaIce.findContacts!(sim, method="")

sim = deepcopy(sim_copy)
sim.ice_floes[1].fixed = true
sim.ice_floes[2].fixed = true
SeaIce.findContacts!(sim)
@test 0 == length(sim.overlaps)
@test 0 == length(sim.contact_pairs)


info("Testing if interact(...) removes contacts correctly")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)
SeaIce.interact!(sim)
SeaIce.findContacts!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

info("Testing findContactsOceanGrid(...)")
sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
SeaIce.sortIceFloesInOceanGrid!(sim)
SeaIce.findContactsOceanGrid!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.ice_floes[1].fixed = true
SeaIce.sortIceFloesInOceanGrid!(sim)
SeaIce.findContactsOceanGrid!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.ice_floes[1].fixed = true
sim.ice_floes[2].fixed = true
SeaIce.sortIceFloesInOceanGrid!(sim)
SeaIce.findContactsOceanGrid!(sim)

@test 0 == length(sim.overlaps)
@test 0 == length(sim.contact_pairs)

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
SeaIce.sortIceFloesInOceanGrid!(sim)
SeaIce.findContacts!(sim)

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

@test 1 == length(sim.overlaps)
@test 1 == length(sim.contact_pairs)
@test_approx_eq [1, 2] sim.contact_pairs[1]
@test_approx_eq [-2., 0.] sim.overlaps[1]

@test_throws ErrorException SeaIce.findContacts!(sim, method="")
