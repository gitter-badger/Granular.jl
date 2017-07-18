#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Testing interIceFloePositionVector(...) and findOverlap(...)")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [18., 0.], 10., 1., verbose=false)

position_ij = SeaIce.interIceFloePositionVector(sim, 1, 2)
overlap_ij = SeaIce.findOverlap(sim, 1, 2, position_ij)

@test [-18., 0.] ≈ position_ij
@test -2. ≈ overlap_ij


info("Testing findContactsAllToAll(...)")
sim_copy = deepcopy(sim)
SeaIce.findContactsAllToAll!(sim)


info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)

sim.ice_floes[1].fixed = true
# The contact should be registered in ice floe 1, but not ice floe 2
@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)

@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts

@test_throws ErrorException SeaIce.findContacts!(sim, method="")

sim = deepcopy(sim_copy)
sim.ice_floes[1].fixed = true
sim.ice_floes[2].fixed = true
SeaIce.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 0 == sim.ice_floes[1].n_contacts
@test 0 == sim.ice_floes[2].n_contacts


sim = deepcopy(sim_copy)
SeaIce.disableIceFloe!(sim, 1)
SeaIce.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 0 == sim.ice_floes[1].n_contacts
@test 0 == sim.ice_floes[2].n_contacts


sim = deepcopy(sim_copy)
SeaIce.disableIceFloe!(sim, 1)
SeaIce.disableIceFloe!(sim, 2)
SeaIce.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 0 == sim.ice_floes[1].n_contacts
@test 0 == sim.ice_floes[2].n_contacts

info("Testing if interact(...) removes contacts correctly")
sim = deepcopy(sim_copy)
SeaIce.findContacts!(sim)
SeaIce.interact!(sim)
SeaIce.findContacts!(sim)

@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts


info("Testing findContactsGrid(...)")
sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
SeaIce.sortIceFloesInGrid!(sim, sim.ocean)
SeaIce.findContactsInGrid!(sim, sim.ocean)

@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts


sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.ice_floes[1].fixed = true
SeaIce.sortIceFloesInGrid!(sim, sim.ocean)
SeaIce.findContactsInGrid!(sim, sim.ocean)

@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts


sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.ice_floes[1].fixed = true
sim.ice_floes[2].fixed = true
SeaIce.sortIceFloesInGrid!(sim, sim.ocean)
SeaIce.findContactsInGrid!(sim, sim.ocean)

for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 0 == sim.ice_floes[1].n_contacts
@test 0 == sim.ice_floes[2].n_contacts

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
SeaIce.sortIceFloesInGrid!(sim, sim.ocean)
SeaIce.findContacts!(sim)

@test 2 == sim.ice_floes[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.ice_floes[1].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.ice_floes[2].contacts[ic]
    @test [0., 0.] ≈ sim.ice_floes[2].contact_parallel_displacement[ic]
end
@test 1 == sim.ice_floes[1].n_contacts
@test 1 == sim.ice_floes[2].n_contacts

@test_throws ErrorException SeaIce.findContacts!(sim, method="")

info("Testing contact registration with multiple contacts")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [2., 2.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [4., 2.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [6., 2.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [2., 4.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [4., 4.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [6., 4.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [2., 6.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [4., 6.], 1.01, 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [6., 6.], 1.01, 1., verbose=false)
sim.ocean = SeaIce.createRegularOceanGrid([4, 4, 2], [8., 8., 2.])
SeaIce.sortIceFloesInGrid!(sim, sim.ocean)
SeaIce.findContacts!(sim)
@test 2 == sim.ice_floes[1].n_contacts
@test 3 == sim.ice_floes[2].n_contacts
@test 2 == sim.ice_floes[3].n_contacts
@test 3 == sim.ice_floes[4].n_contacts
@test 4 == sim.ice_floes[5].n_contacts
@test 3 == sim.ice_floes[6].n_contacts
@test 2 == sim.ice_floes[7].n_contacts
@test 3 == sim.ice_floes[8].n_contacts
@test 2 == sim.ice_floes[9].n_contacts
SeaIce.interact!(sim)
SeaIce.interact!(sim)
SeaIce.interact!(sim)
SeaIce.interact!(sim)
@test 2 == sim.ice_floes[1].n_contacts
@test 3 == sim.ice_floes[2].n_contacts
@test 2 == sim.ice_floes[3].n_contacts
@test 3 == sim.ice_floes[4].n_contacts
@test 4 == sim.ice_floes[5].n_contacts
@test 3 == sim.ice_floes[6].n_contacts
@test 2 == sim.ice_floes[7].n_contacts
@test 3 == sim.ice_floes[8].n_contacts
@test 2 == sim.ice_floes[9].n_contacts
for i=1:9
    sim.ice_floes[i].contact_radius = 0.99
end
SeaIce.interact!(sim)
for i=1:9
    @test sim.ice_floes[i].n_contacts == 0
end
