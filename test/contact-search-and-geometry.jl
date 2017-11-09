#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Testing interGrainPositionVector(...) and findOverlap(...)")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.addGrainCylindrical!(sim, [18., 0.], 10., 1., verbose=false)

position_ij = Granular.interGrainPositionVector(sim, 1, 2)
overlap_ij = Granular.findOverlap(sim, 1, 2, position_ij)

@test [-18., 0.] ≈ position_ij
@test -2. ≈ overlap_ij


info("Testing findContactsAllToAll(...)")
sim_copy = deepcopy(sim)
Granular.findContactsAllToAll!(sim)


info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
Granular.findContacts!(sim)

sim.grains[1].fixed = true
# The contact should be registered in ice floe 1, but not ice floe 2
@test 2 == sim.grains[1].contacts[1]
@test [-18., 0.] ≈ sim.grains[1].position_vector[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
Granular.findContacts!(sim)

@test 2 == sim.grains[1].contacts[1]
@test [-18., 0.] ≈ sim.grains[1].position_vector[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts

@test_throws ErrorException Granular.findContacts!(sim, method="")

sim = deepcopy(sim_copy)
sim.grains[1].fixed = true
sim.grains[2].fixed = true
Granular.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 0 == sim.grains[1].n_contacts
@test 0 == sim.grains[2].n_contacts


sim = deepcopy(sim_copy)
Granular.disableGrain!(sim, 1)
Granular.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 0 == sim.grains[1].n_contacts
@test 0 == sim.grains[2].n_contacts


sim = deepcopy(sim_copy)
Granular.disableGrain!(sim, 1)
Granular.disableGrain!(sim, 2)
Granular.findContacts!(sim)
for ic=1:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 0 == sim.grains[1].n_contacts
@test 0 == sim.grains[2].n_contacts

info("Testing if interact(...) removes contacts correctly")
sim = deepcopy(sim_copy)
Granular.findContacts!(sim)
Granular.interact!(sim)
Granular.findContacts!(sim)

@test 2 == sim.grains[1].contacts[1]
@test [-18., 0.] ≈ sim.grains[1].position_vector[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts


info("Testing findContactsGrid(...)")
sim = deepcopy(sim_copy)
sim.ocean = Granular.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
Granular.sortGrainsInGrid!(sim, sim.ocean)
Granular.findContactsInGrid!(sim, sim.ocean)

@test 2 == sim.grains[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts


sim = deepcopy(sim_copy)
sim.ocean = Granular.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.grains[1].fixed = true
Granular.sortGrainsInGrid!(sim, sim.ocean)
Granular.findContactsInGrid!(sim, sim.ocean)

@test 2 == sim.grains[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts


sim = deepcopy(sim_copy)
sim.ocean = Granular.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
sim.grains[1].fixed = true
sim.grains[2].fixed = true
Granular.sortGrainsInGrid!(sim, sim.ocean)
Granular.findContactsInGrid!(sim, sim.ocean)

for ic=1:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 0 == sim.grains[1].n_contacts
@test 0 == sim.grains[2].n_contacts

info("Testing findContacts(...)")
sim = deepcopy(sim_copy)
sim.ocean = Granular.createRegularOceanGrid([4, 4, 2], [80., 80., 2.])
Granular.sortGrainsInGrid!(sim, sim.ocean)
Granular.findContacts!(sim)

@test 2 == sim.grains[1].contacts[1]
for ic=2:sim.Nc_max
    @test 0 == sim.grains[1].contacts[ic]
    @test [0., 0.] ≈ sim.grains[1].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[1].contact_parallel_displacement[ic]
end
for ic=1:sim.Nc_max
    @test 0 == sim.grains[2].contacts[ic]
    @test [0., 0.] ≈ sim.grains[2].position_vector[ic]
    @test [0., 0.] ≈ sim.grains[2].contact_parallel_displacement[ic]
end
@test 1 == sim.grains[1].n_contacts
@test 1 == sim.grains[2].n_contacts

@test_throws ErrorException Granular.findContacts!(sim, method="")

info("Testing contact registration with multiple contacts")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [2., 2.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [4., 2.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [6., 2.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [2., 4.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [4., 4.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [6., 4.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [2., 6.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [4., 6.], 1.01, 1., verbose=false)
Granular.addGrainCylindrical!(sim, [6., 6.], 1.01, 1., verbose=false)
sim.ocean = Granular.createRegularOceanGrid([4, 4, 2], [8., 8., 2.])
Granular.sortGrainsInGrid!(sim, sim.ocean)
Granular.findContacts!(sim)
@test 2 == sim.grains[1].n_contacts
@test 3 == sim.grains[2].n_contacts
@test 2 == sim.grains[3].n_contacts
@test 3 == sim.grains[4].n_contacts
@test 4 == sim.grains[5].n_contacts
@test 3 == sim.grains[6].n_contacts
@test 2 == sim.grains[7].n_contacts
@test 3 == sim.grains[8].n_contacts
@test 2 == sim.grains[9].n_contacts
Granular.interact!(sim)
Granular.interact!(sim)
Granular.interact!(sim)
Granular.interact!(sim)
@test 2 == sim.grains[1].n_contacts
@test 3 == sim.grains[2].n_contacts
@test 2 == sim.grains[3].n_contacts
@test 3 == sim.grains[4].n_contacts
@test 4 == sim.grains[5].n_contacts
@test 3 == sim.grains[6].n_contacts
@test 2 == sim.grains[7].n_contacts
@test 3 == sim.grains[8].n_contacts
@test 2 == sim.grains[9].n_contacts
for i=1:9
    sim.grains[i].contact_radius = 0.99
end
Granular.interact!(sim)
for i=1:9
    @test sim.grains[i].n_contacts == 0
end
