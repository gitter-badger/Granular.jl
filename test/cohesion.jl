#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical ice floes 

info("#### $(basename(@__FILE__)) ####")

verbose=false

sim_init = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical(sim_init, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim_init, [18., 0.], 10., 1., verbose=verbose)
sim_init.ice_floes[1].youngs_modulus = 1e-5  # repulsion is negligible
sim_init.ice_floes[2].youngs_modulus = 1e-5  # repulsion is negligible
SeaIce.setTimeStep!(sim_init, verbose=verbose)

info("# Check contact age scheme")
sim = deepcopy(sim_init)
SeaIce.setTotalTime!(sim, 10.)
sim.time_step = 1.
SeaIce.run!(sim, verbose=verbose)
@test sim.ice_floes[1].contact_age[1] ≈ sim.time

info("# Check if cohesion increases with time")
sim = SeaIce.createSimulation(id="cohesion")
SeaIce.addIceFloeCylindrical(sim, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.5, 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].youngs_modulus = 1e-5  # repulsion is negligible
sim.ice_floes[2].youngs_modulus = 1e-5  # repulsion is negligible
SeaIce.setTimeStep!(sim)
SeaIce.setTotalTime!(sim, 24.*60.*60.)
sim.file_time_step = 60.
# let the contact age for a while
while sim.time_total*.9 > sim.time
    SeaIce.run!(sim, single_step=true, verbose=verbose)
end
#sim_init.ice_floes[1].youngs_modulus = 2e7
#sim_init.ice_floes[2].youngs_modulus = 2e7
#SeaIce.setTimeStep!(sim)
