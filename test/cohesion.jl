#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical ice floes 

info("#### $(basename(@__FILE__)) ####")

verbose=false

sim_init = SeaIce.createSimulation()
SeaIce.addIceFloeCylindrical!(sim_init, [0., 0.], 10., 1.)
SeaIce.addIceFloeCylindrical!(sim_init, [18., 0.], 10., 1.)
sim_init.ice_floes[1].youngs_modulus = 1e-5  # repulsion is negligible
sim_init.ice_floes[2].youngs_modulus = 1e-5  # repulsion is negligible
SeaIce.setTimeStep!(sim_init, verbose=verbose)

info("# Check contact age scheme")
sim = deepcopy(sim_init)
SeaIce.setTotalTime!(sim, 10.)
sim.time_step = 1.
SeaIce.run!(sim, verbose=verbose)
@test sim.ice_floes[1].contact_age[1] ≈ sim.time

info("# Check if bonds add tensile strength")
sim = SeaIce.createSimulation(id="cohesion")
SeaIce.addIceFloeCylindrical!(sim, [0., 0.], 10., 1., tensile_strength=500e3)
SeaIce.addIceFloeCylindrical!(sim, [20.1, 0.], 10., 1., tensile_strength=500e3)
sim.ice_floes[1].lin_vel[1] = 0.1
SeaIce.setTimeStep!(sim)
SeaIce.setTotalTime!(sim, 10.)
SeaIce.run!(sim, verbose=verbose)
@test sim.ice_floes[1].lin_vel[1] > 0.
@test sim.ice_floes[1].lin_vel[2] ≈ 0.
@test sim.ice_floes[2].lin_vel[1] > 0.
@test sim.ice_floes[2].lin_vel[2] ≈ 0.
@test sim.ice_floes[1].ang_vel ≈ 0.
@test sim.ice_floes[2].ang_vel ≈ 0.
