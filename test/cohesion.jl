#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical grains 

info("#### $(basename(@__FILE__)) ####")

verbose=false

sim_init = Granular.createSimulation()
Granular.addGrainCylindrical!(sim_init, [0., 0.], 10., 1.)
Granular.addGrainCylindrical!(sim_init, [18., 0.], 10., 1.)
sim_init.grains[1].youngs_modulus = 1e-5  # repulsion is negligible
sim_init.grains[2].youngs_modulus = 1e-5  # repulsion is negligible
Granular.setTimeStep!(sim_init, verbose=verbose)

info("# Check contact age scheme")
sim = deepcopy(sim_init)
Granular.setTotalTime!(sim, 10.)
sim.time_step = 1.
Granular.run!(sim, verbose=verbose)
@test sim.grains[1].contact_age[1] ≈ sim.time

info("# Check if bonds add tensile strength")
sim = Granular.createSimulation(id="cohesion")
Granular.addGrainCylindrical!(sim, [0., 0.], 10., 1., tensile_strength=500e3)
Granular.addGrainCylindrical!(sim, [20.1, 0.], 10., 1., tensile_strength=500e3)
sim.grains[1].lin_vel[1] = 0.1
Granular.setTimeStep!(sim)
Granular.setTotalTime!(sim, 10.)
Granular.run!(sim, verbose=verbose)
@test sim.grains[1].lin_vel[1] > 0.
@test sim.grains[1].lin_vel[2] ≈ 0.
@test sim.grains[2].lin_vel[1] > 0.
@test sim.grains[2].lin_vel[2] ≈ 0.
@test sim.grains[1].ang_vel ≈ 0.
@test sim.grains[2].ang_vel ≈ 0.
