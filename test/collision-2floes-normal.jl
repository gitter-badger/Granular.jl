#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical grains 

info("#### $(basename(@__FILE__)) ####")

verbose=false

info("# One ice floe fixed")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
Granular.addGrainCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.grains[1].lin_vel[1] = 0.1
sim.grains[2].fixed = true

E_kin_lin_init = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_init = Granular.totalGrainKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

Granular.setTotalTime!(sim, 10.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Three-term Taylor",
            verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("# Ice floes free to move")

sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
Granular.addGrainCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.grains[1].lin_vel[1] = 0.1

E_kin_lin_init = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_init = Granular.totalGrainKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

Granular.setTotalTime!(sim, 10.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Three-term Taylor",
            verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init ≈ E_kin_lin_final atol=E_kin_lin_init*tol
@test E_kin_rot_init ≈ E_kin_rot_final


info("# Adding contact-normal viscosity")
info("# One ice floe fixed")
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
Granular.addGrainCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.grains[1].lin_vel[1] = 0.1
sim.grains[1].contact_viscosity_normal = 1e4
sim.grains[2].contact_viscosity_normal = 1e4
sim.grains[2].fixed = true

E_kin_lin_init = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_init = Granular.totalGrainKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

Granular.setTotalTime!(sim, 10.0)
sim_init = deepcopy(sim)


info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init > E_kin_lin_final
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Three-term Taylor",
            verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init > E_kin_lin_final
@test E_kin_rot_init ≈ E_kin_rot_final


info("# Ice floes free to move")

sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [0., 0.], 10., 1., verbose=verbose)
Granular.addGrainCylindrical!(sim, [20.05, 0.], 10., 1., verbose=verbose)
sim.grains[1].lin_vel[1] = 0.1
sim.grains[1].contact_viscosity_normal = 1e4
sim.grains[2].contact_viscosity_normal = 1e4

E_kin_lin_init = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_init = Granular.totalGrainKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

Granular.setTotalTime!(sim, 10.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
Granular.run!(sim, temporal_integration_method="Two-term Taylor", verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init > E_kin_lin_final
@test E_kin_rot_init ≈ E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
Granular.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
Granular.run!(sim, temporal_integration_method="Three-term Taylor",
            verbose=verbose)

E_kin_lin_final = Granular.totalGrainKineticTranslationalEnergy(sim)
E_kin_rot_final = Granular.totalGrainKineticRotationalEnergy(sim)
@test E_kin_lin_init > E_kin_lin_final
@test E_kin_rot_init ≈ E_kin_rot_final
