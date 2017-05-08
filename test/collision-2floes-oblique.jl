#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical ice floes 

info("#### $(basename(@__FILE__)) ####")

verbose=false

info("## Contact-normal elasticity only")
info("# One ice floe fixed")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 10.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19., 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1
sim.ice_floes[2].fixed = true

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
#sim.file_time_step = 1.
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor", 
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final

info("# Ice floes free to move")

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 10.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.0, 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor", 
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.02
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor", 
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.01
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="None", verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("## Contact-normal elasticity and tangential viscosity and friction")
info("Testing gamma_t = 0.")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 10.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.0, 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos ≈ 0.
@test sim.ice_floes[1].ang_vel ≈ 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final

info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.02
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos ≈ 0.
@test sim.ice_floes[1].ang_vel ≈ 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final


info("## Contact-normal elasticity and tangential viscosity and friction")
info("# One ice floe fixed")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 10.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19., 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1
sim.ice_floes[2].fixed = true
sim.ice_floes[1].contact_viscosity_tangential = 1e4
sim.ice_floes[2].contact_viscosity_tangential = 1e4

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
#sim.file_time_step = 1.
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional",
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
@test sim.ice_floes[2].ang_vel ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol

info("mu_d = 0.")
sim = deepcopy(sim_init)
sim.ice_floes[1].contact_dynamic_friction = 0.
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.02
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos ≈ 0.
@test sim.ice_floes[1].ang_vel ≈ 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init E_kin_lin_final E_kin_lin_init*tol
@test_approx_eq E_kin_rot_init E_kin_rot_final

info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.06
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
@test sim.ice_floes[2].ang_vel ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.07
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos ≈ 0.
@test sim.ice_floes[2].ang_vel ≈ 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol

info("# Ice floes free to move")

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 10.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.0, 0.], 10., 1., verbose=verbose)
sim.ice_floes[1].lin_vel[1] = 0.1
sim.ice_floes[1].contact_viscosity_tangential = 1e4
sim.ice_floes[2].contact_viscosity_tangential = 1e4

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos > 0.
@test sim.ice_floes[2].ang_vel > 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 

info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.05
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.04
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos > 0.
@test sim.ice_floes[2].ang_vel > 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 


info("# Ice floes free to move, mirrored")

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.0, 10.], 10., 1., verbose=verbose)
sim.ice_floes[2].lin_vel[1] = -0.1
sim.ice_floes[1].contact_viscosity_tangential = 1e4
sim.ice_floes[2].contact_viscosity_tangential = 1e4

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos > 0.
@test sim.ice_floes[2].ang_vel > 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 

info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.05
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.04
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos > 0.
@test sim.ice_floes[1].ang_vel > 0.
@test sim.ice_floes[2].ang_pos > 0.
@test sim.ice_floes[2].ang_vel > 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 


info("# Ice floes free to move, mirrored #2")

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [0., 0.], 10., 1., verbose=verbose)
SeaIce.addIceFloeCylindrical(sim, [19.0, -10.], 10., 1., verbose=verbose)
sim.ice_floes[2].lin_vel[1] = -0.1
sim.ice_floes[1].contact_viscosity_tangential = 1e4
sim.ice_floes[2].contact_viscosity_tangential = 1e4

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

# With decreasing timestep (epsilon towards 0), the explicit integration scheme 
# should become more correct

SeaIce.setTotalTime!(sim, 30.0)
sim_init = deepcopy(sim)

info("Testing kinetic energy conservation with Two-term Taylor scheme")
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.2
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos < 0.
@test sim.ice_floes[1].ang_vel < 0.
@test sim.ice_floes[2].ang_pos < 0.
@test sim.ice_floes[2].ang_vel < 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 

info("Testing kinetic energy conservation with Two-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.007)
tol = 0.05
info("Relative tolerance: $(tol*100.)%")
SeaIce.run!(sim, temporal_integration_method="Two-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 


info("Testing kinetic energy conservation with Three-term Taylor scheme")
sim = deepcopy(sim_init)
SeaIce.setTimeStep!(sim, epsilon=0.07)
tol = 0.04
info("Relative tolerance: $(tol*100.)% with time step: $(sim.time_step)")
SeaIce.run!(sim, temporal_integration_method="Three-term Taylor",
            contact_tangential_rheology="Linear Viscous Frictional", 
            verbose=verbose)

@test sim.ice_floes[1].ang_pos < 0.
@test sim.ice_floes[1].ang_vel < 0.
@test sim.ice_floes[2].ang_pos < 0.
@test sim.ice_floes[2].ang_vel < 0.
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test_approx_eq_eps E_kin_lin_init+E_kin_rot_init E_kin_lin_final+E_kin_rot_final E_kin_lin_init*tol 
