#!/usr/bin/env julia

# Check if atmosphere-specific functions and grid operations are functioning 
# correctly

info("#### $(basename(@__FILE__)) ####")

info("Testing regular grid generation")
sim = SeaIce.createSimulation()
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([6, 6, 6], [1., 1., 1.])
@test size(sim.atmosphere.xq) == (7, 7)
@test size(sim.atmosphere.yq) == (7, 7)
@test size(sim.atmosphere.xh) == (6, 6)
@test size(sim.atmosphere.yh) == (6, 6)
@test sim.atmosphere.xq[1, :, 1] ≈ zeros(7)
@test sim.atmosphere.xq[4, :, 1] ≈ .5*ones(7)
@test sim.atmosphere.xq[end, :, 1] ≈ 1.*ones(7)
@test sim.atmosphere.yq[:, 1, 1] ≈ zeros(7)
@test sim.atmosphere.yq[:, 4, 1] ≈ .5*ones(7)
@test sim.atmosphere.yq[:, end, 1] ≈ 1.*ones(7)
@test size(sim.atmosphere.u) == (7, 7, 6, 1)
@test size(sim.atmosphere.v) == (7, 7, 6, 1)
@test sim.atmosphere.u ≈ zeros(7, 7, 6, 1)
@test sim.atmosphere.v ≈ zeros(7, 7, 6, 1)

info("Testing velocity drag interaction (static atmosphere)")
SeaIce.addIceFloeCylindrical!(sim, [.5, .5], .25, .1)
SeaIce.setTotalTime!(sim, 5.)
SeaIce.setTimeStep!(sim)
sim_init = deepcopy(sim)
sim.ice_floes[1].lin_vel[1] = 0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test E_kin_rot_init ≈ E_kin_rot_final  # no rotation before or after
@test E_kin_lin_init > E_kin_lin_final  # linear velocity lost due to atmos drag
@test sim.ice_floes[1].atmosphere_stress[1] < 0.
@test sim.ice_floes[1].atmosphere_stress[2] ≈ 0.

info("Testing velocity drag interaction (static ice floe)")
sim = deepcopy(sim_init)
sim.atmosphere.v[:, :, 1, 1] = 0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test E_kin_rot_init ≈ E_kin_rot_final  # no rotation before or after
@test E_kin_lin_init < E_kin_lin_final  # linear vel. gained due to atmos drag
@test sim.ice_floes[1].atmosphere_stress[1] ≈ 0.
@test sim.ice_floes[1].atmosphere_stress[2] > 0.

info("Testing vortex interaction (static atmosphere)")
sim = deepcopy(sim_init)
sim.ice_floes[1].ang_vel = 0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test E_kin_rot_init > E_kin_rot_final  # energy lost to atmosphere
@test sim.ice_floes[1].ang_vel > 0.     # check angular velocity orientation
@test sim.ice_floes[1].ang_pos > 0.     # check angular position orientation
@test E_kin_lin_init ≈ E_kin_lin_final  # no linear velocity gained

info("Testing vortex interaction (static ice floe)")
sim = deepcopy(sim_init)
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([1, 1, 1], [1., 1., 1.])
sim.ice_floes[1].lin_pos[1] = 0.5
sim.ice_floes[1].lin_pos[2] = 0.5
sim.atmosphere.v[1, :, 1, 1] = -0.1
sim.atmosphere.v[2, :, 1, 1] = 0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test sim.ice_floes[1].ang_vel > 0.     # check angular velocity orientation
@test sim.ice_floes[1].ang_pos > 0.     # check angular position orientation
@test E_kin_rot_init < E_kin_rot_final  # rotation after due to atm vortex
@test E_kin_lin_init ≈ E_kin_lin_final  # no linear velocity gained

sim = deepcopy(sim_init)
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([1, 1, 1], [1., 1., 1.])
sim.ice_floes[1].lin_pos[1] = 0.5
sim.ice_floes[1].lin_pos[2] = 0.5
sim.atmosphere.v[1, :, 1, 1] = 0.1
sim.atmosphere.v[2, :, 1, 1] = -0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test sim.ice_floes[1].ang_vel < 0.     # check angular velocity orientation
@test sim.ice_floes[1].ang_pos < 0.     # check angular position orientation
@test E_kin_rot_init < E_kin_rot_final  # rotation after due to atm vortex
@test E_kin_lin_init ≈ E_kin_lin_final  # no linear velocity gained

sim = deepcopy(sim_init)
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([1, 1, 1], [1., 1., 1.])
sim.ice_floes[1].lin_pos[1] = 0.5
sim.ice_floes[1].lin_pos[2] = 0.5
sim.atmosphere.u[:, 1, 1, 1] = -0.1
sim.atmosphere.u[:, 2, 1, 1] = 0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test sim.ice_floes[1].ang_vel < 0.     # check angular velocity orientation
@test sim.ice_floes[1].ang_pos < 0.     # check angular position orientation
@test E_kin_rot_init < E_kin_rot_final  # rotation after due to atm vortex
@test E_kin_lin_init ≈ E_kin_lin_final  # no linear velocity gained

sim = deepcopy(sim_init)
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([1, 1, 1], [1., 1., 1.])
sim.ice_floes[1].lin_pos[1] = 0.5
sim.ice_floes[1].lin_pos[2] = 0.5
sim.atmosphere.u[:, 1, 1, 1] = 0.1
sim.atmosphere.u[:, 2, 1, 1] = -0.1
E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
SeaIce.run!(sim, verbose=false)
E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)
@test sim.ice_floes[1].ang_vel > 0.     # check angular velocity orientation
@test sim.ice_floes[1].ang_pos > 0.     # check angular position orientation
@test E_kin_rot_init < E_kin_rot_final  # rotation after due to atm vortex
@test E_kin_lin_init ≈ E_kin_lin_final  # no linear velocity gained

sim = SeaIce.createSimulation()
sim.atmosphere = SeaIce.createRegularAtmosphereGrid([6, 6, 6], [1., 1., 1.])
sim2 = SeaIce.createSimulation()
sim2.atmosphere = SeaIce.createRegularAtmosphereGrid([6, 6, 6], [1., 1., 1.])
SeaIce.compareAtmospheres(sim.atmosphere, sim2.atmosphere)
