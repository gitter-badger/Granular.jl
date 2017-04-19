#!/usr/bin/env julia

# Check for conservation of kinetic energy (=momentum) during a normal collision 
# between two ice cylindrical ice floes 

push!(LOAD_PATH, "../src/")

import Base.Test
import SeaIce

sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1.)
SeaIce.addIceFloeCylindrical(sim, [20., 0.], 10., 1.)
sim.ice_floes[1].lin_vel[1] = 1.

E_kin_lin_init = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_init = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

SeaIce.setTimeStep!(sim)
SeaIce.setTotalTime!(sim, 1.0)
SeaIce.run!(sim)

E_kin_lin_final = SeaIce.totalIceFloeKineticTranslationalEnergy(sim)
E_kin_rot_final = SeaIce.totalIceFloeKineticRotationalEnergy(sim)

Base.Test.@test_approx_eq E_kin_lin_init E_kin_lin_final
Base.Test.@test_approx_eq E_kin_rot_init E_kin_rot_final
