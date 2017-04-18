"""
Update the ice floe kinematic parameters using a temporal integration scheme,
the current force and torque balance, and gravitational acceleration.
"""
function updateIceFloeKinematics!(simulation::Simulation;
                                  method::String = "Three-term Taylor")

    if method == "Two-term Taylor"
        for ice_floe in simulation.ice_floes
            updateIceFloeKinematicsTwoTermTaylor(ice_floe, simulation)
        end
    elseif method == "Three-term Taylor"
        for ice_floe in simulation.ice_floes
            updateIceFloeKinematicsThreeTermTaylor(ice_floe, simulation)
        end
    else
        error("Unknown integration method '$method'")
    end
end

function updateIceFloeKinematicsTwoTermTaylor(icefloe::IceFloeCylindrical,
                                              simulation::Simulation)
    icefloe.lin_acc = icefloe.force/icefloe.mass
    icefloe.ang_acc = icefloe.torque/icefloe.moment_of_inertia

    icefloe.lin_pos +=
        icefloe.lin_vel * simulation.time_step +
        0.5*icefloe.lin_acc * simulation.time_step^2.0
    icefloe.ang_pos +=
        icefloe.ang_vel * simulation.time_step +
        0.5*icefloe.ang_acc * simulation.time_step^2.0

    icefloe.lin_vel += icefloe.lin_acc * simulation.time_step
    icefloe.ang_vel += icefloe.ang_acc * simulation.time_step
end

function updateIceFloeKinematicsThreeTermTaylor(icefloe::IceFloeCylindrical,
                                                simulation::Simulation)
    
    if simulation.time_iteration == 0
        lin_acc_0 = zeros(3)
        ang_acc_0 = zeros(3)
    else
        lin_acc_0 = icefloe.lin_acc
        ang_acc_0 = icefloe.ang_acc
    end

    icefloe.lin_acc = icefloe.force/icefloe.mass
    icefloe.ang_acc = icefloe.torque/icefloe.moment_of_inertia

    # Temporal gradient in acceleration using backwards differences
    d_lin_acc_dt::vector = (icefloe.lin_acc - lin_acc_0)/simulation.time_step
    d_ang_acc_dt::vector = (icefloe.ang_acc - ang_acc_0)/simulation.time_step

    icefloe.lin_pos +=
        icefloe.lin_vel * simulation.time_step +
        0.5*icefloe.lin_acc * simulation.time_step^2. +
        1./6. * d_lin_acc_dt * simulation.time_step^3.
    icefloe.ang_pos +=
        icefloe.ang_vel * simulation.time_step +
        0.5 * icefloe.ang_acc * simulation.time_step^2. +
        1./6. * d_ang_acc_dt * simulation.time_step^3.

    icefloe.lin_vel += icefloe.lin_acc * simulation.time_step +
        0.5 * d_lin_acc_dt * simulation.time_step^2.
    icefloe.ang_vel += icefloe.ang_acc * simulation.time_step +
        0.5 * d_ang_acc_dt * simulation.time_step^2.
end

