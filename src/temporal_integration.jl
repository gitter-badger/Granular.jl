export updateIceFloeKinematics!
"""
    updateIceFloeKinematics!(simulation::Simulation[,
                             method::String = "Three-term Taylor"])

Update the ice floe kinematic parameters using a temporal integration scheme,
the current force and torque balance, and gravitational acceleration.

# Arguments
* `simulation::Simulation`: update the ice floe positions in this object 
    according to temporal integration of length `simulation.time_step`.
* `method::String = "Three-term Taylor"`: the integration method to use.  
    Available methods are "Two-term Taylor" and "Three-term Taylor".  The 
    three-term Taylor expansion is slightly more computationally expensive than 
    the two-term Taylor expansion, but offers an order-of-magnitude increase in 
    precision of ice floe positions.  The two-term expansion can obtain similar 
    precision if the time step is 1/10 the length.
"""
function updateIceFloeKinematics!(simulation::Simulation;
                                  method::String = "Three-term Taylor")

    if method == "Two-term Taylor"
        for ice_floe in simulation.ice_floes
            if !ice_floe.enabled
                continue
            end
            updateIceFloeKinematicsTwoTermTaylor!(ice_floe, simulation)
        end
    elseif method == "Three-term Taylor"
        for ice_floe in simulation.ice_floes
            if !ice_floe.enabled
                continue
            end
            updateIceFloeKinematicsThreeTermTaylor!(ice_floe, simulation)
        end
    else
        error("Unknown integration method '$method'")
    end
end

export updateIceFloeKinematicsTwoTermTaylor!
"""
Use a two-term Taylor expansion for integrating the kinematic degrees of freedom 
for an `icefloe`.
"""
function updateIceFloeKinematicsTwoTermTaylor!(icefloe::IceFloeCylindrical,
                                               simulation::Simulation)
    icefloe.lin_acc = icefloe.force/icefloe.mass
    icefloe.ang_acc = icefloe.torque/icefloe.moment_of_inertia

    if icefloe.fixed
        fill!(icefloe.lin_acc, 0.)
        icefloe.ang_acc = 0.
    elseif !icefloe.rotating
        icefloe.ang_acc = 0.
    end

    icefloe.lin_pos +=
        icefloe.lin_vel * simulation.time_step +
        0.5*icefloe.lin_acc * simulation.time_step^2.0
    icefloe.ang_pos +=
        icefloe.ang_vel * simulation.time_step +
        0.5*icefloe.ang_acc * simulation.time_step^2.0

    icefloe.lin_vel += icefloe.lin_acc * simulation.time_step
    icefloe.ang_vel += icefloe.ang_acc * simulation.time_step
end

export updateIceFloeKinematicsThreeTermTaylor!
"""
Use a three-term Taylor expansion for integrating the kinematic degrees of 
freedom for an `icefloe`.
"""
function updateIceFloeKinematicsThreeTermTaylor!(icefloe::IceFloeCylindrical,
                                                 simulation::Simulation)

    if simulation.time_iteration == 0
        lin_acc_0 = zeros(2)
        ang_acc_0 = 0.
    else
        lin_acc_0 = icefloe.lin_acc
        ang_acc_0 = icefloe.ang_acc
    end

    icefloe.lin_acc = icefloe.force/icefloe.mass
    icefloe.ang_acc = icefloe.torque/icefloe.moment_of_inertia

    if icefloe.fixed
        fill!(icefloe.lin_acc, 0.)
        icefloe.ang_acc = 0.
    elseif !icefloe.rotating
        icefloe.ang_acc = 0.
    end

    # Temporal gradient in acceleration using backwards differences
    d_lin_acc_dt = (icefloe.lin_acc - lin_acc_0)/simulation.time_step
    d_ang_acc_dt = (icefloe.ang_acc - ang_acc_0)/simulation.time_step

    icefloe.lin_pos +=
        icefloe.lin_vel * simulation.time_step +
        0.5 * icefloe.lin_acc * simulation.time_step^2. +
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
