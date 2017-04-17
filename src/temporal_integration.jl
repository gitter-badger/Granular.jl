"""
Update the ice floe kinematic parameters using a temporal integration scheme,
the current force and torque balance, and gravitational acceleration.
"""
function updateIceFloeKinematics(i::Integer;
                                 method::String = "Three-term Taylor")

    if method == "Two-term Taylor"
        updateIceFloeKinematicsTwoTermTaylor(i)
    elseif method == "Three-term Taylor"
        updateIceFloeKinematicsThreeTermTaylor(i)
    else
        error("Unknown integration method '$method'")
    end
end

function updateIceFloeKinematicsTwoTermTaylor(i::Integer)

    # Translational components
    g_acceleration[i] = g_force[i]::vector / g_mass[i]::vector

    velocity_new::vector =
        g_velocity[i]::vector +
        g_acceleration[i]::vector * g_time_step::float

    g_position[i] = g_position[i]::vector +
        g_velocity[i]::vector * g_time_step::float +
        0.5 * g_acceleration[i]::vector * g_time_step::float^2.0

    g_velocity[i] = velocity_new

    # Rotational components
    g_rotational_acceleration[i] =
        g_torque[i]::vector / g_rotational_inertia[i]::vector

    rotational_velocity_new::vector =
        g_rotational_velocity[i]::vector +
        g_rotational_acceleration[i]::vector * g_time_step::float

    g_rotational_position[i] = g_rotational_position[i]::vector +
        g_rotational_velocity[i]::vector * g_time_step::float +
        0.5 * g_rotational_acceleration[i]::Array{float, 1} *
        g_time_step::float^2.0

    g_rotational_velocity[i] = rotational_velocity_new
end

function updateIceFloeKinematicsThreeTermTaylor(i::Integer)
    
    if g_time_iteration == 0
        acceleration_0 = zeros(3)
        rotational_acceleration_0 = zeros(3)
    else
        acceleration_0 = g_acceleration[i]::vector
        rotational_acceleration_0 = g_rotational_acceleration[i]::vector
    end

    # Temporal gradient in acceleration using backwards differences
    dacceleration_dt::vector =
            (g_acceleration[i]::vector - acceleration_0::vector) /
            g_time_step::float

    drotational_acceleration_dt::vector =
            (g_rotational_acceleration[i]::vector -
            rotational_acceleration_0::vector) /
            g_time_step::float

    # Translational components
    g_acceleration[i] = g_force[i]::vector / g_mass[i]::float

    velocity_new::vector =
        g_velocity[i]::vector +
        g_acceleration[i]::vector * g_time_step::float +
        0.5 * dacceleration_dt::vector * g_time_step::float^2.0

    g_position[i] = g_position[i]::vector +
        g_velocity[i]::vector * g_time_step::float +
        0.5 * g_acceleration[i]::vector * g_time_step::float^2.0 +
        1.0 / 6.0 * dacceleration_dt * g_time_step::float^3.0

    g_velocity[i] = velocity_new

    # Rotational components
    g_rotational_acceleration[i] =
        g_torque[i]::vector / g_rotational_inertia[i]::float

    rotational_velocity_new::vector = g_rotational_velocity[i] +
        g_rotational_acceleration[i] * g_time_step::float +
        0.5 * drotational_acceleration_dt * g_time_step::float^2.0

    g_rotational_position[i] = g_rotational_position[i]::vector +
        g_rotational_velocity[i] * g_time_step::float +
        0.5 * g_rotational_acceleration[i] * g_time_step::float^2.0 +
        1.0 / 6.0 * drotational_acceleration_dt * g_time_step::float^3.0

    g_rotational_velocity[i] = rotational_velocity_new
end

