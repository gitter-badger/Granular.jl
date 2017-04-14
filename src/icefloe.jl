## Manage icefloes in the model

# Use immutable composite type for efficiency
immutable IceFloe

    # Material properties
    density::float

    # Geometrical parameters
    thickness::float
    contact_radius::float
    areal_radius::float

    # Linear kinematic degrees of freedom along horizontal plane
    lin_pos::vector
    lin_vel::vector
    lin_acc::vector
    force::vector

    # Angular kinematic degrees of freedom for vertical rotation around center
    ang_pos::float
    ang_vel::float
    ang_acc::float
    torque::float

    # Kinematic constraint flags
    fixed::bool
    rotating::bool
end
    

"""
Adds a grain to the simulation. Example:

    SeaIce.addIceFloeCylindrical([1.0, 2.0, 3.0], 1.0)
"""
function addIceFloeCylindrical(
    position::vector,
    radius::float = 1.0;
    velocity::vector = [0., 0., 0.],
    acceleration::vector = [0., 0., 0.],
    force::vector = [0., 0., 0.],
    rotational_position::vector = [0., 0., 0.],
    rotational_velocity::vector = [0., 0., 0.],
    rotational_acceleration::vector = [0., 0., 0.],
    torque::vector = [0., 0., 0.],
    density::float = 934.,
    verbose::Bool = true)

    # Check input values
    if radius <= 0.0
        error("Radius must be greater than 0.0 (radius = $radius m)")
    elseif density <= 0.0
        error("Density must be greater than 0.0 (density = $density
            kg/m^3)")
    end

    # Save grain in global arrays
    push!(g_radius, radius)

    push!(g_position, position)
    push!(g_velocity, velocity)
    push!(g_acceleration, acceleration)
    push!(g_force, force)

    push!(g_rotational_position, rotational_position)
    push!(g_rotational_velocity, rotational_velocity)
    push!(g_rotational_acceleration, rotational_acceleration)
    push!(g_torque, torque)

    push!(g_density, density)
    volume = 4.0 / 3.0 * pi * radius^3.0
    push!(g_volume, volume)
    push!(g_mass, density*volume)
    push!(g_rotational_inertia, 2.0 / 5.0 * density * volume * radius^2.0)

    if verbose
        info("Added IceFloe $(length(g_radius))")
    end

    min_position = position - radius
    max_position = position + radius

    # Update world boundaries
    for i = 1:3
        if min_position[i] < g_origo[i]
            g_origo[i] = min_position[i]
        end

        if max_position[i] > g_world_size[i]
            g_world_size[i] = max_position[i]
        end
    end
end

function removeIceFloe(i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    delete!(g_radius, i)

    delete!(g_position, i)
    delete!(g_velocity, i)
    delete!(g_acceleration, i)
    delete!(g_force, i)

    delete!(g_rotational_position, i)
    delete!(g_rotational_velocity, i)
    delete!(g_rotational_acceleration, i)
    delete!(g_torque, i)

    delete!(g_density, i)
    delete!(g_volume, i)
    delete!(g_mass, i)
    delete!(g_rotational_inertia, i)
end
