## Manage icefloes in the model

# Use immutable composite type for efficiency
immutable IceFloeCylindrical

    # Material properties
    density::float

    # Geometrical parameters
    thickness::float
    contact_radius::float
    areal_radius::float
    surface_area::float
    volume::float
    mass::float
    moment_of_inertia::float

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
    fixed::Bool
    rotating::Bool

    # Rheological parameters
    contact_stiffness_normal::float
    contact_stiffness_tangential::float
    contact_viscosity_normal::float
    contact_viscosity_tangential::float
    contact_dynamic_friction::float

end
    

"""
Adds a grain to the simulation. Example:

    SeaIce.addIceFloeCylindrical([1.0, 2.0, 3.0], 1.0)
"""
function addIceFloeCylindrical(
    lin_pos::vector,
    contact_radius::float,
    thickness::float;
    areal_radius = false,
    lin_vel::vector = [0., 0.],
    lin_acc::vector = [0., 0.],
    force::vector = [0., 0.],
    ang_pos::float = 0.,
    ang_vel::float = 0.,
    ang_acc::float = 0.,
    torque::float = 0.,
    density::float = 934.,
    contact_stiffness_normal::float = 1.e6,
    contact_stiffness_tangential::float = 1.e6,
    contact_viscosity_normal::float = 0.,
    contact_viscosity_tangential::float = 0.,
    fixed::Bool = false,
    rotating::Bool = true,
    verbose::Bool = true)

    # Check input values
    if contact_radius <= 0.0
        error("Radius must be greater than 0.0 (radius = $contact_radius m)")
    elseif density <= 0.0
        error("Density must be greater than 0.0 (density = $density
            kg/m^3)")
    end

    if !areal_radius
        areal_radius = contact_radius
    end


    # Create icefloe object with placeholder values for surface area, volume, 
    # mass, and moment of inertia.
    icefloe = IceFloeCylindrical(density=density,
                                 thickness=thickness,
                                 contact_radius=contact_radius,
                                 areal_radius=areal_radius,
                                 surface_area=1.0,
                                 volume=1.0,
                                 mass=1.0,
                                 moment_of_inertia=1.0,

                                 lin_pos=lin_pos,
                                 lin_vel=lin_vel,
                                 lin_acc=lin_acc,
                                 force=force,

                                 ang_pos=ang_pos,
                                 ang_vel=ang_vel,
                                 ang_acc=ang_acc,
                                 torque=torque,

                                 fixed=fixed,
                                 rotating=rotating,

                                 contact_stiffness_normal=
                                     contact_stiffness_normal,
                                 contact_stiffness_tangential=
                                     contact_stiffness_tangential,
                                 contact_viscosity_normal=
                                     contact_viscosity_normal,
                                 contact_viscosity_tangential=
                                     contact_viscosity_tangential
                                )

    # Overwrite previous placeholder values
    icefloe.surface_area = iceFloeSurfaceArea(icefloe)
    icefloe.volume = iceFloeVolume(icefloe)
    icefloe.mass = iceFloeMass(icefloe)
    icefloe.moment_of_inertia = iceFloeMomentOfInertia(icefloe)
    
    # Append icefloe to global icefloe array
    push!(g_ice_floes, icefloe)

    if verbose
        info("Added IceFloe $(length(g_ice_floes))")
    end

end

function removeIceFloe(i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    delete!(g_ice_floes, i)
end

function iceFloeSurfaceArea(icefloe::IceFloeCylindrical)
    return pi*icefloe.areal_radius^2.
end

function iceFloeVolume(icefloe::IceFloeCylindrical)
    return iceFloeSurfaceArea(icefloe)*icefloe.thickness
end

function iceFloeMass(icefloe::IceFloeCylindrical)
    return iceFloeVolume(icefloe)*icefloe.density
end

function iceFloeMomentOfInertia(icefloe::IceFloeCylindrical)
    return 0.5*iceFloeMass(icefloe)*icefloe.radius^2.
end
