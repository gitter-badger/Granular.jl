## Manage icefloes in the model

"""
Adds a grain to the simulation. Example:

    SeaIce.addIceFloeCylindrical([1.0, 2.0, 3.0], 1.0)
"""
function addIceFloeCylindrical(simulation::Simulation,
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
                               contact_static_friction::float = 0.4,
                               contact_dynamic_friction::float = 0.4,
                               fixed::Bool = false,
                               rotating::Bool = true,
                               verbose::Bool = true)

    # Check input values
    if length(lin_pos) != 2
        error("Linear position must be a two-element array (lin_pos = ",
              "$lin_pos)")
    end
    if length(lin_vel) != 2
        error("Linear velocity must be a two-element array (lin_vel = ",
              "$lin_vel)")
    end
    if length(lin_acc) != 2
        error("Linear acceleration must be a two-element array (lin_acc = ",
              "$lin_acc)")
    end
    if length(ang_pos) != 1
        error("Angular position must be a scalar (ang_pos = $ang_pos)")
    end
    if length(ang_vel) != 1
        error("Angular velocity must be a scalar (ang_vel = $ang_vel)")
    end
    if length(ang_acc) != 1
        error("Angular acceleration must be a scalar (ang_acc = $ang_acc)")
    end
    if contact_radius <= 0.0
        error("Radius must be greater than 0.0 (radius = $contact_radius m)")
    end
    if density <= 0.0
        error("Density must be greater than 0.0 (density = $density kg/m^3)")
    end

    if !areal_radius
        areal_radius = contact_radius
    end


    # Create icefloe object with placeholder values for surface area, volume, 
    # mass, and moment of inertia.
    icefloe = IceFloeCylindrical(density,

                                 thickness,
                                 contact_radius,
                                 areal_radius,
                                 1.0,  # surface_area
                                 1.0,  # volume
                                 1.0,  # mass
                                 1.0,  # moment_of_inertia
                                 lin_pos,
                                 lin_vel,
                                 lin_acc,
                                 force,

                                 ang_pos,
                                 ang_vel,
                                 ang_acc,
                                 torque,

                                 fixed,
                                 rotating,

                                 contact_stiffness_normal,
                                 contact_stiffness_tangential,
                                 contact_viscosity_normal,
                                 contact_viscosity_tangential,
                                 contact_static_friction,
                                 contact_dynamic_friction
                                )

    # Overwrite previous placeholder values
    icefloe.surface_area = iceFloeSurfaceArea(icefloe)
    icefloe.volume = iceFloeVolume(icefloe)
    icefloe.mass = iceFloeMass(icefloe)
    icefloe.moment_of_inertia = iceFloeMomentOfInertia(icefloe)

    # Add to simulation object
    addIceFloe!(simulation, icefloe, verbose)
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
    return 0.5*iceFloeMass(icefloe)*icefloe.areal_radius^2.
end
