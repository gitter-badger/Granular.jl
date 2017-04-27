## Manage icefloes in the model

export addIceFloeCylindrical
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
                               pressure::float = 0.,
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
                                 1.0,  # circumreference
                                 1.0,  # horizontal_surface_area
                                 1.0,  # side_surface_area
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
                                 contact_dynamic_friction,

                                 pressure
                                )

    # Overwrite previous placeholder values
    icefloe.circumreference = iceFloeCircumreference(icefloe)
    icefloe.horizontal_surface_area = iceFloeHorizontalSurfaceArea(icefloe)
    icefloe.side_surface_area = iceFloeSideSurfaceArea(icefloe)
    icefloe.volume = iceFloeVolume(icefloe)
    icefloe.mass = iceFloeMass(icefloe)
    icefloe.moment_of_inertia = iceFloeMomentOfInertia(icefloe)

    # Add to simulation object
    addIceFloe!(simulation, icefloe, verbose)
end

export iceFloeCircumreference
"Returns the circumreference of the ice floe"
function iceFloeCircumreference(icefloe::IceFloeCylindrical)
    return pi*icefloe.areal_radius*2.
end

export iceFloeHorizontalSurfaceArea
"Returns the top or bottom (horizontal) surface area of the ice floe"
function iceFloeHorizontalSurfaceArea(icefloe::IceFloeCylindrical)
    return pi*icefloe.areal_radius^2.
end

export iceFloeSideSurfaceArea
"Returns the surface area of the ice-floe sides"
function iceFloeSideSurfaceArea(icefloe::IceFloeCylindrical)
    return iceFloeCircumreference(icefloe)*icefloe.thickness
end

export iceFloeVolume
"Returns the volume of the ice floe"
function iceFloeVolume(icefloe::IceFloeCylindrical)
    return iceFloeHorizontalSurfaceArea(icefloe)*icefloe.thickness
end

export iceFloeMass
"Returns the mass of the ice floe"
function iceFloeMass(icefloe::IceFloeCylindrical)
    return iceFloeVolume(icefloe)*icefloe.density
end

export iceFloeMomentOfInertia
"Returns the moment of inertia of the ice floe"
function iceFloeMomentOfInertia(icefloe::IceFloeCylindrical)
    return 0.5*iceFloeMass(icefloe)*icefloe.areal_radius^2.
end

export convertIceFloeDataToArrays
"""
Gathers all ice-floe parameters from the `IceFloeCylindrical` type in continuous 
arrays in an `IceFloeArrays` structure.
"""
function convertIceFloeDataToArrays(simulation::Simulation)

    ifarr = IceFloeArrays(
                          Array(Float64, length(simulation.ice_floes)),

                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),

                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),

                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),

                          Array(Int, length(simulation.ice_floes)),
                          Array(Int, length(simulation.ice_floes)),

                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),
                          Array(Float64, length(simulation.ice_floes)),

                          Array(Float64, length(simulation.ice_floes))
                         )

    # fill arrays
    for i=1:length(simulation.ice_floes)
        ifarr.density[i] = simulation.ice_floes[i].density

        ifarr.thickness[i] = simulation.ice_floes[i].thickness
        ifarr.contact_radius[i] = simulation.ice_floes[i].contact_radius
        ifarr.areal_radius[i] = simulation.ice_floes[i].areal_radius
        ifarr.circumreference[i] = simulation.ice_floes[i].circumreference
        ifarr.horizontal_surface_area[i] =
            simulation.ice_floes[i].horizontal_surface_area
        ifarr.side_surface_area[i] = simulation.ice_floes[i].side_surface_area
        ifarr.volume[i] = simulation.ice_floes[i].volume
        ifarr.mass[i] = simulation.ice_floes[i].mass
        ifarr.moment_of_inertia[i] = simulation.ice_floes[i].moment_of_inertia

        ifarr.lin_pos[1:2, i] = simulation.ice_floes[i].lin_pos
        ifarr.lin_vel[1:2, i] = simulation.ice_floes[i].lin_vel
        ifarr.lin_acc[1:2, i] = simulation.ice_floes[i].lin_acc
        ifarr.force[1:2, i] = simulation.ice_floes[i].force

        ifarr.ang_pos[3, i] = simulation.ice_floes[i].ang_pos
        ifarr.ang_vel[3, i] = simulation.ice_floes[i].ang_vel
        ifarr.ang_acc[3, i] = simulation.ice_floes[i].ang_acc
        ifarr.torque[3, i] = simulation.ice_floes[i].torque

        ifarr.fixed[i] = Int(simulation.ice_floes[i].fixed)
        ifarr.rotating[i] = Int(simulation.ice_floes[i].rotating)

        ifarr.contact_stiffness_normal[i] = 
            simulation.ice_floes[i].contact_stiffness_normal
        ifarr.contact_stiffness_tangential[i] = 
            simulation.ice_floes[i].contact_stiffness_tangential
        ifarr.contact_viscosity_normal[i] = 
            simulation.ice_floes[i].contact_viscosity_normal
        ifarr.contact_viscosity_tangential[i] = 
            simulation.ice_floes[i].contact_viscosity_tangential
        ifarr.contact_static_friction[i] = 
            simulation.ice_floes[i].contact_static_friction
        ifarr.contact_dynamic_friction[i] = 
            simulation.ice_floes[i].contact_dynamic_friction

        ifarr.pressure[i] = simulation.ice_floes[i].pressure
    end

    return ifarr
end

export iceFloeKineticTranslationalEnergy
"Returns the translational kinetic energy of the ice floe"
function iceFloeKineticTranslationalEnergy(icefloe::IceFloeCylindrical)
    return 0.5*icefloe.mass*norm(icefloe.lin_vel)^2.
end

export totalIceFloeKineticTranslationalEnergy
"""
Returns the sum of translational kinetic energies of all ice floes in a 
simulation
"""
function totalIceFloeKineticTranslationalEnergy(simulation::Simulation)
    E_sum = 0.
    for icefloe in simulation.ice_floes
        E_sum += iceFloeKineticTranslationalEnergy(icefloe)
    end
    return E_sum
end

export iceFloeKineticRotationalEnergy
"Returns the rotational kinetic energy of the ice floe"
function iceFloeKineticRotationalEnergy(icefloe::IceFloeCylindrical)
    return 0.5*icefloe.moment_of_inertia*norm(icefloe.ang_vel)^2.
end

export totalIceFloeKineticRotationalEnergy
"""
Returns the sum of rotational kinetic energies of all ice floes in a simulation
"""
function totalIceFloeKineticRotationalEnergy(simulation::Simulation)
    E_sum = 0.
    for icefloe in simulation.ice_floes
        E_sum += iceFloeKineticRotationalEnergy(icefloe)
    end
    return E_sum
end
