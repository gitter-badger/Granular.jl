## Manage icefloes in the model
import PyPlot

export addIceFloeCylindrical!
"""
    function addIceFloeCylindrical!(simulation, lin_pos, contact_radius,
                                    thickness[, areal_radius, lin_vel, lin_acc,
                                    force, ang_pos, ang_vel, ang_acc, torque,
                                    density, contact_stiffness_normal,
                                    contact_stiffness_tangential,
                                    contact_viscosity_normal,
                                    contact_viscosity_tangential,
                                    contact_static_friction,
                                    contact_dynamic_friction,
                                    youngs_modulus, poissons_ratio,
                                    tensile_strength, tensile_heal_rate,
                                    compressive_strength_prefactor,
                                    ocean_drag_coeff_vert,
                                    ocean_drag_coeff_horiz,
                                    atmosphere_drag_coeff_vert,
                                    atmosphere_drag_coeff_horiz,
                                    pressure, fixed, rotating, enabled, verbose,
                                    ocean_grid_pos, atmosphere_grid_pos,
                                    n_contact, granular_stress, ocean_stress,
                                    atmosphere_stress])

Creates and adds a cylindrical icefloe to a simulation. Most of the arguments 
are optional, and come with default values.  The only required arguments are 
`simulation`, `lin_pos`, `contact_radius`, and `thickness`.

# Arguments
- `simulation::Simulation`: the simulation object where the ice floe should be
    added to.
- `lin_pos::Vector{Float64}`: linear position of ice-floe center [m].
- `contact_radius::Float64`: ice-floe radius for granular interaction [m].
- `thickness::Float64`: ice-floe thickness [m].
- `areal_radius = false`: ice-floe radius for determining sea-ice concentration
    [m].
- `lin_vel::Vector{Float64} = [0., 0.]`: linear velocity [m/s].
- `lin_acc::Vector{Float64} = [0., 0.]`: linear acceleration [m/s^2].
- `force::Vector{Float64} = [0., 0.]`: linear force balance [N].
- `ang_pos::Float64 = 0.`: angular position around its center vertical axis
    [rad].
- `ang_vel::Float64 = 0.`: angular velocity around its center vertical axis
    [rad/s].
- `ang_acc::Float64 = 0.`: angular acceleration around its center vertical axis
    [rad/s^2].
- `torque::Float64 = 0.`: torque around its center vertical axis [N*m]
- `density::Float64 = 934.`: ice-floe mean density [kg/m^3].
- `contact_stiffness_normal::Float64 = 1e7`: contact-normal stiffness [N/m];
    overridden if `youngs_modulus` is set to a positive value.
- `contact_stiffness_tangential::Float64 = 0.`: contact-tangential stiffness
    [N/m]; overridden if `youngs_modulus` is set to a positive value.
- `contact_viscosity_normal::Float64 = 0.`: contact-normal viscosity [N/m/s].
- `contact_viscosity_tangential::Float64 = 0.`: contact-tangential viscosity
    [N/m/s].
- `contact_static_friction::Float64 = 0.4`: contact static Coulomb frictional
    coefficient [-].
- `contact_dynamic_friction::Float64 = 0.4`: contact dynamic Coulomb frictional
    coefficient [-].
- `youngs_modulus::Float64 = 2e7`: elastic modulus [Pa]; overrides any value
    set for `k_n`.
- `poissons_ratio::Float64 = 0.185`: Poisson's ratio, used to determine the
    contact-tangential stiffness from `youngs_modulus` [-].
- `tensile_strength::Float64 = 0.`: contact-tensile (cohesive) bond strength
    [Pa].
- `tensile_heal_rate::Float64 = 0.`: rate at which contact-tensile bond strength
    is obtained [1/s].
- `compressive_strength_prefactor::Float64 = 1285e3`: maximum compressive
    strength on granular contact (not currently enforced) [m*Pa].
- `ocean_drag_coeff_vert::Float64 = 0.85`: vertical drag coefficient for ocean
    against ice-floe sides [-].
- `ocean_drag_coeff_horiz::Float64 = 5e-4`: horizontal drag coefficient for
    ocean against ice-floe bottom [-].
- `atmosphere_drag_coeff_vert::Float64 = 0.4`: vertical drag coefficient for
    atmosphere against ice-floe sides [-].
- `atmosphere_drag_coeff_horiz::Float64 = 2.5e-4`: horizontal drag coefficient
    for atmosphere against ice-floe bottom [-].
- `pressure::Float64 = 0.`: current compressive stress on ice floe [Pa].
- `fixed::Bool = false`: ice floe is fixed in space.
- `rotating::Bool = true`: ice floe is allowed to rotate.
- `enabled::Bool = true`: ice floe interacts with other ice floes.
- `verbose::Bool = true`: display diagnostic information during the function
    call.
- `ocean_grid_pos::Array{Int, 1} = [0, 0]`: position of ice floe in the ocean
    grid.
- `atmosphere_grid_pos::Array{Int, 1} = [0, 0]`: position of ice floe in the
    atmosphere grid.
- `n_contacts::Int = 0`: number of contacts with other ice floes.
- `granular_stress::Vector{Float64} = [0., 0.]`: resultant stress on ice floe
    from granular interactions [Pa].
- `ocean_stress::Vector{Float64} = [0., 0.]`: resultant stress on ice floe from
    ocean drag [Pa].
- `atmosphere_stress::Vector{Float64} = [0., 0.]`: resultant stress on ice floe
    from atmosphere drag [Pa].

# Examples
The most basic example adds a new ice floe to the simulation `sim`, with a 
center at `[1., 2., 3.]`, a radius of `1.` meter, and a thickness of `0.5` 
meter.

    SeaIce.addIceFloeCylindrical(sim, [1., 2., 3.], 1., .5)
"""
function addIceFloeCylindrical!(simulation::Simulation,
                                lin_pos::Vector{Float64},
                                contact_radius::Float64,
                                thickness::Float64;
                                areal_radius = false,
                                lin_vel::Vector{Float64} = [0., 0.],
                                lin_acc::Vector{Float64} = [0., 0.],
                                force::Vector{Float64} = [0., 0.],
                                ang_pos::Float64 = 0.,
                                ang_vel::Float64 = 0.,
                                ang_acc::Float64 = 0.,
                                torque::Float64 = 0.,
                                density::Float64 = 934.,
                                contact_stiffness_normal::Float64 = 1e7,
                                contact_stiffness_tangential::Float64 = 0.,
                                contact_viscosity_normal::Float64 = 0.,
                                contact_viscosity_tangential::Float64 = 0.,
                                contact_static_friction::Float64 = 0.4,
                                contact_dynamic_friction::Float64 = 0.4,
                                youngs_modulus::Float64 = 2e7,
                                #youngs_modulus::Float64 = 2e9,  # Hopkins 2004
                                poissons_ratio::Float64 = 0.185,  # Hopkins 2004
                                #tensile_strength::Float64 = 500e3,  # Hopkins2004
                                tensile_strength::Float64 = 0.,
                                tensile_heal_rate::Float64 = 0.,
                                compressive_strength_prefactor::Float64 = 1285e3,  
                                    # Hopkins 2004
                                ocean_drag_coeff_vert::Float64 = 0.85, # H&C 2011
                                ocean_drag_coeff_horiz::Float64 = 5e-4, # H&C 2011
                                atmosphere_drag_coeff_vert::Float64 = 0.4,
                                    # H&C 2011
                                atmosphere_drag_coeff_horiz::Float64 = 2.5e-4,
                                    # H&C2011
                                pressure::Float64 = 0.,
                                fixed::Bool = false,
                                rotating::Bool = true,
                                enabled::Bool = true,
                                verbose::Bool = true,
                                ocean_grid_pos::Array{Int, 1} = [0, 0],
                                atmosphere_grid_pos::Array{Int, 1} = [0, 0],
                                n_contacts::Int = 0,
                                granular_stress::Vector{Float64} = [0., 0.],
                                ocean_stress::Vector{Float64} = [0., 0.],
                                atmosphere_stress::Vector{Float64} = [0., 0.])

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
    if contact_radius <= 0.0
        error("Radius must be greater than 0.0 (radius = $contact_radius m)")
    end
    if density <= 0.0
        error("Density must be greater than 0.0 (density = $density kg/m^3)")
    end

    if !areal_radius
        areal_radius = contact_radius
    end

    contacts::Array{Int, 1} = zeros(Int, simulation.Nc_max)
    contact_parallel_displacement =
        Vector{Vector{Float64}}(simulation.Nc_max)
        contact_age::Vector{Float64} = zeros(Float64, simulation.Nc_max)
    for i=1:simulation.Nc_max
        contact_parallel_displacement[i] = zeros(2)
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
                                 enabled,

                                 contact_stiffness_normal,
                                 contact_stiffness_tangential,
                                 contact_viscosity_normal,
                                 contact_viscosity_tangential,
                                 contact_static_friction,
                                 contact_dynamic_friction,

                                 youngs_modulus,
                                 poissons_ratio,
                                 tensile_strength,
                                 tensile_heal_rate,
                                 compressive_strength_prefactor,

                                 ocean_drag_coeff_vert,
                                 ocean_drag_coeff_horiz,
                                 atmosphere_drag_coeff_vert,
                                 atmosphere_drag_coeff_horiz,

                                 pressure,
                                 n_contacts,
                                 ocean_grid_pos,
                                 atmosphere_grid_pos,
                                 contacts,
                                 contact_parallel_displacement,
                                 contact_age,

                                 granular_stress,
                                 ocean_stress,
                                 atmosphere_stress
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
    nothing
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
                          Array{Float64}(length(simulation.ice_floes)),

                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),

                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),

                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),

                          Array{Int}(length(simulation.ice_floes)),
                          Array{Int}(length(simulation.ice_floes)),
                          Array{Int}(length(simulation.ice_floes)),

                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),

                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),

                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Float64}(length(simulation.ice_floes)),

                          Array{Float64}(length(simulation.ice_floes)),
                          Array{Int}(length(simulation.ice_floes)),

                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
                          zeros(Float64, 3, length(simulation.ice_floes)),
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
        ifarr.enabled[i] = Int(simulation.ice_floes[i].enabled)

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

        ifarr.youngs_modulus[i] = simulation.ice_floes[i].youngs_modulus
        ifarr.poissons_ratio[i] = simulation.ice_floes[i].poissons_ratio
        ifarr.tensile_strength[i] = simulation.ice_floes[i].tensile_strength
        ifarr.compressive_strength_prefactor[i] = 
            simulation.ice_floes[i].compressive_strength_prefactor

        ifarr.ocean_drag_coeff_vert[i] = 
            simulation.ice_floes[i].ocean_drag_coeff_vert
        ifarr.ocean_drag_coeff_horiz[i] = 
            simulation.ice_floes[i].ocean_drag_coeff_horiz
        ifarr.atmosphere_drag_coeff_vert[i] = 
            simulation.ice_floes[i].atmosphere_drag_coeff_vert
        ifarr.atmosphere_drag_coeff_horiz[i] = 
            simulation.ice_floes[i].atmosphere_drag_coeff_horiz

        ifarr.pressure[i] = simulation.ice_floes[i].pressure
        ifarr.n_contacts[i] = simulation.ice_floes[i].n_contacts

        ifarr.granular_stress[1:2, i] = simulation.ice_floes[i].granular_stress
        ifarr.ocean_stress[1:2, i] = simulation.ice_floes[i].ocean_stress
        ifarr.atmosphere_stress[1:2, i] =
            simulation.ice_floes[i].atmosphere_stress
    end

    return ifarr
end

function deleteIceFloeArrays!(ifarr::IceFloeArrays)
    f1 = zeros(1)
    f2 = zeros(1,1)
    i1 = zeros(Int, 1)

    ifarr.density = f1

    ifarr.thickness = f1
    ifarr.contact_radius = f1
    ifarr.areal_radius = f1
    ifarr.circumreference = f1
    ifarr.horizontal_surface_area = f1
    ifarr.side_surface_area = f1
    ifarr.volume = f1
    ifarr.mass = f1
    ifarr.moment_of_inertia = f1

    ifarr.lin_pos = f2
    ifarr.lin_vel = f2
    ifarr.lin_acc = f2
    ifarr.force = f2

    ifarr.ang_pos = f2
    ifarr.ang_vel = f2
    ifarr.ang_acc = f2
    ifarr.torque = f2

    ifarr.fixed = i1
    ifarr.rotating = i1
    ifarr.enabled = i1

    ifarr.contact_stiffness_normal = f1
    ifarr.contact_stiffness_tangential = f1
    ifarr.contact_viscosity_normal = f1
    ifarr.contact_viscosity_tangential = f1
    ifarr.contact_static_friction = f1
    ifarr.contact_dynamic_friction = f1

    ifarr.youngs_modulus = f1
    ifarr.poissons_ratio = f1
    ifarr.tensile_strength = f1
    ifarr.compressive_strength_prefactor = f1

    ifarr.ocean_drag_coeff_vert = f1
    ifarr.ocean_drag_coeff_horiz = f1
    ifarr.atmosphere_drag_coeff_vert = f1
    ifarr.atmosphere_drag_coeff_horiz = f1

    ifarr.pressure = f1
    ifarr.n_contacts = i1

    ifarr.granular_stress = f2
    ifarr.ocean_stress = f2
    ifarr.atmosphere_stress = f2

    gc()
    nothing
end

export printIceFloeInfo
"""
    printIceFloeInfo(icefloe::IceFloeCylindrical)

Prints the contents of an ice floe to stdout in a formatted style.
"""
function printIceFloeInfo(f::IceFloeCylindrical)
    println("  density:                 $(f.density) kg/m^3")
    println("  thickness:               $(f.thickness) m")
    println("  contact_radius:          $(f.contact_radius) m")
    println("  areal_radius:            $(f.areal_radius) m")
    println("  circumreference:         $(f.circumreference) m")
    println("  horizontal_surface_area: $(f.horizontal_surface_area) m^2")
    println("  side_surface_area:       $(f.side_surface_area) m^2")
    println("  volume:                  $(f.volume) m^3")
    println("  mass:                    $(f.mass) kg")
    println("  moment_of_inertia:       $(f.moment_of_inertia) kg*m^2\n")

    println("  lin_pos: $(f.lin_pos) m")
    println("  lin_vel: $(f.lin_vel) m/s")
    println("  lin_acc: $(f.lin_acc) m/s^2")
    println("  force:   $(f.force) N\n")

    println("  ang_pos: $(f.ang_pos) rad")
    println("  ang_vel: $(f.ang_vel) rad/s")
    println("  ang_acc: $(f.ang_acc) rad/s^2")
    println("  torque:  $(f.torque) N*m\n")

    println("  fixed:    $(f.fixed)")
    println("  rotating: $(f.rotating)")
    println("  enabled:  $(f.enabled)\n")

    println("  k_n:     $(f.contact_stiffness_normal) N/m")
    println("  k_t:     $(f.contact_stiffness_tangential) N/m")
    println("  γ_n: $(f.contact_viscosity_normal) N/(m/s)")
    println("  γ_t: $(f.contact_viscosity_tangential) N/(m/s)")
    println("  μ_s:    $(f.contact_static_friction)")
    println("  μ_d:    $(f.contact_dynamic_friction)\n")

    println("  E:      $(f.youngs_modulus) Pa")
    println("  ν:      $(f.poissons_ratio)")
    println("  σ_t:    $(f.tensile_strength) Pa")
    println("  c(σ_c): $(f.compressive_strength_prefactor) m^0.5 Pa\n")

    println("  c_o_v:  $(f.ocean_drag_coeff_vert)")
    println("  c_o_h:  $(f.ocean_drag_coeff_horiz)")
    println("  c_a_v:  $(f.atmosphere_drag_coeff_vert)")
    println("  c_a_h:  $(f.atmosphere_drag_coeff_horiz)\n")

    println("  pressure:   $(f.pressure) Pa")
    println("  n_contacts: $(f.n_contacts)")
    println("  contacts:   $(f.contacts)")
    println("  delta_t:    $(f.contact_parallel_displacement)\n")

    println("  granular_stress:   $(f.granular_stress) Pa")
    println("  ocean_stress:      $(f.ocean_stress) Pa")
    println("  atmosphere_stress: $(f.atmosphere_stress) Pa")
    nothing
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

export compareIceFloes
"""
    compareIceFloes(if1::IceFloeCylindrical, if2::IceFloeCylindrical)

Compare values of two ice floe objects using the `Base.Test` framework.
"""
function compareIceFloes(if1::IceFloeCylindrical, if2::IceFloeCylindrical)

    Base.Test.@test if1.density ≈ if2.density
    Base.Test.@test if1.thickness ≈ if2.thickness
    Base.Test.@test if1.contact_radius ≈
        if2.contact_radius
    Base.Test.@test if1.areal_radius ≈ if2.areal_radius
    Base.Test.@test if1.circumreference ≈
        if2.circumreference
    Base.Test.@test if1.horizontal_surface_area ≈ if2.horizontal_surface_area
    Base.Test.@test if1.side_surface_area ≈ if2.side_surface_area
    Base.Test.@test if1.volume ≈ if2.volume
    Base.Test.@test if1.mass ≈ if2.mass
    Base.Test.@test if1.moment_of_inertia ≈ if2.moment_of_inertia

    Base.Test.@test if1.lin_pos ≈ if2.lin_pos
    Base.Test.@test if1.lin_vel ≈ if2.lin_vel
    Base.Test.@test if1.lin_acc ≈ if2.lin_acc
    Base.Test.@test if1.force ≈ if2.force

    Base.Test.@test if1.ang_pos ≈ if2.ang_pos
    Base.Test.@test if1.ang_vel ≈ if2.ang_vel
    Base.Test.@test if1.ang_acc ≈ if2.ang_acc
    Base.Test.@test if1.torque ≈ if2.torque

    Base.Test.@test if1.fixed == if2.fixed
    Base.Test.@test if1.rotating == if2.rotating
    Base.Test.@test if1.enabled == if2.enabled

    Base.Test.@test if1.contact_stiffness_normal ≈ if2.contact_stiffness_normal
    Base.Test.@test if1.contact_stiffness_tangential ≈ 
        if2.contact_stiffness_tangential
    Base.Test.@test if1.contact_viscosity_normal ≈ if2.contact_viscosity_normal
    Base.Test.@test if1.contact_viscosity_tangential ≈ 
        if2.contact_viscosity_tangential
    Base.Test.@test if1.contact_static_friction ≈ if2.contact_static_friction
    Base.Test.@test if1.contact_dynamic_friction ≈ if2.contact_dynamic_friction

    Base.Test.@test if1.youngs_modulus ≈ if2.youngs_modulus
    Base.Test.@test if1.poissons_ratio ≈ if2.poissons_ratio
    Base.Test.@test if1.tensile_strength ≈ if2.tensile_strength
    Base.Test.@test if1.tensile_heal_rate ≈ if2.tensile_heal_rate
    Base.Test.@test if1.compressive_strength_prefactor ≈
        if2.compressive_strength_prefactor

    Base.Test.@test if1.ocean_drag_coeff_vert ≈ if2.ocean_drag_coeff_vert
    Base.Test.@test if1.ocean_drag_coeff_horiz ≈ if2.ocean_drag_coeff_horiz
    Base.Test.@test if1.atmosphere_drag_coeff_vert ≈ 
        if2.atmosphere_drag_coeff_vert
    Base.Test.@test if1.atmosphere_drag_coeff_horiz ≈ 
        if2.atmosphere_drag_coeff_horiz

    Base.Test.@test if1.pressure ≈ if2.pressure
    Base.Test.@test if1.n_contacts == if2.n_contacts
    Base.Test.@test if1.ocean_grid_pos == if2.ocean_grid_pos
    Base.Test.@test if1.contacts == if2.contacts
    Base.Test.@test if1.contact_parallel_displacement == 
        if2.contact_parallel_displacement
    Base.Test.@test if1.contact_age ≈ if2.contact_age

    Base.Test.@test if1.granular_stress ≈ if2.granular_stress
    Base.Test.@test if1.ocean_stress ≈ if2.ocean_stress
    Base.Test.@test if1.atmosphere_stress ≈ if2.atmosphere_stress
    nothing
end

export plotIceFloeSizeDistribution
"""
    plotIceFloeSizeDistribution(simulation, [filename_postfix], [nbins],
                                [size_type], [figsize], [filetype])

Plot the ice-floe size distribution as a histogram and save it to the disk.  The 
plot is saved accoring to the simulation id, the optional `filename_postfix` 
string, and the `filetype`, and is written to the current folder.

# Arguments
- `simulation::Simulation`: the simulation object containing the ice floes.
- `filename_postfix::String`: optional string for the output filename.
- `nbins::Int`: number of bins in the histogram (default = 12).
- `size_type::String`: specify whether to use the `contact` or `areal` radius 
    for the ice-floe size.  The default is `contact`.
- `figsize::Tuple`: the output figure size in inches (default = (6,4).
- `filetype::String`: the output file type (default = "png").
- `verbose::String`: show output file as info message in stdout (default = 
    true).
- `skip_fixed::Bool`: ommit ice floes that are fixed in space from the size 
    distribution (default = true).
- `logy::Bool`: plot y-axis in log scale.
"""
function plotIceFloeSizeDistribution(simulation::Simulation;
                                     filename_postfix::String = "",
                                     nbins::Int=12,
                                     size_type::String = "contact",
                                     figsize::Tuple = (6,4),
                                     filetype::String = "png",
                                     verbose::Bool = true,
                                     skip_fixed::Bool = true,
                                     log_y::Bool = true)

    diameters = Float64[]
    for i=1:length(simulation.ice_floes)
        if simulation.ice_floes[i].fixed && skip_fixed
            continue
        end
        if size_type == "contact"
            push!(diameters, simulation.ice_floes[i].contact_radius*2.)
        elseif size_type == "areal"
            push!(diameters, simulation.ice_floes[i].areal_radius*2.)
        else
            error("size_type '$size_type' not understood")
        end
    end
    PyPlot.pygui(false)
    PyPlot.figure(figsize=figsize)
    PyPlot.plt[:hist](diameters, nbins, log=log_y)
    PyPlot.xlabel("Diameter [m]")
    PyPlot.ylabel("Count [-]")
    filename = string(simulation.id * filename_postfix * 
                      "-ice-floe-size-distribution." * filetype)
    PyPlot.savefig(filename)
    if verbose
        info(filename)
    end
end
