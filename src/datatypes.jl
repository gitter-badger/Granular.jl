# Define floating point data types
const float = Float64
const vector = Array{Float64, 1}

## Particle composite types
type IceFloeCylindrical

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
    contact_static_friction::float
    contact_dynamic_friction::float
end

# Type for gathering data from ice floe objects into single arrays
type IceFloeArrays

    # Material properties
    density

    # Geometrical parameters
    thickness
    contact_radius
    areal_radius
    surface_area
    volume
    mass
    moment_of_inertia

    # Linear kinematic degrees of freedom along horizontal plane
    lin_pos
    lin_vel
    lin_acc
    force

    # Angular kinematic degrees of freedom for vertical rotation around center
    ang_pos
    ang_vel
    ang_acc
    torque

    # Kinematic constraint flags
    fixed
    rotating

    # Rheological parameters
    contact_stiffness_normal
    contact_stiffness_tangential
    contact_viscosity_normal
    contact_viscosity_tangential
    contact_static_friction
    contact_dynamic_friction
end

## Top-level simulation type

# Simulation-scope data
type Simulation
    id::String

    time_iteration::Int
    time::Float64
    time_total::Float64
    time_step::Float64
    file_time_step::Float64   # 0.0: no output files
    file_number::Int

    ice_floes::Array{IceFloeCylindrical, 1}
    contact_pairs::Array{Array{Int, 1}, 1}
    overlaps::Array{Array{Float64, 1}, 1}
end

