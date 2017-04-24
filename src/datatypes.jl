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

    pressure::float
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

    pressure
end

#=
Type containing all relevant data from MOM6 NetCDF file.  The ocean grid is a 
staggered of Arakawa-C type, with north-east convention centered on the 
h-points.

    q(i-1,  j) --- v(  i,  j) --- q(  i,  j)
         |                             |
         |                             |
         |                             |
         |                             |
    u(i-1,  j)     h(  i,  j)     u(  i,  j)
         |                             |
         |                             |
         |                             |
         |                             |
    q(i-1,j-1) --- v(  i,j-1) --- q(  i,j-1)

Source: 
https://mom6.readthedocs.io/en/latest/api/generated/pages/Horizontal_indexing.html

# Fields
* `input_file::String`: path to input NetCDF file
* `time::Array{Float64, 1}`: time in days
* `xq::Array{Float64, 1}`: nominal longitude of q-points [degrees_E]
* `yq::Array{Float64, 1}`: nominal latitude of q-points [degrees_N]
* `xh::Array{Float64, 1}`: nominal longitude of h-points [degrees_E]
* `yh::Array{Float64, 1}`: nominal latitude of h-points [degrees_N]
* `zl::Array{Float64, 1}`: layer target potential density [kg m^-3]
* `zi::Array{Float64, 1}`: interface target potential density [kg m^-3]
* `u::Array{Float64, Int}`: zonal velocity (positive towards west) [m/s], 
    dimensions correspond to placement in `[xq, yh, zl, time]`.
* `v::Array{Float64, Int}`: meridional velocity (positive towards north) [m/s], 
    dimensions correspond to placement in `[xh, yq, zl, time]`.
* `h::Array{Float64, Int}`: layer thickness [m], dimensions correspond to 
    placement in `[xh, yh, zl, time]`.
* `e::Array{Float64, Int}`: interface height relative to mean sea level [m],  
    dimensions correspond to placement in `[xh, yq, zi, time]`.
=#
type Ocean
    input_file::String

    time::Array{Float64, 1}

    # q-point (cell corner) positions
    xq::Array{Float64, 1}
    yq::Array{Float64, 1}

    # h-point (cell center) positions
    xh::Array{Float64, 1}
    yh::Array{Float64, 1}

    # Vertical positions
    zl::Array{Float64, 1}
    zi::Array{Float64, 1}
    
    # Field values
    u::Array{Float64, 4}
    v::Array{Float64, 4}
    h::Array{Float64, 4}
    e::Array{Float64, 4}
end

# Top-level simulation type
type Simulation
    id::String

    time_iteration::Int
    time::Float64
    time_total::Float64
    time_step::Float64
    file_time_step::Float64
    file_number::Int

    ice_floes::Array{IceFloeCylindrical, 1}
    contact_pairs::Array{Array{Int, 1}, 1}
    overlaps::Array{Array{Float64, 1}, 1}

    ocean::Ocean
end
