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
    circumreference::float
    horizontal_surface_area::float
    side_surface_area::float
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
    enabled::Bool

    # Rheological parameters
    contact_stiffness_normal::float
    contact_stiffness_tangential::float
    contact_viscosity_normal::float
    contact_viscosity_tangential::float
    contact_static_friction::float
    contact_dynamic_friction::float

    youngs_modulus::float
    poissons_ratio::float
    tensile_strength::float
    tensile_heal_rate::float
    compressive_strength_prefactor::float

    # Ocean/atmosphere interaction parameters
    ocean_drag_coeff_vert::float
    ocean_drag_coeff_horiz::float
    atmosphere_drag_coeff_vert::float
    atmosphere_drag_coeff_horiz::float

    # Interaction
    pressure::float
    n_contacts::Int
    ocean_grid_pos::Array{Int, 1}
    atmosphere_grid_pos::Array{Int, 1}
    contacts::Array{Int, 1}
    contact_parallel_displacement::Array{Array{Float64, 1}, 1}
    contact_age::Array{Float64, 1}

    granular_stress::vector
    ocean_stress::vector
    atmosphere_stress::vector
end

# Type for gathering data from ice floe objects into single arrays
type IceFloeArrays

    # Material properties
    density

    # Geometrical parameters
    thickness
    contact_radius
    areal_radius
    circumreference
    horizontal_surface_area
    side_surface_area
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
    enabled

    # Rheological parameters
    contact_stiffness_normal
    contact_stiffness_tangential
    contact_viscosity_normal
    contact_viscosity_tangential
    contact_static_friction
    contact_dynamic_friction

    youngs_modulus
    poissons_ratio
    tensile_strength
    compressive_strength_prefactor

    ocean_drag_coeff_vert
    ocean_drag_coeff_horiz
    atmosphere_drag_coeff_vert
    atmosphere_drag_coeff_horiz

    pressure
    n_contacts

    granular_stress
    ocean_stress
    atmosphere_stress
end

#=
Type containing all relevant data from MOM6 NetCDF files.  The ocean grid is a 
staggered of Arakawa-B type, with south-west convention centered on the 
h-points.  During read, the velocities are interpolated to the cell corners 
(q-points).

    q(  i,j+1) ------------------ q(i+1,j+1)
         |                             |
         |                             |
         |                             |
         |                             |
         |         h(  i,  j)          |
         |                             |
         |                             |
         |                             |
         |                             |
    q(  i,  j) ------------------ q(i+1,  j)

# Fields
* `input_file::String`: path to input NetCDF file
* `time::Array{Float64, 1}`: time in days
* `xq::Array{Float64, 2}`: nominal longitude of q-points [degrees_E]
* `yq::Array{Float64, 2}`: nominal latitude of q-points [degrees_N]
* `xh::Array{Float64, 2}`: nominal longitude of h-points [degrees_E]
* `yh::Array{Float64, 2}`: nominal latitude of h-points [degrees_N]
* `zl::Array{Float64, 1}`: layer target potential density [kg m^-3]
* `zi::Array{Float64, 1}`: interface target potential density [kg m^-3]
* `u::Array{Float64, Int}`: zonal velocity (positive towards west) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `v::Array{Float64, Int}`: meridional velocity (positive towards north) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `h::Array{Float64, Int}`: layer thickness [m], dimensions correspond to 
    placement in `[xh, yh, zl, time]`.
* `e::Array{Float64, Int}`: interface height relative to mean sea level [m],  
    dimensions correspond to placement in `[xh, yh, zi, time]`.
* `ice_floe_list::Array{Float64, Int}`: indexes of ice floes contained in the 
    ocean grid cells.
=#
type Ocean
    input_file::Any

    time::Array{Float64, 1}

    # q-point (cell corner) positions
    xq::Array{Float64, 2}
    yq::Array{Float64, 2}

    # h-point (cell center) positions
    xh::Array{Float64, 2}
    yh::Array{Float64, 2}

    # Vertical positions
    zl::Array{Float64, 1}
    zi::Array{Float64, 1}
    
    # Field values
    u::Array{Float64, 4}
    v::Array{Float64, 4}
    h::Array{Float64, 4}
    e::Array{Float64, 4}

    ice_floe_list::Array{Array{Int, 1}, 2}
end

#=
The atmosphere grid is a staggered of Arakawa-B type, with south-west convention 
centered on the h-points.  During read, the velocities are interpolated to the 
cell corners (q-points).

    q(  i,j+1) ------------------ q(i+1,j+1)
         |                             |
         |                             |
         |                             |
         |                             |
         |         h(  i,  j)          |
         |                             |
         |                             |
         |                             |
         |                             |
    q(  i,  j) ------------------ q(i+1,  j)

# Fields
* `input_file::String`: path to input NetCDF file
* `time::Array{Float64, 1}`: time in days
* `xq::Array{Float64, 1}`: nominal longitude of q-points [degrees_E]
* `yq::Array{Float64, 1}`: nominal latitude of q-points [degrees_N]
* `xh::Array{Float64, 1}`: nominal longitude of h-points [degrees_E]
* `yh::Array{Float64, 1}`: nominal latitude of h-points [degrees_N]
* `zl::Array{Float64, 1}`: vertical position [m]
* `u::Array{Float64, Int}`: zonal velocity (positive towards west) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `v::Array{Float64, Int}`: meridional velocity (positive towards north) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `ice_floe_list::Array{Float64, Int}`: interface height relative to mean sea 
    level [m],  dimensions correspond to placement in `[xh, yh, zi, time]`.
=#
type Atmosphere
    input_file::Any

    time::Array{Float64, 1}

    # q-point (cell corner) positions
    xq::Array{Float64, 2}
    yq::Array{Float64, 2}

    # h-point (cell center) positions
    xh::Array{Float64, 2}
    yh::Array{Float64, 2}

    # Vertical positions
    zl::Array{Float64, 1}
    
    # Field values
    u::Array{Float64, 4}
    v::Array{Float64, 4}

    ice_floe_list::Array{Array{Int, 1}, 2}

    # If true the grid positions are identical to the ocean grid
    collocated_with_ocean_grid::Bool
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
    file_time_since_output_file::Float64

    ice_floes::Array{IceFloeCylindrical, 1}

    ocean::Ocean
    atmosphere::Atmosphere
end
