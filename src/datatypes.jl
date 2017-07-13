## Particle composite types
type IceFloeCylindrical

    # Material properties
    density::Float64

    # Geometrical parameters
    thickness::Float64
    contact_radius::Float64
    areal_radius::Float64
    circumreference::Float64
    horizontal_surface_area::Float64
    side_surface_area::Float64
    volume::Float64
    mass::Float64
    moment_of_inertia::Float64

    # Linear kinematic degrees of freedom along horizontal plane
    lin_pos::Vector{Float64}
    lin_vel::Vector{Float64}
    lin_acc::Vector{Float64}
    force::Vector{Float64}

    # Angular kinematic degrees of freedom for vertical rotation around center
    ang_pos::Float64
    ang_vel::Float64
    ang_acc::Float64
    torque::Float64

    # Kinematic constraint flags
    fixed::Bool
    rotating::Bool
    enabled::Bool

    # Rheological parameters
    contact_stiffness_normal::Float64
    contact_stiffness_tangential::Float64
    contact_viscosity_normal::Float64
    contact_viscosity_tangential::Float64
    contact_static_friction::Float64
    contact_dynamic_friction::Float64

    youngs_modulus::Float64
    poissons_ratio::Float64
    tensile_strength::Float64
    tensile_heal_rate::Float64
    compressive_strength_prefactor::Float64

    # Ocean/atmosphere interaction parameters
    ocean_drag_coeff_vert::Float64
    ocean_drag_coeff_horiz::Float64
    atmosphere_drag_coeff_vert::Float64
    atmosphere_drag_coeff_horiz::Float64

    # Interaction
    pressure::Float64
    n_contacts::Int
    ocean_grid_pos::Vector{Int}
    atmosphere_grid_pos::Vector{Int}
    contacts::Vector{Int}
    contact_parallel_displacement::Vector{Vector{Float64}}
    contact_age::Vector{Float64}

    granular_stress::Vector{Float64}
    ocean_stress::Vector{Float64}
    atmosphere_stress::Vector{Float64}
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

    time::Vector{Float64}

    # q-point (cell corner) positions
    xq::Array{Float64, 2}
    yq::Array{Float64, 2}

    # h-point (cell center) positions
    xh::Array{Float64, 2}
    yh::Array{Float64, 2}

    # Vertical positions
    zl::Vector{Float64}
    zi::Vector{Float64}
    
    # Field values
    u::Array{Float64, 4}
    v::Array{Float64, 4}
    h::Array{Float64, 4}
    e::Array{Float64, 4}

    ice_floe_list::Array{Vector{Int}, 2}
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
* `time::Vector{Float64}`: time in days
* `xq::Array{Float64, 2}`: nominal longitude of q-points [degrees_E]
* `yq::Array{Float64, 2}`: nominal latitude of q-points [degrees_N]
* `xh::Array{Float64, 2}`: nominal longitude of h-points [degrees_E]
* `yh::Array{Float64, 2}`: nominal latitude of h-points [degrees_N]
* `zl::Vector{Float64}`: vertical position [m]
* `u::Array{Float64, Int}`: zonal velocity (positive towards west) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `v::Array{Float64, Int}`: meridional velocity (positive towards north) [m/s], 
    dimensions correspond to placement in `[xq, yq, zl, time]`.
* `ice_floe_list::Array{Float64, Int}`: interface height relative to mean sea 
    level [m],  dimensions correspond to placement in `[xh, yh, zi, time]`.
=#
type Atmosphere
    input_file::Any

    time::Vector{Float64}

    # q-point (cell corner) positions
    xq::Array{Float64, 2}
    yq::Array{Float64, 2}

    # h-point (cell center) positions
    xh::Array{Float64, 2}
    yh::Array{Float64, 2}

    # Vertical positions
    zl::Vector{Float64}
    
    # Field values
    u::Array{Float64, 4}
    v::Array{Float64, 4}

    ice_floe_list::Array{Vector{Int}, 2}

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

    ice_floes::Vector{IceFloeCylindrical}

    ocean::Ocean
    atmosphere::Atmosphere

    Nc_max::Int
end
