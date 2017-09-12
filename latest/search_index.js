var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SeaIce.jl-1",
    "page": "Home",
    "title": "SeaIce.jl",
    "category": "section",
    "text": "A Julia package for particle-based simulation of sea-ice dynamics.SeaIce.jl is a flexible and computationally efficient 2d implementation of the discrete element method, made for simulating sea ice in a Lagrangian manner.  Sea-ice floes are represented as particles, which can be forced by ocean and atmospheric velocity fields.  The ice floes interact through elasto-viscous-frictional contact rheologies and obtain time-dependent tensile strength.The source code for SeaIce.jl is hosted on Github.See the Public API Index for the complete list of documented functions and types."
},

{
    "location": "index.html#Author-1",
    "page": "Home",
    "title": "Author",
    "category": "section",
    "text": "Anders Damsgaard, Geophysical Fluid Dynamics Laboratory, Princeton University."
},

{
    "location": "index.html#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "SeaIce.jl is licensed under the GPLv3; see LICENSE for the full license text."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"man/installation.md\",\n    \"man/simple_example.md\",\n]\nDepth = 1"
},

{
    "location": "index.html#Library-Outline-1",
    "page": "Home",
    "title": "Library Outline",
    "category": "section",
    "text": "Pages = [\n    \"lib/public.md\",\n    \"lib/internals.md\",\n]\nDepth = 1"
},

{
    "location": "man/installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "man/installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "SeaIce.jl can be installed directly from the Julia shell by:julia> Pkg.clone(\"git://github.com/anders-dc/SeaIce.jl.git\")This will install the contents of this repository in the folder  ~/.julia/v$(JULIA_VERSION)/SeaIce and install its requirements.  The package JLD  is used for model restarts and is recommended but not required, and is thus not  automatically installed.You can run the package tests, which are contained in the test/ directory, with the following command:julia> Pkg.test(\"SeaIce\")The package can be updated from this repository using:julia> Pkg.update(\"SeaIce\")"
},

{
    "location": "man/simple_example.html#",
    "page": "A simple example",
    "title": "A simple example",
    "category": "page",
    "text": ""
},

{
    "location": "man/simple_example.html#A-simple-example-1",
    "page": "A simple example",
    "title": "A simple example",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public.html#",
    "page": "Public API",
    "title": "Public API",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public.html#Public-API-documentation-1",
    "page": "Public API",
    "title": "Public API documentation",
    "category": "section",
    "text": "Documentation for SeaIce.jl's public interface.See Package-internal documentation for internal package docs."
},

{
    "location": "lib/public.html#main-index-1",
    "page": "Public API",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public.html#SeaIce",
    "page": "Public API",
    "title": "SeaIce",
    "category": "Module",
    "text": "SeaIce.jl\n\nOffline sea-ice dynamics simulator module.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.addAtmosphereDrag!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.addAtmosphereDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between atmosphere and all  ice floes.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.addIceFloe!",
    "page": "Public API",
    "title": "SeaIce.addIceFloe!",
    "category": "Function",
    "text": "addIceFloe!(simulation::Simulation,\n            icefloe::IceFloeCylindrical,\n            verbose::Bool = False)\n\nAdd an icefloe to the simulation object.  If verbose is true, a short  confirmation message will be printed to stdout.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.addIceFloeCylindrical!-Tuple{SeaIce.Simulation,Array{Float64,1},Float64,Float64}",
    "page": "Public API",
    "title": "SeaIce.addIceFloeCylindrical!",
    "category": "Method",
    "text": "function addIceFloeCylindrical!(simulation, lin_pos, contact_radius,\n                                thickness[, areal_radius, lin_vel, lin_acc,\n                                force, ang_pos, ang_vel, ang_acc, torque,\n                                density, contact_stiffness_normal,\n                                contact_stiffness_tangential,\n                                contact_viscosity_normal,\n                                contact_viscosity_tangential,\n                                contact_static_friction,\n                                contact_dynamic_friction,\n                                youngs_modulus, poissons_ratio,\n                                tensile_strength, tensile_heal_rate,\n                                compressive_strength_prefactor,\n                                ocean_drag_coeff_vert,\n                                ocean_drag_coeff_horiz,\n                                atmosphere_drag_coeff_vert,\n                                atmosphere_drag_coeff_horiz,\n                                pressure, fixed, rotating, enabled, verbose,\n                                ocean_grid_pos, atmosphere_grid_pos,\n                                n_contact, granular_stress, ocean_stress,\n                                atmosphere_stress])\n\nCreates and adds a cylindrical icefloe to a simulation. Most of the arguments  are optional, and come with default values.  The only required arguments are  simulation, lin_pos, contact_radius, and thickness.\n\nArguments\n\nsimulation::Simulation: the simulation object where the ice floe should be   added to.\nlin_pos::Vector{Float64}: linear position of ice-floe center [m].\ncontact_radius::Float64: ice-floe radius for granular interaction [m].\nthickness::Float64: ice-floe thickness [m].\nareal_radius = false: ice-floe radius for determining sea-ice concentration   [m].\nlin_vel::Vector{Float64} = [0., 0.]: linear velocity [m/s].\nlin_acc::Vector{Float64} = [0., 0.]: linear acceleration [m/s^2].\nforce::Vector{Float64} = [0., 0.]: linear force balance [N].\nang_pos::Float64 = 0.: angular position around its center vertical axis   [rad].\nang_vel::Float64 = 0.: angular velocity around its center vertical axis   [rad/s].\nang_acc::Float64 = 0.: angular acceleration around its center vertical axis   [rad/s^2].\ntorque::Float64 = 0.: torque around its center vertical axis [N*m]\ndensity::Float64 = 934.: ice-floe mean density [kg/m^3].\ncontact_stiffness_normal::Float64 = 1e7: contact-normal stiffness [N/m];   overridden if youngs_modulus is set to a positive value.\ncontact_stiffness_tangential::Float64 = 0.: contact-tangential stiffness   [N/m]; overridden if youngs_modulus is set to a positive value.\ncontact_viscosity_normal::Float64 = 0.: contact-normal viscosity [N/m/s].\ncontact_viscosity_tangential::Float64 = 0.: contact-tangential viscosity   [N/m/s].\ncontact_static_friction::Float64 = 0.4: contact static Coulomb frictional   coefficient [-].\ncontact_dynamic_friction::Float64 = 0.4: contact dynamic Coulomb frictional   coefficient [-].\nyoungs_modulus::Float64 = 2e7: elastic modulus [Pa]; overrides any value   set for k_n.\npoissons_ratio::Float64 = 0.185: Poisson's ratio, used to determine the   contact-tangential stiffness from youngs_modulus [-].\ntensile_strength::Float64 = 0.: contact-tensile (cohesive) bond strength   [Pa].\ntensile_heal_rate::Float64 = 0.: rate at which contact-tensile bond strength   is obtained [1/s].\ncompressive_strength_prefactor::Float64 = 1285e3: maximum compressive   strength on granular contact (not currently enforced) [m*Pa].\nocean_drag_coeff_vert::Float64 = 0.85: vertical drag coefficient for ocean   against ice-floe sides [-].\nocean_drag_coeff_horiz::Float64 = 5e-4: horizontal drag coefficient for   ocean against ice-floe bottom [-].\natmosphere_drag_coeff_vert::Float64 = 0.4: vertical drag coefficient for   atmosphere against ice-floe sides [-].\natmosphere_drag_coeff_horiz::Float64 = 2.5e-4: horizontal drag coefficient   for atmosphere against ice-floe bottom [-].\npressure::Float64 = 0.: current compressive stress on ice floe [Pa].\nfixed::Bool = false: ice floe is fixed in space.\nrotating::Bool = true: ice floe is allowed to rotate.\nenabled::Bool = true: ice floe interacts with other ice floes.\nverbose::Bool = true: display diagnostic information during the function   call.\nocean_grid_pos::Array{Int, 1} = [0, 0]: position of ice floe in the ocean   grid.\natmosphere_grid_pos::Array{Int, 1} = [0, 0]: position of ice floe in the   atmosphere grid.\nn_contacts::Int = 0: number of contacts with other ice floes.\ngranular_stress::Vector{Float64} = [0., 0.]: resultant stress on ice floe   from granular interactions [Pa].\nocean_stress::Vector{Float64} = [0., 0.]: resultant stress on ice floe from   ocean drag [Pa].\natmosphere_stress::Vector{Float64} = [0., 0.]: resultant stress on ice floe   from atmosphere drag [Pa].\n\nExamples\n\nThe most basic example adds a new ice floe to the simulation sim, with a  center at [1., 2.], a radius of 1. meter, and a thickness of 0.5  meter:\n\nSeaIce.addIceFloeCylindrical!(sim, [1., 2.], 1., .5)\n\nThe following example will create a ice floe with tensile strength (cohesion), and a velocity of 0.5 m/s towards -x:\n\nSeaIce.addIceFloeCylindrical!(sim, [4., 2.], 1., .5,\n                              tensile_strength = 200e3,\n                              lin_vel = [-.5, 0.])\n\nFixed ice floes are useful for creating walls or coasts, and loops are useful for creating regular arrangements:\n\nfor i=1:5\n    SeaIce.addIceFloeCylindrical!(sim, [i*2., 0., 3.], 1., .5, fixed=true)\nend\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.addOceanDrag!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.addOceanDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between ocean and all ice  floes.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.applyAtmosphereDragToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64,Float64}",
    "page": "Public API",
    "title": "SeaIce.applyAtmosphereDragToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between atmosphere and a single  ice floe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.applyAtmosphereVorticityToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64}",
    "page": "Public API",
    "title": "SeaIce.applyAtmosphereVorticityToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between atmosphere and a  single ice floe.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.applyOceanDragToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64,Float64}",
    "page": "Public API",
    "title": "SeaIce.applyOceanDragToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between ocean and a single ice  floe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.applyOceanVorticityToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64}",
    "page": "Public API",
    "title": "SeaIce.applyOceanVorticityToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between ocean and a  single ice floe.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.areaOfQuadrilateral-NTuple{4,Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.areaOfQuadrilateral",
    "category": "Method",
    "text": "Returns the area of a quadrilateral with corner coordinates a, b, c, and  d.  Corners a and c should be opposite of each other, the same must be  true for b and d.  This is true if the four corners are passed as arguments  in a \"clockwise\" or \"counter-clockwise\" manner.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.areaOfTriangle-Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.areaOfTriangle",
    "category": "Method",
    "text": "Returns the area of an triangle with corner coordinates a, b, and c.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.checkAndAddContact!-Tuple{SeaIce.Simulation,Int64,Int64}",
    "page": "Public API",
    "title": "SeaIce.checkAndAddContact!",
    "category": "Method",
    "text": "checkAndAddContact!(simulation, i, j)\n\nCheck for contact between two ice floes and register the interaction in the  simulation object.  The indexes of the two ice floes is stored in  simulation.contact_pairs as [i, j].  The overlap vector is parallel to a  straight line connecting the ice floe centers, points away from ice floe i and  towards j, and is stored in simulation.overlaps.  A zero-length vector is  written to simulation.contact_parallel_displacement.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\ni::Int: index of the first ice floe.\nj::Int: index of the second ice floe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.checkTimeParameters-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.checkTimeParameters",
    "category": "Method",
    "text": "Checks if simulation temporal parameters are of reasonable values.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.compareAtmospheres-Tuple{SeaIce.Atmosphere,SeaIce.Atmosphere}",
    "page": "Public API",
    "title": "SeaIce.compareAtmospheres",
    "category": "Method",
    "text": "compareAtmospheres(atmosphere1::atmosphere, atmosphere2::atmosphere)\n\nCompare values of two atmosphere objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.compareIceFloes-Tuple{SeaIce.IceFloeCylindrical,SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.compareIceFloes",
    "category": "Method",
    "text": "compareIceFloes(if1::IceFloeCylindrical, if2::IceFloeCylindrical)\n\nCompare values of two ice floe objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.compareOceans-Tuple{SeaIce.Ocean,SeaIce.Ocean}",
    "page": "Public API",
    "title": "SeaIce.compareOceans",
    "category": "Method",
    "text": "compareOceans(ocean1::Ocean, ocean2::Ocean)\n\nCompare values of two Ocean objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.compareSimulations-Tuple{SeaIce.Simulation,SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.compareSimulations",
    "category": "Method",
    "text": "compareSimulations(sim1::Simulation, sim2::Simulation)\n\nCompare values of two Simulation objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.conformalQuadrilateralCoordinates-NTuple{5,Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.conformalQuadrilateralCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional coordinates [x_tilde, y_tilde] of a point p  within a quadrilateral with corner coordinates A, B, C, and D. Points must be ordered in counter-clockwise order, starting from south-west  corner.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.convertIceFloeDataToArrays-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.convertIceFloeDataToArrays",
    "category": "Method",
    "text": "Gathers all ice-floe parameters from the IceFloeCylindrical type in continuous  arrays in an IceFloeArrays structure.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.createEmptyAtmosphere-Tuple{}",
    "page": "Public API",
    "title": "SeaIce.createEmptyAtmosphere",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.createEmptyOcean-Tuple{}",
    "page": "Public API",
    "title": "SeaIce.createEmptyOcean",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.createRegularAtmosphereGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.createRegularAtmosphereGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Atmosphere grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the atmosphere spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.createRegularOceanGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.createRegularOceanGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Ocean grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the ocean spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.createSimulation-Tuple{}",
    "page": "Public API",
    "title": "SeaIce.createSimulation",
    "category": "Method",
    "text": "createSimulation([id::String=\"unnamed\",\n                  time_iteration::Int=0,\n                  time::Float64=0.0,\n                  time_total::Float64=-1.,\n                  time_step::Float64=-1.,\n                  file_time_step::Float64=-1.,\n                  file_number::Int=0,\n                  ice_floes=Array{IceFloeCylindrical, 1}[],\n                  ocean::Ocean,\n                  atmosphere::Atmosphere)\n\nCreate a simulation object containing all relevant variables such as temporal  parameters, and lists of ice floes and contacts.\n\nThe parameter id is used to uniquely identify the simulation when it is  written to disk.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.disableIceFloe!-Tuple{SeaIce.Simulation,Int64}",
    "page": "Public API",
    "title": "SeaIce.disableIceFloe!",
    "category": "Method",
    "text": "Disable ice floe with index i in the simulation object.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.disableOutputFiles!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.disableOutputFiles!",
    "category": "Method",
    "text": "Disables the write of output files to disk during a simulation.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findCellContainingPoint-Tuple{Any,Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.findCellContainingPoint",
    "category": "Method",
    "text": "findCellContainingPoint(grid, point[, method])\n\nReturns the i, j index of the grid cell containing the point. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns the coordinates of the cell.  If no match is  found the function returns (0,0).\n\nArguments\n\ngrid::Any: grid object containing ocean or atmosphere data.\npoint::Vector{Float64}: two-dimensional vector of point to check.\nmethod::String: approach to use for determining if point is inside cell or    not, can be \"Conformal\" (default) or \"Areal\".\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findContacts!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.findContacts!",
    "category": "Method",
    "text": "findContacts!(simulation[, method])\n\nTop-level function to perform an inter-ice floe contact search, based on ice  floe linear positions and contact radii.\n\nThe simplest contact search algorithm (method=\"all to all\") is the most  computationally expensive (O(n^2)).  The method \"ocean grid\" bins the ice floes  into their corresponding cells on the ocean grid and searches for contacts only  within the vicinity.  When this method is applied, it is assumed that the  contact_radius values of the ice floes are smaller than half the cell size.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\nmethod::String: the contact-search method to apply.  Valid options are \"all    to all\" and \"ocean grid\".\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findContactsAllToAll!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.findContactsAllToAll!",
    "category": "Method",
    "text": "findContactsAllToAll!(simulation)\n\nPerform an O(n^2) all-to-all contact search between all ice floes in the  simulation object.  Contacts between fixed ice floes are ignored.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findContactsInGrid!-Tuple{SeaIce.Simulation,Any}",
    "page": "Public API",
    "title": "SeaIce.findContactsInGrid!",
    "category": "Method",
    "text": "findContactsInGrid!(simulation)\n\nPerform an O(n*log(n)) cell-based contact search between all ice floes in the  simulation object.  Contacts between fixed or disabled ice floes are ignored.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findEmptyPositionInGridCell-Tuple{SeaIce.Simulation,Any,Int64,Int64,Float64}",
    "page": "Public API",
    "title": "SeaIce.findEmptyPositionInGridCell",
    "category": "Method",
    "text": "Attempt locate an empty spot for an ice floe with radius r with center  coordinates in a specified grid cell (i, j) without overlapping any other  ice floes in that cell or the neighboring cells.  This function will stop  attempting after n_iter iterations, each with randomly generated positions.\n\nThis function assumes that existing ice floes have been binned according to the  grid (e.g., using sortIceFloesInGrid()).\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findLargestIceFloeStiffness-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.findLargestIceFloeStiffness",
    "category": "Method",
    "text": "Finds the largest elastic stiffness of all ice floes in a simulation.  Used to  determine the optimal time step length.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.findSmallestIceFloeMass-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.findSmallestIceFloeMass",
    "category": "Method",
    "text": "Finds the smallest mass of all ice floes in a simulation.  Used to determine  the optimal time step length.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.getCellCenterCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Public API",
    "title": "SeaIce.getCellCenterCoordinates",
    "category": "Method",
    "text": "getCellCenterCoordinates(grid, i, j)\n\nReturns grid center coordinates (h-point).\n\nArguments\n\nxh::Array{Float64, 2}: nominal longitude of h-points [degrees_E]\nyh::Array{Float64, 2}: nominal latitude of h-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.getCellCornerCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Public API",
    "title": "SeaIce.getCellCornerCoordinates",
    "category": "Method",
    "text": "getCellCornerCoordinates(xq, yq, i, j)\n\nReturns grid corner coordinates in the following order (south-west corner,  south-east corner, north-east corner, north-west corner).\n\nArguments\n\nxq::Array{Float64, 2}: nominal longitude of q-points [degrees_E]\nyq::Array{Float64, 2}: nominal latitude of q-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.getNonDimensionalCellCoordinates-Tuple{Any,Int64,Int64,Array{Float64,1}}",
    "page": "Public API",
    "title": "SeaIce.getNonDimensionalCellCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional conformal mapped coordinates for point point in  cell i,j, based off the coordinates in the grid.\n\nThis function is a wrapper for getCellCornerCoordinates() and  conformalQuadrilateralCoordinates().\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.harmonicMean-Tuple{Number,Number}",
    "page": "Public API",
    "title": "SeaIce.harmonicMean",
    "category": "Method",
    "text": "harmonicMean(a, b)\n\nReturns the harmonic mean of two numbers a::Number and b::Number.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeCircumreference-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeCircumreference",
    "category": "Method",
    "text": "Returns the circumreference of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeHorizontalSurfaceArea-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeHorizontalSurfaceArea",
    "category": "Method",
    "text": "Returns the top or bottom (horizontal) surface area of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeKineticRotationalEnergy-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeKineticRotationalEnergy",
    "category": "Method",
    "text": "Returns the rotational kinetic energy of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeKineticTranslationalEnergy-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeKineticTranslationalEnergy",
    "category": "Method",
    "text": "Returns the translational kinetic energy of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeMass-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeMass",
    "category": "Method",
    "text": "Returns the mass of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeMomentOfInertia-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeMomentOfInertia",
    "category": "Method",
    "text": "Returns the moment of inertia of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeSideSurfaceArea-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeSideSurfaceArea",
    "category": "Method",
    "text": "Returns the surface area of the ice-floe sides\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.iceFloeVolume-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.iceFloeVolume",
    "category": "Method",
    "text": "Returns the volume of the ice floe\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.incrementCurrentTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Public API",
    "title": "SeaIce.incrementCurrentTime!",
    "category": "Method",
    "text": "incrementCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interIceFloePositionVector-Tuple{SeaIce.Simulation,Int64,Int64}",
    "page": "Public API",
    "title": "SeaIce.interIceFloePositionVector",
    "category": "Method",
    "text": "interIceFloePositionVector(simulation, i, j)\n\nReturns a vector pointing from ice floe i to ice floe j in the  simulation.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\ni::Int: index of the first ice floe.\nj::Int: index of the second ice floe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interact!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.interact!",
    "category": "Method",
    "text": "interact!(simulation::Simulation)\n\nResolve mechanical interaction between all particle pairs.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interactIceFloes!-Tuple{SeaIce.Simulation,Int64,Int64,Int64}",
    "page": "Public API",
    "title": "SeaIce.interactIceFloes!",
    "category": "Method",
    "text": "interactIceFloes!(simulation::Simulation, i::Int, j::Int, ic::Int)\n\nResolve an grain-to-grain interaction using a prescibed contact law.  This  function adds the compressive force of the interaction to the ice floe  pressure field of mean compressive stress on the ice floe sides.\n\nThe function uses the macroscopic contact-stiffness parameterization based on  Young's modulus and Poisson's ratio if Young's modulus is a positive value.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interpolateAtmosphereState-Tuple{SeaIce.Atmosphere,Float64}",
    "page": "Public API",
    "title": "SeaIce.interpolateAtmosphereState",
    "category": "Method",
    "text": "Atmosphere data is containted in Atmosphere type at discrete times  (Atmosphere.time).  This function performs linear interpolation between time  steps to get the approximate atmosphere state at any point in time.  If the  Atmosphere data set only contains a single time step, values from that time  are returned.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interpolateAtmosphereVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Public API",
    "title": "SeaIce.interpolateAtmosphereVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interpolateOceanState-Tuple{SeaIce.Ocean,Float64}",
    "page": "Public API",
    "title": "SeaIce.interpolateOceanState",
    "category": "Method",
    "text": "Ocean data is containted in Ocean type at discrete times (Ocean.time).  This  function performs linear interpolation between time steps to get the approximate  ocean state at any point in time.  If the Ocean data set only contains a  single time step, values from that time are returned.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.interpolateOceanVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Public API",
    "title": "SeaIce.interpolateOceanVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.isPointInCell",
    "page": "Public API",
    "title": "SeaIce.isPointInCell",
    "category": "Function",
    "text": "Check if a 2d point is contained inside a cell from the supplied grid. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns true or false.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.plotIceFloeSizeDistribution-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.plotIceFloeSizeDistribution",
    "category": "Method",
    "text": "plotIceFloeSizeDistribution(simulation, [filename_postfix], [nbins],\n                            [size_type], [figsize], [filetype])\n\nPlot the ice-floe size distribution as a histogram and save it to the disk.  The  plot is saved accoring to the simulation id, the optional filename_postfix  string, and the filetype, and is written to the current folder.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\nfilename_postfix::String: optional string for the output filename.\nnbins::Int: number of bins in the histogram (default = 12).\nsize_type::String: specify whether to use the contact or areal radius    for the ice-floe size.  The default is contact.\nfigsize::Tuple: the output figure size in inches (default = (6,4).\nfiletype::String: the output file type (default = \"png\").\nverbose::String: show output file as info message in stdout (default =    true).\nskip_fixed::Bool: ommit ice floes that are fixed in space from the size    distribution (default = true).\nlogy::Bool: plot y-axis in log scale.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.printIceFloeInfo-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Public API",
    "title": "SeaIce.printIceFloeInfo",
    "category": "Method",
    "text": "printIceFloeInfo(icefloe::IceFloeCylindrical)\n\nPrints the contents of an ice floe to stdout in a formatted style.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.printMemoryUsage-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.printMemoryUsage",
    "category": "Method",
    "text": "printMemoryUsage(sim::Simulation)\n\nShows the memory footprint of the simulation object.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.randpower",
    "page": "Public API",
    "title": "SeaIce.randpower",
    "category": "Function",
    "text": "randpower([nvals], [distribution_power], [min_val], [max_val])\n\nReturns one or more random numbers from a power-law probability distribution.\n\nArguments\n\ndims::Any: the dimensions of random values (default = 1)\ndistribution_power::Number: the distribution power (default = 1.)\nmin_val::Number: the lower bound of the distribution range (default = 0.)\nmax_val::Number: the upper bound of the distribution range (default = 1.)\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.readOceanGridNetCDF-Tuple{String}",
    "page": "Public API",
    "title": "SeaIce.readOceanGridNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean supergrid information generated by MOM6 (e.g.   ocean_hrid.nc) from disk and return as Ocean data structure.  This file is  located in the simulation INPUT/ subdirectory.\n\nReturns\n\nxh::Array{Float64, 2}: Longitude for cell centers [deg]\nyh::Array{Float64, 2}: Latitude for cell centers [deg]\nxq::Array{Float64, 2}: Longitude for cell corners [deg]\nyq::Array{Float64, 2}: Latitude for cell corners [deg]\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.readOceanNetCDF-Tuple{String,String}",
    "page": "Public API",
    "title": "SeaIce.readOceanNetCDF",
    "category": "Method",
    "text": "Read ocean NetCDF files generated by MOM6 from disk and return as Ocean data  structure.\n\nArguments\n\nvelocity_file::String: Path to NetCDF file containing ocean velocities,    etc., (e.g. prog__####_###.nc).\ngrid_file::String: Path to NetCDF file containing ocean super-grid    information (typically INPUT/ocean_hgrid.nc).\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.readOceanStateNetCDF-Tuple{String}",
    "page": "Public API",
    "title": "SeaIce.readOceanStateNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean state generated by MOM6 (e.g.  prog__####_###.nc  or ########.ocean_month.nc) from disk and return time stamps, velocity fields,  layer thicknesses, interface heights, and vertical coordinates.\n\nReturns\n\ntime::Vector{Float64}: Time [s]\nu::Array{Float64, 2}: Cell corner zonal velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nv::Array{Float64, 2}: Cell corner meridional velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nh::Array{Float64, 2}: layer thickness [m], dimensions correspond to    placement in [xh, yh, zl, time]\ne::Array{Float64, 2}: interface height relative to mean sea level [m],     dimensions correspond to placement in [xh, yh, zi, time]\nzl::Vector{Float64}: layer target potential density [kg m^-3]\nzi::Vector{Float64}: interface target potential density [kg m^-3]\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.readSimulation",
    "page": "Public API",
    "title": "SeaIce.readSimulation",
    "category": "Function",
    "text": "readSimulation(filename::String=\"\";\n               verbose::Bool=true)\n\nRead all content from Simulation from disk in JDL format.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.readSimulationStatus-Tuple{String}",
    "page": "Public API",
    "title": "SeaIce.readSimulationStatus",
    "category": "Method",
    "text": "readSimulationStatus(filename::String;\n                     folder::String=\".\",\n                     verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.removeSimulationFiles-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.removeSimulationFiles",
    "category": "Method",
    "text": "removeSimulationFiles(simulation[, folder])\n\nRemove all simulation output files from the specified folder.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.reportSimulationTimeToStdout-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.reportSimulationTimeToStdout",
    "category": "Method",
    "text": "Prints the current simulation time and total time to standard out\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.run!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.run!",
    "category": "Method",
    "text": "run!(simulation[,\n     verbose::Bool = true,\n     status_interval = 100.,\n     show_file_output = true,\n     single_step = false,\n     temporal_integration_method = \"Three-term Taylor\"],\n     write_jld = false)\n\nRun the simulation through time until simulation.time equals or exceeds  simulatim.time_total.  This function requires that all ice floes are added to  the simulation and that the length of the computational time step is adjusted  accordingly.\n\nThe function will search for contacts, determine the force balance on each ice  floe, and integrate all kinematic degrees of freedom accordingly.  The temporal  integration is explicit and of length simulation.time_step.  This function  will write VTK files to disk in the intervals simulation.file_time_step by the  function writeVTK.  If this value is negative, no output files will be written  to disk.\n\nArguments\n\nsimulation::Simulation: the simulation to run (object is modified)\nverbose::Bool=true: show verbose information during the time loop\nstatus_interval::Bool=true: show verbose information during the time loop\nshow_file_output::Bool=true: report to stdout when output file is written\nsingle_step::Bool=false: run simulation for a single time step only.  If    this causes simulation.time to exceed simulation.time_total, the latter    is increased accordingly.\ntemporal_integration_method::String=\"Three-term Taylor\": type of integration    method to use.  See updateIceFloeKinematics for details.\nwrite_jld::Bool=false: write simulation state to disk as JLD files (see    SeaIce.writeSimulation(...) whenever saving VTK output.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.setCurrentTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Public API",
    "title": "SeaIce.setCurrentTime!",
    "category": "Method",
    "text": "setCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.setOutputFileInterval!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Public API",
    "title": "SeaIce.setOutputFileInterval!",
    "category": "Method",
    "text": "setOutputFileInterval!(simulation::Simulation, t::Float64)\n\nSets the simulation-time interval between output files are written to disk.  If  this value is zero or negative, no output files will be written.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.setTimeStep!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.setTimeStep!",
    "category": "Method",
    "text": "Find the computational time step length suitable given the grain radii, contact stiffnesses, and grain density. Uses the scheme by Radjaii et al. 2011.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.setTotalTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Public API",
    "title": "SeaIce.setTotalTime!",
    "category": "Method",
    "text": "setTotalTime!(simulation::Simulation, t::Float64)\n\nSets the total simulation time of the simulation object to t, with parameter  value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.sortIceFloesInGrid!-Tuple{SeaIce.Simulation,Any}",
    "page": "Public API",
    "title": "SeaIce.sortIceFloesInGrid!",
    "category": "Method",
    "text": "Find ice-floe positions in grid, based on their center positions.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.status",
    "page": "Public API",
    "title": "SeaIce.status",
    "category": "Function",
    "text": "Shows the status of all simulations with output files written under the  specified folder, which is the current working directory by default.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.totalIceFloeKineticRotationalEnergy-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.totalIceFloeKineticRotationalEnergy",
    "category": "Method",
    "text": "Returns the sum of rotational kinetic energies of all ice floes in a simulation\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.totalIceFloeKineticTranslationalEnergy-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.totalIceFloeKineticTranslationalEnergy",
    "category": "Method",
    "text": "Returns the sum of translational kinetic energies of all ice floes in a  simulation\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.updateIceFloeKinematics!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.updateIceFloeKinematics!",
    "category": "Method",
    "text": "updateIceFloeKinematics!(simulation::Simulation[,\n                         method::String = \"Three-term Taylor\"])\n\nUpdate the ice floe kinematic parameters using a temporal integration scheme, the current force and torque balance, and gravitational acceleration.\n\nArguments\n\nsimulation::Simulation: update the ice floe positions in this object    according to temporal integration of length simulation.time_step.\nmethod::String = \"Three-term Taylor\": the integration method to use.     Available methods are \"Two-term Taylor\" and \"Three-term Taylor\".  The    three-term Taylor expansion is slightly more computationally expensive than    the two-term Taylor expansion, but offers an order-of-magnitude increase in    precision of ice floe positions.  The two-term expansion can obtain similar    precision if the time step is 1/10 the length.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.updateIceFloeKinematicsThreeTermTaylor!-Tuple{SeaIce.IceFloeCylindrical,SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.updateIceFloeKinematicsThreeTermTaylor!",
    "category": "Method",
    "text": "Use a three-term Taylor expansion for integrating the kinematic degrees of  freedom for an icefloe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.updateIceFloeKinematicsTwoTermTaylor!-Tuple{SeaIce.IceFloeCylindrical,SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.updateIceFloeKinematicsTwoTermTaylor!",
    "category": "Method",
    "text": "Use a two-term Taylor expansion for integrating the kinematic degrees of freedom  for an icefloe.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeIceFloeInteractionVTK-Tuple{SeaIce.Simulation,String}",
    "page": "Public API",
    "title": "SeaIce.writeIceFloeInteractionVTK",
    "category": "Method",
    "text": "writeIceFloeInteractionVTK(simulation::Simulation,\n                           filename::String;\n                           verbose::Bool=false)\n\nSaves ice-floe interactions to .vtp files for visualization with VTK, for  example in Paraview.  Convert Cell Data to Point Data and use with Tube filter.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeIceFloeVTK-Tuple{SeaIce.Simulation,String}",
    "page": "Public API",
    "title": "SeaIce.writeIceFloeVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ice floes in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.  This function is called by  writeVTK().\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeParaviewStateFile-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.writeParaviewStateFile",
    "category": "Method",
    "text": "Create a Paraview State File (.pvsm) for the simulation, which reads simulation  output VTK files and applies appropriate glyph filters to the data.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeSimulation-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.writeSimulation",
    "category": "Method",
    "text": "writeSimulation(simulation::Simulation;\n                     filename::String=\"\",\n                     folder::String=\".\",\n                     verbose::Bool=true)\n\nWrite all content from Simulation to disk in JDL format.  If the filename  parameter is not specified, it will be saved to a subdirectory under the current  directory named after the simulation identifier simulation.id.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeSimulationStatus-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.writeSimulationStatus",
    "category": "Method",
    "text": "writeSimulationStatus(simulation::Simulation;\n                      folder::String=\".\",\n                      verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.writeVTK-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.writeVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ice floes in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.\n\nIf the simulation contains an Ocean data structure, it's contents will be  written to separate .vtu files.  This can be disabled by setting the argument  ocean=false.  The same is true for the atmosphere.\n\nThe VTK files will be saved in a subfolder named after the simulation.\n\n\n\n"
},

{
    "location": "lib/public.html#SeaIce.zeroForcesAndTorques!-Tuple{SeaIce.Simulation}",
    "page": "Public API",
    "title": "SeaIce.zeroForcesAndTorques!",
    "category": "Method",
    "text": "Sets the force and torque values of all ice floes to zero.\n\n\n\n"
},

{
    "location": "lib/public.html#Public-Interface-1",
    "page": "Public API",
    "title": "Public Interface",
    "category": "section",
    "text": "Modules = [SeaIce]\nPublic = true\nPrivate = false"
},

{
    "location": "lib/internals.html#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals.html#Package-internal-documentation-1",
    "page": "Internals",
    "title": "Package-internal documentation",
    "category": "section",
    "text": "This page lists all the documented internals of the SeaIce module."
},

{
    "location": "lib/internals.html#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals.html#SeaIce.bilinearInterpolation!-Tuple{Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Int64,Int64}",
    "page": "Internals",
    "title": "SeaIce.bilinearInterpolation!",
    "category": "Method",
    "text": "bilinearInterpolation(field, x_tilde, y_tilde, i, j, k, it)\n\nUse bilinear interpolation to interpolate a staggered grid to an arbitrary  position in a cell.  Assumes south-west convention, i.e. (i,j) is located at the  south-west (-x, -y)-facing corner.\n\nArguments\n\nfield::Array{Float64, 4}: a scalar field to interpolate from\nx_tilde::Float64: x point position [0;1]\ny_tilde::Float64: y point position [0;1]\ni::Int: i-index of cell containing point\nj::Int: j-index of scalar field to interpolate from\nit::Int: time step from scalar field to interpolate from\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.copyGridSortingInfo!-Tuple{SeaIce.Ocean,SeaIce.Atmosphere,Array{SeaIce.IceFloeCylindrical,1}}",
    "page": "Internals",
    "title": "SeaIce.copyGridSortingInfo!",
    "category": "Method",
    "text": "Copy ice floe related information from ocean to atmosphere grid.  This is useful  when the two grids are of identical geometry, meaning only only one sorting  phase is necessary.\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.curl",
    "page": "Internals",
    "title": "SeaIce.curl",
    "category": "Function",
    "text": "curl(grid, x_tilde, y_tilde, i, j, k, it)\n\nUse bilinear interpolation to interpolate curl value for a staggered velocity  grid to an arbitrary position in a cell.  Assumes south-west convention, i.e.   (i,j) is located at the south-west (-x, -y)-facing corner.\n\nArguments\n\ngrid::Any: grid for which to determine curl\nx_tilde::Float64: x point position [0;1]\ny_tilde::Float64: y point position [0;1]\ni::Int: i-index of cell containing point\nj::Int: j-index of scalar field to interpolate from\nit::Int: time step from scalar field to interpolate from\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.findOverlap-Tuple{SeaIce.Simulation,Int64,Int64,Array{Float64,1}}",
    "page": "Internals",
    "title": "SeaIce.findOverlap",
    "category": "Method",
    "text": "position_ij is the inter-grain position vector, and can be found with interIceFloePositionVector().\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.writeGridVTK-Tuple{Any,String}",
    "page": "Internals",
    "title": "SeaIce.writeGridVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ocean data in the simulation in a  structured grid (file type .vts).  These files can be read by ParaView and can  be visualized by applying a Glyph filter.  This function is called by  writeVTK().\n\n\n\n"
},

{
    "location": "lib/internals.html#Internal-Interface-1",
    "page": "Internals",
    "title": "Internal Interface",
    "category": "section",
    "text": "Modules = [SeaIce]\nPublic = false\nPrivate = true"
},

{
    "location": "lib/internals.html#",
    "page": "Package-internal documentation",
    "title": "Package-internal documentation",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals.html#Package-internal-documentation-1",
    "page": "Package-internal documentation",
    "title": "Package-internal documentation",
    "category": "section",
    "text": "This page lists all the documented internals of the SeaIce module."
},

{
    "location": "lib/internals.html#Index-1",
    "page": "Package-internal documentation",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals.html#SeaIce.bilinearInterpolation!-Tuple{Array{Float64,1},Array{Float64,2},Array{Float64,2},Float64,Float64,Int64,Int64}",
    "page": "Package-internal documentation",
    "title": "SeaIce.bilinearInterpolation!",
    "category": "Method",
    "text": "bilinearInterpolation(field, x_tilde, y_tilde, i, j, k, it)\n\nUse bilinear interpolation to interpolate a staggered grid to an arbitrary  position in a cell.  Assumes south-west convention, i.e. (i,j) is located at the  south-west (-x, -y)-facing corner.\n\nArguments\n\nfield::Array{Float64, 4}: a scalar field to interpolate from\nx_tilde::Float64: x point position [0;1]\ny_tilde::Float64: y point position [0;1]\ni::Int: i-index of cell containing point\nj::Int: j-index of scalar field to interpolate from\nit::Int: time step from scalar field to interpolate from\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.copyGridSortingInfo!-Tuple{SeaIce.Ocean,SeaIce.Atmosphere,Array{SeaIce.IceFloeCylindrical,1}}",
    "page": "Package-internal documentation",
    "title": "SeaIce.copyGridSortingInfo!",
    "category": "Method",
    "text": "Copy ice floe related information from ocean to atmosphere grid.  This is useful  when the two grids are of identical geometry, meaning only only one sorting  phase is necessary.\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.curl",
    "page": "Package-internal documentation",
    "title": "SeaIce.curl",
    "category": "Function",
    "text": "curl(grid, x_tilde, y_tilde, i, j, k, it)\n\nUse bilinear interpolation to interpolate curl value for a staggered velocity  grid to an arbitrary position in a cell.  Assumes south-west convention, i.e.   (i,j) is located at the south-west (-x, -y)-facing corner.\n\nArguments\n\ngrid::Any: grid for which to determine curl\nx_tilde::Float64: x point position [0;1]\ny_tilde::Float64: y point position [0;1]\ni::Int: i-index of cell containing point\nj::Int: j-index of scalar field to interpolate from\nit::Int: time step from scalar field to interpolate from\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.findOverlap-Tuple{SeaIce.Simulation,Int64,Int64,Array{Float64,1}}",
    "page": "Package-internal documentation",
    "title": "SeaIce.findOverlap",
    "category": "Method",
    "text": "position_ij is the inter-grain position vector, and can be found with interIceFloePositionVector().\n\n\n\n"
},

{
    "location": "lib/internals.html#SeaIce.writeGridVTK-Tuple{Any,String}",
    "page": "Package-internal documentation",
    "title": "SeaIce.writeGridVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ocean data in the simulation in a  structured grid (file type .vts).  These files can be read by ParaView and can  be visualized by applying a Glyph filter.  This function is called by  writeVTK().\n\n\n\n"
},

{
    "location": "lib/internals.html#Internal-Interface-1",
    "page": "Package-internal documentation",
    "title": "Internal Interface",
    "category": "section",
    "text": "Modules = [SeaIce]\nPublic = false\nPrivate = true"
},

]}
