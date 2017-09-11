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
    "text": "Package for particle-based simulation of sea-ice dynamics"
},

{
    "location": "index.html#Package-features-1",
    "page": "Home",
    "title": "Package features",
    "category": "section",
    "text": "Flexible and computationally efficient 2d implementation of the discrete element method.  The particles represent sea-ice floes, which can be forced by ocean and atmospheric velocity fields.  The ice floes can interact through elasto-viscous-frictional contact rheologies and obtain time-dependent tensile strength."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"installation.md\",\n    \"module.md\"\n]\nDepth = 1"
},

{
    "location": "index.html#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "SeaIce.jl can be installed directly from the Julia shell by:Pkg.clone(\"git://github.com/anders-dc/SeaIce.jl.git\")This will install the contents of this repository in the folder  ~/.julia/v$(JULIA_VERSION)/SeaIce, and install the packages specified as  requirements.  The package JLD  is used for model restarts and is recommended but not required, and thus is not  automatically installed.Import the package contents into the current Julia session or script with:import SeaIceThis will import all functions and data types in the SeaIce namespace.  You  can run the package tests, which are contained in the test/ directory,  with the following command:Pkg.test(\"SeaIce\")The package can be updated from this repository using:Pkg.update(\"SeaIce\")"
},

{
    "location": "module.html#",
    "page": "Modules, constants, types, functions, and macros",
    "title": "Modules, constants, types, functions, and macros",
    "category": "page",
    "text": ""
},

{
    "location": "module.html#SeaIce",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce",
    "category": "Module",
    "text": "SeaIce.jl\n\nOffline sea-ice dynamics simulator by Anders Damsgaard, www.adamsgaard.dk\n\n\n\n"
},

{
    "location": "module.html#SeaIce.addAtmosphereDrag!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.addAtmosphereDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between atmosphere and all  ice floes.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.addIceFloe!",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.addIceFloe!",
    "category": "Function",
    "text": "addIceFloe!(simulation::Simulation,\n            icefloe::IceFloeCylindrical,\n            verbose::Bool = False)\n\nAdd an icefloe to the simulation object.  If verbose is true, a short  confirmation message will be printed to stdout.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.addIceFloeCylindrical!-Tuple{SeaIce.Simulation,Array{Float64,1},Float64,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.addIceFloeCylindrical!",
    "category": "Method",
    "text": "Adds a grain to the simulation. Example:\n\nSeaIce.addIceFloeCylindrical([1.0, 2.0, 3.0], 1.0)\n\n\n\n"
},

{
    "location": "module.html#SeaIce.addOceanDrag!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.addOceanDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between ocean and all ice  floes.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.applyAtmosphereDragToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.applyAtmosphereDragToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between atmosphere and a single  ice floe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.applyAtmosphereVorticityToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.applyAtmosphereVorticityToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between atmosphere and a  single ice floe.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.applyOceanDragToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.applyOceanDragToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between ocean and a single ice  floe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.applyOceanVorticityToIceFloe!-Tuple{SeaIce.IceFloeCylindrical,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.applyOceanVorticityToIceFloe!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between ocean and a  single ice floe.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.areaOfQuadrilateral-NTuple{4,Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.areaOfQuadrilateral",
    "category": "Method",
    "text": "Returns the area of a quadrilateral with corner coordinates a, b, c, and  d.  Corners a and c should be opposite of each other, the same must be  true for b and d.  This is true if the four corners are passed as arguments  in a \"clockwise\" or \"counter-clockwise\" manner.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.areaOfTriangle-Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.areaOfTriangle",
    "category": "Method",
    "text": "Returns the area of an triangle with corner coordinates a, b, and c.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.checkAndAddContact!-Tuple{SeaIce.Simulation,Int64,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.checkAndAddContact!",
    "category": "Method",
    "text": "checkAndAddContact!(simulation, i, j)\n\nCheck for contact between two ice floes and register the interaction in the  simulation object.  The indexes of the two ice floes is stored in  simulation.contact_pairs as [i, j].  The overlap vector is parallel to a  straight line connecting the ice floe centers, points away from ice floe i and  towards j, and is stored in simulation.overlaps.  A zero-length vector is  written to simulation.contact_parallel_displacement.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\ni::Int: index of the first ice floe.\nj::Int: index of the second ice floe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.checkTimeParameters-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.checkTimeParameters",
    "category": "Method",
    "text": "Checks if simulation temporal parameters are of reasonable values.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.compareAtmospheres-Tuple{SeaIce.Atmosphere,SeaIce.Atmosphere}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.compareAtmospheres",
    "category": "Method",
    "text": "compareAtmospheres(atmosphere1::atmosphere, atmosphere2::atmosphere)\n\nCompare values of two atmosphere objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.compareIceFloes-Tuple{SeaIce.IceFloeCylindrical,SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.compareIceFloes",
    "category": "Method",
    "text": "compareIceFloes(if1::IceFloeCylindrical, if2::IceFloeCylindrical)\n\nCompare values of two ice floe objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.compareOceans-Tuple{SeaIce.Ocean,SeaIce.Ocean}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.compareOceans",
    "category": "Method",
    "text": "compareOceans(ocean1::Ocean, ocean2::Ocean)\n\nCompare values of two Ocean objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.compareSimulations-Tuple{SeaIce.Simulation,SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.compareSimulations",
    "category": "Method",
    "text": "compareSimulations(sim1::Simulation, sim2::Simulation)\n\nCompare values of two Simulation objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.conformalQuadrilateralCoordinates-NTuple{5,Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.conformalQuadrilateralCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional coordinates [x_tilde, y_tilde] of a point p  within a quadrilateral with corner coordinates A, B, C, and D. Points must be ordered in counter-clockwise order, starting from south-west  corner.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.convertIceFloeDataToArrays-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.convertIceFloeDataToArrays",
    "category": "Method",
    "text": "Gathers all ice-floe parameters from the IceFloeCylindrical type in continuous  arrays in an IceFloeArrays structure.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.createEmptyAtmosphere-Tuple{}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.createEmptyAtmosphere",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.createEmptyOcean-Tuple{}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.createEmptyOcean",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.createRegularAtmosphereGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.createRegularAtmosphereGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Atmosphere grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the atmosphere spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.createRegularOceanGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.createRegularOceanGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Ocean grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the ocean spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.createSimulation-Tuple{}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.createSimulation",
    "category": "Method",
    "text": "createSimulation([id::String=\"unnamed\",\n                  time_iteration::Int=0,\n                  time::Float64=0.0,\n                  time_total::Float64=-1.,\n                  time_step::Float64=-1.,\n                  file_time_step::Float64=-1.,\n                  file_number::Int=0,\n                  ice_floes=Array{IceFloeCylindrical, 1}[],\n                  ocean::Ocean,\n                  atmosphere::Atmosphere)\n\nCreate a simulation object containing all relevant variables such as temporal  parameters, and lists of ice floes and contacts.\n\nThe parameter id is used to uniquely identify the simulation when it is  written to disk.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.disableIceFloe!-Tuple{SeaIce.Simulation,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.disableIceFloe!",
    "category": "Method",
    "text": "Disable ice floe with index i in the simulation object.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.disableOutputFiles!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.disableOutputFiles!",
    "category": "Method",
    "text": "Disables the write of output files to disk during a simulation.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findCellContainingPoint-Tuple{Any,Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findCellContainingPoint",
    "category": "Method",
    "text": "findCellContainingPoint(grid, point[, method])\n\nReturns the i, j index of the grid cell containing the point. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns the coordinates of the cell.  If no match is  found the function returns (0,0).\n\nArguments\n\ngrid::Any: grid object containing ocean or atmosphere data.\npoint::Vector{Float64}: two-dimensional vector of point to check.\nmethod::String: approach to use for determining if point is inside cell or    not, can be \"Conformal\" (default) or \"Areal\".\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findContacts!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findContacts!",
    "category": "Method",
    "text": "findContacts!(simulation[, method])\n\nTop-level function to perform an inter-ice floe contact search, based on ice  floe linear positions and contact radii.\n\nThe simplest contact search algorithm (method=\"all to all\") is the most  computationally expensive (O(n^2)).  The method \"ocean grid\" bins the ice floes  into their corresponding cells on the ocean grid and searches for contacts only  within the vicinity.  When this method is applied, it is assumed that the  contact_radius values of the ice floes are smaller than half the cell size.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\nmethod::String: the contact-search method to apply.  Valid options are \"all    to all\" and \"ocean grid\".\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findContactsAllToAll!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findContactsAllToAll!",
    "category": "Method",
    "text": "findContactsAllToAll!(simulation)\n\nPerform an O(n^2) all-to-all contact search between all ice floes in the  simulation object.  Contacts between fixed ice floes are ignored.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findContactsInGrid!-Tuple{SeaIce.Simulation,Any}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findContactsInGrid!",
    "category": "Method",
    "text": "findContactsInGrid!(simulation)\n\nPerform an O(n*log(n)) cell-based contact search between all ice floes in the  simulation object.  Contacts between fixed or disabled ice floes are ignored.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findEmptyPositionInGridCell-Tuple{SeaIce.Simulation,Any,Int64,Int64,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findEmptyPositionInGridCell",
    "category": "Method",
    "text": "Attempt locate an empty spot for an ice floe with radius r with center  coordinates in a specified grid cell (i, j) without overlapping any other  ice floes in that cell or the neighboring cells.  This function will stop  attempting after n_iter iterations, each with randomly generated positions.\n\nThis function assumes that existing ice floes have been binned according to the  grid (e.g., using sortIceFloesInGrid()).\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findLargestIceFloeStiffness-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findLargestIceFloeStiffness",
    "category": "Method",
    "text": "Finds the largest elastic stiffness of all ice floes in a simulation.  Used to  determine the optimal time step length.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.findSmallestIceFloeMass-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.findSmallestIceFloeMass",
    "category": "Method",
    "text": "Finds the smallest mass of all ice floes in a simulation.  Used to determine  the optimal time step length.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.getCellCenterCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.getCellCenterCoordinates",
    "category": "Method",
    "text": "getCellCenterCoordinates(grid, i, j)\n\nReturns grid center coordinates (h-point).\n\nArguments\n\nxh::Array{Float64, 2}: nominal longitude of h-points [degrees_E]\nyh::Array{Float64, 2}: nominal latitude of h-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.getCellCornerCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.getCellCornerCoordinates",
    "category": "Method",
    "text": "getCellCornerCoordinates(xq, yq, i, j)\n\nReturns grid corner coordinates in the following order (south-west corner,  south-east corner, north-east corner, north-west corner).\n\nArguments\n\nxq::Array{Float64, 2}: nominal longitude of q-points [degrees_E]\nyq::Array{Float64, 2}: nominal latitude of q-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.getNonDimensionalCellCoordinates-Tuple{Any,Int64,Int64,Array{Float64,1}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.getNonDimensionalCellCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional conformal mapped coordinates for point point in  cell i,j, based off the coordinates in the grid.\n\nThis function is a wrapper for getCellCornerCoordinates() and  conformalQuadrilateralCoordinates().\n\n\n\n"
},

{
    "location": "module.html#SeaIce.harmonicMean-Tuple{Number,Number}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.harmonicMean",
    "category": "Method",
    "text": "harmonicMean(a, b)\n\nReturns the harmonic mean of two numbers a::Number and b::Number.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeCircumreference-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeCircumreference",
    "category": "Method",
    "text": "Returns the circumreference of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeHorizontalSurfaceArea-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeHorizontalSurfaceArea",
    "category": "Method",
    "text": "Returns the top or bottom (horizontal) surface area of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeKineticRotationalEnergy-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeKineticRotationalEnergy",
    "category": "Method",
    "text": "Returns the rotational kinetic energy of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeKineticTranslationalEnergy-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeKineticTranslationalEnergy",
    "category": "Method",
    "text": "Returns the translational kinetic energy of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeMass-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeMass",
    "category": "Method",
    "text": "Returns the mass of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeMomentOfInertia-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeMomentOfInertia",
    "category": "Method",
    "text": "Returns the moment of inertia of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeSideSurfaceArea-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeSideSurfaceArea",
    "category": "Method",
    "text": "Returns the surface area of the ice-floe sides\n\n\n\n"
},

{
    "location": "module.html#SeaIce.iceFloeVolume-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.iceFloeVolume",
    "category": "Method",
    "text": "Returns the volume of the ice floe\n\n\n\n"
},

{
    "location": "module.html#SeaIce.incrementCurrentTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.incrementCurrentTime!",
    "category": "Method",
    "text": "incrementCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interIceFloePositionVector-Tuple{SeaIce.Simulation,Int64,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interIceFloePositionVector",
    "category": "Method",
    "text": "interIceFloePositionVector(simulation, i, j)\n\nReturns a vector pointing from ice floe i to ice floe j in the  simulation.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\ni::Int: index of the first ice floe.\nj::Int: index of the second ice floe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interact!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interact!",
    "category": "Method",
    "text": "interact!(simulation::Simulation)\n\nResolve mechanical interaction between all particle pairs.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interactIceFloes!-Tuple{SeaIce.Simulation,Int64,Int64,Int64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interactIceFloes!",
    "category": "Method",
    "text": "interactIceFloes!(simulation::Simulation, i::Int, j::Int, ic::Int)\n\nResolve an grain-to-grain interaction using a prescibed contact law.  This  function adds the compressive force of the interaction to the ice floe  pressure field of mean compressive stress on the ice floe sides.\n\nThe function uses the macroscopic contact-stiffness parameterization based on  Young's modulus and Poisson's ratio if Young's modulus is a positive value.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interpolateAtmosphereState-Tuple{SeaIce.Atmosphere,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interpolateAtmosphereState",
    "category": "Method",
    "text": "Atmosphere data is containted in Atmosphere type at discrete times  (Atmosphere.time).  This function performs linear interpolation between time  steps to get the approximate atmosphere state at any point in time.  If the  Atmosphere data set only contains a single time step, values from that time  are returned.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interpolateAtmosphereVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interpolateAtmosphereVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interpolateOceanState-Tuple{SeaIce.Ocean,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interpolateOceanState",
    "category": "Method",
    "text": "Ocean data is containted in Ocean type at discrete times (Ocean.time).  This  function performs linear interpolation between time steps to get the approximate  ocean state at any point in time.  If the Ocean data set only contains a  single time step, values from that time are returned.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.interpolateOceanVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.interpolateOceanVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.isPointInCell",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.isPointInCell",
    "category": "Function",
    "text": "Check if a 2d point is contained inside a cell from the supplied grid. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns true or false.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.plotIceFloeSizeDistribution-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.plotIceFloeSizeDistribution",
    "category": "Method",
    "text": "plotIceFloeSizeDistribution(simulation, [filename_postfix], [nbins],\n                            [size_type], [figsize], [filetype])\n\nPlot the ice-floe size distribution as a histogram and save it to the disk.  The  plot is saved accoring to the simulation id, the optional filename_postfix  string, and the filetype, and is written to the current folder.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the ice floes.\nfilename_postfix::String: optional string for the output filename.\nnbins::Int: number of bins in the histogram (default = 12).\nsize_type::String: specify whether to use the contact or areal radius    for the ice-floe size.  The default is contact.\nfigsize::Tuple: the output figure size in inches (default = (6,4).\nfiletype::String: the output file type (default = \"png\").\nverbose::String: show output file as info message in stdout (default =    true).\nskip_fixed::Bool: ommit ice floes that are fixed in space from the size    distribution (default = true).\nlogy::Bool: plot y-axis in log scale.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.printIceFloeInfo-Tuple{SeaIce.IceFloeCylindrical}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.printIceFloeInfo",
    "category": "Method",
    "text": "printIceFloeInfo(icefloe::IceFloeCylindrical)\n\nPrints the contents of an ice floe to stdout in a formatted style.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.printMemoryUsage-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.printMemoryUsage",
    "category": "Method",
    "text": "printMemoryUsage(sim::Simulation)\n\nShows the memory footprint of the simulation object.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.randpower",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.randpower",
    "category": "Function",
    "text": "randpower([nvals], [distribution_power], [min_val], [max_val])\n\nReturns one or more random numbers from a power-law probability distribution.\n\nArguments\n\ndims::Any: the dimensions of random values (default = 1)\ndistribution_power::Number: the distribution power (default = 1.)\nmin_val::Number: the lower bound of the distribution range (default = 0.)\nmax_val::Number: the upper bound of the distribution range (default = 1.)\n\n\n\n"
},

{
    "location": "module.html#SeaIce.readOceanGridNetCDF-Tuple{String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.readOceanGridNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean supergrid information generated by MOM6 (e.g.   ocean_hrid.nc) from disk and return as Ocean data structure.  This file is  located in the simulation INPUT/ subdirectory.\n\nReturns\n\nxh::Array{Float64, 2}: Longitude for cell centers [deg]\nyh::Array{Float64, 2}: Latitude for cell centers [deg]\nxq::Array{Float64, 2}: Longitude for cell corners [deg]\nyq::Array{Float64, 2}: Latitude for cell corners [deg]\n\n\n\n"
},

{
    "location": "module.html#SeaIce.readOceanNetCDF-Tuple{String,String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.readOceanNetCDF",
    "category": "Method",
    "text": "Read ocean NetCDF files generated by MOM6 from disk and return as Ocean data  structure.\n\nArguments\n\nvelocity_file::String: Path to NetCDF file containing ocean velocities,    etc., (e.g. prog__####_###.nc).\ngrid_file::String: Path to NetCDF file containing ocean super-grid    information (typically INPUT/ocean_hgrid.nc).\n\n\n\n"
},

{
    "location": "module.html#SeaIce.readOceanStateNetCDF-Tuple{String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.readOceanStateNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean state generated by MOM6 (e.g.  prog__####_###.nc  or ########.ocean_month.nc) from disk and return time stamps, velocity fields,  layer thicknesses, interface heights, and vertical coordinates.\n\nReturns\n\ntime::Vector{Float64}: Time [s]\nu::Array{Float64, 2}: Cell corner zonal velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nv::Array{Float64, 2}: Cell corner meridional velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nh::Array{Float64, 2}: layer thickness [m], dimensions correspond to    placement in [xh, yh, zl, time]\ne::Array{Float64, 2}: interface height relative to mean sea level [m],     dimensions correspond to placement in [xh, yh, zi, time]\nzl::Vector{Float64}: layer target potential density [kg m^-3]\nzi::Vector{Float64}: interface target potential density [kg m^-3]\n\n\n\n"
},

{
    "location": "module.html#SeaIce.readSimulation",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.readSimulation",
    "category": "Function",
    "text": "readSimulation(filename::String=\"\";\n               verbose::Bool=true)\n\nRead all content from Simulation from disk in JDL format.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.readSimulationStatus-Tuple{String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.readSimulationStatus",
    "category": "Method",
    "text": "readSimulationStatus(filename::String;\n                     folder::String=\".\",\n                     verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.removeSimulationFiles-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.removeSimulationFiles",
    "category": "Method",
    "text": "removeSimulationFiles(simulation[, folder])\n\nRemove all simulation output files from the specified folder.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.reportSimulationTimeToStdout-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.reportSimulationTimeToStdout",
    "category": "Method",
    "text": "Prints the current simulation time and total time to standard out\n\n\n\n"
},

{
    "location": "module.html#SeaIce.run!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.run!",
    "category": "Method",
    "text": "run!(simulation[,\n     verbose::Bool = true,\n     status_interval = 100.,\n     show_file_output = true,\n     single_step = false,\n     temporal_integration_method = \"Three-term Taylor\"],\n     write_jld = false)\n\nRun the simulation through time until simulation.time equals or exceeds  simulatim.time_total.  This function requires that all ice floes are added to  the simulation and that the length of the computational time step is adjusted  accordingly.\n\nThe function will search for contacts, determine the force balance on each ice  floe, and integrate all kinematic degrees of freedom accordingly.  The temporal  integration is explicit and of length simulation.time_step.  This function  will write VTK files to disk in the intervals simulation.file_time_step by the  function writeVTK.  If this value is negative, no output files will be written  to disk.\n\nArguments\n\nsimulation::Simulation: the simulation to run (object is modified)\nverbose::Bool=true: show verbose information during the time loop\nstatus_interval::Bool=true: show verbose information during the time loop\nshow_file_output::Bool=true: report to stdout when output file is written\nsingle_step::Bool=false: run simulation for a single time step only.  If    this causes simulation.time to exceed simulation.time_total, the latter    is increased accordingly.\ntemporal_integration_method::String=\"Three-term Taylor\": type of integration    method to use.  See updateIceFloeKinematics for details.\nwrite_jld::Bool=false: write simulation state to disk as JLD files (see    SeaIce.writeSimulation(...) whenever saving VTK output.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.setCurrentTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.setCurrentTime!",
    "category": "Method",
    "text": "setCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.setOutputFileInterval!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.setOutputFileInterval!",
    "category": "Method",
    "text": "setOutputFileInterval!(simulation::Simulation, t::Float64)\n\nSets the simulation-time interval between output files are written to disk.  If  this value is zero or negative, no output files will be written.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.setTimeStep!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.setTimeStep!",
    "category": "Method",
    "text": "Find the computational time step length suitable given the grain radii, contact stiffnesses, and grain density. Uses the scheme by Radjaii et al. 2011.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.setTotalTime!-Tuple{SeaIce.Simulation,Float64}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.setTotalTime!",
    "category": "Method",
    "text": "setTotalTime!(simulation::Simulation, t::Float64)\n\nSets the total simulation time of the simulation object to t, with parameter  value checks.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.sortIceFloesInGrid!-Tuple{SeaIce.Simulation,Any}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.sortIceFloesInGrid!",
    "category": "Method",
    "text": "Find ice-floe positions in grid, based on their center positions.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.status",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.status",
    "category": "Function",
    "text": "Shows the status of all simulations with output files written under the  specified folder, which is the current working directory by default.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.totalIceFloeKineticRotationalEnergy-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.totalIceFloeKineticRotationalEnergy",
    "category": "Method",
    "text": "Returns the sum of rotational kinetic energies of all ice floes in a simulation\n\n\n\n"
},

{
    "location": "module.html#SeaIce.totalIceFloeKineticTranslationalEnergy-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.totalIceFloeKineticTranslationalEnergy",
    "category": "Method",
    "text": "Returns the sum of translational kinetic energies of all ice floes in a  simulation\n\n\n\n"
},

{
    "location": "module.html#SeaIce.updateIceFloeKinematics!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.updateIceFloeKinematics!",
    "category": "Method",
    "text": "updateIceFloeKinematics!(simulation::Simulation[,\n                         method::String = \"Three-term Taylor\"])\n\nUpdate the ice floe kinematic parameters using a temporal integration scheme, the current force and torque balance, and gravitational acceleration.\n\nArguments\n\nsimulation::Simulation: update the ice floe positions in this object    according to temporal integration of length simulation.time_step.\nmethod::String = \"Three-term Taylor\": the integration method to use.     Available methods are \"Two-term Taylor\" and \"Three-term Taylor\".  The    three-term Taylor expansion is slightly more computationally expensive than    the two-term Taylor expansion, but offers an order-of-magnitude increase in    precision of ice floe positions.  The two-term expansion can obtain similar    precision if the time step is 1/10 the length.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.updateIceFloeKinematicsThreeTermTaylor!-Tuple{SeaIce.IceFloeCylindrical,SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.updateIceFloeKinematicsThreeTermTaylor!",
    "category": "Method",
    "text": "Use a three-term Taylor expansion for integrating the kinematic degrees of  freedom for an icefloe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.updateIceFloeKinematicsTwoTermTaylor!-Tuple{SeaIce.IceFloeCylindrical,SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.updateIceFloeKinematicsTwoTermTaylor!",
    "category": "Method",
    "text": "Use a two-term Taylor expansion for integrating the kinematic degrees of freedom  for an icefloe.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeIceFloeInteractionVTK-Tuple{SeaIce.Simulation,String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeIceFloeInteractionVTK",
    "category": "Method",
    "text": "writeIceFloeInteractionVTK(simulation::Simulation,\n                           filename::String;\n                           verbose::Bool=false)\n\nSaves ice-floe interactions to .vtp files for visualization with VTK, for  example in Paraview.  Convert Cell Data to Point Data and use with Tube filter.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeIceFloeVTK-Tuple{SeaIce.Simulation,String}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeIceFloeVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ice floes in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.  This function is called by  writeVTK().\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeParaviewStateFile-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeParaviewStateFile",
    "category": "Method",
    "text": "Create a Paraview State File (.pvsm) for the simulation, which reads simulation  output VTK files and applies appropriate glyph filters to the data.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeSimulation-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeSimulation",
    "category": "Method",
    "text": "writeSimulation(simulation::Simulation;\n                     filename::String=\"\",\n                     folder::String=\".\",\n                     verbose::Bool=true)\n\nWrite all content from Simulation to disk in JDL format.  If the filename  parameter is not specified, it will be saved to a subdirectory under the current  directory named after the simulation identifier simulation.id.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeSimulationStatus-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeSimulationStatus",
    "category": "Method",
    "text": "writeSimulationStatus(simulation::Simulation;\n                      folder::String=\".\",\n                      verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.writeVTK-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.writeVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all ice floes in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.\n\nIf the simulation contains an Ocean data structure, it's contents will be  written to separate .vtu files.  This can be disabled by setting the argument  ocean=false.  The same is true for the atmosphere.\n\nThe VTK files will be saved in a subfolder named after the simulation.\n\n\n\n"
},

{
    "location": "module.html#SeaIce.zeroForcesAndTorques!-Tuple{SeaIce.Simulation}",
    "page": "Modules, constants, types, functions, and macros",
    "title": "SeaIce.zeroForcesAndTorques!",
    "category": "Method",
    "text": "Sets the force and torque values of all ice floes to zero.\n\n\n\n"
},

{
    "location": "module.html#Modules,-constants,-types,-functions,-and-macros-1",
    "page": "Modules, constants, types, functions, and macros",
    "title": "Modules, constants, types, functions, and macros",
    "category": "section",
    "text": "Modules = [SeaIce]\nPrivate = false"
},

]}
