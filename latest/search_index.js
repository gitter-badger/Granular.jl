var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Granular.jl-1",
    "page": "Home",
    "title": "Granular.jl",
    "category": "section",
    "text": "A Julia package for particle-based simulation of  granular dynamics.Granular.jl is a flexible and computationally efficient 2d implementation of  the discrete element method.  Grains are represented as particles, which can be  forced by drag in fluid grids.  The grains interact through  elasto-viscous-frictional contact rheologies and can obtain time-dependent  tensile strength.The source code for Granular.jl is hosted on Github.See the Public API Index for the complete list of documented functions and types."
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
    "text": "Granular.jl is licensed under the GPLv3; see LICENSE for the full license text."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"man/installation.md\",\n    \"man/package_contents.md\",\n    \"man/methods.md\",\n    \"man/getting_started.md\",\n]\nDepth = 1"
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
    "text": ""
},

{
    "location": "man/installation.html#Prerequisites-1",
    "page": "Installation",
    "title": "Prerequisites",
    "category": "section",
    "text": "Granular.jl is written as a package for the Julia programming  language, which is a computationally efficient, yet  high-level language. Julia also includes a very useful package manager which  makes it easy to install packages and their requirements, as well as convenient  updating features."
},

{
    "location": "man/installation.html#Installing-Julia-1",
    "page": "Installation",
    "title": "Installing Julia",
    "category": "section",
    "text": "If you do not have Julia installed, download the current release from the  official Julia download page, or using your  system package manager (e.g. brew cask install julia on macOS with the  Homebrew package manager).  Afterwards, the program julia  can be launched from the terminal."
},

{
    "location": "man/installation.html#Installing-Paraview-1",
    "page": "Installation",
    "title": "Installing Paraview",
    "category": "section",
    "text": "The core visualization functionality of Granular.jl is based on VTK and  ParaView.  The most recent stable release can be downloaded from the ParaView  downloads page.  Alternatively, on macOS  with Homebrew, Paraview can be installed from the terminal with brew cask  install paraview."
},

{
    "location": "man/installation.html#Stable-installation-(recommended)-1",
    "page": "Installation",
    "title": "Stable installation (recommended)",
    "category": "section",
    "text": "The latest stable release of Granular.jl can be installed directly from the  Julia shell by:julia> Pkg.add(\"Granular\")This will install the contents of this repository in the folder  ~/.julia/v$(JULIA_VERSION)/Granular and install its requirements.  The  package JLD is used for model restarts and  is recommended but not required, and is thus not automatically installed."
},

{
    "location": "man/installation.html#Development-installation-1",
    "page": "Installation",
    "title": "Development installation",
    "category": "section",
    "text": "If desired, the current developmental version of the Granular.jl Github  repository can be installed with the  command:julia> Pkg.clone(\"git://github.com/anders-dc/Granular.jl\")Please note: The developmental version is considered unstable and should only  be used over the stable version if there is a compelling reason to do so."
},

{
    "location": "man/installation.html#Keeping-the-package-up-to-date-1",
    "page": "Installation",
    "title": "Keeping the package up to date",
    "category": "section",
    "text": "With the Pkg.update() command, Julia checks and updates all installed  packages to their newest versions."
},

{
    "location": "man/installation.html#Package-tests-1",
    "page": "Installation",
    "title": "Package tests",
    "category": "section",
    "text": "The Granular.jl package contains many tests that verify that the functionality  works as intended.  The extent of test coverage of the source code is monitored  and published with CodeCov.The package tests are during development continuously run with  Travis-CI for Mac (latest stable  release) and Linux (Ubuntu stable (trusty)), and  AppVeyor for Windows.The test scripts are contained in the test/ directory, can be run locally  with the following command:julia> Pkg.test(\"Granular\")In case any of these tests fail, please open a Github  Issue describing the problems  so further investigation and diagnosis can follow."
},

{
    "location": "man/package_contents.html#",
    "page": "Package contents",
    "title": "Package contents",
    "category": "page",
    "text": ""
},

{
    "location": "man/package_contents.html#Package-contents-1",
    "page": "Package contents",
    "title": "Package contents",
    "category": "section",
    "text": "This package follows the official  guidelines  for Julia package layout and contents. "
},

{
    "location": "man/package_contents.html#File-locations-after-installation-1",
    "page": "Package contents",
    "title": "File locations after installation",
    "category": "section",
    "text": "After installation, the package contents will be installed inside the hidden  ~/.julia/ folder in the home directory.  The path can be printed from inside  the julia shell by the command:julia> Pkg.dir(\"Granular\")\n\"/Users/ad/.julia/v0.7/Granular\"The above output will be different for different platforms and Julia versions.  In order to open this directory on macOS, run the following command:julia> run(`open $(Pkg.dir(\"Granular\"))`)On Linux, use the following command:julia> run(`xdg-open $(Pkg.dir(\"Granular\"))`)The above commands will open the directory containing all of the Granular.jl  components. The main component of Granular.jl is the source code contained in  the src/ directory.  The docs/  directory contains the documentation source via Markdown files.  The online  documentation is generated from these files via  Documenter.jl by the  docs/make.jl  script.  The documentation consists of manual pages, as well as auto-generated  API reference that is parsed from the documentation of the Granular.jl source  code (src/  directory)."
},

{
    "location": "man/package_contents.html#Example-scripts-1",
    "page": "Package contents",
    "title": "Example scripts",
    "category": "section",
    "text": "The examples  directory contains several annotated examples, which are useful for getting  started with the Granular.jl package and for demonstrating some of its  features.  The examples are generally heavily annotated with comments to  illustrate the purpose of the included commands.The examples can be run by either navigating to the examples directory from the  command line, and launching them with a command like julia -e logo.jl, or  directly from the julia shell with:julia> include(\"$(Pkg.dir(\"Granular\"))/examples/logo.jl\")It is recommended that the source code of the examples is inspected beforehand."
},

{
    "location": "man/methods.html#",
    "page": "Computational methods",
    "title": "Computational methods",
    "category": "page",
    "text": ""
},

{
    "location": "man/methods.html#Computational-methods-1",
    "page": "Computational methods",
    "title": "Computational methods",
    "category": "section",
    "text": ""
},

{
    "location": "man/getting_started.html#",
    "page": "Getting started",
    "title": "Getting started",
    "category": "page",
    "text": ""
},

{
    "location": "man/getting_started.html#Getting-started-1",
    "page": "Getting started",
    "title": "Getting started",
    "category": "section",
    "text": "In the following, two simple examples are presented using some of the core  commands of Granular.jl.  For more examples, see the scripts in the examples/  directory.The relevant functions are all contained in the Granular module, which can be  imported with import Granular at the top of your script.  Note: As per  Julia conventions, functions that contain an exclamation mark (!) modify the  values of the arguments.Any of the functions called below are documented in the source code, and this  documentation can be found in the Public API Index in the  online documentation, or simply from the Julia shell by typing ?<function  name>.  An example:julia> ?Granular.fitGridToGrains!\n\n  fitGridToGrains!(simulation, grid[, padding])\n\n  Fit the ocean or atmosphere grid for a simulation to the current grains and their positions.\n\n     Arguments\n    ≡≡≡≡≡≡≡≡≡≡≡\n\n    •    simulation::Simulation: simulation object to manipulate.\n\n    •    grid::Any: Ocean or Atmosphere grid to manipulate.\n\n    •    padding::Real: optional padding around edges [m].\n\n    •    verbose::Bool: show grid information when function completes.\n"
},

{
    "location": "man/getting_started.html#Collision-between-two-particles-1",
    "page": "Getting started",
    "title": "Collision between two particles",
    "category": "section",
    "text": "For this simple example (example/two-grains.jl), we will create two grains,  where one of the grains is bumping in to the other.As the first command, we import all the Granular.jl functionality:julia> import Granular"
},

{
    "location": "man/getting_started.html#Simulation-setup-1",
    "page": "Getting started",
    "title": "Simulation setup",
    "category": "section",
    "text": "Next, we create our simulation object which holds all information on the  simulated grains.  This object can be called whatever is appropriate.  In this  documentation, we will use the name sim:julia> sim = Granular.createSimulation(id=\"two-grains\")\nGranular.Simulation(\"two-grains\", 0, 0.0, -1.0, -1.0, -1.0, 0, 0.0, \nGranular.GrainCylindrical[], Granular.Ocean(false, [0.0], [0.0], [0.0], [0.0], \n[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], Array{Int64,1}[#undef], 1, 1, \n1, 1), Granular.Atmosphere(false, [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], \n[0.0], [0.0], Array{Int64,1}[#undef], 1, 1, 1, 1, false), 16)We will be presented with some output about the contents of the sim  simulation object.  This is of minor importance as of now, and can safely be  ignored.During the above createSimulation call, the id argument is optional, but is  used to name simulation output files that are written to the disk.  It is good  practice to use the same name as used for the simulation script file."
},

{
    "location": "man/getting_started.html#Adding-grains-one-by-one-1",
    "page": "Getting started",
    "title": "Adding grains one by one",
    "category": "section",
    "text": "We have now created a simulation object, which will be used during all of the  following.  Next, we add grains to this object.  The first grain is cylinder  shaped, placed at the x-y position (0,0) m, has a radius of 0.1 m, and a  thickness of 0.05 m (along z).  As this call modifies the sim object, the  function contains an exclamation mark (!).  For further information regarding  this call, see the reference in the Public API Index.julia> Granular.addGrainCylindrical!(sim, [0.0, 0.0], 0.1, 0.05)\nINFO: Added Grain 1Let's add another grain, placed at some distance from the first grain:julia> Granular.addGrainCylindrical!(sim, [0.5, 0.0], 0.1, 0.05)\nINFO: Added Grain 2We now want to prescribe a linear (not rotational or angular) velocity to the  first grain, to make it bump into the second grain.The simulation object sim contains an array of all grains that are added to  it.  We can directly inspect the grains and their values from the simulation  object.  Let's take a look at the default value of the linear velocity, called  lin_vel:julia> sim.grains[1].lin_vel\n2-element Array{Float64, 1}:\n 0.0\n 0.0The default value is a 0,0 vector.  Similarly, we can modify the properties of  the first grain directly with the following call:julia> sim.grains[1].lin_vel = [1.0, 0.0]\n2-element Array{Float64, 1}:\n 1.0\n 0.0The first grain (index 1 in sim.grains) now has a positive velocity along x  with the value of 1.0 meter per second."
},

{
    "location": "man/getting_started.html#Setting-temporal-parameters-for-the-simulation-1",
    "page": "Getting started",
    "title": "Setting temporal parameters for the simulation",
    "category": "section",
    "text": "Before we can start the simulation, we need to tell the code vital information,  like what time step to use, how often to write output files to the disk, and  for how long to run the simulation.  To set the computational time step, we  call the following:julia> Granular.setTimeStep!(sim)\nINFO: Time step length t=8.478741928617433e-5 sBased on the elastic stiffness of the grains and their size, a suitable time  step is automatically determined.The two grains are 0.3 meters apart, as the centers are placed 0.5 meter from  each other and each grain has a diameter of 0.1 m.  With a linear velocity of  1.0 m/s, the collision should occur after 0.3 seconds.  For that reason it  seems reasonable to run the simulation for a total of 1.0 seconds, and produce  an output file every 0.05 seconds.  We will later use the output files for  visualization purposes.julia> Granular.setOutputFileInterval!(sim, 0.05)\n\njulia> Granular.setTotalTime!(sim, 1.0)"
},

{
    "location": "man/getting_started.html#Running-the-simulation-1",
    "page": "Getting started",
    "title": "Running the simulation",
    "category": "section",
    "text": "We are now ready to run the simulation.  For illustrative purposes, let us  compare the kinetic energy in the granular system before and after the  collision.  For now, we save the initial value using the following call:julia> Granular.totalGrainKineticTranslationalEnergy(sim)\n0.7335618846132168The above value is the total translational (not angular) kinetic energy in  Joules before the simulation is started.For running the simulation, we have two choices; we can either run the entire  simulation length with a single call, which steps time until the total time is  reached and generates output files along the way.  Alternatively, we can run  the simulation for a single time step a time, and inspect the progress or do  other modifications along the way.Here, we will run the entire simulation in one go, and afterwards visualize the  grains from their output files using ParaView.julia> Granular.run!(sim)\n\nINFO: Output file: ./two-grains/two-grains.grains.1.vtu\nINFO: wrote status to ./two-grains/two-grains.status.txt\n  t = 0.04239370964308682/1.0 s\nINFO: Output file: ./two-grains/two-grains.grains.2.vtu\nINFO: wrote status to ./two-grains/two-grains.status.txt\n\n...\n\nINFO: Output file: ./two-grains/two-grains.grains.20.vtu\nINFO: wrote status to ./two-grains/two-grains.status.txt\n  t = 0.9920128056483803/1.0 s\nINFO: ./two-grains/two-grains.py written, execute with 'pvpython /Users/ad/code/Granular-ext/two-grains/two-grains.py'\nINFO: wrote status to ./two-grains/two-grains.status.txt\n  t = 1.0000676104805686/1.0 sThe script informs us of the simulation progress, and notifies us as output  files are generated (this can be disabled by passing verbose=false to the  run!() command).  Finally, it tells us that it has generated a ParaView  python file for visualization.Before going further, we are interested in getting an immediate idea of how the  collision went.  We print the new velocities with the following commands:julia> sim.grains[1].lin_vel\n2-element Array{Float64, 1}:\n 7.58343e-5\n 0.0\n\njulia> sim.grains[2].lin_vel\n2-element Array{Float64, 1}:\n 0.999924\n 0.0The first grain has transferred effectively all of its kinetic energy to the  second grain during the cause of the simulation.  The total kinetic energy now  is the following:julia> Granular.totalGrainKineticTranslationalEnergy(sim)\n0.7334506347624973The before and after values are reasonably close (to less than 0.1 percent),  which is what can be expected given the computation accuracy in the algorithm."
},

{
    "location": "man/getting_started.html#Visualizing-the-output-1",
    "page": "Getting started",
    "title": "Visualizing the output",
    "category": "section",
    "text": "To visualize the output we open ParaView.  The  output files of the simulation are written using the VTK (visualization  toolkit) format, which is natively supported by ParaView.While the .vtu files produced during the simulation can be opened with  ParaView and visualized manually using Glyph filters, the simplest and  fastest way to visualize the data is to use the Python script generated for the  simulation by Granular.jl.Open ParaView and open the Python Shell, found under the menu Tools > Python  Shell.  In the pop-up dialog we select Run Script, which opens yet another  dialog prompting us to locate the visualization script (two-grains.py, in our  example).  We locate this file, which is placed under the directory from where  we launched the julia session with the commands above.After selecting the two-grains/two-grains.py script, we can close the Python  Shell window to inspect our simulation.  Press the Play symbol in the top  toolbar, and see what happens!Alternatively, you can color the grains using different parameters, such as  velocity, number of contacts, etc.  These can be selected by changing the  chosen parameter under the Glyph1 object in the Pipeline Browser on the  left, and selecting a different field for Coloring.  Press the Apply button  to see the changes in effect."
},

{
    "location": "man/getting_started.html#Exercises-1",
    "page": "Getting started",
    "title": "Exercises",
    "category": "section",
    "text": "To gain more familiarity with the simulation procedure, I suggest experimenting  with the following:What effect does the grain size have on the time step?\nTry to make an oblique collision by placing one of the grains at a different    y position.\nWhat happens if the second grains is set to be fixed in space    (sim.grains[2].fixed = true)?\nHow is the relationship between total kinetic energy before and after    affected by the choice of time step length?  Try setting different time    step values, e.g. with sim.time_step = 0.1234 and rerun the simulation."
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
    "text": "Documentation for Granular.jl's public interface.See Package-internal documentation for internal package docs."
},

{
    "location": "lib/public.html#main-index-1",
    "page": "Public API",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public.html#Granular",
    "page": "Public API",
    "title": "Granular",
    "category": "Module",
    "text": "Granular.jl\n\nOffline granular dynamics simulator module.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.addAtmosphereDrag!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.addAtmosphereDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between atmosphere and all  grains.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.addBodyForce!-Tuple{Granular.GrainCylindrical,Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.addBodyForce!",
    "category": "Method",
    "text": "setBodyForce!(grain, force)\n\nAdd to the value of the external body force on a grain.\n\nArguments\n\ngrain::GrainCylindrical: the grain to set the body force for.\nforce::Vector{Float64}: a vector of force [N]\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.addGrain!",
    "page": "Public API",
    "title": "Granular.addGrain!",
    "category": "Function",
    "text": "addGrain!(simulation::Simulation,\n            grain::GrainCylindrical,\n            verbose::Bool = False)\n\nAdd an grain to the simulation object.  If verbose is true, a short  confirmation message will be printed to stdout.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.addGrainCylindrical!-Tuple{Granular.Simulation,Array{Float64,1},Float64,Float64}",
    "page": "Public API",
    "title": "Granular.addGrainCylindrical!",
    "category": "Method",
    "text": "function addGrainCylindrical!(simulation, lin_pos, contact_radius,\n                                thickness[, areal_radius, lin_vel, lin_acc,\n                                force, ang_pos, ang_vel, ang_acc, torque,\n                                density, contact_stiffness_normal,\n                                contact_stiffness_tangential,\n                                contact_viscosity_normal,\n                                contact_viscosity_tangential,\n                                contact_static_friction,\n                                contact_dynamic_friction,\n                                youngs_modulus, poissons_ratio,\n                                tensile_strength, tensile_heal_rate,\n                                compressive_strength_prefactor,\n                                ocean_drag_coeff_vert,\n                                ocean_drag_coeff_horiz,\n                                atmosphere_drag_coeff_vert,\n                                atmosphere_drag_coeff_horiz,\n                                pressure, fixed, rotating, enabled, verbose,\n                                ocean_grid_pos, atmosphere_grid_pos,\n                                n_contact, granular_stress, ocean_stress,\n                                atmosphere_stress])\n\nCreates and adds a cylindrical grain to a simulation. Most of the arguments  are optional, and come with default values.  The only required arguments are  simulation, lin_pos, contact_radius, and thickness.\n\nArguments\n\nsimulation::Simulation: the simulation object where the grain should be   added to.\nlin_pos::Vector{Float64}: linear position of grain center [m].\ncontact_radius::Float64: grain radius for granular interaction [m].\nthickness::Float64: grain thickness [m].\nareal_radius = false: grain radius for determining sea-ice concentration   [m].\nlin_vel::Vector{Float64} = [0., 0.]: linear velocity [m/s].\nlin_acc::Vector{Float64} = [0., 0.]: linear acceleration [m/s^2].\nforce::Vector{Float64} = [0., 0.]: linear force balance [N].\nang_pos::Float64 = 0.: angular position around its center vertical axis   [rad].\nang_vel::Float64 = 0.: angular velocity around its center vertical axis   [rad/s].\nang_acc::Float64 = 0.: angular acceleration around its center vertical axis   [rad/s^2].\ntorque::Float64 = 0.: torque around its center vertical axis [N*m]\ndensity::Float64 = 934.: grain mean density [kg/m^3].\ncontact_stiffness_normal::Float64 = 1e7: contact-normal stiffness [N/m];   overridden if youngs_modulus is set to a positive value.\ncontact_stiffness_tangential::Float64 = 0.: contact-tangential stiffness   [N/m]; overridden if youngs_modulus is set to a positive value.\ncontact_viscosity_normal::Float64 = 0.: contact-normal viscosity [N/m/s].\ncontact_viscosity_tangential::Float64 = 0.: contact-tangential viscosity   [N/m/s].\ncontact_static_friction::Float64 = 0.4: contact static Coulomb frictional   coefficient [-].\ncontact_dynamic_friction::Float64 = 0.4: contact dynamic Coulomb frictional   coefficient [-].\nyoungs_modulus::Float64 = 2e7: elastic modulus [Pa]; overrides any value   set for contact_stiffness_normal.\npoissons_ratio::Float64 = 0.185: Poisson's ratio, used to determine the   contact-tangential stiffness from youngs_modulus [-].\ntensile_strength::Float64 = 0.: contact-tensile (cohesive) bond strength   [Pa].\ntensile_heal_rate::Float64 = 0.: rate at which contact-tensile bond strength   is obtained [1/s].\ncompressive_strength_prefactor::Float64 = 1285e3: maximum compressive   strength on granular contact (not currently enforced) [m*Pa].\nocean_drag_coeff_vert::Float64 = 0.85: vertical drag coefficient for ocean   against grain sides [-].\nocean_drag_coeff_horiz::Float64 = 5e-4: horizontal drag coefficient for   ocean against grain bottom [-].\natmosphere_drag_coeff_vert::Float64 = 0.4: vertical drag coefficient for   atmosphere against grain sides [-].\natmosphere_drag_coeff_horiz::Float64 = 2.5e-4: horizontal drag coefficient   for atmosphere against grain bottom [-].\npressure::Float64 = 0.: current compressive stress on grain [Pa].\nfixed::Bool = false: grain is fixed to a constant velocity (e.g. zero).\nrotating::Bool = true: grain is allowed to rotate.\nenabled::Bool = true: grain interacts with other grains.\nverbose::Bool = true: display diagnostic information during the function   call.\nocean_grid_pos::Array{Int, 1} = [0, 0]: position of grain in the ocean   grid.\natmosphere_grid_pos::Array{Int, 1} = [0, 0]: position of grain in the   atmosphere grid.\nn_contacts::Int = 0: number of contacts with other grains.\ngranular_stress::Vector{Float64} = [0., 0.]: resultant stress on grain   from granular interactions [Pa].\nocean_stress::Vector{Float64} = [0., 0.]: resultant stress on grain from   ocean drag [Pa].\natmosphere_stress::Vector{Float64} = [0., 0.]: resultant stress on grain   from atmosphere drag [Pa].\n\nExamples\n\nThe most basic example adds a new grain to the simulation sim, with a  center at [1., 2.], a radius of 1. meter, and a thickness of 0.5  meter:\n\nGranular.addGrainCylindrical!(sim, [1., 2.], 1., .5)\n\nThe following example will create a grain with tensile strength (cohesion), and a velocity of 0.5 m/s towards -x:\n\nGranular.addGrainCylindrical!(sim, [4., 2.], 1., .5,\n                              tensile_strength = 200e3,\n                              lin_vel = [-.5, 0.])\n\nFixed grains are useful for creating walls or coasts, and loops are useful for creating regular arrangements:\n\nfor i=1:5\n    Granular.addGrainCylindrical!(sim, [i*2., 0., 3.], 1., .5, fixed=true)\nend\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.addOceanDrag!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.addOceanDrag!",
    "category": "Method",
    "text": "Add drag from linear and angular velocity difference between ocean and all ice  floes.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.applyAtmosphereDragToGrain!-Tuple{Granular.GrainCylindrical,Float64,Float64}",
    "page": "Public API",
    "title": "Granular.applyAtmosphereDragToGrain!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between atmosphere and a single  grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.applyAtmosphereVorticityToGrain!-Tuple{Granular.GrainCylindrical,Float64}",
    "page": "Public API",
    "title": "Granular.applyAtmosphereVorticityToGrain!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between atmosphere and a  single grain.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.applyOceanDragToGrain!-Tuple{Granular.GrainCylindrical,Float64,Float64}",
    "page": "Public API",
    "title": "Granular.applyOceanDragToGrain!",
    "category": "Method",
    "text": "Add Stokes-type drag from velocity difference between ocean and a single ice  floe.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.applyOceanVorticityToGrain!-Tuple{Granular.GrainCylindrical,Float64}",
    "page": "Public API",
    "title": "Granular.applyOceanVorticityToGrain!",
    "category": "Method",
    "text": "Add Stokes-type torque from angular velocity difference between ocean and a  single grain.  See Eq. 9.28 in \"Introduction to Fluid Mechanics\" by Nakayama  and Boucher, 1999.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.areaOfQuadrilateral-NTuple{4,Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.areaOfQuadrilateral",
    "category": "Method",
    "text": "Returns the area of a quadrilateral with corner coordinates a, b, c, and  d.  Corners a and c should be opposite of each other, the same must be  true for b and d.  This is true if the four corners are passed as arguments  in a \"clockwise\" or \"counter-clockwise\" manner.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.areaOfTriangle-Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.areaOfTriangle",
    "category": "Method",
    "text": "Returns the area of an triangle with corner coordinates a, b, and c.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.checkAndAddContact!",
    "page": "Public API",
    "title": "Granular.checkAndAddContact!",
    "category": "Function",
    "text": "checkAndAddContact!(simulation, i, j)\n\nCheck for contact between two grains and register the interaction in the  simulation object.  The indexes of the two grains is stored in  simulation.contact_pairs as [i, j].  The overlap vector is parallel to a  straight line connecting the grain centers, points away from grain i and  towards j, and is stored in simulation.overlaps.  A zero-length vector is  written to simulation.contact_parallel_displacement.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the grains.\ni::Int: index of the first grain.\nj::Int: index of the second grain.\ndistance_Modifier::Vector{Float64}: vector modifying percieved   inter-particle distance, which is used for contact search across periodic   boundaries.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.checkTimeParameters-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.checkTimeParameters",
    "category": "Method",
    "text": "Checks if simulation temporal parameters are of reasonable values.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.compareAtmospheres-Tuple{Granular.Atmosphere,Granular.Atmosphere}",
    "page": "Public API",
    "title": "Granular.compareAtmospheres",
    "category": "Method",
    "text": "compareAtmospheres(atmosphere1::atmosphere, atmosphere2::atmosphere)\n\nCompare values of two atmosphere objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.compareGrains-Tuple{Granular.GrainCylindrical,Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.compareGrains",
    "category": "Method",
    "text": "compareGrains(if1::GrainCylindrical, if2::GrainCylindrical)\n\nCompare values of two grain objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.compareOceans-Tuple{Granular.Ocean,Granular.Ocean}",
    "page": "Public API",
    "title": "Granular.compareOceans",
    "category": "Method",
    "text": "compareOceans(ocean1::Ocean, ocean2::Ocean)\n\nCompare values of two Ocean objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.compareSimulations-Tuple{Granular.Simulation,Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.compareSimulations",
    "category": "Method",
    "text": "compareSimulations(sim1::Simulation, sim2::Simulation)\n\nCompare values of two Simulation objects using the Base.Test framework.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.conformalQuadrilateralCoordinates-NTuple{5,Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.conformalQuadrilateralCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional coordinates [x_tilde, y_tilde] of a point p  within a quadrilateral with corner coordinates A, B, C, and D. Points must be ordered in counter-clockwise order, starting from south-west  corner.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.convertGrainDataToArrays-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.convertGrainDataToArrays",
    "category": "Method",
    "text": "Gathers all grain parameters from the GrainCylindrical type in continuous  arrays in an GrainArrays structure.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.createEmptyAtmosphere-Tuple{}",
    "page": "Public API",
    "title": "Granular.createEmptyAtmosphere",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.createEmptyOcean-Tuple{}",
    "page": "Public API",
    "title": "Granular.createEmptyOcean",
    "category": "Method",
    "text": "Returns empty ocean type for initialization purposes.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.createRegularAtmosphereGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.createRegularAtmosphereGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Atmosphere grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the atmosphere spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.createRegularOceanGrid-Tuple{Array{Int64,1},Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.createRegularOceanGrid",
    "category": "Method",
    "text": "Initialize and return a regular, Cartesian Ocean grid with n[1] by n[2]  cells in the horizontal dimension, and n[3] vertical cells.  The cell corner  and center coordinates will be set according to the grid spatial dimensions  L[1], L[2], and L[3].  The grid u, v, h, and e fields will contain  one 4-th dimension matrix per time step.  Sea surface will be at z=0. with  the ocean spanning z<0..  Vertical indexing starts with k=0 at the sea  surface, and increases downwards.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.createSimulation-Tuple{}",
    "page": "Public API",
    "title": "Granular.createSimulation",
    "category": "Method",
    "text": "createSimulation([id::String=\"unnamed\",\n                  time_iteration::Int=0,\n                  time::Float64=0.0,\n                  time_total::Float64=-1.,\n                  time_step::Float64=-1.,\n                  file_time_step::Float64=-1.,\n                  file_number::Int=0,\n                  grains=Array{GrainCylindrical, 1}[],\n                  ocean::Ocean,\n                  atmosphere::Atmosphere)\n\nCreate a simulation object containing all relevant variables such as temporal  parameters, and lists of grains and contacts. The parameter id is used to uniquely identify the simulation when it is written to disk.\n\nArguments\n\nid::String=\"unnamed\":\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.disableAtmosphereDrag!-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.disableAtmosphereDrag!",
    "category": "Method",
    "text": "disableAtmosphereDrag!(grain)\n\nDisable atmosphere-caused drag on the grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.disableGrain!-Tuple{Granular.Simulation,Int64}",
    "page": "Public API",
    "title": "Granular.disableGrain!",
    "category": "Method",
    "text": "Disable grain with index i in the simulation object.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.disableOceanDrag!-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.disableOceanDrag!",
    "category": "Method",
    "text": "disableOceanDrag!(grain)\n\nDisable ocean-caused drag on the grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.disableOutputFiles!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.disableOutputFiles!",
    "category": "Method",
    "text": "Disables the write of output files to disk during a simulation.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.enableAtmosphereDrag!-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.enableAtmosphereDrag!",
    "category": "Method",
    "text": "enableAtmosphereDrag!(grain)\n\nEnable atmosphere-caused drag on the grain, with values by Hunke and Comeau (2011).\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.enableOceanDrag!-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.enableOceanDrag!",
    "category": "Method",
    "text": "enableOceanDrag!(grain)\n\nEnable ocean-caused drag on the grain, with values by Hunke and Comeau (2011).\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findCellContainingPoint",
    "page": "Public API",
    "title": "Granular.findCellContainingPoint",
    "category": "Function",
    "text": "findCellContainingPoint(grid, point[, method])\n\nReturns the i, j index of the grid cell containing the point. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns the coordinates of the cell.  If no match is  found the function returns (0,0).\n\nArguments\n\ngrid::Any: grid object containing ocean or atmosphere data.\npoint::Vector{Float64}: two-dimensional vector of point to check.\nmethod::String: approach to use for determining if point is inside cell or    not, can be \"Conformal\" (default) or \"Areal\".\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findContacts!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.findContacts!",
    "category": "Method",
    "text": "findContacts!(simulation[, method])\n\nTop-level function to perform an inter-grain contact search, based on grain  linear positions and contact radii.\n\nThe simplest contact search algorithm (method=\"all to all\") is the most  computationally expensive (O(n^2)).  The method \"ocean grid\" bins the grains  into their corresponding cells on the ocean grid and searches for contacts only  within the vicinity.  When this method is applied, it is assumed that the  contact_radius values of the grains are smaller than half the cell size.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the grains.\nmethod::String: the contact-search method to apply.  Valid options are \"all    to all\" and \"ocean grid\".\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findContactsAllToAll!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.findContactsAllToAll!",
    "category": "Method",
    "text": "findContactsAllToAll!(simulation)\n\nPerform an O(n^2) all-to-all contact search between all grains in the  simulation object.  Contacts between fixed grains are ignored.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findContactsInGrid!-Tuple{Granular.Simulation,Any}",
    "page": "Public API",
    "title": "Granular.findContactsInGrid!",
    "category": "Method",
    "text": "findContactsInGrid!(simulation)\n\nPerform an O(n*log(n)) cell-based contact search between all grains in the  simulation object.  Contacts between fixed or disabled grains are ignored.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findEmptyPositionInGridCell-Tuple{Granular.Simulation,Any,Int64,Int64,Float64}",
    "page": "Public API",
    "title": "Granular.findEmptyPositionInGridCell",
    "category": "Method",
    "text": "Attempt locate an empty spot for an grain with radius r with center  coordinates in a specified grid cell (i, j) without overlapping any other  grains in that cell or the neighboring cells.  This function will stop  attempting after n_iter iterations, each with randomly generated positions.\n\nThis function assumes that existing grains have been binned according to the  grid (e.g., using sortGrainsInGrid()).\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findLargestGrainStiffness-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.findLargestGrainStiffness",
    "category": "Method",
    "text": "Finds the largest elastic stiffness of all grains in a simulation.  Used to  determine the optimal time step length.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.findSmallestGrainMass-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.findSmallestGrainMass",
    "category": "Method",
    "text": "Finds the smallest mass of all grains in a simulation.  Used to determine  the optimal time step length.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.fitGridToGrains!-Tuple{Granular.Simulation,Any}",
    "page": "Public API",
    "title": "Granular.fitGridToGrains!",
    "category": "Method",
    "text": "fitGridToGrains!(simulation, grid[, padding])\n\nFit the ocean or atmosphere grid for a simulation to the current grains and their positions.\n\nArguments\n\nsimulation::Simulation: simulation object to manipulate.\ngrid::Any: Ocean or Atmosphere grid to manipulate.\npadding::Real: optional padding around edges [m].\nverbose::Bool: show grid information when function completes.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.getCellCenterCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Public API",
    "title": "Granular.getCellCenterCoordinates",
    "category": "Method",
    "text": "getCellCenterCoordinates(grid, i, j)\n\nReturns grid center coordinates (h-point).\n\nArguments\n\nxh::Array{Float64, 2}: nominal longitude of h-points [degrees_E]\nyh::Array{Float64, 2}: nominal latitude of h-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.getCellCornerCoordinates-Tuple{Array{Float64,2},Array{Float64,2},Int64,Int64}",
    "page": "Public API",
    "title": "Granular.getCellCornerCoordinates",
    "category": "Method",
    "text": "getCellCornerCoordinates(xq, yq, i, j)\n\nReturns grid-cell corner coordinates in the following order (south-west corner,  south-east corner, north-east corner, north-west corner).\n\nArguments\n\nxq::Array{Float64, 2}: nominal longitude of q-points [degrees_E]\nyq::Array{Float64, 2}: nominal latitude of q-points [degrees_N]\ni::Int: x-index of cell.\nj::Int: y-index of cell.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.getGridCornerCoordinates-Tuple{Array{Float64,2},Array{Float64,2}}",
    "page": "Public API",
    "title": "Granular.getGridCornerCoordinates",
    "category": "Method",
    "text": "getGridCornerCoordinates(xq, yq)\n\nReturns grid corner coordinates in the following order (south-west corner,  south-east corner, north-east corner, north-west corner).\n\nArguments\n\nxq::Array{Float64, 2}: nominal longitude of q-points [degrees_E]\nyq::Array{Float64, 2}: nominal latitude of q-points [degrees_N]\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.getNonDimensionalCellCoordinates-Tuple{Any,Int64,Int64,Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.getNonDimensionalCellCoordinates",
    "category": "Method",
    "text": "Returns the non-dimensional conformal mapped coordinates for point point in  cell i,j, based off the coordinates in the grid.\n\nThis function is a wrapper for getCellCornerCoordinates() and  conformalQuadrilateralCoordinates().\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainCircumreference-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainCircumreference",
    "category": "Method",
    "text": "Returns the circumreference of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainHorizontalSurfaceArea-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainHorizontalSurfaceArea",
    "category": "Method",
    "text": "Returns the top or bottom (horizontal) surface area of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainKineticRotationalEnergy-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainKineticRotationalEnergy",
    "category": "Method",
    "text": "Returns the rotational kinetic energy of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainKineticTranslationalEnergy-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainKineticTranslationalEnergy",
    "category": "Method",
    "text": "Returns the translational kinetic energy of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainMass-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainMass",
    "category": "Method",
    "text": "Returns the mass of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainMomentOfInertia-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainMomentOfInertia",
    "category": "Method",
    "text": "Returns the moment of inertia of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainSideSurfaceArea-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainSideSurfaceArea",
    "category": "Method",
    "text": "Returns the surface area of the grain sides\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.grainVolume-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.grainVolume",
    "category": "Method",
    "text": "Returns the volume of the grain\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.harmonicMean-Tuple{Number,Number}",
    "page": "Public API",
    "title": "Granular.harmonicMean",
    "category": "Method",
    "text": "harmonicMean(a, b)\n\nReturns the harmonic mean of two numbers a::Number and b::Number.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.incrementCurrentTime!-Tuple{Granular.Simulation,Float64}",
    "page": "Public API",
    "title": "Granular.incrementCurrentTime!",
    "category": "Method",
    "text": "incrementCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interGrainPositionVector-Tuple{Granular.Simulation,Int64,Int64}",
    "page": "Public API",
    "title": "Granular.interGrainPositionVector",
    "category": "Method",
    "text": "interGrainPositionVector(simulation, i, j)\n\nReturns a vector pointing from grain i to grain j in the  simulation.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the grains.\ni::Int: index of the first grain.\nj::Int: index of the second grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interact!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.interact!",
    "category": "Method",
    "text": "interact!(simulation::Simulation)\n\nResolve mechanical interaction between all particle pairs.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interactGrains!-Tuple{Granular.Simulation,Int64,Int64,Int64}",
    "page": "Public API",
    "title": "Granular.interactGrains!",
    "category": "Method",
    "text": "interactGrains!(simulation::Simulation, i::Int, j::Int, ic::Int)\n\nResolve an grain-to-grain interaction using a prescibed contact law.  This  function adds the compressive force of the interaction to the grain  pressure field of mean compressive stress on the grain sides.\n\nThe function uses the macroscopic contact-stiffness parameterization based on  Young's modulus and Poisson's ratio if Young's modulus is a positive value.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interpolateAtmosphereState-Tuple{Granular.Atmosphere,Float64}",
    "page": "Public API",
    "title": "Granular.interpolateAtmosphereState",
    "category": "Method",
    "text": "Atmosphere data is containted in Atmosphere type at discrete times  (Atmosphere.time).  This function performs linear interpolation between time  steps to get the approximate atmosphere state at any point in time.  If the  Atmosphere data set only contains a single time step, values from that time  are returned.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interpolateAtmosphereVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Public API",
    "title": "Granular.interpolateAtmosphereVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interpolateOceanState-Tuple{Granular.Ocean,Float64}",
    "page": "Public API",
    "title": "Granular.interpolateOceanState",
    "category": "Method",
    "text": "Ocean data is containted in Ocean type at discrete times (Ocean.time).  This  function performs linear interpolation between time steps to get the approximate  ocean state at any point in time.  If the Ocean data set only contains a  single time step, values from that time are returned.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.interpolateOceanVelocitiesToCorners-Tuple{Array{Float64,4},Array{Float64,4}}",
    "page": "Public API",
    "title": "Granular.interpolateOceanVelocitiesToCorners",
    "category": "Method",
    "text": "Convert gridded data from Arakawa-C type (decomposed velocities at faces) to  Arakawa-B type (velocities at corners) through interpolation.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.isPointInCell",
    "page": "Public API",
    "title": "Granular.isPointInCell",
    "category": "Function",
    "text": "Check if a 2d point is contained inside a cell from the supplied grid. The function uses either an area-based approach (method = \"Area\"), or a  conformal mapping approach (method = \"Conformal\").  The area-based approach is  more robust.  This function returns true or false.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.isPointInGrid",
    "page": "Public API",
    "title": "Granular.isPointInGrid",
    "category": "Function",
    "text": "Check if a 2d point is contained inside the grid.  The function uses either an area-based approach (method = \"Area\"), or a conformal mapping approach (method = \"Conformal\").  The area-based approach is more robust.  This function returns true or false.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.plotGrainSizeDistribution-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.plotGrainSizeDistribution",
    "category": "Method",
    "text": "plotGrainSizeDistribution(simulation, [filename_postfix], [nbins],\n                            [size_type], [figsize], [filetype])\n\nPlot the grain size distribution as a histogram and save it to the disk.  The  plot is saved accoring to the simulation id, the optional filename_postfix  string, and the filetype, and is written to the current folder.\n\nArguments\n\nsimulation::Simulation: the simulation object containing the grains.\nfilename_postfix::String: optional string for the output filename.\nnbins::Int: number of bins in the histogram (default = 12).\nsize_type::String: specify whether to use the contact or areal radius    for the grain size.  The default is contact.\nfigsize::Tuple: the output figure size in inches (default = (6,4).\nfiletype::String: the output file type (default = \"png\").\nverbose::String: show output file as info message in stdout (default =    true).\nskip_fixed::Bool: ommit grains that are fixed in space from the size    distribution (default = true).\nlogy::Bool: plot y-axis in log scale.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.poissonDiscSampling-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.poissonDiscSampling",
    "category": "Method",
    "text": "Generate disc packing in 2D using Poisson disc sampling with O(N) complexity, as described by Robert Bridson (2007).\n\nArguments\n\nsimulation::Simulation: simulation object where grains are inserted.\nradius_max::Real: largest grain radius to use.\nradius_min::Real: smallest grain radius to use.\nsample_limit::Integer=30: number of points to sample around each grain   before giving up.\nmax_padding_factor::Real=2.: this factor scales the padding to use during ice   floe generation in numbers of grain diameters.\nverbose::Bool=true: show diagnostic information to stdout.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.printGrainInfo-Tuple{Granular.GrainCylindrical}",
    "page": "Public API",
    "title": "Granular.printGrainInfo",
    "category": "Method",
    "text": "printGrainInfo(grain::GrainCylindrical)\n\nPrints the contents of an grain to stdout in a formatted style.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.printMemoryUsage-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.printMemoryUsage",
    "category": "Method",
    "text": "printMemoryUsage(sim::Simulation)\n\nShows the memory footprint of the simulation object.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.randpower",
    "page": "Public API",
    "title": "Granular.randpower",
    "category": "Function",
    "text": "randpower([nvals], [distribution_power], [min_val], [max_val])\n\nReturns one or more random numbers from a power-law probability distribution.\n\nArguments\n\ndims::Any: the dimensions of random values (default = 1)\ndistribution_power::Number: the distribution power (default = 1.)\nmin_val::Number: the lower bound of the distribution range (default = 0.)\nmax_val::Number: the upper bound of the distribution range (default = 1.)\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.readOceanGridNetCDF-Tuple{String}",
    "page": "Public API",
    "title": "Granular.readOceanGridNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean supergrid information generated by MOM6 (e.g.   ocean_hrid.nc) from disk and return as Ocean data structure.  This file is  located in the simulation INPUT/ subdirectory.\n\nReturns\n\nxh::Array{Float64, 2}: Longitude for cell centers [deg]\nyh::Array{Float64, 2}: Latitude for cell centers [deg]\nxq::Array{Float64, 2}: Longitude for cell corners [deg]\nyq::Array{Float64, 2}: Latitude for cell corners [deg]\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.readOceanNetCDF-Tuple{String,String}",
    "page": "Public API",
    "title": "Granular.readOceanNetCDF",
    "category": "Method",
    "text": "Read ocean NetCDF files generated by MOM6 from disk and return as Ocean data  structure.\n\nArguments\n\nvelocity_file::String: Path to NetCDF file containing ocean velocities,    etc., (e.g. prog__####_###.nc).\ngrid_file::String: Path to NetCDF file containing ocean super-grid    information (typically INPUT/ocean_hgrid.nc).\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.readOceanStateNetCDF-Tuple{String}",
    "page": "Public API",
    "title": "Granular.readOceanStateNetCDF",
    "category": "Method",
    "text": "Read NetCDF file with ocean state generated by MOM6 (e.g.  prog__####_###.nc  or ########.ocean_month.nc) from disk and return time stamps, velocity fields,  layer thicknesses, interface heights, and vertical coordinates.\n\nReturns\n\ntime::Vector{Float64}: Time [s]\nu::Array{Float64, 2}: Cell corner zonal velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nv::Array{Float64, 2}: Cell corner meridional velocity [m/s],   dimensions correspond to placement in [xq, yq, zl, time]\nh::Array{Float64, 2}: layer thickness [m], dimensions correspond to    placement in [xh, yh, zl, time]\ne::Array{Float64, 2}: interface height relative to mean sea level [m],     dimensions correspond to placement in [xh, yh, zi, time]\nzl::Vector{Float64}: layer target potential density [kg m^-3]\nzi::Vector{Float64}: interface target potential density [kg m^-3]\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.readSimulation",
    "page": "Public API",
    "title": "Granular.readSimulation",
    "category": "Function",
    "text": "readSimulation(filename::String=\"\";\n               verbose::Bool=true)\n\nRead all content from Simulation from disk in JDL format.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.readSimulationStatus-Tuple{String}",
    "page": "Public API",
    "title": "Granular.readSimulationStatus",
    "category": "Method",
    "text": "readSimulationStatus(filename::String;\n                     folder::String=\".\",\n                     verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.reflectGrainsFromImpermeableBoundaries!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.reflectGrainsFromImpermeableBoundaries!",
    "category": "Method",
    "text": "reflectGrainsFromImpermeableBoundaries!(simulation::Simulation)\n\nIf the ocean or atmosphere grids are impermeable, reflect grain trajectories by reversing the velocity vectors normal to the boundary.  This function is to be called after temporal integration of the grain positions.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.regularPacking!-Tuple{Granular.Simulation,Array{Int64,1},Real,Real}",
    "page": "Public API",
    "title": "Granular.regularPacking!",
    "category": "Method",
    "text": "regularPacking!(simulation, n, r_min, r_max[, padding_factor,\n                size_distribution, size_distribution_parameter])\n\nCreate a grid-based regular packing with grain numbers along each axis specified by the n vector.\n\nArguments\n\nsimulation::Simulation: simulation object where the grains are inserted,   preferably not containing prior grains.\nn::Vector{Integer}: 2-element vector determining number of grains along the   x and y axes.\nr_min::Real: minimum desired grain radius.\nr_max::Real: maximum desired grain radius.\npadding_factor::Real: percentage-wise padding around each grain to allow for   random perturbations to grain position.\nsize_distribution::String: grain-size distribution to sample. Valid values   are \"powerlaw\" and \"uniform\".\nsize_distribution_parameter::Real: parameter to pass to the grain-size   distribution generating function.\nseed::Integer: seed value to the pseudo-random number generator.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.removeSimulationFiles-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.removeSimulationFiles",
    "category": "Method",
    "text": "removeSimulationFiles(simulation[, folder])\n\nRemove all simulation output files from the specified folder.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.render-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.render",
    "category": "Method",
    "text": "render(simulation[, pvpython, images, animation])\n\nWrapper function which calls writeParaviewPythonScript(...) and executes it from the shell using the supplied pvpython argument.\n\nArguments\n\nsimulation::Simulation: simulation object containing the grain data.\npvpython::String: path to the pvpython executable to use.  By default, the   script uses the pvpython in the system PATH.\nimages::Bool: render images to disk (default: true)\nanimation::Bool: render animation to disk (default: false)\ntrim::Bool: trim images in animated sequence (default: true)\nreverse::Bool: if images=true additionally render reverse-animated gif   (default: false)\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.reportGridBoundaryConditions-Tuple{Any}",
    "page": "Public API",
    "title": "Granular.reportGridBoundaryConditions",
    "category": "Method",
    "text": "reportGridBoundaryConditions(grid)\n\nReport the boundary conditions for the grid to the console.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.reportSimulationTimeToStdout-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.reportSimulationTimeToStdout",
    "category": "Method",
    "text": "Prints the current simulation time and total time to standard out\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.run!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.run!",
    "category": "Method",
    "text": "run!(simulation[,\n     verbose::Bool = true,\n     status_interval = 100.,\n     show_file_output = true,\n     single_step = false,\n     temporal_integration_method = \"Three-term Taylor\"],\n     write_jld = false)\n\nRun the simulation through time until simulation.time equals or exceeds  simulatim.time_total.  This function requires that all grains are added to  the simulation and that the length of the computational time step is adjusted  accordingly.\n\nThe function will search for contacts, determine the force balance on each ice  floe, and integrate all kinematic degrees of freedom accordingly.  The temporal  integration is explicit and of length simulation.time_step.  This function  will write VTK files to disk in the intervals simulation.file_time_step by the  function writeVTK.  If this value is negative, no output files will be written  to disk.\n\nArguments\n\nsimulation::Simulation: the simulation to run (object is modified)\nverbose::Bool=true: show verbose information during the time loop\nstatus_interval::Bool=true: show verbose information during the time loop\nshow_file_output::Bool=true: report to stdout when output file is written\nsingle_step::Bool=false: run simulation for a single time step only.  If    this causes simulation.time to exceed simulation.time_total, the latter    is increased accordingly.\ntemporal_integration_method::String=\"Three-term Taylor\": type of integration    method to use.  See updateGrainKinematics for details.\nwrite_jld::Bool=false: write simulation state to disk as JLD files (see    Granular.writeSimulation(...) whenever saving VTK output.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setBodyForce!-Tuple{Granular.GrainCylindrical,Array{Float64,1}}",
    "page": "Public API",
    "title": "Granular.setBodyForce!",
    "category": "Method",
    "text": "setBodyForce!(grain, force)\n\nSet the value of the external body force on a grain.\n\nArguments\n\ngrain::GrainCylindrical: the grain to set the body force for.\nforce::Vector{Float64}: a vector of force [N]\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setCurrentTime!-Tuple{Granular.Simulation,Float64}",
    "page": "Public API",
    "title": "Granular.setCurrentTime!",
    "category": "Method",
    "text": "setCurrentTime!(simulation::Simulation, t::Float64)\n\nSets the current simulation time of the simulation object to t, with  parameter value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setGridBoundaryConditions!",
    "page": "Public API",
    "title": "Granular.setGridBoundaryConditions!",
    "category": "Function",
    "text": "setGridBoundaryConditions!(grid, grid_face, mode)\n\nSet boundary conditions for the granular phase at the edges of Ocean or Atmosphere grids.  The target boundary can be selected through the grid_face argument, or the same boundary condition can be applied to all grid boundaries at once.\n\nWhen the center coordinate of grains crosses an inactive boundary (mode = \"inactive\"), the grain is disabled (GrainCylindrical.enabled = false).  This keeps the grain in memory, but stops it from moving or interacting with other grains.  By default, all boundaries are inactive.\n\nIf the center coordinate of a grain crosses a periodic boundary (mode = periodic), the grain is repositioned to the opposite side of the model domain. Grains can interact mechanically across the periodic boundary.\n\nArguments\n\ngrid::Any: Ocean or Atmosphere grid to apply the boundary condition to.\ngrid_face::String: Grid face to apply the boundary condition to.  Valid   values are any combination and sequence of \"west\" (-x), \"south\" (-y),   \"east\" (+x), \"north\" (+y).  The values may be delimited in any way.   Also, and by default, all boundaries can be selected with \"all\" (-x, -y,   +x, +y), which overrides any other face selection.\nmode::String: Boundary behavior, accepted values are \"inactive\",   \"periodic\", and \"impermeable\".  You cannot specify more than one mode at   a time, so if several modes are desired as boundary conditions for the grid,   several calls to this function should be made.\nverbose::Bool: Confirm boundary conditions by reporting values to console.\n\nExamples\n\nSet all boundaries for the ocean grid to be periodic:\n\nsetGridBoundaryConditions!(ocean, \"periodic\", \"all\")\n\nSet the south-north boundaries to be inactive, but the west-east boundaries to be periodic:\n\nsetGridBoundaryConditions!(ocean, \"inactive\", \"south north\")\nsetGridBoundaryConditions!(ocean, \"periodic\", \"west east\")\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setOutputFileInterval!-Tuple{Granular.Simulation,Float64}",
    "page": "Public API",
    "title": "Granular.setOutputFileInterval!",
    "category": "Method",
    "text": "setOutputFileInterval!(simulation::Simulation, t::Float64)\n\nSets the simulation-time interval between output files are written to disk.  If  this value is zero or negative, no output files will be written.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setTimeStep!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.setTimeStep!",
    "category": "Method",
    "text": "Find the computational time step length suitable given the grain radii, contact stiffnesses, and grain density. Uses the scheme by Radjaii et al. 2011.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.setTotalTime!-Tuple{Granular.Simulation,Float64}",
    "page": "Public API",
    "title": "Granular.setTotalTime!",
    "category": "Method",
    "text": "setTotalTime!(simulation::Simulation, t::Float64)\n\nSets the total simulation time of the simulation object to t, with parameter  value checks.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.sortGrainsInGrid!-Tuple{Granular.Simulation,Any}",
    "page": "Public API",
    "title": "Granular.sortGrainsInGrid!",
    "category": "Method",
    "text": "Find grain positions in grid, based on their center positions.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.status",
    "page": "Public API",
    "title": "Granular.status",
    "category": "Function",
    "text": "Shows the status of all simulations with output files written under the  specified folder, which is the current working directory by default.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.totalGrainKineticRotationalEnergy-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.totalGrainKineticRotationalEnergy",
    "category": "Method",
    "text": "totalGrainKineticRotationalEnergy(simulation)\n\nReturns the sum of rotational kinetic energies of all grains in a simulation\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.totalGrainKineticTranslationalEnergy-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.totalGrainKineticTranslationalEnergy",
    "category": "Method",
    "text": "totalGrainKineticTranslationalEnergy(simulation)\n\nReturns the sum of translational kinetic energies of all grains in a  simulation\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.updateGrainKinematics!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.updateGrainKinematics!",
    "category": "Method",
    "text": "updateGrainKinematics!(simulation::Simulation[,\n                         method::String = \"Three-term Taylor\"])\n\nUpdate the grain kinematic parameters using a temporal integration scheme, the current force and torque balance, and gravitational acceleration.  If the simulation contains a grid with periodic boundaries, affected grain positions are adjusted accordingly.\n\nArguments\n\nsimulation::Simulation: update the grain positions in this object    according to temporal integration of length simulation.time_step.\nmethod::String = \"Three-term Taylor\": the integration method to use.     Available methods are \"Two-term Taylor\" and \"Three-term Taylor\".  The    three-term Taylor expansion is slightly more computationally expensive than    the two-term Taylor expansion, but offers an order-of-magnitude increase in    precision of grain positions.  The two-term expansion can obtain similar    precision if the time step is 1/10 the length.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.updateGrainKinematicsThreeTermTaylor!-Tuple{Granular.GrainCylindrical,Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.updateGrainKinematicsThreeTermTaylor!",
    "category": "Method",
    "text": "Use a three-term Taylor expansion for integrating the kinematic degrees of  freedom for an grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.updateGrainKinematicsTwoTermTaylor!-Tuple{Granular.GrainCylindrical,Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.updateGrainKinematicsTwoTermTaylor!",
    "category": "Method",
    "text": "Use a two-term Taylor expansion for integrating the kinematic degrees of freedom  for an grain.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeGrainInteractionVTK-Tuple{Granular.Simulation,String}",
    "page": "Public API",
    "title": "Granular.writeGrainInteractionVTK",
    "category": "Method",
    "text": "writeGrainInteractionVTK(simulation::Simulation,\n                           filename::String;\n                           verbose::Bool=false)\n\nSaves grain interactions to .vtp files for visualization with VTK, for  example in Paraview.  Convert Cell Data to Point Data and use with Tube filter.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeGrainVTK-Tuple{Granular.Simulation,String}",
    "page": "Public API",
    "title": "Granular.writeGrainVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all grains in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.  This function is called by  writeVTK().\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeParaviewPythonScript-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.writeParaviewPythonScript",
    "category": "Method",
    "text": "function writeParaviewPythonScript(simulation,                                    [filename, folder, vtk_folder, verbose])\n\nCreate a \".py\" script for visualizing the simulation VTK files in Paraview. The script can be run from the command line with pvpython (bundled with Paraview), or from the interactive Python shell inside Paraview.\n\nArguments\n\nsimulation::Simulation: input simulation file containing the data.\nfilename::String: output file name for the Python script. At its default   (blank) value, the script is named after the simulation id (simulation.id).\nfolder::String: output directory, current directory the default.\nvtk_folder::String: directory containing the VTK output files, by default   points to the full system path equivalent to \"./<simulation.id>/\".\nsave_animation::Bool: make the generated script immediately save a rendered   animation to disk when the \".py\" script is called.\nverbose::Bool: show diagnostic information during\n\nfunction call, on by     default.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeSimulation-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.writeSimulation",
    "category": "Method",
    "text": "writeSimulation(simulation::Simulation;\n                     filename::String=\"\",\n                     folder::String=\".\",\n                     verbose::Bool=true)\n\nWrite all content from Simulation to disk in JDL format.  If the filename  parameter is not specified, it will be saved to a subdirectory under the current  directory named after the simulation identifier simulation.id.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeSimulationStatus-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.writeSimulationStatus",
    "category": "Method",
    "text": "writeSimulationStatus(simulation::Simulation;\n                      folder::String=\".\",\n                      verbose::Bool=false)\n\nWrite current simulation status to disk in a minimal txt file.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.writeVTK-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.writeVTK",
    "category": "Method",
    "text": "Write a VTK file to disk containing all grains in the simulation in an  unstructured mesh (file type .vtu).  These files can be read by ParaView and  can be visualized by applying a Glyph filter.\n\nIf the simulation contains an Ocean data structure, it's contents will be  written to separate .vtu files.  This can be disabled by setting the argument  ocean=false.  The same is true for the atmosphere.\n\nThe VTK files will be saved in a subfolder named after the simulation.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.zeroForcesAndTorques!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.zeroForcesAndTorques!",
    "category": "Method",
    "text": "Sets the force and torque values of all grains to zero.\n\n\n\n"
},

{
    "location": "lib/public.html#Granular.zeroKinematics!-Tuple{Granular.Simulation}",
    "page": "Public API",
    "title": "Granular.zeroKinematics!",
    "category": "Method",
    "text": "zeroKinematics!(simulation)\n\nSet all grain forces, torques, accelerations, and velocities (linear and rotational) to zero in order to get rid of all kinetic energy.\n\n\n\n"
},

{
    "location": "lib/public.html#Public-Interface-1",
    "page": "Public API",
    "title": "Public Interface",
    "category": "section",
    "text": "Modules = [Granular]\nPublic = true\nPrivate = false"
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
    "text": "This page lists all the documented internals of the Granular module."
},

{
    "location": "lib/internals.html#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals.html#Internal-Interface-1",
    "page": "Internals",
    "title": "Internal Interface",
    "category": "section",
    "text": "`autodocs Modules = Granular Public = false Private = true a"
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
    "text": "This page lists all the documented internals of the Granular module."
},

{
    "location": "lib/internals.html#Index-1",
    "page": "Package-internal documentation",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module.Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals.html#Internal-Interface-1",
    "page": "Package-internal documentation",
    "title": "Internal Interface",
    "category": "section",
    "text": "`autodocs Modules = Granular Public = false Private = true a"
},

]}
