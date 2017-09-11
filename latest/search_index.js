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
    "text": "Flexible and computationally efficient 2d implementation of the discrete element method.  The particles represent sea-ice floes, which can be forced by ocean and velocity fields.  The ice floes can interact through elasto-viscous-frictional contact rheologies and obtain time-dependent tensile strength."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"installation.md\"\n]\nDepth = 1"
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

]}
