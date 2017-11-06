#/usr/bin/env julia
import Granular

################################################################################
#### Step 1: Create a loose granular assemblage and let it settle at -y        #
################################################################################
sim = Granular.createSimulation(id="shear-init")

# Generate 10 grains along x and 100 grains along y, with radii between 0.2 and
# 1.0 m.
Granular.regularPacking!(sim, [10, 100], 0.2, 1.0)

# Create a grid for contact searching spanning the extent of the grains
Granular.fitGridToGrains!(sim, sim.ocean)

# Make the ocean grid drag grains uniformly towards -y
sim.ocean.v[:, :, 1, 1] = -5.0

# Make the top and bottom boundaries impermeable, and the side boundaries
# periodic, which will come in handy during shear
Granular.setGridBoundaryConditions!(sim.ocean, "impermeable", "north south")
Granular.setGridBoundaryConditions!(sim.ocean, "periodic", "east west")

# Automatically set the computational time step based on grain sizes and
# properties
Granular.setTimeStep!(sim)

# Set the total simulation time for this step [s]
Granular.setTotalTime!(sim, 30.)

# Set the interval in model time between simulation files [s]
Granular.setOutputFileInterval!(sim, .2)

# Visualize the grain-size distribution
Granular.plotGrainSizeDistribution(sim)

# Start the simulation
Granular.run!(sim)

# Try to render the simulation if `pvpython` is installed on the system
Granular.render(sim)










################################################################################
#### Step 2: Consolidate the previous output under a constant normal stress    #
################################################################################

















################################################################################
#### Step 3: Shear the consolidated assemblage with a constant velocity        #
################################################################################




















