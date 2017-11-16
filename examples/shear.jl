#/usr/bin/env julia
import Granular
import JLD
ENV["MPLBACKEND"] = "Agg"
import PyPlot

# Common simulation identifier
const id_prefix = "test0"

################################################################################
#### Step 1: Create a loose granular assemblage and let it settle at -y        #
################################################################################
sim = Granular.createSimulation(id="$(id_prefix)-init")

# Generate 10 grains along x and 50 grains along y, with radii between 0.2 and
# 1.0 m.
const nx = 10
const ny = 50
const r_min = 0.2
const r_max = 1.0
Granular.regularPacking!(sim, [nx, ny], r_min, r_max)

# Create a grid for contact searching spanning the extent of the grains
Granular.fitGridToGrains!(sim, sim.ocean)

# Make the top and bottom boundaries impermeable, and the side boundaries
# periodic, which will come in handy during shear
Granular.setGridBoundaryConditions!(sim.ocean, "impermeable", "north south")
Granular.setGridBoundaryConditions!(sim.ocean, "periodic", "east west")

# Add gravitational acceleration to all grains and disable ocean-grid drag
const g = [0., -9.8]
for grain in sim.grains
    Granular.addBodyForce!(grain, grain.mass*g)
    Granular.disableOceanDrag!(grain)
end

# Automatically set the computational time step based on grain sizes and
# properties
Granular.setTimeStep!(sim)

# Set the total simulation time for this step [s]
Granular.setTotalTime!(sim, 30.)

# Set the interval in model time between simulation files [s]
Granular.setOutputFileInterval!(sim, .2)

# Visualize the grain-size distribution
#Granular.plotGrainSizeDistribution(sim)

# Start the simulation
Granular.run!(sim)

# Try to render the simulation if `pvpython` is installed on the system
Granular.render(sim, trim=false)

# Save the simulation state to disk in case we need to reuse the current state
# This step requires the JLD package (Pkg.add("JLD"))
Granular.writeSimulation(sim)

# Also copy the simulation in memory, in case we want to loop over different
# normal stresses below:
sim_init = deepcopy(sim)


################################################################################
#### Step 2: Consolidate the previous output under a constant normal stress    #
################################################################################

# Select a normal stress for the consolidation [Pa]
N = 10e3

# Rename the simulation so we don't overwrite output from the previous step
sim.id = "$(id_prefix)-cons-N$(N)Pa"

# Set all linear and rotational velocities to zero
Granular.zeroKinematics!(sim)

# Add a dynamic wall to the top which adds a normal stress downwards.  The
# normal of this wall is downwards, and we place it at the top of the granular
# assemblage
y_top = -Inf
for grain in grains
    if y_top < grain.lin_pos[2] + grain.contact_radius
        y_top = grain.lin_pos[2] + grain.contact_radius
    end
end
Granular.addWallLinearFrictionless!(sim, normal=[0., 1.], pos=y_top,
                                    bc="normal stress", normal_stress=-N)

# Resize the grid to span the current state
Granular.fitGridToGrains!(sim, sim.ocean)

# Lock the grains at the very bottom so that the lower boundary is rough
y_bot = Inf
for grain in grains
    if y_bot > grain.lin_pos[2] - grain.contact_radius
        y_bot = grain.lin_pos[2] - grain.contact_radius
    end
end
const fixed_thickness = 2. * r_max
for grain in grains
    if grain.lin_pos[2] <= fixed_thickness
        grain.fixed = true  # set x and y acceleration to zero
    end
end

# Set current time to zero and reset output file counter
Granular.resetTime!(sim)

# Set the simulation time to run the consolidation for
Granular.setTotalTime!(sim, 5.0)

# Run the consolidation experiment
Granular.run!(sim)

# Try to render the simulation if `pvpython` is installed on the system
Granular.render(sim, trim=false)

# Save the simulation state to disk in case we need to reuse the consolidated
# state (e.g. different shear velocities below)
Granular.writeSimulation(sim)

# Also copy the simulation in memory, in case we want to loop over different
# normal stresses below:
sim_cons = deepcopy(sim)


################################################################################
#### Step 3: Shear the consolidated assemblage with a constant velocity        #
################################################################################

# Select a shear velocity for the consolidation [m/s]
const vel_shear = 0.1

# Rename the simulation so we don't overwrite output from the previous step
sim.id = "$(id_prefix)-shear-N$(N)Pa-vel_shear$(vel_shear)m-s"

# Set all linear and rotational velocities to zero
Granular.zeroKinematics!(sim)

# Prescribe the shear velocity to the uppermost grains
for grain in grains
    if grain.lin_pos[2] <= fixed_thickness

        # do not allow changes in velocity
        grain.fixed = true
        
        # allow free up/down movement to permit dilation, which partially
        # overrides the `fixed` flag
        grain.allow_y_acc = true

        grain.lin_vel[1] = vel_shear
    end
end

# Set current time to zero and reset output file counter
Granular.resetTime!(sim)

# Set the simulation time to run the shear experiment for
Granular.setTotalTime!(sim, 15.0)

# Run the shear experiment
time = Float64[]
shear_stress = Float64[]
shear_strain = Float64[]
dilation = Float64[]
const thickness_initial = sim.walls[1].pos - y_bot
x_min = +Inf
x_max = -Inf
for grain in grains
    if x_min > grain.lin_pos[1] - grain.contact_radius
        x_min = grain.lin_pos[1] - grain.contact_radius
    end
    if x_max < grain.lin_pos[1] + grain.contact_radius
        x_max = grain.lin_pos[1] + grain.contact_radius
    end
end
const surface_area = (x_max - x_min)
while sim.time < sim.time_total

    for i=1:100  # run for 100 steps before measuring shear stress and dilation
        Granular.run!(sim, single_step=true)
    end

    append!(time, sim.time)

    # Determine the current shear stress
    for grain in sim.grains
        if grain.fixed && grain.allow_y_acc
            shear_force += -grain.force[1]
        end
    end
    append!(shear_stress, shear_force/surface_area)

    # Determine the current shear strain
    append!(shear_strain, sim.time*vel_shear/thickness_initial)

    # Determine the current dilation
    append!(dilation, (sim.walls[1].pos - y_bot)/thickness_initial)

end

# Try to render the simulation if `pvpython` is installed on the system
Granular.render(sim, trim=false)

# Save the simulation state to disk in case we need to reuse the sheared state
Granular.writeSimulation(sim)

# Plot time vs. shear stress and dilation
PyPlot.subplot(211)
PyPlot.plot(time, shear_stress)
PyPlot.subplot(212)
PyPlot.plot(time, dilation)
PyPlot.xlabel("Time [s]")
PyPlot.ylabel("Shear stress [Pa]")
PyPlot.savefig(sim.id * "-time_vs_shear-stress.pdf")
PyPlot.clf()

# Plot shear strain vs. shear stress and dilation
PyPlot.subplot(211)
PyPlot.plot(time, shear_stress)
PyPlot.subplot(212)
PyPlot.plot(shear_strain, dilation)
PyPlot.xlabel("Shear strain [-]")
PyPlot.ylabel("Shear stress [Pa]")
PyPlot.savefig(sim.id * "-shear-strain_vs_shear-stress.pdf")
PyPlot.clf()

# Plot time vs. shear strain (boring when the shear experiment is rate controlled)
PyPlot.plot(time, shear_strain)
PyPlot.xlabel("Time [s]")
PyPlot.ylabel("Shear strain [-]")
PyPlot.savefig(sim.id * "-time_vs_shear-strain.pdf")
PyPlot.clf()
