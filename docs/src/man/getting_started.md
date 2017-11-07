# Getting started
If you are new to the [Julia](https://julialang.org) programming language, the 
official manual has a useful guide to [getting started with 
Julia](https://docs.julialang.org/en/latest/manual/getting-started/).

In the following, two simple examples are presented using some of the core 
commands of Granular.jl.  For more examples, see the scripts in the `examples/` 
directory.

The relevant functions are all contained in the `Granular` module, which can be 
imported with `import Granular` at the top of your script.  *Note:* As per 
Julia conventions, functions that contain an exclamation mark (!) modify the 
values of the arguments.

Any of the functions called below are documented in the source code, and this 
documentation can be found in the [Public API Index](@ref main-index) in the 
online documentation, or simply from the Julia shell by typing `?<function 
name>`.  An example:

```julia-repl
julia> ?Granular.fitGridToGrains!

  fitGridToGrains!(simulation, grid[, padding])

  Fit the ocean or atmosphere grid for a simulation to the current grains and their positions.

     Arguments
    ≡≡≡≡≡≡≡≡≡≡≡

    •    simulation::Simulation: simulation object to manipulate.

    •    grid::Any: Ocean or Atmosphere grid to manipulate.

    •    padding::Real: optional padding around edges [m].

    •    verbose::Bool: show grid information when function completes.

```

## Collision between two particles
For this simple example (`example/two-grains.jl`), we will create two grains, 
where one of the grains is bumping in to the other.

As the first command, we import all the Granular.jl functionality:

```julia-repl
julia> import Granular
```

### Simulation setup
Next, we create our simulation object which holds all information on the 
simulated grains.  This object can be called whatever is appropriate.  In this 
documentation, we will use the name `sim`:

```julia-repl
julia> sim = Granular.createSimulation(id="two-grains")
Granular.Simulation("two-grains", 0, 0.0, -1.0, -1.0, -1.0, 0, 0.0, 
Granular.GrainCylindrical[], Granular.Ocean(false, [0.0], [0.0], [0.0], [0.0], 
[0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], Array{Int64,1}[#undef], 1, 1, 
1, 1), Granular.Atmosphere(false, [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], 
[0.0], [0.0], Array{Int64,1}[#undef], 1, 1, 1, 1, false), 16)
```

We will be presented with some output about the contents of the `sim` 
simulation object.  This is of minor importance as of now, and can safely be 
ignored.

During the above `createSimulation` call, the `id` argument is optional, but is 
used to name simulation output files that are written to the disk.  It is good 
practice to use the same name as used for the simulation script file.

### Adding grains one by one
We have now created a simulation object, which will be used during all of the 
following.  Next, we add grains to this object.  The first grain is cylinder 
shaped, placed at the x-y position (0,0) m, has a radius of 0.1 m, and a 
thickness of 0.05 m (along z).  As this call modifies the `sim` object, the 
function contains an exclamation mark (!).  For further information regarding 
this call, see the reference in the [Public API Index](@ref main-index).

```julia-repl
julia> Granular.addGrainCylindrical!(sim, [0.0, 0.0], 0.1, 0.05)
INFO: Added Grain 1
```

Let's add another grain, placed at some distance from the first grain:

```julia-repl
julia> Granular.addGrainCylindrical!(sim, [0.5, 0.0], 0.1, 0.05)
INFO: Added Grain 2
```

We now want to prescribe a linear (not rotational or angular) velocity to the 
first grain, to make it bump into the second grain.

The simulation object `sim` contains an array of all grains that are added to 
it.  We can directly inspect the grains and their values from the simulation 
object.  Let's take a look at the default value of the linear velocity, called 
`lin_vel`:

```julia-repl
julia> sim.grains[1].lin_vel
2-element Array{Float64, 1}:
 0.0
 0.0
```

The default value is a 0,0 vector.  Similarly, we can modify the properties of 
the first grain directly with the following call:

```julia-repl
julia> sim.grains[1].lin_vel = [1.0, 0.0]
2-element Array{Float64, 1}:
 1.0
 0.0
```

The first grain (index 1 in `sim.grains`) now has a positive velocity along `x` 
with the value of 1.0 meter per second.

### Setting temporal parameters for the simulation
Before we can start the simulation, we need to tell the code vital information, 
like what time step to use, how often to write output files to the disk, and 
for how long to run the simulation.  To set the computational time step, we 
call the following:

```julia-repl
julia> Granular.setTimeStep!(sim)
INFO: Time step length t=8.478741928617433e-5 s
```

Based on the elastic stiffness of the grains and their size, a suitable time 
step is automatically determined.

The two grains are 0.3 meters apart, as the centers are placed 0.5 meter from 
each other and each grain has a diameter of 0.1 m.  With a linear velocity of 
1.0 m/s, the collision should occur after 0.3 seconds.  For that reason it 
seems reasonable to run the simulation for a total of 1.0 seconds, and produce 
an output file every 0.05 seconds.  We will later use the output files for 
visualization purposes.

```julia-repl
julia> Granular.setOutputFileInterval!(sim, 0.05)

julia> Granular.setTotalTime!(sim, 1.0)
```

### Running the simulation
We are now ready to run the simulation.  For illustrative purposes, let us 
compare the kinetic energy in the granular system before and after the 
collision.  For now, we save the initial value using the following call:

```julia-repl
julia> Granular.totalGrainKineticTranslationalEnergy(sim)
0.7335618846132168
```

The above value is the total translational (not angular) kinetic energy in 
Joules before the simulation is started.

For running the simulation, we have two choices; we can either run the entire 
simulation length with a single call, which steps time until the total time is 
reached and generates output files along the way.  Alternatively, we can run 
the simulation for a single time step a time, and inspect the progress or do 
other modifications along the way.

Here, we will run the entire simulation in one go, and afterwards visualize the 
grains from their output files using ParaView.

```julia-repl
julia> Granular.run!(sim)

INFO: Output file: ./two-grains/two-grains.grains.1.vtu
INFO: wrote status to ./two-grains/two-grains.status.txt
  t = 0.04239370964308682/1.0 s
INFO: Output file: ./two-grains/two-grains.grains.2.vtu
INFO: wrote status to ./two-grains/two-grains.status.txt

...

INFO: Output file: ./two-grains/two-grains.grains.20.vtu
INFO: wrote status to ./two-grains/two-grains.status.txt
  t = 0.9920128056483803/1.0 s
INFO: ./two-grains/two-grains.py written, execute with 'pvpython /Users/ad/code/Granular-ext/two-grains/two-grains.py'
INFO: wrote status to ./two-grains/two-grains.status.txt
  t = 1.0000676104805686/1.0 s
```

The script informs us of the simulation progress, and notifies us as output 
files are generated (this can be disabled by passing `verbose=false` to the 
`run!()` command).  Finally, it tells us that it has generated a ParaView 
python file for visualization.

Before going further, we are interested in getting an immediate idea of how the 
collision went.  We print the new velocities with the following commands:

```julia-repl
julia> sim.grains[1].lin_vel
2-element Array{Float64, 1}:
 7.58343e-5
 0.0

julia> sim.grains[2].lin_vel
2-element Array{Float64, 1}:
 0.999924
 0.0
```

The first grain has transferred effectively all of its kinetic energy to the 
second grain during the cause of the simulation.  The total kinetic energy now 
is the following:

```julia-repl
julia> Granular.totalGrainKineticTranslationalEnergy(sim)
0.7334506347624973
```

The before and after values are reasonably close (to less than 0.1 percent), 
which is what can be expected given the computation accuracy in the algorithm.

### Visualizing the output
To visualize the output we open [ParaView](https://www.paraview.org).  The 
output files of the simulation are written using the VTK (visualization 
toolkit) format, which is natively supported by ParaView.

While the `.vtu` files produced during the simulation can be opened with 
ParaView and visualized manually using *Glyph* filters, the simplest and 
fastest way to visualize the data is to use the Python script generated for the 
simulation by Granular.jl.

Open ParaView and open the *Python Shell*, found under the menu *Tools > Python 
Shell*.  In the pop-up dialog we select *Run Script*, which opens yet another 
dialog prompting us to locate the visualization script (`two-grains.py`, in our 
example).  We locate this file, which is placed under the directory from where 
we launched the `julia` session with the commands above.

After selecting the `two-grains/two-grains.py` script, we can close the *Python 
Shell* window to inspect our simulation.  Press the *Play* symbol in the top 
toolbar, and see what happens!

Alternatively, you can color the grains using different parameters, such as 
velocity, number of contacts, etc.  These can be selected by changing the 
chosen parameter under the *Glyph1* object in the *Pipeline Browser* on the 
left, and selecting a different field for *Coloring*.  Press the *Apply* button 
to see the changes in effect.

### Exercises
To gain more familiarity with the simulation procedure, I suggest experimenting 
with the following:

- What effect does the grain size have on the time step?
- Try to make an oblique collision by placing one of the grains at a different 
    `y` position.
- What happens if the second grains is set to be fixed in space 
    (`sim.grains[2].fixed = true`)?
- How is the relationship between total kinetic energy before and after 
    affected by the choice of time step length?  Try setting different time 
    step values, e.g. with `sim.time_step = 0.1234` and rerun the simulation.

