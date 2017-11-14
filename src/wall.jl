## Manage dynamic walls in the model

export addWallLinearFrictionless!
"""
    function addWallLinear!(simulation, normal, pos[, bc, mass, thickness, 
                            normal_stress, vel, force, verbose])

Creates and adds a linear (flat) and frictionless dynamic wall to a grain to a
simulation. Most of the arguments are optional, and come with default values.
The only required arguments are 
`simulation`, `normal`, `pos`, and `bc`.

# Arguments
* `simulation::Simulation`: the simulation object where the wall should be
    added to.
* `normal::Vector{Float64}`: 2d vector denoting the normal to the wall [m].
* `pos::Float64`: position along axis parallel to the normal vector [m].
* `bc::String="fixed"`: boundary condition, possible values are `"fixed"`
    (default), `"normal stress"`, or `"velocity"`.
* `mass::Float64=NaN`: wall mass, which is used if wall boundary conditions
    differs from `bc="fixed"`.  If the parameter is left to its default value,
    the wall mass is set to be equal the total mass of grains in the simulation.
    Units: [kg]
* `thickness::Float64=NaN`: wall thickness, which is used for determining wall
    surface area.  If the parameter is left to its default value, the wall
    thickness is set to be equal to the thickest grain in the simulation.
    Units: [m].
* `normal_stress::Float64=0.`: the wall normal stress when `bc == "normal
    stress"` [Pa].
* `vel::Float64=0.`: the wall velocity along the `normal` vector.  If the
    wall boundary condition is `bc = "velocity"` the wall will move according to
    this constant value.  If `bc = "normal stress"` the velocity will be a free
    parameter. Units: [m/s]
* `force::Float64=0.`: sum of normal forces on the wall from interaction with
    grains [N].
* `verbose::Bool=true`: show verbose information during function call.

# Examples
The most basic example adds a new fixed wall to the simulation `sim`, with a 
wall-face normal of `[1., 0.]` (wall along *y* and normal to *x*), a position of
`1.5` meter:

```julia
Granular.addWallLinearFrictionless!(sim, [1., 0.], 1.5)
```

The following example creates a wall with a velocity of 0.5 m/s towards *-y*:

```julia
Granular.addWallLinearFrictionless!(sim, [0., 1.], 1.5,
                                    bc="velocity",
                                    vel=-0.5)
```

To create a wall parallel to the *x* axis with a constant normal stress of 100
kPa:

```julia
Granular.addWallLinearFrictionless!(sim, [1., 0.], 3.5,
                                    bc="normal stress",
                                    normal_stress=100e3)
```
"""
function addWallLinearFrictionless!(simulation::Simulation,
                                    normal::Vector{Float64},
                                    pos::Float64;
                                    bc::String = "fixed",
                                    mass::Float64 = NaN,
                                    thickness::Float64 = NaN,
                                    normal_stress::Float64 = 0.,
                                    vel::Float64 = 0.,
                                    force::Float64 = 0.,
                                    verbose::Bool=true)

    # Check input values
    if length(normal) != 2
        error("Wall normal must be a two-element array (normal = ",
              "$normal)")
    end

    if !(normal ≈ [1., 0.]) && !(normal ≈ [0., 1.])
        error("Currently only walls with normals orthogonal to the " *
              "coordinate system are allowed, i.e. normals parallel to the " *
              "x or y axes.  Accepted values for `normal` " *
              "are [1., 0.] and [0., 1.].  The passed normal was $normal")
    end

    # if not set, set wall mass to equal the mass of all grains.
    if isnan(mass)
        if length(simulation.grains < 1)
            error("If wall mass is not specified, walls should be added " *
                  "after grains have been added to the simulation.")
        end
        mass = 0.
        for grain in simulation.grains
            mass += grain.mass
        end
        info("Setting wall mass to total grain mass: $mass kg")
    end

    # if not set, set wall thickness to equal largest grain thickness
    if isnan(thickness)
        if length(simulation.grains < 1)
            error("If wall thickness is not specified, walls should be added " *
                  "after grains have been added to the simulation.")
        end
        thickness = -Inf
        for grain in simulation.grains
            if grain.thickess > thickness
                thickness = grain.thickness
            end
        end
        info("Setting wall thickness to largest grain thickness: $thickness m")
    end

    # Create wall object
    wall = WallLinearFrictionless(normal,
                                   pos,
                                   bc,
                                   mass,
                                   thickness,
                                   normal_stress,
                                   vel,
                                   force)

    # Add to simulation object
    addWall!(simulation, wall, verbose)
    nothing
end

