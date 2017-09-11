# SeaIce.jl

*A [Julia](https://julialang.org) package for particle-based simulation of sea-ice dynamics.*

## Package features

Flexible and computationally efficient 2d implementation of the discrete element method.  The particles represent sea-ice floes, which can be forced by ocean and atmospheric velocity fields.  The ice floes can interact through elasto-viscous-frictional contact rheologies and obtain time-dependent tensile strength.

See the [Public API Index](@ref main-index) for the complete list of documented functions and types.

---

### Author
[Anders Damsgaard](https://adamsgaard.dk), Geophysical Fluid Dynamics Laboratory, Princeton University.

### License
SeaIce.jl is licensed under the GPLv3; see [LICENSE](https://github.com/anders-dc/SeaIce.jl/blob/master/LICENSE.md) for the full license text.

## Manual Outline

```@contents
Pages = [
    "man/installation.md",
]
Depth = 1
```

## Library Outline
```@contents
Pages = [
    "lib/public.md",
    "lib/internals.md",
]
Depth = 1
```
