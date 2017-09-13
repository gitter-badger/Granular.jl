#!/usr/bin/env julia

import SeaIce

const verbose = true

const text = "SeaIce.jl"

const forcing = "gyres"
#const forcing = "down"
#const forcing = "convergent"

# Font created with `figlet` and the font 'pebbles'.  If figlet is not installed 
# on your system, use the string below:
#logo_string =
#""".oOOOo.               ooOoOOo                     o 
#o     o                  O                     O O  
#O.                       o                       o  
# `OOoo.                  O                       O  
#      `O .oOo. .oOoO'    o    .oOo  .oOo.     'o o  
#       o OooO' O   o     O    O     OooO'      O O  
#O.    .O O     o   O     O    o     O     oO   o o  
# `oooO'  `OoO' `OoO'o ooOOoOo `OoO' `OoO' Oo   O Oo 
#                                               o    
#                                             oO'    """

logo_string = readstring(`figlet -f pebbles "$text"`)

const dx = 1.
const dy = dx

const logo_string_split = split(logo_string, '\n')

const ny = length(logo_string_split)
maxwidth = 0
for i=1:ny
    if maxwidth < length(logo_string_split[i])
        maxwidth = length(logo_string_split[i])
    end
end
const nx = maxwidth + 1

const Lx = nx*dx
const Ly = ny*dy

x = 0.
y = 0.
r = 0.
c = ' '
h = .5
const youngs_modulus = 2e6

sim = SeaIce.createSimulation(id="logo")

print(logo_string)
info("nx = $nx, ny = $ny")

for iy=1:length(logo_string_split)
    for ix=1:length(logo_string_split[iy])

        c = logo_string_split[iy][ix]

        if c == ' '
            continue
        elseif c == 'O'
            x = ix*dx - .5*dx
            y = Ly - (iy*dy - .5*dy)
            r = .5*dx
        elseif c == 'o'
            x = ix*dx - .5*dx
            y = Ly - (iy*dy - .33*dy)
            r = .33*dx
        elseif c == 'o'
            x = ix*dx - .5*dx
            y = Ly - (iy*dy - .25*dy)
            r = .25*dx
        elseif c == '\''
            x = ix*dx - .75*dx
            y = Ly - (iy*dy - .75*dy)
            r = .25*dx
        elseif c == '`'
            x = ix*dx - .25*dx
            y = Ly - (iy*dy - .75*dy)
            r = .25*dx
        end

        if r > 0.
            SeaIce.addIceFloeCylindrical!(sim, [x + dx, y - dy], r, h,
                                          tensile_strength=200e3,
                                          youngs_modulus=youngs_modulus,
                                          verbose=verbose)
        end
        r = -1.
    end
end

# set ocean forcing
sim.ocean = SeaIce.createRegularOceanGrid([nx, ny, 1], [Lx, Ly, 1.],
                                          name="logo_ocean")

if forcing == "gyres"
    epsilon = 0.25  # amplitude of periodic oscillations
    t = 0.
    a = epsilon*sin(2.*pi*t)
    b = 1. - 2.*epsilon*sin(2.*pi*t)
    for i=1:size(sim.ocean.u, 1)
        for j=1:size(sim.ocean.u, 2)

            x = sim.ocean.xq[i, j]/(Lx*.5)  # x in [0;2]
            y = sim.ocean.yq[i, j]/Ly       # y in [0;1]

            f = a*x^2. + b*x
            df_dx = 2.*a*x + b

            sim.ocean.u[i, j, 1, 1] = -pi/10.*sin(pi*f)*cos(pi*y) * 2e1
            sim.ocean.v[i, j, 1, 1] = pi/10.*cos(pi*f)*sin(pi*y)*df_dx * 2e1
        end
    end

elseif forcing == "down"
    Base.Random.srand(1)
    sim.ocean.u[:, :, 1, 1] = (Base.Random.rand(nx+1, ny+1) - .5)*.1
    sim.ocean.v[:, :, 1, 1] = -5.

elseif forcing == "convergent"
    Base.Random.srand(1)
    sim.ocean.u[:, :, 1, 1] = (Base.Random.rand(nx+1, ny+1) - .5)*.1
    for j=1:size(sim.ocean.u, 2)
        sim.ocean.v[:, j, 1, 1] = -(j/ny - .5)*10.
    end

else
    error("Forcing not understood")
end

# Initialize confining walls, which are ice floes that are fixed in space
r = dx/4.

## N-S wall segments
for y in linspace(r, Ly-r, Int(round((Ly - 2.*r)/(r*2))))
    SeaIce.addIceFloeCylindrical!(sim, [r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    SeaIce.addIceFloeCylindrical!(sim, [Lx-r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end

## E-W wall segments
for x in linspace(3.*r, Lx-3.*r, Int(round((Lx - 6.*r)/(r*2))))
    SeaIce.addIceFloeCylindrical!(sim, [x, r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    SeaIce.addIceFloeCylindrical!(sim, [x, Ly-r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end


# Finalize setup and start simulation
SeaIce.setTimeStep!(sim, verbose=verbose)

SeaIce.setTotalTime!(sim, 5.)
SeaIce.setOutputFileInterval!(sim, .1)

SeaIce.removeSimulationFiles(sim)
SeaIce.run!(sim, verbose=verbose)
