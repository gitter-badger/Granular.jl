#!/usr/bin/env julia

import Granular
import FileIO
import Colors

const verbose = true

const img_file = "aadcroft.png"

img = FileIO.load(img_file)

# resize the image if it is too large, preceed with lopass to avoid antialias
max_pixels = 100^2
if size(img, 1)*size(img, 2) > max_pixels
    cp(img_file, "backup-" * img_file, remove_destination=true)
    run(`convert $(img_file) -resize "$(max_pixels)@>" $(img_file)`)
    img = FileIO.load(img_file)
end

const img_bw = Colors.Gray.(img)

const forcing = "gyres"
#const forcing = "down"
#const forcing = "convergent"
#const forcing = "sandpile"

const dx = 1.
const dy = dx

const nx = size(img_bw, 2) + 1
const ny = size(img_bw, 1) + 1

Lx = nx*dx
if forcing == "sandpile"
    Lx *= 1.5
end
const Ly = ny*dy

const youngs_modulus = 2e7
const tensile_strength = 0e3
const h = .5

sim = Granular.createSimulation(id="image")

info("nx = $nx, ny = $ny")

for iy=1:size(img_bw, 1)
    for ix=1:size(img_bw, 2)

        x = ix*dx - dx
        if forcing == "sandpile"
            x += Lx/6.
        end
        y = Ly - (iy*dy - dy)
        r = .5*dx*((1. - Float64(img_bw[iy, ix])))

        if r > .1*dx
            Granular.addGrainCylindrical!(sim, [x + dx, y - dy], r, h,
                                          tensile_strength=tensile_strength,
                                          youngs_modulus=youngs_modulus,
                                          verbose=verbose)
        end
    end
end

# set ocean forcing
sim.ocean = Granular.createRegularOceanGrid([nx, ny, 1], [Lx, Ly, 1.],
                                          name="image_ocean")

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

            sim.ocean.u[i, j, 1, 1] = -pi/10.*sin(pi*f)*cos(pi*y) * 4e1
            sim.ocean.v[i, j, 1, 1] = pi/10.*cos(pi*f)*sin(pi*y)*df_dx * 4e1
        end
    end

elseif forcing == "down" || forcing == "sandpile"
    srand(1)
    sim.ocean.u[:, :, 1, 1] = (rand(nx+1, ny+1) - .5)*.1
    sim.ocean.v[:, :, 1, 1] = -Ly/5.

elseif forcing == "convergent"
    srand(1)
    sim.ocean.u[:, :, 1, 1] = (rand(nx+1, ny+1) - .5)*.1
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
    Granular.addGrainCylindrical!(sim, [r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    Granular.addGrainCylindrical!(sim, [Lx-r, y], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end

## E-W wall segments
for x in linspace(3.*r, Lx-3.*r, Int(round((Lx - 6.*r)/(r*2))))
    Granular.addGrainCylindrical!(sim, [x, r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
    Granular.addGrainCylindrical!(sim, [x, Ly-r], r, h, fixed=true,
                                  youngs_modulus=youngs_modulus,
                                  verbose=false)
end


# Finalize setup and start simulation
Granular.setTimeStep!(sim, verbose=verbose)

Granular.setTotalTime!(sim, 5.)
Granular.setOutputFileInterval!(sim, .1)

Granular.removeSimulationFiles(sim)
Granular.run!(sim, verbose=verbose)
