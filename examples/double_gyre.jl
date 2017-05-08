#!/usr/bin/env julia
import SeaIce

sim = SeaIce.createSimulation(id="double_gyre")

# Initialize ocean
L = [100e3, 50e3, 1e3]
Ly_constriction = 20e3
#n = [750, 500, 2]  # high resolution
n = [30, 15, 2]  # intermedite resolution
#n = [8, 5, 2]  # coarse resolution
sim.ocean = SeaIce.createRegularOceanGrid(n, L, name="double_gyre")

epsilon = 0.25  # amplitude of periodic oscillations
t = 0.
a = epsilon*sin(2.*pi*t)
b = 1. - 2.*epsilon*sin(2.*pi*t)
for i=1:size(sim.ocean.u, 1)
    for j=1:size(sim.ocean.u, 2)

        x = sim.ocean.xq[i, j]/(L[1]*.5)  # x in [0;2]
        y = sim.ocean.yq[i, j]/L[2]       # y in [0;1]

        f = a*x^2. + b*x
        df_dx = 2.*a*x + b

        sim.ocean.u[i, j, 1, 1] = -pi/10.*sin(pi*f)*cos(pi*y) * 1e1
        sim.ocean.v[i, j, 1, 1] = pi/10.*cos(pi*f)*sin(pi*y)*df_dx * 1e1
    end
end

# Initialize confining walls, which are ice floes that are fixed in space
r = minimum(L[1:2]/n[1:2])/2.
h = 1.

## N-S wall segments
for y in linspace(r, L[2]-r, Int(round((L[2] - 2.*r)/(r*2))))
    SeaIce.addIceFloeCylindrical(sim, [r, y], r, h, fixed=true,
                                 verbose=false)
    SeaIce.addIceFloeCylindrical(sim, [L[1]-r, y], r, h, fixed=true,
                                 verbose=false)
end

## E-W wall segments
for x in linspace(3.*r, L[1]-3.*r, Int(round((L[1] - 6.*r)/(r*2))))
    SeaIce.addIceFloeCylindrical(sim, [x, r], r, h, fixed=true,
                                 verbose=false)
    SeaIce.addIceFloeCylindrical(sim, [x, L[2]-r], r, h, fixed=true,
                                 verbose=false)
end

n_walls = length(sim.ice_floes)
info("added $(n_walls) fixed ice floes as walls")



# Initialize ice floes everywhere
floe_padding = .5*r
noise_amplitude = .8*floe_padding
Base.Random.srand(1)
for y in (4.*r + noise_amplitude):(2.*r + floe_padding):(L[2] - 4.*r - 
                                                         noise_amplitude)
                                                         
    for x in (4.*r + noise_amplitude):(2.*r + floe_padding):(L[1] - 4.*r - 
                                                             noise_amplitude)
        #if iy % 2 == 0
            #x += 1.5*r
        #end

        x_ = x + noise_amplitude*(0.5 - Base.Random.rand())
        y_ = y + noise_amplitude*(0.5 - Base.Random.rand())

        SeaIce.addIceFloeCylindrical(sim, [x_, y_], r, h, verbose=false)
    end
end
n = length(sim.ice_floes) - n_walls
info("added $(n) ice floes")

# Remove old simulation files
SeaIce.removeSimulationFiles(sim)

k_n = 1e6  # N/m
gamma_t = 1e7  # N/(m/s)
mu_d = 0.7
rotating = false
for i=1:length(sim.ice_floes)
    sim.ice_floes[i].contact_stiffness_normal = k_n
    sim.ice_floes[i].contact_stiffness_tangential = k_n
    sim.ice_floes[i].contact_viscosity_tangential = gamma_t
    sim.ice_floes[i].contact_dynamic_friction = mu_d
    sim.ice_floes[i].rotating = rotating
end

# Set temporal parameters
SeaIce.setTotalTime!(sim, 12.*60.*60.)
SeaIce.setOutputFileInterval!(sim, 60.)
SeaIce.setTimeStep!(sim)

SeaIce.run!(sim, status_interval=1,
            contact_tangential_rheology="Linear Viscous Frictional")
