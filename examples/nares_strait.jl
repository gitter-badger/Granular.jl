#!/usr/bin/env julia
import SeaIce

sim = SeaIce.createSimulation(id="nares_strait")
n = [25, 25, 2]

#sim = SeaIce.createSimulation(id="nares_strait_coarse")
#n = [6, 6, 2]

# Initialize ocean
Lx = 50.e3
Lx_constriction = 10e3
L = [Lx, Lx*1.5, 1e3]
Ly_constriction = 20e3
sim.ocean = SeaIce.createRegularOceanGrid(n, L, name="poiseuille_flow")
sim.ocean.v[:, :, 1, 1] = 1e-8*((sim.ocean.xq - Lx/2.).^2 - Lx^2./4.)

# Initialize confining walls, which are ice floes that are fixed in space
r = minimum(L[1:2]/n[1:2])/2.
h = 1.

## N-S segments
for y in linspace((L[2] - Ly_constriction)/2.,
                  Ly_constriction + (L[2] - Ly_constriction)/2., 
                  Int(round(Ly_constriction/(r*2))))
    SeaIce.addIceFloeCylindrical(sim, [(Lx - Lx_constriction)/2., y], r, h, 
                                 fixed=true, verbose=false)
end
for y in linspace((L[2] - Ly_constriction)/2.,
                  Ly_constriction + (L[2] - Ly_constriction)/2., 
                  Int(round(Ly_constriction/(r*2))))
    SeaIce.addIceFloeCylindrical(sim,
                                 [Lx_constriction + (L[1] - Lx_constriction)/2., 
                                  y], r, h, fixed=true, verbose=false)
end

dx = 2.*r*sin(atan((Lx - Lx_constriction)/(L[2] - Ly_constriction)))

## NW diagonal
x = r:dx:((Lx - Lx_constriction)/2.)
y = linspace(L[2] - r, (L[2] - Ly_constriction)/2. + Ly_constriction + r, 
             length(x))
for i in 1:length(x)
    SeaIce.addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, 
                                 verbose=false)
end

## NE diagonal
x = (L[1] - r):(-dx):((Lx - Lx_constriction)/2. + Lx_constriction)
y = linspace(L[2] - r, (L[2] - Ly_constriction)/2. + Ly_constriction + r, 
             length(x))
for i in 1:length(x)
    SeaIce.addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, 
                                 verbose=false)
end

## SW diagonal
x = r:dx:((Lx - Lx_constriction)/2.)
y = linspace(r, (L[2] - Ly_constriction)/2. - r, length(x))
for i in 1:length(x)
    SeaIce.addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, 
                                 verbose=false)
end

## SE diagonal
x = (L[1] - r):(-dx):((Lx - Lx_constriction)/2. + Lx_constriction)
y = linspace(r, (L[2] - Ly_constriction)/2. - r, length(x))
for i in 1:length(x)
    SeaIce.addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, 
                                 verbose=false)
end

n_walls = length(sim.ice_floes)
info("added $(n_walls) fixed ice floes as walls")

# Initialize ice floes in wedge north of the constriction
iy = 1
dy = sqrt((2.*r)^2. - dx^2.)
spacing_to_boundaries = 4.*r
floe_padding = .5*r
noise_amplitude = floe_padding
Base.Random.srand(1)
for y in (L[2] - r - noise_amplitude):(-(2.*r + floe_padding)):((L[2] - 
    Ly_constriction)/2. + Ly_constriction)
    for x in (r + noise_amplitude):(2.*r + floe_padding):(Lx - r - 
                                                          noise_amplitude)
        if iy % 2 == 0
            x += 1.5*r
        end

        x_ = x + noise_amplitude*(0.5 - Base.Random.rand())
        y_ = y + noise_amplitude*(0.5 - Base.Random.rand())

        if y_ < -dy/dx*x_ + L[2] + spacing_to_boundaries
            continue
        end
            
        if y_ < dy/dx*x_ + (L[2] - dy/dx*Lx) + spacing_to_boundaries
            continue
        end
            
        SeaIce.addIceFloeCylindrical(sim, [x_, y_], r, h, verbose=false)
    end
    iy += 1
end
n = length(sim.ice_floes) - n_walls
info("added $(n) ice floes")

# Remove old simulation files
SeaIce.removeSimulationFiles(sim)

k_n = 1e6  # N/m
k_t = k_n
gamma_t = 1e7  # N/(m/s)
mu_d = 0.7
rotating = true
for i=1:length(sim.ice_floes)
    sim.ice_floes[i].contact_stiffness_normal = k_n
    sim.ice_floes[i].contact_stiffness_tangential = k_t
    sim.ice_floes[i].contact_viscosity_tangential = gamma_t
    sim.ice_floes[i].contact_dynamic_friction = mu_d
    sim.ice_floes[i].rotating = rotating
end

# Set temporal parameters
SeaIce.setTotalTime!(sim, 12.*60.*60.)
SeaIce.setOutputFileInterval!(sim, 60.)
SeaIce.setTimeStep!(sim)

# Run simulation for 10 time steps, then add new icefloes from the top
while sim.time < sim.time_total
    for it=1:10
        SeaIce.run!(sim, status_interval=1, single_step=true,
                    contact_tangential_rheology="Linear Viscous Frictional")
    end
    for i=1:size(sim.ocean.xh, 1)
        if sim.ocean.ice_floe_list[i, end] == []

            x, y = SeaIce.getCellCenterCoordinates(sim.ocean, i, 
                                                   size(sim.ocean.xh, 2))

            # Enable for high mass flux
            SeaIce.addIceFloeCylindrical(sim, [x-r, y-r], r, h, verbose=false,
                    contact_stiffness_normal=k_n,
                    contact_stiffness_normal=k_n,
                    contact_viscosity_tangential=gamma_t,
                    contact_dynamic_friction = mu_d,
                    rotating=rotating)
            SeaIce.addIceFloeCylindrical(sim, [x+r, y-r], r, h, verbose=false,
                    contact_stiffness_normal=k_n,
                    contact_stiffness_normal=k_n,
                    contact_viscosity_tangential=gamma_t,
                    contact_dynamic_friction = mu_d,
                    rotating=rotating)
            SeaIce.addIceFloeCylindrical(sim, [x+r, y+r], r, h, verbose=false,
                    contact_stiffness_normal=k_n,
                    contact_stiffness_normal=k_n,
                    contact_viscosity_tangential=gamma_t,
                    contact_dynamic_friction = mu_d,
                    rotating=rotating)
            SeaIce.addIceFloeCylindrical(sim, [x-r, y+r], r, h, verbose=false,
                    contact_stiffness_normal=k_n,
                    contact_stiffness_normal=k_n,
                    contact_viscosity_tangential=gamma_t,
                    contact_dynamic_friction = mu_d,
                    rotating=rotating)

            # Enable for low mass flux
            #x += noise_amplitude*(0.5 - Base.Random.rand())
            #y += noise_amplitude*(0.5 - Base.Random.rand())
            #SeaIce.addIceFloeCylindrical(sim, [x, y], r, h, verbose=false)
        end
    end
end
