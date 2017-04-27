#!/usr/bin/env julia
using SeaIce

sim = createSimulation(id="nares_strait")

# Initialize ocean
Lx = 50.e3
Lx_constriction = Lx*.33
L = [Lx, Lx, 1e3]
Ly_constriction = L[2]/2.
n = [100, 100, 2]
sim.ocean = createRegularOceanGrid(n, L, name="poiseuille_flow")
sim.ocean.v[:, :, 1, 1] = 1e-8*((sim.ocean.xq - Lx/2.).^2 - Lx^2./4.)

# Initialize confining walls, which are ice floes that are fixed in space
r = .5e3
h = 1.

## N-S segments
for y in linspace((Lx - Lx_constriction)/2.,
                  Ly_constriction + (L[2] - Ly_constriction)/2., 
                  Int(floor(Ly_constriction/(r*2))))
    addIceFloeCylindrical(sim, [(Lx - Lx_constriction)/2., y], r, h, fixed=true, 
verbose=false)
end
for y in linspace((L[2] - Ly_constriction)/2.,
                  Ly_constriction + (L[2] - Ly_constriction)/2., 
                  Int(floor(Ly_constriction/(r*2))))
    addIceFloeCylindrical(sim, [Ly_constriction + (L[2] - Ly_constriction)/2., 
y], r, h, fixed=true, verbose=false)
end

dx = r*cos(deg2rad(atan(Ly_constriction/Lx_constriction)))*2.
dy = r*cos(deg2rad(atan(Lx_constriction/Ly_constriction)))*2.

## NW diagonal
x = r:dx:((Lx - Lx_constriction)/2.)
#y = (L[2] - r):(-dy):((L[2] - Ly_constriction)/2. + Ly_constriction)
y = linspace(L[2] - r, (L[2] - Ly_constriction)/2. + Ly_constriction, length(x))
for i in 1:length(x)
    addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, verbose=false)
end

## NE diagonal
x = (L[1] - r):(-dx):((Lx - Lx_constriction)/2. + Lx_constriction)
#y = (L[2] - r):(-dy):((L[2] - Ly_constriction)/2. + Ly_constriction)
y = linspace(L[2] - r, (L[2] - Ly_constriction)/2. + Ly_constriction, length(x))
for i in 1:length(x)
    addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, verbose=false)
end

## SW diagonal
x = r:dx:((Lx - Lx_constriction)/2.)
#y = (r):(dy):((L[2] - Ly_constriction)/2.)
y = linspace(r, (L[2] - Ly_constriction)/2., length(x))
for i in 1:length(x)
    addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, verbose=false)
end

## SE diagonal
x = (L[1] - r):(-dx):((Lx - Lx_constriction)/2. + Lx_constriction)
#y = (r):(dy):((L[2] - Ly_constriction)/2.)
y = linspace(r, (L[2] - Ly_constriction)/2., length(x))
for i in 1:length(x)
    addIceFloeCylindrical(sim, [x[i], y[i]], r, h, fixed=true, verbose=false)
end



info("added $(length(sim.ice_floes)) fixed ice floes as walls")
writeVTK(sim)


# Initialize ice floes


# Run temporal loop
