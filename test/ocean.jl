#!/usr/bin/env julia

# Check if ocean-specific functions and grid operations are functioning 
# correctly

info("#### $(basename(@__FILE__)) ####")

info("Testing regular grid generation")
ocean = SeaIce.createRegularOceanGrid([6, 6, 6], [1., 1., 1.])
@test size(ocean.xq) == (7, 7)
@test size(ocean.yq) == (7, 7)
@test size(ocean.xh) == (6, 6)
@test size(ocean.yh) == (6, 6)
@test ocean.xq[1, :, 1] ≈ zeros(7)
@test ocean.xq[4, :, 1] ≈ .5*ones(7)
@test ocean.xq[end, :, 1] ≈ 1.*ones(7)
@test ocean.yq[:, 1, 1] ≈ zeros(7)
@test ocean.yq[:, 4, 1] ≈ .5*ones(7)
@test ocean.yq[:, end, 1] ≈ 1.*ones(7)
@test size(ocean.u) == (7, 7, 6, 1)
@test size(ocean.v) == (7, 7, 6, 1)
@test size(ocean.h) == (7, 7, 6, 1)
@test size(ocean.e) == (7, 7, 6, 1)
@test ocean.u ≈ zeros(7, 7, 6, 1)
@test ocean.v ≈ zeros(7, 7, 6, 1)
@test ocean.h ≈ zeros(7, 7, 6, 1)
@test ocean.e ≈ zeros(7, 7, 6, 1)
