#!/usr/bin/env julia

# Check if NetCDF files are read correctly from the disk.

import Base.Test
import SeaIce

info("#### $(basename(@__FILE__)) ####")

info("Testing dimensions of content read from prog__0001_006.nc")
ocean = SeaIce.readOceanNetCDF("prog__0001_006.nc")
@test length(ocean.xq) == 44
@test length(ocean.xh) == 44
@test length(ocean.yq) == 40
@test length(ocean.yh) == 40
@test ocean.time ≈ [5., 10.]
@test size(ocean.u) == (44, 40, 2, 2)
@test size(ocean.v) == (44, 40, 2, 2)
@test size(ocean.h) == (44, 40, 2, 2)
@test size(ocean.e) == (44, 40, 3, 2)
