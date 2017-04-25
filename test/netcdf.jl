#!/usr/bin/env julia

# Check if NetCDF files are read correctly from the disk.

info("#### $(basename(@__FILE__)) ####")

info("Testing dimensions of content read from Baltic test case")
ocean = SeaIce.readOceanNetCDF("Baltic/00010101.ocean_month.nc",
                               "Baltic/ocean_hgrid.nc")
@test ocean.time/(24.*60.*60.) ≈ [.5, 1.5, 2.5, 3.5, 4.5]
@test size(ocean.xq) == (24, 15)
@test size(ocean.yq) == (24, 15)
@test size(ocean.xh) == (23, 14)
@test size(ocean.yh) == (23, 14)
@test size(ocean.u) == (24, 15, 63, 5)
@test size(ocean.v) == (24, 15, 63, 5)
@test size(ocean.h) == (23, 14, 63, 5)
@test size(ocean.e) == (23, 14, 64, 5)
