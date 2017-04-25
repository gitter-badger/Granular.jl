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

info("Testing ocean state interpolation")
@test_throws ErrorException SeaIce.findContacts!(ocean, time=0.)
@test_throws ErrorException SeaIce.findContacts!(ocean, time=1.e34)
u1, v1, h1, e1 = SeaIce.interpolateOceanState(ocean, ocean.time[1])
u2, v2, h2, e2 = SeaIce.interpolateOceanState(ocean, ocean.time[2])
u1_5, v1_5, h1_5, e1_5 = SeaIce.interpolateOceanState(ocean,
    ocean.time[1] + (ocean.time[2] - ocean.time[1])/2.)
@test_approx_eq u1 ocean.u[:, :, :, 1]
@test_approx_eq v1 ocean.v[:, :, :, 1]
@test_approx_eq h1 ocean.h[:, :, :, 1]
@test_approx_eq e1 ocean.e[:, :, :, 1]
@test_approx_eq u2 ocean.u[:, :, :, 2]
@test_approx_eq v2 ocean.v[:, :, :, 2]
@test_approx_eq h2 ocean.h[:, :, :, 2]
@test_approx_eq e2 ocean.e[:, :, :, 2]
@test_approx_eq u1_5 (ocean.u[:, :, :, 1] + ocean.u[:, :, :, 2])/2.
@test_approx_eq v1_5 (ocean.v[:, :, :, 1] + ocean.v[:, :, :, 2])/2.
@test_approx_eq h1_5 (ocean.h[:, :, :, 1] + ocean.h[:, :, :, 2])/2.
@test_approx_eq e1_5 (ocean.e[:, :, :, 1] + ocean.e[:, :, :, 2])/2.
