#!/usr/bin/env julia

# Check if NetCDF files are read correctly from the disk.

info("#### $(basename(@__FILE__)) ####")

Test.@test_throws ErrorException Granular.readOceanStateNetCDF("nonexistentfile")
Test.@test_throws ErrorException Granular.readOceanGridNetCDF("nonexistentfile")

info("Testing dimensions of content read from Baltic test case")
ocean = Granular.readOceanNetCDF("Baltic/00010101.ocean_month.nc",
                               "Baltic/ocean_hgrid.nc")
Test.@test ocean.time/(24.*60.*60.) ≈ [.5, 1.5, 2.5, 3.5, 4.5]
Test.@test size(ocean.xq) == (24, 15)
Test.@test size(ocean.yq) == (24, 15)
Test.@test size(ocean.xh) == (23, 14)
Test.@test size(ocean.yh) == (23, 14)
Test.@test size(ocean.u) == (24, 15, 63, 5)
Test.@test size(ocean.v) == (24, 15, 63, 5)
Test.@test size(ocean.h) == (23, 14, 63, 5)
Test.@test size(ocean.e) == (23, 14, 64, 5)

info("Testing ocean state interpolation")
Test.@test_throws ErrorException Granular.interpolateOceanState(ocean, time=0.)
Test.@test_throws ErrorException Granular.interpolateOceanState(ocean, time=1.e34)
u1, v1, h1, e1 = Granular.interpolateOceanState(ocean, ocean.time[1])
u2, v2, h2, e2 = Granular.interpolateOceanState(ocean, ocean.time[2])
Test.@test_throws ErrorException Granular.interpolateOceanState(ocean, -1.)
u1_5, v1_5, h1_5, e1_5 = Granular.interpolateOceanState(ocean,
    ocean.time[1] + (ocean.time[2] - ocean.time[1])/2.)
Test.@test u1 ≈ ocean.u[:, :, :, 1]
Test.@test v1 ≈ ocean.v[:, :, :, 1]
Test.@test h1 ≈ ocean.h[:, :, :, 1]
Test.@test e1 ≈ ocean.e[:, :, :, 1]
Test.@test u2 ≈ ocean.u[:, :, :, 2]
Test.@test v2 ≈ ocean.v[:, :, :, 2]
Test.@test h2 ≈ ocean.h[:, :, :, 2]
Test.@test e2 ≈ ocean.e[:, :, :, 2]
Test.@test u1_5 ≈ (ocean.u[:, :, :, 1] + ocean.u[:, :, :, 2])/2.
Test.@test v1_5 ≈ (ocean.v[:, :, :, 1] + ocean.v[:, :, :, 2])/2.
Test.@test h1_5 ≈ (ocean.h[:, :, :, 1] + ocean.h[:, :, :, 2])/2.
Test.@test e1_5 ≈ (ocean.e[:, :, :, 1] + ocean.e[:, :, :, 2])/2.
