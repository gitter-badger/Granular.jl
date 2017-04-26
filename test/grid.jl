#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

ocean = SeaIce.readOceanNetCDF("Baltic/00010101.ocean_month.nc",
                               "Baltic/ocean_hgrid.nc")

info("Testing area-determination methods")
@test SeaIce.areaOfTriangle([0., 0.], [1., 0.], [0., 1.]) ≈ .5
@test SeaIce.areaOfTriangle([1., 0.], [0., 1.], [0., 0.]) ≈ .5
@test SeaIce.areaOfQuadrilateral([1., 0.], [0., 1.], [0., 0.], [1., 1.]) ≈ 1.

info("Testing area-based cell content determination")
@test SeaIce.isPointInCell(ocean, 2, 2, [6.5, 53.5]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.5]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.0, 53.5]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.7]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.9]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.99999]) == true
@test SeaIce.isPointInCell(ocean, 2, 2, [7.5, 53.5]) == false
@test SeaIce.isPointInCell(ocean, 2, 2, [0.0, 53.5]) == false

info("Testing conformal mapping methods")
@test SeaIce.conformalQuadrilateralCoordinates([0., 0.],
                                               [5., 0.],
                                               [5., 3.],
                                               [0., 3.],
                                               [2.5, 1.5]) ≈ [0.5, 0.5]
@test SeaIce.conformalQuadrilateralCoordinates([0., 0.],
                                               [5., 0.],
                                               [5., 3.],
                                               [0., 3.],
                                               [7.5, 1.5]) ≈ [1.5, 0.5]
@test SeaIce.conformalQuadrilateralCoordinates([0., 0.],
                                               [5., 0.],
                                               [5., 3.],
                                               [0., 3.],
                                               [7.5,-1.5]) ≈ [1.5,-0.5]
info("Checking cell content using conformal mapping methods")
@test SeaIce.isPointInCell(ocean, 2, 2, [6.4, 53.4], method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.5], method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.0, 53.5], method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.7], method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.9], method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [6.1, 53.99999],
                           method="Conformal") == true
@test SeaIce.isPointInCell(ocean, 2, 2, [7.5, 53.5],
                           method="Conformal") == false
@test SeaIce.isPointInCell(ocean, 2, 2, [0.0, 53.5],
                           method="Conformal") == false
