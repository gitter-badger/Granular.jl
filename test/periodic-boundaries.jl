#!/usr/bin/env julia

info("#### $(basename(@__FILE__)) ####")

verbose=false
i = 0

info("Testing assignment and reporting of grid boundary conditions")
ocean = Granular.createEmptyOcean()

Test.@test ocean.bc_west == 1
Test.@test ocean.bc_east == 1
Test.@test ocean.bc_north == 1
Test.@test ocean.bc_south == 1

const originalSTDOUT = STDOUT
(out_r, out_w) = redirect_stdout()
Granular.reportGridBoundaryConditions(ocean)
close(out_w)
redirect_stdout(originalSTDOUT)
output = convert(String, readavailable(out_r))
Test.@test output == """West  (-x): inactive\t(1)
East  (+x): inactive\t(1)
South (-y): inactive\t(1)
North (+y): inactive\t(1)
"""

(out_r, out_w) = redirect_stdout()
Granular.setGridBoundaryConditions!(ocean, "periodic", "south, west", verbose=true)
close(out_w)
redirect_stdout(originalSTDOUT)
output = convert(String, readavailable(out_r))
Test.@test output == """West  (-x): periodic\t(2)
East  (+x): inactive\t(1)
South (-y): periodic\t(2)
North (+y): inactive\t(1)
"""
Test.@test ocean.bc_west == 2
Test.@test ocean.bc_east == 1
Test.@test ocean.bc_north == 1
Test.@test ocean.bc_south == 2

(out_r, out_w) = redirect_stdout()
Granular.setGridBoundaryConditions!(ocean, "inactive", "all", verbose=false)
close(out_w)
redirect_stdout(originalSTDOUT)
output = convert(String, readavailable(out_r))
Test.@test output == ""
Test.@test ocean.bc_west == 1
Test.@test ocean.bc_east == 1
Test.@test ocean.bc_north == 1
Test.@test ocean.bc_south == 1

(out_r, out_w) = redirect_stdout()
Granular.setGridBoundaryConditions!(ocean, "periodic", "all")
close(out_w)
output = convert(String, readavailable(out_r))
redirect_stdout(originalSTDOUT)
Test.@test output == """West  (-x): periodic\t(2)
East  (+x): periodic\t(2)
South (-y): periodic\t(2)
North (+y): periodic\t(2)
"""
Test.@test ocean.bc_west == 2
Test.@test ocean.bc_east == 2
Test.@test ocean.bc_north == 2
Test.@test ocean.bc_south == 2

(out_r, out_w) = redirect_stdout()
Granular.setGridBoundaryConditions!(ocean, "inactive")
close(out_w)
output = convert(String, readavailable(out_r))
redirect_stdout(originalSTDOUT)
Test.@test output == """West  (-x): inactive\t(1)
East  (+x): inactive\t(1)
South (-y): inactive\t(1)
North (+y): inactive\t(1)
"""
Test.@test ocean.bc_west == 1
Test.@test ocean.bc_east == 1
Test.@test ocean.bc_north == 1
Test.@test ocean.bc_south == 1

(out_r, out_w) = redirect_stdout()
Granular.setGridBoundaryConditions!(ocean, "periodic")
close(out_w)
output = convert(String, readavailable(out_r))
redirect_stdout(originalSTDOUT)
Test.@test output == """West  (-x): periodic\t(2)
East  (+x): periodic\t(2)
South (-y): periodic\t(2)
North (+y): periodic\t(2)
"""
Test.@test ocean.bc_west == 2
Test.@test ocean.bc_east == 2
Test.@test ocean.bc_north == 2
Test.@test ocean.bc_south == 2

Test.@test_throws ErrorException Granular.setGridBoundaryConditions!(ocean,
                                                                     "periodic",
                                                                     "asdf")

Test.@test_throws ErrorException Granular.setGridBoundaryConditions!(ocean,
                                                                     "asdf")
