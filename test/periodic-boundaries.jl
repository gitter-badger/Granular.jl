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


info("Testing granular interaction across periodic boundaries")
sim = Granular.createSimulation()
sim.ocean = Granular.createRegularOceanGrid([5, 5, 2], [1., 1., 1.])
Granular.setGridBoundaryConditions!(sim.ocean, "periodic")
Granular.addGrainCylindrical!(sim, [0.1, 0.5], 0.11, 0.1, verbose=false)
Granular.addGrainCylindrical!(sim, [0.9, 0.5], 0.11, 0.1, verbose=false)

# there should be an error if all-to-all contact search is used
Test.@test_throws ErrorException Granular.findContacts!(sim)
Test.@test_throws ErrorException Granular.findContacts!(sim, method="all to all")
Test.@test_throws ErrorException Granular.findContactsAllToAll!(sim)

Granular.sortGrainsInGrid!(sim, sim.ocean, verbose=false)
Granular.findContacts!(sim, method="ocean grid")
Test.@test 2 == sim.grains[1].contacts[1]
Test.@test 1 == sim.grains[1].n_contacts
Test.@test 1 == sim.grains[2].n_contacts

