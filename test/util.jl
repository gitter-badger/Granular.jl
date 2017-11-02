#!/usr/bin/env julia

info("#### $(basename(@__FILE__)) ####")

info("Testing power-law RNG")

Test.@test 1 == length(Granular.randpower())
Test.@test () == size(Granular.randpower())
Test.@test 1 == length(Granular.randpower(1))
Test.@test () == size(Granular.randpower(1))
Test.@test 4 == length(Granular.randpower((2,2)))
Test.@test (2,2) == size(Granular.randpower((2,2)))
Test.@test 5 == length(Granular.randpower(5))
Test.@test (5,) == size(Granular.randpower(5))

srand(1)
for i=1:10^5
    Test.@test 0. <= Granular.randpower() <= 1.
    Test.@test 0. <= Granular.randpower(1, 1., 0., 1.) <= 1.
    Test.@test 0. <= Granular.randpower(1, 1., 0., .1) <= .1
    Test.@test 5. <= Granular.randpower(1, 1., 5., 6.) <= 6.
    Test.@test 0. <= minimum(Granular.randpower((2,2), 1., 0., 1.))
    Test.@test 1. >= maximum(Granular.randpower((2,2), 1., 0., 1.))
    Test.@test 0. <= minimum(Granular.randpower(5, 1., 0., 1.))
    Test.@test 1. >= minimum(Granular.randpower(5, 1., 0., 1.))
end
