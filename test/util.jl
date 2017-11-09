#!/usr/bin/env julia

info("#### $(basename(@__FILE__)) ####")

info("Testing power-law RNG")

@test 1 == length(Granular.randpower())
@test () == size(Granular.randpower())
@test 1 == length(Granular.randpower(1))
@test () == size(Granular.randpower(1))
@test 4 == length(Granular.randpower((2,2)))
@test (2,2) == size(Granular.randpower((2,2)))
@test 5 == length(Granular.randpower(5))
@test (5,) == size(Granular.randpower(5))

srand(1)
for i=1:10^5
    @test 0. <= Granular.randpower() <= 1.
    @test 0. <= Granular.randpower(1, 1., 0., 1.) <= 1.
    @test 0. <= Granular.randpower(1, 1., 0., .1) <= .1
    @test 5. <= Granular.randpower(1, 1., 5., 6.) <= 6.
    @test 0. <= minimum(Granular.randpower((2,2), 1., 0., 1.))
    @test 1. >= maximum(Granular.randpower((2,2), 1., 0., 1.))
    @test 0. <= minimum(Granular.randpower(5, 1., 0., 1.))
    @test 1. >= minimum(Granular.randpower(5, 1., 0., 1.))
end
