#!/usr/bin/env julia
import SeaIce
using Base.Test

info("#### $(basename(@__FILE__)) ####")

info("Testing power-law RNG")

@test 1 == length(SeaIce.randpower())
@test (1,) == size(SeaIce.randpower())
@test 1 == length(SeaIce.randpower(1))
@test (1,) == size(SeaIce.randpower(1))
@test 4 == length(SeaIce.randpower((2,2)))
@test (2,2) == size(SeaIce.randpower((2,2)))
@test 5 == length(SeaIce.randpower(5))
@test (5,) == size(SeaIce.randpower(5))

Base.Random.srand(1)
for i=1:10^5
    @test 0. <= SeaIce.randpower()[1] <= 1.
    @test 0. <= SeaIce.randpower(1, 1., 0., 1.)[1] <= 1.
    @test 0. <= SeaIce.randpower(1, 1., 0., .1)[1] <= .1
    @test 0. <= minimum(SeaIce.randpower((2,2), 1., 0., 1.))
    @test 1. >= maximum(SeaIce.randpower((2,2), 1., 0., 1.))
    @test 0. <= minimum(SeaIce.randpower(5, 1., 0., 1.))
    @test 1. >= minimum(SeaIce.randpower(5, 1., 0., 1.))
end
