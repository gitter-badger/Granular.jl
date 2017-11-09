## Functions for creating grain packings

export regularPacking!
"""

    regularPacking!(simulation, n, r_min, r_max[, padding_factor,
                    size_distribution, size_distribution_parameter])

Create a grid-based regular packing with grain numbers along each axis specified
by the `n` vector.

# Arguments
* `simulation::Simulation`: simulation object where the grains are inserted,
    preferably not containing prior grains.
* `n::Vector{Integer}`: 2-element vector determining number of grains along the
    `x` and `y` axes.
* `r_min::Real`: minimum desired grain radius.
* `r_max::Real`: maximum desired grain radius.
* `padding_factor::Real`: percentage-wise padding around each grain to allow for
    random perturbations to grain position.
* `size_distribution::String`: grain-size distribution to sample. Valid values
    are "powerlaw" and "uniform".
* `size_distribution_parameter::Real`: parameter to pass to the grain-size
    distribution generating function.
* `seed::Integer`: seed value to the pseudo-random number generator.
"""
function regularPacking!(simulation::Simulation,
                         n::Vector{Int},
                         r_min::Real,
                         r_max::Real;
                         padding_factor::Real=.1,
                         size_distribution::String="powerlaw",
                         size_distribution_parameter::Real=-1.8,
                         seed::Integer=1)

    r_rand = 0.
    pos = zeros(2)
    h = .5   # disc tickness
    dx = r_max * 2. * (1. + padding_factor)  # cell size
    dx_padding = r_max * 2. * padding_factor
    srand(seed)

    for iy in 1:n[2]
        for ix in 1:n[1]

            if size_distribution == "powerlaw"
                r_rand = Granular.randpower(1, size_distribution_parameter,
                                            r_min, r_max)
            elseif size_distribution == "uniform"
                r_rand = rand()*(r_max - r_min) + r_min
            end

            # Determine position from grid index and sample randomly from within
            # padding
            pos = [ix*dx - .5*dx, iy*dx - .5*dx] .+
                rand(2) .* dx_padding .- .5*dx_padding

            addGrainCylindrical!(simulation, pos, r_rand, h, verbose=false)
        end
    end

end

"""
Return random point in spherical annulus between `(r_i + r_j)` and `2.*(r_i +
r_j)` around `x_i`.  Note: there is slightly higher point density towards (r_i +
r_j).
"""
function generateNeighboringPoint(x_i::Vector, r_i::Real, r_j::Real,
                                  max_padding_factor::Real)

    R = rand() * (r_i + r_j) * max_padding_factor + 2. * (r_i + r_j)
    T = rand() * 2. * pi
    return [x_i[1] + R * sin(T), x_i[2] + R * cos(T)]
end

export poissonDiscSampling
"""
Generate disc packing in 2D using Poisson disc sampling with O(N) complexity, as
described by [Robert Bridson (2007)](http://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf).

# Arguments
* `simulation::Simulation`: simulation object where grains are inserted.
* `radius_max::Real`: largest grain radius to use.
* `radius_min::Real`: smallest grain radius to use.
* `sample_limit::Integer=30`: number of points to sample around each grain
    before giving up.
* `max_padding_factor::Real=2.`: this factor scales the padding to use during ice
    floe generation in numbers of grain diameters.
* `verbose::Bool=true`: show diagnostic information to stdout.

"""
function poissonDiscSampling(simulation::Simulation;
                             radius_max::Real=.1,
                             radius_min::Real=.1,
                             sample_limit::Integer=30,
                             max_padding_factor::Real=2.,
                             thickness::Real=1.,
                             verbose::Bool=true)

    active_list = Int[]  # list of points to originate search from

    # Step 0: Use existing `grid` (ocean or atmosphere) for contact-search grid
    if typeof(simulation.ocean.input_file) != Bool
        grid = simulation.ocean
    elseif typeof(simulation.atmosphere.input_file) != Bool
        grid = simulation.atmosphere
    else
        error("poissonDiscSampling requires an ocean or atmosphere grid")
    end
    # save grid boundaries
    sw, se, ne, nw = getGridCornerCoordinates(grid.xq, grid.yq)
    width_x = se[1] - sw[1]  # assume regular grid
    width_y = nw[2] - sw[2]  # assume regular grid

    # Step 1: If grid is empty: select random initial sample and save its index
    # to the background grid. Else: Make all grains active for search
    if isempty(simulation.grains)
        r = rand()*(radius_max - radius_min) + radius_min
        x0 = rand(2).*[width_x, width_y] + sw
        addGrainCylindrical!(simulation, x0, r, thickness)
        sortGrainsInGrid!(simulation, grid)
        push!(active_list, 1)
    else
        for i=1:length(simulation.grains)
            push!(active_list, i)
        end
    end

    # Step 2: While the active list is not empty, choose a random index `i` from
    # it.  Generate up to `k` points chosen uniformly from the spherical annulus
    # between radius `2*(r_i+r_j)` and `max_padding_factor*2*(r_i+r_j)` around
    # `x_i`.
    i = 0; j = 0; x_i = zeros(2); x_j = zeros(2); r_i = 0.; r_j = 0.
    while !isempty(active_list)

        i = rand(1:length(active_list))
        deleteat!(active_list, i)

        x_i = simulation.grains[i].lin_pos
        r_i = simulation.grains[i].contact_radius
        r_j = rand()*(radius_max - radius_min) + radius_min

        for j=1:sample_limit

            x_j = generateNeighboringPoint(x_i, r_i, r_j, max_padding_factor)
            if !(isPointInGrid(grid, x_j))
                continue
            end

            # if it doesn't collide, add it to the active list
            if !checkCollisions(grid, x_j, dx, position_list, radius_list)
                push!(position_list, x_j)
                push!(radius_list, r_j)
                push!(active_list, length(radius_list))
            end
        end
        println("Active points: $(length(active_list))")
    end
    if verbose
        info("Generated $(length(radius_list)) points")
    end
    return position_list, radius_list
end
