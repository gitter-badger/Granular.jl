
"""
Use bilinear interpolation to interpolate a staggered grid to an arbitrary 
position in a cell.  Assumes north-east convention, i.e. (i,j) is located at the 
north-east corner.

# Arguments
* `field::Array{Float64, 4}`: a scalar field to interpolate from
* `xi::float`: relative x position in cell [-], must be in `[0., 1.]`
* `yj::float`: relative y position in cell [-], must be in `[0., 1.]`
* `i::Int`: i-index of cell containing point
* `j::Int`: j-index of cell containing point
* `grid_type::String="Arakawa A"`: grid system for `field`
"""
function bilinearInterpolation(field::Array{Float64, 4},
                               xi::float,
                               yj::float,
                               i::Int,
                               j::Int;
                               grid_type::String="Arakawa A")

    if xi < 0. || xi > 1. || yj < 0. || yj > 1.
        error("relative coordinates outside bounds ($(xi), $(yj))")
    end

    if grid_type == "Arakawa A"
        return (field[i,j]*xi + field[i-1,j]*(1. - xi))*yi +
            (field[i,j-1]*xi + field[i-1,j-1]*(1. - xi))*(1. - yi)
    else
        error("grid type not understood.")
    end
end

"""
Find ice-floe positions in ocean grid, based on their center positions.
"""
function sortIceFloesInOceanGrid!(simulation::Simulation, verbose=true)

    # TODO: initialize empty ice_floe_list before appending to list

    for idx in 1:length(simulation.ice_floes)

        if cellContainsIceFloe(simulation.ocean, i, j,
                               simulation.ice_floes[idx])

            # add cell to ice floe
            simulation.ice_floes[idx].ocean_grid_pos = [i, j]

            # add ice floe to cell
            push!(simulation.ice_floe_list[i, j], idx)
        end
    end
end

"""
Returns the `i`, `j` index of the ocean grid cell containing the `point`.
"""
function findCellContainingPoint(ocean::Ocean, i::Int, j::Int, 
                                 point::Array{float, 2})



    return i, j
end


"""
Check if a 2d point is contained inside a cell from the ocean grid.
The function uses either an area-based approach (`method = "Area"`), or a 
conformal mapping approach (`method = "Conformal"`).  The area-based approach is 
more robust.  This function returns `true` or `false`.
"""
function isPointInCell(ocean::Ocean, i::Int, j::Int, point::Array{float, 1};
                      method::String="Area")

    sw, se, ne, nw = getCellCornerCoordinates(ocean, i, j)

    if method == "Area"
        if areaOfQuadrilateral(sw, se, ne, nw) ≈
            areaOfTriangle(point, sw, se) +
            areaOfTriangle(point, se, ne) +
            areaOfTriangle(point, ne, nw) +
            areaOfTriangle(point, nw, sw)
            return true
        else
            return false
        end

    elseif method == "Conformal"
        x_tilde, y_tilde = conformalQuadrilateralCoordinates(sw, se, ne, nw,
                                                             point)
        if x_tilde >= 0. && x_tilde <= 1. && y_tilde >= 0. && y_tilde <= 1.
            return true
        else
            return false
        end
    else
        error("method not understood")
    end
end

"""
Returns ocean-grid corner coordinates in the following order (south-west corner, 
south-east corner, north-east corner, north-west corner).
"""
function getCellCornerCoordinates(ocean::Ocean, i::Int, j::Int)
    sw = [ocean.xq[i-1, j-1], ocean.yq[i-1, j-1]]
    se = [ocean.xq[  i, j-1], ocean.yq[  i, j-1]]
    ne = [ocean.xq[  i,   j], ocean.yq[  i,   j]]
    nw = [ocean.xq[i-1,   j], ocean.yq[i-1,   j]]
    return sw, se, ne, nw
end

"Returns the area of an triangle with corner coordinates `a`, `b`, and `c`."
function areaOfTriangle(a::Array{float, 1},
                        b::Array{float, 1},
                        c::Array{float, 1})
    return abs(
               (a[1]*(b[2] - c[2]) +
                b[1]*(c[2] - a[2]) +
                c[1]*(a[2] - b[2]))/2.
              )
end

"""
Returns the area of a quadrilateral with corner coordinates `a`, `b`, `c`, and 
`d`.  Corners `a` and `c` should be opposite of each other, the same must be 
true for `b` and `d`.  This is true if the four corners are passed as arguments 
in a "clockwise" or "counter-clockwise" manner.
"""
function areaOfQuadrilateral(a::Array{float, 1},
                             b::Array{float, 1},
                             c::Array{float, 1},
                             d::Array{float, 1})
    return areaOfTriangle(a, b, c) + areaOfTriangle(c, d, a)
end

"""
Returns the non-dimensional coordinates `[x_tilde, y_tilde]` of a point `p` 
within a quadrilateral with corner coordinates `A`, `B`, `C`, and `D`.
Points must be ordered in counter-clockwise order, starting from south-west 
corner.
"""
function conformalQuadrilateralCoordinates(A::Array{float, 1},
                                           B::Array{float, 1},
                                           C::Array{float, 1},
                                           D::Array{float, 1},
                                           p::Array{float, 1})

    if !(A[1] < B[1] && B[2] < C[2] && C[1] > D[1])
        error("corner coordinates are not passed in the correct order")
    end
    alpha = B[1] - A[1]
    delta = B[2] - A[2]
    beta = D[1] - A[1]
    epsilon = D[2] - A[2]
    gamma = C[1] - A[1] - (alpha + beta)
    kappa = C[2] - A[2] - (delta + epsilon)
    a = kappa*beta - gamma*epsilon
    dx = p[1] - A[1]
    dy = p[2] - A[2]
    b = (delta*beta - alpha*epsilon) - (kappa*dx - gamma*dy)
    c = alpha*dy - delta*dx
    if abs(a) > 0.
        d = b^2./4. - a*c
        if d >= 0.
            yy1 = -(b/2. + sqrt(d))/a
            yy2 = -(b/2. - sqrt(d))/a
            if abs(yy1 - .5) < abs(yy2 - .5)
                y_tilde = yy1
            else
                y_tilde = yy2
            end
        else
            error("could not perform conformal mapping\n",
                  "A = $(A), B = $(B), C = $(C), D = $(D), point = $(p),\n",
                  "alpha = $(alpha), beta = $(beta), gamma = $(gamma), ",
                  "delta = $(delta), epsilon = $(epsilon), kappa = $(kappa)")
        end
    else
        if !(b ≈ 0.)
            y_tilde = -c/b
        else
            y_tilde = 0.
        end
    end
    a = alpha + gamma*y_tilde
    b = delta + kappa*y_tilde
    if !(a ≈ 0.)
        x_tilde = (dx - beta*y_tilde)/a
    elseif !(b ≈ 0.)
        x_tilde = (dy - epsilon*y_tilde)/b
    else
        error("could not determine non-dimensional position in quadrilateral ",
              "(a = 0. and b = 0.)\n",
              "A = $(A), B = $(B), C = $(C), D = $(D), point = $(p),\n",
              "alpha = $(alpha), beta = $(beta), gamma = $(gamma), ",
              "delta = $(delta), epsilon = $(epsilon), kappa = $(kappa)")
    end
    return [x_tilde, y_tilde]
end
