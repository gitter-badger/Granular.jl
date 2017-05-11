"""
    bilinearInterpolation(field, x_tilde, y_tilde, i, j, k, it)

Use bilinear interpolation to interpolate a staggered grid to an arbitrary 
position in a cell.  Assumes south-west convention, i.e. (i,j) is located at the 
south-west (-x, -y)-facing corner.

# Arguments
* `field::Array{Float64, 4}`: a scalar field to interpolate from
* `x_tilde::float`: x point position [0;1]
* `y_tilde::float`: y point position [0;1]
* `i::Int`: i-index of cell containing point
* `j::Int`: j-index of scalar field to interpolate from
* `it::Int`: time step from scalar field to interpolate from
"""
function bilinearInterpolation(field::Array{Float64, 4},
                               x_tilde::Float64,
                               y_tilde::Float64,
                               i::Int,
                               j::Int,
                               k::Int,
                               it::Int)

    if x_tilde < 0. || x_tilde > 1. || y_tilde < 0. || y_tilde > 1.
        error("relative coordinates outside bounds ($(x_tilde), $(y_tilde))")
    end

    return (field[i+1, j+1, k, it]*x_tilde +
            field[i, j+1, k, it]*(1. - x_tilde))*y_tilde +
           (field[i+1, j, k, it]*x_tilde +
            field[i, j, k, it]*(1. - x_tilde))*(1. - y_tilde)
end

"""
    curl(ocean, x_tilde, y_tilde, i, j, k, it)

Use bilinear interpolation to interpolate curl value for a staggered velocity 
grid to an arbitrary position in a cell.  Assumes south-west convention, i.e.  
(i,j) is located at the south-west (-x, -y)-facing corner.

# Arguments
* `ocean::Ocean`: grid for which to determine curl
* `x_tilde::float`: x point position [0;1]
* `y_tilde::float`: y point position [0;1]
* `i::Int`: i-index of cell containing point
* `j::Int`: j-index of scalar field to interpolate from
* `it::Int`: time step from scalar field to interpolate from
"""
function curl(ocean::Ocean,
              x_tilde::Float64,
              y_tilde::Float64,
              i::Int,
              j::Int,
              k::Int,
              it::Int)

    sw, se, ne, nw = getCellCornerCoordinates(ocean, i, j)
    sw_se = norm(sw - se)
    se_ne = norm(se - ne)
    nw_ne = norm(nw - ne)
    sw_nw = norm(sw - nw)

    return (
    ((ocean.v[i+1, j  , k,it] - ocean.v[i  , j  , k,it])/sw_se*(1. - y_tilde) +
     ((ocean.v[i+1, j+1, k,it] - ocean.v[i  , j+1, k,it])/nw_ne)*y_tilde) -
    ((ocean.u[i  , j+1, k,it] - ocean.u[i  , j  , k,it])/sw_nw*(1. - x_tilde) +
     ((ocean.u[i+1, j+1, k,it] - ocean.u[i+1, j  , k,it])/se_ne)*x_tilde))
end

export sortIceFloesInOceanGrid!
"""
Find ice-floe positions in ocean grid, based on their center positions.
"""
function sortIceFloesInOceanGrid!(simulation::Simulation; verbose=false)

    simulation.ocean.ice_floe_list =
        Array{Array{Int, 1}}(size(simulation.ocean.xh, 1), 
                             size(simulation.ocean.xh, 2))
    for i=1:size(simulation.ocean.xh, 1)
        for j=1:size(simulation.ocean.xh, 2)
            simulation.ocean.ice_floe_list[i, j] = Int[]
        end
    end

    for idx in 1:length(simulation.ice_floes)

        i, j = findCellContainingPoint(simulation.ocean,
                                       simulation.ice_floes[idx].lin_pos)

        # remove ice floe if it is outside of the grid
        if i == 0 && j == 0
            disableIceFloe!(simulation, idx)
            continue
        end

        # add cell to ice floe
        simulation.ice_floes[idx].ocean_grid_pos = [i, j]

        # add ice floe to cell
        push!(simulation.ocean.ice_floe_list[i, j], idx)
    end
end

export findCellContainingPoint
"""
    findCellContainingPoint(ocean, point[, method])

Returns the `i`, `j` index of the ocean grid cell containing the `point`.
The function uses either an area-based approach (`method = "Area"`), or a 
conformal mapping approach (`method = "Conformal"`).  The area-based approach is 
more robust.  This function returns the coordinates of the cell.  If no match is 
found the function returns `(0,0)`.

# Arguments
* `ocean::Ocean`: object containing ocean data.
* `point::Array{float, 1}`: two-dimensional vector of point to check.
* `method::String`: approach to use for determining if point is inside cell or 
    not, can be "Conformal" (default) or "Areal".
"""
function findCellContainingPoint(ocean::Ocean, point::Array{float, 1};
                                 method::String="Conformal")
    for i=1:size(ocean.xh, 1)
        for j=1:size(ocean.yh, 2)
            if isPointInCell(ocean, i, j, point, method=method)
                return i, j
            end
        end
    end
    return 0, 0
end

export getNonDimensionalCellCoordinates
"""
Returns the non-dimensional conformal mapped coordinates for point `point` in 
cell `i,j`, based off the coordinates in the `ocean` grid.

This function is a wrapper for `getCellCornerCoordinates()` and 
`conformalQuadrilateralCoordinates()`.
"""
function getNonDimensionalCellCoordinates(ocean::Ocean, i::Int, j::Int,
                                          point::Array{float, 1})

    sw, se, ne, nw = getCellCornerCoordinates(ocean, i, j)
    x_tilde, y_tilde = conformalQuadrilateralCoordinates(sw, se, ne, nw, point)
    return [x_tilde, y_tilde]
end

export isPointInCell
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

export getCellCornerCoordinates
"""
    getCellCornerCoordinates(ocean, i, j)

Returns ocean-grid corner coordinates in the following order (south-west corner, 
south-east corner, north-east corner, north-west corner).

# Arguments
* `ocean::Ocean`: ocean object containing grid.
* `i::Int`: x-index of cell.
* `j::Int`: y-index of cell.
"""
function getCellCornerCoordinates(ocean::Ocean, i::Int, j::Int)
    sw = [ocean.xq[  i,   j], ocean.yq[  i,   j]]
    se = [ocean.xq[i+1,   j], ocean.yq[i+1,   j]]
    ne = [ocean.xq[i+1, j+1], ocean.yq[i+1, j+1]]
    nw = [ocean.xq[  i, j+1], ocean.yq[  i, j+1]]
    return sw, se, ne, nw
end

export getCellCenterCoordinates
"""
    getCellCenterCoordinates(ocean, i, j)

Returns ocean-grid center coordinates (h-point).

# Arguments
* `ocean::Ocean`: ocean object containing grid.
* `i::Int`: x-index of cell.
* `j::Int`: y-index of cell.
"""
function getCellCenterCoordinates(ocean::Ocean, i::Int, j::Int)
    return [ocean.xh[i, j], ocean.yh[i, j]]
end

export areaOfTriangle
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

export areaOfQuadrilateral
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

export conformalQuadrilateralCoordinates
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
