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
function bilinearInterpolation(field::Array{Float64, 2},
                               x_tilde::Float64,
                               y_tilde::Float64,
                               i::Int,
                               j::Int)

    if x_tilde < 0. || x_tilde > 1. || y_tilde < 0. || y_tilde > 1.
        error("relative coordinates outside bounds ($(x_tilde), $(y_tilde))")
    end

    return (field[i+1, j+1]*x_tilde + field[i, j+1]*(1. - x_tilde))*y_tilde +
           (field[i+1, j]*x_tilde + field[i, j]*(1. - x_tilde))*(1. - y_tilde)
end

"""
    curl(grid, x_tilde, y_tilde, i, j, k, it)

Use bilinear interpolation to interpolate curl value for a staggered velocity 
grid to an arbitrary position in a cell.  Assumes south-west convention, i.e.  
(i,j) is located at the south-west (-x, -y)-facing corner.

# Arguments
* `grid::Any`: grid for which to determine curl
* `x_tilde::float`: x point position [0;1]
* `y_tilde::float`: y point position [0;1]
* `i::Int`: i-index of cell containing point
* `j::Int`: j-index of scalar field to interpolate from
* `it::Int`: time step from scalar field to interpolate from
"""
function curl(grid::Any,
              x_tilde::Float64,
              y_tilde::Float64,
              i::Int,
              j::Int,
              k::Int,
              it::Int)

    sw, se, ne, nw = getCellCornerCoordinates(grid, i, j)
    sw_se = norm(sw - se)
    se_ne = norm(se - ne)
    nw_ne = norm(nw - ne)
    sw_nw = norm(sw - nw)

    return (
    ((grid.v[i+1, j  , k,it] - grid.v[i  , j  , k,it])/sw_se*(1. - y_tilde) +
     ((grid.v[i+1, j+1, k,it] - grid.v[i  , j+1, k,it])/nw_ne)*y_tilde) -
    ((grid.u[i  , j+1, k,it] - grid.u[i  , j  , k,it])/sw_nw*(1. - x_tilde) +
     ((grid.u[i+1, j+1, k,it] - grid.u[i+1, j  , k,it])/se_ne)*x_tilde))
end

export sortIceFloesInGrid!
"""
Find ice-floe positions in grid, based on their center positions.
"""
function sortIceFloesInGrid!(simulation::Simulation, grid::Any; verbose=false)

    grid.ice_floe_list =
        Array{Array{Int, 1}}(size(grid.xh, 1), size(grid.xh, 2))
    for i=1:size(grid.xh, 1)
        for j=1:size(grid.xh, 2)
            grid.ice_floe_list[i, j] = Int[]
        end
    end

    for idx in 1:length(simulation.ice_floes)

        if !simulation.ice_floes[idx].enabled
            continue
        end

        # After first iteration, check if ice floe is in same cell before 
        # traversing entire grid
        if typeof(grid) == Ocean
            i_old, j_old = simulation.ice_floes[idx].ocean_grid_pos
        elseif typeof(grid) == Atmosphere
            i_old, j_old = simulation.ice_floes[idx].atmosphere_grid_pos
        else
            error("grid type not understood.")
        end
        if simulation.time > 0. &&
            i_old > 0 && j_old > 0 &&
            isPointInCell(grid, i_old, j_old,
                         simulation.ice_floes[idx].lin_pos)
            i = i_old
            j = j_old

        else

            # Search for point in 8 neighboring cells
            nx = size(grid.xh, 1)
            ny = size(grid.xh, 2)
            found = false
            for i_rel=-1:1
                for j_rel=-1:1
                    if i_rel == 0 && j_rel == 0
                        continue  # cell previously searched
                    end
                    i_t = max(min(i_old + i_rel, nx), 1)
                    j_t = max(min(j_old + j_rel, ny), 1)
                    
                    if isPointInCell(grid, i_t, j_t,
                                     simulation.ice_floes[idx].lin_pos)
                        i = i_t
                        j = j_t
                        found = true
                        break
                    end
                end
                if found
                    break
                end
            end

            if !found
                i, j = findCellContainingPoint(grid,
                                               simulation.ice_floes[idx].lin_pos)
            end

            # remove ice floe if it is outside of the grid
            if i == 0 && j == 0
                disableIceFloe!(simulation, idx)
                continue
            end

            # add cell to ice floe
            if typeof(grid) == Ocean
                simulation.ice_floes[idx].ocean_grid_pos = [i, j]
            elseif typeof(grid) == Atmosphere
                simulation.ice_floes[idx].atmosphere_grid_pos = [i, j]
            else
                error("grid type not understood.")
            end
        end

        # add ice floe to cell
        push!(grid.ice_floe_list[i, j], idx)
    end
end

export findCellContainingPoint
"""
    findCellContainingPoint(grid, point[, method])

Returns the `i`, `j` index of the grid cell containing the `point`.
The function uses either an area-based approach (`method = "Area"`), or a 
conformal mapping approach (`method = "Conformal"`).  The area-based approach is 
more robust.  This function returns the coordinates of the cell.  If no match is 
found the function returns `(0,0)`.

# Arguments
* `grid::Any`: grid object containing ocean or atmosphere data.
* `point::Array{float, 1}`: two-dimensional vector of point to check.
* `method::String`: approach to use for determining if point is inside cell or 
    not, can be "Conformal" (default) or "Areal".
"""
function findCellContainingPoint(grid::Any, point::Array{float, 1};
                                 method::String="Conformal")

    for i=1:size(grid.xh, 1)
        for j=1:size(grid.yh, 2)
            if isPointInCell(grid, i, j, point, method=method)
                return i, j
            end
        end
    end
    return 0, 0
end

export getNonDimensionalCellCoordinates
"""
Returns the non-dimensional conformal mapped coordinates for point `point` in 
cell `i,j`, based off the coordinates in the grid.

This function is a wrapper for `getCellCornerCoordinates()` and 
`conformalQuadrilateralCoordinates()`.
"""
function getNonDimensionalCellCoordinates(grid::Any, i::Int, j::Int,
                                          point::Array{float, 1})

    sw, se, ne, nw = getCellCornerCoordinates(grid, i, j)
    x_tilde, y_tilde = conformalQuadrilateralCoordinates(sw, se, ne, nw, point)
    return [x_tilde, y_tilde]
end

export isPointInCell
"""
Check if a 2d point is contained inside a cell from the supplied grid.
The function uses either an area-based approach (`method = "Area"`), or a 
conformal mapping approach (`method = "Conformal"`).  The area-based approach is 
more robust.  This function returns `true` or `false`.
"""
function isPointInCell(grid::Any, i::Int, j::Int, point::Array{float, 1};
                       method::String="Conformal")

    sw, se, ne, nw = getCellCornerCoordinates(grid, i, j)

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
    getCellCornerCoordinates(grid, i, j)

Returns grid corner coordinates in the following order (south-west corner, 
south-east corner, north-east corner, north-west corner).

# Arguments
* `grid::Any`: grid object (Ocean or Atmosphere) containing grid.
* `i::Int`: x-index of cell.
* `j::Int`: y-index of cell.
"""
function getCellCornerCoordinates(grid::Any, i::Int, j::Int)
    sw = [grid.xq[  i,   j], grid.yq[  i,   j]]
    se = [grid.xq[i+1,   j], grid.yq[i+1,   j]]
    ne = [grid.xq[i+1, j+1], grid.yq[i+1, j+1]]
    nw = [grid.xq[  i, j+1], grid.yq[  i, j+1]]
    return sw, se, ne, nw
end

export getCellCenterCoordinates
"""
    getCellCenterCoordinates(grid, i, j)

Returns grid center coordinates (h-point).

# Arguments
* `grid::Any`: grid object containing grid.
* `i::Int`: x-index of cell.
* `j::Int`: y-index of cell.
"""
function getCellCenterCoordinates(grid::Any, i::Int, j::Int)
    return [grid.xh[i, j], grid.yh[i, j]]
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

export findEmptyPositionInGridCell
"""
Attempt locate an empty spot for an ice floe with radius `r` with center 
coordinates in a specified grid cell (`i`, `j`) without overlapping any other 
ice floes in that cell or the neighboring cells.  This function will stop 
attempting after `n_iter` iterations, each with randomly generated positions.

This function assumes that existing ice floes have been binned according to the 
grid (e.g., using `sortIceFloesInGrid()`).
"""
function findEmptyPositionInGridCell(simulation::Simulation,
                                     grid::Any,
                                     i::Int,
                                     j::Int,
                                     r::float;
                                     n_iter::Int = 10,
                                     seed::Int = 1,
                                     verbose::Bool = false)
    overlap_found = false
    i_iter = 0
    pos = [NaN NaN]

    nx, ny = size(grid.xh)

    for i_iter=1:n_iter

        overlap_found = false
        Base.Random.srand(i*j*seed*i_iter)
        # generate random candidate position
        x_tilde = Base.Random.rand()
        y_tilde = Base.Random.rand()
        pos = [bilinearInterpolation(grid.xq, x_tilde, y_tilde, i, j)
               bilinearInterpolation(grid.yq, x_tilde, y_tilde, i, j)]
        if verbose
            info("trying poisition $pos in cell $i,$j")
        end

        # search for contacts in current and eight neighboring cells
        for i_neighbor_corr=[0 -1 1]
            for j_neighbor_corr=[0 -1 1]

                # target cell index
                it = i + i_neighbor_corr
                jt = j + j_neighbor_corr

                # do not search outside grid boundaries
                if it < 1 || it > nx || jt < 1 || jt > ny
                    continue
                end

                # traverse list of ice floes in the target cell and check 
                # for overlaps
                for icefloe_idx in grid.ice_floe_list[it, jt]
                    overlap = norm(simulation.ice_floes[icefloe_idx].lin_pos - 
                                   pos) -
                        (simulation.ice_floes[icefloe_idx].contact_radius + r)

                    if overlap < 0.
                        if verbose
                            info("overlap with $icefloe_idx in cell $i,$j")
                        end
                        overlap_found = true
                        break
                    end
                end
            end
            if overlap_found == true
                break
            end
        end
        if overlap_found == false
            break
        end
    end
    if verbose && overlap_found == false
        info("Found position $pos in cell $i,$j after $i_iter iterations")
    elseif verbose && overlap_found
        info("Free position not found in cell $i,$j")
    end

    if overlap_found == false
        if isnan(pos[1]) || isnan(pos[2])
            error("fatal error: could not determine free position in cell")
        end
        return pos
    else
        if verbose
            warn("could not insert an ice floe into " *
                 "$(typeof(grid)) grid cell ($i, $j)")
        end
        return false
    end
end
