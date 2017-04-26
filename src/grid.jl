
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
Check if a 2d point is contained inside a cell from the ocean grid.  Returns 
`true`/`false`.
"""
function isPointInCell(ocean::Ocean, i::Int, j::Int, point::Array{float, 1})

    sw, nw, se, ne = getCellCornerCoordinates(ocean, i, j)

    if areaOfQuadrilateral(sw, nw, se, ne) ≈
        areaOfTriangle(point, sw, se) +
        areaOfTriangle(point, se, ne) +
        areaOfTriangle(point, ne, nw) +
        areaOfTriangle(point, nw, sw)
        return true
    else
        return false
    end
end

"""
Returns ocean-grid corner coordinates in the following order (south-west corner, 
north-west corner, south-east corner, north-east corner).
"""
function getCellCornerCoordinates(ocean::Ocean, i::Int, j::Int)
    sw = [ocean.xq[i-1, j-1], ocean.yq[i-1, j-1]]
    nw = [ocean.xq[i-1,   j], ocean.yq[i-1,   j]]
    se = [ocean.xq[  i, j-1], ocean.yq[  i, j-1]]
    ne = [ocean.xq[  i,   j], ocean.yq[  i,   j]]
    return sw, nw, se, ne
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

