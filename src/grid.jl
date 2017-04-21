function fitWorldSize(margin::vector = zeros(3))

    if length(g_radius) < 1
        error("At least a single grain must be added before setting world
        size.")
    end

    minpos = [Inf, Inf, Inf]
    maxpos = [-Inf, -Inf, -Inf]
    for i = 1:length(g_radius)

        min_position = g_position[i]::vector - g_radius[i]::float
        max_position = g_position[i]::vector + g_radius[i]::float

        for i = 1:3
            if min_position[i] < g_origo[i]
                g_origo[i] = min_position[i]
            end

            if max_position[i] > g_world_size[i]
                g_world_size[i] = max_position[i]
            end
        end
    end

    g_origo[:] = minpos[:] - margin[:]
    g_world_size[:] = maxpos[:] + margin[:]
end

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
