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

