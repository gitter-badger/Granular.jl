## Contact mapping
function findContacts!(simulation::Simulation,
                       method::String = "all to all")

    if method == "all to all"
        findContactsAllToAll(simulation)
    else
        error("Unknown contact search method '$method'")
    end
end

function interIceFloePositionVector(i::Integer, j::Integer)
    return g_position[j]::vector - g_position[i]::vector
end

"""
position_ij is the inter-grain position vector, and can be found with
interIceFloePositionVector().
"""
function findOverlap(i::Integer, j::Integer, position_ij::vector)
    return norm(position_ij) - (g_radius[i]::float + g_radius[j]::float)
end

function findContactsAllToAll(simulation::Simulation)

    for i = 1:length(simulation.ice_floes)

        # Check contacts with other grains
        for j = 1:length(simulation.ice_floes)
            if i < j

                # Inter-grain position vector and grain overlap
                position_ij = interIceFloePositionVector(i, j)
                overlap_ij = findOverlap(i, j, position_ij)

                # Check if grains overlap (overlap when negative)
                if overlap_ij < 0.0
                    push!(g_contact_pairs, [i, j])
                    #push!(g_positions, position_ij)
                    push!(g_overlaps, overlap_ij)
                end
            end
        end
    end # `i` grain loop end
end

"""
Check contacts with world boundaries (walls). Contacts are stored in
g_wall_contacts with the format [i, w], where i is the grain number and w is the
wall number. The walls are orthorectangular. The wall at negative x is called
-1, and 1 at positive x. For negative y it is called -2, and 2 at positive y.
For negative z it is called -3, and 3 at positive z.
"""
function findIceFloeWallContacts()

    for i = 1:length(g_radius)

        # Overlap with origo
        if g_position[i][1]::float - g_radius[i]::float - g_origo[1] < 0.0
            push!(g_wall_contacts, [i, -1])
        end

        if g_position[i][2]::float - g_radius[i]::float - g_origo[2] < 0.0
            push!(g_wall_contacts, [i, -2])
        end

        if g_position[i][3]::float - g_radius[i]::float - g_origo[3] < 0.0
            push!(g_wall_contacts, [i, -3])
        end

        # Overlap with world_size
        if g_world_size[1] - g_position[i][1]::float - g_radius[i]::float < 0.0
            push!(g_wall_contacts, [i, 1])
        end

        if g_world_size[2] - g_position[i][2]::float - g_radius[i]::float < 0.0
            push!(g_wall_contacts, [i, 2])
        end

        if g_world_size[3] - g_position[i][3]::float - g_radius[i]::float < 0.0
            push!(g_wall_contacts, [i, 3])
        end

    end
end
