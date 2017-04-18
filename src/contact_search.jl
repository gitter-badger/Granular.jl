## Contact mapping
function findContacts!(simulation::Simulation,
                       method::String = "all to all")

    if method == "all to all"
        findContactsAllToAll!(simulation)
    else
        error("Unknown contact search method '$method'")
    end
end

function interIceFloePositionVector(simulation::Simulation,
                                    i::Integer, j::Integer)
    return simulation.ice_floes[j].lin_pos
    - simulation.ice_floes[i].lin_pos
end

"""
position_ij is the inter-grain position vector, and can be found with
interIceFloePositionVector().
"""
function findOverlap(simulation::Simulation, i::Integer, j::Integer, 
                     position_ij::vector)
    return norm(position_ij) - (simulation.ice_floes[i].contact_radius + 
                                simulation.ice_floes[j].contact_radius)
end

function findContactsAllToAll!(simulation::Simulation)

    for i = 1:length(simulation.ice_floes)

        # Check contacts with other grains
        for j = 1:length(simulation.ice_floes)
            if i < j

                # Inter-grain position vector and grain overlap
                position_ij = interIceFloePositionVector(simulation, i, j)
                overlap_ij = findOverlap(simulation, i, j, position_ij)

                # Check if grains overlap (overlap when negative)
                if overlap_ij < 0.0
                    push!(simulation.contact_pairs, [i, j])
                    #push!(g_positions, position_ij)
                    push!(simulation.overlaps, overlap_ij)
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
function findIceFloeWallContacts(simulation::Simulation)

    for i = 1:length(g_radius)

        # Overlap with origo
        if (simulation.ice_floes[i].lin_pos[1] -
            simulation.ice_floes[i].contact_radius - simulation.origo[1]
            < 0.0)
            push!(simulation.wall_contacts, [i, -1])
        end

        if (simulation.ice_floes[i].lin_pos[2] - 
            simulation.ice_floes[i].contact_radius - simulation.origo[2] < 0.0)
            push!(simulation.ice_floes[i].wall_contacts, [i, -2])
        end

        # Overlap with world_size
        if (simulation.world_size[1] - simulation.ice_floes[i].lin_pos[1] - 
            simulation.ice_floes[i].contact_radius < 0.0)
            push!(simulation.ice_floes[i].wall_contacts, [i, 1])
        end

        if (simulation.world_size[2] - simulation.ice_floes[i].lin_pos[2] - 
            simulation.ice_floes[i].contact_radius < 0.0)
            push!(simulation.ice_floes[i].wall_contacts, [i, 2])
        end
    end
end
