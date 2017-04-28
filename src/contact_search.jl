## Contact mapping
export findContacts
"""
Top-level function to perform an inter-ice floe contact search, based on ice 
floe linear positions and contact radii.

The simplest contact search algorithm (`method="all to all"`) is the most 
computationally expensive (O(n^2)).
"""
function findContacts!(simulation::Simulation,
                       method::String = "all to all")

    if method == "all to all"
        findContactsAllToAll!(simulation)
    else
        error("Unknown contact search method '$method'")
    end
end

export interIceFloePositionVector
function interIceFloePositionVector(simulation::Simulation,
                                    i::Integer, j::Integer)
    return simulation.ice_floes[j].lin_pos - simulation.ice_floes[i].lin_pos
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

export findContactsAllToAll
"""
Perform an O(n^2) all-to-all contact search between all ice floes in the 
`simulation` object.  Contacts between fixed ice floes are ignored.
"""
function findContactsAllToAll!(simulation::Simulation)

    for i = 1:length(simulation.ice_floes)

        # Check contacts with other grains
        for j = 1:length(simulation.ice_floes)
            if i < j

                if simulation.ice_floes[i].fixed &&
                    simulation.ice_floes[j].fixed
                    continue
                end

                # Inter-grain position vector and grain overlap
                position_ij = interIceFloePositionVector(simulation, i, j)
                overlap_ij = findOverlap(simulation, i, j, position_ij)

                # Check if grains overlap (overlap when negative)
                if overlap_ij < 0.0
                    push!(simulation.contact_pairs, [i, j])
                    push!(simulation.overlaps, 
                          overlap_ij*position_ij/norm(position_ij))
                end
            end
        end
    end
end
