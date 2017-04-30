## Contact mapping
export findContacts!
"""
    findContacts!(simulation[, method])
    
Top-level function to perform an inter-ice floe contact search, based on ice 
floe linear positions and contact radii.

The simplest contact search algorithm (`method="all to all"`) is the most 
computationally expensive (O(n^2)).  The method "ocean grid" bins the ice floes 
into their corresponding cells on the ocean grid and searches for contacts only 
within the vicinity.  When this method is applied, it is assumed that the 
`contact_radius` values of the ice floes are *smaller than half the cell size*.

# Arguments
* `simulation::Simulation`: the simulation object containing the ice floes.
* `method::String`: the contact-search method to apply.  Valid options are "all 
    to all" and "ocean grid".
"""
function findContacts!(simulation::Simulation,
                       method::String = "all to all")

    if method == "all to all"
        findContactsAllToAll!(simulation)
    elseif method == "ocean grid"
        findContactsOceanGrid!(simulation)
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

export findContactsAllToAll!
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
