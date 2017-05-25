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
function findContacts!(simulation::Simulation;
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
"""
    interIceFloePositionVector(simulation, i, j)

Returns a `vector` pointing from ice floe `i` to ice floe `j` in the 
`simulation`.

# Arguments
* `simulation::Simulation`: the simulation object containing the ice floes.
* `i::Int`: index of the first ice floe.
* `j::Int`: index of the second ice floe.
"""
function interIceFloePositionVector(simulation::Simulation,
                                    i::Integer, j::Integer)
    return simulation.ice_floes[i].lin_pos - simulation.ice_floes[j].lin_pos
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
    findContactsAllToAll!(simulation)

Perform an O(n^2) all-to-all contact search between all ice floes in the 
`simulation` object.  Contacts between fixed ice floes are ignored.
"""
function findContactsAllToAll!(simulation::Simulation)

    for i = 1:length(simulation.ice_floes)

        # Check contacts with other grains
        for j = 1:length(simulation.ice_floes)
            checkAndAddContact!(simulation, i, j)
        end
    end
end

export findContactsOceanGrid!
"""
    findContactsOceanGrid!(simulation)

Perform an O(n*log(n)) cell-based contact search between all ice floes in the 
`simulation` object.  Contacts between fixed ice floes are ignored.
"""
function findContactsOceanGrid!(simulation::Simulation)

    for idx_i = 1:length(simulation.ice_floes)

        grid_pos = simulation.ice_floes[idx_i].ocean_grid_pos
        nx, ny = size(simulation.ocean.xh)

        for i=(grid_pos[1] - 1):(grid_pos[1] + 1)
            for j=(grid_pos[2] - 1):(grid_pos[2] + 1)

                # only check for contacts within grid boundaries
                if i < 1 || i > nx || j < 1 || j > ny
                    continue
                end

                for idx_j in simulation.ocean.ice_floe_list[i, j]
                    checkAndAddContact!(simulation, idx_i, idx_j)
                end
            end
        end
    end
end

export addContact!
"""
    checkAndAddContact!(simulation, i, j)

Check for contact between two ice floes and register the interaction in the 
`simulation` object.  The indexes of the two ice floes is stored in 
`simulation.contact_pairs` as `[i, j]`.  The overlap vector is parallel to a 
straight line connecting the ice floe centers, points away from ice floe `i` and 
towards `j`, and is stored in `simulation.overlaps`.  A zero-length vector is 
written to `simulation.contact_parallel_displacement`.

# Arguments
* `simulation::Simulation`: the simulation object containing the ice floes.
* `i::Int`: index of the first ice floe.
* `j::Int`: index of the second ice floe.
"""
function checkAndAddContact!(sim::Simulation, i::Int, j::Int)
    if i < j

        if (sim.ice_floes[i].fixed && sim.ice_floes[j].fixed) ||
            !sim.ice_floes[i].enabled || !sim.ice_floes[j].enabled
            return
        end

        # Inter-grain position vector and grain overlap
        position_ij = interIceFloePositionVector(sim, i, j)
        overlap_ij = findOverlap(sim, i, j, position_ij)

        # Check if grains overlap (overlap when negative)
        if overlap_ij < 0.
            for ic=1:(Nc_max + 1)
                if ic == (Nc_max + 1)
                    error("contact $i-$j exceeds max. number of contacts " *
                          "(Nc_max = $Nc_max) for ice floe $i")

                else
                    if sim.ice_floes[i].contacts[ic] == j
                        break  # contact already registered

                    elseif sim.ice_floes[i].contacts[ic] == 0  # empty
                        sim.ice_floes[i].n_contacts += 1  # register new contact
                        sim.ice_floes[j].n_contacts += 1
                        sim.ice_floes[i].contacts[ic] = j
                        sim.ice_floes[i].contact_parallel_displacement[ic] =
                            zeros(2)
                        sim.ice_floes[i].contact_age[ic] = 0.
                        break
                    end
                end
            end
        end
    end
end
