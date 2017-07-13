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
        findContactsInGrid!(simulation, simulation.ocean)

    elseif method == "atmosphere grid"
        findContactsInGrid!(simulation, simulation.atmosphere)

    else
        error("Unknown contact search method '$method'")
    end
    nothing
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
                                    i::Int, j::Int)
    @inbounds return simulation.ice_floes[i].lin_pos - 
    simulation.ice_floes[j].lin_pos
end

"""
position_ij is the inter-grain position vector, and can be found with
interIceFloePositionVector().
"""
function findOverlap(simulation::Simulation, i::Int, j::Int, 
                     position_ij::Vector{Float64})
    @inbounds return norm(position_ij) - (simulation.ice_floes[i].contact_radius 
                                + simulation.ice_floes[j].contact_radius)
end

export findContactsAllToAll!
"""
    findContactsAllToAll!(simulation)

Perform an O(n^2) all-to-all contact search between all ice floes in the 
`simulation` object.  Contacts between fixed ice floes are ignored.
"""
function findContactsAllToAll!(simulation::Simulation)

    @inbounds for i = 1:length(simulation.ice_floes)

        # Check contacts with other grains
        for j = 1:length(simulation.ice_floes)
            checkAndAddContact!(simulation, i, j)
        end
    end
    nothing
end

export findContactsInGrid!
"""
    findContactsInGrid!(simulation)

Perform an O(n*log(n)) cell-based contact search between all ice floes in the 
`simulation` object.  Contacts between fixed or disabled ice floes are ignored.
"""
function findContactsInGrid!(simulation::Simulation, grid::Any)

    for idx_i = 1:length(simulation.ice_floes)

        if typeof(grid) == Ocean
            grid_pos = simulation.ice_floes[idx_i].ocean_grid_pos
        elseif typeof(grid) == Atmosphere
            grid_pos = simulation.ice_floes[idx_i].atmosphere_grid_pos
        else
            error("grid type not understood")
        end
        nx, ny = size(grid.xh)

        for i=(grid_pos[1] - 1):(grid_pos[1] + 1)
            for j=(grid_pos[2] - 1):(grid_pos[2] + 1)

                # only check for contacts within grid boundaries
                if i < 1 || i > nx || j < 1 || j > ny
                    continue
                end

                @inbounds for idx_j in grid.ice_floe_list[i, j]
                    checkAndAddContact!(simulation, idx_i, idx_j)
                end
            end
        end
    end
    nothing
end

export checkAndAddContact!
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

        @inbounds if (sim.ice_floes[i].fixed && sim.ice_floes[j].fixed) ||
            !sim.ice_floes[i].enabled || !sim.ice_floes[j].enabled
            return
        end

        # Inter-grain position vector and grain overlap
        position_ij = interIceFloePositionVector(sim, i, j)
        overlap_ij = findOverlap(sim, i, j, position_ij)

        # Check if grains overlap (overlap when negative)
        if overlap_ij < 0.
            for ic=1:(sim.Nc_max + 1)
                if ic == (sim.Nc_max + 1)
                    error("contact $i-$j exceeds max. number of contacts " *
                          "(sim.Nc_max = $(sim.Nc_max)) for ice floe $i")

                else
                    @inbounds if sim.ice_floes[i].contacts[ic] == j
                        break  # contact already registered

                    elseif sim.ice_floes[i].contacts[ic] == 0  # empty
                        @inbounds sim.ice_floes[i].n_contacts += 1
                        @inbounds sim.ice_floes[j].n_contacts += 1
                        @inbounds sim.ice_floes[i].contacts[ic] = j
                        @inbounds fill!(sim.ice_floes[i].
                              contact_parallel_displacement[ic] , 0.)
                        @inbounds sim.ice_floes[i].contact_age[ic] = 0.
                        break
                    end
                end
            end
        end
    end
    nothing
end
