## Contact mapping
export findContacts!
"""
    findContacts!(simulation[, method])
    
Top-level function to perform an inter-grain contact search, based on grain 
linear positions and contact radii.

The simplest contact search algorithm (`method="all to all"`) is the most 
computationally expensive (O(n^2)).  The method "ocean grid" bins the grains 
into their corresponding cells on the ocean grid and searches for contacts only 
within the vicinity.  When this method is applied, it is assumed that the 
`contact_radius` values of the grains are *smaller than half the cell size*.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
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

export interGrainPositionVector
"""
    interGrainPositionVector(simulation, i, j)

Returns a `vector` pointing from grain `i` to grain `j` in the 
`simulation`.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
* `i::Int`: index of the first grain.
* `j::Int`: index of the second grain.
"""
function interGrainPositionVector(simulation::Simulation,
                                    i::Int, j::Int)
    @inbounds return simulation.grains[i].lin_pos - 
    simulation.grains[j].lin_pos
end

"""
position_ij is the inter-grain position vector, and can be found with
interGrainPositionVector().
"""
function findOverlap(simulation::Simulation, i::Int, j::Int, 
                     position_ij::Vector{Float64})
    @inbounds return norm(position_ij) - (simulation.grains[i].contact_radius 
                                + simulation.grains[j].contact_radius)
end

export findContactsAllToAll!
"""
    findContactsAllToAll!(simulation)

Perform an O(n^2) all-to-all contact search between all grains in the 
`simulation` object.  Contacts between fixed grains are ignored.
"""
function findContactsAllToAll!(simulation::Simulation)

    if simulation.ocean.bc_west > 1 ||
        simulation.ocean.bc_east > 1 ||
        simulation.ocean.bc_north > 1 ||
        simulation.ocean.bc_south > 1
        error("Ocean boundary conditions to not work with all-to-all contact " *
              "search")
    end
    if simulation.atmosphere.bc_west > 1 ||
        simulation.atmosphere.bc_east > 1 ||
        simulation.atmosphere.bc_north > 1 ||
        simulation.atmosphere.bc_south > 1
        error("Atmopshere boundary conditions to not work with all-to-all " *
              "contact search")
    end

    @inbounds for i = 1:length(simulation.grains)

        # Check contacts with other grains
        for j = 1:length(simulation.grains)
            checkAndAddContact!(simulation, i, j)
        end
    end
    nothing
end

export findContactsInGrid!
"""
    findContactsInGrid!(simulation)

Perform an O(n*log(n)) cell-based contact search between all grains in the 
`simulation` object.  Contacts between fixed or disabled grains are ignored.
"""
function findContactsInGrid!(simulation::Simulation, grid::Any)

    distance_modifier = [0., 0.]
    i_corrected = 0
    j_corrected = 0

    for idx_i = 1:length(simulation.grains)

        if typeof(grid) == Ocean
            grid_pos = simulation.grains[idx_i].ocean_grid_pos
        elseif typeof(grid) == Atmosphere
            grid_pos = simulation.grains[idx_i].atmosphere_grid_pos
        else
            error("Grid type not understood")
        end
        nx, ny = size(grid.xh)

        for i=(grid_pos[1] - 1):(grid_pos[1] + 1)
            for j=(grid_pos[2] - 1):(grid_pos[2] + 1)

                # i and j are not corrected for periodic boundaries
                i_corrected = i
                j_corrected = j

                # vector for correcting inter-particle distance in case of
                # boundary periodicity
                distance_modifier .= [0., 0.]

                # only check for contacts within grid boundaries, and wrap
                # around if they are periodic
                if i < 1 || i > nx || j < 1 || j > ny

                    if i < 1 && grid.bc_west == 2  # periodic -x
                        distance_modifier[1] = grid.xq[end] - grid.xq[1]
                        i_corrected = nx
                    elseif i > nx && grid.bc_east == 2  # periodic +x
                        distance_modifier[1] = -(grid.xq[end] - grid.xq[1])
                        i_corrected = 1
                    end

                    if j < 1 && grid.bc_south == 2  # periodic -y
                        distance_modifier[2] = grid.yq[end] - grid.yq[1]
                        j_corrected = ny
                    elseif j > ny && grid.bc_north == 2  # periodic +y
                        distance_modifier[2] = -(grid.yq[end] - grid.yq[1])
                        j_corrected = 1
                    end

                    # skip iteration if target still falls outside grid after
                    # periodicity correction
                    if i_corrected < 1 || i_corrected > nx ||
                        j_corrected < 1 || j_corrected > ny
                        continue
                    end
                end

                @inbounds for idx_j in grid.grain_list[i_corrected, j_corrected]
                    checkAndAddContact!(simulation, idx_i, idx_j,
                                        distance_modifier)
                end
            end
        end
    end
    nothing
end

export checkAndAddContact!
"""
    checkAndAddContact!(simulation, i, j)

Check for contact between two grains and register the interaction in the 
`simulation` object.  The indexes of the two grains is stored in 
`simulation.contact_pairs` as `[i, j]`.  The overlap vector is parallel to a 
straight line connecting the grain centers, points away from grain `i` and 
towards `j`, and is stored in `simulation.overlaps`.  A zero-length vector is 
written to `simulation.contact_parallel_displacement`.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
* `i::Int`: index of the first grain.
* `j::Int`: index of the second grain.
* `distance_Modifier::Vector{Float64}`: vector modifying percieved
    inter-particle distance, which is used for contact search across periodic
    boundaries.
"""
function checkAndAddContact!(sim::Simulation, i::Int, j::Int,
                             distance_modifier::Vector{Float64} = [0., 0.])
    if i < j

        @inbounds if (sim.grains[i].fixed && sim.grains[j].fixed) ||
            !sim.grains[i].enabled || !sim.grains[j].enabled
            return
        end

        # Inter-grain position vector and grain overlap
        position_ij = interGrainPositionVector(sim, i, j) + distance_modifier
        overlap_ij = findOverlap(sim, i, j, position_ij)

        contact_found = false

        # Check if contact is already registered, and update position if so
        for ic=1:sim.Nc_max
            @inbounds if sim.grains[i].contacts[ic] == j
                contact_found = true
                @inbounds sim.grains[i].position_vector[ic] .= position_ij
                nothing  # contact already registered
            end
        end

        # Check if grains overlap (overlap when negative)
        if overlap_ij < 0.

            # Register as new contact in first empty position
            if !contact_found

                for ic=1:(sim.Nc_max + 1)

                    # Test if this contact exceeds the number of contacts
                    if ic == (sim.Nc_max + 1)
                        for ic=1:sim.Nc_max
                            warn("grains[$i].contacts[$ic] = " *
                                 "$(sim.grains[i].contacts[ic])")
                            warn("grains[$i].contact_age[$ic] = " *
                                 "$(sim.grains[i].contact_age[ic])")
                        end
                        error("contact $i-$j exceeds max. number of contacts " *
                              "(sim.Nc_max = $(sim.Nc_max)) for grain $i")
                    end

                    # Register as new contact
                    @inbounds if sim.grains[i].contacts[ic] == 0  # empty
                        @inbounds sim.grains[i].n_contacts += 1
                        @inbounds sim.grains[j].n_contacts += 1
                        @inbounds sim.grains[i].contacts[ic] = j
                        @inbounds sim.grains[i].position_vector[ic] .=
                            position_ij
                        @inbounds fill!(sim.grains[i].
                              contact_parallel_displacement[ic] , 0.)
                        @inbounds sim.grains[i].contact_age[ic] = 0.
                        break
                    end
                end
            end
        end
    end
    nothing
end
