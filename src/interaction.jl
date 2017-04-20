## Interaction functions

"""
Resolve mechanical interaction between all grain pairs and walls.
"""
function interact!(simulation::Simulation)

    # IceFloe to grain collisions
    while !isempty(simulation.contact_pairs)
        contact_pair = pop!(simulation.contact_pairs)
        overlap_vector = pop!(simulation.overlaps)
        interactIceFloes!(simulation, contact_pair[1], contact_pair[2],
                          overlap_vector)
    end

    # IceFloe to wall collisions
    while !isempty(simulation.wall_contacts)
        contact_pair = pop!(simulation.wall_contacts)
        interactIceFloeWall!(simulation, contact_pair[1], contact_pair[2])
    end
end

"""
Resolve an grain-to-grain interaction using a prescibed contact law.
"""
function interactIceFloes!(simulation::Simulation,
                           i::Integer, j::Integer,
                           overlap_vector::Array{Float64, 1};
                           contact_normal::String = "LinearElastic")

    force = zeros(2)

    if contact_normal == "None"
        # do nothing

    elseif contact_normal == "LinearElastic"
        force = interactNormalLinearElastic(simulation, i, j, overlap_vector)

    else
        error("Unknown contact_normal interaction model '$contact_normal'")
    end

    simulation.ice_floes[i].force += force;
    simulation.ice_floes[j].force -= force;

end

function interactNormalLinearElastic(simulation::Simulation,
                                     i::Integer, j::Integer,
                                     overlap_vector::vector)

    k_n_i = simulation.ice_floes[i].contact_stiffness_normal
    k_n_j = simulation.ice_floes[j].contact_stiffness_normal
    k_n_harmonic_mean = 2.*k_n_i*k_n_j/(k_n_i + k_n_j)

    return k_n_harmonic_mean * overlap_vector
end
