## Interaction functions

"""
Resolve mechanical interaction between all particle pairs.
"""
function interact!(simulation::Simulation)

    # IceFloe to grain collisions
    while !isempty(simulation.contact_pairs)
        contact_pair = pop!(simulation.contact_pairs)
        overlap_vector = pop!(simulation.overlaps)
        interactIceFloes!(simulation, contact_pair[1], contact_pair[2],
                          overlap_vector)
    end
end

"""
Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the ice floe 
`pressure` field of mean compressive stress on the ice floe sides.
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

    simulation.ice_floes[i].pressure += 
        norm(force)/simulation.ice_floes[i].side_surface_area;
    simulation.ice_floes[j].pressure += 
        norm(force)/simulation.ice_floes[j].side_surface_area;
end

"""
Resolves linear-elastic interaction between two ice floes in the contact-normal 
direction.
"""
function interactNormalLinearElastic(simulation::Simulation,
                                     i::Integer, j::Integer,
                                     overlap_vector::vector)

    k_n_i = simulation.ice_floes[i].contact_stiffness_normal
    k_n_j = simulation.ice_floes[j].contact_stiffness_normal
    k_n_harmonic_mean = 2.*k_n_i*k_n_j/(k_n_i + k_n_j)

    return k_n_harmonic_mean * overlap_vector
end
