## Interaction functions

export interact!
"""
Resolve mechanical interaction between all particle pairs.
"""
function interact!(simulation::Simulation)

    # IceFloe to grain collisions
    while !isempty(simulation.contact_pairs)
        contact_pair = pop!(simulation.contact_pairs)
        overlap_vector = pop!(simulation.overlaps)
        contact_parallel_displacement = 
            pop!(simulation.contact_parallel_displacement)
        interactIceFloes!(simulation, contact_pair[1], contact_pair[2],
                          overlap_vector, contact_parallel_displacement)
    end
end

export interactIceFloes!
"""
Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the ice floe 
`pressure` field of mean compressive stress on the ice floe sides.
"""
function interactIceFloes!(simulation::Simulation,
                           i::Integer, j::Integer,
                           overlap_vector::Array{Float64, 1},
                           contact_parallel_displacement::Array{Float64, 1};
                           contact_normal::String = "LinearElasticNoTangential")

    force = zeros(2)

    if contact_normal == "None"
        return nothing

    elseif contact_normal == "LinearElasticNoTangential" ||
        contact_normal == "LinearElastic"
        force_n = interactNormalLinearElastic(simulation, i, j, overlap_vector)

        if contact_normal == "LinearElastic"
            force_t, torque = interactTangentialLinearElastic(simulation, i, j,
                                                              overlap_vector,
                                                              contact_parallel_displacement)
        end

    else
        error("Unknown contact_normal interaction model '$contact_normal'")
    end

    simulation.ice_floes[i].force += force_n;
    simulation.ice_floes[j].force -= force_n;

    simulation.ice_floes[i].pressure += 
        norm(force_n)/simulation.ice_floes[i].side_surface_area;
    simulation.ice_floes[j].pressure += 
        norm(force_n)/simulation.ice_floes[j].side_surface_area;
end

export interactNormalLinearElastic
"""
Resolves linear-elastic interaction between two ice floes in the contact-normal 
direction.
"""
function interactNormalLinearElastic(simulation::Simulation,
                                     i::Integer, j::Integer,
                                     overlap_vector::vector)

    k_n_harmonic_mean = 
        harmonicMean(simulation.ice_floes[i].contact_stiffness_normal,
                     simulation.ice_floes[j].contact_stiffness_normal)

    return k_n_harmonic_mean * overlap_vector
end

export interactTangentialLinearElastic
"""
Resolves linear-elastic interaction between two ice floes in the 
contact-parallel (tangential) direction.
"""
function interactTangentialLinearElastic(simulation::Simulation,
                                     i::Integer, j::Integer,
                                     overlap_vector::vector,
                                     contact_parallel_displacement::vector)

    k_t_harmonic_mean = 
        harmonicMean(simulation.ice_floes[i].contact_stiffness_tangential,
                     simulation.ice_floes[j].contact_stiffness_tangential)

    # correct displacement for contact rotation
    n = overlap_vector/norm(overlap_vector)
    contact_parallel_displacement -= (n * dot(n, contact_parallel_displacement))


    force_t = k_t_harmonic_mean * contact_parallel_displacement

    return force, torque, contact_parallel_displacement

end

function harmonicMean(a::Any, b::Any)
    return 2.*a*b/(a + b)
end
