## Interaction functions

export interact!
"""
Resolve mechanical interaction between all particle pairs.
"""
function interact!(simulation::Simulation;
                   contact_normal_rheology::String = "Linear Elastic",
                   contact_tangential_rheology::String = "None")

    # IceFloe to grain collisions
    while !isempty(simulation.contact_pairs)
        contact_pair = pop!(simulation.contact_pairs)
        overlap_vector = pop!(simulation.overlaps)
        contact_parallel_displacement = 
            pop!(simulation.contact_parallel_displacement)
        interactIceFloes!(simulation, contact_pair[1], contact_pair[2],
                          overlap_vector, contact_parallel_displacement,
                          contact_normal_rheology=contact_normal_rheology,
                          contact_tangential_rheology=
                              contact_tangential_rheology)
    end
end

export interactIceFloes!
"""
Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the ice floe 
`pressure` field of mean compressive stress on the ice floe sides.
"""
function interactIceFloes!(simulation::Simulation,
                           i::Int, j::Int,
                           overlap_vector::Array{Float64, 1},
                           contact_parallel_displacement::Array{Float64, 1};
                           contact_normal_rheology::String = "Linear Elastic",
                           contact_tangential_rheology::String = "None")

    force_n = zeros(2)
    force_t = zeros(2)

    if contact_normal_rheology == "None"
        # do nothing
    elseif contact_normal_rheology == "Linear Elastic"
        force_n = interactNormalLinearElastic(simulation, i, j, overlap_vector)
    else
        error("Unknown contact_normal_rheology '$contact_normal_rheology'")
    end

    if contact_tangential_rheology == "None"
        # do nothing
    elseif contact_tangential_rheology == "Linear Viscous Frictional"
        force_t = interactTangentialLinearViscousFrictional(simulation, i, j,
                                                            overlap_vector,
                                                            force_n)
    else
        error("Unknown contact_tangential_rheology ", 
              "'$contact_tangential_rheology'")
    end

    simulation.ice_floes[i].force += force_n + force_t;
    simulation.ice_floes[j].force -= force_n + force_t;

    if norm(force_t) > 0.
        torque = findTorque(simulation, overlap_vector, force_t, i, j)
        simulation.ice_floes[i].torque += torque
        simulation.ice_floes[j].torque += torque
    end

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
                                     i::Int, j::Int,
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
function interactTangentialLinearViscousFrictional(simulation::Simulation,
                                                   i::Int, j::Integer,
                                                   overlap_vector::vector,
                                                   force_n::vector)

    contact_parallel_velocity = findContactParallelVelocity(simulation, i, j, 
                                                            overlap_vector)

    if norm(contact_parallel_velocity) ≈ 0.
        return [0., 0.]
    end

    gamma_t_harmonic_mean = harmonicMean(
                     simulation.ice_floes[i].contact_viscosity_tangential,
                     simulation.ice_floes[j].contact_viscosity_tangential)

    if norm(gamma_t_harmonic_mean) ≈ 0.
        return [0., 0.]
    end

    mu_d_minimum = min(simulation.ice_floes[i].contact_dynamic_friction,
                       simulation.ice_floes[j].contact_dynamic_friction)

    return -min(norm(gamma_t_harmonic_mean*contact_parallel_velocity), 
                mu_d_minimum*norm(force_n))*
        contact_parallel_velocity/norm(contact_parallel_velocity)
end

function harmonicMean(a::Any, b::Any)
    hm = 2.*a*b/(a + b)
    if isnan(hm)
        return 0.
    end
    return hm
end

function findTorque(simulation::Simulation, overlap_vector::vector, 
                    force_t::vector, i::Int, j::Int)
    n = overlap_vector/norm(overlap_vector)
    return -findEffectiveRadius(simulation, i, j, overlap_vector)*
        (n[1]*force_t[2] - n[2]*force_t[1])
end

function findEffectiveRadius(simulation::Simulation, i::Int, j::Int,
                             overlap_vector::vector)
    return harmonicMean(simulation.ice_floes[i].contact_radius,
                        simulation.ice_floes[j].contact_radius)
                        - norm(overlap_vector)/2.
end

function findContactParallelVelocity(simulation::Simulation, i::Int, j::Int,
                                     overlap_vector::vector)
    return simulation.ice_floes[i].lin_vel -
        simulation.ice_floes[j].lin_vel +
        findEffectiveRadius(simulation, i, j, overlap_vector)*
        (simulation.ice_floes[i].ang_vel + simulation.ice_floes[j].ang_vel)
end
