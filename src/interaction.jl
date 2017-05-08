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
        overlap = pop!(simulation.overlaps)
        contact_parallel_displacement = 
            pop!(simulation.contact_parallel_displacement)
        interactIceFloes!(simulation, contact_pair[1], contact_pair[2],
                          overlap, contact_parallel_displacement,
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
                           overlap::vector,
                           contact_parallel_displacement::vector;
                           contact_normal_rheology::String = "Linear Elastic",
                           contact_tangential_rheology::String = "None")

    force_n = zeros(2)
    force_t = zeros(2)

    if contact_normal_rheology == "None"
        # do nothing
    elseif contact_normal_rheology == "Linear Elastic"
        force_n = interactNormalLinearElastic(simulation, i, j, overlap)
    else
        error("Unknown contact_normal_rheology '$contact_normal_rheology'")
    end

    if contact_tangential_rheology == "None"
        # do nothing
    elseif contact_tangential_rheology == "Linear Viscous Frictional"
        force_t = interactTangentialLinearViscousFrictional(simulation, i, j,
                                                            overlap,
                                                            force_n)
    else
        error("Unknown contact_tangential_rheology ", 
              "'$contact_tangential_rheology'")
    end

    simulation.ice_floes[i].force += force_n + force_t;
    simulation.ice_floes[j].force -= force_n + force_t;

    if norm(force_t) > 0.
        torque = -findTorque(simulation, overlap, force_t, i, j)
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
                                     overlap::vector)

    k_n_harmonic_mean = 
        harmonicMean(simulation.ice_floes[i].contact_stiffness_normal,
                     simulation.ice_floes[j].contact_stiffness_normal)

    return -k_n_harmonic_mean*overlap
end

export interactTangentialLinearElastic
"""
    interactTangentialLinearViscousFrictional(simulation, i, j,
                                              overlap, force_n)

Resolves linear-viscous interaction between two ice floes in the contact- 
parallel (tangential) direction, with a frictional limit dependent on the 
Coulomb criterion.

# Arguments
* `simulation::Simulation`: the simulation object containing the ice floes.
* `i::Int`: index of the first ice floe.
* `j::Int`: index of the second ice floe.
* `overlap::vector`: two-dimensional vector pointing from i to j.
* `force_n::vector`: normal force from the interaction.
"""
function interactTangentialLinearViscousFrictional(simulation::Simulation,
                                                   i::Int, j::Integer,
                                                   overlap::vector,
                                                   force_n::vector)

    """
    contact_parallel_velocity = findContactParallelVelocity(simulation, i, j, 
                                                            overlap)



    if contact_parallel_velocity ≈ 0.
        return [0., 0.]
    end
    gamma_t_harmonic_mean = harmonicMean(
                     simulation.ice_floes[i].contact_viscosity_tangential,
                     simulation.ice_floes[j].contact_viscosity_tangential)
    mu_d_minimum = min(simulation.ice_floes[i].contact_dynamic_friction,
                       simulation.ice_floes[j].contact_dynamic_friction)

    if gamma_t_harmonic_mean ≈ 0. || mu_d_minimum ≈ 0.
        return [0., 0.]
    end

    force_t = abs(gamma_t_harmonic_mean*contact_parallel_velocity)
    if force_t > mu_d_minimum*norm(force_n)
        force_t = mu_d_minimum*norm(force_n)
    end
    if contact_parallel_velocity > 0.
        force_t = -force_t
    end

    return force_t*contactParallelVector(contactNormalVector(overlap))
    """

    p = simulation.ice_floes[i].lin_pos - simulation.ice_floes[j].lin_pos

    r_i = simulation.ice_floes[i].contact_radius
    r_j = simulation.ice_floes[j].contact_radius

    dist = norm(p)
    dn = dist - (r_i + r_j)

    if dn < 0.
        n = p/dist
        t = [-n[2], n[1]]

        vel_lin = simulation.ice_floes[i].lin_vel -
            simulation.ice_floes[j].lin_vel

        vel_n = dot(vel_lin, n)
        vel_t = dot(vel_lin, t) -
            harmonicMean(r_i, r_j)*(simulation.ice_floes[i].ang_vel +
                                    simulation.ice_floes[j].ang_vel)

        #force_n = -kn * dn - nu * vn;

        gamma_t_harmonic_mean = harmonicMean(
                     simulation.ice_floes[i].contact_viscosity_tangential,
                     simulation.ice_floes[j].contact_viscosity_tangential)

        force_t = abs(gamma_t_harmonic_mean * vel_t)

        mu_d_minimum = min(simulation.ice_floes[i].contact_dynamic_friction,
                           simulation.ice_floes[j].contact_dynamic_friction)

        if force_t > mu_d_minimum*norm(force_n)
            force_t = mu_d_minimum*norm(force_n)
        end
        if vel_t > 0.
            force_t = -force_t
        end

        return force_t*t

    end

end

function harmonicMean(a::Any, b::Any)
    hm = 2.*a*b/(a + b)
    if isnan(hm)
        return 0.
    end
    return hm
end

function findTorque(simulation::Simulation, overlap::vector, force_t::vector, 
                    i::Int, j::Int)
    return -findEffectiveRadius(simulation, i, j, overlap)*norm(force_t)
end

function findEffectiveRadius(simulation::Simulation, i::Int, j::Int,
                             overlap::vector)
    return harmonicMean(simulation.ice_floes[i].contact_radius,
                        simulation.ice_floes[j].contact_radius) -
                        norm(overlap)/2.
end

function findContactNormalVelocity(simulation::Simulation, i::Int, j::Int,
                                   overlap::vector)
    n = contactNormalVector(overlap)
    v_ij = simulation.ice_floes[i].lin_vel - simulation.ice_floes[j].lin_vel
    return v_ij[1]*n[1] + v_ij[2]*n[2]
end

function findContactParallelVelocity(simulation::Simulation, i::Int, j::Int,
                                     overlap::vector)
    v_ij = simulation.ice_floes[i].lin_vel - simulation.ice_floes[j].lin_vel
    n = contactNormalVector(overlap)
    t = contactParallelVector(n)
    return (v_ij[1]*t[1] + v_ij[2]*t[2] -
        findEffectiveRadius(simulation, i, j, overlap)*
        (simulation.ice_floes[i].ang_vel + simulation.ice_floes[j].ang_vel))
end

function contactNormalVector(overlap::vector)
    return overlap/norm(overlap)
end

function contactParallelVector(n::vector)
    return [-n[2], n[1]]
end
