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

    force_n = 0.
    force_t = 0.

    p = simulation.ice_floes[i].lin_pos - simulation.ice_floes[j].lin_pos
    dist = norm(p)

    r_i = simulation.ice_floes[i].contact_radius
    r_j = simulation.ice_floes[j].contact_radius

    dn = dist - (r_i + r_j)

    if dn < 0.  # Contact (this should always occur)

        # Local axes
        n = p/dist
        t = [-n[2], n[1]]

        # Contact kinematics
        vel_lin = simulation.ice_floes[i].lin_vel -
            simulation.ice_floes[j].lin_vel
        vel_n = dot(vel_lin, n)
        vel_t = dot(vel_lin, t) -
            harmonicMean(r_i, r_j)*(simulation.ice_floes[i].ang_vel +
                                    simulation.ice_floes[j].ang_vel)

        R_ij = harmonicMean(simulation.ice_floes[i].contact_radius,
                            simulation.ice_floes[j].contact_radius) -
               norm(overlap)/2.

        # Contact mechanical parameters
        k_n_harmonic_mean = 
            harmonicMean(simulation.ice_floes[i].contact_stiffness_normal,
                         simulation.ice_floes[j].contact_stiffness_normal)

        k_t_harmonic_mean = 
            harmonicMean(simulation.ice_floes[i].contact_stiffness_tangential,
                         simulation.ice_floes[j].contact_stiffness_tangential)

        gamma_n_harmonic_mean = harmonicMean(
                     simulation.ice_floes[i].contact_viscosity_normal,
                     simulation.ice_floes[j].contact_viscosity_normal)

        gamma_t_harmonic_mean = harmonicMean(
                     simulation.ice_floes[i].contact_viscosity_tangential,
                     simulation.ice_floes[j].contact_viscosity_tangential)

        mu_d_minimum = min(simulation.ice_floes[i].contact_dynamic_friction,
                           simulation.ice_floes[j].contact_dynamic_friction)

        # Determine contact forces
        if contact_normal_rheology == "None"
            force_n = 0.
        elseif contact_normal_rheology == "Linear Elastic"
            force_n = -k_n_harmonic_mean*dn
        elseif contact_normal_rheology == "Linear Viscous Elastic"
            force_n = -k_n_harmonic_mean*dn - gamma_n_harmonic_mean*vel_n
            if force_n < 0.
                force_n = 0.
            end
        else
            error("unknown contact_normal_rheology '$contact_normal_rheology'")
        end

        if contact_tangential_rheology == "None"
            # do nothing
        elseif contact_tangential_rheology == "Linear Elastic Frictional"
            error("not yet implemented")
        elseif contact_tangential_rheology == "Linear Viscous Frictional"
            force_t = abs(gamma_t_harmonic_mean * vel_t)
            if force_t > mu_d_minimum*norm(force_n)
                force_t = mu_d_minimum*norm(force_n)
            end
            if vel_t > 0.
                force_t = -force_t
            end
        else
            error("unknown contact_tangential_rheology ", 
                  "'$contact_tangential_rheology'")
        end

    else
        error("function called to process non-existent contact between ice " *
              "floes $i and $j")
    end

    simulation.ice_floes[i].force += force_n*n + force_t*t;
    simulation.ice_floes[j].force -= force_n*n + force_t*t;

    simulation.ice_floes[i].torque += -force_t*R_ij
    simulation.ice_floes[j].torque += -force_t*R_ij

    simulation.ice_floes[i].pressure += 
        force_n/simulation.ice_floes[i].side_surface_area;
    simulation.ice_floes[j].pressure += 
        force_n/simulation.ice_floes[j].side_surface_area;
end

function harmonicMean(a::Any, b::Any)
    hm = 2.*a*b/(a + b)
    if isnan(hm)
        return 0.
    end
    return hm
end
