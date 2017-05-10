## Interaction functions

export interact!
"""
    interact!(simulation::Simulation)

Resolve mechanical interaction between all particle pairs.
"""
function interact!(simulation::Simulation)
    for i=1:Int(ceil(length(simulation.ice_floes)/2.))  # i <= Int(N/2)
        for ic=1:Nc_max

            j = simulation.ice_floes[i].contacts[ic]

            if j == 0
                break  # end of contact list reached
            end

            if norm(simulation.ice_floes[i].lin_pos - 
                    simulation.ice_floes[j].lin_pos) - 
                (simulation.ice_floes[i].contact_radius + 
                 simulation.ice_floes[j].contact_radius) > 0.

                simulation.ice_floes[i].contacts[ic] = 0  # remove contact
                simulation.ice_floes[i].n_contacts -= 1
            else
                interactIceFloes!(simulation, i, j, ic)
            end
        end
    end
end

export interactIceFloes!
"""
Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the ice floe 
`pressure` field of mean compressive stress on the ice floe sides.
"""
function interactIceFloes!(simulation::Simulation, i::Int, j::Int, ic::Int)

    force_n = 0.  # Contact-normal force
    force_t = 0.  # Contact-parallel (tangential) force

    # Inter-position vector
    p = simulation.ice_floes[i].lin_pos - simulation.ice_floes[j].lin_pos
    dist = norm(p)

    r_i = simulation.ice_floes[i].contact_radius
    r_j = simulation.ice_floes[j].contact_radius

    # Floe distance
    delta_n = dist - (r_i + r_j)

    if delta_n < 0.  # Contact (this should always occur)

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

        # Correct old tangential displacement for contact rotation and add new
        delta_t0 =simulation.ice_floes[i].contact_parallel_displacement[ic]
        delta_t = dot(t, delta_t0 - (n*dot(n, delta_t0))) +
            vel_t*simulation.time_step

        # Effective radius
        R_ij = harmonicMean(simulation.ice_floes[i].contact_radius,
                            simulation.ice_floes[j].contact_radius) - 
            abs(delta_n)/2.

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
        if k_n_harmonic_mean ≈ 0. && gamma_n_harmonic_mean ≈ 0.
            force_n = 0.

        elseif k_n_harmonic_mean > 0. && gamma_n_harmonic_mean ≈ 0.
            force_n = -k_n_harmonic_mean*delta_n

        elseif k_n_harmonic_mean > 0. && gamma_n_harmonic_mean > 0.
            force_n = -k_n_harmonic_mean*delta_n - gamma_n_harmonic_mean*vel_n
            if force_n < 0.
                force_n = 0.
            end

        else
            error("unknown contact_normal_rheology (k_n = $k_n_harmonic_mean," *
                  " gamma_n = $gamma_n_harmonic_mean")
        end

        if k_t_harmonic_mean ≈ 0. && gamma_t_harmonic_mean ≈ 0.
            # do nothing

        elseif k_t_harmonic_mean ≈ 0. && gamma_t_harmonic_mean > 0.
            force_t = abs(gamma_t_harmonic_mean * vel_t)
            if force_t > mu_d_minimum*abs(force_n)
                force_t = mu_d_minimum*abs(force_n)
            end
            if vel_t > 0.
                force_t = -force_t
            end

        elseif k_t_harmonic_mean > 0.

                force_t = -k_t_harmonic_mean*delta_t -
                    gamma_t_harmonic_mean*vel_t

            if abs(force_t) > mu_d_minimum*abs(force_n)
                force_t = mu_d_minimum*abs(force_n)*force_t/abs(force_t)
                delta_t = (-force_t - gamma_t_harmonic_mean*vel_t)/
                    k_t_harmonic_mean
            end

        else
            error("unknown contact_tangential_rheology (k_t = " *
                  "$k_t_harmonic_mean, gamma_t = $gamma_t_harmonic_mean")
        end

    else
        error("function called to process non-existent contact between ice " *
              "floes $i and $j")
    end
    simulation.ice_floes[i].contact_parallel_displacement[ic] = delta_t*t

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
