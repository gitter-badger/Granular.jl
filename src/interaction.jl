## Interaction functions

export interact!
"""
    interact!(simulation::Simulation)

Resolve mechanical interaction between all particle pairs.
"""
function interact!(simulation::Simulation)
    for i=1:length(simulation.ice_floes)
        for ic=1:simulation.Nc_max

            @inbounds j = simulation.ice_floes[i].contacts[ic]

            if i > j  # skip i > j and j == 0
                continue
            end

            """
            if norm(simulation.ice_floes[i].lin_pos - 
                    simulation.ice_floes[j].lin_pos) - 
                (simulation.ice_floes[i].contact_radius + 
                 simulation.ice_floes[j].contact_radius) > 0.

                simulation.ice_floes[i].contacts[ic] = 0  # remove contact
                simulation.ice_floes[i].n_contacts -= 1
                simulation.ice_floes[j].n_contacts -= 1
            else
            """
            interactIceFloes!(simulation, i, j, ic)
            #end
        end
    end

    for i=1:length(simulation.ice_floes)
        @inbounds simulation.ice_floes[i].granular_stress = 
            simulation.ice_floes[i].force/
            simulation.ice_floes[i].horizontal_surface_area
    end
    nothing
end

export interactIceFloes!
"""
    interactIceFloes!(simulation::Simulation, i::Int, j::Int, ic::Int)

Resolve an grain-to-grain interaction using a prescibed contact law.  This 
function adds the compressive force of the interaction to the ice floe 
`pressure` field of mean compressive stress on the ice floe sides.

The function uses the macroscopic contact-stiffness parameterization based on 
Young's modulus and Poisson's ratio if Young's modulus is a positive value.
"""
function interactIceFloes!(simulation::Simulation, i::Int, j::Int, ic::Int)

    if !simulation.ice_floes[i].enabled || !simulation.ice_floes[j].enabled
        return
    end

    force_n = 0.  # Contact-normal force
    force_t = 0.  # Contact-parallel (tangential) force

    # Inter-position vector
    p = simulation.ice_floes[i].lin_pos - simulation.ice_floes[j].lin_pos
    dist = norm(p)

    r_i = simulation.ice_floes[i].contact_radius
    r_j = simulation.ice_floes[j].contact_radius

    # Floe distance; <0: compression, >0: tension
    δ_n = dist - (r_i + r_j)

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
    δ_t0 =simulation.ice_floes[i].contact_parallel_displacement[ic]
    δ_t = dot(t, δ_t0 - (n*dot(n, δ_t0))) + vel_t*simulation.time_step

    # Effective radius
    R_ij = harmonicMean(r_i, r_j)

    # Contact area
    A_ij = R_ij*min(simulation.ice_floes[i].thickness, 
                    simulation.ice_floes[j].thickness)

    # Contact mechanical parameters
    if simulation.ice_floes[i].youngs_modulus > 0. &&
        simulation.ice_floes[j].youngs_modulus > 0.

        E = harmonicMean(simulation.ice_floes[i].youngs_modulus,
                         simulation.ice_floes[j].youngs_modulus)
        ν = harmonicMean(simulation.ice_floes[i].poissons_ratio,
                         simulation.ice_floes[j].poissons_ratio)

        # Effective normal and tangential stiffness
        k_n = E*A_ij/R_ij
        #k_t = k_n*ν   # Kneib et al 2016
        k_t = k_n*2.*(1. - ν^2.)/((2. - ν)*(1. + ν))  # Obermayr et al 2011

    else  # Micromechanical parameterization: k_n and k_t set explicitly
        k_n = harmonicMean(simulation.ice_floes[i].contact_stiffness_normal,
                         simulation.ice_floes[j].contact_stiffness_normal)

        k_t = harmonicMean(simulation.ice_floes[i].contact_stiffness_tangential,
                           simulation.ice_floes[j].contact_stiffness_tangential)
    end

    γ_n = harmonicMean(simulation.ice_floes[i].contact_viscosity_normal,
                           simulation.ice_floes[j].contact_viscosity_normal)

    γ_t = harmonicMean(simulation.ice_floes[i].contact_viscosity_tangential,
                       simulation.ice_floes[j].contact_viscosity_tangential)

    μ_d_minimum = min(simulation.ice_floes[i].contact_dynamic_friction,
                       simulation.ice_floes[j].contact_dynamic_friction)

    # Determine contact forces
    if k_n ≈ 0. && γ_n ≈ 0.
        force_n = 0.

    elseif k_n > 0. && γ_n ≈ 0.
        force_n = -k_n*δ_n

    elseif k_n > 0. && γ_n > 0.
        force_n = -k_n*δ_n - γ_n*vel_n
        if force_n < 0.
            force_n = 0.
        end

    else
        error("unknown contact_normal_rheology (k_n = $k_n, γ_n = $γ_n")
    end

    # Contact tensile strength increases linearly with contact age until tensile 
    # stress exceeds tensile strength
    if δ_n > 0.

        # linearly increase tensile strength with time until max. value
        tensile_strength = min(simulation.ice_floes[i].contact_age[ic]/
                               max(simulation.ice_floes[i].tensile_heal_rate,
                                   1e-12), 1.)*
                               simulation.ice_floes[i].tensile_strength

        # break bond
        if abs(force_n) >= tensile_strength*A_ij
            force_n = 0.
            force_t = 0.
            simulation.ice_floes[i].contacts[ic] = 0  # remove contact
            simulation.ice_floes[i].n_contacts -= 1
            simulation.ice_floes[j].n_contacts -= 1
        end
    end

    if k_t ≈ 0. && γ_t ≈ 0.
        # do nothing

    elseif k_t ≈ 0. && γ_t > 0.
        force_t = abs(γ_t * vel_t)
        if force_t > μ_d_minimum*abs(force_n)
            force_t = μ_d_minimum*abs(force_n)
        end
        if vel_t > 0.
            force_t = -force_t
        end

    elseif k_t > 0.

        force_t = -k_t*δ_t - γ_t*vel_t

        if abs(force_t) > μ_d_minimum*abs(force_n)
            force_t = μ_d_minimum*abs(force_n)*force_t/abs(force_t)
            δ_t = (-force_t - γ_t*vel_t)/k_t
        end

    else
        error("unknown contact_tangential_rheology (k_t = $k_t, γ_t = $γ_t")
    end

    simulation.ice_floes[i].contact_parallel_displacement[ic] = δ_t*t
    simulation.ice_floes[i].contact_age[ic] += simulation.time_step

    simulation.ice_floes[i].force += force_n*n + force_t*t;
    simulation.ice_floes[j].force -= force_n*n + force_t*t;

    simulation.ice_floes[i].torque += -force_t*R_ij
    simulation.ice_floes[j].torque += -force_t*R_ij

    simulation.ice_floes[i].pressure += 
        force_n/simulation.ice_floes[i].side_surface_area;
    simulation.ice_floes[j].pressure += 
        force_n/simulation.ice_floes[j].side_surface_area;
    nothing
end

function harmonicMean(a::Any, b::Any)
    if a ≈ 0. && b ≈ 0
        return 0.
    else
        return 2.*a*b/(a + b)
    end
end
