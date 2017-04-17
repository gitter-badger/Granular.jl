## Set temporal parameters
function setTotalTime!(simulation::Simulation, t::float)
    if t <= 0.0
        error("Total time should be a positive value (t = $t s)")
    end
    simulation.time_total = t
end

function setCurrentTime!(simulation::Simulation, t::float)
    if t <= 0.0
        error("Current time should be a positive value (t = $t s)")
    end
    simulation.time = t
end

function incrementCurrentTime!(simulation::Simulation, t::float)
    if t <= 0.0
        error("Current time increment should be a positive value (t = $t s)")
    end
    simulation.time = simulation.time + t
end

function setOutputFileInterval!(simulation::Simulation, t::float)
    if t <= 0.0
        info("The output file interval should be a positive value (t = $t s)")
    end
    simulation.file_time_step = t
end

function disableOutputFiles!(simulation::Simulation)
    simulation.file_time_step = 0.0
end

function checkTimeParameters(simulation::Simulation)
    if simulation.time_total <= 0.0 || simulation.time_total <= simulation.time
        error("Total time should be positive and larger than current time.\n",
            "  t_total = ", simulation.time_total, " s, t = ", simulation.time, 
            " s")
    end
    if simulation.time_step <= 0.0
        error("Time step should be positive (t = ", simulation.time_step, ")")
    end
end

function findSmallestIceFloeMass(simulation::Simulation)
    m_min = 1e20
    i_min = -1
    for i in length(simulation.ice_floes)
        icefloe = simulation.ice_floes[i]
        if icefloe.mass < m_min
            m_min = icefloe.mass
            i_min = i
        end
    end
    return m_min, i_min
end

function findLargestIceFloeStiffness(simulation::Simulation)
    k_n_max = 0.
    k_t_max = 0.
    i_n_max = -1
    i_t_max = -1
    for i in length(simulation.ice_floes)
        icefloe = simulation.ice_floes[i]
        if icefloe.contact_stiffness_normal > k_n_max
            k_n_max = icefloe.contact_stiffness_normal
            i_n_max = i
        end
        if icefloe.contact_stiffness_tangential > k_t_max
            k_t_max = icefloe.contact_stiffness_tangential
            i_t_max = i
        end
    end
    return k_n_max, k_t_max, i_n_max, i_t_max
end

"""
Find the computational time step length suitable given the grain radii, contact
stiffnesses, and grain density. Uses the scheme by Radjaii et al. 2011.
"""
function setTimeStep!(simulation::Simulation;
                      epsilon::float=0.07, verbose::Bool=true)

    if length(simulation.ice_floes) < 1
        error("At least 1 grain is needed to find the computational time step.")
    end

    k_n_max, k_t_max, _, _ = findLargestIceFloeStiffness(simulation)
    m_min, _ = findSmallestIceFloeMass(simulation)

    simulation.time_step = epsilon/(sqrt(maximum([k_n_max, k_t_max])/m_min))

    if simulation.time_step <= 1.0e-20
        error("Time step too small or negative (", simulation.time_step, " s)")
    end

    if verbose
        info("Time step length t=",  simulation.time_step, " s")
    end
end

"""
Perform temporal integration for all grains.
"""
function updateKinematics!(simulation::Simulation)
    for i = 1:length(simulation.ice_floes)
        updateIceFloeKinematics(simulation.ice_floes[i])
    end
end
