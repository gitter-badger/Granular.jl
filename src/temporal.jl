## Set temporal parameters
function setTotalTime(t::float)
    if t <= 0.0
        error("Total time should be a positive value (t = $t s)")
    end
    global g_time_total = t
end

function setCurrentTime(t::float)
    if t <= 0.0
        error("Current time should be a positive value (t = $t s)")
    end
    global g_time = t
end

function incrementCurrentTime(t::float)
    if t <= 0.0
        error("Current time increment should be a positive value (t = $t s)")
    end
    global g_time = g_time::float + t
end

function setOutputFileInterval(t::float)
    if t <= 0.0
        info("The output file interval should be a positive value (t = $t s)")
    end
    global g_file_time_step = t
end

function disableOutputFiles()
    global g_file_time_step = 0.0
end

function checkTimeParameters()
    if g_time_total <= 0.0 || g_time_total <= g_time
        error("Total time should be positive and larger than current time. " *
            "t_total = $g_time_total s, t = $g_time s")
    end
    if g_time_step <= 0.0
        error("Time step should be positive (Δt = $g_time_step)")
    end
end


"""
Find the computational time step length suitable given the grain radii, contact
stiffnesses, and grain density. Uses the scheme by Radjaii et al. 2011.
"""
function findTimeStep(epsilon::float=0.07, verbose::Bool=true)

    if length(g_radius) < 1
        error("At least 1 grain is needed to find the computational time step.")
    end

    global g_time_step = epsilon / (sqrt(
        maximum([g_contact_stiffness_normal::float,
        g_contact_stiffness_tangential::float]) /
        minimum(g_mass)))

    if g_time_step <= 1.0e-20
        error("Time step too small or negative ($g_time_step s)")
    end

    if verbose
        info("Time step length Δt = $g_time_step s")
    end
end


"""
Perform temporal integration for all grains.
"""
function updateKinematics()
    for i = 1:length(g_radius)
        updateIceFloeKinematics(i)
    end
end

