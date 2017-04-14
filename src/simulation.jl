## General simulation functions

function id(identifier::String)
    g_simulation_id = identifier
end

function id()
    return g_simulation_id
end


function run(verbose::Bool = true,
    status_interval = 100.,
    show_file_output = true)

    checkTimeParameters()
    if g_file_time_step > 0.0
        writeVTK(verbose = show_file_output)
    end

    time_since_output_file = 0.0

    while g_time <= g_time_total

        if g_file_time_step > 0.0 && time_since_output_file >= g_file_time_step
            writeVTK(verbose = show_file_output)
            time_since_output_file = 0.0
        end

        if verbose && g_time_iteration % status_interval == 0
            print("\r  t = $g_time/$g_time_total s            ")
        end

        findContacts()
        interact()
        updateKinematics()

        # Update time variables
        global g_time_iteration = g_time_iteration::Integer + 1
        incrementCurrentTime(g_time_step::float)
        time_since_output_file = time_since_output_file + g_time_step::float

    end

    if g_file_time_step > 0.0
        writeVTK(verbose = show_file_output)
    end
end
