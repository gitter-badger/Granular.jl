## General simulation types and functions

# Simulation-scope data
type simulation

    time_iteration::Int
    time::Float64
    time_total::Float64
    time_step::Float64
    file_time_step::Float64   # 0.0: no output files
    file_number::Int

    gravitational_acceleration::vector

    id::String

    ice_floes::Array{IceFloeCylindrical, 1}
    contact_pairs::Array{Integer, 1}
    wall_contacts::Array{Integer, 1}

    # default values
    simulation(time_iteration=0.0,
              time=0.0,
              origo=[0., 0.],
              file_number=0,
              ice_floes=Array{IceFloeCylindrical, 1}[],
              contact_pairs=Array{Integer, 1}[],
              wall_contacts=Array{Integer, 1}[]) = new(simulation)
end

function id(simulation::simulation, identifier::String)
    simulation.id = identifier
end

function id(simulation)
    return simulation.id
end


function run(simulation::simulation,
             verbose::Bool = true,
             status_interval = 100.,
             show_file_output = true)

    checkTimeParameters(simulation)
    if simulation.file_time_step > 0.0
        writeVTK(simulation, verbose = show_file_output)
    end

    time_since_output_file = 0.0

    while simulation.time <= simulation.time_total

        if simulation.file_time_step > 0.0 &&
            simulation.time_since_output_file >= simulation.file_time_step

            writeVTK(simulation=simulation, verbose=show_file_output)
            time_since_output_file = 0.0
        end

        if verbose && simulation.time_iteration % status_interval == 0
            print("\r  t = ", simulation.time, '/', simulation.time_total,
                  " s            ")
        end

        findContacts()
        interact()
        updateKinematics()

        # Update time variables
        simulation.time_iteration += 1
        simulation.time += simulation.time_step
        time_since_output_file = time_since_output_file + simulation.time_step
    end

    if simulation.file_time_step > 0.0
        writeVTK(simulation, verbose=show_file_output)
    end
end
