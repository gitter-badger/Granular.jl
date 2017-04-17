## General simulation types and functions

# Simulation-scope data
type Simulation
    id::String

    time_iteration::Int
    time::Float64
    time_total::Float64
    time_step::Float64
    file_time_step::Float64   # 0.0: no output files
    file_number::Int

    gravitational_acceleration::vector

    ice_floes::Array{IceFloeCylindrical, 1}
    contact_pairs::Array{Integer, 1}
    wall_contacts::Array{Integer, 1}
end

function createSimulation(;id::String="unnamed",
                          time_iteration::Int=0,
                          time::Float64=0.0,
                          time_total::Float64=-1.,
                          time_step::Float64=-1.,
                          file_time_step::Float64=-1.,
                          file_number::Int=0,
                          gravitational_acceleration::vector=[0., 0.],
                          ice_floes=Array{IceFloeCylindrical, 1}[],
                          contact_pairs=Array{Integer, 1}[],
                          wall_contacts=Array{Integer, 1}[])

    return Simulation(id,
                      time_iteration,
                      time,
                      time_total,
                      time_step,
                      file_time_step,
                      file_number,
                      gravitational_acceleration,
                      ice_floes,
                      contact_pairs,
                      wall_contacts)
end

function id!(simulation::Simulation, identifier::String)
    simulation.id = identifier
end

function id(simulation::Simulation)
    return simulation.id
end


function run!(simulation::Simulation,
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

function addIceFloe!(simulation::Simulation, icefloe::IceFloeCylindrical)
    # Append icefloe to global icefloe array
    push!(simulation.ice_floes, icefloe)

    if verbose
        info("Added IceFloe $(length(g_ice_floes))")
    end
end

function removeIceFloe!(simulation::Simulation, i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    delete!(simulation.ice_floes, i)
end

