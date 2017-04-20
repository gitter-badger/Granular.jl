## General simulation functions

function createSimulation(;id::String="unnamed",
                          time_iteration::Int=0,
                          time::Float64=0.0,
                          time_total::Float64=-1.,
                          time_step::Float64=-1.,
                          file_time_step::Float64=-1.,
                          file_number::Int=0,
                          ice_floes=Array{IceFloeCylindrical, 1}[],
                          contact_pairs=Array{Int64, 1}[],
                          overlaps=Array{Array{Float64, 1}, 1}[],
                          wall_contacts=Array{Int64, 1}[])

    return Simulation(id,
                      time_iteration,
                      time,
                      time_total,
                      time_step,
                      file_time_step,
                      file_number,
                      ice_floes,
                      contact_pairs,
                      overlaps,
                      wall_contacts)
end

function run!(simulation::Simulation;
              verbose::Bool = true,
              status_interval = 100.,
              show_file_output = true,
              single_step=false,
              temporal_integration_method="Three-term Taylor")

    if single_step && simulation.time >= simulation.time_total
        simulation.time_total += simulation.time_step
    end

    checkTimeParameters(simulation)
    if simulation.file_time_step > 0.0
        writeVTK(simulation, verbose=show_file_output)
    end

    time_since_output_file = 0.0

    while simulation.time <= simulation.time_total

        if simulation.file_time_step > 0.0 &&
            time_since_output_file >= simulation.file_time_step

            writeVTK(simulation, verbose=show_file_output)
            time_since_output_file = 0.0
        end

        if verbose && simulation.time_iteration % status_interval == 0
            reportSimulationTimeToStdout(simulation)
        end

        zeroForcesAndTorques!(simulation)
        findContacts!(simulation)
        interact!(simulation)
        updateIceFloeKinematics!(simulation, method=temporal_integration_method)

        # Update time variables
        simulation.time_iteration += 1
        simulation.time += simulation.time_step
        time_since_output_file = time_since_output_file + simulation.time_step

        if single_step
            if verbose
                println("Current time: ", simulation.time)
            end
            return
        end
    end
    if verbose
        reportSimulationTimeToStdout(simulation)
        println()
    end
end

function addIceFloe!(simulation::Simulation,
                     icefloe::IceFloeCylindrical,
                     verbose::Bool = False)
    # Append icefloe to global icefloe array
    push!(simulation.ice_floes, icefloe)

    if verbose
        info("Added IceFloe $(length(simulation.ice_floes))")
    end
end

function removeIceFloe!(simulation::Simulation, i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    delete!(simulation.ice_floes, i)
end

function zeroForcesAndTorques!(simulation::Simulation)
    for icefloe in simulation.ice_floes
        icefloe.force = zeros(2)
        icefloe.torque = 0.
    end
end

function reportSimulationTimeToStdout(simulation::Simulation)
    print("\r  t = ", simulation.time, '/', simulation.time_total,
          " s            ")
end
