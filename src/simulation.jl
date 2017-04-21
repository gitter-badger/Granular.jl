## General simulation functions

"""
    createSimulation([id::String="unnamed",
                      time_iteration::Int=0,
                      time::Float64=0.0,
                      time_total::Float64=-1.,
                      time_step::Float64=-1.,
                      file_time_step::Float64=-1.,
                      file_number::Int=0,
                      ice_floes=Array{IceFloeCylindrical, 1}[],
                      contact_pairs=Array{Int64, 1}[],
                      overlaps=Array{Array{Float64, 1}, 1}[],

Create a simulation object containing all relevant variables such as temporal 
parameters, and lists of ice floes and contacts.

The parameter `id` is used to uniquely identify the simulation when it is 
written to disk.
"""
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
                          ocean::Ocean=createEmptyOcean())

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
                      ocean)
end

"""
    run!(simulation[,
         verbose::Bool = true,
         status_interval = 100.,
         show_file_output = true,
         single_step=false,
         temporal_integration_method="Three-term Taylor"])

Run the `simulation` through time until `simulation.time` equals or exceeds 
`simulatim.time_total`.  This function requires that all ice floes are added to 
the simulation and that the length of the computational time step is adjusted 
accordingly.

The function will search for contacts, determine the force balance on each ice 
floe, and integrate all kinematic degrees of freedom accordingly.  The temporal 
integration is explicit and of length `simulation.time_step`.  This function 
will write VTK files to disk in the intervals `simulation.file_time_step` by the 
function `writeVTK`.  If this value is negative, no output files will be written 
to disk.

# Arguments
* `simulation::Simulation`: the simulation to run (object is modified)
* `verbose::Bool=true`: show verbose information during the time loop
* `status_interval::Bool=true`: show verbose information during the time loop
* `show_file_output::Bool=true`: report to stdout when output file is written
* `single_step::Bool=false`: run simulation for a single time step only.  If 
    this causes `simulation.time` to exceed `simulation.time_total`, the latter 
    is increased accordingly.
* `temporal_integration_method::String="Three-term Taylor"`: type of integration 
    method to use.  See `updateIceFloeKinematics` for details.
"""
function run!(simulation::Simulation;
              verbose::Bool=true,
              status_interval::Int=100,
              show_file_output::Bool=true,
              single_step::Bool=false,
              temporal_integration_method::String="Three-term Taylor")

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
        incrementCurrentTime!(simulation, simulation.time_step)
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

"Add an `icefloe` to the `simulation` object.  If `verbose` is true, a short 
confirmation message will be printed to stdout`."
function addIceFloe!(simulation::Simulation,
                     icefloe::IceFloeCylindrical,
                     verbose::Bool = False)
    # Append icefloe to global icefloe array
    push!(simulation.ice_floes, icefloe)

    if verbose
        info("Added IceFloe $(length(simulation.ice_floes))")
    end
end

"Remove ice floe with index `i` from the `simulation` object."
function removeIceFloe!(simulation::Simulation, i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    delete!(simulation.ice_floes, i)
end

"Sets the `force` and `torque` values of all ice floes to zero."
function zeroForcesAndTorques!(simulation::Simulation)
    for icefloe in simulation.ice_floes
        icefloe.force = zeros(2)
        icefloe.torque = 0.
    end
end

"Prints the current simulation time and total time to standard out"
function reportSimulationTimeToStdout(simulation::Simulation)
    print("\r  t = ", simulation.time, '/', simulation.time_total,
          " s            ")
end
