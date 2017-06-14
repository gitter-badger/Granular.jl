## General simulation functions

export createSimulation
"""
    createSimulation([id::String="unnamed",
                      time_iteration::Int=0,
                      time::Float64=0.0,
                      time_total::Float64=-1.,
                      time_step::Float64=-1.,
                      file_time_step::Float64=-1.,
                      file_number::Int=0,
                      ice_floes=Array{IceFloeCylindrical, 1}[],
                      ocean::Ocean,
                      atmosphere::Atmosphere)

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
                          file_time_since_output_file::Float64=0.,
                          ice_floes=Array{IceFloeCylindrical, 1}[],
                          ocean::Ocean=createEmptyOcean(),
                          atmosphere::Atmosphere=createEmptyAtmosphere())

    return Simulation(id,
                      time_iteration,
                      time,
                      time_total,
                      time_step,
                      file_time_step,
                      file_number,
                      file_time_since_output_file,
                      ice_floes,
                      ocean,
                      atmosphere)
end

export run!
"""
    run!(simulation[,
         verbose::Bool = true,
         status_interval = 100.,
         show_file_output = true,
         single_step = false,
         temporal_integration_method = "Three-term Taylor"],
         write_jld = false)

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
* `write_jld::Bool=false`: write simulation state to disk as JLD files (see 
    `SeaIce.writeSimulation(...)` whenever saving VTK output.
"""
function run!(simulation::Simulation;
              verbose::Bool=true,
              status_interval::Int=100,
              show_file_output::Bool=true,
              single_step::Bool=false,
              temporal_integration_method::String="Three-term Taylor",
              write_jld::Bool=false)

    if single_step && simulation.time >= simulation.time_total
        simulation.time_total += simulation.time_step
    end

    checkTimeParameters(simulation)

    # number of time steps between output files
    n_file_time_step = Int(ceil(simulation.file_time_step/simulation.time_step))

    while simulation.time <= simulation.time_total

        if simulation.file_time_step > 0.0 &&
            simulation.time_iteration % n_file_time_step == 0

            if show_file_output
                println()
            end
            if write_jld
                writeSimulation(simulation, verbose=show_file_output)
            end
            writeVTK(simulation, verbose=show_file_output)
            writeSimulationStatus(simulation, verbose=show_file_output)
            simulation.file_time_since_output_file = 0.0
        end

        if verbose && simulation.time_iteration % status_interval == 0
            reportSimulationTimeToStdout(simulation)
        end

        zeroForcesAndTorques!(simulation)

        if typeof(simulation.atmosphere.input_file) != Bool
            sortIceFloesInGrid!(simulation, simulation.atmosphere)
        end

        if typeof(simulation.ocean.input_file) != Bool
            sortIceFloesInGrid!(simulation, simulation.ocean)
            findContacts!(simulation, method="ocean grid")

        elseif typeof(simulation.atmosphere.input_file) != Bool
            findContacts!(simulation, method="atmosphere grid")

        else
            findContacts!(simulation, method="all to all")
        end

        interact!(simulation)

        if typeof(simulation.ocean.input_file) != Bool
            addOceanDrag!(simulation)
        end

        if typeof(simulation.atmosphere.input_file) != Bool
            addAtmosphereDrag!(simulation)
        end

        updateIceFloeKinematics!(simulation, method=temporal_integration_method)

        # Update time variables
        simulation.time_iteration += 1
        incrementCurrentTime!(simulation, simulation.time_step)

        if single_step
            return
        end
    end
    if verbose
        reportSimulationTimeToStdout(simulation)
        println()
    end
end

export addIceFloe!
"""
    addIceFloe!(simulation::Simulation,
                icefloe::IceFloeCylindrical,
                verbose::Bool = False)

Add an `icefloe` to the `simulation` object.  If `verbose` is true, a short 
confirmation message will be printed to stdout.
"""
function addIceFloe!(simulation::Simulation,
                     icefloe::IceFloeCylindrical,
                     verbose::Bool = False)
    push!(simulation.ice_floes, icefloe)

    if verbose
        info("Added IceFloe $(length(simulation.ice_floes))")
    end
end

export disableIceFloe!
"Disable ice floe with index `i` in the `simulation` object."
function disableIceFloe!(simulation::Simulation, i::Integer)
    if i < 1
        error("Index must be greater than 0 (i = $i)")
    end

    simulation.ice_floes[i].enabled = false
end

export zeroForcesAndTorques!
"Sets the `force` and `torque` values of all ice floes to zero."
function zeroForcesAndTorques!(simulation::Simulation)
    for icefloe in simulation.ice_floes
        icefloe.force = zeros(2)
        icefloe.torque = 0.
        icefloe.pressure = 0.
    end
end

export reportSimulationTimeToStdout
"Prints the current simulation time and total time to standard out"
function reportSimulationTimeToStdout(simulation::Simulation)
    print("\r  t = ", simulation.time, '/', simulation.time_total,
          " s            ")
end

export compareSimulations
"""
    compareSimulations(sim1::Simulation, sim2::Simulation)

Compare values of two `Simulation` objects using the `Base.Test` framework.
"""
function compareSimulations(sim1::Simulation, sim2::Simulation)

    Base.Test.@test sim1.id == sim2.id

    Base.Test.@test sim1.time_iteration == sim2.time_iteration
    Base.Test.@test sim1.time ≈ sim2.time
    Base.Test.@test sim1.time_total ≈ sim2.time_total
    Base.Test.@test sim1.time_step ≈ sim2.time_step
    Base.Test.@test sim1.file_time_step ≈ sim2.file_time_step
    Base.Test.@test sim1.file_number == sim2.file_number
    Base.Test.@test sim1.file_time_since_output_file ≈ 
        sim2.file_time_since_output_file

    for i=1:length(sim1.ice_floes)
        compareIceFloes(sim1.ice_floes[i], sim2.ice_floes[i])
    end
    compareOceans(sim1.ocean, sim2.ocean)
    compareAtmospheres(sim1.atmosphere, sim2.atmosphere)
end
