#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical!(sim, [18., 0.], 10., 1., verbose=false)
sim.ocean = SeaIce.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
SeaIce.findContacts!(sim, method="all to all")
SeaIce.writeVTK(sim, verbose=false)

cmd_post = ""
if is_linux()
    cmd = "sha256sum"
elseif is_apple()
    cmd = ["shasum", "-a", "256"]
elseif is_windows()
    info("checksum verification not yet implemented on Windows")
    exit()
    cmd = ["powershell", "-Command", "\"Get-FileHash", "-Algorithm", "SHA256"]
    cmd_post = "\""
else
    error("checksum verification of VTK file not supported on this platform")
end

icefloepath = "test/test.icefloes.1.vtu"
icefloechecksum = 
"c75ffde29fbdd80161dafd524e690fbcbae2136d4f68c29f725d2d2454c6a162  " *
icefloepath * "\n"

icefloeinteractionpath = "test/test.icefloe-interaction.1.vtp"
icefloeinteractionchecksum = 
"881598f8f7279ece4301f6c94cb8f9146eb695f8c710edb446f49c1f7a061b84  " *
icefloeinteractionpath * "\n"

oceanpath = "test/test.ocean.1.vts"
oceanchecksum =
"d56ffb109841a803f2b2b94c74c87f7a497237204841d557d2b1043694d51f0d  " *
oceanpath * "\n"

@test readstring(`$(cmd) $(icefloepath)$(cmd_post)`) == icefloechecksum
@test readstring(`$(cmd) $(icefloeinteractionpath)$(cmd_post)`) == 
    icefloeinteractionchecksum
@test readstring(`$(cmd) $(oceanpath)$(cmd_post)`) == oceanchecksum

SeaIce.removeSimulationFiles(sim)

info("Testing VTK write during run!()")
SeaIce.setOutputFileInterval!(sim, 1e-9)
SeaIce.setTotalTime!(sim, 1.5)
SeaIce.setTimeStep!(sim)
sim.file_number = 0
SeaIce.run!(sim, single_step=true)
@test SeaIce.readSimulationStatus(sim.id) == 1
SeaIce.setOutputFileInterval!(sim, 0.1)
SeaIce.run!(sim)

SeaIce.status()

info("Testing generation of Paraview Python script")
SeaIce.writeParaviewPythonScript(sim,
                                 save_animation=true,
                                 save_images=false)
@test isfile("$(sim.id)/$(sim.id).py") && filesize("$(sim.id)/$(sim.id).py") > 0

info("Testing Paraview rendering if `pvpython` is present")
try
    cd(sim.id)
    run(`pvpython $(sim.id).py`)
    cd("..")
catch return_signal
    if !isa(return_signal, Base.UVError)
        @test isfile("$(sim.id)/$(sim.id).avi")
    end
end

SeaIce.writeParaviewPythonScript(sim,
                                 save_animation=false,
                                 save_images=true)
try
    cd(sim.id)
    run(`pvpython $(sim.id).py`)
    cd("..")
catch return_signal
    if !isa(return_signal, Base.UVError)
        @test isfile("$(sim.id)/$(sim.id).0000.png")
        @test isfile("$(sim.id)/$(sim.id).0014.png")
        SeaIce.render(sim)
        @test isfile("$(sim.id)/$(sim.id).0001.png")
    end
end

@test readstring(`$(cmd) $(icefloepath)$(cmd_post)`) == icefloechecksum
@test readstring(`$(cmd) $(icefloeinteractionpath)$(cmd_post)`) == 
    icefloeinteractionchecksum
@test readstring(`$(cmd) $(oceanpath)$(cmd_post)`) == oceanchecksum

SeaIce.removeSimulationFiles(sim)
