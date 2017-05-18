#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)
sim.ocean = SeaIce.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
SeaIce.writeVTK(sim, verbose=false)

cmd_post = ""
if Base.is_linux()
    cmd = "sha256sum"
elseif Base.is_apple()
    cmd = ["shasum", "-a", "256"]
elseif Base.is_windows()
    info("checksum verification not yet implemented on Windows")
    exit()
    cmd = ["powershell", "-Command", "\"Get-FileHash", "-Algorithm", "SHA256"]
    cmd_post = "\""
else
    error("checksum verification of VTK file not supported on this platform")
end

icefloechecksum = 
"88daceb1b99c519154b1acdcf8f340967794c552c74ea70c4af8954d8af5296a  " *
"test.icefloes.1.vtu\n"

oceanchecksum =
"d56ffb109841a803f2b2b94c74c87f7a497237204841d557d2b1043694d51f0d  " *
"test.ocean.1.vts\n"

@test readstring(`$(cmd) test.icefloes.1.vtu$(cmd_post)`) == icefloechecksum
@test readstring(`$(cmd) test.ocean.1.vts$(cmd_post)`) == oceanchecksum

SeaIce.removeSimulationFiles(sim)

info("Testing VTK write during run!()")
SeaIce.setOutputFileInterval!(sim, 1e-9)
SeaIce.setTotalTime!(sim, 1.5)
SeaIce.setTimeStep!(sim)
sim.file_number = 0
SeaIce.run!(sim, single_step=true)
SeaIce.run!(sim, single_step=true)

@test readstring(`$(cmd) test.icefloes.1.vtu$(cmd_post)`) == icefloechecksum
@test readstring(`$(cmd) test.ocean.1.vts$(cmd_post)`) == oceanchecksum

SeaIce.removeSimulationFiles(sim)
