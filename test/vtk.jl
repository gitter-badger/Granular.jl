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
    cmd = ["powershell", "-Command", "\"Get-FileHash", "-Algorithm", "SHA256"]
    cmd_post = "\""
else
    error("checksum verification of VTK file not supported on this platform")
end

@test readstring(`$(cmd) test.icefloes.1.vtu$(cmd_post)`) == 
"88daceb1b99c519154b1acdcf8f340967794c552c74ea70c4af8954d8af5296a  " *
"test.icefloes.1.vtu\n"

@test readstring(`$(cmd) test.ocean.1.vts$(cmd_post)`) == 
"f0117e414c4e71a0c55980f63865eb03b6c597fa2546983258b8a57eb4ff2a25  " * 
"test.ocean.1.vts\n"

SeaIce.removeSimulationFiles(sim)
