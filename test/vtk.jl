#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)
sim.ocean = SeaIce.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
SeaIce.writeVTK(sim, verbose=false)

if Base.is_linux()
    cmd = "sha256sum"
elseif Base.is_apple()
    cmd = ["shasum", "-a", "256"]
else
    error("checksum verification of VTK file not supported on this platform")
end
@test readstring(`$(cmd) test.icefloes.1.vtu`) == "72f4e4b854d7e92afd8cde0b79a4af6a29e49714b751ffc30a4ff3867f44b50  test.icefloes.1.vtu\n"
@test readstring(`$(cmd) test.ocean.1.vts`) == "f0117e414c4e71a0c55980f63865eb03b6c597fa2546983258b8a57eb4ff2a25  test.ocean.1.vts\n"

rm("test.icefloes.1.vtu")
rm("test.ocean.1.vts")
