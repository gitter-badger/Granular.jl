#!/usr/bin/env julia

# Check the contact search and geometry of a two-particle interaction

info("#### $(basename(@__FILE__)) ####")

info("Writing simple simulation to VTK file")
sim = SeaIce.createSimulation(id="test")
SeaIce.addIceFloeCylindrical(sim, [ 0., 0.], 10., 1., verbose=false)
SeaIce.addIceFloeCylindrical(sim, [18., 0.], 10., 1., verbose=false)
SeaIce.writeVTK(sim)

if Base.is_linux()
    checksum = readstring(`sha256sum test.1.vtu`)
elseif Base.is_apple()
    checksum = readstring(`shasum -a 256 test.1.vtu`)
else
    warn("checksum verification of VTK file not supported on this platform")
end
rm("test.1.vtu")
@test checksum == "1c0c2bdd265abdda22ef3727e7cac829e2321462d494be2e23364653f9529c87  test.1.vtu\n"
