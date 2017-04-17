## IO functions

function writeVTK(simulation::Simulation;
                  folder::String=".",
                  verbose::Bool=true)

    simulation.file_number += 1
    filename = string(folder, "/", simulation.id, ".", simulation.file_number, 
                      ".vtu")

    if verbose
        println("Output file: $filename")
    end

    f = open(filename, "w")

    write(f, "<?xml version=\"1.0\"?>\n" * # XML header
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " *
        "byte_order=\"LittleEndian\">\n" * # VTK header
        "  <UnstructuredGrid>\n" *
        "    <Piece NumberOfPoints=\"$(length(simulation.ice_floes))\" " *
        "NumberOfCells=\"0\">\n")

    # Coordinates for each point (positions)
    write(f, "      <Points>\n" *
        "        <DataArray name=\"Position [m]\" type=\"Float32\" " *
        "NumberOfComponents=\"3\" format=\"ascii\">\n" *
        "          ")

    for i=1:length(simulation.ice_floes)
        #=write(f, "%f %f %f " %
            g_position[i][1]::float,
            g_position[i][2]::float,
            g_position[i][3]::float)=#
            write(f, simulation.ice_floes[i].lin_pos[1], "\n",
                  simulation.ice_floes[i].lin_pos[2], "\n")
    end

    write(f, "\n" *
        "        </DataArray>\n" *
        "      </Points>\n")

    ### Data attributes
    write(f, "      <PointData Scalars=\"Diameter [m]\" " *
        "Vectors=\"vector\">\n")

    # Radii
    write(f, "        <DataArray type=\"Float32\" Name=\"Diameter (areal)\" " *
        "format=\"ascii\">\n" *
        "          ")
    for i=1:length(simulation.ice_floes)
        #write(f, "%f " % g_radius[i]::float*2.0)
        write(f, "$(simulation.ice_floes[i].areal_radius*2.) ")
    end
    write(f, "\n" *
        "        </DataArray>\n")

    # Footer
    write(f, "      </PointData>\n" *
        "      <Cells>\n" *
        "        <DataArray type=\"Int32\" Name=\"connectivity\" " *
        "format=\"ascii\">\n" *
        "        </DataArray>\n" *
        "        <DataArray type=\"Int32\" Name=\"offsets\" " *
        "format=\"ascii\">\n" *
        "        </DataArray>\n" *
        "        <DataArray type=\"UInt8\" Name=\"types\" " *
        "format=\"ascii\">\n" *
        "        </DataArray>\n" *
        "      </Cells>\n" *
        "    </Piece>\n" *
        "  </UnstructuredGrid>\n" *
        "</VTKFile>")

    close(f)

    println()
end
