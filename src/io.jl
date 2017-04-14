## IO functions

function writeVTK(folder::String = ".", verbose::Bool = true)

    global g_file_number += 1
    filename = "$(folder)/$(g_simulation_id).$(g_file_number).t=$(g_time)s.vtu"

    if verbose
        println("Output file: $filename")
    end

    f = open(filename, "w")

    write(f, "<?xml version=\"1.0\"?>\n" * # XML header
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " *
        "byte_order=\"LittleEndian\">\n" * # VTK header
        "  <UnstructuredGrid>\n" *
        "    <Piece NumberOfPoints=\"$(length(g_radius))\" " *
        "NumberOfCells=\"0\">\n")

    # Coordinates for each point (positions)
    write(f, "      <Points>\n" *
        "        <DataArray name=\"Position [m]\" type=\"Float32\" " *
        "NumberOfComponents=\"3\" format=\"ascii\">\n" *
        "          ")

    for i=1:length(g_radius)
        #=write(f, "%f %f %f " %
            g_position[i][1]::float,
            g_position[i][2]::float,
            g_position[i][3]::float)=#
        write(f, "$(g_position[i][1]::float)
            $(g_position[i][2]::float)
            $(g_position[i][3]::float) ")
    end

    write(f, "\n" *
        "        </DataArray>\n" *
        "      </Points>\n")

    ### Data attributes
    write(f, "      <PointData Scalars=\"Diameter [m]\" " *
        "Vectors=\"vector\">\n")

    # Radii
    write(f, "        <DataArray type=\"Float32\" Name=\"Diameter\" " *
        "format=\"ascii\">\n" *
        "          ")
    for i=1:length(g_radius)
        #write(f, "%f " % g_radius[i]::float*2.0)
        write(f, "$(g_radius[i]::float*2.0) ")
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
