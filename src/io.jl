import WriteVTK
import NetCDF
hasJLD = false
if typeof(Pkg.installed("JLD")) == VersionNumber
    import JLD
    hasJLD = true
end

## IO functions

export writeSimulation
"""
    writeSimulation(simulation::Simulation;
                         filename::String="",
                         folder::String=".",
                         verbose::Bool=true)

Write all content from `Simulation` to disk in JDL format.  If the `filename` 
parameter is not specified, it will be saved to a subdirectory under the current 
directory named after the simulation identifier `simulation.id`.
"""
function writeSimulation(simulation::Simulation;
                         filename::String="",
                         folder::String=".",
                         verbose::Bool=true)
    if !hasJLD
        warn("Package JLD not found. Simulation save/read not supported. " * 
             "Please install JLD and its " *
             "requirements with `Pkg.add(\"JLD\")`.")
    else
        if filename == ""
            folder = folder * "/" * simulation.id
            mkpath(folder)
            filename = string(folder, "/", simulation.id, ".",
                              simulation.file_number, ".jld")
        end

        JLD.save(filename, "simulation", simulation)

        if verbose
            info("simulation written to $filename")
        end
    end
    nothing
end

export readSimulation
"""
    readSimulation(filename::String="";
                   verbose::Bool=true)

Read all content from `Simulation` from disk in JDL format.
"""
function readSimulation(filename::String="";
                         verbose::Bool=true)
    if !hasJLD
        warn("Package JLD not found. Simulation save/read not supported. " * 
             "Please install JLD and its " *
             "requirements with `Pkg.add(\"JLD\")`.")
        nothing
    else
        if verbose
            info("reading simulation from $filename")
        end
        return JLD.load(filename, "simulation")
    end
end

export writeSimulationStatus
"""
    writeSimulationStatus(simulation::Simulation;
                          folder::String=".",
                          verbose::Bool=false)

Write current simulation status to disk in a minimal txt file.
"""
function writeSimulationStatus(simulation::Simulation;
                               folder::String=".",
                               verbose::Bool=false)
    folder = folder * "/" * simulation.id
    mkpath(folder)
    filename = string(folder, "/", simulation.id, ".status.txt")

    writedlm(filename, [simulation.time
                        simulation.time/simulation.time_total*100.
                        float(simulation.file_number)])
    if verbose
        info("wrote status to $filename")
    end
    nothing
end

export readSimulationStatus
"""
    readSimulationStatus(filename::String;
                         folder::String=".",
                         verbose::Bool=false)

Write current simulation status to disk in a minimal txt file.
"""
function readSimulationStatus(simulation_id::String;
                              folder::String=".",
                              verbose::Bool=true)

    folder = folder * "/" * simulation_id
    filename = string(folder, "/", simulation_id, ".status.txt")

    data = readdlm(filename)
    if verbose
        info("$simulation_id:\n" *
             "  time:             $(data[1]) s\n" *
             "  complete:         $(data[2])%\n" *
             "  last output file: $(Int(round(data[3])))\n")
    end
    return data[3]
end

export status
"""
Shows the status of all simulations with output files written under the 
specified `folder`, which is the current working directory by default.
"""
function status(folder::String=".";
                loop::Bool=false,
                t_int::Int=10,
                colored_output::Bool=true,
                write_header::Bool=true)

    if colored_output
        id_color_complete = :green
        id_color_in_progress = :yellow
        time_color = :white
        percentage_color = :blue
        lastfile_color = :cyan
    else
        id_color_complete = :default
        id_color_in_progress = :default
        time_color = :default
        percentage_color = :default
        lastfile_color = :default
    end

    repeat = true
    while repeat

        status_files = String[]
        println(Dates.format(DateTime(now()), Dates.RFC1123Format))

        for (root, dirs, files) in walkdir(folder, follow_symlinks=false)

            for file in files
                if contains(file, ".status.txt")
                    push!(status_files, joinpath(root, file))
                end
            end
        end

        if length(status_files) > 0
            if write_header
                println("--------------------------------------" * 
                        "--------------------------------------")
                print_with_color(:default, "simulation folder \t")
                print_with_color(time_color, "      time \t")
                print_with_color(percentage_color, "      completed  ")
                print_with_color(lastfile_color, "last file \n")
                println("--------------------------------------" * 
                        "--------------------------------------")
            end

            for file in status_files
                data = readdlm(file)
                id = replace(file, ".status.txt", "")
                id = replace(id, "./", "")
                id = replace(id, r".*/", "")
                time_s = @sprintf "%6.2fs" data[1]
                time_h = @sprintf "%5.1fh" data[1]/(60.*60.)
                percentage = @sprintf "%3.0f%%" data[2]
                lastfile = @sprintf "%5d" data[3]
                if data[2] < 99.
                    print_with_color(id_color_in_progress, "$id \t")
                else
                    print_with_color(id_color_complete, "$id \t")
                end
                print_with_color(time_color, "$time_s ($time_h) \t")
                print_with_color(percentage_color, "$percentage \t")
                print_with_color(lastfile_color, "$lastfile \n")
            end
            if write_header
                println("--------------------------------------" * 
                        "--------------------------------------")
            end
        else
            warn("no simulations found in $(pwd())/$folder")
        end

        if loop && t_int > 0
            sleep(t_int)
        end
        if !loop
            repeat = false
        end
    end
    nothing
end

export writeVTK
"""
Write a VTK file to disk containing all ice floes in the `simulation` in an 
unstructured mesh (file type `.vtu`).  These files can be read by ParaView and 
can be visualized by applying a *Glyph* filter.

If the simulation contains an `Ocean` data structure, it's contents will be 
written to separate `.vtu` files.  This can be disabled by setting the argument 
`ocean=false`.  The same is true for the atmosphere.

The VTK files will be saved in a subfolder named after the simulation.
"""
function writeVTK(simulation::Simulation;
                  folder::String=".",
                  verbose::Bool=true,
                  ocean::Bool=true,
                  atmosphere::Bool=true)

    simulation.file_number += 1
    folder = folder * "/" * simulation.id
    mkpath(folder)

    filename = string(folder, "/", simulation.id, ".icefloes.", 
                      simulation.file_number)
    writeIceFloeVTK(simulation, filename, verbose=verbose)

    filename = string(folder, "/", simulation.id, ".icefloe-interaction.", 
                      simulation.file_number)
    writeIceFloeInteractionVTK(simulation, filename, verbose=verbose)

    if typeof(simulation.ocean.input_file) != Bool && ocean
        filename = string(folder, "/", simulation.id, ".ocean.", 
                        simulation.file_number)
        writeGridVTK(simulation.ocean, filename, verbose=verbose)
    end

    if typeof(simulation.atmosphere.input_file) != Bool && atmosphere
        filename = string(folder, "/", simulation.id, ".atmosphere.", 
                        simulation.file_number)
        writeGridVTK(simulation.atmosphere, filename, verbose=verbose)
    end
    nothing
end

export writeIceFloeVTK
"""
Write a VTK file to disk containing all ice floes in the `simulation` in an 
unstructured mesh (file type `.vtu`).  These files can be read by ParaView and 
can be visualized by applying a *Glyph* filter.  This function is called by 
`writeVTK()`.
"""
function writeIceFloeVTK(simulation::Simulation,
                         filename::String;
                         verbose::Bool=false)

    ifarr = convertIceFloeDataToArrays(simulation)
    
    # add arrays to VTK file
    vtkfile = WriteVTK.vtk_grid(filename, ifarr.lin_pos, WriteVTK.MeshCell[])

    WriteVTK.vtk_point_data(vtkfile, ifarr.density, "Density [kg m^-3]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.thickness, "Thickness [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_radius*2.,
                            "Diameter (contact) [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.areal_radius*2.,
                            "Diameter (areal) [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.circumreference,
                            "Circumreference  [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.horizontal_surface_area,
                            "Horizontal surface area [m^2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.side_surface_area,
                            "Side surface area [m^2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.volume, "Volume [m^3]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.mass, "Mass [kg]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.moment_of_inertia,
                            "Moment of inertia [kg m^2]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.lin_vel, "Linear velocity [m s^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.lin_acc,
                            "Linear acceleration [m s^-2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.force, "Sum of forces [N]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_pos, "Angular position [rad]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_vel,
                            "Angular velocity [rad s^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_acc,
                            "Angular acceleration [rad s^-2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.torque, "Sum of torques [N*m]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.fixed, "Fixed in space [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.rotating, "Free to rotate [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.enabled, "Enabled [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_stiffness_normal,
                            "Contact stiffness (normal) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_stiffness_tangential,
                            "Contact stiffness (tangential) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_viscosity_normal,
                            "Contact viscosity (normal) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_viscosity_tangential,
                            "Contact viscosity (tangential) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_static_friction,
                            "Contact friction (static) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_dynamic_friction,
                            "Contact friction (dynamic) [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.youngs_modulus,
                            "Young's modulus [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.poissons_ratio,
                            "Poisson's ratio [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.tensile_strength,
                            "Tensile strength [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.compressive_strength_prefactor,
                            "Compressive strength prefactor [m^0.5 Pa]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_drag_coeff_vert,
                            "Ocean drag coefficient (vertical) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_drag_coeff_horiz,
                            "Ocean drag coefficient (horizontal) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_drag_coeff_vert,
                            "Atmosphere drag coefficient (vertical) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_drag_coeff_horiz,
                            "Atmosphere drag coefficient (horizontal) [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.pressure,
                            "Contact pressure [Pa]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.n_contacts,
                            "Number of contacts [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.granular_stress,
                            "Granular stress [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_stress,
                            "Ocean stress [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_stress,
                            "Atmosphere stress [Pa]")

    deleteIceFloeArrays!(ifarr)
    ifarr = 0
    gc()

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    end
    nothing
end

export writeIceFloeInteractionVTK
"""
    writeIceFloeInteractionVTK(simulation::Simulation,
                               filename::String;
                               verbose::Bool=false)

Saves ice-floe interactions to `.vtp` files for visualization with VTK, for 
example in Paraview.  Convert Cell Data to Point Data and use with Tube filter.
"""
function writeIceFloeInteractionVTK(simulation::Simulation,
                                    filename::String;
                                    verbose::Bool=false)

    i1 = Int64[]
    i2 = Int64[]
    inter_particle_vector = Vector{Float64}[]
    force = Float64[]
    effective_radius = Float64[]
    contact_area = Float64[]
    contact_stiffness = Float64[]
    tensile_stress = Float64[]
    shear_displacement = Vector{Float64}[]
    contact_age = Float64[]
    for i=1:length(simulation.ice_floes)
        for ic=1:simulation.Nc_max
            if simulation.ice_floes[i].contacts[ic] > 0
                j = simulation.ice_floes[i].contacts[ic]

                if !simulation.ice_floes[i].enabled ||
                    !simulation.ice_floes[j].enabled
                    continue
                end

                p = simulation.ice_floes[i].lin_pos -
                    simulation.ice_floes[j].lin_pos
                dist = norm(p)

                r_i = simulation.ice_floes[i].contact_radius
                r_j = simulation.ice_floes[j].contact_radius
                δ_n = dist - (r_i + r_j)
                R_ij = harmonicMean(r_i, r_j)

                if simulation.ice_floes[i].youngs_modulus > 0. &&
                    simulation.ice_floes[j].youngs_modulus > 0.
                    E_ij = harmonicMean(simulation.ice_floes[i].
                                        youngs_modulus,
                                        simulation.ice_floes[j].
                                        youngs_modulus)
                    A_ij = R_ij*min(simulation.ice_floes[i].thickness, 
                                    simulation.ice_floes[j].thickness)
                    k_n = E_ij*A_ij/R_ij
                else
                    k_n = harmonicMean(simulation.ice_floes[i].
                                       contact_stiffness_normal,
                                       simulation.ice_floes[j].
                                       contact_stiffness_normal)
                end

                
                push!(i1, i)
                push!(i2, j)
                push!(inter_particle_vector, p)

                push!(force, k_n*δ_n)
                push!(effective_radius, R_ij)
                push!(contact_area, A_ij)
                push!(contact_stiffness, k_n)
                push!(tensile_stress, k_n*δ_n/A_ij)

                push!(shear_displacement, simulation.ice_floes[i].
                      contact_parallel_displacement[ic])

                push!(contact_age, simulation.ice_floes[i].contact_age[ic])
            end
        end
    end

    # Insert a piece for each ice floe interaction using ice floe positions as 
    # coordinates and connect them with lines by referencing their indexes.
    open(filename * ".vtp", "w") do f
        write(f, "<?xml version=\"1.0\"?>\n")
        write(f, "<VTKFile type=\"PolyData\" version=\"0.1\" " *
              "byte_order=\"LittleEndian\">\n")
        write(f, "  <PolyData>\n")
        write(f, "    <Piece " *
              "NumberOfPoints=\"$(length(simulation.ice_floes))\" " *
              "NumberOfVerts=\"0\" " *
              "NumberOfLines=\"$(length(i1))\" " *
              "NumberOfStrips=\"0\" " *
              "NumberOfPolys=\"0\">\n")
        write(f, "      <PointData>\n")
        write(f, "      </PointData>\n")
        write(f, "      <CellData>\n")

        # Write values associated to each line
        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Inter-particle vector [m]\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(inter_particle_vector[i][1]) ")
            write(f, "$(inter_particle_vector[i][2]) ")
            write(f, "0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Shear displacement [m]\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(shear_displacement[i][1]) ")
            @inbounds write(f, "$(shear_displacement[i][2]) ")
            write(f, "0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" Name=\"Force [N]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(force[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Effective radius [m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(effective_radius[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact area [m^2]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_area[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact stiffness [N/m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_stiffness[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Tensile stress [Pa]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(tensile_stress[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact age [s]\" NumberOfComponents=\"1\" 
        format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_age[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </CellData>\n")
        write(f, "      <Points>\n")

        # Write line endpoints (ice floe centers)
        #write(f, "        <DataArray Name=\"Position [m]\" type=\"Float32\" " *
        write(f, "        <DataArray type=\"Float32\" Name=\"Points\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i in simulation.ice_floes
            @inbounds write(f, "$(i.lin_pos[1]) $(i.lin_pos[2]) 0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Points>\n")
        write(f, "      <Verts>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Verts>\n")
        write(f, "      <Lines>\n")

        # Write contact connectivity by referring to point indexes
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(i1[i] - 1) $(i2[i] - 1) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        
        # Write 0-indexed offset for the connectivity array for the end of each 
        # cell
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$((i - 1)*2 + 2) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </Lines>\n")
        write(f, "      <Strips>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Strips>\n")
        write(f, "      <Polys>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Polys>\n")
        write(f, "    </Piece>\n")
        write(f, "  </PolyData>\n")
        write(f, "</VTKFile>\n")
    end
    nothing
end

export writeOceanVTK
"""
Write a VTK file to disk containing all ocean data in the `simulation` in a 
structured grid (file type `.vts`).  These files can be read by ParaView and can 
be visualized by applying a *Glyph* filter.  This function is called by 
`writeVTK()`.
"""
function writeGridVTK(grid::Any,
                      filename::String;
                      verbose::Bool=false)
    
    # make each coordinate array three-dimensional
    xq = similar(grid.u[:,:,:,1])
    yq = similar(grid.u[:,:,:,1])
    zq = similar(grid.u[:,:,:,1])

    for iz=1:size(xq, 3)
        @inbounds xq[:,:,iz] = grid.xq
        @inbounds yq[:,:,iz] = grid.yq
    end
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            @inbounds zq[ix,iy,:] = grid.zl
        end
    end

    # add arrays to VTK file
    vtkfile = WriteVTK.vtk_grid(filename, xq, yq, zq)

    WriteVTK.vtk_point_data(vtkfile, grid.u[:, :, :, 1],
                            "u: Zonal velocity [m/s]")
    WriteVTK.vtk_point_data(vtkfile, grid.v[:, :, :, 1],
                            "v: Meridional velocity [m/s]")
    # write velocities as 3d vector
    vel = zeros(3, size(xq, 1), size(xq, 2), size(xq, 3))
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            for iz=1:size(xq, 3)
                @inbounds vel[1, ix, iy, iz] = grid.u[ix, iy, iz, 1]
                @inbounds vel[2, ix, iy, iz] = grid.v[ix, iy, iz, 1]
            end
        end
    end
    
    WriteVTK.vtk_point_data(vtkfile, vel, "Velocity vector [m/s]")

    if typeof(grid) == Ocean
        WriteVTK.vtk_point_data(vtkfile, grid.h[:, :, :, 1],
                                "h: Layer thickness [m]")
        WriteVTK.vtk_point_data(vtkfile, grid.e[:, :, :, 1],
                                "e: Relative interface height [m]")
    end

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    end
    nothing
end

export writeParaviewPythonScript
"""
function writeParaviewPythonScript(simulation,
                                   [filename, folder, vtk_folder, verbose])

Create a `".py"` script for visualizing the simulation VTK files in Paraview.
The script can be run from the command line with `pvpython` (bundled with
Paraview), or from the interactive Python shell inside Paraview.

# Arguments
* `simulation::Simulation`: input simulation file containing the data.
* `filename::String`: output file name for the Python script. At its default
    (blank) value, the script is named after the simulation id (`simulation.id`).
* `folder::String`: output directory, current directory the default.
* `vtk_folder::String`: directory containing the VTK output files, by default
    `"./<simulation.id>/"`.
* `save_animation::Bool`: make the generated script immediately save a rendered
    animation to disk when the `".py"` script is called.
* `verbose::Bool`: show diagnostic information during
function call, on by
    default.
"""
function writeParaviewPythonScript(simulation::Simulation;
                                   filename::String="",
                                   folder::String=".",
                                   vtk_folder::String="",
                                   save_animation::Bool=true,
                                   save_images::Bool=false,
                                   width::Integer=1920,
                                   height::Integer=1080,
                                   framerate::Integer=10,
                                   ice_floes_color_scheme::String="X Ray",
                                   verbose::Bool=true)
    if filename == ""
        mkpath(folder)
        filename = string(folder, "/", simulation.id, ".py")
    end
    if vtk_folder == ""
        vtk_folder = "."
    end

    open(filename, "w") do f
        write(f, """from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
FileName=[""")
        for i=1:simulation.file_number
            write(f, "'$(vtk_folder)/$(simulation.id).icefloes.$(i).vtu', ")
        end
        write(f, """]
imageicefloes = XMLUnstructuredGridReader(FileName=FileName)

imageicefloes.PointArrayStatus = [
'Density [kg m^-3]',
'Thickness [m]',
'Diameter (contact) [m]',
'Diameter (areal) [m]',
'Circumreference  [m]',
'Horizontal surface area [m^2]',
'Side surface area [m^2]',
'Volume [m^3]',
'Mass [kg]',
'Moment of inertia [kg m^2]',
'Linear velocity [m s^-1]',
'Linear acceleration [m s^-2]',
'Sum of forces [N]',
'Angular position [rad]',
'Angular velocity [rad s^-1]',
'Angular acceleration [rad s^-2]',
'Sum of torques [N*m]',
'Fixed in space [-]',
'Free to rotate [-]',
'Enabled [-]',
'Contact stiffness (normal) [N m^-1]',
'Contact stiffness (tangential) [N m^-1]',
'Contact viscosity (normal) [N m^-1 s]',
'Contact viscosity (tangential) [N m^-1 s]',
'Contact friction (static) [-]',
'Contact friction (dynamic) [-]',
"Young's modulus [Pa]",
"Poisson's ratio [-]",
'Tensile strength [Pa]'
'Compressive strength prefactor [m^0.5 Pa]',
'Ocean drag coefficient (vertical) [-]',
'Ocean drag coefficient (horizontal) [-]',
'Atmosphere drag coefficient (vertical) [-]',
'Atmosphere drag coefficient (horizontal) [-]',
'Contact pressure [Pa]',
'Number of contacts [-]',
'Granular stress [Pa]',
'Ocean stress [Pa]',
'Atmosphere stress [Pa]']

animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2478, 1570]

# show data in view
imageicefloesDisplay = Show(imageicefloes, renderView1)
# trace defaults for the display properties.
imageicefloesDisplay.Representation = 'Surface'
imageicefloesDisplay.AmbientColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.ColorArrayName = [None, '']
imageicefloesDisplay.OSPRayScaleArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
imageicefloesDisplay.SelectOrientationVectors = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.ScaleFactor = 6.050000000000001
imageicefloesDisplay.SelectScaleArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.GlyphType = 'Arrow'
imageicefloesDisplay.GlyphTableIndexArray = 'Angular acceleration [rad s^-2]'
imageicefloesDisplay.DataAxesGrid = 'GridAxesRepresentation'
imageicefloesDisplay.PolarAxes = 'PolarAxesRepresentation'
imageicefloesDisplay.ScalarOpacityUnitDistance = 64.20669746996803
imageicefloesDisplay.GaussianRadius = 3.0250000000000004
imageicefloesDisplay.SetScaleArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imageicefloesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
imageicefloesDisplay.OpacityArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imageicefloesDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
imageicefloesDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
imageicefloesDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
imageicefloesDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(Input=imageicefloes,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
glyph1.Vectors = ['POINTS', 'Angular acceleration [rad s^-2]']
glyph1.ScaleFactor = 6.050000000000001
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.Scalars = ['POINTS', 'Diameter (areal) [m]']
glyph1.Vectors = ['POINTS', 'Angular position [rad]']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0
glyph1.GlyphMode = 'All Points'

# get color transfer function/color map for 'Diameterarealm'
diameterarealmLUT = GetColorTransferFunction('Diameterarealm')

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.LookupTable = diameterarealmLUT
glyph1Display.OSPRayScaleArray = 'Diameter (areal) [m]'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 6.1000000000000005
glyph1Display.SelectScaleArray = 'Diameter (areal) [m]'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Diameter (areal) [m]'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 3.0500000000000003
glyph1Display.SetScaleArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# rescale color and/or opacity maps used to exactly fit the current data range
glyph1Display.RescaleTransferFunctionToDataRange(False, True)

# get opacity transfer function/opacity map for 'Diameterarealm'
diameterarealmPWF = GetOpacityTransferFunction('Diameterarealm')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
diameterarealmLUT.ApplyPreset('$(ice_floes_color_scheme)', True)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
""")
        if save_animation
            write(f, """
SaveAnimation('$(folder)/$(simulation.id).avi', renderView1,
ImageResolution=[$(width), $(height)],
FrameRate=$(framerate),
FrameWindow=[0, $(simulation.file_number)])
""")
        end

        if save_images
            write(f, """
SaveAnimation('$(folder)/$(simulation.id).png', renderView1,
ImageResolution=[$(width), $(height)],
FrameRate=$(framerate),
FrameWindow=[0, $(simulation.file_number)])
""")
        end
    end
    if verbose
        info("$(filename) written, execute with " *
             "`cd $folder && pvpython $(simulation.id).py`")
    end
end

export render
"""
    render(simulation[, pvpython, images, animation])

Wrapper function which calls `writeParaviewPythonScript(...)` and executes it
from the shell using the supplied `pvpython` argument.

# Arguments
* `simulation::Simulation`: simulation object containing the ice-floe data.
* `pvpython::String`: path to the `pvpython` executable to use.  By default, the
    script uses the pvpython in the system PATH.
* `images::Bool`: render images to disk (default: true)
* `animation::Bool`: render animation to disk (default: false)
"""
function render(simulation::Simulation; pvpython::String="pvpython",
                images::Bool=true,
                animation::Bool=false)

    writeParaviewPythonScript(simulation, save_animation=animation,
                              save_images=images, verbose=false)
    try
        cd(simulation.id)
        run(`$(pvpython) $(simulation.id).py`)

        # if available, use imagemagick to create gif from images
        if images
            try
                run(`convert -trim +repage -delay 10 -transparent-color white 
                    -loop 0 $(simulation.id)*.png $(simulation.id).gif`)
            catch return_signal
                if isa(return_signal, Base.UVError)
                    error("skipping gif merge since `convert` was not found.")
                end
            end
        end
    catch return_signal
        if isa(return_signal, Base.UVError)
            error("`pvpython` was not found.")
        end
    end
end

export removeSimulationFiles
"""
    removeSimulationFiles(simulation[, folder])

Remove all simulation output files from the specified folder.
"""
function removeSimulationFiles(simulation::Simulation; folder::String=".")
    folder = folder * "/" * simulation.id
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtu"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtp"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vts"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).status.txt"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.jld"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).py"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).avi"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.png"`)
    run(`bash -c "rm -rf $(folder)"`)
    nothing
end
