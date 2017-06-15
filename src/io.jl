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

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    else
        return nothing
    end
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
    inter_particle_vector = Array{float, 1}[]
    force = float[]
    effective_radius = float[]
    contact_area = float[]
    contact_stiffness = float[]
    tensile_stress = float[]
    shear_displacement = Array{float, 1}[]
    contact_age = float[]
    for i=1:length(simulation.ice_floes)
        for ic=1:Nc_max
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

                if simulation.ice_floes[i].youngs_modulus > 0. &&
                    simulation.ice_floes[j].youngs_modulus > 0.
                    R_ij = harmonicMean(r_i, r_j)
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
            write(f, "$(shear_displacement[i][1]) ")
            write(f, "$(shear_displacement[i][2]) ")
            write(f, "0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" Name=\"Force [N]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(force[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Effective radius [m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(effective_radius[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact area [m^2]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(contact_area[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact stiffness [N/m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(contact_stiffness[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Tensile stress [Pa]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(tensile_stress[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact age [s]\" NumberOfComponents=\"1\" 
        format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(contact_age[i]) ")
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
            write(f, "$(i.lin_pos[1]) $(i.lin_pos[2]) 0.0 ")
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
            write(f, "$(i1[i] - 1) $(i2[i] - 1) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        
        # Write 0-indexed offset for the connectivity array for the end of each 
        # cell
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$((i - 1)*2 + 2) ")
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
        xq[:,:,iz] = grid.xq
        yq[:,:,iz] = grid.yq
    end
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            zq[ix,iy,:] = grid.zl
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
                vel[1, ix, iy, iz] = grid.u[ix, iy, iz, 1]
                vel[2, ix, iy, iz] = grid.v[ix, iy, iz, 1]
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
    else
        return nothing
    end
end

export writeParaviewStateFile
"""
Create a Paraview State File (.pvsm) for the simulation, which reads simulation 
output VTK files and applies appropriate glyph filters to the data.
"""
function writeParaviewStateFile(simulation::Simulation;
                                filename::String="",
                                folder::String=".",
                                vtk_folder::String=".",  # maybe expand to full path?
                                verbose::Bool=true)
    if filename == ""
        folder = folder * "/" * simulation.id
        mkpath(folder)
        filename = string(folder, "/", simulation.id, ".pvsm")
    end

    open(filename, "w") do f
        write(f, """<ParaView>
  <ServerManagerState version="4.4.0">
    <Proxy group="animation" type="AnimationScene" id="263" servers="16">
      <Property name="AnimationTime" id="263.AnimationTime" number_of_elements="1">
        <Element index="0" value="410"/>
      </Property>
      <Property name="Caching" id="263.Caching" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="263.Caching.bool"/>
      </Property>
      <Property name="Cues" id="263.Cues" number_of_elements="1">
        <Proxy value="265"/>
        <Domain name="groups" id="263.Cues.groups"/>
      </Property>
      <Property name="Duration" id="263.Duration" number_of_elements="1">
        <Element index="0" value="10"/>
      </Property>
      <Property name="EndTime" id="263.EndTime" number_of_elements="1">""")
        write(f, """<Element index="0" value="$(simulation.file_number)"/>""")
        write(f, """
      </Property>
      <Property name="FramesPerTimestep" id="263.FramesPerTimestep" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="263.FramesPerTimestep.range"/>
      </Property>
      <Property name="GoToFirst" id="263.GoToFirst"/>
      <Property name="GoToLast" id="263.GoToLast"/>
      <Property name="GoToNext" id="263.GoToNext"/>
      <Property name="GoToPrevious" id="263.GoToPrevious"/>
      <Property name="LockEndTime" id="263.LockEndTime" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="263.LockEndTime.bool"/>
      </Property>
      <Property name="LockStartTime" id="263.LockStartTime" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="263.LockStartTime.bool"/>
      </Property>
      <Property name="Loop" id="263.Loop" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="263.Loop.bool"/>
      </Property>
      <Property name="NumberOfFrames" id="263.NumberOfFrames" number_of_elements="1">
        <Element index="0" value="10"/>
        <Domain name="range" id="263.NumberOfFrames.range"/>
      </Property>
      <Property name="Play" id="263.Play"/>
      <Property name="PlayMode" id="263.PlayMode" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="enum" id="263.PlayMode.enum">
          <Entry value="0" text="Sequence"/>
          <Entry value="1" text="Real Time"/>
          <Entry value="2" text="Snap To TimeSteps"/>
        </Domain>
      </Property>
      <Property name="StartTime" id="263.StartTime" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
      <Property name="Stop" id="263.Stop"/>
      <Property name="TimeKeeper" id="263.TimeKeeper" number_of_elements="1">
        <Proxy value="259"/>
      </Property>
      <Property name="ViewModules" id="263.ViewModules" number_of_elements="2">
        <Proxy value="6135"/>
        <Proxy value="6135"/>
        <Domain name="groups" id="263.ViewModules.groups"/>
      </Property>
    </Proxy>
    <Proxy group="animation" type="TimeAnimationCue" id="265" servers="16">
      <Property name="AnimatedDomainName" id="265.AnimatedDomainName" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="AnimatedElement" id="265.AnimatedElement" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
      <Property name="AnimatedPropertyName" id="265.AnimatedPropertyName" number_of_elements="1">
        <Element index="0" value="Time"/>
      </Property>
      <Property name="AnimatedProxy" id="265.AnimatedProxy" number_of_elements="1">
        <Proxy value="259"/>
        <Domain name="groups" id="265.AnimatedProxy.groups"/>
      </Property>
      <Property name="Enabled" id="265.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="265.Enabled.bool"/>
      </Property>
      <Property name="EndTime" id="265.EndTime" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="KeyFrames" id="265.KeyFrames">
        <Domain name="groups" id="265.KeyFrames.groups"/>
      </Property>
      <Property name="LastAddedKeyFrameIndex" id="265.LastAddedKeyFrameIndex"/>
      <Property name="StartTime" id="265.StartTime" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
      <Property name="TimeMode" id="265.TimeMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="265.TimeMode.enum">
          <Entry value="0" text="Normalized"/>
          <Entry value="1" text="Relative"/>
        </Domain>
      </Property>
      <Property name="UseAnimationTime" id="265.UseAnimationTime" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="265.UseAnimationTime.bool"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="ViewLayout" id="6136" servers="16">
      <Layout number_of_elements="1">
        <Item direction="0" fraction="0.5" view="6135"/>
      </Layout>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5797" servers="21">
      <Property name="AboveRangeColor" id="5797.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5797.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5797.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5797.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5797.Annotations"/>
      <Property name="BelowRangeColor" id="5797.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5797.Build"/>
      <Property name="ColorSpace" id="5797.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5797.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5797.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5797.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5797.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5797.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5797.IndexedColors"/>
      <Property name="IndexedLookup" id="5797.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5797.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5797.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5797.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5797.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5797.RGBPoints" number_of_elements="12">
        <Element index="0" value="0.00025"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="0.0002500012500125"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="0.000250002500025"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5797.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5797.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5796"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5797.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5797.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5797.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5797.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5797.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5797.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5797.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5797.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5797.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5797.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5797.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5229" servers="21">
      <Property name="AboveRangeColor" id="5229.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5229.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5229.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5229.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5229.Annotations"/>
      <Property name="BelowRangeColor" id="5229.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5229.Build"/>
      <Property name="ColorSpace" id="5229.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5229.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5229.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5229.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5229.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5229.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5229.IndexedColors"/>
      <Property name="IndexedLookup" id="5229.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5229.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5229.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5229.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5229.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5229.RGBPoints" number_of_elements="12">
        <Element index="0" value="44"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="1364.87667828827"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="2685.75335652961"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5229.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5229.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5228"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5229.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5229.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5229.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5229.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5229.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5229.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5229.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5229.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5229.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5229.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5229.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5261" servers="21">
      <Property name="AboveRangeColor" id="5261.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5261.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5261.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5261.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5261.Annotations"/>
      <Property name="BelowRangeColor" id="5261.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5261.Build"/>
      <Property name="ColorSpace" id="5261.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5261.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5261.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5261.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5261.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5261.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5261.IndexedColors"/>
      <Property name="IndexedLookup" id="5261.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5261.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5261.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5261.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5261.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5261.RGBPoints" number_of_elements="12">
        <Element index="0" value="0"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="0.50000500005"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="1.0000100001"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5261.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5261.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5260"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5261.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5261.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5261.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5261.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5261.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5261.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5261.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5261.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5261.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5261.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5261.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5253" servers="21">
      <Property name="AboveRangeColor" id="5253.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5253.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5253.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5253.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5253.Annotations"/>
      <Property name="BelowRangeColor" id="5253.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5253.Build"/>
      <Property name="ColorSpace" id="5253.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5253.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5253.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5253.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5253.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5253.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5253.IndexedColors"/>
      <Property name="IndexedLookup" id="5253.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5253.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5253.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5253.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5253.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5253.RGBPoints" number_of_elements="12">
        <Element index="0" value="0"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="0.5"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="1"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5253.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5253.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5252"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5253.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5253.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5253.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5253.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5253.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5253.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5253.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5253.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5253.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5253.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5253.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5245" servers="21">
      <Property name="AboveRangeColor" id="5245.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5245.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5245.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5245.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5245.Annotations"/>
      <Property name="BelowRangeColor" id="5245.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5245.Build"/>
      <Property name="ColorSpace" id="5245.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5245.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5245.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5245.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5245.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5245.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5245.IndexedColors"/>
      <Property name="IndexedLookup" id="5245.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5245.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5245.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5245.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5245.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5245.RGBPoints" number_of_elements="12">
        <Element index="0" value="0"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="5e-17"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="1e-16"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5245.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5245.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5244"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5245.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5245.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5245.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5245.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5245.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5245.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5245.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5245.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5245.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5245.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5245.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5740" servers="21">
      <Property name="AboveRangeColor" id="5740.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5740.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5740.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5740.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5740.Annotations"/>
      <Property name="BelowRangeColor" id="5740.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5740.Build"/>
      <Property name="ColorSpace" id="5740.ColorSpace" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5740.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5740.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5740.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5740.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5740.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5740.IndexedColors"/>
      <Property name="IndexedLookup" id="5740.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5740.LockScalarRange" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5740.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5740.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5740.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5740.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5740.RGBPoints" number_of_elements="16">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0.999"/>
        <Element index="5" value="0"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0.501960784314"/>
        <Element index="8" value="1.998"/>
        <Element index="9" value="0"/>
        <Element index="10" value="0.501960784314"/>
        <Element index="11" value="1"/>
        <Element index="12" value="3"/>
        <Element index="13" value="1"/>
        <Element index="14" value="1"/>
        <Element index="15" value="1"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5740.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5740.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5739"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5740.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5740.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5740.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5740.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5740.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5740.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5740.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5740.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5740.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5740.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5740.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5237" servers="21">
      <Property name="AboveRangeColor" id="5237.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5237.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5237.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5237.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5237.Annotations"/>
      <Property name="BelowRangeColor" id="5237.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5237.Build"/>
      <Property name="ColorSpace" id="5237.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5237.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5237.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5237.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5237.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5237.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5237.IndexedColors"/>
      <Property name="IndexedLookup" id="5237.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5237.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5237.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5237.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5237.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5237.RGBPoints" number_of_elements="12">
        <Element index="0" value="0"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="2"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="4"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5237.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5237.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5236"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5237.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5237.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5237.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5237.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5237.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5237.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5237.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5237.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5237.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5237.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5237.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="6072" servers="21">
      <Property name="AboveRangeColor" id="6072.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="6072.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="6072.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6072.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="6072.Annotations"/>
      <Property name="BelowRangeColor" id="6072.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="6072.Build"/>
      <Property name="ColorSpace" id="6072.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="6072.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="6072.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6072.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="6072.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="6072.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="6072.IndexedColors"/>
      <Property name="IndexedLookup" id="6072.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="6072.LockScalarRange" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6072.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="6072.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="6072.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="6072.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="6072.RGBPoints" number_of_elements="12">
        <Element index="0" value="-500000"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="0"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="500000"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="6072.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="6072.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="6071"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="6072.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6072.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="6072.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="6072.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="6072.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="6072.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6072.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="6072.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6072.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="6072.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6072.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5299" servers="21">
      <Property name="AboveRangeColor" id="5299.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5299.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5299.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5299.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5299.Annotations"/>
      <Property name="BelowRangeColor" id="5299.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5299.Build"/>
      <Property name="ColorSpace" id="5299.ColorSpace" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5299.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5299.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5299.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5299.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5299.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5299.IndexedColors"/>
      <Property name="IndexedLookup" id="5299.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5299.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5299.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5299.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5299.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5299.RGBPoints" number_of_elements="16">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0.665553313626219"/>
        <Element index="5" value="0.501960784314"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="1.33110662725244"/>
        <Element index="9" value="1"/>
        <Element index="10" value="0.501960784314"/>
        <Element index="11" value="0"/>
        <Element index="12" value="1.99865859947813"/>
        <Element index="13" value="1"/>
        <Element index="14" value="1"/>
        <Element index="15" value="1"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5299.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5299.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5298"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5299.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5299.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5299.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5299.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5299.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5299.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5299.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5299.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5299.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5299.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5299.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5269" servers="21">
      <Property name="AboveRangeColor" id="5269.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5269.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5269.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5269.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5269.Annotations"/>
      <Property name="BelowRangeColor" id="5269.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5269.Build"/>
      <Property name="ColorSpace" id="5269.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="5269.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5269.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5269.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5269.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5269.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5269.IndexedColors"/>
      <Property name="IndexedLookup" id="5269.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5269.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5269.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5269.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5269.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5269.RGBPoints" number_of_elements="12">
        <Element index="0" value="0"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="5e-17"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="1e-16"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5269.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5269.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5268"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5269.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5269.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5269.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5269.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5269.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5269.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5269.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5269.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5269.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5269.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5269.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="5888" servers="21">
      <Property name="AboveRangeColor" id="5888.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="5888.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="5888.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5888.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="5888.Annotations"/>
      <Property name="BelowRangeColor" id="5888.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="5888.Build"/>
      <Property name="ColorSpace" id="5888.ColorSpace" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="enum" id="5888.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="5888.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5888.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="5888.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="5888.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="5888.IndexedColors"/>
      <Property name="IndexedLookup" id="5888.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="5888.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="5888.NanColor" number_of_elements="3">
        <Element index="0" value="0.25"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="5888.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="5888.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="5888.RGBPoints" number_of_elements="84">
        <Element index="0" value="0"/>
        <Element index="1" value="0.960784"/>
        <Element index="2" value="1"/>
        <Element index="3" value="0.980392"/>
        <Element index="4" value="5e-18"/>
        <Element index="5" value="0.815686"/>
        <Element index="6" value="0.960784"/>
        <Element index="7" value="0.913725"/>
        <Element index="8" value="1e-17"/>
        <Element index="9" value="0.670588"/>
        <Element index="10" value="0.929412"/>
        <Element index="11" value="0.870588"/>
        <Element index="12" value="1.5e-17"/>
        <Element index="13" value="0.556863"/>
        <Element index="14" value="0.901961"/>
        <Element index="15" value="0.843137"/>
        <Element index="16" value="2e-17"/>
        <Element index="17" value="0.478431"/>
        <Element index="18" value="0.870588"/>
        <Element index="19" value="0.823529"/>
        <Element index="20" value="2.5e-17"/>
        <Element index="21" value="0.439216"/>
        <Element index="22" value="0.831373"/>
        <Element index="23" value="0.803922"/>
        <Element index="24" value="3e-17"/>
        <Element index="25" value="0.4"/>
        <Element index="26" value="0.8"/>
        <Element index="27" value="0.788235"/>
        <Element index="28" value="3.5e-17"/>
        <Element index="29" value="0.376471"/>
        <Element index="30" value="0.768627"/>
        <Element index="31" value="0.768627"/>
        <Element index="32" value="4e-17"/>
        <Element index="33" value="0.34902"/>
        <Element index="34" value="0.709804"/>
        <Element index="35" value="0.729412"/>
        <Element index="36" value="4.5e-17"/>
        <Element index="37" value="0.32549"/>
        <Element index="38" value="0.654902"/>
        <Element index="39" value="0.690196"/>
        <Element index="40" value="5e-17"/>
        <Element index="41" value="0.301961"/>
        <Element index="42" value="0.607843"/>
        <Element index="43" value="0.658824"/>
        <Element index="44" value="5.5e-17"/>
        <Element index="45" value="0.247059"/>
        <Element index="46" value="0.545098"/>
        <Element index="47" value="0.619608"/>
        <Element index="48" value="6e-17"/>
        <Element index="49" value="0.239216"/>
        <Element index="50" value="0.494118"/>
        <Element index="51" value="0.580392"/>
        <Element index="52" value="6.5e-17"/>
        <Element index="53" value="0.227451"/>
        <Element index="54" value="0.439216"/>
        <Element index="55" value="0.541176"/>
        <Element index="56" value="7e-17"/>
        <Element index="57" value="0.227451"/>
        <Element index="58" value="0.403922"/>
        <Element index="59" value="0.521569"/>
        <Element index="60" value="7.5e-17"/>
        <Element index="61" value="0.231373"/>
        <Element index="62" value="0.368627"/>
        <Element index="63" value="0.501961"/>
        <Element index="64" value="8e-17"/>
        <Element index="65" value="0.227451"/>
        <Element index="66" value="0.321569"/>
        <Element index="67" value="0.470588"/>
        <Element index="68" value="8.5e-17"/>
        <Element index="69" value="0.219608"/>
        <Element index="70" value="0.282353"/>
        <Element index="71" value="0.439216"/>
        <Element index="72" value="9e-17"/>
        <Element index="73" value="0.192157"/>
        <Element index="74" value="0.235294"/>
        <Element index="75" value="0.4"/>
        <Element index="76" value="9.5e-17"/>
        <Element index="77" value="0.160784"/>
        <Element index="78" value="0.184314"/>
        <Element index="79" value="0.34902"/>
        <Element index="80" value="1e-16"/>
        <Element index="81" value="0.133333"/>
        <Element index="82" value="0.12549"/>
        <Element index="83" value="0.301961"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="5888.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5888.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5887"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5888.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5888.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="5888.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="5888.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="5888.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="5888.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5888.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="5888.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5888.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="5888.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5888.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="lookup_tables" type="PVLookupTable" id="6693" servers="21">
      <Property name="AboveRangeColor" id="6693.AboveRangeColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
      <Property name="ActiveAnnotatedValues" id="6693.ActiveAnnotatedValues"/>
      <Property name="AllowDuplicateScalars" id="6693.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6693.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Annotations" id="6693.Annotations"/>
      <Property name="BelowRangeColor" id="6693.BelowRangeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Build" id="6693.Build"/>
      <Property name="ColorSpace" id="6693.ColorSpace" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="enum" id="6693.ColorSpace.enum">
          <Entry value="0" text="RGB"/>
          <Entry value="1" text="HSV"/>
          <Entry value="2" text="Lab"/>
          <Entry value="3" text="Diverging"/>
        </Domain>
      </Property>
      <Property name="Discretize" id="6693.Discretize" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6693.Discretize.bool"/>
      </Property>
      <Property name="EnableOpacityMapping" id="6693.EnableOpacityMapping" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.EnableOpacityMapping.bool"/>
      </Property>
      <Property name="HSVWrap" id="6693.HSVWrap" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.HSVWrap.bool"/>
      </Property>
      <Property name="IndexedColors" id="6693.IndexedColors"/>
      <Property name="IndexedLookup" id="6693.IndexedLookup" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.IndexedLookup.bool"/>
      </Property>
      <Property name="LockScalarRange" id="6693.LockScalarRange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.LockScalarRange.bool"/>
      </Property>
      <Property name="NanColor" id="6693.NanColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="NumberOfTableValues" id="6693.NumberOfTableValues" number_of_elements="1">
        <Element index="0" value="256"/>
        <Domain name="range" id="6693.NumberOfTableValues.range"/>
      </Property>
      <Property name="RGBPoints" id="6693.RGBPoints" number_of_elements="12">
        <Element index="0" value="-30"/>
        <Element index="1" value="0.231373"/>
        <Element index="2" value="0.298039"/>
        <Element index="3" value="0.752941"/>
        <Element index="4" value="-15"/>
        <Element index="5" value="0.865003"/>
        <Element index="6" value="0.865003"/>
        <Element index="7" value="0.865003"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0.705882"/>
        <Element index="10" value="0.0156863"/>
        <Element index="11" value="0.14902"/>
      </Property>
      <Property name="RescaleOnVisibilityChange" id="6693.RescaleOnVisibilityChange" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.RescaleOnVisibilityChange.bool"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="6693.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="6692"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="6693.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6693.ScalarRangeInitialized.bool"/>
      </Property>
      <Property name="ShowIndexedColorActiveValues" id="6693.ShowIndexedColorActiveValues" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.ShowIndexedColorActiveValues.bool"/>
      </Property>
      <Property name="UseAboveRangeColor" id="6693.UseAboveRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.UseAboveRangeColor.bool"/>
      </Property>
      <Property name="UseBelowRangeColor" id="6693.UseBelowRangeColor" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.UseBelowRangeColor.bool"/>
      </Property>
      <Property name="UseLogScale" id="6693.UseLogScale" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6693.UseLogScale.bool"/>
      </Property>
      <Property name="VectorComponent" id="6693.VectorComponent" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6693.VectorComponent.range"/>
      </Property>
      <Property name="VectorMode" id="6693.VectorMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6693.VectorMode.enum">
          <Entry value="0" text="Magnitude"/>
          <Entry value="1" text="Component"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5796" servers="21">
      <Property name="AllowDuplicateScalars" id="5796.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5796.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5796.Points" number_of_elements="8">
        <Element index="0" value="0.00025"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0.000250002500025"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5796.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5796.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5228" servers="21">
      <Property name="AllowDuplicateScalars" id="5228.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5228.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5228.Points" number_of_elements="8">
        <Element index="0" value="44"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="2685.75335652961"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5228.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5228.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5260" servers="21">
      <Property name="AllowDuplicateScalars" id="5260.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5260.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5260.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1.0000100001"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5260.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5260.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5252" servers="21">
      <Property name="AllowDuplicateScalars" id="5252.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5252.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5252.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5252.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5252.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5244" servers="21">
      <Property name="AllowDuplicateScalars" id="5244.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5244.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5244.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1e-16"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5244.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5244.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5739" servers="21">
      <Property name="AllowDuplicateScalars" id="5739.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5739.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5739.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="3"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5739.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5739.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5236" servers="21">
      <Property name="AllowDuplicateScalars" id="5236.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5236.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5236.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="4"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5236.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5236.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="6071" servers="21">
      <Property name="AllowDuplicateScalars" id="6071.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6071.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="6071.Points" number_of_elements="8">
        <Element index="0" value="-500000"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="500000"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="6071.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6071.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5298" servers="21">
      <Property name="AllowDuplicateScalars" id="5298.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5298.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5298.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1.99865859947813"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5298.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5298.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5268" servers="21">
      <Property name="AllowDuplicateScalars" id="5268.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5268.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5268.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1e-16"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5268.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5268.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="5887" servers="21">
      <Property name="AllowDuplicateScalars" id="5887.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5887.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="5887.Points" number_of_elements="8">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="1e-16"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="5887.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5887.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="piecewise_functions" type="PiecewiseFunction" id="6692" servers="21">
      <Property name="AllowDuplicateScalars" id="6692.AllowDuplicateScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6692.AllowDuplicateScalars.bool"/>
      </Property>
      <Property name="Points" id="6692.Points" number_of_elements="8">
        <Element index="0" value="-30"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0.5"/>
        <Element index="7" value="0"/>
      </Property>
      <Property name="ScalarRangeInitialized" id="6692.ScalarRangeInitialized" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6692.ScalarRangeInitialized.bool"/>
      </Property>
    </Proxy>
    <Proxy group="annotations" type="GridAxes3DActor" id="5227" servers="21">
      <Property name="AxesToLabel" id="5227.AxesToLabel" number_of_elements="1">
        <Element index="0" value="63"/>
        <Domain name="range" id="5227.AxesToLabel.range"/>
      </Property>
      <Property name="DataPosition" id="5227.DataPosition" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.DataPosition.range"/>
      </Property>
      <Property name="DataScale" id="5227.DataScale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5227.DataScale.range"/>
      </Property>
      <Property name="FacesToRender" id="5227.FacesToRender" number_of_elements="1">
        <Element index="0" value="63"/>
        <Domain name="range" id="5227.FacesToRender.range"/>
      </Property>
      <Property name="LabelUniqueEdgesOnly" id="5227.LabelUniqueEdgesOnly" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5227.LabelUniqueEdgesOnly.bool"/>
      </Property>
      <Property name="ModelBounds" id="5227.ModelBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="ModelTransformMatrix" id="5227.ModelTransformMatrix" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
      <Property name="ShowEdges" id="5227.ShowEdges" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5227.ShowEdges.bool"/>
      </Property>
      <Property name="ShowGrid" id="5227.ShowGrid" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ShowGrid.bool"/>
      </Property>
      <Property name="ShowTicks" id="5227.ShowTicks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5227.ShowTicks.bool"/>
      </Property>
      <Property name="UseModelTransform" id="5227.UseModelTransform" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5227.UseModelTransform.range"/>
      </Property>
      <Property name="Visibility" id="5227.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5227.Visibility.bool"/>
      </Property>
      <Property name="XAxisLabels" id="5227.XAxisLabels">
        <Domain name="scalar_range" id="5227.XAxisLabels.scalar_range"/>
      </Property>
      <Property name="XAxisNotation" id="5227.XAxisNotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.XAxisNotation.enum">
          <Entry value="0" text="Mixed"/>
          <Entry value="1" text="Scientific"/>
          <Entry value="2" text="Fixed"/>
        </Domain>
      </Property>
      <Property name="XAxisPrecision" id="5227.XAxisPrecision" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5227.XAxisPrecision.range"/>
      </Property>
      <Property name="XAxisUseCustomLabels" id="5227.XAxisUseCustomLabels" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XAxisUseCustomLabels.bool"/>
      </Property>
      <Property name="XTitle" id="5227.XTitle" number_of_elements="1">
        <Element index="0" value="X Axis"/>
      </Property>
      <Property name="YAxisLabels" id="5227.YAxisLabels">
        <Domain name="scalar_range" id="5227.YAxisLabels.scalar_range"/>
      </Property>
      <Property name="YAxisNotation" id="5227.YAxisNotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.YAxisNotation.enum">
          <Entry value="0" text="Mixed"/>
          <Entry value="1" text="Scientific"/>
          <Entry value="2" text="Fixed"/>
        </Domain>
      </Property>
      <Property name="YAxisPrecision" id="5227.YAxisPrecision" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5227.YAxisPrecision.range"/>
      </Property>
      <Property name="YAxisUseCustomLabels" id="5227.YAxisUseCustomLabels" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YAxisUseCustomLabels.bool"/>
      </Property>
      <Property name="YTitle" id="5227.YTitle" number_of_elements="1">
        <Element index="0" value="Y Axis"/>
      </Property>
      <Property name="ZAxisLabels" id="5227.ZAxisLabels">
        <Domain name="scalar_range" id="5227.ZAxisLabels.scalar_range"/>
      </Property>
      <Property name="ZAxisNotation" id="5227.ZAxisNotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.ZAxisNotation.enum">
          <Entry value="0" text="Mixed"/>
          <Entry value="1" text="Scientific"/>
          <Entry value="2" text="Fixed"/>
        </Domain>
      </Property>
      <Property name="ZAxisPrecision" id="5227.ZAxisPrecision" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5227.ZAxisPrecision.range"/>
      </Property>
      <Property name="ZAxisUseCustomLabels" id="5227.ZAxisUseCustomLabels" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZAxisUseCustomLabels.bool"/>
      </Property>
      <Property name="ZTitle" id="5227.ZTitle" number_of_elements="1">
        <Element index="0" value="Z Axis"/>
      </Property>
      <Property name="CullBackface" id="5227.CullBackface" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.CullBackface.bool"/>
      </Property>
      <Property name="CullFrontface" id="5227.CullFrontface" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5227.CullFrontface.bool"/>
      </Property>
      <Property name="GridColor" id="5227.GridColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.GridColor.range"/>
      </Property>
      <Property name="XLabelBold" id="5227.XLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XLabelBold.bool"/>
      </Property>
      <Property name="XLabelColor" id="5227.XLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.XLabelColor.range"/>
      </Property>
      <Property name="XLabelFontFamily" id="5227.XLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.XLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="XLabelFontSize" id="5227.XLabelFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.XLabelFontSize.range"/>
      </Property>
      <Property name="XLabelItalic" id="5227.XLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XLabelItalic.bool"/>
      </Property>
      <Property name="XLabelOpacity" id="5227.XLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.XLabelOpacity.range"/>
      </Property>
      <Property name="XLabelShadow" id="5227.XLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XLabelShadow.bool"/>
      </Property>
      <Property name="XTitleBold" id="5227.XTitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XTitleBold.bool"/>
      </Property>
      <Property name="XTitleColor" id="5227.XTitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.XTitleColor.range"/>
      </Property>
      <Property name="XTitleFontFamily" id="5227.XTitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.XTitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="XTitleFontSize" id="5227.XTitleFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.XTitleFontSize.range"/>
      </Property>
      <Property name="XTitleItalic" id="5227.XTitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XTitleItalic.bool"/>
      </Property>
      <Property name="XTitleOpacity" id="5227.XTitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.XTitleOpacity.range"/>
      </Property>
      <Property name="XTitleShadow" id="5227.XTitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.XTitleShadow.bool"/>
      </Property>
      <Property name="YLabelBold" id="5227.YLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YLabelBold.bool"/>
      </Property>
      <Property name="YLabelColor" id="5227.YLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.YLabelColor.range"/>
      </Property>
      <Property name="YLabelFontFamily" id="5227.YLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.YLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="YLabelFontSize" id="5227.YLabelFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.YLabelFontSize.range"/>
      </Property>
      <Property name="YLabelItalic" id="5227.YLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YLabelItalic.bool"/>
      </Property>
      <Property name="YLabelOpacity" id="5227.YLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.YLabelOpacity.range"/>
      </Property>
      <Property name="YLabelShadow" id="5227.YLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YLabelShadow.bool"/>
      </Property>
      <Property name="YTitleBold" id="5227.YTitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YTitleBold.bool"/>
      </Property>
      <Property name="YTitleColor" id="5227.YTitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.YTitleColor.range"/>
      </Property>
      <Property name="YTitleFontFamily" id="5227.YTitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.YTitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="YTitleFontSize" id="5227.YTitleFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.YTitleFontSize.range"/>
      </Property>
      <Property name="YTitleItalic" id="5227.YTitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YTitleItalic.bool"/>
      </Property>
      <Property name="YTitleOpacity" id="5227.YTitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.YTitleOpacity.range"/>
      </Property>
      <Property name="YTitleShadow" id="5227.YTitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.YTitleShadow.bool"/>
      </Property>
      <Property name="ZLabelBold" id="5227.ZLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZLabelBold.bool"/>
      </Property>
      <Property name="ZLabelColor" id="5227.ZLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.ZLabelColor.range"/>
      </Property>
      <Property name="ZLabelFontFamily" id="5227.ZLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.ZLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="ZLabelFontSize" id="5227.ZLabelFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.ZLabelFontSize.range"/>
      </Property>
      <Property name="ZLabelItalic" id="5227.ZLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZLabelItalic.bool"/>
      </Property>
      <Property name="ZLabelOpacity" id="5227.ZLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.ZLabelOpacity.range"/>
      </Property>
      <Property name="ZLabelShadow" id="5227.ZLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZLabelShadow.bool"/>
      </Property>
      <Property name="ZTitleBold" id="5227.ZTitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZTitleBold.bool"/>
      </Property>
      <Property name="ZTitleColor" id="5227.ZTitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5227.ZTitleColor.range"/>
      </Property>
      <Property name="ZTitleFontFamily" id="5227.ZTitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5227.ZTitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="ZTitleFontSize" id="5227.ZTitleFontSize" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5227.ZTitleFontSize.range"/>
      </Property>
      <Property name="ZTitleItalic" id="5227.ZTitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZTitleItalic.bool"/>
      </Property>
      <Property name="ZTitleOpacity" id="5227.ZTitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5227.ZTitleOpacity.range"/>
      </Property>
      <Property name="ZTitleShadow" id="5227.ZTitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5227.ZTitleShadow.bool"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6182" servers="16">
      <Property name="Source" id="6182.Source" number_of_elements="1">
        <Proxy value="4765"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6183" servers="16">
      <Property name="Source" id="6183.Source" number_of_elements="1">
        <Proxy value="4787"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6184" servers="16">
      <Property name="Source" id="6184.Source" number_of_elements="1">
        <Proxy value="4809"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6185" servers="16">
      <Property name="Source" id="6185.Source" number_of_elements="1">
        <Proxy value="4831"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="4842" servers="1">
      <Property name="Position" id="4842.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4842.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="4842.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="4842.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4842.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="4842.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="4842.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="4842.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="4842.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6186" servers="16">
      <Property name="Source" id="6186.Source" number_of_elements="1">
        <Proxy value="4920"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="4843" servers="1">
      <Property name="Invert" id="4843.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="4843.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="4843.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="4843.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="4843.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4843.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="4843.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="4843.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="4843.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="4843.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="4843.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4843.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="4854" servers="1">
      <Property name="Capping" id="4854.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="4854.Capping.bool"/>
      </Property>
      <Property name="Center" id="4854.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4854.Center.range"/>
      </Property>
      <Property name="Direction" id="4854.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4854.Direction.range"/>
      </Property>
      <Property name="Height" id="4854.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4854.Height.range"/>
      </Property>
      <Property name="Radius" id="4854.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4854.Radius.range"/>
      </Property>
      <Property name="Resolution" id="4854.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4854.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="4865" servers="1">
      <Property name="Center" id="4865.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4865.Center.range"/>
      </Property>
      <Property name="XLength" id="4865.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4865.XLength.range"/>
      </Property>
      <Property name="YLength" id="4865.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4865.YLength.range"/>
      </Property>
      <Property name="ZLength" id="4865.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4865.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="4876" servers="1">
      <Property name="Capping" id="4876.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="4876.Capping.bool"/>
      </Property>
      <Property name="Center" id="4876.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4876.Center.range"/>
      </Property>
      <Property name="Height" id="4876.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4876.Height.range"/>
      </Property>
      <Property name="Radius" id="4876.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4876.Radius.range"/>
      </Property>
      <Property name="Resolution" id="4876.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4876.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="4887" servers="1">
      <Property name="Point1" id="4887.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4887.Point1.range"/>
      </Property>
      <Property name="Point2" id="4887.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4887.Point2.range"/>
      </Property>
      <Property name="Resolution" id="4887.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4887.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="4898" servers="1">
      <Property name="Center" id="4898.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4898.Center.range"/>
      </Property>
      <Property name="EndPhi" id="4898.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="4898.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="4898.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="4898.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="4898.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="4898.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="4898.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4898.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="4898.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="4898.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="4898.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="4898.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="4898.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="4898.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="4909" servers="1">
      <Property name="Center" id="4909.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4909.Center.range"/>
      </Property>
      <Property name="Filled" id="4909.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="4909.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="4909.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="4909.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="4931" servers="1">
      <Property name="Position" id="4931.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4931.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="4931.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="4931.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4931.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="4931.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="4931.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="4931.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="4931.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6187" servers="16">
      <Property name="Source" id="6187.Source" number_of_elements="1">
        <Proxy value="5009"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="4932" servers="1">
      <Property name="Invert" id="4932.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="4932.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="4932.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="4932.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="4932.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4932.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="4932.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="4932.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="4932.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="4932.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="4932.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4932.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="4943" servers="1">
      <Property name="Capping" id="4943.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="4943.Capping.bool"/>
      </Property>
      <Property name="Center" id="4943.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4943.Center.range"/>
      </Property>
      <Property name="Direction" id="4943.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4943.Direction.range"/>
      </Property>
      <Property name="Height" id="4943.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4943.Height.range"/>
      </Property>
      <Property name="Radius" id="4943.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4943.Radius.range"/>
      </Property>
      <Property name="Resolution" id="4943.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4943.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="4954" servers="1">
      <Property name="Center" id="4954.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4954.Center.range"/>
      </Property>
      <Property name="XLength" id="4954.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4954.XLength.range"/>
      </Property>
      <Property name="YLength" id="4954.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4954.YLength.range"/>
      </Property>
      <Property name="ZLength" id="4954.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4954.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="4965" servers="1">
      <Property name="Capping" id="4965.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="4965.Capping.bool"/>
      </Property>
      <Property name="Center" id="4965.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4965.Center.range"/>
      </Property>
      <Property name="Height" id="4965.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="4965.Height.range"/>
      </Property>
      <Property name="Radius" id="4965.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4965.Radius.range"/>
      </Property>
      <Property name="Resolution" id="4965.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4965.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="4976" servers="1">
      <Property name="Point1" id="4976.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4976.Point1.range"/>
      </Property>
      <Property name="Point2" id="4976.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4976.Point2.range"/>
      </Property>
      <Property name="Resolution" id="4976.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="4976.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="4987" servers="1">
      <Property name="Center" id="4987.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4987.Center.range"/>
      </Property>
      <Property name="EndPhi" id="4987.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="4987.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="4987.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="4987.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="4987.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="4987.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="4987.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="4987.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="4987.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="4987.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="4987.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="4987.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="4987.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="4987.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="4998" servers="1">
      <Property name="Center" id="4998.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="4998.Center.range"/>
      </Property>
      <Property name="Filled" id="4998.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="4998.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="4998.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="4998.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="5020" servers="1">
      <Property name="Position" id="5020.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5020.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="5020.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="5020.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5020.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="5020.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="5020.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5020.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="5020.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6188" servers="16">
      <Property name="Source" id="6188.Source" number_of_elements="1">
        <Proxy value="5098"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="5021" servers="1">
      <Property name="Invert" id="5021.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5021.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="5021.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="5021.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="5021.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5021.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="5021.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="5021.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="5021.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="5021.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="5021.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5021.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="5032" servers="1">
      <Property name="Capping" id="5032.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5032.Capping.bool"/>
      </Property>
      <Property name="Center" id="5032.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5032.Center.range"/>
      </Property>
      <Property name="Direction" id="5032.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5032.Direction.range"/>
      </Property>
      <Property name="Height" id="5032.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5032.Height.range"/>
      </Property>
      <Property name="Radius" id="5032.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5032.Radius.range"/>
      </Property>
      <Property name="Resolution" id="5032.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5032.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="5043" servers="1">
      <Property name="Center" id="5043.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5043.Center.range"/>
      </Property>
      <Property name="XLength" id="5043.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5043.XLength.range"/>
      </Property>
      <Property name="YLength" id="5043.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5043.YLength.range"/>
      </Property>
      <Property name="ZLength" id="5043.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5043.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="5054" servers="1">
      <Property name="Capping" id="5054.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5054.Capping.bool"/>
      </Property>
      <Property name="Center" id="5054.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5054.Center.range"/>
      </Property>
      <Property name="Height" id="5054.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5054.Height.range"/>
      </Property>
      <Property name="Radius" id="5054.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5054.Radius.range"/>
      </Property>
      <Property name="Resolution" id="5054.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5054.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="5065" servers="1">
      <Property name="Point1" id="5065.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5065.Point1.range"/>
      </Property>
      <Property name="Point2" id="5065.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5065.Point2.range"/>
      </Property>
      <Property name="Resolution" id="5065.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5065.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="5076" servers="1">
      <Property name="Center" id="5076.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5076.Center.range"/>
      </Property>
      <Property name="EndPhi" id="5076.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="5076.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="5076.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="5076.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="5076.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="5076.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="5076.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5076.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="5076.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5076.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="5076.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5076.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="5076.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="5076.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="5087" servers="1">
      <Property name="Center" id="5087.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5087.Center.range"/>
      </Property>
      <Property name="Filled" id="5087.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5087.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="5087.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="5087.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="5109" servers="1">
      <Property name="Position" id="5109.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5109.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="5109.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="5109.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5109.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="5109.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="5109.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5109.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="5109.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6189" servers="16">
      <Property name="Source" id="6189.Source" number_of_elements="1">
        <Proxy value="5187"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="5110" servers="1">
      <Property name="Invert" id="5110.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5110.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="5110.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="5110.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="5110.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5110.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="5110.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="5110.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="5110.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="5110.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="5110.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5110.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="5121" servers="1">
      <Property name="Capping" id="5121.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5121.Capping.bool"/>
      </Property>
      <Property name="Center" id="5121.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5121.Center.range"/>
      </Property>
      <Property name="Direction" id="5121.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5121.Direction.range"/>
      </Property>
      <Property name="Height" id="5121.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5121.Height.range"/>
      </Property>
      <Property name="Radius" id="5121.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5121.Radius.range"/>
      </Property>
      <Property name="Resolution" id="5121.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5121.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="5132" servers="1">
      <Property name="Center" id="5132.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5132.Center.range"/>
      </Property>
      <Property name="XLength" id="5132.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5132.XLength.range"/>
      </Property>
      <Property name="YLength" id="5132.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5132.YLength.range"/>
      </Property>
      <Property name="ZLength" id="5132.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5132.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="5143" servers="1">
      <Property name="Capping" id="5143.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5143.Capping.bool"/>
      </Property>
      <Property name="Center" id="5143.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5143.Center.range"/>
      </Property>
      <Property name="Height" id="5143.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5143.Height.range"/>
      </Property>
      <Property name="Radius" id="5143.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5143.Radius.range"/>
      </Property>
      <Property name="Resolution" id="5143.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5143.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="5154" servers="1">
      <Property name="Point1" id="5154.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5154.Point1.range"/>
      </Property>
      <Property name="Point2" id="5154.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5154.Point2.range"/>
      </Property>
      <Property name="Resolution" id="5154.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="5154.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="5165" servers="1">
      <Property name="Center" id="5165.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5165.Center.range"/>
      </Property>
      <Property name="EndPhi" id="5165.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="5165.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="5165.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="5165.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="5165.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="5165.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="5165.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="5165.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="5165.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5165.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="5165.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5165.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="5165.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="5165.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="5176" servers="1">
      <Property name="Center" id="5176.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5176.Center.range"/>
      </Property>
      <Property name="Filled" id="5176.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5176.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="5176.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="5176.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6190" servers="16">
      <Property name="Source" id="6190.Source" number_of_elements="1">
        <Proxy value="5198"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6191" servers="16">
      <Property name="Source" id="6191.Source" number_of_elements="1">
        <Proxy value="5209"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6337" servers="16">
      <Property name="Source" id="6337.Source" number_of_elements="1">
        <Proxy value="6326"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform3" id="6325" servers="1">
      <Property name="Position" id="6325.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6325.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="6325.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="6325.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6325.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="6325.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="6325.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6325.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="6325.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="6525" servers="1">
      <Property name="Position" id="6525.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6525.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="6525.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="6525.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6525.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="6525.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="6525.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6525.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="6525.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6614" servers="16">
      <Property name="Source" id="6614.Source" number_of_elements="1">
        <Proxy value="6603"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="6526" servers="1">
      <Property name="Invert" id="6526.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6526.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="6526.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="6526.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="6526.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6526.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="6526.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="6526.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="6526.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="6526.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="6526.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6526.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="6537" servers="1">
      <Property name="Capping" id="6537.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6537.Capping.bool"/>
      </Property>
      <Property name="Center" id="6537.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6537.Center.range"/>
      </Property>
      <Property name="Direction" id="6537.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6537.Direction.range"/>
      </Property>
      <Property name="Height" id="6537.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6537.Height.range"/>
      </Property>
      <Property name="Radius" id="6537.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6537.Radius.range"/>
      </Property>
      <Property name="Resolution" id="6537.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6537.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="6548" servers="1">
      <Property name="Center" id="6548.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6548.Center.range"/>
      </Property>
      <Property name="XLength" id="6548.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6548.XLength.range"/>
      </Property>
      <Property name="YLength" id="6548.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6548.YLength.range"/>
      </Property>
      <Property name="ZLength" id="6548.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6548.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="6559" servers="1">
      <Property name="Capping" id="6559.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6559.Capping.bool"/>
      </Property>
      <Property name="Center" id="6559.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6559.Center.range"/>
      </Property>
      <Property name="Height" id="6559.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6559.Height.range"/>
      </Property>
      <Property name="Radius" id="6559.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6559.Radius.range"/>
      </Property>
      <Property name="Resolution" id="6559.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6559.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="6570" servers="1">
      <Property name="Point1" id="6570.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6570.Point1.range"/>
      </Property>
      <Property name="Point2" id="6570.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6570.Point2.range"/>
      </Property>
      <Property name="Resolution" id="6570.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6570.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="6581" servers="1">
      <Property name="Center" id="6581.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6581.Center.range"/>
      </Property>
      <Property name="EndPhi" id="6581.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="6581.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="6581.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="6581.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="6581.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="6581.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="6581.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6581.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="6581.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6581.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="6581.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6581.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="6581.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="6581.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="6592" servers="1">
      <Property name="Center" id="6592.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6592.Center.range"/>
      </Property>
      <Property name="Filled" id="6592.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6592.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="6592.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="6592.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6712" servers="16">
      <Property name="Source" id="6712.Source" number_of_elements="1">
        <Proxy value="6701"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform3" id="6700" servers="1">
      <Property name="Position" id="6700.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6700.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="6700.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="6700.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6700.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="6700.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="6700.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6700.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="6700.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6805" servers="16">
      <Property name="Source" id="6805.Source" number_of_elements="1">
        <Proxy value="6794"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6894" servers="16">
      <Property name="Source" id="6894.Source" number_of_elements="1">
        <Proxy value="6883"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="6984" servers="16">
      <Property name="Source" id="6984.Source" number_of_elements="1">
        <Proxy value="6973"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform3" id="6972" servers="1">
      <Property name="Position" id="6972.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6972.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="6972.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="6972.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6972.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="6972.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="6972.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6972.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="6972.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="7082" servers="1">
      <Property name="Position" id="7082.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7082.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="7082.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="7082.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7082.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="7082.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="7082.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7082.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="7082.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="7171" servers="16">
      <Property name="Source" id="7171.Source" number_of_elements="1">
        <Proxy value="7160"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="7083" servers="1">
      <Property name="Invert" id="7083.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7083.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="7083.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="7083.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="7083.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7083.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="7083.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="7083.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="7083.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="7083.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="7083.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7083.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="7094" servers="1">
      <Property name="Capping" id="7094.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7094.Capping.bool"/>
      </Property>
      <Property name="Center" id="7094.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7094.Center.range"/>
      </Property>
      <Property name="Direction" id="7094.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7094.Direction.range"/>
      </Property>
      <Property name="Height" id="7094.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7094.Height.range"/>
      </Property>
      <Property name="Radius" id="7094.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7094.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7094.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7094.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="7105" servers="1">
      <Property name="Center" id="7105.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7105.Center.range"/>
      </Property>
      <Property name="XLength" id="7105.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7105.XLength.range"/>
      </Property>
      <Property name="YLength" id="7105.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7105.YLength.range"/>
      </Property>
      <Property name="ZLength" id="7105.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7105.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="7116" servers="1">
      <Property name="Capping" id="7116.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7116.Capping.bool"/>
      </Property>
      <Property name="Center" id="7116.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7116.Center.range"/>
      </Property>
      <Property name="Height" id="7116.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7116.Height.range"/>
      </Property>
      <Property name="Radius" id="7116.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7116.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7116.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7116.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="7127" servers="1">
      <Property name="Point1" id="7127.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7127.Point1.range"/>
      </Property>
      <Property name="Point2" id="7127.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7127.Point2.range"/>
      </Property>
      <Property name="Resolution" id="7127.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7127.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="7138" servers="1">
      <Property name="Center" id="7138.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7138.Center.range"/>
      </Property>
      <Property name="EndPhi" id="7138.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="7138.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="7138.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="7138.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="7138.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7138.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="7138.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7138.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="7138.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7138.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="7138.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7138.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="7138.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7138.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="7149" servers="1">
      <Property name="Center" id="7149.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7149.Center.range"/>
      </Property>
      <Property name="Filled" id="7149.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7149.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="7149.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="7149.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="7249" servers="1">
      <Property name="Position" id="7249.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7249.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="7249.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="7249.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7249.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="7249.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="7249.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7249.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="7249.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="7338" servers="16">
      <Property name="Source" id="7338.Source" number_of_elements="1">
        <Proxy value="7327"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="7250" servers="1">
      <Property name="Invert" id="7250.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7250.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="7250.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="7250.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="7250.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7250.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="7250.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="7250.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="7250.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="7250.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="7250.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7250.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="7261" servers="1">
      <Property name="Capping" id="7261.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7261.Capping.bool"/>
      </Property>
      <Property name="Center" id="7261.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7261.Center.range"/>
      </Property>
      <Property name="Direction" id="7261.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7261.Direction.range"/>
      </Property>
      <Property name="Height" id="7261.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7261.Height.range"/>
      </Property>
      <Property name="Radius" id="7261.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7261.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7261.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7261.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="7272" servers="1">
      <Property name="Center" id="7272.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7272.Center.range"/>
      </Property>
      <Property name="XLength" id="7272.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7272.XLength.range"/>
      </Property>
      <Property name="YLength" id="7272.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7272.YLength.range"/>
      </Property>
      <Property name="ZLength" id="7272.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7272.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="7283" servers="1">
      <Property name="Capping" id="7283.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7283.Capping.bool"/>
      </Property>
      <Property name="Center" id="7283.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7283.Center.range"/>
      </Property>
      <Property name="Height" id="7283.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7283.Height.range"/>
      </Property>
      <Property name="Radius" id="7283.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7283.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7283.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7283.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="7294" servers="1">
      <Property name="Point1" id="7294.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7294.Point1.range"/>
      </Property>
      <Property name="Point2" id="7294.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7294.Point2.range"/>
      </Property>
      <Property name="Resolution" id="7294.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7294.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="7305" servers="1">
      <Property name="Center" id="7305.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7305.Center.range"/>
      </Property>
      <Property name="EndPhi" id="7305.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="7305.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="7305.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="7305.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="7305.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7305.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="7305.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7305.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="7305.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7305.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="7305.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7305.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="7305.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7305.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="7316" servers="1">
      <Property name="Center" id="7316.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7316.Center.range"/>
      </Property>
      <Property name="Filled" id="7316.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7316.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="7316.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="7316.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="7428" servers="16">
      <Property name="Source" id="7428.Source" number_of_elements="1">
        <Proxy value="7417"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform3" id="7416" servers="1">
      <Property name="Position" id="7416.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7416.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="7416.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="7416.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7416.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="7416.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="7416.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7416.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="7416.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="extended_sources" type="Transform2" id="7526" servers="1">
      <Property name="Position" id="7526.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7526.Position.range"/>
      </Property>
      <Property name="PositionInfo" id="7526.PositionInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Rotation" id="7526.Rotation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7526.Rotation.range"/>
      </Property>
      <Property name="RotationInfo" id="7526.RotationInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Scale" id="7526.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7526.Scale.range"/>
      </Property>
      <Property name="ScaleInfo" id="7526.ScaleInfo" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="misc" type="RepresentationAnimationHelper" id="7615" servers="16">
      <Property name="Source" id="7615.Source" number_of_elements="1">
        <Proxy value="7604"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ArrowSource" id="7527" servers="1">
      <Property name="Invert" id="7527.Invert" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7527.Invert.bool"/>
      </Property>
      <Property name="ShaftRadius" id="7527.ShaftRadius" number_of_elements="1">
        <Element index="0" value="0.03"/>
        <Domain name="range" id="7527.ShaftRadius.range"/>
      </Property>
      <Property name="ShaftResolution" id="7527.ShaftResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7527.ShaftResolution.range"/>
      </Property>
      <Property name="TipLength" id="7527.TipLength" number_of_elements="1">
        <Element index="0" value="0.35"/>
        <Domain name="range" id="7527.TipLength.range"/>
      </Property>
      <Property name="TipRadius" id="7527.TipRadius" number_of_elements="1">
        <Element index="0" value="0.1"/>
        <Domain name="range" id="7527.TipRadius.range"/>
      </Property>
      <Property name="TipResolution" id="7527.TipResolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7527.TipResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="ConeSource" id="7538" servers="1">
      <Property name="Capping" id="7538.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7538.Capping.bool"/>
      </Property>
      <Property name="Center" id="7538.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7538.Center.range"/>
      </Property>
      <Property name="Direction" id="7538.Direction" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7538.Direction.range"/>
      </Property>
      <Property name="Height" id="7538.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7538.Height.range"/>
      </Property>
      <Property name="Radius" id="7538.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7538.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7538.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7538.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CubeSource" id="7549" servers="1">
      <Property name="Center" id="7549.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7549.Center.range"/>
      </Property>
      <Property name="XLength" id="7549.XLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7549.XLength.range"/>
      </Property>
      <Property name="YLength" id="7549.YLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7549.YLength.range"/>
      </Property>
      <Property name="ZLength" id="7549.ZLength" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7549.ZLength.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="CylinderSource" id="7560" servers="1">
      <Property name="Capping" id="7560.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7560.Capping.bool"/>
      </Property>
      <Property name="Center" id="7560.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7560.Center.range"/>
      </Property>
      <Property name="Height" id="7560.Height" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7560.Height.range"/>
      </Property>
      <Property name="Radius" id="7560.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7560.Radius.range"/>
      </Property>
      <Property name="Resolution" id="7560.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7560.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="LineSource" id="7571" servers="1">
      <Property name="Point1" id="7571.Point1" number_of_elements="3">
        <Element index="0" value="-0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7571.Point1.range"/>
      </Property>
      <Property name="Point2" id="7571.Point2" number_of_elements="3">
        <Element index="0" value="0.5"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7571.Point2.range"/>
      </Property>
      <Property name="Resolution" id="7571.Resolution" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="7571.Resolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="SphereSource" id="7582" servers="1">
      <Property name="Center" id="7582.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7582.Center.range"/>
      </Property>
      <Property name="EndPhi" id="7582.EndPhi" number_of_elements="1">
        <Element index="0" value="180"/>
        <Domain name="range" id="7582.EndPhi.range"/>
      </Property>
      <Property name="EndTheta" id="7582.EndTheta" number_of_elements="1">
        <Element index="0" value="360"/>
        <Domain name="range" id="7582.EndTheta.range"/>
      </Property>
      <Property name="PhiResolution" id="7582.PhiResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7582.PhiResolution.range"/>
      </Property>
      <Property name="Radius" id="7582.Radius" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7582.Radius.range"/>
      </Property>
      <Property name="StartPhi" id="7582.StartPhi" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7582.StartPhi.range"/>
      </Property>
      <Property name="StartTheta" id="7582.StartTheta" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7582.StartTheta.range"/>
      </Property>
      <Property name="ThetaResolution" id="7582.ThetaResolution" number_of_elements="1">
        <Element index="0" value="8"/>
        <Domain name="range" id="7582.ThetaResolution.range"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="GlyphSource2D" id="7593" servers="1">
      <Property name="Center" id="7593.Center" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7593.Center.range"/>
      </Property>
      <Property name="Filled" id="7593.Filled" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7593.Filled.bool"/>
      </Property>
      <Property name="GlyphType" id="7593.GlyphType" number_of_elements="1">
        <Element index="0" value="9"/>
        <Domain name="enum" id="7593.GlyphType.enum">
          <Entry value="1" text="Vertex"/>
          <Entry value="2" text="Dash"/>
          <Entry value="3" text="Cross"/>
          <Entry value="4" text="ThickCross"/>
          <Entry value="5" text="Triangle"/>
          <Entry value="6" text="Square"/>
          <Entry value="7" text="Circle"/>
          <Entry value="8" text="Diamond"/>
          <Entry value="9" text="Arrow"/>
          <Entry value="10" text="ThickArrow"/>
          <Entry value="11" text="HookedArrow"/>
          <Entry value="12" text="EdgeArrow"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="5623" servers="21">
      <Property name="CubeAxesVisibility" id="5623.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="5623.Input" number_of_elements="1">
        <Proxy value="4787" output_port="0"/>
        <Domain name="input_array_cell" id="5623.Input.input_array_cell"/>
        <Domain name="input_array_point" id="5623.Input.input_array_point"/>
        <Domain name="input_type" id="5623.Input.input_type"/>
      </Property>
      <Property name="Representation" id="5623.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="5623.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="5623.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="5623.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Contact age [s]"/>
        <Domain name="array_list" id="5623.SelectionCellFieldDataArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="5623.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalPointIds"/>
        <Domain name="array_list" id="5623.SelectionPointFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionVisibility" id="5623.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="5623.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="5623.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5623.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="5623.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="5623.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="5623.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="5623.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="5623.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="5623.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="5623.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="5623.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="5623.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5623.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="5623.BlockOpacity"/>
      <Property name="BlockVisibility" id="5623.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="5623.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="5623.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value=""/>
        <Element index="4" value=""/>
        <Domain name="array_list" id="5623.ColorArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
        <Domain name="field_list" id="5623.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="5623.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="5623.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5623.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="5623.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5623.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="5623.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="5623.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="5623.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="5623.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="5623.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="5623.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="5623.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="5623.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="5623.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="5623.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="5623.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="5623.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="5623.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="5623.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="5623.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="5623.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="5623.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="5623.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="5623.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="5623.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="5623.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="5623.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="5623.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="5623.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="5623.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="5623.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="5623.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="5623.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="5623.DataBounds" number_of_elements="6">
        <Element index="0" value="6.10282454773138e-316"/>
        <Element index="1" value="0"/>
        <Element index="2" value="6.10282099045873e-316"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="6.10282296672131e-316"/>
        <Element index="5" value="4.44659081257122e-323"/>
      </Property>
      <Property name="Diffuse" id="5623.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="5623.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="5623.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="5623.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="5623.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="5623.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5623.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="5623.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="5623.LookupTable">
        <Domain name="groups" id="5623.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="5623.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="5623.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="5623.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="5623.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="5623.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.Opacity.range"/>
      </Property>
      <Property name="Orient" id="5623.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="5623.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="5623.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="5623.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="5623.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="5623.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5623.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="5623.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5623.PointSize.range"/>
      </Property>
      <Property name="Position" id="5623.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.Position.range"/>
      </Property>
      <Property name="Scale" id="5623.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="5623.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7432.5"/>
        <Domain name="bounds" id="5623.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="5623.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="5623.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="5623.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="5623.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="5623.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="5623.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="5623.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="5623.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="5623.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="5623.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="5623.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5623.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="5623.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="5623.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="5623.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="5623.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="5623.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="5623.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="5623.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="5623.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5623.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="5623.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="5623.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="5623.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5623.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="5623.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="5623.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5623.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="5623.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="5623.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="5623.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5623.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="5623.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5623.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="5623.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="5623.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="5623.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5623.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="5623.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5623.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="5623.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="5623.Source">
        <Domain name="input_type" id="5623.Source.input_type"/>
      </Property>
      <Property name="Specular" id="5623.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5623.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="5623.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5623.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="5623.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="5623.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="5623.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="5623.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="5623.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="5623.Texture">
        <Domain name="groups" id="5623.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="5623.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="5623.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5623.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="5623.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="7682" servers="21">
      <Property name="CubeAxesVisibility" id="7682.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="7682.Input" number_of_elements="1">
        <Proxy value="7604" output_port="0"/>
        <Domain name="input_array_cell" id="7682.Input.input_array_cell"/>
        <Domain name="input_array_point" id="7682.Input.input_array_point"/>
        <Domain name="input_type" id="7682.Input.input_type"/>
      </Property>
      <Property name="Representation" id="7682.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="7682.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="7682.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="7682.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="7682.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="7682.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="e: Relative interface height [m]"/>
        <Domain name="array_list" id="7682.SelectionPointFieldDataArrayName.array_list">
          <String text="GlyphVector"/>
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="7682.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="7682.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="7682.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7682.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="7682.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="7682.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="7682.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="7682.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="7682.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="7682.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="7682.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="7682.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="7682.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="7682.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="7682.BlockOpacity"/>
      <Property name="BlockVisibility" id="7682.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="7682.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="7682.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="7682.ColorArrayName.array_list">
          <String text="GlyphVector"/>
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="7682.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="7682.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="7682.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7682.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="7682.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7682.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="7682.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="7682.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="7682.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="7682.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="7682.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="7682.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="7682.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="7682.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="7682.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="7682.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="7682.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="7682.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="7682.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="7682.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="7682.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="7682.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="7682.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="7682.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="7682.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="7682.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="7682.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="7682.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="7682.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="7682.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="7682.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="7682.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="7682.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="7682.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="7682.DataBounds" number_of_elements="6">
        <Element index="0" value="1.92137910742298e-314"/>
        <Element index="1" value="2.10077583305562e-312"/>
        <Element index="2" value="1.43279037293961e-322"/>
        <Element index="3" value="1.3570418091293e-316"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1.37113018983359e-316"/>
      </Property>
      <Property name="Diffuse" id="7682.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="7682.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="7682.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="7682.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="7682.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="7682.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7682.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="7682.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="7682.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="7682.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="7682.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="7682.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="7682.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="7682.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="7682.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.Opacity.range"/>
      </Property>
      <Property name="Orient" id="7682.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="7682.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="7682.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="7682.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="7682.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="7682.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7682.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="7682.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7682.PointSize.range"/>
      </Property>
      <Property name="Position" id="7682.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.Position.range"/>
      </Property>
      <Property name="Scale" id="7682.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="7682.ScaleFactor" number_of_elements="1">
        <Element index="0" value="8073.48125"/>
        <Domain name="bounds" id="7682.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7682.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7682.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7682.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="7682.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="7682.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="7682.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="7682.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="7682.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="7682.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="7682.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="7682.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7682.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="7682.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="7682.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="7682.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="7682.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="7682.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="7682.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="7682.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="7682.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7682.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="7682.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="7682.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="7682.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7682.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="7682.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="7682.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7682.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="7682.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="7682.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="7682.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7682.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="7682.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7682.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="7682.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="7682.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="7682.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="7682.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="7682.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7682.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="7682.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="7682.Source">
        <Domain name="input_type" id="7682.Source.input_type"/>
      </Property>
      <Property name="Specular" id="7682.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7682.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="7682.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7682.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="7682.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="7682.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="7682.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="7682.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="7682.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="7682.Texture">
        <Domain name="groups" id="7682.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="7682.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="7682.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7682.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="7682.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="7238" servers="21">
      <Property name="CubeAxesVisibility" id="7238.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="7238.Input" number_of_elements="1">
        <Proxy value="7160" output_port="0"/>
        <Domain name="input_array_cell" id="7238.Input.input_array_cell"/>
        <Domain name="input_array_point" id="7238.Input.input_array_point"/>
        <Domain name="input_type" id="7238.Input.input_type"/>
      </Property>
      <Property name="Representation" id="7238.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="7238.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="7238.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="7238.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="7238.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="7238.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Diameter (areal) [m]"/>
        <Domain name="array_list" id="7238.SelectionPointFieldDataArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="GlyphVector"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Normals"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="7238.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="7238.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="7238.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7238.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="7238.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="7238.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="7238.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="7238.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="7238.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="7238.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="7238.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="7238.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="7238.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="7238.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="7238.BlockOpacity"/>
      <Property name="BlockVisibility" id="7238.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="7238.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="7238.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Linear velocity [m s^-1]"/>
        <Domain name="array_list" id="7238.ColorArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="GlyphVector"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Normals"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="7238.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="7238.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="7238.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7238.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="7238.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7238.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="7238.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="7238.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="7238.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="7238.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="7238.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="7238.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="7238.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="7238.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="7238.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="7238.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="7238.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="7238.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="7238.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="7238.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="7238.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="7238.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="7238.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="7238.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="7238.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="7238.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="7238.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="7238.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="7238.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="7238.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="7238.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="7238.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="7238.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="7238.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="7238.DataBounds" number_of_elements="6">
        <Element index="0" value="3.95252516672997e-322"/>
        <Element index="1" value="1.36291059754737e-315"/>
        <Element index="2" value="1.36283474858942e-315"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="1.36283447191265e-315"/>
        <Element index="5" value="4.64476894223412e-317"/>
      </Property>
      <Property name="Diffuse" id="7238.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="7238.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="7238.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="7238.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="7238.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="7238.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7238.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="7238.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="7238.LookupTable" number_of_elements="1">
        <Proxy value="5740"/>
        <Domain name="groups" id="7238.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="7238.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="7238.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="7238.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="7238.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="7238.Opacity" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="7238.Opacity.range"/>
      </Property>
      <Property name="Orient" id="7238.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="7238.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="7238.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="7238.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="7238.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="7238.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7238.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="7238.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7238.PointSize.range"/>
      </Property>
      <Property name="Position" id="7238.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.Position.range"/>
      </Property>
      <Property name="Scale" id="7238.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="7238.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7498.30772294998"/>
        <Domain name="bounds" id="7238.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7238.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7238.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7238.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="7238.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="7238.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="7238.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="7238.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="7238.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="7238.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="7238.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="7238.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7238.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="7238.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="7238.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="7238.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="7238.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="7238.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="7238.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="7238.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="7238.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7238.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="7238.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="7238.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="7238.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7238.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="7238.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="7238.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7238.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="7238.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="7238.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="7238.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7238.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="7238.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7238.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="7238.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="7238.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="7238.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="7238.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="7238.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7238.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="7238.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="7238.Source">
        <Domain name="input_type" id="7238.Source.input_type"/>
      </Property>
      <Property name="Specular" id="7238.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7238.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="7238.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7238.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="7238.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="7238.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="7238.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="7238.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="7238.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="7238.Texture">
        <Domain name="groups" id="7238.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="7238.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="7238.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7238.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="7238.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="7405" servers="21">
      <Property name="CubeAxesVisibility" id="7405.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="7405.Input" number_of_elements="1">
        <Proxy value="7327" output_port="0"/>
        <Domain name="input_array_cell" id="7405.Input.input_array_cell"/>
        <Domain name="input_array_point" id="7405.Input.input_array_point"/>
        <Domain name="input_type" id="7405.Input.input_type"/>
      </Property>
      <Property name="Representation" id="7405.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="7405.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="7405.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="7405.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="7405.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="7405.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Atmosphere drag coefficient (horizontal) [-]"/>
        <Domain name="array_list" id="7405.SelectionPointFieldDataArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="GlyphVector"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="7405.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="7405.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="7405.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7405.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="7405.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="7405.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="7405.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="7405.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="7405.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="7405.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="7405.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="7405.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="7405.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="7405.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="7405.BlockOpacity"/>
      <Property name="BlockVisibility" id="7405.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="7405.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="7405.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Linear velocity [m s^-1]"/>
        <Domain name="array_list" id="7405.ColorArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="GlyphVector"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="7405.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="7405.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="7405.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7405.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="7405.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7405.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="7405.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="7405.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="7405.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="7405.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="7405.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="7405.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="7405.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="7405.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="7405.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="7405.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="7405.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="7405.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="7405.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="7405.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="7405.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="7405.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="7405.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="7405.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="7405.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="7405.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="7405.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="7405.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="7405.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="7405.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="7405.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="7405.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="7405.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="7405.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="7405.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1.19488140101286e-315"/>
        <Element index="3" value="1.19488155911387e-315"/>
        <Element index="4" value="3.49897290384771e-320"/>
        <Element index="5" value="5.73829384359551e+194"/>
      </Property>
      <Property name="Diffuse" id="7405.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="7405.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="7405.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="7405.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="7405.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="7405.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7405.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="7405.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="7405.LookupTable" number_of_elements="1">
        <Proxy value="5740"/>
        <Domain name="groups" id="7405.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="7405.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="7405.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="7405.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="7405.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="7405.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.Opacity.range"/>
      </Property>
      <Property name="Orient" id="7405.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="7405.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="7405.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="7405.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="7405.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="7405.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7405.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="7405.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7405.PointSize.range"/>
      </Property>
      <Property name="Position" id="7405.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.Position.range"/>
      </Property>
      <Property name="Scale" id="7405.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="7405.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7432.5"/>
        <Domain name="bounds" id="7405.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7405.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7405.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7405.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="7405.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="7405.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="7405.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="7405.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="7405.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="7405.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="7405.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="7405.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7405.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="7405.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="7405.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="7405.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="7405.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="7405.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="7405.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="7405.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="7405.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7405.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="7405.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="7405.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="7405.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7405.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="7405.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="7405.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7405.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="7405.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="7405.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="7405.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7405.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="7405.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7405.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="7405.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="7405.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="7405.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="7405.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="7405.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7405.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="7405.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="7405.Source">
        <Domain name="input_type" id="7405.Source.input_type"/>
      </Property>
      <Property name="Specular" id="7405.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7405.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="7405.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7405.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="7405.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="7405.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="7405.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="7405.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="7405.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="7405.Texture">
        <Domain name="groups" id="7405.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="7405.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="7405.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7405.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="7405.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="6783" servers="21">
      <Property name="CubeAxesVisibility" id="6783.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="6783.Input" number_of_elements="1">
        <Proxy value="6701" output_port="0"/>
        <Domain name="input_array_cell" id="6783.Input.input_array_cell"/>
        <Domain name="input_array_point" id="6783.Input.input_array_point"/>
        <Domain name="input_type" id="6783.Input.input_type"/>
      </Property>
      <Property name="Representation" id="6783.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="6783.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="6783.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="6783.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Contact age [s]"/>
        <Domain name="array_list" id="6783.SelectionCellFieldDataArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="6783.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalPointIds"/>
        <Domain name="array_list" id="6783.SelectionPointFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionVisibility" id="6783.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="6783.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="6783.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6783.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="6783.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="6783.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="6783.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="6783.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="6783.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="6783.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="6783.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="6783.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="6783.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="6783.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="6783.BlockOpacity"/>
      <Property name="BlockVisibility" id="6783.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="6783.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="6783.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value=""/>
        <Element index="4" value=""/>
        <Domain name="array_list" id="6783.ColorArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
        <Domain name="field_list" id="6783.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="6783.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="6783.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6783.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="6783.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6783.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="6783.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="6783.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="6783.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="6783.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="6783.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="6783.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="6783.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="6783.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="6783.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="6783.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="6783.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="6783.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="6783.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="6783.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="6783.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="6783.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="6783.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="6783.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="6783.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="6783.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="6783.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="6783.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="6783.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="6783.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="6783.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="6783.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="6783.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="6783.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="6783.DataBounds" number_of_elements="6">
        <Element index="0" value="0.0078125"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1.3487992131466e-321"/>
        <Element index="3" value="1.15850394038835e-315"/>
        <Element index="4" value="6.92299120110413e-310"/>
        <Element index="5" value="0.0078125"/>
      </Property>
      <Property name="Diffuse" id="6783.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="6783.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="6783.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="6783.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="6783.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="6783.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6783.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="6783.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="6783.LookupTable">
        <Domain name="groups" id="6783.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="6783.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="6783.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="6783.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="6783.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="6783.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.Opacity.range"/>
      </Property>
      <Property name="Orient" id="6783.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="6783.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="6783.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="6783.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="6783.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="6783.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6783.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="6783.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6783.PointSize.range"/>
      </Property>
      <Property name="Position" id="6783.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.Position.range"/>
      </Property>
      <Property name="Scale" id="6783.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="6783.ScaleFactor" number_of_elements="1">
        <Element index="0" value="3650.19360351562"/>
        <Domain name="bounds" id="6783.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6783.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6783.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6783.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="6783.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="6783.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="6783.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="6783.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="6783.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="6783.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="6783.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="6783.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6783.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="6783.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="6783.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="6783.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="6783.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="6783.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="6783.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="6783.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="6783.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6783.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="6783.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="6783.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="6783.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6783.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="6783.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="6783.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6783.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="6783.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="6783.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="6783.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6783.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="6783.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6783.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="6783.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="6783.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="6783.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6783.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="6783.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6783.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="6783.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="6783.Source">
        <Domain name="input_type" id="6783.Source.input_type"/>
      </Property>
      <Property name="Specular" id="6783.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6783.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="6783.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6783.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="6783.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="6783.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="6783.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="6783.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="6783.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="6783.Texture">
        <Domain name="groups" id="6783.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="6783.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="6783.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6783.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="6783.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="6872" servers="21">
      <Property name="CubeAxesVisibility" id="6872.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="6872.Input" number_of_elements="1">
        <Proxy value="6794" output_port="0"/>
        <Domain name="input_array_cell" id="6872.Input.input_array_cell"/>
        <Domain name="input_array_point" id="6872.Input.input_array_point"/>
        <Domain name="input_type" id="6872.Input.input_type"/>
      </Property>
      <Property name="Representation" id="6872.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="6872.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="6872.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="6872.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="6872.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="6872.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Contact age [s]"/>
        <Domain name="array_list" id="6872.SelectionPointFieldDataArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="6872.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="6872.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="6872.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6872.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="6872.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="6872.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="6872.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="6872.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="6872.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="6872.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="6872.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="6872.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="6872.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="6872.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="6872.BlockOpacity"/>
      <Property name="BlockVisibility" id="6872.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="6872.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="6872.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value=""/>
        <Element index="4" value=""/>
        <Domain name="array_list" id="6872.ColorArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
        <Domain name="field_list" id="6872.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="6872.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="6872.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6872.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="6872.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6872.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="6872.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="6872.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="6872.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="6872.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="6872.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="6872.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="6872.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="6872.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="6872.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="6872.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="6872.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="6872.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="6872.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="6872.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="6872.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="6872.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="6872.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="6872.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="6872.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="6872.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="6872.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="6872.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="6872.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="6872.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="6872.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="6872.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="6872.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="6872.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="6872.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="3.33948851337015e-319"/>
        <Element index="2" value="1.17167584908223e-315"/>
        <Element index="3" value="8.29005889303207e-317"/>
        <Element index="4" value="7.95445689804407e-322"/>
        <Element index="5" value="1.15144283322849e-315"/>
      </Property>
      <Property name="Diffuse" id="6872.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="6872.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="6872.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="6872.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="6872.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="6872.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6872.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="6872.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="6872.LookupTable">
        <Domain name="groups" id="6872.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="6872.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="6872.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="6872.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="6872.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="6872.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.Opacity.range"/>
      </Property>
      <Property name="Orient" id="6872.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="6872.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="6872.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="6872.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="6872.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="6872.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6872.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="6872.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6872.PointSize.range"/>
      </Property>
      <Property name="Position" id="6872.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.Position.range"/>
      </Property>
      <Property name="Scale" id="6872.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="6872.ScaleFactor" number_of_elements="1">
        <Element index="0" value="3650.19360351562"/>
        <Domain name="bounds" id="6872.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6872.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6872.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6872.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="6872.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="6872.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="6872.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="6872.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="6872.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="6872.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="6872.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="6872.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6872.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="6872.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="6872.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="6872.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="6872.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="6872.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="6872.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="6872.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="6872.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6872.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="6872.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="6872.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="6872.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6872.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="6872.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="6872.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6872.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="6872.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="6872.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="6872.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6872.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="6872.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6872.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="6872.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="6872.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="6872.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6872.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="6872.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6872.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="6872.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="6872.Source">
        <Domain name="input_type" id="6872.Source.input_type"/>
      </Property>
      <Property name="Specular" id="6872.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6872.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="6872.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6872.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="6872.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="6872.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="6872.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="6872.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="6872.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="6872.Texture">
        <Domain name="groups" id="6872.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="6872.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="6872.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6872.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="6872.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="6961" servers="21">
      <Property name="CubeAxesVisibility" id="6961.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="6961.Input" number_of_elements="1">
        <Proxy value="6883" output_port="0"/>
        <Domain name="input_array_cell" id="6961.Input.input_array_cell"/>
        <Domain name="input_array_point" id="6961.Input.input_array_point"/>
        <Domain name="input_type" id="6961.Input.input_type"/>
      </Property>
      <Property name="Representation" id="6961.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="6961.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="6961.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="6961.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="6961.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="6961.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Contact age [s]"/>
        <Domain name="array_list" id="6961.SelectionPointFieldDataArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
          <String text="TubeNormals"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="6961.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="6961.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="6961.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6961.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="6961.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="6961.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="6961.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="6961.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="6961.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="6961.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="6961.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="6961.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="6961.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="6961.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="6961.BlockOpacity"/>
      <Property name="BlockVisibility" id="6961.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="6961.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="6961.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Tensile stress [Pa]"/>
        <Domain name="array_list" id="6961.ColorArrayName.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Tensile stress [Pa]"/>
          <String text="TubeNormals"/>
        </Domain>
        <Domain name="field_list" id="6961.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="6961.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="6961.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6961.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="6961.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6961.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="6961.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="6961.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="6961.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="6961.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="6961.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="6961.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="6961.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="6961.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="6961.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="6961.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="6961.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="6961.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="6961.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="6961.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="6961.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="6961.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="6961.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="6961.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="6961.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="6961.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="6961.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="6961.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="6961.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="6961.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="6961.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="6961.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="6961.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="6961.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="6961.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="8.43208122494927e-316"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="0"/>
        <Element index="5" value="0"/>
      </Property>
      <Property name="Diffuse" id="6961.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="6961.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="6961.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="6961.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="6961.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="6961.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6961.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="6961.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="6961.LookupTable" number_of_elements="1">
        <Proxy value="6072"/>
        <Domain name="groups" id="6961.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="6961.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="6961.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="6961.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="6961.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="6961.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.Opacity.range"/>
      </Property>
      <Property name="Orient" id="6961.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="6961.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="6961.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="6961.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="6961.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="6961.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6961.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="6961.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6961.PointSize.range"/>
      </Property>
      <Property name="Position" id="6961.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.Position.range"/>
      </Property>
      <Property name="Scale" id="6961.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="6961.ScaleFactor" number_of_elements="1">
        <Element index="0" value="3685.41860351563"/>
        <Domain name="bounds" id="6961.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6961.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6961.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6961.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="6961.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="6961.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="6961.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="6961.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="6961.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="6961.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="6961.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="6961.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6961.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="6961.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="6961.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="6961.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="6961.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="6961.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="6961.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="6961.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="6961.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6961.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="6961.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="6961.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="6961.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6961.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="6961.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="6961.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6961.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="6961.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="6961.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="6961.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6961.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="6961.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6961.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="6961.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="6961.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="6961.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6961.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="6961.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6961.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="6961.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="6961.Source">
        <Domain name="input_type" id="6961.Source.input_type"/>
      </Property>
      <Property name="Specular" id="6961.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6961.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="6961.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6961.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="6961.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="6961.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="6961.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="6961.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="6961.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="6961.Texture">
        <Domain name="groups" id="6961.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="6961.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="6961.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6961.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="6961.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="GeometryRepresentation" id="6681" servers="21">
      <Property name="CubeAxesVisibility" id="6681.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="6681.Input" number_of_elements="1">
        <Proxy value="6603" output_port="0"/>
        <Domain name="input_array_cell" id="6681.Input.input_array_cell"/>
        <Domain name="input_array_point" id="6681.Input.input_array_point"/>
        <Domain name="input_type" id="6681.Input.input_type"/>
      </Property>
      <Property name="Representation" id="6681.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="6681.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="6681.RepresentationTypesInfo" number_of_elements="6">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="6681.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="6681.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="6681.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
        <Domain name="array_list" id="6681.SelectionPointFieldDataArrayName.array_list">
          <String text="GlyphVector"/>
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="6681.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="6681.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="6681.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6681.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="6681.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="6681.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="6681.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="6681.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="6681.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="6681.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="6681.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="6681.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="6681.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="6681.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="6681.BlockOpacity"/>
      <Property name="BlockVisibility" id="6681.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="6681.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="6681.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="6681.ColorArrayName.array_list">
          <String text="GlyphVector"/>
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="6681.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="6681.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="6681.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6681.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="6681.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6681.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="6681.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="6681.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="6681.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="6681.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="6681.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="6681.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="6681.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="6681.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="6681.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="6681.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="6681.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="6681.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="6681.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="6681.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="6681.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="6681.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="6681.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="6681.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="6681.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="6681.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="6681.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="6681.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="6681.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="6681.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="6681.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="6681.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="6681.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="6681.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="6681.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1.33511968839482e-306"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="4.00758928069504e-270"/>
      </Property>
      <Property name="Diffuse" id="6681.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="6681.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="6681.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="6681.EdgeColor.range"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="6681.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="6681.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6681.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="6681.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="6681.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="6681.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="6681.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="6681.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="6681.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="6681.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="6681.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.Opacity.range"/>
      </Property>
      <Property name="Orient" id="6681.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="6681.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="6681.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="6681.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="6681.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="6681.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6681.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="6681.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6681.PointSize.range"/>
      </Property>
      <Property name="Position" id="6681.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.Position.range"/>
      </Property>
      <Property name="Scale" id="6681.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="6681.ScaleFactor" number_of_elements="1">
        <Element index="0" value="30000"/>
        <Domain name="bounds" id="6681.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6681.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6681.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6681.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="6681.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.Scaling.bool"/>
      </Property>
      <Property name="SelectMaskArray" id="6681.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="6681.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="6681.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="6681.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="6681.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="6681.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="6681.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6681.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="6681.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="6681.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="6681.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="6681.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="6681.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="6681.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="6681.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="6681.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6681.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="6681.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="6681.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="6681.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6681.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="6681.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="6681.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6681.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="6681.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="6681.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="6681.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6681.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="6681.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6681.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="6681.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="6681.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="6681.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6681.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="6681.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6681.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="6681.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="6681.Source">
        <Domain name="input_type" id="6681.Source.input_type"/>
      </Property>
      <Property name="Specular" id="6681.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6681.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="6681.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6681.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="6681.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="6681.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="6681.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="6681.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="6681.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="6681.Texture">
        <Domain name="groups" id="6681.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="6681.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="6681.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6681.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="6681.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="StructuredGridRepresentation" id="5360" servers="21">
      <Property name="CubeAxesVisibility" id="5360.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="5360.Input" number_of_elements="1">
        <Proxy value="4765" output_port="0"/>
        <Domain name="input_array_point" id="5360.Input.input_array_point"/>
        <Domain name="input_type" id="5360.Input.input_type"/>
      </Property>
      <Property name="Representation" id="5360.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="5360.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="5360.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="5360.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5360.SelectionCellFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="5360.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5360.SelectionPointFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="5360.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="5360.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="5360.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5360.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="5360.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="5360.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="5360.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="5360.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="5360.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="5360.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="5360.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="5360.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="5360.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5360.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="5360.BlockOpacity"/>
      <Property name="BlockVisibility" id="5360.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="5360.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="5360.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5360.ColorArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
        <Domain name="field_list" id="5360.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="5360.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="5360.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5360.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="5360.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5360.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="5360.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="5360.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="5360.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="5360.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="5360.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="5360.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="5360.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="5360.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="5360.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="5360.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="5360.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="5360.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="5360.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="5360.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="5360.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="5360.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="5360.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="5360.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="5360.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="5360.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="5360.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="5360.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="5360.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="5360.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="5360.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="5360.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="5360.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="5360.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="5360.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="5.21444076216655e-316"/>
        <Element index="2" value="0"/>
        <Element index="3" value="8.48798316336702e-314"/>
        <Element index="4" value="5.21444945772192e-316"/>
        <Element index="5" value="8.48798316386109e-314"/>
      </Property>
      <Property name="Diffuse" id="5360.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="5360.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="5360.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="5360.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="5360.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="5360.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="5360.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="5360.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5360.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="5360.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="5360.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="5360.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="5360.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="5360.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="5360.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="5360.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="5360.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.Opacity.range"/>
      </Property>
      <Property name="Orient" id="5360.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="5360.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="5360.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="5360.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="5360.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="5360.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="5360.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5360.PointSize.range"/>
      </Property>
      <Property name="Position" id="5360.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5360.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5298"/>
        <Domain name="groups" id="5360.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="5360.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="17359.7055400195"/>
        <Domain name="bounds" id="5360.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="5360.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="5360.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7500"/>
        <Domain name="bounds" id="5360.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="5360.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="5360.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="5360.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="5360.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="5360.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="5360.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="5360.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="5360.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="5360.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="5360.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="5360.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="5360.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="5360.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5360.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="5360.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="5360.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="5360.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="5360.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="5360.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="5360.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="5360.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="5360.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5360.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="5360.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="5360.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="5360.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5360.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="5360.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="5360.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5360.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="5360.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="5360.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="5360.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5360.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="5360.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5360.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="5360.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="5360.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="5360.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5360.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="5360.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5360.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="5360.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="5360.Source">
        <Domain name="input_type" id="5360.Source.input_type"/>
      </Property>
      <Property name="Specular" id="5360.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5360.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="5360.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5360.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="5360.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="5360.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="5360.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="5360.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="5360.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="5360.Texture">
        <Domain name="groups" id="5360.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="5360.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="5360.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5360.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UseDataParititions" id="5360.UseDataParititions" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5360.UseDataParititions.bool"/>
      </Property>
      <Property name="UserTransform" id="5360.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="StructuredGridRepresentation" id="5546" servers="21">
      <Property name="CubeAxesVisibility" id="5546.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="5546.Input" number_of_elements="1">
        <Proxy value="4831" output_port="0"/>
        <Domain name="input_array_point" id="5546.Input.input_array_point"/>
        <Domain name="input_type" id="5546.Input.input_type"/>
      </Property>
      <Property name="Representation" id="5546.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="5546.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="5546.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="5546.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5546.SelectionCellFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="5546.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5546.SelectionPointFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="5546.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="5546.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="5546.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5546.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="5546.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="5546.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="5546.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="5546.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="5546.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="5546.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="5546.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="5546.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="5546.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5546.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="5546.BlockOpacity"/>
      <Property name="BlockVisibility" id="5546.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="5546.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="5546.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="5546.ColorArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
        <Domain name="field_list" id="5546.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="5546.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="5546.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5546.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="5546.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5546.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="5546.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="5546.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="5546.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="5546.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="5546.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="5546.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="5546.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="5546.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="5546.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="5546.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="5546.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="5546.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="5546.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="5546.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="5546.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="5546.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="5546.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="5546.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="5546.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="5546.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="5546.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="5546.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="5546.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="5546.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="5546.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="5546.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="5546.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="5546.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="5546.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="6.10341782175891e-316"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="4.94065645841247e-324"/>
        <Element index="5" value="1.97626258336499e-321"/>
      </Property>
      <Property name="Diffuse" id="5546.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="5546.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="5546.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="5546.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="5546.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="5546.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="5546.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="5546.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5546.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="5546.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="5546.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="5546.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="5546.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="5546.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="5546.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="5546.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="5546.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.Opacity.range"/>
      </Property>
      <Property name="Orient" id="5546.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="5546.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="5546.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="5546.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="5546.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="5546.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="5546.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5546.PointSize.range"/>
      </Property>
      <Property name="Position" id="5546.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5546.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5298"/>
        <Domain name="groups" id="5546.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="5546.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="17359.7055400195"/>
        <Domain name="bounds" id="5546.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="5546.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="5546.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7500"/>
        <Domain name="bounds" id="5546.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="5546.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="5546.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="5546.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="5546.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="5546.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="5546.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="5546.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="5546.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="5546.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="5546.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="5546.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="5546.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="5546.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5546.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="5546.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="5546.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="5546.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="5546.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="5546.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="5546.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="5546.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="5546.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5546.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="5546.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="5546.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="5546.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5546.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="5546.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="5546.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5546.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="5546.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="5546.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="5546.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5546.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="5546.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5546.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="5546.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="5546.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="5546.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5546.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="5546.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5546.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="5546.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="5546.Source">
        <Domain name="input_type" id="5546.Source.input_type"/>
      </Property>
      <Property name="Specular" id="5546.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5546.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="5546.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5546.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="5546.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="5546.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="5546.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="5546.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="5546.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="5546.Texture">
        <Domain name="groups" id="5546.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="5546.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="5546.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5546.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UseDataParititions" id="5546.UseDataParititions" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5546.UseDataParititions.bool"/>
      </Property>
      <Property name="UserTransform" id="5546.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="StructuredGridRepresentation" id="6424" servers="21">
      <Property name="CubeAxesVisibility" id="6424.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="6424.Input" number_of_elements="1">
        <Proxy value="6326" output_port="0"/>
        <Domain name="input_array_point" id="6424.Input.input_array_point"/>
        <Domain name="input_type" id="6424.Input.input_type"/>
      </Property>
      <Property name="Representation" id="6424.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="6424.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="6424.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="6424.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="6424.SelectionCellFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="6424.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="6424.SelectionPointFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="6424.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="6424.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="6424.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6424.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="6424.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="6424.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="6424.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="6424.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="6424.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="6424.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="6424.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="6424.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="6424.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="6424.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="6424.BlockOpacity"/>
      <Property name="BlockVisibility" id="6424.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="6424.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="6424.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="6424.ColorArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="6424.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="6424.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="6424.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6424.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="6424.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6424.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="6424.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="6424.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="6424.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="6424.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="6424.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="6424.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="6424.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="6424.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="6424.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="6424.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="6424.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="6424.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="6424.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="6424.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="6424.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="6424.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="6424.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="6424.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="6424.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="6424.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="6424.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="6424.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="6424.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="6424.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="6424.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="6424.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="6424.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="6424.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="6424.DataBounds" number_of_elements="6">
        <Element index="0" value="-5.4861292803319e+303"/>
        <Element index="1" value="-5.4861292803319e+303"/>
        <Element index="2" value="-5.4861292803319e+303"/>
        <Element index="3" value="-5.4861292803319e+303"/>
        <Element index="4" value="-5.4861292803319e+303"/>
        <Element index="5" value="-5.4861292803319e+303"/>
      </Property>
      <Property name="Diffuse" id="6424.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="6424.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="6424.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="6424.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="6424.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="6424.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="6424.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="6424.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6424.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="6424.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="6424.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="6424.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="6424.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="6424.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="6424.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="6424.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="6424.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.Opacity.range"/>
      </Property>
      <Property name="Orient" id="6424.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="6424.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="6424.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="6424.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="6424.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="6424.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="6424.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6424.PointSize.range"/>
      </Property>
      <Property name="Position" id="6424.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="6424.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5298"/>
        <Domain name="groups" id="6424.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="6424.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="17359.7055400195"/>
        <Domain name="bounds" id="6424.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="6424.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="6424.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7500"/>
        <Domain name="bounds" id="6424.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6424.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6424.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6424.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="6424.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="6424.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="6424.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="6424.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="6424.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="6424.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="6424.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="6424.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="6424.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="6424.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6424.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="6424.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="6424.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="6424.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="6424.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="6424.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="6424.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="6424.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="6424.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6424.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="6424.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="6424.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="6424.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6424.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="6424.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="6424.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="6424.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="6424.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="6424.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="6424.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6424.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="6424.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6424.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="6424.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="6424.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="6424.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6424.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="6424.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6424.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="6424.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="6424.Source">
        <Domain name="input_type" id="6424.Source.input_type"/>
      </Property>
      <Property name="Specular" id="6424.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6424.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="6424.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6424.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="6424.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="6424.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="6424.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="6424.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="6424.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="6424.Texture">
        <Domain name="groups" id="6424.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="6424.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="6424.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6424.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UseDataParititions" id="6424.UseDataParititions" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6424.UseDataParititions.bool"/>
      </Property>
      <Property name="UserTransform" id="6424.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="StructuredGridRepresentation" id="7515" servers="21">
      <Property name="CubeAxesVisibility" id="7515.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="7515.Input" number_of_elements="1">
        <Proxy value="7417" output_port="0"/>
        <Domain name="input_array_point" id="7515.Input.input_array_point"/>
        <Domain name="input_type" id="7515.Input.input_type"/>
      </Property>
      <Property name="Representation" id="7515.Representation" number_of_elements="1">
        <Element index="0" value="Surface With Edges"/>
        <Domain name="list" id="7515.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="7515.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="7515.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="7515.SelectionCellFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="7515.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="7515.SelectionPointFieldDataArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="7515.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="7515.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="7515.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7515.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="7515.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="7515.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="7515.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="7515.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="7515.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="7515.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="7515.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="7515.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="7515.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="7515.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="7515.BlockOpacity"/>
      <Property name="BlockVisibility" id="7515.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="7515.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="7515.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="7515.ColorArrayName.array_list">
          <String text="Velocity vector [m/s]"/>
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="cellNormals"/>
        </Domain>
        <Domain name="field_list" id="7515.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="7515.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="7515.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7515.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="7515.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7515.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="7515.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="7515.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="7515.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="7515.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="7515.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="7515.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="7515.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="7515.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="7515.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="7515.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="7515.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="7515.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="7515.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="7515.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="7515.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="7515.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="7515.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="7515.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="7515.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="7515.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="7515.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="7515.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="7515.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="7515.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="7515.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="7515.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="7515.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="7515.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="7515.DataBounds" number_of_elements="6">
        <Element index="0" value="1.35587423319505e-318"/>
        <Element index="1" value="1.35595531042681e-311"/>
        <Element index="2" value="1.39106676629983e-315"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="6.36598737289582e-314"/>
        <Element index="5" value="1.37113018983359e-316"/>
      </Property>
      <Property name="Diffuse" id="7515.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="7515.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="7515.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="7515.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="7515.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="7515.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="7515.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="7515.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7515.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="7515.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="7515.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="7515.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="7515.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="7515.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="7515.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="7515.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="7515.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.Opacity.range"/>
      </Property>
      <Property name="Orient" id="7515.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="7515.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="7515.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="7515.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="7515.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="7515.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="7515.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7515.PointSize.range"/>
      </Property>
      <Property name="Position" id="7515.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="7515.ScalarOpacityFunction" number_of_elements="1">
        <Proxy value="5298"/>
        <Domain name="groups" id="7515.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="7515.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="17359.7055400195"/>
        <Domain name="bounds" id="7515.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="7515.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="7515.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7500"/>
        <Domain name="bounds" id="7515.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7515.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7515.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7515.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="7515.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="7515.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="7515.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="7515.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="7515.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="7515.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="7515.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="7515.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="7515.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="7515.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7515.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="7515.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="7515.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="7515.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="7515.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="7515.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="7515.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="7515.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="7515.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7515.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="7515.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="7515.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="7515.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7515.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="7515.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="7515.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7515.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="7515.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="7515.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="7515.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7515.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="7515.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7515.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="7515.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="7515.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="7515.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="7515.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="7515.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7515.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="7515.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="7515.Source">
        <Domain name="input_type" id="7515.Source.input_type"/>
      </Property>
      <Property name="Specular" id="7515.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7515.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="7515.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7515.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="7515.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="7515.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="7515.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="7515.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="7515.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="7515.Texture">
        <Domain name="groups" id="7515.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="7515.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="7515.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7515.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UseDataParititions" id="7515.UseDataParititions" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7515.UseDataParititions.bool"/>
      </Property>
      <Property name="UserTransform" id="7515.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="UnstructuredGridRepresentation" id="5453" servers="21">
      <Property name="CubeAxesVisibility" id="5453.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="5453.Input" number_of_elements="1">
        <Proxy value="4809" output_port="0"/>
        <Domain name="input_array_cell" id="5453.Input.input_array_cell"/>
        <Domain name="input_array_point" id="5453.Input.input_array_point"/>
        <Domain name="input_type" id="5453.Input.input_type"/>
      </Property>
      <Property name="Representation" id="5453.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="5453.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="5453.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="5453.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="5453.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="5453.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Angular acceleration [rad s^-2]"/>
        <Domain name="array_list" id="5453.SelectionPointFieldDataArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="5453.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="5453.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="5453.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5453.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="5453.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="5453.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="5453.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="5453.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="5453.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="5453.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="5453.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="5453.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="5453.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="5453.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="5453.BlockOpacity"/>
      <Property name="BlockVisibility" id="5453.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="5453.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="5453.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value=""/>
        <Element index="4" value=""/>
        <Domain name="array_list" id="5453.ColorArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
        <Domain name="field_list" id="5453.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="5453.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="5453.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5453.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="5453.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5453.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="5453.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="5453.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="5453.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="5453.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="5453.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="5453.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="5453.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="5453.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="5453.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="5453.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="5453.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="5453.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="5453.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="5453.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="5453.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="5453.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="5453.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="5453.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="5453.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="5453.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="5453.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="5453.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="5453.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="5453.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="5453.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="5453.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="5453.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="5453.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="5453.DataBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="5.75407783742264e-316"/>
        <Element index="3" value="8.48798316386109e-314"/>
        <Element index="4" value="4.94065645841247e-324"/>
        <Element index="5" value="1.97626258336499e-321"/>
      </Property>
      <Property name="Diffuse" id="5453.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="5453.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="5453.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="5453.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="5453.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="5453.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="5453.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="5453.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5453.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="5453.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="5453.LookupTable">
        <Domain name="groups" id="5453.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="5453.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="5453.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="5453.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="5453.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="5453.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.Opacity.range"/>
      </Property>
      <Property name="Orient" id="5453.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="5453.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="5453.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="5453.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="5453.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="5453.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5453.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="5453.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5453.PointSize.range"/>
      </Property>
      <Property name="Position" id="5453.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="5453.ScalarOpacityFunction">
        <Domain name="groups" id="5453.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="5453.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="89202.9217570815"/>
        <Domain name="bounds" id="5453.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="5453.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="5453.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7432.5"/>
        <Domain name="bounds" id="5453.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="5453.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="5453.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="5453.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="5453.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="5453.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="5453.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="5453.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="5453.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="5453.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="5453.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="5453.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="5453.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="5453.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5453.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="5453.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="5453.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="5453.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="5453.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="5453.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="5453.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="5453.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="5453.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="5453.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="5453.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="5453.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="5453.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5453.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="5453.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="5453.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="5453.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="5453.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="5453.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="5453.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5453.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="5453.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5453.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="5453.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="5453.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="5453.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5453.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="5453.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5453.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="5453.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="5453.Source">
        <Domain name="input_type" id="5453.Source.input_type"/>
      </Property>
      <Property name="Specular" id="5453.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="5453.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="5453.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="5453.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="5453.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="5453.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="5453.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="5453.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="5453.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="5453.Texture">
        <Domain name="groups" id="5453.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="5453.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="5453.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5453.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="5453.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="UnstructuredGridRepresentation" id="7071" servers="21">
      <Property name="CubeAxesVisibility" id="7071.CubeAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.CubeAxesVisibility.bool"/>
      </Property>
      <Property name="Input" id="7071.Input" number_of_elements="1">
        <Proxy value="6973" output_port="0"/>
        <Domain name="input_array_cell" id="7071.Input.input_array_cell"/>
        <Domain name="input_array_point" id="7071.Input.input_array_point"/>
        <Domain name="input_type" id="7071.Input.input_type"/>
      </Property>
      <Property name="Representation" id="7071.Representation" number_of_elements="1">
        <Element index="0" value="Surface"/>
        <Domain name="list" id="7071.Representation.list">
          <String text="3D Glyphs"/>
          <String text="Outline"/>
          <String text="Points"/>
          <String text="Surface"/>
          <String text="Surface With Edges"/>
          <String text="Volume"/>
          <String text="Wireframe"/>
        </Domain>
      </Property>
      <Property name="RepresentationTypesInfo" id="7071.RepresentationTypesInfo" number_of_elements="7">
        <Element index="0" value="3D Glyphs"/>
        <Element index="1" value="Outline"/>
        <Element index="2" value="Points"/>
        <Element index="3" value="Surface"/>
        <Element index="4" value="Surface With Edges"/>
        <Element index="5" value="Volume"/>
        <Element index="6" value="Wireframe"/>
      </Property>
      <Property name="SelectionCellFieldDataArrayName" id="7071.SelectionCellFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="vtkOriginalCellIds"/>
        <Domain name="array_list" id="7071.SelectionCellFieldDataArrayName.array_list"/>
      </Property>
      <Property name="SelectionPointFieldDataArrayName" id="7071.SelectionPointFieldDataArrayName" number_of_elements="1">
        <Element index="0" value="Angular acceleration [rad s^-2]"/>
        <Domain name="array_list" id="7071.SelectionPointFieldDataArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectionVisibility" id="7071.SelectionVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.SelectionVisibility.bool"/>
      </Property>
      <Property name="Visibility" id="7071.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.Visibility.bool"/>
      </Property>
      <Property name="Ambient" id="7071.Ambient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7071.Ambient.range"/>
      </Property>
      <Property name="AmbientColor" id="7071.AmbientColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.AmbientColor.range"/>
      </Property>
      <Property name="AxesOrigin" id="7071.AxesOrigin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.AxesOrigin.range"/>
      </Property>
      <Property name="BackfaceAmbientColor" id="7071.BackfaceAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.BackfaceAmbientColor.range"/>
      </Property>
      <Property name="BackfaceDiffuseColor" id="7071.BackfaceDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.BackfaceDiffuseColor.range"/>
      </Property>
      <Property name="BackfaceOpacity" id="7071.BackfaceOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.BackfaceOpacity.range"/>
      </Property>
      <Property name="BackfaceRepresentation" id="7071.BackfaceRepresentation" number_of_elements="1">
        <Element index="0" value="400"/>
        <Domain name="enum" id="7071.BackfaceRepresentation.enum">
          <Entry value="400" text="Follow Frontface"/>
          <Entry value="401" text="Cull Backface"/>
          <Entry value="402" text="Cull Frontface"/>
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
          <Entry value="3" text="Surface With Edges"/>
        </Domain>
      </Property>
      <Property name="BlockColor" id="7071.BlockColor"/>
      <Property name="BlockColorsDistinctValues" id="7071.BlockColorsDistinctValues" number_of_elements="1">
        <Element index="0" value="12"/>
        <Domain name="range" id="7071.BlockColorsDistinctValues.range"/>
      </Property>
      <Property name="BlockOpacity" id="7071.BlockOpacity"/>
      <Property name="BlockVisibility" id="7071.BlockVisibility"/>
      <Property name="CenterStickyAxes" id="7071.CenterStickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.CenterStickyAxes.bool"/>
      </Property>
      <Property name="ColorArrayName" id="7071.ColorArrayName" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value=""/>
        <Element index="4" value=""/>
        <Domain name="array_list" id="7071.ColorArrayName.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
        <Domain name="field_list" id="7071.ColorArrayName.field_list">
          <Entry value="0" text="Point Data"/>
          <Entry value="1" text="Cell Data"/>
          <Entry value="2" text="Field Data"/>
          <Entry value="4" text="Vertex Data"/>
          <Entry value="5" text="Edge Data"/>
          <Entry value="6" text="Row Data"/>
        </Domain>
      </Property>
      <Property name="CubeAxesColor" id="7071.CubeAxesColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.CubeAxesColor.range"/>
      </Property>
      <Property name="CubeAxesCornerOffset" id="7071.CubeAxesCornerOffset" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7071.CubeAxesCornerOffset.range"/>
      </Property>
      <Property name="CubeAxesFlyMode" id="7071.CubeAxesFlyMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7071.CubeAxesFlyMode.enum">
          <Entry value="0" text="Outer Edges"/>
          <Entry value="1" text="Closest Triad"/>
          <Entry value="2" text="Furthest Triad"/>
          <Entry value="3" text="Static Triad"/>
          <Entry value="4" text="Static Edges"/>
        </Domain>
      </Property>
      <Property name="CubeAxesGridLineLocation" id="7071.CubeAxesGridLineLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.CubeAxesGridLineLocation.enum">
          <Entry value="0" text="All Faces"/>
          <Entry value="1" text="Closest Faces"/>
          <Entry value="2" text="Furthest Faces"/>
        </Domain>
      </Property>
      <Property name="CubeAxesInertia" id="7071.CubeAxesInertia" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.CubeAxesInertia.range"/>
      </Property>
      <Property name="CubeAxesTickLocation" id="7071.CubeAxesTickLocation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.CubeAxesTickLocation.enum">
          <Entry value="0" text="Inside"/>
          <Entry value="1" text="Outside"/>
          <Entry value="2" text="Both"/>
        </Domain>
      </Property>
      <Property name="CubeAxesUseDefaultXTitle" id="7071.CubeAxesUseDefaultXTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesUseDefaultXTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultYTitle" id="7071.CubeAxesUseDefaultYTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesUseDefaultYTitle.bool"/>
      </Property>
      <Property name="CubeAxesUseDefaultZTitle" id="7071.CubeAxesUseDefaultZTitle" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesUseDefaultZTitle.bool"/>
      </Property>
      <Property name="CubeAxesXAxisMinorTickVisibility" id="7071.CubeAxesXAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesXAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisTickVisibility" id="7071.CubeAxesXAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesXAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXAxisVisibility" id="7071.CubeAxesXAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesXAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesXGridLines" id="7071.CubeAxesXGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.CubeAxesXGridLines.bool"/>
      </Property>
      <Property name="CubeAxesXLabelFormat" id="7071.CubeAxesXLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesXTitle" id="7071.CubeAxesXTitle" number_of_elements="1">
        <Element index="0" value="X-Axis"/>
      </Property>
      <Property name="CubeAxesYAxisMinorTickVisibility" id="7071.CubeAxesYAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesYAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisTickVisibility" id="7071.CubeAxesYAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesYAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYAxisVisibility" id="7071.CubeAxesYAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesYAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesYGridLines" id="7071.CubeAxesYGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.CubeAxesYGridLines.bool"/>
      </Property>
      <Property name="CubeAxesYLabelFormat" id="7071.CubeAxesYLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesYTitle" id="7071.CubeAxesYTitle" number_of_elements="1">
        <Element index="0" value="Y-Axis"/>
      </Property>
      <Property name="CubeAxesZAxisMinorTickVisibility" id="7071.CubeAxesZAxisMinorTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesZAxisMinorTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisTickVisibility" id="7071.CubeAxesZAxisTickVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesZAxisTickVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZAxisVisibility" id="7071.CubeAxesZAxisVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.CubeAxesZAxisVisibility.bool"/>
      </Property>
      <Property name="CubeAxesZGridLines" id="7071.CubeAxesZGridLines" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.CubeAxesZGridLines.bool"/>
      </Property>
      <Property name="CubeAxesZLabelFormat" id="7071.CubeAxesZLabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="CubeAxesZTitle" id="7071.CubeAxesZTitle" number_of_elements="1">
        <Element index="0" value="Z-Axis"/>
      </Property>
      <Property name="CustomBounds" id="7071.CustomBounds" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomBoundsActive" id="7071.CustomBoundsActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CustomRange" id="7071.CustomRange" number_of_elements="6">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Element index="3" value="1"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="CustomRangeActive" id="7071.CustomRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="DataBounds" id="7071.DataBounds" number_of_elements="6">
        <Element index="0" value="8.40731984053726e-316"/>
        <Element index="1" value="7.95961415288162e-316"/>
        <Element index="2" value="2.42092166462211e-322"/>
        <Element index="3" value="7.95962521995209e-316"/>
        <Element index="4" value="7.95961336237659e-316"/>
        <Element index="5" value="6.9229924411334e-310"/>
      </Property>
      <Property name="Diffuse" id="7071.Diffuse" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.Diffuse.range"/>
      </Property>
      <Property name="DiffuseColor" id="7071.DiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.DiffuseColor.range"/>
      </Property>
      <Property name="EdgeColor" id="7071.EdgeColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.5"/>
        <Domain name="range" id="7071.EdgeColor.range"/>
      </Property>
      <Property name="ExtractedBlockIndex" id="7071.ExtractedBlockIndex" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="tree" id="7071.ExtractedBlockIndex.tree"/>
      </Property>
      <Property name="InterpolateScalarsBeforeMapping" id="7071.InterpolateScalarsBeforeMapping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.InterpolateScalarsBeforeMapping.bool"/>
      </Property>
      <Property name="Interpolation" id="7071.Interpolation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7071.Interpolation.enum">
          <Entry value="0" text="Flat"/>
          <Entry value="1" text="Gouraud"/>
        </Domain>
      </Property>
      <Property name="LineWidth" id="7071.LineWidth" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.LineWidth.range"/>
      </Property>
      <Property name="LookupTable" id="7071.LookupTable">
        <Domain name="groups" id="7071.LookupTable.groups"/>
      </Property>
      <Property name="MapScalars" id="7071.MapScalars" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.MapScalars.bool"/>
      </Property>
      <Property name="Masking" id="7071.Masking" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.Masking.bool"/>
      </Property>
      <Property name="MeshVisibility" id="7071.MeshVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.MeshVisibility.bool"/>
      </Property>
      <Property name="NonlinearSubdivisionLevel" id="7071.NonlinearSubdivisionLevel" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.NonlinearSubdivisionLevel.range"/>
      </Property>
      <Property name="Opacity" id="7071.Opacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.Opacity.range"/>
      </Property>
      <Property name="Orient" id="7071.Orient" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.Orient.bool"/>
      </Property>
      <Property name="Orientation" id="7071.Orientation" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.Orientation.range"/>
      </Property>
      <Property name="OrientationMode" id="7071.OrientationMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.OrientationMode.enum">
          <Entry value="0" text="Direction"/>
          <Entry value="1" text="Rotation"/>
        </Domain>
      </Property>
      <Property name="Origin" id="7071.Origin" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.Origin.range"/>
      </Property>
      <Property name="OriginalBoundsRangeActive" id="7071.OriginalBoundsRangeActive" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="Pickable" id="7071.Pickable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7071.Pickable.bool"/>
      </Property>
      <Property name="PointSize" id="7071.PointSize" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7071.PointSize.range"/>
      </Property>
      <Property name="Position" id="7071.Position" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.Position.range"/>
      </Property>
      <Property name="ScalarOpacityFunction" id="7071.ScalarOpacityFunction">
        <Domain name="groups" id="7071.ScalarOpacityFunction.groups"/>
      </Property>
      <Property name="ScalarOpacityUnitDistance" id="7071.ScalarOpacityUnitDistance" number_of_elements="1">
        <Element index="0" value="89202.9217570815"/>
        <Domain name="bounds" id="7071.ScalarOpacityUnitDistance.bounds"/>
      </Property>
      <Property name="Scale" id="7071.Scale" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.Scale.range"/>
      </Property>
      <Property name="ScaleFactor" id="7071.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7432.5"/>
        <Domain name="bounds" id="7071.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7071.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7071.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7071.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.ScaleMode.enum">
          <Entry value="0" text="No Data Scaling Off"/>
          <Entry value="1" text="Magnitude"/>
          <Entry value="2" text="Vector Components"/>
        </Domain>
      </Property>
      <Property name="Scaling" id="7071.Scaling" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.Scaling.bool"/>
      </Property>
      <Property name="SelectMapper" id="7071.SelectMapper" number_of_elements="1">
        <Element index="0" value="Projected tetra"/>
        <Domain name="list" id="7071.SelectMapper.list">
          <String text="Projected tetra"/>
          <String text="HAVS"/>
          <String text="Z sweep"/>
          <String text="Bunyk ray cast"/>
        </Domain>
      </Property>
      <Property name="SelectMaskArray" id="7071.SelectMaskArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectOrientationVectors" id="7071.SelectOrientationVectors" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectScaleArray" id="7071.SelectScaleArray" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelBold" id="7071.SelectionCellLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionCellLabelBold.bool"/>
      </Property>
      <Property name="SelectionCellLabelColor" id="7071.SelectionCellLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.SelectionCellLabelColor.range"/>
      </Property>
      <Property name="SelectionCellLabelFontFamily" id="7071.SelectionCellLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.SelectionCellLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelFontSize" id="7071.SelectionCellLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7071.SelectionCellLabelFontSize.range"/>
      </Property>
      <Property name="SelectionCellLabelFormat" id="7071.SelectionCellLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionCellLabelItalic" id="7071.SelectionCellLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionCellLabelItalic.bool"/>
      </Property>
      <Property name="SelectionCellLabelJustification" id="7071.SelectionCellLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.SelectionCellLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionCellLabelOpacity" id="7071.SelectionCellLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.SelectionCellLabelOpacity.range"/>
      </Property>
      <Property name="SelectionCellLabelShadow" id="7071.SelectionCellLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionCellLabelShadow.bool"/>
      </Property>
      <Property name="SelectionCellLabelVisibility" id="7071.SelectionCellLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionCellLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionColor" id="7071.SelectionColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.SelectionColor.range"/>
      </Property>
      <Property name="SelectionLineWidth" id="7071.SelectionLineWidth" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="7071.SelectionLineWidth.range"/>
      </Property>
      <Property name="SelectionOpacity" id="7071.SelectionOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.SelectionOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelBold" id="7071.SelectionPointLabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionPointLabelBold.bool"/>
      </Property>
      <Property name="SelectionPointLabelColor" id="7071.SelectionPointLabelColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="7071.SelectionPointLabelColor.range"/>
      </Property>
      <Property name="SelectionPointLabelFontFamily" id="7071.SelectionPointLabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.SelectionPointLabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelFontSize" id="7071.SelectionPointLabelFontSize" number_of_elements="1">
        <Element index="0" value="18"/>
        <Domain name="range" id="7071.SelectionPointLabelFontSize.range"/>
      </Property>
      <Property name="SelectionPointLabelFormat" id="7071.SelectionPointLabelFormat" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="SelectionPointLabelItalic" id="7071.SelectionPointLabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionPointLabelItalic.bool"/>
      </Property>
      <Property name="SelectionPointLabelJustification" id="7071.SelectionPointLabelJustification" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7071.SelectionPointLabelJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Center"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="SelectionPointLabelOpacity" id="7071.SelectionPointLabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7071.SelectionPointLabelOpacity.range"/>
      </Property>
      <Property name="SelectionPointLabelShadow" id="7071.SelectionPointLabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionPointLabelShadow.bool"/>
      </Property>
      <Property name="SelectionPointLabelVisibility" id="7071.SelectionPointLabelVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionPointLabelVisibility.bool"/>
      </Property>
      <Property name="SelectionPointSize" id="7071.SelectionPointSize" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="7071.SelectionPointSize.range"/>
      </Property>
      <Property name="SelectionRepresentation" id="7071.SelectionRepresentation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7071.SelectionRepresentation.enum">
          <Entry value="0" text="Points"/>
          <Entry value="1" text="Wireframe"/>
          <Entry value="2" text="Surface"/>
        </Domain>
      </Property>
      <Property name="SelectionUseOutline" id="7071.SelectionUseOutline" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SelectionUseOutline.bool"/>
      </Property>
      <Property name="Source" id="7071.Source">
        <Domain name="input_type" id="7071.Source.input_type"/>
      </Property>
      <Property name="Specular" id="7071.Specular" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="7071.Specular.range"/>
      </Property>
      <Property name="SpecularColor" id="7071.SpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="7071.SpecularColor.range"/>
      </Property>
      <Property name="SpecularPower" id="7071.SpecularPower" number_of_elements="1">
        <Element index="0" value="100"/>
        <Domain name="range" id="7071.SpecularPower.range"/>
      </Property>
      <Property name="StaticMode" id="7071.StaticMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.StaticMode.bool"/>
      </Property>
      <Property name="StickyAxes" id="7071.StickyAxes" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.StickyAxes.bool"/>
      </Property>
      <Property name="SuppressLOD" id="7071.SuppressLOD" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.SuppressLOD.bool"/>
      </Property>
      <Property name="Texture" id="7071.Texture">
        <Domain name="groups" id="7071.Texture.groups"/>
      </Property>
      <Property name="Triangulate" id="7071.Triangulate" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.Triangulate.bool"/>
      </Property>
      <Property name="UseAxesOrigin" id="7071.UseAxesOrigin" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="7071.UseAxesOrigin.bool"/>
      </Property>
      <Property name="UserTransform" id="7071.UserTransform" number_of_elements="16">
        <Element index="0" value="1"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Element index="3" value="0"/>
        <Element index="4" value="0"/>
        <Element index="5" value="1"/>
        <Element index="6" value="0"/>
        <Element index="7" value="0"/>
        <Element index="8" value="0"/>
        <Element index="9" value="0"/>
        <Element index="10" value="1"/>
        <Element index="11" value="0"/>
        <Element index="12" value="0"/>
        <Element index="13" value="0"/>
        <Element index="14" value="0"/>
        <Element index="15" value="1"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5235" servers="21">
      <Property name="Enabled" id="5235.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5235.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5235.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5235.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5235.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5235.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5235.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5235.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5235.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5235.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5235.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5235.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5235.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5235.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5235.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5235.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5235.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5235.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5235.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5235.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5235.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5235.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5235.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5235.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5235.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5235.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5235.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5235.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5235.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5235.LookupTable" number_of_elements="1">
        <Proxy value="5229"/>
        <Domain name="groups" id="5235.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5235.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5235.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5235.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5235.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5235.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5235.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5235.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5235.Position.range"/>
      </Property>
      <Property name="Position2" id="5235.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5235.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5235.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5235.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5235.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5235.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5235.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5235.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5235.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5235.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5235.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5235.Title" number_of_elements="1">
        <Element index="0" value="Diameter (areal) [m]"/>
      </Property>
      <Property name="TitleBold" id="5235.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5235.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5235.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5235.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5235.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5235.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5235.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5235.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5235.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5235.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5235.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5235.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5235.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5235.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5235.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5894" servers="21">
      <Property name="Enabled" id="5894.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5894.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5894.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5894.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5894.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5894.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5894.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5894.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5894.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5894.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5894.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5894.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5894.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5894.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5894.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5894.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5894.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5894.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5894.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5894.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5894.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5894.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5894.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5894.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5894.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5894.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5894.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5894.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5894.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5894.LookupTable" number_of_elements="1">
        <Proxy value="5888"/>
        <Domain name="groups" id="5894.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5894.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5894.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5894.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5894.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5894.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5894.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5894.Position" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.52"/>
        <Domain name="range" id="5894.Position.range"/>
      </Property>
      <Property name="Position2" id="5894.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5894.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5894.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5894.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5894.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5894.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5894.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5894.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5894.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5894.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5894.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5894.Title" number_of_elements="1">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
      </Property>
      <Property name="TitleBold" id="5894.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5894.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5894.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5894.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5894.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5894.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5894.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5894.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5894.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5894.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5894.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5894.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5894.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5894.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5894.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="6133" servers="21">
      <Property name="Enabled" id="6133.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="6133.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="6133.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="6133.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="6133.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="6133.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="6133.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="6133.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="6133.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="6133.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="6133.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="6133.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="6133.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="6133.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="6133.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="6133.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="6133.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="6133.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="6133.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6133.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="6133.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6133.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="6133.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="6133.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="6133.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="6133.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="6133.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6133.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="6133.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="6133.LookupTable" number_of_elements="1">
        <Proxy value="6072"/>
        <Domain name="groups" id="6133.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="6133.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="6133.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6133.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="6133.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6133.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="6133.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="6133.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.52"/>
        <Domain name="range" id="6133.Position.range"/>
      </Property>
      <Property name="Position2" id="6133.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="6133.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="6133.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="6133.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.52"/>
      </Property>
      <Property name="RangeLabelFormat" id="6133.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="6133.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="6133.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="6133.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="6133.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6133.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="6133.Title" number_of_elements="1">
        <Element index="0" value="Tensile stress [Pa]"/>
      </Property>
      <Property name="TitleBold" id="6133.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="6133.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6133.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="6133.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6133.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="6133.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="6133.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="6133.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="6133.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6133.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="6133.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6133.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="6133.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6133.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="6133.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6133.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="6699" servers="21">
      <Property name="Enabled" id="6699.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="6699.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="6699.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="6699.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="6699.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="6699.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="6699.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="6699.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="6699.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="6699.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="6699.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="6699.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="6699.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="6699.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="6699.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="6699.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="6699.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="6699.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="6699.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6699.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="6699.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6699.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="6699.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="6699.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="6699.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="6699.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="6699.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6699.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="6699.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="6699.LookupTable" number_of_elements="1">
        <Proxy value="6693"/>
        <Domain name="groups" id="6699.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="6699.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="6699.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="6699.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="6699.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6699.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="6699.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="6699.Position" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.52"/>
        <Domain name="range" id="6699.Position.range"/>
      </Property>
      <Property name="Position2" id="6699.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="6699.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="6699.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="6699.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="6699.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="6699.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="6699.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6699.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="6699.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="6699.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6699.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="6699.Title" number_of_elements="1">
        <Element index="0" value="v: Meridional velocity [m/s]"/>
      </Property>
      <Property name="TitleBold" id="6699.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="6699.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="6699.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="6699.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6699.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="6699.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="6699.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="6699.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="6699.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6699.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="6699.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6699.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="6699.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="6699.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6699.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5243" servers="21">
      <Property name="Enabled" id="5243.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5243.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5243.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5243.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5243.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5243.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5243.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5243.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5243.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5243.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5243.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5243.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5243.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5243.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5243.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5243.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5243.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5243.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5243.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5243.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5243.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5243.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5243.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5243.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5243.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5243.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5243.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5243.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5243.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5243.LookupTable" number_of_elements="1">
        <Proxy value="5237"/>
        <Domain name="groups" id="5243.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5243.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5243.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5243.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5243.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5243.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5243.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5243.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5243.Position.range"/>
      </Property>
      <Property name="Position2" id="5243.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5243.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5243.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5243.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5243.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5243.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5243.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5243.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5243.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5243.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5243.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5243.Title" number_of_elements="1">
        <Element index="0" value="Number of contacts [-]"/>
      </Property>
      <Property name="TitleBold" id="5243.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5243.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5243.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5243.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5243.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5243.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5243.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5243.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5243.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5243.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5243.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5243.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5243.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5243.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5243.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5251" servers="21">
      <Property name="Enabled" id="5251.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5251.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5251.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5251.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5251.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5251.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5251.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5251.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5251.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5251.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5251.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5251.ComponentTitle" number_of_elements="1">
        <Element index="0" value="Magnitude"/>
      </Property>
      <Property name="DrawAnnotations" id="5251.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5251.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5251.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5251.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5251.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5251.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5251.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5251.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5251.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5251.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5251.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5251.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5251.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5251.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5251.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5251.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5251.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5251.LookupTable" number_of_elements="1">
        <Proxy value="5245"/>
        <Domain name="groups" id="5251.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5251.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5251.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5251.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5251.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5251.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5251.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5251.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5251.Position.range"/>
      </Property>
      <Property name="Position2" id="5251.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5251.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5251.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5251.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5251.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5251.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5251.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5251.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5251.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5251.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5251.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5251.Title" number_of_elements="1">
        <Element index="0" value="Linear acceleration [m s^-2]"/>
      </Property>
      <Property name="TitleBold" id="5251.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5251.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5251.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5251.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5251.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5251.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5251.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5251.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5251.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5251.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5251.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5251.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5251.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5251.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5251.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5259" servers="21">
      <Property name="Enabled" id="5259.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5259.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5259.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5259.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5259.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5259.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5259.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5259.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5259.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5259.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5259.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5259.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5259.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5259.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5259.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5259.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5259.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5259.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5259.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5259.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5259.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5259.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5259.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5259.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5259.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5259.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5259.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5259.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5259.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5259.LookupTable" number_of_elements="1">
        <Proxy value="5253"/>
        <Domain name="groups" id="5259.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5259.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5259.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5259.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5259.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5259.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5259.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5259.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5259.Position.range"/>
      </Property>
      <Property name="Position2" id="5259.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5259.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5259.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5259.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5259.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5259.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5259.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5259.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5259.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5259.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5259.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5259.Title" number_of_elements="1">
        <Element index="0" value="Fixed in space [-]"/>
      </Property>
      <Property name="TitleBold" id="5259.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5259.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5259.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5259.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5259.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5259.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5259.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5259.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5259.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5259.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5259.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5259.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5259.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5259.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5259.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5267" servers="21">
      <Property name="Enabled" id="5267.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5267.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5267.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5267.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5267.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5267.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5267.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5267.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5267.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5267.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5267.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5267.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5267.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5267.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5267.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5267.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5267.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5267.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5267.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5267.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5267.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5267.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5267.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5267.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5267.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5267.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5267.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5267.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5267.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5267.LookupTable" number_of_elements="1">
        <Proxy value="5261"/>
        <Domain name="groups" id="5267.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5267.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5267.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5267.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5267.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5267.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5267.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5267.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5267.Position.range"/>
      </Property>
      <Property name="Position2" id="5267.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5267.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5267.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5267.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5267.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5267.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5267.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5267.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5267.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5267.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5267.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5267.Title" number_of_elements="1">
        <Element index="0" value="Enabled [-]"/>
      </Property>
      <Property name="TitleBold" id="5267.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5267.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5267.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5267.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5267.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5267.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5267.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5267.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5267.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5267.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5267.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5267.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5267.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5267.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5267.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5275" servers="21">
      <Property name="Enabled" id="5275.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5275.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5275.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5275.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5275.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5275.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5275.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5275.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5275.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5275.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5275.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5275.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5275.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5275.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5275.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5275.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5275.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5275.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5275.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5275.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5275.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5275.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5275.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5275.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5275.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5275.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5275.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5275.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5275.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5275.LookupTable" number_of_elements="1">
        <Proxy value="5269"/>
        <Domain name="groups" id="5275.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5275.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5275.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5275.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5275.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5275.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5275.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5275.Position" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5275.Position.range"/>
      </Property>
      <Property name="Position2" id="5275.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5275.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5275.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5275.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.52"/>
      </Property>
      <Property name="RangeLabelFormat" id="5275.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5275.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5275.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5275.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5275.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5275.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5275.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5275.Title" number_of_elements="1">
        <Element index="0" value="e: Relative interface height [m]"/>
      </Property>
      <Property name="TitleBold" id="5275.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5275.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5275.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5275.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5275.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5275.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5275.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5275.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5275.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5275.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5275.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5275.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5275.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5275.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5275.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5639" servers="21">
      <Property name="Enabled" id="5639.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5639.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5639.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5639.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5639.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5639.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5639.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5639.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5639.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5639.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5639.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5639.ComponentTitle" number_of_elements="1">
        <Element index="0" value="Magnitude"/>
      </Property>
      <Property name="DrawAnnotations" id="5639.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5639.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5639.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5639.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5639.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5639.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5639.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5639.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5639.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5639.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5639.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5639.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5639.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5639.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5639.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5639.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5639.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5639.LookupTable" number_of_elements="1">
        <Proxy value="5299"/>
        <Domain name="groups" id="5639.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5639.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5639.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5639.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5639.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5639.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5639.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5639.Position" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5639.Position.range"/>
      </Property>
      <Property name="Position2" id="5639.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5639.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5639.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5639.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5639.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5639.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5639.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5639.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5639.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5639.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5639.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5639.Title" number_of_elements="1">
        <Element index="0" value="Velocity vector [m/s]"/>
      </Property>
      <Property name="TitleBold" id="5639.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5639.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5639.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5639.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5639.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5639.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5639.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5639.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5639.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5639.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5639.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5639.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5639.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5639.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5639.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5803" servers="21">
      <Property name="Enabled" id="5803.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5803.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5803.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5803.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5803.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5803.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5803.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5803.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5803.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5803.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5803.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5803.ComponentTitle" number_of_elements="1">
        <Element index="0" value=""/>
      </Property>
      <Property name="DrawAnnotations" id="5803.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5803.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5803.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5803.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5803.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5803.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5803.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5803.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5803.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5803.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5803.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5803.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5803.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5803.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5803.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5803.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5803.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5803.LookupTable" number_of_elements="1">
        <Proxy value="5797"/>
        <Domain name="groups" id="5803.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5803.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5803.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5803.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5803.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5803.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5803.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5803.Position" number_of_elements="2">
        <Element index="0" value="0.1285"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5803.Position.range"/>
      </Property>
      <Property name="Position2" id="5803.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5803.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5803.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5803.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.52"/>
      </Property>
      <Property name="RangeLabelFormat" id="5803.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5803.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5803.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5803.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5803.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5803.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5803.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5803.Title" number_of_elements="1">
        <Element index="0" value="Atmosphere drag coefficient (horizontal) [-]"/>
      </Property>
      <Property name="TitleBold" id="5803.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5803.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5803.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5803.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5803.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5803.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5803.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5803.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5803.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5803.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5803.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5803.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5803.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5803.Visibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5803.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="representations" type="ScalarBarWidgetRepresentation" id="5809" servers="21">
      <Property name="Enabled" id="5809.Enabled" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.Enabled.bool"/>
      </Property>
      <Property name="LockPosition" id="5809.LockPosition" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.LockPosition.bool"/>
      </Property>
      <Property name="UseNonCompositedRenderer" id="5809.UseNonCompositedRenderer" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="AddRangeAnnotations" id="5809.AddRangeAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.AddRangeAnnotations.bool"/>
      </Property>
      <Property name="AddRangeLabels" id="5809.AddRangeLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.AddRangeLabels.bool"/>
      </Property>
      <Property name="AspectRatio" id="5809.AspectRatio" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="5809.AspectRatio.range"/>
      </Property>
      <Property name="AutoOrient" id="5809.AutoOrient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.AutoOrient.bool"/>
      </Property>
      <Property name="AutoOrientInfo" id="5809.AutoOrientInfo" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.AutoOrientInfo.bool"/>
      </Property>
      <Property name="AutomaticAnnotations" id="5809.AutomaticAnnotations" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.AutomaticAnnotations.bool"/>
      </Property>
      <Property name="AutomaticLabelFormat" id="5809.AutomaticLabelFormat" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.AutomaticLabelFormat.bool"/>
      </Property>
      <Property name="ComponentTitle" id="5809.ComponentTitle" number_of_elements="1">
        <Element index="0" value="Magnitude"/>
      </Property>
      <Property name="DrawAnnotations" id="5809.DrawAnnotations" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.DrawAnnotations.bool"/>
      </Property>
      <Property name="DrawNanAnnotation" id="5809.DrawNanAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.DrawNanAnnotation.bool"/>
      </Property>
      <Property name="DrawSubTickMarks" id="5809.DrawSubTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.DrawSubTickMarks.bool"/>
      </Property>
      <Property name="DrawTickLabels" id="5809.DrawTickLabels" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.DrawTickLabels.bool"/>
      </Property>
      <Property name="DrawTickMarks" id="5809.DrawTickMarks" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.DrawTickMarks.bool"/>
      </Property>
      <Property name="LabelBold" id="5809.LabelBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.LabelBold.bool"/>
      </Property>
      <Property name="LabelColor" id="5809.LabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5809.LabelColor.range"/>
      </Property>
      <Property name="LabelFontFamily" id="5809.LabelFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5809.LabelFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="LabelFontSize" id="5809.LabelFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5809.LabelFontSize.range"/>
      </Property>
      <Property name="LabelFormat" id="5809.LabelFormat" number_of_elements="1">
        <Element index="0" value="%-#6.3g"/>
      </Property>
      <Property name="LabelItalic" id="5809.LabelItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.LabelItalic.bool"/>
      </Property>
      <Property name="LabelOpacity" id="5809.LabelOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5809.LabelOpacity.range"/>
      </Property>
      <Property name="LabelShadow" id="5809.LabelShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.LabelShadow.bool"/>
      </Property>
      <Property name="LookupTable" id="5809.LookupTable" number_of_elements="1">
        <Proxy value="5740"/>
        <Domain name="groups" id="5809.LookupTable.groups"/>
      </Property>
      <Property name="NanAnnotation" id="5809.NanAnnotation" number_of_elements="1">
        <Element index="0" value="NaN"/>
      </Property>
      <Property name="NumberOfLabels" id="5809.NumberOfLabels" number_of_elements="1">
        <Element index="0" value="5"/>
        <Domain name="range" id="5809.NumberOfLabels.range"/>
      </Property>
      <Property name="Orientation" id="5809.Orientation" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5809.Orientation.enum">
          <Entry value="0" text="Horizontal"/>
          <Entry value="1" text="Vertical"/>
        </Domain>
      </Property>
      <Property name="OrientationInfo" id="5809.OrientationInfo" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="Position" id="5809.Position" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
        <Domain name="range" id="5809.Position.range"/>
      </Property>
      <Property name="Position2" id="5809.Position2" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
        <Domain name="range" id="5809.Position2.range"/>
      </Property>
      <Property name="Position2Info" id="5809.Position2Info" number_of_elements="2">
        <Element index="0" value="0.12"/>
        <Element index="1" value="0.43"/>
      </Property>
      <Property name="PositionInfo" id="5809.PositionInfo" number_of_elements="2">
        <Element index="0" value="0.85"/>
        <Element index="1" value="0.05"/>
      </Property>
      <Property name="RangeLabelFormat" id="5809.RangeLabelFormat" number_of_elements="1">
        <Element index="0" value="%4.3e"/>
      </Property>
      <Property name="Repositionable" id="5809.Repositionable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.Repositionable.bool"/>
      </Property>
      <Property name="Resizable" id="5809.Resizable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.Resizable.bool"/>
      </Property>
      <Property name="Selectable" id="5809.Selectable" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.Selectable.bool"/>
      </Property>
      <Property name="TextPosition" id="5809.TextPosition" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5809.TextPosition.enum">
          <Entry value="1" text="Ticks right/top, annotations left/bottom"/>
          <Entry value="0" text="Ticks left/bottom, annotations right/top"/>
        </Domain>
      </Property>
      <Property name="Title" id="5809.Title" number_of_elements="1">
        <Element index="0" value="Linear velocity [m s^-1]"/>
      </Property>
      <Property name="TitleBold" id="5809.TitleBold" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.TitleBold.bool"/>
      </Property>
      <Property name="TitleColor" id="5809.TitleColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
        <Domain name="range" id="5809.TitleColor.range"/>
      </Property>
      <Property name="TitleFontFamily" id="5809.TitleFontFamily" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="5809.TitleFontFamily.enum">
          <Entry value="0" text="Arial"/>
          <Entry value="1" text="Courier"/>
          <Entry value="2" text="Times"/>
        </Domain>
      </Property>
      <Property name="TitleFontSize" id="5809.TitleFontSize" number_of_elements="1">
        <Element index="0" value="7"/>
        <Domain name="range" id="5809.TitleFontSize.range"/>
      </Property>
      <Property name="TitleItalic" id="5809.TitleItalic" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.TitleItalic.bool"/>
      </Property>
      <Property name="TitleJustification" id="5809.TitleJustification" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="5809.TitleJustification.enum">
          <Entry value="0" text="Left"/>
          <Entry value="1" text="Centered"/>
          <Entry value="2" text="Right"/>
        </Domain>
      </Property>
      <Property name="TitleOpacity" id="5809.TitleOpacity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="5809.TitleOpacity.range"/>
      </Property>
      <Property name="TitleShadow" id="5809.TitleShadow" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="5809.TitleShadow.bool"/>
      </Property>
      <Property name="Visibility" id="5809.Visibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="5809.Visibility.bool"/>
      </Property>
    </Proxy>
    <Proxy group="filters" type="TransformFilter" id="6326" servers="1">
      <Property name="Input" id="6326.Input" number_of_elements="1">
        <Proxy value="4765" output_port="0"/>
        <Domain name="groups" id="6326.Input.groups"/>
        <Domain name="input_type" id="6326.Input.input_type"/>
      </Property>
      <Property name="Transform" id="6326.Transform" number_of_elements="1">
        <Proxy value="6325"/>
        <Domain name="groups" id="6326.Transform.groups"/>
        <Domain name="proxy_list" id="6326.Transform.proxy_list">
          <Proxy value="6325"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="Glyph" id="6603" servers="1">
      <Property name="GlyphMode" id="6603.GlyphMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6603.GlyphMode.enum">
          <Entry value="0" text="All Points"/>
          <Entry value="1" text="Every Nth Point"/>
          <Entry value="2" text="Uniform Spatial Distribution"/>
        </Domain>
      </Property>
      <Property name="GlyphTransform" id="6603.GlyphTransform" number_of_elements="1">
        <Proxy value="6525"/>
        <Domain name="proxy_list" id="6603.GlyphTransform.proxy_list">
          <Proxy value="6525"/>
        </Domain>
      </Property>
      <Property name="Input" id="6603.Input" number_of_elements="1">
        <Proxy value="6326" output_port="0"/>
        <Domain name="groups" id="6603.Input.groups"/>
        <Domain name="input_array1" id="6603.Input.input_array1"/>
        <Domain name="input_array2" id="6603.Input.input_array2"/>
        <Domain name="input_type" id="6603.Input.input_type"/>
      </Property>
      <Property name="MaximumNumberOfSamplePoints" id="6603.MaximumNumberOfSamplePoints" number_of_elements="1">
        <Element index="0" value="5000"/>
        <Domain name="range" id="6603.MaximumNumberOfSamplePoints.range"/>
      </Property>
      <Property name="Orient" id="6603.Orient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6603.Orient.bool"/>
      </Property>
      <Property name="Scalars" id="6603.Scalars" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="u: Zonal velocity [m/s]"/>
        <Domain name="array_list" id="6603.Scalars.array_list">
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="ScaleFactor" id="6603.ScaleFactor" number_of_elements="1">
        <Element index="0" value="500"/>
        <Domain name="bounds" id="6603.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="6603.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="6603.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="6603.ScaleMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6603.ScaleMode.enum">
          <Entry value="0" text="scalar"/>
          <Entry value="1" text="vector"/>
          <Entry value="2" text="vector_components"/>
          <Entry value="3" text="off"/>
        </Domain>
      </Property>
      <Property name="Seed" id="6603.Seed" number_of_elements="1">
        <Element index="0" value="10339"/>
        <Domain name="range" id="6603.Seed.range"/>
      </Property>
      <Property name="Source" id="6603.Source" number_of_elements="1">
        <Proxy value="6526" output_port="0"/>
        <Domain name="groups" id="6603.Source.groups"/>
        <Domain name="input_type" id="6603.Source.input_type"/>
        <Domain name="proxy_list" id="6603.Source.proxy_list">
          <Proxy value="6526"/>
          <Proxy value="6537"/>
          <Proxy value="6548"/>
          <Proxy value="6559"/>
          <Proxy value="6570"/>
          <Proxy value="6581"/>
          <Proxy value="6592"/>
        </Domain>
      </Property>
      <Property name="Stride" id="6603.Stride" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6603.Stride.range"/>
      </Property>
      <Property name="Vectors" id="6603.Vectors" number_of_elements="5">
        <Element index="0" value="1"/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="6603.Vectors.array_list">
          <String text="Velocity vector [m/s]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="CellDataToPointData" id="6794" servers="1">
      <Property name="Input" id="6794.Input" number_of_elements="1">
        <Proxy value="6701" output_port="0"/>
        <Domain name="groups" id="6794.Input.groups"/>
        <Domain name="input_array" id="6794.Input.input_array"/>
        <Domain name="input_type" id="6794.Input.input_type"/>
      </Property>
      <Property name="PassCellData" id="6794.PassCellData" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6794.PassCellData.bool"/>
      </Property>
      <Property name="PieceInvariant" id="6794.PieceInvariant" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6794.PieceInvariant.bool"/>
      </Property>
    </Proxy>
    <Proxy group="filters" type="TransformFilter" id="6701" servers="1">
      <Property name="Input" id="6701.Input" number_of_elements="1">
        <Proxy value="4787" output_port="0"/>
        <Domain name="groups" id="6701.Input.groups"/>
        <Domain name="input_type" id="6701.Input.input_type"/>
      </Property>
      <Property name="Transform" id="6701.Transform" number_of_elements="1">
        <Proxy value="6700"/>
        <Domain name="groups" id="6701.Transform.groups"/>
        <Domain name="proxy_list" id="6701.Transform.proxy_list">
          <Proxy value="6700"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="TubeFilter" id="6883" servers="1">
      <Property name="Capping" id="6883.Capping" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6883.Capping.bool"/>
      </Property>
      <Property name="DefaultNormal" id="6883.DefaultNormal" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6883.DefaultNormal.range"/>
      </Property>
      <Property name="Input" id="6883.Input" number_of_elements="1">
        <Proxy value="6794" output_port="0"/>
        <Domain name="groups" id="6883.Input.groups"/>
        <Domain name="input_array1" id="6883.Input.input_array1"/>
        <Domain name="input_array2" id="6883.Input.input_array2"/>
        <Domain name="input_type" id="6883.Input.input_type"/>
      </Property>
      <Property name="NumberOfSides" id="6883.NumberOfSides" number_of_elements="1">
        <Element index="0" value="6"/>
        <Domain name="range" id="6883.NumberOfSides.range"/>
      </Property>
      <Property name="Radius" id="6883.Radius" number_of_elements="1">
        <Element index="0" value="250"/>
        <Domain name="bounds" id="6883.Radius.bounds"/>
      </Property>
      <Property name="RadiusFactor" id="6883.RadiusFactor" number_of_elements="1">
        <Element index="0" value="250"/>
        <Domain name="range" id="6883.RadiusFactor.range"/>
      </Property>
      <Property name="SelectInputScalars" id="6883.SelectInputScalars" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Tensile stress [Pa]"/>
        <Domain name="array_list" id="6883.SelectInputScalars.array_list">
          <String text="Contact age [s]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Effective radius [m]"/>
          <String text="Force [N]"/>
          <String text="Tensile stress [Pa]"/>
        </Domain>
      </Property>
      <Property name="SelectInputVectors" id="6883.SelectInputVectors" number_of_elements="5">
        <Element index="0" value="1"/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Inter-particle vector [m]"/>
        <Domain name="array_list" id="6883.SelectInputVectors.array_list">
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
        </Domain>
      </Property>
      <Property name="UseDefaultNormal" id="6883.UseDefaultNormal" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6883.UseDefaultNormal.bool"/>
      </Property>
      <Property name="VaryRadius" id="6883.VaryRadius" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6883.VaryRadius.enum">
          <Entry value="0" text="Off"/>
          <Entry value="1" text="By Scalar"/>
          <Entry value="2" text="By Vector"/>
          <Entry value="3" text="By Absolute Scalar"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="Glyph" id="7327" servers="1">
      <Property name="GlyphMode" id="7327.GlyphMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7327.GlyphMode.enum">
          <Entry value="0" text="All Points"/>
          <Entry value="1" text="Every Nth Point"/>
          <Entry value="2" text="Uniform Spatial Distribution"/>
        </Domain>
      </Property>
      <Property name="GlyphTransform" id="7327.GlyphTransform" number_of_elements="1">
        <Proxy value="7249"/>
        <Domain name="proxy_list" id="7327.GlyphTransform.proxy_list">
          <Proxy value="7249"/>
        </Domain>
      </Property>
      <Property name="Input" id="7327.Input" number_of_elements="1">
        <Proxy value="6973" output_port="0"/>
        <Domain name="groups" id="7327.Input.groups"/>
        <Domain name="input_array1" id="7327.Input.input_array1"/>
        <Domain name="input_array2" id="7327.Input.input_array2"/>
        <Domain name="input_type" id="7327.Input.input_type"/>
      </Property>
      <Property name="MaximumNumberOfSamplePoints" id="7327.MaximumNumberOfSamplePoints" number_of_elements="1">
        <Element index="0" value="5000"/>
        <Domain name="range" id="7327.MaximumNumberOfSamplePoints.range"/>
      </Property>
      <Property name="Orient" id="7327.Orient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7327.Orient.bool"/>
      </Property>
      <Property name="Scalars" id="7327.Scalars" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Diameter (areal) [m]"/>
        <Domain name="array_list" id="7327.Scalars.array_list">
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="ScaleFactor" id="7327.ScaleFactor" number_of_elements="1">
        <Element index="0" value="3000"/>
        <Domain name="bounds" id="7327.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7327.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7327.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7327.ScaleMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7327.ScaleMode.enum">
          <Entry value="0" text="scalar"/>
          <Entry value="1" text="vector"/>
          <Entry value="2" text="vector_components"/>
          <Entry value="3" text="off"/>
        </Domain>
      </Property>
      <Property name="Seed" id="7327.Seed" number_of_elements="1">
        <Element index="0" value="10339"/>
        <Domain name="range" id="7327.Seed.range"/>
      </Property>
      <Property name="Source" id="7327.Source" number_of_elements="1">
        <Proxy value="7250" output_port="0"/>
        <Domain name="groups" id="7327.Source.groups"/>
        <Domain name="input_type" id="7327.Source.input_type"/>
        <Domain name="proxy_list" id="7327.Source.proxy_list">
          <Proxy value="7250"/>
          <Proxy value="7261"/>
          <Proxy value="7272"/>
          <Proxy value="7283"/>
          <Proxy value="7294"/>
          <Proxy value="7305"/>
          <Proxy value="7316"/>
        </Domain>
      </Property>
      <Property name="Stride" id="7327.Stride" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7327.Stride.range"/>
      </Property>
      <Property name="Vectors" id="7327.Vectors" number_of_elements="5">
        <Element index="0" value="1"/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Linear velocity [m s^-1]"/>
        <Domain name="array_list" id="7327.Vectors.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="Glyph" id="7160" servers="1">
      <Property name="GlyphMode" id="7160.GlyphMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7160.GlyphMode.enum">
          <Entry value="0" text="All Points"/>
          <Entry value="1" text="Every Nth Point"/>
          <Entry value="2" text="Uniform Spatial Distribution"/>
        </Domain>
      </Property>
      <Property name="GlyphTransform" id="7160.GlyphTransform" number_of_elements="1">
        <Proxy value="7082"/>
        <Domain name="proxy_list" id="7160.GlyphTransform.proxy_list">
          <Proxy value="7082"/>
        </Domain>
      </Property>
      <Property name="Input" id="7160.Input" number_of_elements="1">
        <Proxy value="6973" output_port="0"/>
        <Domain name="groups" id="7160.Input.groups"/>
        <Domain name="input_array1" id="7160.Input.input_array1"/>
        <Domain name="input_array2" id="7160.Input.input_array2"/>
        <Domain name="input_type" id="7160.Input.input_type"/>
      </Property>
      <Property name="MaximumNumberOfSamplePoints" id="7160.MaximumNumberOfSamplePoints" number_of_elements="1">
        <Element index="0" value="5000"/>
        <Domain name="range" id="7160.MaximumNumberOfSamplePoints.range"/>
      </Property>
      <Property name="Orient" id="7160.Orient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7160.Orient.bool"/>
      </Property>
      <Property name="Scalars" id="7160.Scalars" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Diameter (areal) [m]"/>
        <Domain name="array_list" id="7160.Scalars.array_list">
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Circumreference  [m]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Density [kg m^-3]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Enabled [-]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Number of contacts [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Thickness [m]"/>
          <String text="Volume [m^3]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
        </Domain>
      </Property>
      <Property name="ScaleFactor" id="7160.ScaleFactor" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bounds" id="7160.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7160.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7160.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7160.ScaleMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7160.ScaleMode.enum">
          <Entry value="0" text="scalar"/>
          <Entry value="1" text="vector"/>
          <Entry value="2" text="vector_components"/>
          <Entry value="3" text="off"/>
        </Domain>
      </Property>
      <Property name="Seed" id="7160.Seed" number_of_elements="1">
        <Element index="0" value="10339"/>
        <Domain name="range" id="7160.Seed.range"/>
      </Property>
      <Property name="Source" id="7160.Source" number_of_elements="1">
        <Proxy value="7138" output_port="0"/>
        <Domain name="groups" id="7160.Source.groups"/>
        <Domain name="input_type" id="7160.Source.input_type"/>
        <Domain name="proxy_list" id="7160.Source.proxy_list">
          <Proxy value="7083"/>
          <Proxy value="7094"/>
          <Proxy value="7105"/>
          <Proxy value="7116"/>
          <Proxy value="7127"/>
          <Proxy value="7138"/>
          <Proxy value="7149"/>
        </Domain>
      </Property>
      <Property name="Stride" id="7160.Stride" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7160.Stride.range"/>
      </Property>
      <Property name="Vectors" id="7160.Vectors" number_of_elements="5">
        <Element index="0" value="1"/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Angular position [rad]"/>
        <Domain name="array_list" id="7160.Vectors.array_list">
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Atmosphere stress [Pa]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Sum of forces [N]"/>
          <String text="Sum of torques [N*m]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="TransformFilter" id="6973" servers="1">
      <Property name="Input" id="6973.Input" number_of_elements="1">
        <Proxy value="4809" output_port="0"/>
        <Domain name="groups" id="6973.Input.groups"/>
        <Domain name="input_type" id="6973.Input.input_type"/>
      </Property>
      <Property name="Transform" id="6973.Transform" number_of_elements="1">
        <Proxy value="6972"/>
        <Domain name="groups" id="6973.Transform.groups"/>
        <Domain name="proxy_list" id="6973.Transform.proxy_list">
          <Proxy value="6972"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="TransformFilter" id="7417" servers="1">
      <Property name="Input" id="7417.Input" number_of_elements="1">
        <Proxy value="4831" output_port="0"/>
        <Domain name="groups" id="7417.Input.groups"/>
        <Domain name="input_type" id="7417.Input.input_type"/>
      </Property>
      <Property name="Transform" id="7417.Transform" number_of_elements="1">
        <Proxy value="7416"/>
        <Domain name="groups" id="7417.Transform.groups"/>
        <Domain name="proxy_list" id="7417.Transform.proxy_list">
          <Proxy value="7416"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="filters" type="Glyph" id="7604" servers="1">
      <Property name="GlyphMode" id="7604.GlyphMode" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="7604.GlyphMode.enum">
          <Entry value="0" text="All Points"/>
          <Entry value="1" text="Every Nth Point"/>
          <Entry value="2" text="Uniform Spatial Distribution"/>
        </Domain>
      </Property>
      <Property name="GlyphTransform" id="7604.GlyphTransform" number_of_elements="1">
        <Proxy value="7526"/>
        <Domain name="proxy_list" id="7604.GlyphTransform.proxy_list">
          <Proxy value="7526"/>
        </Domain>
      </Property>
      <Property name="Input" id="7604.Input" number_of_elements="1">
        <Proxy value="7417" output_port="0"/>
        <Domain name="groups" id="7604.Input.groups"/>
        <Domain name="input_array1" id="7604.Input.input_array1"/>
        <Domain name="input_array2" id="7604.Input.input_array2"/>
        <Domain name="input_type" id="7604.Input.input_type"/>
      </Property>
      <Property name="MaximumNumberOfSamplePoints" id="7604.MaximumNumberOfSamplePoints" number_of_elements="1">
        <Element index="0" value="5000"/>
        <Domain name="range" id="7604.MaximumNumberOfSamplePoints.range"/>
      </Property>
      <Property name="Orient" id="7604.Orient" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="7604.Orient.bool"/>
      </Property>
      <Property name="Scalars" id="7604.Scalars" number_of_elements="5">
        <Element index="0" value=""/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="e: Relative interface height [m]"/>
        <Domain name="array_list" id="7604.Scalars.array_list">
          <String text="e: Relative interface height [m]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
        </Domain>
      </Property>
      <Property name="ScaleFactor" id="7604.ScaleFactor" number_of_elements="1">
        <Element index="0" value="7500"/>
        <Domain name="bounds" id="7604.ScaleFactor.bounds"/>
        <Domain name="scalar_range" id="7604.ScaleFactor.scalar_range"/>
        <Domain name="vector_range" id="7604.ScaleFactor.vector_range"/>
      </Property>
      <Property name="ScaleMode" id="7604.ScaleMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="7604.ScaleMode.enum">
          <Entry value="0" text="scalar"/>
          <Entry value="1" text="vector"/>
          <Entry value="2" text="vector_components"/>
          <Entry value="3" text="off"/>
        </Domain>
      </Property>
      <Property name="Seed" id="7604.Seed" number_of_elements="1">
        <Element index="0" value="10339"/>
        <Domain name="range" id="7604.Seed.range"/>
      </Property>
      <Property name="Source" id="7604.Source" number_of_elements="1">
        <Proxy value="7527" output_port="0"/>
        <Domain name="groups" id="7604.Source.groups"/>
        <Domain name="input_type" id="7604.Source.input_type"/>
        <Domain name="proxy_list" id="7604.Source.proxy_list">
          <Proxy value="7527"/>
          <Proxy value="7538"/>
          <Proxy value="7549"/>
          <Proxy value="7560"/>
          <Proxy value="7571"/>
          <Proxy value="7582"/>
          <Proxy value="7593"/>
        </Domain>
      </Property>
      <Property name="Stride" id="7604.Stride" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="7604.Stride.range"/>
      </Property>
      <Property name="Vectors" id="7604.Vectors" number_of_elements="5">
        <Element index="0" value="1"/>
        <Element index="1" value=""/>
        <Element index="2" value=""/>
        <Element index="3" value="0"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Domain name="array_list" id="7604.Vectors.array_list">
          <String text="Velocity vector [m/s]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="sources" type="XMLStructuredGridReader" id="4765" servers="1">
""")
        write(f, """<Property name="FileName" id="4765.FileName" """ *
              """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" """ *
                  """value="$(vtk_folder)/""" * 
                  """$(simulation.id).atmosphere.$(i).vts"/>\n""")
        end
        write(f, """<Domain name="files" id="4765.FileName.files"/>
      </Property>
      <Property name="FileNameInfo" id="4765.FileNameInfo" """ * 
      """number_of_elements="1">\n""")
        write(f, """<Element index="0" """ *
              """value="$(vtk_folder)/$(simulation.id).atmosphere.1.vts"/>
      </Property>
      <Property name="TimestepValues" id="4765.TimestepValues" """ * 
      """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" value="$(i - 1)"/>\n""")
        end
        write(f, """ </Property>
      <Property name="CellArrayInfo" id="4765.CellArrayInfo"/>
      <Property name="CellArrayStatus" id="4765.CellArrayStatus">
        <Domain name="array_list" id="4765.CellArrayStatus.array_list"/>
      </Property>
      <Property name="PointArrayInfo" id="4765.PointArrayInfo" number_of_elements="6">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="v: Meridional velocity [m/s]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Element index="5" value="1"/>
      </Property>
      <Property name="PointArrayStatus" id="4765.PointArrayStatus" number_of_elements="6">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="v: Meridional velocity [m/s]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Element index="5" value="1"/>
        <Domain name="array_list" id="4765.PointArrayStatus.array_list">
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="Velocity vector [m/s]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="sources" type="XMLPolyDataReader" id="4787" servers="1">""")
        write(f, """<Property name="FileName" id="4787.FileName" """ *
              """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" """ *
                  """value="$(vtk_folder)/""" * 
                  """$(simulation.id).icefloe-interaction.$(i).vtp"/>\n""")
        end
        write(f, """<Domain name="files" id="4787.FileName.files"/>
      </Property>
      <Property name="FileNameInfo" id="4787.FileNameInfo" """ * 
      """number_of_elements="1">\n""")
        write(f, """<Element index="0" 
              value="$(vtk_folder)/$(simulation.id).icefloe-interaction.1.vtp"/>
      </Property>
      <Property name="TimestepValues" id="4787.TimestepValues" """ * 
      """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" value="$(i - 1)"/>\n""")
        end
        write(f, """ </Property>
      <Property name="CellArrayInfo" id="4787.CellArrayInfo" number_of_elements="16">
        <Element index="0" value="Inter-particle vector [m]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="Shear displacement [m]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Force [N]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="Effective radius [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="Contact area [m^2]"/>
        <Element index="9" value="1"/>
        <Element index="10" value="Contact stiffness [N/m]"/>
        <Element index="11" value="1"/>
        <Element index="12" value="Tensile stress [Pa]"/>
        <Element index="13" value="1"/>
        <Element index="14" value="Contact age [s]"/>
        <Element index="15" value="1"/>
      </Property>
      <Property name="CellArrayStatus" id="4787.CellArrayStatus" number_of_elements="16">
        <Element index="0" value="Inter-particle vector [m]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="Shear displacement [m]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Force [N]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="Effective radius [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="Contact area [m^2]"/>
        <Element index="9" value="1"/>
        <Element index="10" value="Contact stiffness [N/m]"/>
        <Element index="11" value="1"/>
        <Element index="12" value="Tensile stress [Pa]"/>
        <Element index="13" value="1"/>
        <Element index="14" value="Contact age [s]"/>
        <Element index="15" value="1"/>
        <Domain name="array_list" id="4787.CellArrayStatus.array_list">
          <String text="Inter-particle vector [m]"/>
          <String text="Shear displacement [m]"/>
          <String text="Force [N]"/>
          <String text="Effective radius [m]"/>
          <String text="Contact area [m^2]"/>
          <String text="Contact stiffness [N/m]"/>
          <String text="Tensile stress [Pa]"/>
          <String text="Contact age [s]"/>
        </Domain>
      </Property>
      <Property name="PointArrayInfo" id="4787.PointArrayInfo"/>
      <Property name="PointArrayStatus" id="4787.PointArrayStatus">
        <Domain name="array_list" id="4787.PointArrayStatus.array_list"/>
      </Property>
    </Proxy>
    <Proxy group="sources" type="XMLUnstructuredGridReader" id="4809" """ * 
    """servers="1">\n""")
        write(f, """<Property name="FileName" id="4809.FileName" """ *
              """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" """ *
                  """value="$(vtk_folder)/""" * 
                  """$(simulation.id).icefloes.$(i).vtu"/>\n""")
        end
        write(f, """<Domain name="files" id="4809.FileName.files"/>
      </Property>
      <Property name="FileNameInfo" id="4809.FileNameInfo" """ * 
      """number_of_elements="1">\n""")
        write(f, """<Element index="0" """ *
              """value="$(vtk_folder)/$(simulation.id).icefloes.1.vtu"/>
      </Property>
      <Property name="TimestepValues" id="4809.TimestepValues" """ * 
      """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" value="$(i - 1)"/>\n""")
        end
        write(f, """ </Property>
      <Property name="CellArrayInfo" id="4809.CellArrayInfo"/>
      <Property name="CellArrayStatus" id="4809.CellArrayStatus">
        <Domain name="array_list" id="4809.CellArrayStatus.array_list"/>
      </Property>
      <Property name="PointArrayInfo" id="4809.PointArrayInfo" number_of_elements="78">
        <Element index="0" value="Density [kg m^-3]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="Thickness [m]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Diameter (contact) [m]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="Diameter (areal) [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="Circumreference  [m]"/>
        <Element index="9" value="1"/>
        <Element index="10" value="Horizontal surface area [m^2]"/>
        <Element index="11" value="1"/>
        <Element index="12" value="Side surface area [m^2]"/>
        <Element index="13" value="1"/>
        <Element index="14" value="Volume [m^3]"/>
        <Element index="15" value="1"/>
        <Element index="16" value="Mass [kg]"/>
        <Element index="17" value="1"/>
        <Element index="18" value="Moment of inertia [kg m^2]"/>
        <Element index="19" value="1"/>
        <Element index="20" value="Linear velocity [m s^-1]"/>
        <Element index="21" value="1"/>
        <Element index="22" value="Linear acceleration [m s^-2]"/>
        <Element index="23" value="1"/>
        <Element index="24" value="Sum of forces [N]"/>
        <Element index="25" value="1"/>
        <Element index="26" value="Angular position [rad]"/>
        <Element index="27" value="1"/>
        <Element index="28" value="Angular velocity [rad s^-1]"/>
        <Element index="29" value="1"/>
        <Element index="30" value="Angular acceleration [rad s^-2]"/>
        <Element index="31" value="1"/>
        <Element index="32" value="Sum of torques [N*m]"/>
        <Element index="33" value="1"/>
        <Element index="34" value="Fixed in space [-]"/>
        <Element index="35" value="1"/>
        <Element index="36" value="Free to rotate [-]"/>
        <Element index="37" value="1"/>
        <Element index="38" value="Enabled [-]"/>
        <Element index="39" value="1"/>
        <Element index="40" value="Contact stiffness (normal) [N m^-1]"/>
        <Element index="41" value="1"/>
        <Element index="42" value="Contact stiffness (tangential) [N m^-1]"/>
        <Element index="43" value="1"/>
        <Element index="44" value="Contact viscosity (normal) [N m^-1 s]"/>
        <Element index="45" value="1"/>
        <Element index="46" value="Contact viscosity (tangential) [N m^-1 s]"/>
        <Element index="47" value="1"/>
        <Element index="48" value="Contact friction (static) [-]"/>
        <Element index="49" value="1"/>
        <Element index="50" value="Contact friction (dynamic) [-]"/>
        <Element index="51" value="1"/>
        <Element index="52" value="Young&#x27;s modulus [Pa]"/>
        <Element index="53" value="1"/>
        <Element index="54" value="Poisson&#x27;s ratio [-]"/>
        <Element index="55" value="1"/>
        <Element index="56" value="Tensile strength [Pa]"/>
        <Element index="57" value="1"/>
        <Element index="58" value="Compressive strength prefactor [m^0.5 Pa]"/>
        <Element index="59" value="1"/>
        <Element index="60" value="Ocean drag coefficient (vertical) [-]"/>
        <Element index="61" value="1"/>
        <Element index="62" value="Ocean drag coefficient (horizontal) [-]"/>
        <Element index="63" value="1"/>
        <Element index="64" value="Atmosphere drag coefficient (vertical) [-]"/>
        <Element index="65" value="1"/>
        <Element index="66" value="Atmosphere drag coefficient (horizontal) [-]"/>
        <Element index="67" value="1"/>
        <Element index="68" value="Contact pressure [Pa]"/>
        <Element index="69" value="1"/>
        <Element index="70" value="Number of contacts [-]"/>
        <Element index="71" value="1"/>
        <Element index="72" value="Granular stress [Pa]"/>
        <Element index="73" value="1"/>
        <Element index="74" value="Ocean stress [Pa]"/>
        <Element index="75" value="1"/>
        <Element index="76" value="Atmosphere stress [Pa]"/>
        <Element index="77" value="1"/>
      </Property>
      <Property name="PointArrayStatus" id="4809.PointArrayStatus" number_of_elements="78">
        <Element index="0" value="Density [kg m^-3]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="Thickness [m]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Diameter (contact) [m]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="Diameter (areal) [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="Circumreference  [m]"/>
        <Element index="9" value="1"/>
        <Element index="10" value="Horizontal surface area [m^2]"/>
        <Element index="11" value="1"/>
        <Element index="12" value="Side surface area [m^2]"/>
        <Element index="13" value="1"/>
        <Element index="14" value="Volume [m^3]"/>
        <Element index="15" value="1"/>
        <Element index="16" value="Mass [kg]"/>
        <Element index="17" value="1"/>
        <Element index="18" value="Moment of inertia [kg m^2]"/>
        <Element index="19" value="1"/>
        <Element index="20" value="Linear velocity [m s^-1]"/>
        <Element index="21" value="1"/>
        <Element index="22" value="Linear acceleration [m s^-2]"/>
        <Element index="23" value="1"/>
        <Element index="24" value="Sum of forces [N]"/>
        <Element index="25" value="1"/>
        <Element index="26" value="Angular position [rad]"/>
        <Element index="27" value="1"/>
        <Element index="28" value="Angular velocity [rad s^-1]"/>
        <Element index="29" value="1"/>
        <Element index="30" value="Angular acceleration [rad s^-2]"/>
        <Element index="31" value="1"/>
        <Element index="32" value="Sum of torques [N*m]"/>
        <Element index="33" value="1"/>
        <Element index="34" value="Fixed in space [-]"/>
        <Element index="35" value="1"/>
        <Element index="36" value="Free to rotate [-]"/>
        <Element index="37" value="1"/>
        <Element index="38" value="Enabled [-]"/>
        <Element index="39" value="1"/>
        <Element index="40" value="Contact stiffness (normal) [N m^-1]"/>
        <Element index="41" value="1"/>
        <Element index="42" value="Contact stiffness (tangential) [N m^-1]"/>
        <Element index="43" value="1"/>
        <Element index="44" value="Contact viscosity (normal) [N m^-1 s]"/>
        <Element index="45" value="1"/>
        <Element index="46" value="Contact viscosity (tangential) [N m^-1 s]"/>
        <Element index="47" value="1"/>
        <Element index="48" value="Contact friction (static) [-]"/>
        <Element index="49" value="1"/>
        <Element index="50" value="Contact friction (dynamic) [-]"/>
        <Element index="51" value="1"/>
        <Element index="52" value="Young&#x27;s modulus [Pa]"/>
        <Element index="53" value="1"/>
        <Element index="54" value="Poisson&#x27;s ratio [-]"/>
        <Element index="55" value="1"/>
        <Element index="56" value="Tensile strength [Pa]"/>
        <Element index="57" value="1"/>
        <Element index="58" value="Compressive strength prefactor [m^0.5 Pa]"/>
        <Element index="59" value="1"/>
        <Element index="60" value="Ocean drag coefficient (vertical) [-]"/>
        <Element index="61" value="1"/>
        <Element index="62" value="Ocean drag coefficient (horizontal) [-]"/>
        <Element index="63" value="1"/>
        <Element index="64" value="Atmosphere drag coefficient (vertical) [-]"/>
        <Element index="65" value="1"/>
        <Element index="66" value="Atmosphere drag coefficient (horizontal) [-]"/>
        <Element index="67" value="1"/>
        <Element index="68" value="Contact pressure [Pa]"/>
        <Element index="69" value="1"/>
        <Element index="70" value="Number of contacts [-]"/>
        <Element index="71" value="1"/>
        <Element index="72" value="Granular stress [Pa]"/>
        <Element index="73" value="1"/>
        <Element index="74" value="Ocean stress [Pa]"/>
        <Element index="75" value="1"/>
        <Element index="76" value="Atmosphere stress [Pa]"/>
        <Element index="77" value="1"/>
        <Domain name="array_list" id="4809.PointArrayStatus.array_list">
          <String text="Density [kg m^-3]"/>
          <String text="Thickness [m]"/>
          <String text="Diameter (contact) [m]"/>
          <String text="Diameter (areal) [m]"/>
          <String text="Circumreference  [m]"/>
          <String text="Horizontal surface area [m^2]"/>
          <String text="Side surface area [m^2]"/>
          <String text="Volume [m^3]"/>
          <String text="Mass [kg]"/>
          <String text="Moment of inertia [kg m^2]"/>
          <String text="Linear velocity [m s^-1]"/>
          <String text="Linear acceleration [m s^-2]"/>
          <String text="Sum of forces [N]"/>
          <String text="Angular position [rad]"/>
          <String text="Angular velocity [rad s^-1]"/>
          <String text="Angular acceleration [rad s^-2]"/>
          <String text="Sum of torques [N*m]"/>
          <String text="Fixed in space [-]"/>
          <String text="Free to rotate [-]"/>
          <String text="Enabled [-]"/>
          <String text="Contact stiffness (normal) [N m^-1]"/>
          <String text="Contact stiffness (tangential) [N m^-1]"/>
          <String text="Contact viscosity (normal) [N m^-1 s]"/>
          <String text="Contact viscosity (tangential) [N m^-1 s]"/>
          <String text="Contact friction (static) [-]"/>
          <String text="Contact friction (dynamic) [-]"/>
          <String text="Young&#x27;s modulus [Pa]"/>
          <String text="Poisson&#x27;s ratio [-]"/>
          <String text="Tensile strength [Pa]"/>
          <String text="Compressive strength prefactor [m^0.5 Pa]"/>
          <String text="Ocean drag coefficient (vertical) [-]"/>
          <String text="Ocean drag coefficient (horizontal) [-]"/>
          <String text="Atmosphere drag coefficient (vertical) [-]"/>
          <String text="Atmosphere drag coefficient (horizontal) [-]"/>
          <String text="Contact pressure [Pa]"/>
          <String text="Number of contacts [-]"/>
          <String text="Granular stress [Pa]"/>
          <String text="Ocean stress [Pa]"/>
          <String text="Atmosphere stress [Pa]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="sources" type="XMLStructuredGridReader" id="4831" """ * 
    """servers="1">\n""")

        write(f, """<Property name="FileName" id="4831.FileName" """ *
              """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" """ *
                  """value="$(vtk_folder)/""" * 
                  """$(simulation.id).ocean.$(i).vts"/>\n""")
        end
        write(f, """<Domain name="files" id="4831.FileName.files"/>
      </Property>
      <Property name="FileNameInfo" id="4831.FileNameInfo" """ * 
      """number_of_elements="1">\n""")
        write(f, """<Element index="0" """ *
              """value="$(vtk_folder)/$(simulation.id).ocean.1.vts"/>
      </Property>
      <Property name="TimestepValues" id="4831.TimestepValues" """ * 
      """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" value="$(i - 1)"/>\n""")
        end
        write(f, """</Property>
      <Property name="CellArrayInfo" id="4831.CellArrayInfo"/>
      <Property name="CellArrayStatus" id="4831.CellArrayStatus">
        <Domain name="array_list" id="4831.CellArrayStatus.array_list"/>
      </Property>
      <Property name="PointArrayInfo" id="4831.PointArrayInfo" number_of_elements="10">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="v: Meridional velocity [m/s]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="h: Layer thickness [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="e: Relative interface height [m]"/>
        <Element index="9" value="1"/>
      </Property>
      <Property name="PointArrayStatus" id="4831.PointArrayStatus" number_of_elements="10">
        <Element index="0" value="u: Zonal velocity [m/s]"/>
        <Element index="1" value="1"/>
        <Element index="2" value="v: Meridional velocity [m/s]"/>
        <Element index="3" value="1"/>
        <Element index="4" value="Velocity vector [m/s]"/>
        <Element index="5" value="1"/>
        <Element index="6" value="h: Layer thickness [m]"/>
        <Element index="7" value="1"/>
        <Element index="8" value="e: Relative interface height [m]"/>
        <Element index="9" value="1"/>
        <Domain name="array_list" id="4831.PointArrayStatus.array_list">
          <String text="u: Zonal velocity [m/s]"/>
          <String text="v: Meridional velocity [m/s]"/>
          <String text="Velocity vector [m/s]"/>
          <String text="h: Layer thickness [m]"/>
          <String text="e: Relative interface height [m]"/>
        </Domain>
      </Property>
    </Proxy>
    <Proxy group="misc" type="TimeKeeper" id="259" servers="16">
      <Property name="SuppressedTimeSources" id="259.SuppressedTimeSources"/>
      <Property name="Time" id="259.Time" number_of_elements="1">
        <Element index="0" value="410"/>
        <Domain name="range" id="259.Time.range"/>
      </Property>
      <Property name="TimeLabel" id="259.TimeLabel" number_of_elements="1">
        <Element index="0" value="Time"/>
      </Property>
      <Property name="TimeRange" id="259.TimeRange" number_of_elements="2">
        <Element index="0" value="0"/>
        <Element index="1" value="$(simulation.file_number)"/>
      </Property>
      <Property name="TimeSources" id="259.TimeSources" number_of_elements="14">
        <Proxy value="4765"/>
        <Proxy value="4787"/>
        <Proxy value="4809"/>
        <Proxy value="4831"/>
        <Proxy value="6326"/>
        <Proxy value="6603"/>
        <Proxy value="6701"/>
        <Proxy value="6794"/>
        <Proxy value="6883"/>
        <Proxy value="6973"/>
        <Proxy value="7160"/>
        <Proxy value="7327"/>
        <Proxy value="7417"/>
        <Proxy value="7604"/>
        </Property>\n""")
        write(f, """<Property name="TimestepValues" id="259.TimestepValues" """*
            """number_of_elements="$(simulation.file_number)">\n""")
        for i=1:simulation.file_number
            write(f, """<Element index="$(i - 1)" value="$(i - 1)"/>\n""")
        end
        write(f, """ </Property>
      <Property name="Views" id="259.Views" number_of_elements="1">
        <Proxy value="6135"/>
      </Property>
    </Proxy>
    <Proxy group="views" type="RenderView" id="6135" servers="21">
      <Property name="AlphaBitPlanes" id="6135.AlphaBitPlanes" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.AlphaBitPlanes.bool"/>
      </Property>
      <Property name="AxesGrid" id="6135.AxesGrid" number_of_elements="1">
        <Proxy value="5227"/>
        <Domain name="proxy_list" id="6135.AxesGrid.proxy_list">
          <Proxy value="5227"/>
        </Domain>
      </Property>
      <Property name="BackLightAzimuth" id="6135.BackLightAzimuth" number_of_elements="1">
        <Element index="0" value="110"/>
        <Domain name="range" id="6135.BackLightAzimuth.range"/>
      </Property>
      <Property name="BackLightElevation" id="6135.BackLightElevation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6135.BackLightElevation.range"/>
      </Property>
      <Property name="BackLightK:B Ratio" id="6135.BackLightK:B Ratio" number_of_elements="1">
        <Element index="0" value="3.5"/>
        <Domain name="range" id="6135.BackLightK:B Ratio.range"/>
      </Property>
      <Property name="BackLightWarmth" id="6135.BackLightWarmth" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6135.BackLightWarmth.range"/>
      </Property>
      <Property name="Background" id="6135.Background" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6135.Background.range"/>
      </Property>
      <Property name="Background2" id="6135.Background2" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0.165"/>
        <Domain name="range" id="6135.Background2.range"/>
      </Property>
      <Property name="BackgroundTexture" id="6135.BackgroundTexture">
        <Domain name="groups" id="6135.BackgroundTexture.groups"/>
      </Property>
      <Property name="CacheKey" id="6135.CacheKey" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="range" id="6135.CacheKey.range"/>
      </Property>
      <Property name="Camera2DManipulators" id="6135.Camera2DManipulators" number_of_elements="9">
        <Element index="0" value="1"/>
        <Element index="1" value="3"/>
        <Element index="2" value="2"/>
        <Element index="3" value="2"/>
        <Element index="4" value="2"/>
        <Element index="5" value="2"/>
        <Element index="6" value="3"/>
        <Element index="7" value="1"/>
        <Element index="8" value="4"/>
        <Domain name="enum" id="6135.Camera2DManipulators.enum">
          <Entry value="0" text="None"/>
          <Entry value="1" text="Pan"/>
          <Entry value="2" text="Zoom"/>
          <Entry value="3" text="Roll"/>
          <Entry value="4" text="Rotate"/>
        </Domain>
      </Property>
      <Property name="Camera3DManipulators" id="6135.Camera3DManipulators" number_of_elements="9">
        <Element index="0" value="4"/>
        <Element index="1" value="1"/>
        <Element index="2" value="2"/>
        <Element index="3" value="3"/>
        <Element index="4" value="4"/>
        <Element index="5" value="1"/>
        <Element index="6" value="2"/>
        <Element index="7" value="4"/>
        <Element index="8" value="2"/>
        <Domain name="enum" id="6135.Camera3DManipulators.enum">
          <Entry value="0" text="None"/>
          <Entry value="1" text="Pan"/>
          <Entry value="2" text="Zoom"/>
          <Entry value="3" text="Roll"/>
          <Entry value="4" text="Rotate"/>
          <Entry value="5" text="Multi-Rotate"/>
        </Domain>
      </Property>
      <Property name="CenterAxesVisibility" id="6135.CenterAxesVisibility" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.CenterAxesVisibility.bool"/>
      </Property>
      <Property name="CenterOfRotation" id="6135.CenterOfRotation" number_of_elements="3">
        <Element index="0" value="25000"/>
        <Element index="1" value="37500"/>
        <Element index="2" value="-375"/>
      </Property>
      <Property name="CollectGeometryThreshold" id="6135.CollectGeometryThreshold" number_of_elements="1">
        <Element index="0" value="100"/>
      </Property>
      <Property name="CompressorConfig" id="6135.CompressorConfig" number_of_elements="1">
        <Element index="0" value="vtkSquirtCompressor 0 3"/>
      </Property>
      <Property name="DepthPeeling" id="6135.DepthPeeling" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.DepthPeeling.bool"/>
      </Property>
      <Property name="EnableRenderOnInteraction" id="6135.EnableRenderOnInteraction" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.EnableRenderOnInteraction.bool"/>
      </Property>
      <Property name="FillLightAzimuth" id="6135.FillLightAzimuth" number_of_elements="1">
        <Element index="0" value="-10"/>
        <Domain name="range" id="6135.FillLightAzimuth.range"/>
      </Property>
      <Property name="FillLightElevation" id="6135.FillLightElevation" number_of_elements="1">
        <Element index="0" value="-75"/>
        <Domain name="range" id="6135.FillLightElevation.range"/>
      </Property>
      <Property name="FillLightK:F Ratio" id="6135.FillLightK:F Ratio" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="range" id="6135.FillLightK:F Ratio.range"/>
      </Property>
      <Property name="FillLightWarmth" id="6135.FillLightWarmth" number_of_elements="1">
        <Element index="0" value="0.4"/>
        <Domain name="range" id="6135.FillLightWarmth.range"/>
      </Property>
      <Property name="HeadLightK:H Ratio" id="6135.HeadLightK:H Ratio" number_of_elements="1">
        <Element index="0" value="3"/>
        <Domain name="range" id="6135.HeadLightK:H Ratio.range"/>
      </Property>
      <Property name="HeadLightWarmth" id="6135.HeadLightWarmth" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6135.HeadLightWarmth.range"/>
      </Property>
      <Property name="ImageReductionFactor" id="6135.ImageReductionFactor" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6135.ImageReductionFactor.range"/>
      </Property>
      <Property name="InteractionMode" id="6135.InteractionMode" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6135.InteractionMode.enum">
          <Entry value="0" text="3D"/>
          <Entry value="1" text="2D"/>
          <Entry value="2" text="Selection"/>
        </Domain>
      </Property>
      <Property name="KeyLightAzimuth" id="6135.KeyLightAzimuth" number_of_elements="1">
        <Element index="0" value="10"/>
        <Domain name="range" id="6135.KeyLightAzimuth.range"/>
      </Property>
      <Property name="KeyLightElevation" id="6135.KeyLightElevation" number_of_elements="1">
        <Element index="0" value="50"/>
        <Domain name="range" id="6135.KeyLightElevation.range"/>
      </Property>
      <Property name="KeyLightIntensity" id="6135.KeyLightIntensity" number_of_elements="1">
        <Element index="0" value="0.75"/>
        <Domain name="range" id="6135.KeyLightIntensity.range"/>
      </Property>
      <Property name="KeyLightWarmth" id="6135.KeyLightWarmth" number_of_elements="1">
        <Element index="0" value="0.6"/>
        <Domain name="range" id="6135.KeyLightWarmth.range"/>
      </Property>
      <Property name="LODResolution" id="6135.LODResolution" number_of_elements="1">
        <Element index="0" value="0.5"/>
        <Domain name="range" id="6135.LODResolution.range"/>
      </Property>
      <Property name="LODThreshold" id="6135.LODThreshold" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="6135.LODThreshold.range"/>
      </Property>
      <Property name="LightAmbientColor" id="6135.LightAmbientColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6135.LightAmbientColor.range"/>
      </Property>
      <Property name="LightDiffuseColor" id="6135.LightDiffuseColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6135.LightDiffuseColor.range"/>
      </Property>
      <Property name="LightIntensity" id="6135.LightIntensity" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6135.LightIntensity.range"/>
      </Property>
      <Property name="LightSpecularColor" id="6135.LightSpecularColor" number_of_elements="3">
        <Element index="0" value="1"/>
        <Element index="1" value="1"/>
        <Element index="2" value="1"/>
        <Domain name="range" id="6135.LightSpecularColor.range"/>
      </Property>
      <Property name="LightSwitch" id="6135.LightSwitch" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.LightSwitch.bool"/>
      </Property>
      <Property name="LightType" id="6135.LightType" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="enum" id="6135.LightType.enum">
          <Entry value="1" text="HeadLight"/>
          <Entry value="2" text="CameraLight"/>
          <Entry value="3" text="SceneLight"/>
        </Domain>
      </Property>
      <Property name="MaintainLuminance" id="6135.MaintainLuminance" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.MaintainLuminance.bool"/>
      </Property>
      <Property name="MaximumNumberOfPeels" id="6135.MaximumNumberOfPeels" number_of_elements="1">
        <Element index="0" value="4"/>
        <Domain name="range" id="6135.MaximumNumberOfPeels.range"/>
      </Property>
      <Property name="MultiSamples" id="6135.MultiSamples" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
      <Property name="NonInteractiveRenderDelay" id="6135.NonInteractiveRenderDelay" number_of_elements="1">
        <Element index="0" value="0"/>
      </Property>
      <Property name="OrientationAxesInteractivity" id="6135.OrientationAxesInteractivity" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.OrientationAxesInteractivity.bool"/>
      </Property>
      <Property name="OrientationAxesLabelColor" id="6135.OrientationAxesLabelColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="OrientationAxesOutlineColor" id="6135.OrientationAxesOutlineColor" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="0"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="OrientationAxesVisibility" id="6135.OrientationAxesVisibility" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.OrientationAxesVisibility.bool"/>
      </Property>
      <Property name="RemoteRenderThreshold" id="6135.RemoteRenderThreshold" number_of_elements="1">
        <Element index="0" value="20"/>
        <Domain name="range" id="6135.RemoteRenderThreshold.range"/>
      </Property>
      <Property name="Representations" id="6135.Representations" number_of_elements="26">
        <Proxy value="5235"/>
        <Proxy value="5243"/>
        <Proxy value="5251"/>
        <Proxy value="5259"/>
        <Proxy value="5267"/>
        <Proxy value="5275"/>
        <Proxy value="5360"/>
        <Proxy value="5453"/>
        <Proxy value="5546"/>
        <Proxy value="5623"/>
        <Proxy value="5639"/>
        <Proxy value="5803"/>
        <Proxy value="5809"/>
        <Proxy value="5894"/>
        <Proxy value="6133"/>
        <Proxy value="6424"/>
        <Proxy value="6681"/>
        <Proxy value="6699"/>
        <Proxy value="6783"/>
        <Proxy value="6872"/>
        <Proxy value="6961"/>
        <Proxy value="7071"/>
        <Proxy value="7238"/>
        <Proxy value="7405"/>
        <Proxy value="7515"/>
        <Proxy value="7682"/>
      </Property>
      <Property name="RotationFactor" id="6135.RotationFactor" number_of_elements="1">
        <Element index="0" value="1"/>
      </Property>
      <Property name="ServerStereoType" id="6135.ServerStereoType" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6135.ServerStereoType.enum">
          <Entry value="0" text="Same As Client"/>
          <Entry value="1" text="Crystal Eyes"/>
          <Entry value="2" text="Red-Blue"/>
          <Entry value="3" text="Interlaced"/>
          <Entry value="4" text="Left"/>
          <Entry value="5" text="Right"/>
          <Entry value="6" text="Dresden"/>
          <Entry value="7" text="Anaglyph"/>
          <Entry value="8" text="Checkerboard"/>
          <Entry value="9" text="SplitViewportHorizontal"/>
          <Entry value="10" text="None"/>
        </Domain>
      </Property>
      <Property name="ShowAnnotation" id="6135.ShowAnnotation" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.ShowAnnotation.bool"/>
      </Property>
      <Property name="StencilCapable" id="6135.StencilCapable" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.StencilCapable.bool"/>
      </Property>
      <Property name="StereoRender" id="6135.StereoRender" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.StereoRender.bool"/>
      </Property>
      <Property name="StereoType" id="6135.StereoType" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="enum" id="6135.StereoType.enum">
          <Entry value="1" text="Crystal Eyes"/>
          <Entry value="2" text="Red-Blue"/>
          <Entry value="3" text="Interlaced"/>
          <Entry value="4" text="Left"/>
          <Entry value="5" text="Right"/>
          <Entry value="6" text="Dresden"/>
          <Entry value="7" text="Anaglyph"/>
          <Entry value="8" text="Checkerboard"/>
          <Entry value="9" text="SplitViewportHorizontal"/>
          <Entry value="10" text="None"/>
        </Domain>
      </Property>
      <Property name="StillRenderImageReductionFactor" id="6135.StillRenderImageReductionFactor" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="range" id="6135.StillRenderImageReductionFactor.range"/>
      </Property>
      <Property name="UseCache" id="6135.UseCache" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseCache.bool"/>
      </Property>
      <Property name="UseGradientBackground" id="6135.UseGradientBackground" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseGradientBackground.bool"/>
      </Property>
      <Property name="UseInteractiveRenderingForScreenshots" id="6135.UseInteractiveRenderingForScreenshots" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseInteractiveRenderingForScreenshots.bool"/>
      </Property>
      <Property name="UseLight" id="6135.UseLight" number_of_elements="1">
        <Element index="0" value="1"/>
        <Domain name="bool" id="6135.UseLight.bool"/>
      </Property>
      <Property name="UseOffscreenRendering" id="6135.UseOffscreenRendering" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseOffscreenRendering.bool"/>
      </Property>
      <Property name="UseOffscreenRenderingForScreenshots" id="6135.UseOffscreenRenderingForScreenshots" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseOffscreenRenderingForScreenshots.bool"/>
      </Property>
      <Property name="UseOutlineForLODRendering" id="6135.UseOutlineForLODRendering" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseOutlineForLODRendering.bool"/>
      </Property>
      <Property name="UseTexturedBackground" id="6135.UseTexturedBackground" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.UseTexturedBackground.bool"/>
      </Property>
      <Property name="ViewSize" id="6135.ViewSize" number_of_elements="2">
        <Element index="0" value="983"/>
        <Element index="1" value="929"/>
      </Property>
      <Property name="ViewTime" id="6135.ViewTime" number_of_elements="1">
        <Element index="0" value="410"/>
        <Domain name="range" id="6135.ViewTime.range"/>
      </Property>
      <Property name="CameraClippingRange" id="6135.CameraClippingRange" number_of_elements="2">
        <Element index="0" value="169422.800869707"/>
        <Element index="1" value="179651.834875773"/>
      </Property>
      <Property name="CameraClippingRangeInfo" id="6135.CameraClippingRangeInfo" number_of_elements="2">
        <Element index="0" value="169422.800869707"/>
        <Element index="1" value="179651.834875773"/>
      </Property>
      <Property name="CameraFocalPoint" id="6135.CameraFocalPoint" number_of_elements="3">
        <Element index="0" value="28327.3614847311"/>
        <Element index="1" value="47107.885941104"/>
        <Element index="2" value="-375"/>
      </Property>
      <Property name="CameraFocalPointInfo" id="6135.CameraFocalPointInfo" number_of_elements="3">
        <Element index="0" value="28327.3614847311"/>
        <Element index="1" value="47107.885941104"/>
        <Element index="2" value="-375"/>
      </Property>
      <Property name="CameraParallelProjection" id="6135.CameraParallelProjection" number_of_elements="1">
        <Element index="0" value="0"/>
        <Domain name="bool" id="6135.CameraParallelProjection.bool"/>
      </Property>
      <Property name="CameraParallelScale" id="6135.CameraParallelScale" number_of_elements="1">
        <Element index="0" value="54535.8507228272"/>
      </Property>
      <Property name="CameraParallelScaleInfo" id="6135.CameraParallelScaleInfo" number_of_elements="1">
        <Element index="0" value="54535.8507228272"/>
      </Property>
      <Property name="CameraPosition" id="6135.CameraPosition" number_of_elements="3">
        <Element index="0" value="28327.3614847311"/>
        <Element index="1" value="47107.885941104"/>
        <Element index="2" value="173765.782386196"/>
      </Property>
      <Property name="CameraPositionInfo" id="6135.CameraPositionInfo" number_of_elements="3">
        <Element index="0" value="28327.3614847311"/>
        <Element index="1" value="47107.885941104"/>
        <Element index="2" value="173765.782386196"/>
      </Property>
      <Property name="CameraViewAngle" id="6135.CameraViewAngle" number_of_elements="1">
        <Element index="0" value="30"/>
      </Property>
      <Property name="CameraViewUp" id="6135.CameraViewUp" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="CameraViewUpInfo" id="6135.CameraViewUpInfo" number_of_elements="3">
        <Element index="0" value="0"/>
        <Element index="1" value="1"/>
        <Element index="2" value="0"/>
      </Property>
      <Property name="EyeAngle" id="6135.EyeAngle" number_of_elements="1">
        <Element index="0" value="2"/>
        <Domain name="range" id="6135.EyeAngle.range"/>
      </Property>
    </Proxy>
    <ProxyCollection name="animation">
      <Item id="263" name="AnimationScene1"/>
      <Item id="265" name="TimeAnimationCue1"/>
    </ProxyCollection>
    <ProxyCollection name="layouts">
      <Item id="6136" name="ViewLayout1"/>
    </ProxyCollection>
    <ProxyCollection name="lookup_tables">
      <Item id="5797" name="Atmospheredragcoefficienthorizontal.PVLookupTable"/>
      <Item id="5229" name="Diameterarealm.PVLookupTable"/>
      <Item id="5261" name="Enabled.PVLookupTable"/>
      <Item id="5253" name="Fixedinspace.PVLookupTable"/>
      <Item id="5245" name="Linearaccelerationms2.PVLookupTable"/>
      <Item id="5740" name="Linearvelocityms1.PVLookupTable"/>
      <Item id="5237" name="Numberofcontacts.PVLookupTable"/>
      <Item id="6072" name="TensilestressPa.PVLookupTable"/>
      <Item id="5299" name="Velocityvectorms.PVLookupTable"/>
      <Item id="5269" name="eRelativeinterfaceheightm.PVLookupTable"/>
      <Item id="5888" name="uZonalvelocityms.PVLookupTable"/>
      <Item id="6693" name="vMeridionalvelocityms.PVLookupTable"/>
    </ProxyCollection>
    <ProxyCollection name="piecewise_functions">
      <Item id="5796" name="Atmospheredragcoefficienthorizontal.PiecewiseFunction"/>
      <Item id="5228" name="Diameterarealm.PiecewiseFunction"/>
      <Item id="5260" name="Enabled.PiecewiseFunction"/>
      <Item id="5252" name="Fixedinspace.PiecewiseFunction"/>
      <Item id="5244" name="Linearaccelerationms2.PiecewiseFunction"/>
      <Item id="5739" name="Linearvelocityms1.PiecewiseFunction"/>
      <Item id="5236" name="Numberofcontacts.PiecewiseFunction"/>
      <Item id="6071" name="TensilestressPa.PiecewiseFunction"/>
      <Item id="5298" name="Velocityvectorms.PiecewiseFunction"/>
      <Item id="5268" name="eRelativeinterfaceheightm.PiecewiseFunction"/>
      <Item id="5887" name="uZonalvelocityms.PiecewiseFunction"/>
      <Item id="6692" name="vMeridionalvelocityms.PiecewiseFunction"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.11810">
      <Item id="5227" name="AxesGrid"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17038">
      <Item id="6182" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17106">
      <Item id="6183" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17174">
      <Item id="6184" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17198">
      <Item id="6185" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17652">
      <Item id="4842" name="GlyphTransform"/>
      <Item id="6186" name="RepresentationAnimationHelper"/>
      <Item id="4843" name="Source"/>
      <Item id="4854" name="Source"/>
      <Item id="4865" name="Source"/>
      <Item id="4876" name="Source"/>
      <Item id="4887" name="Source"/>
      <Item id="4898" name="Source"/>
      <Item id="4909" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.17821">
      <Item id="4931" name="GlyphTransform"/>
      <Item id="6187" name="RepresentationAnimationHelper"/>
      <Item id="4932" name="Source"/>
      <Item id="4943" name="Source"/>
      <Item id="4954" name="Source"/>
      <Item id="4965" name="Source"/>
      <Item id="4976" name="Source"/>
      <Item id="4987" name="Source"/>
      <Item id="4998" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.18006">
      <Item id="5020" name="GlyphTransform"/>
      <Item id="6188" name="RepresentationAnimationHelper"/>
      <Item id="5021" name="Source"/>
      <Item id="5032" name="Source"/>
      <Item id="5043" name="Source"/>
      <Item id="5054" name="Source"/>
      <Item id="5065" name="Source"/>
      <Item id="5076" name="Source"/>
      <Item id="5087" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.18352">
      <Item id="5109" name="GlyphTransform"/>
      <Item id="6189" name="RepresentationAnimationHelper"/>
      <Item id="5110" name="Source"/>
      <Item id="5121" name="Source"/>
      <Item id="5132" name="Source"/>
      <Item id="5143" name="Source"/>
      <Item id="5154" name="Source"/>
      <Item id="5165" name="Source"/>
      <Item id="5176" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.18443">
      <Item id="6190" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.18532">
      <Item id="6191" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6326">
      <Item id="6337" name="RepresentationAnimationHelper"/>
      <Item id="6325" name="Transform"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6603">
      <Item id="6525" name="GlyphTransform"/>
      <Item id="6614" name="RepresentationAnimationHelper"/>
      <Item id="6526" name="Source"/>
      <Item id="6537" name="Source"/>
      <Item id="6548" name="Source"/>
      <Item id="6559" name="Source"/>
      <Item id="6570" name="Source"/>
      <Item id="6581" name="Source"/>
      <Item id="6592" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6701">
      <Item id="6712" name="RepresentationAnimationHelper"/>
      <Item id="6700" name="Transform"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6794">
      <Item id="6805" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6883">
      <Item id="6894" name="RepresentationAnimationHelper"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.6973">
      <Item id="6984" name="RepresentationAnimationHelper"/>
      <Item id="6972" name="Transform"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.7160">
      <Item id="7082" name="GlyphTransform"/>
      <Item id="7171" name="RepresentationAnimationHelper"/>
      <Item id="7083" name="Source"/>
      <Item id="7094" name="Source"/>
      <Item id="7105" name="Source"/>
      <Item id="7116" name="Source"/>
      <Item id="7127" name="Source"/>
      <Item id="7138" name="Source"/>
      <Item id="7149" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.7327">
      <Item id="7249" name="GlyphTransform"/>
      <Item id="7338" name="RepresentationAnimationHelper"/>
      <Item id="7250" name="Source"/>
      <Item id="7261" name="Source"/>
      <Item id="7272" name="Source"/>
      <Item id="7283" name="Source"/>
      <Item id="7294" name="Source"/>
      <Item id="7305" name="Source"/>
      <Item id="7316" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.7417">
      <Item id="7428" name="RepresentationAnimationHelper"/>
      <Item id="7416" name="Transform"/>
    </ProxyCollection>
    <ProxyCollection name="pq_helper_proxies.7604">
      <Item id="7526" name="GlyphTransform"/>
      <Item id="7615" name="RepresentationAnimationHelper"/>
      <Item id="7527" name="Source"/>
      <Item id="7538" name="Source"/>
      <Item id="7549" name="Source"/>
      <Item id="7560" name="Source"/>
      <Item id="7571" name="Source"/>
      <Item id="7582" name="Source"/>
      <Item id="7593" name="Source"/>
    </ProxyCollection>
    <ProxyCollection name="representations">
      <Item id="5623" name="GeometryRepresentation1"/>
      <Item id="7682" name="GeometryRepresentation2"/>
      <Item id="7238" name="GeometryRepresentation3"/>
      <Item id="7405" name="GeometryRepresentation4"/>
      <Item id="6783" name="GeometryRepresentation5"/>
      <Item id="6872" name="GeometryRepresentation6"/>
      <Item id="6961" name="GeometryRepresentation7"/>
      <Item id="6681" name="GeometryRepresentation8"/>
      <Item id="5360" name="StructuredGridRepresentation1"/>
      <Item id="5546" name="StructuredGridRepresentation2"/>
      <Item id="6424" name="StructuredGridRepresentation3"/>
      <Item id="7515" name="StructuredGridRepresentation4"/>
      <Item id="5453" name="UnstructuredGridRepresentation1"/>
      <Item id="7071" name="UnstructuredGridRepresentation2"/>
    </ProxyCollection>
    <ProxyCollection name="scalar_bars">
      <Item id="5235" name="ScalarBarWidgetRepresentation1"/>
      <Item id="5894" name="ScalarBarWidgetRepresentation10"/>
      <Item id="6133" name="ScalarBarWidgetRepresentation11"/>
      <Item id="6699" name="ScalarBarWidgetRepresentation12"/>
      <Item id="5243" name="ScalarBarWidgetRepresentation2"/>
      <Item id="5251" name="ScalarBarWidgetRepresentation3"/>
      <Item id="5259" name="ScalarBarWidgetRepresentation4"/>
      <Item id="5267" name="ScalarBarWidgetRepresentation5"/>
      <Item id="5275" name="ScalarBarWidgetRepresentation6"/>
      <Item id="5639" name="ScalarBarWidgetRepresentation7"/>
      <Item id="5803" name="ScalarBarWidgetRepresentation8"/>
      <Item id="5809" name="ScalarBarWidgetRepresentation9"/>
    </ProxyCollection>
    <ProxyCollection name="sources">
      <Item id="6326" name="Atmosphere Transform"/>
      <Item id="6603" name="Atmosphere Velocity Vectors"/>
      <Item id="6794" name="CellDatatoPointData"/>
      <Item id="6701" name="Icefloe Interaction Transform"/>
      <Item id="6883" name="Icefloe Interaction Tube"/>
      <Item id="7327" name="Icefloe Velocity Vectors"/>
      <Item id="7160" name="Icefloes"/>
      <Item id="6973" name="Icefloes Transform"/>
      <Item id="7417" name="Ocean Transform"/>
      <Item id="7604" name="Ocean Velocity Vectors"/>\n""")
        write(f, """ <Item id="4765" name="$(simulation.id).atmosphere.*"/>
              <Item id="4787" name="$(simulation.id).icefloe-interaction.*"/>
              <Item id="4809" name="$(simulation.id).icefloes.*"/>
              <Item id="4831" name="$(simulation.id).ocean.*"/>\n""")
        write(f, """ <Item id="259" name="TimeKeeper1"/>
    </ProxyCollection>
    <ProxyCollection name="views">
      <Item id="6135" name="RenderView1"/>
    </ProxyCollection>
    <CustomProxyDefinitions>
      <CustomProxyDefinition name="SeaIce" group="filters">
        <CompoundSourceProxy id="6702" servers="1">
          <Proxy group="filters" type="Glyph" id="6613" servers="1" compound_name="Glyph1">
            <Property name="GlyphMode" id="6613.GlyphMode" number_of_elements="1">
              <Element index="0" value="0"/>
              <Domain name="enum" id="6613.GlyphMode.enum">
                <Entry value="0" text="All Points"/>
                <Entry value="1" text="Every Nth Point"/>
                <Entry value="2" text="Uniform Spatial Distribution"/>
              </Domain>
            </Property>
            <Property name="GlyphTransform" id="6613.GlyphTransform" number_of_elements="1">
              <Proxy value="6535"/>
              <Domain name="proxy_list" id="6613.GlyphTransform.proxy_list">
                <Proxy value="6535"/>
              </Domain>
            </Property>
            <Property name="Input" id="6613.Input" number_of_elements="1">
              <Domain name="groups" id="6613.Input.groups"/>
              <Domain name="input_array1" id="6613.Input.input_array1"/>
              <Domain name="input_array2" id="6613.Input.input_array2"/>
              <Domain name="input_type" id="6613.Input.input_type"/>
            </Property>
            <Property name="MaximumNumberOfSamplePoints" id="6613.MaximumNumberOfSamplePoints" number_of_elements="1">
              <Element index="0" value="5000"/>
              <Domain name="range" id="6613.MaximumNumberOfSamplePoints.range"/>
            </Property>
            <Property name="Orient" id="6613.Orient" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="bool" id="6613.Orient.bool"/>
            </Property>
            <Property name="Scalars" id="6613.Scalars" number_of_elements="5">
              <Element index="0" value=""/>
              <Element index="1" value=""/>
              <Element index="2" value=""/>
              <Element index="3" value="0"/>
              <Element index="4" value="Diameter (areal) [m]"/>
              <Domain name="array_list" id="6613.Scalars.array_list">
                <String text="Contact friction (dynamic) [-]"/>
                <String text="Contact friction (static) [-]"/>
                <String text="Contact stiffness (normal) [N m^-1]"/>
                <String text="Contact stiffness (tangential) [N m^-1]"/>
                <String text="Contact viscosity (normal) [N m^-1 s]"/>
                <String text="Contact viscosity (tangential) [N m^-1 s]"/>
                <String text="Density [kg m^-3]"/>
                <String text="Diameter (areal) [m]"/>
                <String text="Diameter (contact) [m]"/>
                <String text="Fixed in space [-]"/>
                <String text="Free to rotate [-]"/>
                <String text="Mass [kg]"/>
                <String text="Moment of inertia [kg m^2]"/>
                <String text="Surface area [m^2]"/>
                <String text="Thickness [m]"/>
                <String text="Volume [m^3]"/>
              </Domain>
            </Property>
            <Property name="ScaleFactor" id="6613.ScaleFactor" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="bounds" id="6613.ScaleFactor.bounds"/>
              <Domain name="scalar_range" id="6613.ScaleFactor.scalar_range"/>
              <Domain name="vector_range" id="6613.ScaleFactor.vector_range"/>
            </Property>
            <Property name="ScaleMode" id="6613.ScaleMode" number_of_elements="1">
              <Element index="0" value="0"/>
              <Domain name="enum" id="6613.ScaleMode.enum">
                <Entry value="0" text="scalar"/>
                <Entry value="1" text="vector"/>
                <Entry value="2" text="vector_components"/>
                <Entry value="3" text="off"/>
              </Domain>
            </Property>
            <Property name="Seed" id="6613.Seed" number_of_elements="1">
              <Element index="0" value="10339"/>
              <Domain name="range" id="6613.Seed.range"/>
            </Property>
            <Property name="Source" id="6613.Source" number_of_elements="1">
              <Proxy value="6591" output_port="0"/>
              <Domain name="groups" id="6613.Source.groups"/>
              <Domain name="input_type" id="6613.Source.input_type"/>
              <Domain name="proxy_list" id="6613.Source.proxy_list">
                <Proxy value="6536"/>
                <Proxy value="6547"/>
                <Proxy value="6558"/>
                <Proxy value="6569"/>
                <Proxy value="6580"/>
                <Proxy value="6591"/>
                <Proxy value="6602"/>
              </Domain>
            </Property>
            <Property name="Stride" id="6613.Stride" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="range" id="6613.Stride.range"/>
            </Property>
            <Property name="Vectors" id="6613.Vectors" number_of_elements="5">
              <Element index="0" value="1"/>
              <Element index="1" value=""/>
              <Element index="2" value=""/>
              <Element index="3" value="0"/>
              <Element index="4" value="Angular position [rad]"/>
              <Domain name="array_list" id="6613.Vectors.array_list">
                <String text="Angular acceleration [rad s^-2]"/>
                <String text="Angular position [rad]"/>
                <String text="Angular velocity [rad s^-1]"/>
                <String text="Linear acceleration [m s^-2]"/>
                <String text="Linear velocity [m s^-1]"/>
                <String text="Sum of forces [N]"/>
                <String text="Sum of torques [N*m]"/>
              </Domain>
            </Property>
          </Proxy>
          <Proxy group="extended_sources" type="Transform2" id="6535" servers="1" compound_name="auto_6535">
            <Property name="Position" id="6535.Position" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="0"/>
              <Domain name="range" id="6535.Position.range"/>
            </Property>
            <Property name="PositionInfo" id="6535.PositionInfo" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="0"/>
            </Property>
            <Property name="Rotation" id="6535.Rotation" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="0"/>
              <Domain name="range" id="6535.Rotation.range"/>
            </Property>
            <Property name="RotationInfo" id="6535.RotationInfo" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="0"/>
            </Property>
            <Property name="Scale" id="6535.Scale" number_of_elements="3">
              <Element index="0" value="1"/>
              <Element index="1" value="1"/>
              <Element index="2" value="1"/>
              <Domain name="range" id="6535.Scale.range"/>
            </Property>
            <Property name="ScaleInfo" id="6535.ScaleInfo" number_of_elements="3">
              <Element index="0" value="1"/>
              <Element index="1" value="1"/>
              <Element index="2" value="1"/>
            </Property>
          </Proxy>
          <Proxy group="sources" type="SphereSource" id="6591" servers="1" compound_name="auto_6591">
            <Property name="Center" id="6591.Center" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="0"/>
              <Domain name="range" id="6591.Center.range"/>
            </Property>
            <Property name="EndPhi" id="6591.EndPhi" number_of_elements="1">
              <Element index="0" value="180"/>
              <Domain name="range" id="6591.EndPhi.range"/>
            </Property>
            <Property name="EndTheta" id="6591.EndTheta" number_of_elements="1">
              <Element index="0" value="360"/>
              <Domain name="range" id="6591.EndTheta.range"/>
            </Property>
            <Property name="PhiResolution" id="6591.PhiResolution" number_of_elements="1">
              <Element index="0" value="8"/>
              <Domain name="range" id="6591.PhiResolution.range"/>
            </Property>
            <Property name="Radius" id="6591.Radius" number_of_elements="1">
              <Element index="0" value="0.5"/>
              <Domain name="range" id="6591.Radius.range"/>
            </Property>
            <Property name="StartPhi" id="6591.StartPhi" number_of_elements="1">
              <Element index="0" value="0"/>
              <Domain name="range" id="6591.StartPhi.range"/>
            </Property>
            <Property name="StartTheta" id="6591.StartTheta" number_of_elements="1">
              <Element index="0" value="0"/>
              <Domain name="range" id="6591.StartTheta.range"/>
            </Property>
            <Property name="ThetaResolution" id="6591.ThetaResolution" number_of_elements="1">
              <Element index="0" value="8"/>
              <Domain name="range" id="6591.ThetaResolution.range"/>
            </Property>
          </Proxy>
          <ExposedProperties>
            <Property name="GlyphMode" proxy_name="Glyph1" exposed_name="Glyph Mode"/>
            <Property name="GlyphTransform" proxy_name="Glyph1" exposed_name="Glyph Transform"/>
            <Property name="Input" proxy_name="Glyph1" exposed_name="Input"/>
            <Property name="MaximumNumberOfSamplePoints" proxy_name="Glyph1" exposed_name="Maximum Number Of Sample Points"/>
            <Property name="Orient" proxy_name="Glyph1" exposed_name="Orient"/>
            <Property name="Scalars" proxy_name="Glyph1" exposed_name="Scalars"/>
            <Property name="ScaleFactor" proxy_name="Glyph1" exposed_name="Scale Factor"/>
            <Property name="ScaleMode" proxy_name="Glyph1" exposed_name="Scale Mode"/>
            <Property name="Vectors" proxy_name="Glyph1" exposed_name="Vectors"/>
          </ExposedProperties>
          <OutputPort name="Output" proxy="Glyph1" port_index="0"/>
          <Hints>
            <ShowInMenu/>
          </Hints>
        </CompoundSourceProxy>
      </CustomProxyDefinition>
      <CustomProxyDefinition name="SeaIceGranularInteraction" group="filters">
        <CompoundSourceProxy id="6731" servers="1">
          <Proxy group="filters" type="TubeFilter" id="6634" servers="1" compound_name="Tube1">
            <Property name="Capping" id="6634.Capping" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="bool" id="6634.Capping.bool"/>
            </Property>
            <Property name="DefaultNormal" id="6634.DefaultNormal" number_of_elements="3">
              <Element index="0" value="0"/>
              <Element index="1" value="0"/>
              <Element index="2" value="1"/>
              <Domain name="range" id="6634.DefaultNormal.range"/>
            </Property>
            <Property name="Input" id="6634.Input" number_of_elements="1">
              <Domain name="groups" id="6634.Input.groups"/>
              <Domain name="input_array1" id="6634.Input.input_array1"/>
              <Domain name="input_array2" id="6634.Input.input_array2"/>
              <Domain name="input_type" id="6634.Input.input_type"/>
            </Property>
            <Property name="NumberOfSides" id="6634.NumberOfSides" number_of_elements="1">
              <Element index="0" value="6"/>
              <Domain name="range" id="6634.NumberOfSides.range"/>
            </Property>
            <Property name="Radius" id="6634.Radius" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="bounds" id="6634.Radius.bounds"/>
            </Property>
            <Property name="RadiusFactor" id="6634.RadiusFactor" number_of_elements="1">
              <Element index="0" value="250"/>
              <Domain name="range" id="6634.RadiusFactor.range"/>
            </Property>
            <Property name="SelectInputScalars" id="6634.SelectInputScalars" number_of_elements="5">
              <Element index="0" value=""/>
              <Element index="1" value=""/>
              <Element index="2" value=""/>
              <Element index="3" value="0"/>
              <Element index="4" value="Tensile stress [Pa]"/>
              <Domain name="array_list" id="6634.SelectInputScalars.array_list">
                <String text="Contact age [s]"/>
                <String text="Contact area [m^2]"/>
                <String text="Contact stiffness [N/m]"/>
                <String text="Effective radius [m]"/>
                <String text="Force [N]"/>
                <String text="Tensile stress [Pa]"/>
              </Domain>
            </Property>
            <Property name="SelectInputVectors" id="6634.SelectInputVectors" number_of_elements="5">
              <Element index="0" value="1"/>
              <Element index="1" value=""/>
              <Element index="2" value=""/>
              <Element index="3" value="0"/>
              <Element index="4" value="Inter-particle vector [m]"/>
              <Domain name="array_list" id="6634.SelectInputVectors.array_list">
                <String text="Inter-particle vector [m]"/>
                <String text="Shear displacement [m]"/>
              </Domain>
            </Property>
            <Property name="UseDefaultNormal" id="6634.UseDefaultNormal" number_of_elements="1">
              <Element index="0" value="0"/>
              <Domain name="bool" id="6634.UseDefaultNormal.bool"/>
            </Property>
            <Property name="VaryRadius" id="6634.VaryRadius" number_of_elements="1">
              <Element index="0" value="1"/>
              <Domain name="enum" id="6634.VaryRadius.enum">
                <Entry value="0" text="Off"/>
                <Entry value="1" text="By Scalar"/>
                <Entry value="2" text="By Vector"/>
                <Entry value="3" text="By Absolute Scalar"/>
              </Domain>
            </Property>
          </Proxy>
          <ExposedProperties>
            <Property name="Input" proxy_name="Tube1" exposed_name="Input"/>
          </ExposedProperties>
          <OutputPort name="Output" proxy="Tube1" port_index="0"/>
          <Hints>
            <ShowInMenu/>
          </Hints>
        </CompoundSourceProxy>
      </CustomProxyDefinition>
    </CustomProxyDefinitions>
    <Links>
      <GlobalPropertyLink global_name="BackgroundColor" proxy="6135" property="Background"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="5360" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="5453" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="5546" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="5623" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="6424" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="6681" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="6783" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="6872" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="6961" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="7071" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="7238" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="7405" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="7515" property="EdgeColor"/>
      <GlobalPropertyLink global_name="EdgeColor" proxy="7682" property="EdgeColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5227" property="GridColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6135" property="OrientationAxesOutlineColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5360" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5360" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5453" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5453" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5546" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5546" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5623" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="5623" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6424" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6424" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6681" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6681" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6783" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6783" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6872" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6872" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6961" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="6961" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7071" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7071" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7238" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7238" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7405" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7405" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7515" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7515" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7682" property="AmbientColor"/>
      <GlobalPropertyLink global_name="ForegroundColor" proxy="7682" property="CubeAxesColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="5360" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="5453" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="5546" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="5623" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="6424" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="6681" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="6783" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="6872" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="6961" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="7071" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="7238" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="7405" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="7515" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SelectionColor" proxy="7682" property="SelectionColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5360" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5360" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5453" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5453" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5546" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5546" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5623" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="5623" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6424" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6424" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6681" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6681" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6783" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6783" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6872" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6872" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6961" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="6961" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7071" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7071" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7238" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7238" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7405" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7405" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7515" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7515" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7682" property="BackfaceDiffuseColor"/>
      <GlobalPropertyLink global_name="SurfaceColor" proxy="7682" property="DiffuseColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="XLabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="XTitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="YLabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="YTitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="ZLabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5227" property="ZTitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="6135" property="OrientationAxesLabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5235" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5235" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5243" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5243" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5251" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5251" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5259" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5259" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5267" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5267" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5275" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5275" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5639" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5639" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5803" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5803" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5809" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5809" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5894" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="5894" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="6133" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="6133" property="TitleColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="6699" property="LabelColor"/>
      <GlobalPropertyLink global_name="TextAnnotationColor" proxy="6699" property="TitleColor"/>
      <GlobalPropertyLink global_name="UseGradientBackground" proxy="6135" property="UseGradientBackground"/>
    </Links>
  </ServerManagerState>
  </ParaView>""")

    end
    if verbose
        info("Output file: " * filename)
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
    run(`bash -c "rm -rf $(folder)"`)
end
