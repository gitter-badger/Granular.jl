import WriteVTK
import NetCDF

## IO functions

export writeVTK
"""
Write a VTK file to disk containing all ice floes in the `simulation` in an 
unstructured mesh (file type `.vtu`).  These files can be read by ParaView and 
can be visualized by applying a *Glyph* filter.

If the simulation contains an `Ocean` data structure, it's contents will be 
written to separate `.vtu` files.  This can be disabled by setting the argument 
`ocean=false`.
"""
function writeVTK(simulation::Simulation;
                  folder::String=".",
                  verbose::Bool=true,
                  ocean::Bool=true)

    simulation.file_number += 1
    filename = string(folder, "/", simulation.id, ".icefloes.", 
                      simulation.file_number)
    writeIceFloeVTK(simulation, filename, verbose=verbose)

    filename = string(folder, "/", simulation.id, ".icefloe-interaction.", 
                      simulation.file_number)
    writeIceFloeInteractionVTK(simulation, filename, verbose=verbose)

    if typeof(simulation.ocean.input_file) != Bool && ocean
        filename = string(folder, "/", simulation.id, ".ocean.", 
                        simulation.file_number)
        writeOceanVTK(simulation.ocean, filename, verbose=verbose)
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
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmos_drag_coeff_vert,
                            "Atmosphere drag coefficient (vertical) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmos_drag_coeff_horiz,
                            "Atmosphere drag coefficient (horizontal) [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.pressure,
                            "Contact pressure [Pa]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.n_contacts,
                            "Number of contacts [-]")

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    else
        return nothing
    end
end

export writeIceFloeInteractionVTK
function writeIceFloeInteractionVTK(simulation::Simulation,
                                    filename::String;
                                    verbose::Bool=false)

    # Save ice-floe indexes and metrics for all interactions
    i1 = []
    i2 = []
    force = []
    shear_displacement_1 = []
    shear_displacement_2 = []
    contact_age = []
    for i=1:length(simulation.ice_floes)
        for ic=1:Nc_max
            if simulation.ice_floes[i].contacts[ic] > 0
                j = simulation.ice_floes[i].contacts[ic]

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

                
                append!(i1, i)
                append!(i2, j)

                append!(force, k_n*δ_n)

                append!(shear_displacement_1, simulation.ice_floes[i].
                        contact_parallel_displacement[ic][1])
                append!(shear_displacement_2, 
                        simulation.ice_floes[i].
                        contact_parallel_displacement[ic][2])

                append!(contact_age, simulation.ice_floes[i].contact_age[ic])
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
        write(f, "        <DataArray name=\"Force [N]\" type=\"Float32\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(force[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray name=\"Contact age [s]\" " * 
              "type=\"Float32\" NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(contact_age[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </CellData>\n")
        write(f, "      <Points>\n")

        # Write line endpoints (ice floe centers)
        #write(f, "        <DataArray name=\"Position [m]\" type=\"Float32\" " *
        write(f, "        <DataArray name=\"Points\" type=\"Float32\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i in simulation.ice_floes
            write(f, "$(i.lin_pos[1]) $(i.lin_pos[2]) 0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Points>\n")
        write(f, "      <Verts>\n")
        write(f, "        <DataArray name=\"connectivity\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray name=\"offsets\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Verts>\n")
        write(f, "      <Lines>\n")

        # Write contact connectivity by referring to point indexes
        write(f, "        <DataArray name=\"connectivity\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(i1[i] - 1) $(i2[i] - 1) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        
        # Write 0-indexed offset for the connectivity array for the end of each 
        # cell
        write(f, "        <DataArray name=\"offsets\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$((i - 1)*2 + 2) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </Lines>\n")
        write(f, "      <Strips>\n")
        write(f, "        <DataArray name=\"connectivity\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray name=\"offsets\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Strips>\n")
        write(f, "      <Polys>\n")
        write(f, "        <DataArray name=\"connectivity\" type=\"Int64\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray name=\"offsets\" type=\"Int64\" " *
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
function writeOceanVTK(ocean::Ocean,
                       filename::String;
                       verbose::Bool=false)
    
    # make each coordinate array three-dimensional
    xq = similar(ocean.u[:,:,:,1])
    yq = similar(ocean.u[:,:,:,1])
    zq = similar(ocean.u[:,:,:,1])

    for iz=1:size(xq, 3)
        xq[:,:,iz] = ocean.xq
        yq[:,:,iz] = ocean.yq
    end
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            zq[ix,iy,:] = ocean.zl
        end
    end

    # add arrays to VTK file
    vtkfile = WriteVTK.vtk_grid(filename, xq, yq, zq)

    WriteVTK.vtk_point_data(vtkfile, ocean.u[:, :, :, 1],
                            "u: Zonal velocity [m/s]")
    WriteVTK.vtk_point_data(vtkfile, ocean.v[:, :, :, 1],
                            "v: Meridional velocity [m/s]")
    # write velocities as 3d vector
    vel = zeros(3, size(xq, 1), size(xq, 2), size(xq, 3))
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            for iz=1:size(xq, 3)
                vel[1, ix, iy, iz] = ocean.u[ix, iy, iz, 1]
                vel[2, ix, iy, iz] = ocean.v[ix, iy, iz, 1]
            end
        end
    end
    
    WriteVTK.vtk_point_data(vtkfile, vel, "Velocity vector [m/s]")

    WriteVTK.vtk_point_data(vtkfile, ocean.h[:, :, :, 1],
                            "h: Layer thickness [m]")
    WriteVTK.vtk_point_data(vtkfile, ocean.e[:, :, :, 1],
                            "e: Relative interface height [m]")

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    else
        return nothing
    end
end

export removeSimulationFiles
"""
    removeSimulationFiles(simulation[, folder])

Remove all simulation output files from the specified folder.
"""
function removeSimulationFiles(simulation::Simulation; folder::String=".")
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtu"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtp"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vts"`)
end
