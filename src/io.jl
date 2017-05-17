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
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vts"`)
end
