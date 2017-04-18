import WriteVTK  # Install with Pkg.add("WriteVTK")

## IO functions

function writeVTK(simulation::Simulation;
                     folder::String=".",
                     verbose::Bool=false)

    simulation.file_number += 1
    filename = string(folder, "/", simulation.id, ".", simulation.file_number)

    # allocate arrays
    density = Array(Float64, length(simulation.ice_floes))

    thickness = Array(Float64, length(simulation.ice_floes))
    contact_radius = Array(Float64, length(simulation.ice_floes))
    areal_radius = Array(Float64, length(simulation.ice_floes))
    surface_area = Array(Float64, length(simulation.ice_floes))
    volume = Array(Float64, length(simulation.ice_floes))
    mass = Array(Float64, length(simulation.ice_floes))
    moment_of_inertia = Array(Float64, length(simulation.ice_floes))

    lin_pos = Array(Float64, 3, length(simulation.ice_floes))
    lin_vel = Array(Float64, 3, length(simulation.ice_floes))
    lin_acc = Array(Float64, 3, length(simulation.ice_floes))
    force = Array(Float64, 3, length(simulation.ice_floes))

    ang_pos = Array(Float64, length(simulation.ice_floes))
    ang_vel = Array(Float64, length(simulation.ice_floes))
    ang_acc = Array(Float64, length(simulation.ice_floes))
    torque = Array(Float64, length(simulation.ice_floes))

    ang_pos = Array(Float64, 3, length(simulation.ice_floes))
    ang_vel = Array(Float64, 3, length(simulation.ice_floes))
    ang_acc = Array(Float64, 3, length(simulation.ice_floes))
    torque = Array(Float64, 3, length(simulation.ice_floes))

    fixed = Array(Int, length(simulation.ice_floes))
    rotating = Array(Int, length(simulation.ice_floes))

    contact_stiffness_normal = Array(Float64, length(simulation.ice_floes))
    contact_stiffness_tangential = Array(Float64, length(simulation.ice_floes))
    contact_viscosity_normal = Array(Float64, length(simulation.ice_floes))
    contact_viscosity_tangential = Array(Float64, length(simulation.ice_floes))
    contact_static_friction = Array(Float64, length(simulation.ice_floes))
    contact_dynamic_friction = Array(Float64, length(simulation.ice_floes))

    # fill arrays
    for i=1:length(simulation.ice_floes)
        density[i] = simulation.ice_floes[i].density

        thickness[i] = simulation.ice_floes[i].thickness
        contact_radius[i] = simulation.ice_floes[i].contact_radius
        areal_radius[i] = simulation.ice_floes[i].areal_radius
        surface_area[i] = simulation.ice_floes[i].surface_area
        volume[i] = simulation.ice_floes[i].volume
        mass[i] = simulation.ice_floes[i].mass
        moment_of_inertia[i] = simulation.ice_floes[i].moment_of_inertia

        lin_pos[1:2, i] = simulation.ice_floes[i].lin_pos
        lin_vel[1:2, i] = simulation.ice_floes[i].lin_vel
        lin_acc[1:2, i] = simulation.ice_floes[i].lin_acc
        force[1:2, i] = simulation.ice_floes[i].force

        ang_pos[3, i] = simulation.ice_floes[i].ang_pos
        ang_vel[3, i] = simulation.ice_floes[i].ang_vel
        ang_acc[3, i] = simulation.ice_floes[i].ang_acc
        torque[3, i] = simulation.ice_floes[i].torque

        fixed[i] = Int(simulation.ice_floes[i].fixed)
        rotating[i] = Int(simulation.ice_floes[i].rotating)

        contact_stiffness_normal[i] = 
            simulation.ice_floes[i].contact_stiffness_normal
        contact_stiffness_tangential[i] = 
            simulation.ice_floes[i].contact_stiffness_tangential
        contact_viscosity_normal[i] = 
            simulation.ice_floes[i].contact_viscosity_normal
        contact_viscosity_tangential[i] = 
            simulation.ice_floes[i].contact_viscosity_tangential
        contact_static_friction[i] = 
            simulation.ice_floes[i].contact_static_friction
        contact_dynamic_friction[i] = 
            simulation.ice_floes[i].contact_dynamic_friction
    end
    
    # write to disk
    vtkfile = WriteVTK.vtk_grid(filename, lin_pos, WriteVTK.MeshCell[])

    WriteVTK.vtk_point_data(vtkfile, density, "Density [kg m^-3]")

    WriteVTK.vtk_point_data(vtkfile, thickness, "Thickness [m]")
    WriteVTK.vtk_point_data(vtkfile, contact_radius*2.,
                            "Diameter (contact) [m]")
    WriteVTK.vtk_point_data(vtkfile, areal_radius*2., "Diameter (areal) [m]")
    WriteVTK.vtk_point_data(vtkfile, surface_area, "Surface area [m^2]")
    WriteVTK.vtk_point_data(vtkfile, volume, "Volume [m^3]")
    WriteVTK.vtk_point_data(vtkfile, mass, "Mass [kg]")
    WriteVTK.vtk_point_data(vtkfile, moment_of_inertia,
                            "Moment of inertia [kg m^2]")

    WriteVTK.vtk_point_data(vtkfile, lin_vel, "Linear velocity [m s^-1]")
    WriteVTK.vtk_point_data(vtkfile, lin_acc, "Linear acceleration [m s^-2]")
    WriteVTK.vtk_point_data(vtkfile, force, "Sum of forces [N]")

    WriteVTK.vtk_point_data(vtkfile, ang_pos, "Angular position [rad]")
    WriteVTK.vtk_point_data(vtkfile, ang_vel, "Angular velocity [rad s^-1]")
    WriteVTK.vtk_point_data(vtkfile, ang_acc, "Angular acceleration [rad s^-2]")
    WriteVTK.vtk_point_data(vtkfile, torque, "Sum of torques [N*m]")

    WriteVTK.vtk_point_data(vtkfile, fixed, "Fixed in space [-]")
    WriteVTK.vtk_point_data(vtkfile, rotating, "Free to rotate [-]")

    WriteVTK.vtk_point_data(vtkfile, contact_stiffness_normal,
                            "Contact stiffness (normal) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, contact_stiffness_tangential,
                            "Contact stiffness (tangential) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, contact_viscosity_normal,
                            "Contact viscosity (normal) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, contact_viscosity_tangential,
                            "Contact viscosity (tangential) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, contact_static_friction,
                            "Contact friction (static) [-]")
    WriteVTK.vtk_point_data(vtkfile, contact_dynamic_friction,
                            "Contact friction (dynamic) [-]")

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        println("Output file: " * outfiles[1])
    else
        return nothing
    end
end
