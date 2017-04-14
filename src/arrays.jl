## Global arrays for general simulation data
# TODO: Variables should be local, or passed as arguments to functions whenever
# possible
g_position = vector[]
g_velocity = vector[]
g_acceleration = vector[]
g_force = vector[]

g_rotational_position = vector[]
g_rotational_velocity = vector[]
g_rotational_acceleration = vector[]
g_torque = vector[]

g_radius = float[]

g_density = float[]
g_volume = float[]
g_mass = float[]
g_rotational_inertia = float[]

g_time_iteration = 0
g_time = 0.0
g_time_total = 0.0
g_time_step = 0.0
g_file_time_step = 0.0   # <= 0.0: no output files
g_file_number = 0

g_contact_stiffness_normal = 1.16e9
g_contact_stiffness_tangential = 0.0
g_contact_viscosity_normal = 1.16e9
g_contact_viscosity_tangential = 0.0

g_gravitational_acceleration = zeros(3)

g_contact_pairs = Array{Integer, 1}[]
#g_positions = float[]
g_overlaps = float[]
g_wall_contacts = Array{Integer, 1}[]

g_simulation_id = "unnamed"

g_origo = zeros(3)
g_world_size = zeros(3)
