## Global arrays for general simulation data
g_ice_floes = IceFloeCylindrical[]

g_time_iteration = 0
g_time = 0.0
g_time_total = 0.0
g_time_step = 0.0
g_file_time_step = 0.0   # <= 0.0: no output files
g_file_number = 0

g_gravitational_acceleration = zeros(3)

g_contact_pairs = Array{Integer, 1}[]
#g_positions = float[]
g_overlaps = float[]
g_wall_contacts = Array{Integer, 1}[]

g_simulation_id = "unnamed"

g_origo = zeros(3)
g_world_size = zeros(3)
