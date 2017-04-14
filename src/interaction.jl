## Interaction functions

"""
Resolve mechanical interaction between all grain pairs and walls.
"""
function interact()
    contact_pair = Array{Integer, 1}
    #overlap_ij = float

    # IceFloe to grain collisions
    while !isempty(g_contact_pairs)
        contact_pair = pop!(g_contact_pairs)
        interactIceFloes(contact_pair[1], contact_pair[2])
    end

    # IceFloe to wall collisions
    while !isempty(g_wall_contacts)
        contact_pair = pop!(g_wall_contacts)
        interactIceFloeWall(contact_pair[1], contact_pair[2])
    end

    #for k=1:length(g_contact_pairs)
    #end
end

"""
Resolve an grain-to-grain interaction using a prescibed contact law.
"""
function interactIceFloes(i::Integer, j::Integer,
    overlap_vector::vector;
    contact_normal::String = "LinearElastic")

    if contact_normal == "None"
        # do nothing

    elseif contact_normal == "LinearElastic"
        interactNormalLinearViscous(i, j, overlap_vector)

    else
        error("Unknown contact_normal interaction model '$contact_normal'")
    end

end

function interactNormalLinearElastic(i::Integer, j::Integer,
    overlap_vector::vector)

    force = -g_contact_stiffness_normal * overlap_vector

    g_force[i] = g_force[i]::vector + force;
    g_force[j] = g_force[j]::vector - force;
end
