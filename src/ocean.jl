"Returns empty ocean type for initialization purposes."
function createEmptyOcean()
    return Ocean("",
                 zeros(1),
                 zeros(1),
                 zeros(1),
                 zeros(1),
                 zeros(1),
                 zeros(1),
                 zeros(1),
                 zeros(1,1,1,1),
                 zeros(1,1,1,1),
                 zeros(1,1,1,1),
                 zeros(1,1,1,1))
end
