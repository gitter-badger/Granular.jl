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

"""
Read NetCDF file generated by MOM6 (e.g. `prog__####_###.nc`) from disk and 
return as `Ocean` data structure.
"""
function readOceanNetCDF(filename::String)

    if !isfile(filename)
        error("$(filename) could not be opened")
    end
    u_staggered::Array{float, 4} = NetCDF.ncread(filename, "u")
    v_staggered::Array{float, 4} = NetCDF.ncread(filename, "v")
    u, v = convertToColocatedOceanGrid(u_staggered, v_staggered)

    ocean = Ocean(filename,
                  NetCDF.ncread(filename, "Time"),

                  NetCDF.ncread(filename, "xq"),
                  NetCDF.ncread(filename, "yq"),
                  NetCDF.ncread(filename, "xh"),
                  NetCDF.ncread(filename, "yh"),
                  NetCDF.ncread(filename, "zl"),
                  NetCDF.ncread(filename, "zi"),

                  u,
                  v,
                  NetCDF.ncread(filename, "h"),
                  NetCDF.ncread(filename, "e"))
    return ocean
end

"""
Convert gridded data from staggered (Arakawa-C) to collocated grid (Arakawa-A) 
through interpolation.  The new data points are located in the centers of the 
original staggered grid (spatial coordinates `xh` and `yh`).
"""
function convertToColocatedOceanGrid(u_in::Array{float, 4},
                                     v_in::Array{float, 4})
    u = Array{float}(size(u_in))
    v = Array{float}(size(v_in))
    nx = size(u_in)[1]
    ny = size(u_in)[2]
    for i=1:nx
        for j=1:ny
            if j < ny - 1
                u[i, j, :, :] = (u_in[i, j, :, :] + u_in[i, j+1, :, :])/2.
            else
                u[i, j, :, :] = u_in[i, j, :, :]
            end
            if i < nx - 1
                v[i, j, :, :] = (v_in[i, j, :, :] + v_in[i+1, j, :, :])/2.
            else
                v[i, j, :, :] = v_in[i, j, :, :]
            end
        end
    end
    return u, v
end
