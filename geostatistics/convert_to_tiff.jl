# ---------
# convert_to_tiff.jl 
# author: Yudha Styawan
# email: yudhastyawan26@gmail.com 
# free to use
# ---------

# framework
using DataFrames
import ArchGDAL as AG

# IO modules
using CSV

# read data
# data must be delimited in tab
# long	lat	coulomb
# 99.000    	     -1.000   	    0.0011530
# 99.000    	     -0.950   	    0.0011910
# ...
df = DataFrame(CSV.File("path/to/coulomb.dat"; skipto=2, header=[:x, :y, :value], delim="\t"))

# get the length
lenx = length(df[df.x .== df.x[1],:x])
leny = length(df[df.y .== df.y[1],:y])

# reshape to matrix and transpose
ZZ = Matrix(transpose(reshape(df.value, (lenx, leny))))

# Resolution
resx = (df.x[end] - df.x[1]) /lenx
resy = (df.y[end] - df.y[1]) /leny

# Upper-left pixel coordinates
# ul_x = df.x[1] - resx/2
# ul_y = df.y[1] - resy/2

# coordinate information
gt = [df.x[1], resx, 0.0, df.y[1], 0.0, resy]
# gt = [ul_x, resx, 0.0, ul_y, 0.0, resy]

# CRS
crs = AG.toWKT(AG.importPROJ4("+proj=latlong"))

# Save to raster (.tiff)
AG.create(
    "./coulomb.tiff",
    driver = AG.getdriver("GTiff"),
    width=lenx,
    height=leny,
    nbands=1,
    dtype=Float64
) do dataset
    AG.write!(dataset, ZZ, 1)

    AG.setgeotransform!(dataset, gt)
    AG.setproj!(dataset, crs)
end