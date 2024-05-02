# ---------
# kriging_interpolation.jl 
# author: Yudha Styawan
# email: yudhastyawan26@gmail.com 
# free to use
# ---------

# framework
using GeoStats
import ArchGDAL as AG
using Statistics

# IO modules
using CSV
using GeoIO
using Printf

# viz modules
import CairoMakie as Mke
using Makie

# EDIT the parameters here as the kriging output
lon = (99.0, 101.0)
lat = (-1.0, 1.0)
dimen = (1000, 1000)

# read data
# data must be delimited in tab
# long	lat	coulomb
# 99.000    	     -1.000   	    0.0011530
# 99.000    	     -0.950   	    0.0011910
# ...
dtable = georef(CSV.File("/Users/yudhastyawan/Documents/projects/julia/miscs/coulomb_outUSGS.dat"; skipto=2, header=[:x, :y, :value], delim="\t"), ("x", "y"))

# check value distribution
fig = Mke.Figure()
Mke.Axis(fig[1,1], xlabel="value", ylabel="occurence number")
Mke.hist!(fig[1,1], dtable.value, bins=15, color = :red, strokewidth = 1, strokecolor = :black)
fig
save("data_distribution.png", fig)

# generate the empirical variogram and fitted to the theoretical variogram
g = EmpiricalVariogram(dtable, "value", estimator = :cressie)

# possible variograms
γ_str = ["Gaussian", "Spherical", "Exponential", "Cubic"]
γ_func = [GaussianVariogram, SphericalVariogram, ExponentialVariogram, CubicVariogram]
γs1 = [GeoStatsFunctions.fit(v, g) for v in γ_func]
γs2 = [GeoStatsFunctions.fit(v, g; maxrange=2.0, maxsill=1.0) for v in γ_func]
γs = permutedims([γs1 γs2])

# plot empirical and theoretical variograms with the correlation coefficients
fig = Mke.Figure(size=(1000,500))
n = 1
for i in 1:2
    for j in 1:4
        γ = γs[i,j]
        theo = [γ(x) for x in g.abscissa]
        corcoeff = cor(theo, g.ordinate)
        ax = Mke.Axis(fig[i,j], title=@sprintf "cor=%.2f" corcoeff)
        if i == 2
            ax.xlabel = "lags"
        end
        Mke.lines!(fig[i,j], g.abscissa, theo.*10^4, color=:red, label="est")
        Mke.scatter!(fig[i,j], g.abscissa, g.ordinate.*10^4, label="obs")
        Mke.Label(fig[i, j, Mke.Top()], halign = :left, L"\times 10^{-4}")
        Mke.text!(0,3, text="$(γ_str[j])\n├─ sill: $(@sprintf "%.3E" sill(γ))\n├─ range: $(@sprintf "%.3E" range(γ))\n└─ nugget: $(@sprintf "%.3E" nugget(γ))", align=(:left, :top))
        axislegend(position=:rb)
        n = n+1
    end
end
fig
save("variograms.png", fig)

# selected variogram
γ = γs[2,4]

# kriging
model = GeoStatsModels.OrdinaryKriging(γ)
grid = CartesianGrid((lon[1],lat[1],), (lon[2],lat[2],), dims=dimen)
interp = dtable |> InterpolateNeighbors(grid, model, prob=true)

μ_value = [x.μ for x in interp.value]
μ_table = georef((value=μ_value,),grid)

σ_value = [x.σ for x in interp.value]
σ_table = georef((err_std=σ_value,),grid)

# check error STD distribution
fig = Mke.Figure()
Mke.Axis(fig[1,1], xlabel="error σ", ylabel="occurence number")
Mke.hist!(fig[1,1], σ_table.err_std, bins=15, color = :red, strokewidth = 1, strokecolor = :black)
fig
save("error_distribution.png", fig)

# show the distribution after interpolation
fig = Mke.Figure()
crange = (-0.3, 0.3)
cmap = "bwr"
Mke.Axis(fig[1,1], xlabel="longitude", ylabel="latitude")
viz!(fig[1,1], μ_table.geometry, color=μ_table.value, colormap=cmap, colorrange=crange)
cbar(fig[1,2], μ_table.value, colorrange=crange, colormap=cmap)
fig
save("kriging_output.png", fig)

# check observation and estimation values of the nearest points
geom_ctrs = centroid.(μ_table.geometry)
idcs = [argmin([√(sum((x.coords .- y.coords).^2)) for x in geom_ctrs]) for y in dtable.geometry]
μ_nearest = μ_table.value[idcs]
σ_nearest = σ_table.err_std[idcs]

# check observation and calculation correlation
fig = Mke.Figure()
Mke.Axis(fig[1,1], xlabel="closest-point estimated value", ylabel="observed value")
Mke.scatter!(fig[1,1], μ_nearest, dtable.value)
Mke.ablines!(fig[1,1], 0,1, color=:red, label=L"y = x")
axislegend(position=:rb)
fig
save("obs_vs_est.png", fig)

# check error STD distribution
fig = Mke.Figure()
Mke.Axis(fig[1,1], xlabel="error σ", ylabel="occurence number")
Mke.hist!(fig[1,1], σ_nearest, bins=15, color = :red, strokewidth = 1, strokecolor = :black)
fig
save("nearest_error_distribution.png", fig)

# save to tiff
# ---------------
function save_to_tiff(values, outnames)
    Z = reshape(values[:], dimen)

    # Resolution
    resx = (lon[2] - lon[1]) /dimen[1]
    resy = (lat[2] - lat[1]) /dimen[2]

    # Upper-left pixel coordinates
    ul_x = lon[1] - resx/2
    ul_y = lat[1] - resy/2

    # gt = [ul_x, resx, 0.0, ul_y, 0.0, resy]
    gt = [lon[1], resx, 0.0, lat[1], 0.0, resy]

    crs = AG.toWKT(AG.importPROJ4("+proj=latlong"))

    AG.create(
        outnames,
        driver = AG.getdriver("GTiff"),
        width=dimen[1],
        height=dimen[2],
        nbands=1,
        dtype=Float64
    ) do dataset
        AG.write!(dataset, Z, 1)

        AG.setgeotransform!(dataset, gt)
        AG.setproj!(dataset, crs)
    end
end

save_to_tiff(μ_table.value, "output_values.tiff")
save_to_tiff(σ_table.err_std, "output_errstd.tiff")
# ---------------