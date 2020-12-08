using JLD2
using Plots
using Printf
using Oceananigans
using Oceananigans.Grids

filename = "oceananigans_unstable_bickley_jet_Nh512_WENO5_fields.jld2"

file = jldopen(filename)

iterations = parse.(Int, keys(file["timeseries/t"]))
grid = file["serialized/grid"]

xc, yc, zc = nodes((Cell, Cell, Cell), grid)
xω, yω, zω = nodes((Face, Face, Cell), grid)

snapshots = []

for iteration in iterations[1:6:54]

    t = file["timeseries/t/$iteration"]
    ω = file["timeseries/ω/$iteration"][:, :, 1]
    c = file["timeseries/c/$iteration"][:, :, 1]

    kwargs = Dict(
                  :aspectratio => 1,
                  :linewidth => 0,
                  :colorbar => false,
                  :xticks => false,
                  :yticks => false,
                  :clims => (-1, 1),
                  :xlims => (-grid.Lx/2, grid.Lx/2),
                  :ylims => (-grid.Ly/2, grid.Ly/2),
                  :title => @sprintf("t = %.1f", t)
                 )

    #snap = heatmap(xc, yc, clamp.(c, -1, 1)'; color = :thermal, kwargs...)
    snap = heatmap(xω, yω, clamp.(ω, -1, 1)'; color = :balance, kwargs...)

    push!(snapshots, snap)
end

close(file)

plt = plot(snapshots..., size = (1800, 1600))

display(plt)

savefig("vorticity.png")
#savefig("dye.png")
