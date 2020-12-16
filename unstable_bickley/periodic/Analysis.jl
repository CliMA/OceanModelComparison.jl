module Analysis

using Oceananigans
using JLD2
using Plots

using Oceananigans.Grids
using Oceananigans.Fields: offset_data
import Oceananigans.Fields: Field

Field(loc::Tuple, grid::AbstractGrid, raw_data::Array) = Field(loc, CPU(), grid, nothing, offset_data(raw_data, grid, loc))
field_from_file(file, loc, name, iter) = Field(loc, file["serialized/grid"], file["timeseries/$name/$iter"])

CellFileField(file, name, iter)  = Field((Cell, Cell, Cell), file["serialized/grid"], file["timeseries/$name/$iter"])
XFaceFileField(file, name, iter) = Field((Face, Cell, Cell), file["serialized/grid"], file["timeseries/$name/$iter"])
YFaceFileField(file, name, iter) = Field((Cell, Face, Cell), file["serialized/grid"], file["timeseries/$name/$iter"])

max_abs_tracer(file, iter) = maximum(abs, interior(CellFileField(file, :c, iter)))

tracer_variance(file, iter) = 1/2 * sum(interior(CellFileField(file, :c, iter)).^2)

kinetic_energy(file, iter) = 1/2 * (sum(interior(XFaceFileField(file, :u, iter)).^2) + 
                                    sum(interior(YFaceFileField(file, :v, iter)).^2))

function oceananigans_statistics(filename)
    file = jldopen(filename)

    iterations = parse.(Int, keys(file["timeseries/t"]))

    t  = [ file["timeseries/t/$iter"]  for iter in iterations]
    e  = [ kinetic_energy(file, iter)  for iter in iterations]
    χ  = [ tracer_variance(file, iter) for iter in iterations]
    c★ = [ max_abs_tracer(file, iter)  for iter in iterations]

    close(file)

    Δχ = χ ./ χ[1]
    Δe = e ./ e[1]

    return t, Δe, Δχ, c★
end

function plot_resolutions(t, data, resolutions; kwargs...)
    plt = plot(t[1], data[1]; label = "N² = $(resolutions[1])²", kwargs...)
    for i = 2:length(data)
        plot!(plt, t[i], data[i]; label = "N² = $(resolutions[i])²", kwargs...)
    end
    return plt
end

end # module
