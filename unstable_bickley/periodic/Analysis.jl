module Analysis

using JLD2
using Plots

max_abs_tracer(file, iter) = maximum(abs, file["timeseries/c/$iter"])

tracer_variance(file, iter) = 1/2 * sum(file["timeseries/c/$iter"].^2)

kinetic_energy(file, iter) = 1/2 * (sum(file["timeseries/u/$iter"].^2) +
                                    sum(file["timeseries/v/$iter"].^2))

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
