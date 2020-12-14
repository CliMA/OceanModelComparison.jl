using Plots

include("Analysis.jl")

using .Analysis

filenames = [
             "oceananigans_unstable_bickley_jet_Nh32_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh64_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh128_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh256_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh512_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh1024_WENO5_fields.jld2",
             "oceananigans_unstable_bickley_jet_Nh2048_WENO5_fields.jld2"
            ]

resolutions = [32, 64, 128, 256, 512, 1024, 2048]

results = [Analysis.oceananigans_statistics(filename) for filename in filenames]

 t = map(r -> r[1], results)
Δe = map(r -> r[2], results)
Δχ = map(r -> r[3], results)
c★ = map(r -> r[4], results)

kwargs = Dict(
              :xlabel => "Time",
             )

energy_plot   = Analysis.plot_resolutions(t, Δe, resolutions;
                                          ylabel = "Kinetic energy",
                                          legend = :bottomleft,
                                          kwargs...)

max_c_plot    = Analysis.plot_resolutions(t, c★, resolutions;
                                          ylabel = "\$ \\max |c| \$",
                                          legend = :bottomleft,
                                          kwargs...)

variance_plot = Analysis.plot_resolutions(t, Δχ, resolutions;
                                          ylabel = "Tracer variance", 
                                          legend = :bottomleft,
                                          kwargs...)

display(plot(energy_plot, variance_plot, max_c_plot, layout=(3, 1), size=(800, 600)))
