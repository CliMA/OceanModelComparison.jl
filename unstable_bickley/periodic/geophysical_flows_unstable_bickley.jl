# Unstable Bickley jet in GeophysicalFlows

ENV["GKSwstype"] = "nul"
using Printf
using Random
using Plots

using FFTW: rfft, irfft
using CUDA
using JLD2

using FourierFlows
using GeophysicalFlows: TwoDNavierStokes

include("Bickley.jl")

using .Bickley

function run(; dev = CPU(), nx = 128, dt = 1e-2)

    nsubs = ceil(Int, 2 / dt) # number of steps between each plot
    nsteps = 100 * nsubs # total number of steps

    problem = TwoDNavierStokes.Problem(dev; nx = nx,
                                            Lx = 4π,
                                            ny = nx,
                                            Ly = 4π,
                                            dt = dt,
                                            stepper = "FilteredRK4")

    ϵ = 0.1
    k = 0.5
    ℓ = 0.5

    DeviceArray = dev isa CPU ? Array : CuArray

    grid = problem.grid
    x = repeat(reshape(grid.x, grid.nx, 1), 1, grid.ny)
    y = repeat(reshape(grid.y, 1, grid.ny), grid.nx, 1)

    ψᵢ = DeviceArray(@. Bickley.Ψ(y, grid.Ly) + ϵ * Bickley.ψ̃(x, y, k, ℓ))

    ψᵢh = rfft(ψᵢ)
    CUDA.@allowscalar ψᵢh[1, 1] = 0

    ψᵢ = irfft(ψᵢh, grid.nx)

    @. problem.vars.zetah = - grid.Krsq * ψᵢh

    ζᵢ = irfft(problem.vars.zetah, problem.grid.nx)

    TwoDNavierStokes.set_zeta!(problem, ζᵢ)
    TwoDNavierStokes.set_c!(problem, DeviceArray(sin.(y/2)))

    # ## Output

    # We choose folder for outputing `.jld2` files and snapshots (`.png` files).
    filename = "geophysical_flows_unstable_bickley_jet_Nh$nx.jld2"

    isfile(filename) && rm(filename, force=true)

    output = Output(problem, filename,
        (:ζh, p -> Array(p.vars.zetah)),
        (:ζ,  p -> Array(irfft(p.vars.zetah, p.grid.nx))),
        (:ch, p -> Array(p.vars.ch)),
        (:c,  p -> Array(irfft(p.vars.ch, p.grid.nx))),
    )

    saveproblem(output)

    start_time = time_ns()

    e₀ = TwoDNavierStokes.energy(problem)
    χ₀ = TwoDNavierStokes.tracer_variance(problem)

    for j = 0:Int(nsteps/nsubs)

        stepforward!(problem, nsubs)
        TwoDNavierStokes.updatevars!(problem)
        saveoutput(output)

        e = TwoDNavierStokes.energy(problem)
        χ = TwoDNavierStokes.tracer_variance(problem)

        log = @sprintf("Step: %04d, t: %d, Δe/e₀: %.2e, Δχ/χ₀: %.2e, wall_time: %.2f seconds",
                       problem.clock.step,
                       problem.clock.t,
                       e / e₀ - 1,
                       χ / χ₀ - 1,
                       (time_ns() - start_time) * 1e-9,
                      )
        println(log)

    end

    run_time = (time_ns() - start_time) * 1e-9

    return output.path, run_time, nsteps
end

function visualize(filename)

    file = jldopen(filename)

    iterations = parse.(Int, keys(file["snapshots/t"]))

    x = file["grid/x"]
    y = file["grid/y"]

    animation = @animate for (i, iter) in enumerate(iterations)

        @info "Plotting frame $i of $(length(iterations))"

        t = file["snapshots/t/$iter"]
        ζ = file["snapshots/ζ/$iter"]
        c = file["snapshots/c/$iter"]

        kwargs = Dict(
                      :aspectratio => 1,
                      :linewidth => 0,
                      :ticks => nothing,
                      :colorbar => :none,
                      :clims => (-1, 1),
                      :xlims => (-2π, 2π),
                      :ylims => (-2π, 2π),
                     )

        ζ_plot = heatmap(x, y, clamp.(ζ, -1, 1)'; color=:balance, kwargs...)
        c_plot = heatmap(x, y, clamp.(c, -1, 1)'; color=:thermal, kwargs...)

        ω_title = @sprintf("ζ at t = %.1f", t)
        c_title = @sprintf("c at t = %.1f", t)

        plot(ζ_plot, c_plot, title = [ω_title c_title], size=(4000, 2000))
    end

    close(file)

    name = filename[1:end-5]

    gif(animation, name * ".gif", fps=8)

    return nothing
end

cfl = 0.1
for nx in (32, 
           64, 128, 256, 512,
           1024, 2048)
           #32, 64, 128, 256, 512)

    @show dt = cfl * 4π / nx
    filename, run_time, nsteps = run(dev=GPU(), nx=nx, dt=dt)

    DOF = nx^2
    @info "Results for nx = $nx !"
    @show cost = run_time / (nsteps * DOF)

    filename = "geophysical_flows_unstable_bickley_jet_Nh$nx.jld2"
    visualize(filename)
end
