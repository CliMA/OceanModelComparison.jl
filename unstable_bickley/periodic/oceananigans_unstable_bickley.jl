using JLD2

ENV["GKSwstype"] = "nul"
using Plots
using Printf
using Statistics

using Oceananigans
using Oceananigans.Advection
using Oceananigans.AbstractOperations
using Oceananigans.OutputWriters
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.Forcings

include("Bickley.jl")

using .Bickley

make_name(Nh, advection) = "oceananigans_unstable_bickley_jet_Nh$(Nh)_$(typeof(advection).name.wrapper)"

function run(; Nh = 128,
               output_time_interval = 2,
               stop_time = 200,
               arch = CPU(),
               ν = 0,
               advection = WENO5())

    name = make_name(Nh, advection)

    grid = RegularCartesianGrid(size=(Nh, Nh, 1),
                                x = (-2π, 2π), y=(-2π, 2π), z=(0, 1),
                                topology = (Periodic, Periodic, Bounded))

    model = IncompressibleModel(architecture = arch,
                                 timestepper = :RungeKutta3, 
                                   advection = advection,
                                        grid = grid,
                                     tracers = :c,
                                     closure = IsotropicDiffusivity(ν=ν, κ=ν),
                                    buoyancy = nothing)

    # ** Initial conditions **
    #
    # u, v: Large-scale jet + vortical perturbations
    #    c: Sinusoid
    
    # Parameters
    ϵ = 0.1 # perturbation magnitude
    ℓ = 0.5 # Gaussian width
    k = 0.5 # Sinusoidal wavenumber

    # Total initial conditions
    uᵢ(x, y, z) = Bickley.U(y) + ϵ * Bickley.ũ(x, y, ℓ, k)
    vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
    cᵢ(x, y, z) = Bickley.C(y, grid.Ly)

    set!(model, u=uᵢ, v=vᵢ, c=cᵢ)
    
    progress(sim) = @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|u|: %.2f",
                                   sim.model.clock.iteration, sim.model.clock.time,
                                   sim.Δt.Δt, maximum(abs, u.data.parent)))

    wizard = TimeStepWizard(cfl=0.5, Δt=1e-4, max_change=1.1, max_Δt=10.0)

    simulation = Simulation(model, Δt=wizard, stop_time=stop_time,
                            iteration_interval=10, progress=progress)

    # Output: primitive fields + computations
    u, v, w, c = primitives = merge(model.velocities, model.tracers)

    ω   = ComputedField( ∂x(v) - ∂y(u)      )
    ω²  = ComputedField( (∂x(v) - ∂y(u))^2  )
    c²  = ComputedField( c^2                )

    ∇c² = @at (Cell, Cell, Cell) ∂x(c)^2 + ∂y(c)^2
    ∇c² = ComputedField(∇c²)

    computations = (ω=ω, ω²=ω², c²=c², ∇c²=∇c²)
    outputs = merge(primitives, computations)

    save_grid = (file, model) -> file["serialized/grid"] = model.grid

    #JLD2OutputWriter(model, merge(model.velocities, model.tracers, (ω=ω, ∇c²=∇c²)),
    #
    simulation.output_writers[:fields] =
        JLD2OutputWriter(model, merge(model.velocities, model.tracers, (ω=ω,)),
                                schedule = TimeInterval(output_time_interval),
                                init = save_grid,
                                prefix = name * "_fields",
                                force = true)

    #=
    averages = Dict(name => mean(ϕ, dims=(1, 2, 3))
                    for (name, ϕ) in zip(keys(computations), values(computations)))

    simulation.output_writers[:averages] =
        JLD2OutputWriter(model, averages,
                         schedule = TimeInterval(output_time_interval),
                         init = save_grid,
                         prefix = name * "_statistics",
                         force = true)
    =#

    @info "Running a simulation of an unstable Bickley jet with $(Nh)² degrees of freedom..."

    start_time = time_ns()

    run!(simulation)

    run_time = (time_ns() - start_time) * 1e-9

    @info "The simulation with Nh = $Nh ran for $run_time seconds"
    @show cost = run_time / (model.clock.iteration * Nh^2)

    return name
end

function analyze(name)
    @info "Analyzing the results of an unstable Bickley jet simulation..."

    statistics_file = jldopen(name * "_statistics.jld2")

    iterations = parse.(Int, keys(statistics_file["timeseries/t"]))
    grid = statistics_file["serialized/grid"]

    ∇c² = [statistics_file["timeseries/∇c²/$iter"][1, 1, 1] for iter in iterations]
    c²  = [statistics_file["timeseries/c²/$iter"][1, 1, 1] for iter in iterations]
    t   = [statistics_file["timeseries/t/$iter"][1, 1, 1] for iter in iterations]

    close(statistics_file)

    plot(t, 2π / grid.Ly .* ∇c² ./ c²)
    savefig(name * ".png")

    return nothing
end

function visualize(name, contours=false)
    @info "Making a fun movie about an unstable Bickley jet..."

    fields_file = jldopen(name * "_fields.jld2")

    iterations = parse.(Int, keys(fields_file["timeseries/t"]))
    grid = fields_file["serialized/grid"]

    xu, yu, zu = nodes((Face, Cell, Cell), grid)
    xω, yω, zω = nodes((Face, Face, Cell), grid)
    xc, yc, zc = nodes((Cell, Cell, Cell), grid)

    anim = @animate for (i, iteration) in enumerate(iterations)

        @info "    Plotting frame $i from iteration $iteration..."
        
        t = fields_file["timeseries/t/$iteration"]
        ω = fields_file["timeseries/ω/$iteration"][:, :, 1]
        u = fields_file["timeseries/u/$iteration"][:, :, 1]
        c = fields_file["timeseries/c/$iteration"][:, :, 1]

        kwargs = Dict(:xlabel => "x",
                      :ylabel => "y",
                      :aspectratio => 1,
                      :linewidth => 0,
                      :colorbar => true,
                      :clims => (-1, 1),
                      :xlims => (-grid.Lx/2, grid.Lx/2),
                      :ylims => (-grid.Ly/2, grid.Ly/2))

        contours && (kwargs[:levels] = range(-1, 1, length=31))
        plotter = contours ? contourf : heatmap

        u_plot = plotter(xu, yu, clamp.(u, -1, 1)'; color = :balance, kwargs...)
        ω_plot = plotter(xω, yω, clamp.(ω, -1, 1)'; color = :balance, kwargs...)
        c_plot = plotter(xc, yc, clamp.(c, -1, 1)'; color = :thermal, kwargs...)

        u_title = @sprintf("u at t = %.1f", t)
        ω_title = @sprintf("ω at t = %.1f", t)
        c_title = @sprintf("c at t = %.1f", t)

        plot(u_plot, ω_plot, c_plot,
             title = [u_title ω_title c_title],
             layout = (1, 3),
             size = (1600, 400))
    end

    gif(anim, name * ".gif", fps = 8)

    return nothing
end

for Nh in (32, 64, 128, 256, 512, 1024, 2048)
    #name = run(Nh=Nh, arch=GPU(), advection=WENO5())
    name = make_name(Nh, WENO5())
    visualize(name)
end

#=
for Nh in (32, 64, 128, 256, 512, 1024)
    name = run(Nh=Nh, arch=CPU(), advection=UpwindBiasedFifthOrder())
end
=#
