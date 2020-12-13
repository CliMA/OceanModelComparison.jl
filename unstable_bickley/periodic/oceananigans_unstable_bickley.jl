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
    uᵢ(x, y, z) = Bickley.U(y, grid.Ly) + ϵ * Bickley.ũ(x, y, ℓ, k)
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

    ζ   = ComputedField(∂x(v) - ∂y(u))

    save_grid = (file, model) -> file["serialized/grid"] = model.grid

    #=
    simulation.output_writers[:fields] =
        JLD2OutputWriter(model, merge(model.velocities, model.tracers, (ζ=ζ,)),
                                schedule = TimeInterval(output_time_interval),
                                init = save_grid,
                                prefix = name * "_fields",
                                field_slicer = nothing,
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

using Oceananigans.Grids: AbstractGrid
using Oceananigans.Fields: offset_data
import Oceananigans.Fields: Field

Field(loc::Tuple, grid::AbstractGrid, raw_data::Array) = Field(loc, CPU(), grid, nothing, offset_data(raw_data, grid, loc))
field_from_file(file, loc, name, iter) = Field(loc, file["serialized/grid"], file["timeseries/$name/$iter"])

function visualize(name, contours=false)
    @info "Making a fun movie about an unstable Bickley jet..."

    file = jldopen(name * "_fields.jld2")

    iterations = parse.(Int, keys(file["timeseries/t"]))
    grid = file["serialized/grid"]

    xu, yu, zu = nodes((Face, Cell, Cell), grid)
    xζ, yζ, zζ = nodes((Face, Face, Cell), grid)
    xc, yc, zc = nodes((Cell, Cell, Cell), grid)

    anim = @animate for (i, iteration) in enumerate(iterations)

        @info "    Plotting frame $i from iteration $iteration..."
        
        t = file["timeseries/t/$iteration"]

        ζ = field_from_file(file, (Face, Face, Cell), :ζ, iteration)
        c = field_from_file(file, (Cell, Cell, Cell), :c, iteration)

        ζi = interior(ζ)[:, :, 1]
        ci = interior(c)[:, :, 1]

        kwargs = Dict(
                      :aspectratio => 1,
                      :linewidth => 0,
                      :colorbar => :none,
                      :ticks => nothing,
                      :clims => (-1, 1),
                      :xlims => (-grid.Lx/2, grid.Lx/2),
                      :ylims => (-grid.Ly/2, grid.Ly/2)
                     )

        contours && (kwargs[:levels] = range(-1, 1, length=31))
        plotter = contours ? contourf : heatmap

        ζ_plot = plotter(xζ, yζ, clamp.(ζi, -1, 1)'; color = :balance, kwargs...)
        c_plot = plotter(xc, yc, clamp.(ci, -1, 1)'; color = :thermal, kwargs...)

        ζ_title = @sprintf("ζ at t = %.1f", t)
        c_title = @sprintf("c at t = %.1f", t)

        plot(ζ_plot, c_plot, title = [ζ_title c_title], size = (4000, 2000))
    end

    gif(anim, name * ".gif", fps = 8)

    return nothing
end

for Nh in (32, 32, 64, 128, 256, 512)
    name = run(Nh=Nh, arch=CPU(), advection=WENO5())
    name = run(Nh=Nh, arch=CPU(), advection=UpwindBiasedFifthOrder())
    #visualize(name)
end
