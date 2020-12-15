# Unstable Bickley jet

using Printf

ENV["GKSwstype"] = "nul" 
using Plots
using Revise
using CUDA
using ClimateMachine

ClimateMachine.init()

using ClimateMachine.Ocean
using ClimateMachine.Ocean.Domains
using ClimateMachine.Ocean.Fields

using ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
using ClimateMachine.GenericCallbacks: EveryXSimulationTime
using ClimateMachine.GenericCallbacks: EveryXSimulationSteps
using ClimateMachine.Ocean: current_step, Δt, current_time
using ClimateMachine.Ocean: JLD2Writer, OutputTimeSeries, write!
using CLIMAParameters: AbstractEarthParameterSet, Planet

struct NonDimensionalParameters <: AbstractEarthParameterSet end
Planet.grav(::NonDimensionalParameters) = 10
c = sqrt(Planet.grav(NonDimensionalParameters())) # gravity wave speed for unit depth

include("Bickley.jl")
using .Bickley

# Low-p assumption:
effective_node_spacing(Ne, Np, Lx=4π) = Lx / (Ne * (Np + 1))

function run(;
             Ne = 4,
             Np = 4,
             ν = 0,
             time_step = 0.1 * effective_node_spacing(Ne, Np) / c,
             array_type = Array,
             output_time_interval = 2,
             stabilizing_dissipation = nothing,
             stop_time = 200)

    name = @sprintf("climate_machine_unstable_bickley_jet_Ne%d_Np%d_ν%.1e_no_rotation", Ne, Np, ν)

    ClimateMachine.Settings.array_type = array_type

    # Domain

    domain = RectangularDomain(Ne = (Ne, Ne, 1), Np = Np,
                               x = (-2π, 2π), y = (-2π, 2π), z = (0, 1),
                               periodicity = (true, true, false))

    # Physical parameters:
    g = Planet.grav(NonDimensionalParameters())

    # Non-dimensional parameters
    ϵ = 0.1 # Perturbation amplitude
    ℓ = 0.5 # Perturbation width
    k = 0.5 # Perturbation wavenumber

    # Initial conditions: Jet/tracer + perturbations
    uᵢ(x, y, z) = Bickley.U(y, domain.L.y) + ϵ * Bickley.ũ(x, y, ℓ, k)
    vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
    θᵢ(x, y, z) = Bickley.C(y, domain.L.y)

    initial_conditions = InitialConditions(u=uᵢ, v=vᵢ, θ=θᵢ)

    model = Ocean.HydrostaticBoussinesqSuperModel(
        domain = domain,
        time_step = time_step,
        initial_conditions = initial_conditions,
        parameters = NonDimensionalParameters(),
        turbulence_closure = (νʰ = ν, κʰ = ν, νᶻ = ν, κᶻ = ν),
        rusanov_wave_speeds = (cʰ = sqrt(g * domain.L.z), cᶻ = 1e-2),
        stabilizing_dissipation = stabilizing_dissipation,
        coriolis = (f₀ = 0, β = 0),
        buoyancy = (αᵀ = 0,),
        boundary_tags = ((0, 0), (1, 1), (1, 2)),
        boundary_conditions = (OceanBC(Impenetrable(FreeSlip()), Insulating()),
                               OceanBC(Penetrable(FreeSlip()), Insulating()))
    )

    name = @sprintf("climate_machine_unstable_bickley_jet_Ne%d_Np%d_ν%.1e_no_rotation", Ne, Np, ν)

    ClimateMachine.Settings.array_type = array_type

    # Domain

    domain = RectangularDomain(Ne = (Ne, Ne, 1), Np = Np,
                               x = (-2π, 2π), y = (-2π, 2π), z = (0, 1),
                               periodicity = (true, true, false))

    # Physical parameters:
    g = Planet.grav(NonDimensionalParameters())

    # Non-dimensional parameters
    ϵ = 0.1 # Perturbation amplitude
    ℓ = 0.5 # Perturbation width
    k = 0.5 # Perturbation wavenumber

    # Initial conditions: Jet/tracer + perturbations
    uᵢ(x, y, z) = Bickley.U(y, domain.L.y) + ϵ * Bickley.ũ(x, y, ℓ, k)
    vᵢ(x, y, z) = ϵ * Bickley.ṽ(x, y, ℓ, k)
    θᵢ(x, y, z) = Bickley.C(y, domain.L.y)

    initial_conditions = InitialConditions(u=uᵢ, v=vᵢ, θ=θᵢ)

    model = Ocean.HydrostaticBoussinesqSuperModel(
        domain = domain,
        time_step = time_step,
        initial_conditions = initial_conditions,
        parameters = NonDimensionalParameters(),
        turbulence_closure = (νʰ = ν, κʰ = ν, νᶻ = ν, κᶻ = ν),
        rusanov_wave_speeds = (cʰ = sqrt(g * domain.L.z), cᶻ = 1e-2),
        stabilizing_dissipation = stabilizing_dissipation,
        coriolis = (f₀ = 0, β = 0),
        buoyancy = (αᵀ = 0,),
        boundary_tags = ((0, 0), (1, 1), (1, 2)),
        boundary_conditions = (OceanBC(Impenetrable(FreeSlip()), Insulating()),
                               OceanBC(Penetrable(FreeSlip()), Insulating()))
    )    # We prepare a callback that periodically fetches the horizontal velocity and
    # tracer concentration for later animation,

    writer = JLD2Writer(model, filepath = name * ".jld2", overwrite_existing = true)
    cpu_grid = DiscontinuousSpectralElementGrid(domain, array_type=Array)

    start_time = time_ns()

    data_fetcher = EveryXSimulationTime(output_time_interval) do
        write!(writer)

        cpu_data = convert(Array, model.state.realdata)
        u = SpectralElementField(domain, cpu_grid, view(cpu_data, :, 1, :))

        # Print a helpful message
        step = @sprintf("Step: %d", current_step(model))
        time = @sprintf("time: %.2f", current_time(model))
        max_u = @sprintf("max|u|: %.6f", maximum(abs, u))

        elapsed = (time_ns() - start_time) * 1e-9
        wall_time = @sprintf("elapsed wall time: %.2f min", elapsed / 60)  

        isnan(maximum(abs, u)) && error("NaNs.") 

        @info "$step, $time, $max_u, $wall_time"
    end

    # and then run the simulation.

    model.solver_configuration.timeend = stop_time

    total_steps = ceil(Int, stop_time / time_step)
    @info @sprintf("Running a simulation of the instability of the Bickley jet (Δt=%.2e, steps=%d)", time_step, total_steps)

    try
        result = ClimateMachine.invoke!(model.solver_configuration;
                                        user_callbacks = [data_fetcher])
    catch err
        @warn "Simulation ended prematurely because $(sprint(showerror, err))"
    end

    return name
end

function visualize(name, contours=false)

    filepath = name * ".jld2"

    u_timeseries = OutputTimeSeries(:u, filepath)
    v_timeseries = OutputTimeSeries(:v, filepath)
    η_timeseries = OutputTimeSeries(:η, filepath)
    c_timeseries = OutputTimeSeries(:θ, filepath)

    u₀ = u_timeseries[1]
    domain = u₀.domain
    assembled_u₀ = assemble(u₀)
    x = assembled_u₀.x[:, 1, 1]
    y = assembled_u₀.y[1, :, 1]

    times = u_timeseries.times

    animation = @animate for i = 1:length(u_timeseries)

        @info "Plotting frame $i of $(length(u_timeseries))..."

        kwargs = (xlim = domain.x, ylim = domain.y, linewidth = 0, aspectratio = 1)

        u = assemble(u_timeseries[i]).data[:, :, 1]
        v = assemble(v_timeseries[i]).data[:, :, 1]
        η = assemble(η_timeseries[i]).data[:, :, 1]
        c = assemble(c_timeseries[i]).data[:, :, 1]

        if ~isnan(maximum(abs, u))

            s = @. sqrt(u^2 + v^2)

            ηmin = minimum(η)
            ηmax = maximum(η)
            ηlevels = range(ηmin, ηmax, length=31)

            kwargs = Dict(:xlabel => "x",
                          :ylabel => "y",
                          :aspectratio => 1,
                          :linewidth => 0,
                          :colorbar => true,
                          :xlims => (-domain.L.x/2, domain.L.x/2),
                          :ylims => (-domain.L.y/2, domain.L.y/2))

            plotter = contours ? contourf : heatmap

            kwargs[:clims] = (0, 1)
            contours && (kwargs[:levels] = range(0, 1, length=31))

            s_plot = plotter(x, y, clamp.(s, 0, 1)'; color = :thermal, kwargs...)

            kwargs[:clims] = (-1, 1)
            contours && (kwargs[:levels] = range(-1, 1, length=31))

            u_plot = plotter(x, y, clamp.(u, -1, 1)'; color = :balance, kwargs...)
            c_plot = plotter(x, y, clamp.(c, -1, 1)'; color = :thermal, kwargs...)
            #η_plot = plotter(x, y, clamp.(c, -1, 1)'; color = :thermal, kwargs...)

            u_title = @sprintf("u at t = %.2f", times[i])
            s_title = @sprintf("speed at t = %.2f", times[i])
            c_title = @sprintf("c at t = %.2f", times[i])

            plot(u_plot, s_plot, c_plot,
                 title = [u_title s_title c_title],
                 layout = (1, 3),
                 size = (1600, 400))

        end
    end

    gif(animation, name * ".gif", fps = 8)

    return nothing
end

include("StabilizingDissipations.jl")
using .StabilizingDissipations: StabilizingDissipation

Ne = 8
Np = 3

time_step = 0.1 * effective_node_spacing(Ne, Np) / c

test_dissipation = StabilizingDissipation(minimum_node_spacing = effective_node_spacing(Ne, Np),
                                          time_step = time_step,
                                          Δu = 1e-3,
                                          Δθ = 1e-3)

name = run(Ne=Ne, Np=Np, stabilizing_dissipation=test_dissipation)
visualize(name)

#=
for DOF in (32, 64, 128, 256, 512)
    for Np in (2, 3, 4, 5, 6)
        Ne = round(Int, DOF / (Np+1))
        name = run(Ne=16, Np=3, safety=0.1)
        visualize(name)
    end
end

for DOF in (512,)
    for Np in (2, 3, 4, 5, 6)
        Ne = round(Int, DOF / (Np+1))
        name = run(Ne=16, Np=3, safety=0.1, ν=1e-4)
        visualize(name)
    end
end

for DOF in (1024,)
    for Np in (2, 3, 4, 5, 6)
        Ne = round(Int, DOF / (Np+1))
        name = run(Ne=16, Np=3, safety=0.1, ν=1e-5)
        visualize(name)
    end
end
=#
