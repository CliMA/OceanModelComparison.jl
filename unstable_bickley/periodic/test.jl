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

Ne = 4
Np = 7
ν = 0
time_step = 0.1 * effective_node_spacing(Ne, Np) / c
array_type = Array
output_time_interval = 2
stabilizing_dissipation = nothing
stop_time = 200

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
    state_filter_order = Np-1,
    coriolis = (f₀ = 0, β = 0),
    buoyancy = (αᵀ = 0,),
    boundary_tags = ((0, 0), (1, 1), (1, 2)),
    boundary_conditions = (OceanBC(Impenetrable(FreeSlip()), Insulating()),
                           OceanBC(Penetrable(FreeSlip()), Insulating()))
)

# We prepare a callback that periodically fetches the horizontal velocity and
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

model.solver_configuration.timeend = 200

total_steps = ceil(Int, stop_time / time_step)
@info @sprintf("Running a simulation of the instability of the Bickley jet (Δt=%.2e, steps=%d)", time_step, total_steps)

try
    result = ClimateMachine.invoke!(model.solver_configuration;
                                    user_callbacks = [data_fetcher])
catch err
    @warn "Simulation ended prematurely because $(sprint(showerror, err))"
end
