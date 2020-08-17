using Dates
using Logging
using Printf
using LinearAlgebra

using MPI
using StaticArrays
using NCDatasets

using ClimateMachine
ClimateMachine.init(parse_clargs = true)

using ClimateMachine.BalanceLaws: vars_state_conservative, vars_state_auxiliary
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Filters
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates: flattenednames
using ClimateMachine.SplitExplicit01
using ClimateMachine.GenericCallbacks
using ClimateMachine.VTK

using CLIMAParameters
using CLIMAParameters.Planet: grav

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

import ClimateMachine.SplitExplicit01:
    ocean_init_aux!,
    ocean_init_state!,
    ocean_boundary_state!,
    CoastlineFreeSlip,
    CoastlineNoSlip,
    OceanFloorFreeSlip,
    OceanFloorNoSlip,
    OceanSurfaceNoStressNoForcing,
    OceanSurfaceStressNoForcing,
    OceanSurfaceNoStressForcing,
    OceanSurfaceStressForcing

import ClimateMachine.DGMethods:
    update_auxiliary_state!,
    update_auxiliary_state_gradient!,
    vars_state_conservative,
    vars_state_auxiliary,
    VerticalDirection

const ArrayType = ClimateMachine.array_type()

struct SimpleBox{T} <: AbstractOceanProblem
    Lˣ::T
    Lʸ::T
    H ::T
    τₒ::T
    λʳ::T
    θᴱ::T
end

@inline function ocean_boundary_state!(m::OceanModel, p::SimpleBox, bctype, x...)
    if bctype == 1
        ocean_boundary_state!(m, CoastlineNoSlip(), x...)
    elseif bctype == 2
        ocean_boundary_state!(m, OceanFloorNoSlip(), x...)
    elseif bctype == 3
        ocean_boundary_state!(m, OceanSurfaceStressForcing(), x...)
    end
end

@inline function ocean_boundary_state!(m::Continuity3dModel, p::SimpleBox, bctype, x...)
   #if bctype == 1
        ocean_boundary_state!(m, CoastlineNoSlip(), x...)
   #end
end

@inline function ocean_boundary_state!(m::BarotropicModel, p::SimpleBox, bctype, x...)
    return ocean_boundary_state!(m, CoastlineNoSlip(), x...)
end

function ocean_init_state!(p::SimpleBox, Q, A, coords, t)
    @inbounds y = coords[2]
    @inbounds z = coords[3]
    @inbounds H = p.H

    Q.u = @SVector [-0, -0]
    Q.η = -0
    Q.θ = (5 + 4 * cos(y * π / p.Lʸ)) * (1 + z / H)

    return nothing
end

function ocean_init_aux!(m::OceanModel, p::SimpleBox, A, geom)
    FT = eltype(A)
    @inbounds A.y = geom.coord[2]

    # not sure if this is needed but getting weird intialization stuff
    A.w = -0
    A.pkin = -0
    A.wz0 = -0
    A.u_d = @SVector [-0, -0]
    A.ΔGu = @SVector [-0, -0]

    return nothing
end

# A is Filled afer the state
function ocean_init_aux!(m::BarotropicModel, P::SimpleBox, A, geom)
    @inbounds A.y = geom.coord[2]

    A.Gᵁ = @SVector [-0, -0]
    A.U_c = @SVector [-0, -0]
    A.η_c = -0
    A.U_s = @SVector [-0, -0]
    A.η_s = -0
    A.Δu = @SVector [-0, -0]
    A.η_diag = -0
    A.Δη = -0

    return nothing
end

function make_callbacks(step, nout, mpicomm, odesolver, dg_slow, model_slow, Q_slow, dg_fast, model_fast, Q_fast)
    Np = N
    Nᵖ⁺¹ = Np + 1

    ΣNˣ = (Np+1)*Nˣ
    ΣNʸ = (Np+1)*Nʸ
    ΣNᶻ = (Np+1)*Nᶻ

    gnd = reshape(dg_slow.grid.vgeo, (Np+1, Np+1, Np+1, 16, Nᶻ, Nʸ, Nˣ))
    x = gnd[:, :, :, 13, :, :, :] |> Array
    y = gnd[:, :, :, 14, :, :, :] |> Array
    z = gnd[:, :, :, 15, :, :, :] |> Array

    ΔX = (maximum(x) - minimum(x)) / Nˣ
    ΔY = (maximum(y) - minimum(y)) / Nʸ
    ΔZ = (maximum(z) - minimum(z)) / Nᶻ

    ds = NCDataset("simple_box_ivd.nc", "c")

    defDim(ds, "time", Inf)
    defVar(ds, "time", Float64, ("time",))

    defDim(ds, "x", ΣNˣ)
    defDim(ds, "y", ΣNʸ)
    defDim(ds, "z", ΣNᶻ)

    defVar(ds, "x_nodal", Float64, ("x", "y", "z"))
    defVar(ds, "y_nodal", Float64, ("x", "y", "z"))
    defVar(ds, "z_nodal", Float64, ("x", "y", "z"))

    defVar(ds, "u", Float64, ("x", "y", "z", "time"))
    defVar(ds, "v", Float64, ("x", "y", "z", "time"))
    defVar(ds, "η", Float64, ("x", "y", "z", "time"))
    defVar(ds, "θ", Float64, ("x", "y", "z", "time"))

    netcdf_output = GenericCallbacks.EveryXSimulationSteps(1) do (init = false)
        time_index = length(ds["time"]) + 1

        if time_index == 1
            xs = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
            ys = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
            zs = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
        end

        Q3nd = reshape(Q_slow.realdata, (Np+1, Np+1, Np+1, 4, Nᶻ, Nʸ, Nˣ))
        u = Q3nd[:, :, :, 1, :, :, :] |> Array
        v = Q3nd[:, :, :, 2, :, :, :] |> Array
        η = Q3nd[:, :, :, 3, :, :, :] |> Array
        θ = Q3nd[:, :, :, 4, :, :, :] |> Array

        us = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
        vs = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
        ηs = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)
        θs = zeros(ΣNˣ, ΣNʸ, ΣNᶻ)

        for I in 1:Nˣ, J in 1:Nʸ, K in 1:Nᶻ
            I′ = x[:, :, :, K, J, I] ./ ΔX |> maximum |> round |> Int
            J′ = y[:, :, :, K, J, I] ./ ΔY |> maximum |> round |> Int
            K′ = z[:, :, :, K, J, I] ./ ΔZ |> maximum |> round |> Int
            K′ = Nᶻ + K′

            @debug "($I, $J, $K) -> ($I′, $J′, $(Nᶻ+K′))"

            i_elem = (I′-1) * Nᵖ⁺¹ + 1 : I′ * Nᵖ⁺¹
            j_elem = (J′-1) * Nᵖ⁺¹ + 1 : J′ * Nᵖ⁺¹
            k_elem = (K′-1) * Nᵖ⁺¹ + 1 : K′ * Nᵖ⁺¹

            if time_index == 1
                xs[i_elem, j_elem, k_elem] .= x[:, :, :, K, J, I]
                ys[i_elem, j_elem, k_elem] .= y[:, :, :, K, J, I]
                zs[i_elem, j_elem, k_elem] .= z[:, :, :, K, J, I]
            end

            us[i_elem, j_elem, k_elem] .= u[:, :, :, K, J, I]
            vs[i_elem, j_elem, k_elem] .= v[:, :, :, K, J, I]
            ηs[i_elem, j_elem, k_elem] .= η[:, :, :, K, J, I]
            θs[i_elem, j_elem, k_elem] .= θ[:, :, :, K, J, I]
        end

        if time_index == 1
            ds["x_nodal"][:, :, :] = xs
            ds["y_nodal"][:, :, :] = ys
            ds["z_nodal"][:, :, :] = zs
        end

        ds["u"][:, :, :, time_index] = us
        ds["v"][:, :, :, time_index] = vs
        ds["η"][:, :, :, time_index] = ηs
        ds["θ"][:, :, :, time_index] = θs

        sync(ds)
    end

    starttime = Ref(now())
    cbinfo = GenericCallbacks.EveryXWallTimeSeconds(60, mpicomm) do (s = false)
        if s
            starttime[] = now()
        else
            energy = norm(Q_slow)
            @info @sprintf(
                """Update
                simtime = %8.2f / %8.2f
                runtime = %s
                norm(Q) = %.16e""",
                ODESolvers.gettime(odesolver), timeend,
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy)
        end
    end

    return (netcdf_output, cbinfo)
end

#################
# RUN THE TESTS #
#################
FT = Float64

const timeend = 5 * 24 * 3600 # s
const tout = 24 * 3600 # s

const N = 4
const Nˣ = 20
const Nʸ = 20
const Nᶻ = 20
const Lˣ = 4e6  # m
const Lʸ = 4e6  # m
const H = 1000  # m

xrange = range(FT(0); length = Nˣ + 1, stop = Lˣ)
yrange = range(FT(0); length = Nʸ + 1, stop = Lʸ)
zrange = range(FT(-H); length = Nᶻ + 1, stop = 0)

#const cʰ = sqrt(gravity * H)
const cʰ = 1  # typical of ocean internal-wave speed
const cᶻ = 0

#- inverse ratio of additional fast time steps (for weighted average)
#  --> do 1/add more time-steps and average from: 1 - 1/add up to: 1 + 1/add
# e.g., = 1 --> 100% more ; = 2 --> 50% more ; = 3 --> 33% more ...
add_fast_substeps = 2

#- number of Implicit vertical-diffusion sub-time-steps within one model full time-step
# default = 0 : disable implicit vertical diffusion
numImplSteps = 5

#const τₒ = 2e-1  # (Pa = N/m^2)
# since we are using old BC (with factor of 2), take only half:
const τₒ = 1e-1
const λʳ = 10 // 86400 # m/s
const θᴱ = 10    # deg.C

mpicomm = MPI.COMM_WORLD

ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
loglevel = ll == "DEBUG" ? Logging.Debug :
           ll == "WARN"  ? Logging.Warn  :
           ll == "ERROR" ? Logging.Error : Logging.Info
logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
global_logger(ConsoleLogger(logger_stream, loglevel))

brickrange_2D = (xrange, yrange)
topl_2D = BrickTopology(mpicomm, brickrange_2D, periodicity = (false, false))
grid_2D = DiscontinuousSpectralElementGrid(topl_2D, FloatType = FT, DeviceArray = ArrayType,
                                           polynomialorder = N)

brickrange_3D = (xrange, yrange, zrange)
topl_3D = StackedBrickTopology(mpicomm, brickrange_3D; periodicity = (false, false, false),
                               boundary = ((1, 1), (1, 1), (2, 3)))
grid_3D = DiscontinuousSpectralElementGrid(topl_3D, FloatType = FT, DeviceArray = ArrayType, polynomialorder = N)

prob = SimpleBox{FT}(Lˣ, Lʸ, H, τₒ, λʳ, θᴱ)
gravity = grav(param_set)

#- set model time-step:
dt_fast = 240
dt_slow = 5400

nout = ceil(Int64, tout / dt_slow)
dt_slow = tout / nout
numImplSteps > 0 ? ivdc_dt = dt_slow / FT(numImplSteps) : ivdc_dt = dt_slow

model = OceanModel{FT}(prob, grav = gravity, cʰ = cʰ, add_fast_substeps = add_fast_substeps,
                       numImplSteps = numImplSteps, ivdc_dt = ivdc_dt, κᶜ = FT(0.1))

barotropicmodel = BarotropicModel(model)

minΔx = min_node_distance(grid_3D, HorizontalDirection())
minΔz = min_node_distance(grid_3D, VerticalDirection())

#- 2 horiz directions
gravity_max_dT = 1 / ( 2 * sqrt(gravity * H) / minΔx )

#- 2 horiz directions + harmonic visc or diffusion: 2^2 factor in CFL:
viscous_max_dT = 1 / ( 2 * model.νʰ / minΔx^2 + model.νᶻ / minΔz^2 )/ 4
diffusive_max_dT = 1 / ( 2 * model.κʰ / minΔx^2 + model.κᶻ / minΔz^2 )/ 4

@info @sprintf(
    """Update
       Gravity Max-dT = %.1f
       Timestep       = %.1f""",
    gravity_max_dT, dt_fast)

@info @sprintf(
    """Update
   Viscous   Max-dT = %.1f
   Diffusive Max-dT = %.1f
   Timestep      = %.1f""",
    viscous_max_dT, diffusive_max_dT, dt_slow)

dg = OceanDGModel(
    model, grid_3D,
#   CentralNumericalFluxFirstOrder(),
    RusanovNumericalFlux(),
    CentralNumericalFluxSecondOrder(),
    CentralNumericalFluxGradient(),
)

barotropic_dg = DGModel(
    barotropicmodel, grid_2D,
#   CentralNumericalFluxFirstOrder(),
    RusanovNumericalFlux(),
    CentralNumericalFluxSecondOrder(),
    CentralNumericalFluxGradient(),
)

Q_3D = init_ode_state(dg, FT(0); init_on_cpu = true)
Q_2D = init_ode_state(barotropic_dg, FT(0); init_on_cpu = true)

lsrk_ocean = LSRK54CarpenterKennedy(dg, Q_3D, dt = dt_slow, t0 = 0)
lsrk_barotropic = LSRK54CarpenterKennedy(barotropic_dg, Q_2D, dt = dt_fast, t0 = 0)

odesolver = SplitExplicitLSRK2nSolver(lsrk_ocean, lsrk_barotropic)

# Set up State Check call back for config state arrays, called every ntFreq time steps
ntFreq=1
cbcs_dg=ClimateMachine.StateCheck.sccreate(
        [(Q_3D,"oce Q_3D",),
         (dg.state_auxiliary,"oce aux",),
    #    (dg.diffstate,"oce diff",),
    #    (lsrk_ocean.dQ,"oce_dQ",),
    #    (dg.modeldata.tendency_dg.state_auxiliary,"tend Int aux",),
    #    (dg.modeldata.conti3d_Q,"conti3d_Q",),
         (Q_2D,"baro Q_2D",),
         (barotropic_dg.state_auxiliary ,"baro aux",)
    ],
        ntFreq; prec=12)

step = [0, 0]
cbvector = make_callbacks(step, nout, mpicomm, odesolver, dg, model, Q_3D,
      barotropic_dg, barotropicmodel, Q_2D)

eng0 = norm(Q_3D)

@info @sprintf """Starting
norm(Q₀) = %.16e
ArrayType = %s""" eng0 ArrayType

# slow fast state tuple
Qvec = (slow = Q_3D, fast = Q_2D)
# solve!(Qvec, odesolver; timeend = timeend, callbacks = cbvector)
cbv=(cbvector...,cbcs_dg)
solve!(Qvec, odesolver; timeend = timeend, callbacks = cbv)

## Enable the code block below to print table for use in reference value code
## reference value code sits in a file named $(@__FILE__)_refvals.jl. It is hand
## edited using code generated by block below when reference values are updated.
regenRefVals = false
if regenRefVals
    ## Print state statistics in format for use as reference values
    println(
        "# SC ========== Test number ", 1,
        " reference values and precision match template. =======",
    )
    println("# SC ========== $(@__FILE__) test reference values ======================================")
    ClimateMachine.StateCheck.scprintref(cbcs_dg)
    println("# SC ====================================================================================")
end

## Check results against reference if present
checkRefVals = true
if checkRefVals
    include("simple_box_ivd_refvals.jl")
    refDat = (refVals[1], refPrecs[1])
    checkPass = ClimateMachine.StateCheck.scdocheck(cbcs_dg, refDat)
    checkPass ? checkRep = "Pass" : checkRep = "Fail"
    @info @sprintf("""Compare vs RefVals: %s""", checkRep )
end
