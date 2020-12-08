# Unstable Bickley jet in GeophysicalFlows

using FourierFlows, Printf, Random, Plots
 
using Random: seed!
using FFTW: rfft, irfft

import GeophysicalFlows.TwoDNavierStokes
import GeophysicalFlows.TwoDNavierStokes: energy, enstrophy
import GeophysicalFlows: peakedisotropicspectrum

include("Bickley.jl")

using .Bickley

dev = CPU()

n, L  = 128, 4π

    dt = 1e-2  # timestep
nsteps = 4000  # total number of steps
 nsubs = 20    # number of steps between each plot

problem = TwoDNavierStokes.Problem(dev;
                                   nx = n,
                                   Lx = 4π,
                                   ny = n,
                                   Ly = 4π,
                                   dt = dt,
                                   stepper = "FilteredRK4")

ϵ = 0.1
k = 0.5
ℓ = 0.5

grid = problem.grid
x, y = grid.x, grid.y

ψᵢ = @. Bickley.Ψ(y) + ϵ * Bickley.ψ̃(x, y, k, ℓ)

ψᵢh = zeros(Complex{Float64}, grid.nkr, grid.nl)
ζᵢ = zeros(Float64, grid.nx, grid.ny)

mul!(ψᵢh, grid.rfftplan, ψ)

ζh = @. grid.invKrsq * ψh
ζh[1, 1] = 0

mul(ζᵢ

#=
TwoDNavierStokes.set_zeta!(prob, ζ₀)

# Let's plot the initial vorticity field:
heatmap(x, y, vs.zeta,
         aspectratio = 1,
              c = :balance,
           clim = (-40, 40),
          xlims = (-L/2, L/2),
          ylims = (-L/2, L/2),
         xticks = -3:3,
         yticks = -3:3,
         xlabel = "x",
         ylabel = "y",
          title = "initial vorticity",
     framestyle = :box)
            

#=
# ## Diagnostics

# Create Diagnostics -- `energy` and `enstrophy` functions are imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
diags = [E, Z] # A list of Diagnostics types passed to "stepforward!" will  be updated every timestep.
nothing # hide


# ## Output

# We choose folder for outputing `.jld2` files and snapshots (`.png` files).
filepath = "."
plotpath = "./plots_decayingTwoDNavierStokes"
plotname = "snapshots"
filename = joinpath(filepath, "decayingTwoDNavierStokes.jld2")
nothing # hide

# Do some basic file management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end
nothing # hide

# And then create Output
get_sol(prob) = Array(prob.sol) # extracts the Fourier-transformed solution
get_u(prob) = Array(irfft(im*gr.l.*gr.invKrsq.*sol, gr.nx))
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))
saveproblem(out)
nothing # hide


# ## Visualizing the simulation

# We initialize a plot with the vorticity field and the time-series of
# energy and enstrophy diagnostics.

p1 = heatmap(x, y, vs.zeta,
         aspectratio = 1,
                   c = :balance,
                clim = (-40, 40),
               xlims = (-L/2, L/2),
               ylims = (-L/2, L/2),
              xticks = -3:3,
              yticks = -3:3,
              xlabel = "x",
              ylabel = "y",
               title = "vorticity, t="*@sprintf("%.2f", cl.t),
          framestyle = :box)

p2 = plot(2, # this means "a plot with two series"
               label = ["energy E(t)/E(0)" "enstrophy Z(t)/Z(0)"],
              legend = :right,
           linewidth = 2,
               alpha = 0.7,
              xlabel = "t",
               xlims = (0, 41),
               ylims = (0, 1.1))

l = @layout grid(1, 2)
p = plot(p1, p2, layout = l, size = (900, 400))


# ## Time-stepping the `Problem` forward

# We time-step the `Problem` forward in time.

startwalltime = time()

anim = @animate for j = 0:Int(nsteps/nsubs)
    
  log = @sprintf("step: %04d, t: %d, ΔE: %.4f, ΔZ: %.4f, walltime: %.2f min",
      cl.step, cl.t, E.data[E.i]/E.data[1], Z.data[Z.i]/Z.data[1], (time()-startwalltime)/60)
  
  if j%(1000/nsubs)==0; println(log) end  

  p[1][1][:z] = vs.zeta
  p[1][:title] = "vorticity, t="*@sprintf("%.2f", cl.t)
  push!(p[2][1], E.t[E.i], E.data[E.i]/E.data[1])
  push!(p[2][2], Z.t[Z.i], Z.data[Z.i]/Z.data[1])

  stepforward!(prob, diags, nsubs)
  TwoDNavierStokes.updatevars!(prob)  
  
end

mp4(anim, "twodturb.mp4", fps=18)


# Last we save the output.
saveoutput(out)


# ## Radial energy spectrum

# After the simulation is done we plot the radial energy spectrum to illustrate
# how `FourierFlows.radialspectrum` can be used,

E  = @. 0.5*(vs.u^2 + vs.v^2) # energy density
Eh = rfft(E)                  # Fourier transform of energy density
kr, Ehr = FourierFlows.radialspectrum(Eh, gr, refinement=1) # compute radial specturm of `Eh`
nothing # hide

# and we plot it.
plot(kr, abs.(Ehr),
    linewidth = 2,
        alpha = 0.7,
       xlabel = "kᵣ", ylabel = "∫ |Ê| kᵣ dk_θ",
        xlims = (5e-1, gr.nx),
       xscale = :log10, yscale = :log10,
        title = "Radial energy spectrum",
       legend = false)
       
=#
=#
