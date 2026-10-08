# Minimal example: on a grid with a Flat y-direction, the diffusive time-step limit for a
# HorizontalScalarDiffusivity uses Δy = 1 m instead of ignoring the Flat direction.
#
# Run with:  julia --project=. simulations/mwe_flat_y_diffusion_timescale.jl
# Found with Oceananigans 0.113.5 on Julia 1.13.1 (macOS, Apple M1, CPU). Not reported upstream (yet).

using Oceananigans
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures: cell_diffusion_timescale

Δx = 2000.0                                   # horizontal grid spacing (m)
ν = 0.01                                      # horizontal viscosity (m^2/s)

grid = RectilinearGrid(CPU(); size = (32, 5), x = (0, 32Δx), z = (-20, 0),
                       topology = (Bounded, Flat, Bounded))           # y is Flat: there is no y-direction

model = HydrostaticFreeSurfaceModel(grid; closure = HorizontalScalarDiffusivity(ν = ν))

println("Oceananigans ", pkgversion(Oceananigans), ", Julia ", VERSION)
println("grid topology:                         ", Oceananigans.Grids.topology(grid))
println("minimum_xspacing:                      ", minimum_xspacing(grid, Center(), Center(), Center()), " m")
println("minimum_yspacing (Flat direction):     ", minimum_yspacing(grid, Center(), Center(), Center()), " m")
println("cell_diffusion_timescale from model:   ", cell_diffusion_timescale(model), " s")
println("expected Δx^2 / ν (y is Flat):         ", Δx^2 / ν, " s")

# Consequence for the TimeStepWizard: with diffusive_cfl = 0.2 the time step is capped at
# 0.2 * (1 m)^2 / ν = 20 s, although the real horizontal cells are 2000 m wide.
simulation = Simulation(model; Δt = 1000.0, stop_iteration = 1)
wizard = TimeStepWizard(cfl = 0.2, diffusive_cfl = 0.2, max_change = Inf, min_change = 0.0)
wizard(simulation)
println("time step after the wizard (Δt was 1000 s): ", simulation.Δt, " s")
