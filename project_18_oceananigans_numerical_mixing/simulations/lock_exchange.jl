# lock_exchange.jl
#
# One plain function, run_lock_exchange, that builds and runs one 2D (x-z) lock-exchange simulation
# with Oceananigans and saves its output. The run scripts for each step call it with different settings.
#
# Setup follows Ilicak et al. (2012, Ocean Modelling 45-46, Section 3.1):
#   - channel 64 km long and 20 m deep, no rotation, closed (no-flow, free-slip) walls
#   - dense water on the left half, light water on the right half, density difference 5 kg/m^3
#   - zero explicit tracer diffusivity, vertical viscosity 1e-4 m^2/s, lateral viscosity varied
# I use buoyancy b directly instead of temperature: b = -g (rho - rho0) / rho0, with rho0 = 1000 kg/m^3.
# The dense water has b = 0 and the light water b = g * 5 / 1000 = 0.04905 m/s^2.
#
# Units: metres, seconds, m/s^2 for buoyancy.

using Oceananigans
using Oceananigans.Units
using Printf
using NCDatasets                          # loads Oceananigans' NetCDF writer (a package extension)

const GRAVITY = 9.81                    # gravitational acceleration (m/s^2)
const RHO_0 = 1000.0                    # reference density (kg/m^3)
const DELTA_RHO = 5.0                   # density difference between the two water masses (kg/m^3)
const DELTA_B = GRAVITY * DELTA_RHO / RHO_0   # buoyancy difference (m/s^2) = reduced gravity g'
const LENGTH = 64kilometers             # channel length (m)
const DEPTH = 20meters                  # channel depth (m)
const FREE_SURFACE_GRAVITY = 9.81e5    # gravity used ONLY by the free surface (m/s^2), 100 000 x the real value.
                                        # A stiff free surface barely moves (max |eta| ~ 3e-7 m), so the model acts
                                        # like a rigid lid and total buoyancy is conserved to ~1e-8. With the real g,
                                        # buoyancy leaked through the moving surface by ~1e-3 (see TUTORIAL step 1).

"""
    run_lock_exchange(; name, output_dir, kwargs...)

Build and run one lock-exchange simulation. Keyword arguments:

  name                 run name, also the folder name inside output_dir
  output_dir           parent folder for raw output (outside the git repository)
  Δx, Δz               grid spacing in x and z (m)
  tracer_advection     advection scheme for buoyancy, e.g. WENO(order=5), Centered(order=2)
  momentum_advection   advection scheme for velocity
  νh, νz               lateral and vertical viscosity (m^2/s); tracer diffusivity is always zero
  stop_time            model time to stop at (s)
  cfl                  target advective CFL number for adaptive time stepping (ignored if fixed_Δt is set)
  fixed_Δt             a fixed time step (s); `nothing` means adaptive (TimeStepWizard)
  save_interval        time between saved snapshots of b, u, w (s)
  check_every          how often (in iterations) to check for NaN or out-of-range values and log statistics
  labels               a Dict of text labels saved with the run (scheme names etc.)
  b_range_limit        stop if b goes further than b_range_limit x DELTA_B outside the starting range [0, DELTA_B]
                       (default 1); NaN and Inf always stop the run
  consistency_tracer   if true, add a passive tracer c that starts at exactly 1 everywhere; a consistent
                       advection scheme must keep it at exactly 1, so its largest deviation is recorded

Returns a NamedTuple with the run status and timing. Writes into output_dir/name/:
  fields.nc     snapshots of b, u, w (64-bit floats)
  log.csv       statistics every check_every iterations
  run_info.txt  settings, final status, iterations and wall-clock time
"""
function run_lock_exchange(; name,
                             output_dir,
                             Δx = 500.0,
                             Δz = 1.0,
                             tracer_advection = WENO(order=5),
                             momentum_advection = WENO(order=5),
                             νh = 1e-2,
                             νz = 1e-4,
                             stop_time = 17hours,
                             cfl = 0.2,
                             fixed_Δt = nothing,
                             save_interval = 10minutes,
                             check_every = 10,
                             labels = Dict{String, String}(),
                             consistency_tracer = false,
                             b_range_limit = 1.0)

    run_dir = joinpath(output_dir, name)
    mkpath(run_dir)                                            # create the run folder if needed

    # ---- grid: 2D x-z, uniform spacing, closed walls ------------------------------------------
    Nx = round(Int, LENGTH / Δx)                               # number of cells along the channel
    Nz = round(Int, DEPTH / Δz)                                # number of cells in the vertical
    grid = RectilinearGrid(CPU();
                           size = (Nx, Nz),
                           x = (0, LENGTH),
                           z = (-DEPTH, 0),                    # fixed vertical grid (see TUTORIAL step 1)
                           halo = (5, 5),                      # wide enough for 5th-order WENO
                           topology = (Bounded, Flat, Bounded))   # walls at both ends, no y-direction

    # ---- physics: viscosity only, NO tracer diffusivity ----------------------------------------
    closure = (HorizontalScalarDiffusivity(ν = νh, κ = 0),     # lateral viscosity, zero lateral diffusivity
               VerticalScalarDiffusivity(ν = νz, κ = 0))       # vertical viscosity, zero vertical diffusivity

    model = HydrostaticFreeSurfaceModel(grid;
                                        tracers = consistency_tracer ? (:b, :c) : (:b,),   # c: passive check tracer
                                        buoyancy = BuoyancyTracer(),       # density depends on b only (linear)
                                        momentum_advection = momentum_advection,
                                        tracer_advection = tracer_advection,
                                        closure = closure,
                                        free_surface = ImplicitFreeSurface(gravitational_acceleration = FREE_SURFACE_GRAVITY))

    # ---- initial state: dense water left (b = 0), light water right (b = DELTA_B), at rest ------
    initial_b(x, z) = x < LENGTH / 2 ? 0.0 : DELTA_B
    set!(model, b = initial_b)
    if consistency_tracer
        set!(model, c = 1.0)                                   # uniform passive tracer, must stay exactly 1
    end

    # ---- time stepping ------------------------------------------------------------------------
    if fixed_Δt === nothing
        simulation = Simulation(model; Δt = 1.0, stop_time = stop_time)   # start small; the wizard adjusts it
        # The wizard's own viscous check treats the Flat y-direction as 1 m wide, which would cap the step
        # at 0.2 * (1 m)^2 / νh (20 s for νh = 0.01). So I switch that check off and give the wizard the real
        # viscous limit myself, from the actual cell sizes (see TUTORIAL step 1).
        viscous_limit = 0.2 * min(Δx^2 / νh, Δz^2 / νz)       # explicit viscosity is stable below ~ Δ^2 / (2ν)
        conjure_time_step_wizard!(simulation, IterationInterval(10);
                                  cfl = cfl,                               # advective CFL target
                                  diffusive_cfl = Inf,                     # wizard's viscous check off (see above)
                                  max_Δt = viscous_limit,                  # real viscous limit instead
                                  max_change = 1.1)
    else
        simulation = Simulation(model; Δt = fixed_Δt, stop_time = stop_time)   # one fixed time step
    end

    # ---- output: snapshots of b, u, w as 64-bit floats ---------------------------------------
    b = model.tracers.b
    u, v, w = model.velocities
    simulation.output_writers[:fields] = NetCDFWriter(model, (; b, u, w);
                                                      filename = "fields.nc",
                                                      dir = run_dir,
                                                      schedule = TimeInterval(save_interval),
                                                      array_type = Array{Float64},   # full precision for RPE
                                                      overwrite_files = true)        # replace output from an earlier run of the same name

    # ---- health check and log, every check_every iterations ----------------------------------
    cell_volume = Δx * Δz                                      # per metre of width (2D), the same for every cell
    log_file = open(joinpath(run_dir, "log.csv"), "w")
    println(log_file, "iteration,time_s,dt_s,b_min,b_max,b_total_m3_per_s2,u_max_abs,w_max_abs,advective_cfl")
    status = Ref("running")
    c_deviation = Ref(0.0)                                     # largest |c - 1| seen so far (if c is used)

    function check_and_log(sim)
        bi = interior(sim.model.tracers.b)
        ui = interior(sim.model.velocities.u)
        wi = interior(sim.model.velocities.w)
        b_min, b_max = minimum(bi), maximum(bi)
        u_max, w_max = maximum(abs, ui), maximum(abs, wi)
        dt = sim.Δt
        courant = dt * (u_max / Δx + w_max / Δz)              # advective CFL number (upper bound)
        if consistency_tracer
            c_deviation[] = max(c_deviation[], maximum(abs, interior(sim.model.tracers.c) .- 1))
        end
        @printf(log_file, "%d,%.6f,%.6f,%.12e,%.12e,%.15e,%.6e,%.6e,%.6f\n",
                iteration(sim), time(sim), dt, b_min, b_max, sum(bi) * cell_volume, u_max, w_max, courant)
        # stop cleanly if anything is not a number or clearly outside a physical range
        if !(isfinite(b_min) && isfinite(b_max) && isfinite(u_max) && isfinite(w_max))
            status[] = @sprintf("blew up: NaN or Inf at t = %.1f s (iteration %d)", time(sim), iteration(sim))
            sim.running = false
        elseif b_min < -b_range_limit * DELTA_B || b_max > (1 + b_range_limit) * DELTA_B || u_max > 10.0
            status[] = @sprintf("blew up: out of range at t = %.1f s (b in [%.3e, %.3e], max|u| = %.2f m/s)",
                                time(sim), b_min, b_max, u_max)
            sim.running = false
        end
        return nothing
    end
    add_callback!(simulation, check_and_log, IterationInterval(check_every))

    # ---- run ---------------------------------------------------------------------------------
    start = time_ns()
    run!(simulation)
    wall_time = 1e-9 * (time_ns() - start)
    close(log_file)
    if status[] == "running"
        status[] = "completed"
    end

    # ---- record settings and outcome --------------------------------------------------------
    open(joinpath(run_dir, "run_info.txt"), "w") do io
        println(io, "name = ", name)
        println(io, "status = ", status[])
        println(io, "oceananigans_version = ", pkgversion(Oceananigans))
        println(io, "julia_version = ", VERSION)
        println(io, "Nx = ", Nx, "\nNz = ", Nz, "\ndx_m = ", Δx, "\ndz_m = ", Δz)
        println(io, "tracer_advection = ", get(labels, "tracer_advection", summary(tracer_advection)))
        println(io, "momentum_advection = ", get(labels, "momentum_advection", summary(momentum_advection)))
        println(io, "nu_h_m2_s = ", νh, "\nnu_z_m2_s = ", νz, "\ntracer_diffusivity_m2_s = 0")
        println(io, "time_stepping = ", fixed_Δt === nothing ? @sprintf("adaptive, cfl = %.2f", cfl) : @sprintf("fixed, dt = %.4f s", fixed_Δt))
        println(io, "stop_time_s = ", stop_time, "\nmodel_time_s = ", time(simulation))
        println(io, "iterations = ", iteration(simulation))
        println(io, @sprintf("wall_time_s = %.2f", wall_time))
        println(io, "delta_b_m_s2 = ", DELTA_B)
        println(io, "b_range_limit = ", b_range_limit)
        println(io, "free_surface = ", nameof(typeof(model.free_surface)),
                ", gravitational_acceleration = ", model.free_surface.gravitational_acceleration, " m/s^2")
        # an EXPLICIT free surface would need Δt < Δx / sqrt(g H) for its surface gravity waves;
        # the implicit one does not, so this number only shows how far above it the actual steps were
        println(io, @sprintf("explicit_surface_wave_dt_limit_s = %.4f", Δx / sqrt(FREE_SURFACE_GRAVITY * DEPTH)))
        if consistency_tracer
            println(io, @sprintf("consistency_tracer_max_deviation = %.3e", c_deviation[]))
        end
    end

    @info @sprintf("%s: %s after %d iterations, model time %.0f s, wall time %.1f s",
                   name, status[], iteration(simulation), time(simulation), wall_time)
    return (; name, status = status[], iterations = iteration(simulation), model_time = time(simulation), wall_time)
end
