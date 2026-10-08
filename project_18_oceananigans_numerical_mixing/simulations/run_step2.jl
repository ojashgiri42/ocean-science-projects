# run_step2.jl
# Step 2: the baseline lock exchange, following Ilicak et al. (2012): dx = 500 m, dz = 1 m (128 x 20 cells),
# lateral viscosity 1e-2 m^2/s (their base value), WENO (5th order) for buoyancy and momentum, 17 hours.
# Run from the project folder with:   julia --project=. simulations/run_step2.jl

include("lock_exchange.jl")

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo

result = run_lock_exchange(name = "step2_baseline_weno",
                           output_dir = OUTPUT_DIR,
                           Δx = 500.0,                        # horizontal grid spacing (m)
                           Δz = 1.0,                          # vertical grid spacing (m)
                           tracer_advection = WENO(order=5),
                           momentum_advection = WENO(order=5),
                           νh = 1e-2,                         # lateral viscosity (m^2/s), Ilicak's base case
                           νz = 1e-4,                         # vertical viscosity (m^2/s)
                           stop_time = 17hours,               # Ilicak's analysis time (s)
                           cfl = 0.2,                         # adaptive time step, advective CFL target
                           save_interval = 10minutes,         # snapshot interval (s)
                           consistency_tracer = true,         # also check that a uniform tracer stays at 1
                           labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => "WENO(order=5)"))
println(result)
