# run_step1.jl
# Step 1: a tiny, short lock exchange to test the whole chain (Julia run -> NetCDF file -> Python).
# Run from the project folder with:   julia --project=. simulations/run_step1.jl

include("lock_exchange.jl")                       # brings in run_lock_exchange and the constants

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo

result = run_lock_exchange(name = "step1_tiny",
                           output_dir = OUTPUT_DIR,
                           Δx = 2000.0,                       # 32 cells along the channel (m)
                           Δz = 4.0,                          # 5 cells in the vertical (m)
                           tracer_advection = WENO(order=5),
                           momentum_advection = WENO(order=5),
                           νh = 1e-2,                         # Ilicak's base lateral viscosity (m^2/s)
                           νz = 1e-4,                         # vertical viscosity (m^2/s)
                           stop_time = 2hours,                # short test (s)
                           cfl = 0.2,
                           save_interval = 10minutes,
                           labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => "WENO(order=5)"))
println(result)
