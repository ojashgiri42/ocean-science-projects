# run_step4.jl
# Step 4: grid resolution. For the centred (2nd order) and WENO (5th order) tracer schemes, run the lock exchange
# at dx = 2000, 1000, 250 and 125 m, with dz = dx / 500 so every grid has the same cell shape (dx/dz = 500).
# The dx = 500 m runs are reused from steps 2 and 3. Everything else as in the baseline: WENO momentum advection,
# nu_h = 1e-2 m^2/s, nu_z = 1e-4 m^2/s, 17 hours. The centred runs use the wide stopping limit from step 3.
# Run from the project folder with:   julia --project=. simulations/run_step4.jl [first]
# With the argument "first", only the most expensive run (WENO, dx = 125 m) is done, to time it.

include("lock_exchange.jl")

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo

schemes = [("weno5",     WENO(order=5),     "WENO(order=5)",     1.0),     # (short name, scheme, label, b_range_limit)
           ("centered2", Centered(order=2), "Centered(order=2)", 10.0)]
spacings = [125.0, 250.0, 1000.0, 2000.0]                                  # dx (m); dz = dx / 500

only_first = length(ARGS) > 0 && ARGS[1] == "first"

for (short, scheme, label, limit) in schemes, Δx in spacings
    name = "step4_" * short * "_dx" * string(round(Int, Δx))
    if only_first && !(short == "weno5" && Δx == 125.0)
        continue
    end
    if !only_first && short == "weno5" && Δx == 125.0 && isfile(joinpath(OUTPUT_DIR, name, "run_info.txt"))
        println("skipping ", name, " (already done in the timing run)")
        continue
    end
    result = run_lock_exchange(name = name,
                               output_dir = OUTPUT_DIR,
                               Δx = Δx, Δz = Δx / 500,                # same cell shape on every grid (m)
                               tracer_advection = scheme,
                               momentum_advection = WENO(order=5),    # fixed, as in steps 2 and 3
                               νh = 1e-2, νz = 1e-4,                  # viscosities (m^2/s)
                               stop_time = 17hours,
                               cfl = 0.2,
                               save_interval = 10minutes,
                               consistency_tracer = true,
                               b_range_limit = limit,
                               labels = Dict("tracer_advection" => label, "momentum_advection" => "WENO(order=5)"))
    println(result)
end
