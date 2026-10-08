# run_step3.jl
# Step 3: change ONLY the tracer (buoyancy) advection scheme. Everything else as in the step 2 baseline:
# dx = 500 m, dz = 1 m, WENO(order=5) momentum advection, nu_h = 1e-2 m^2/s, nu_z = 1e-4 m^2/s, 17 hours.
# The WENO(order=5) tracer case is the step 2 baseline run, so it is not repeated here.
# Run from the project folder with:   julia --project=. simulations/run_step3.jl

include("lock_exchange.jl")

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo

# (run name, tracer advection scheme, label)
cases = [("step3_tracer_centered2", Centered(order=2),     "Centered(order=2)"),       # 2nd-order centred, no upwinding
         ("step3_tracer_upwind3",   UpwindBiased(order=3), "UpwindBiased(order=3)")]   # 3rd-order upwind-biased

for (name, scheme, label) in cases
    result = run_lock_exchange(name = name,
                               output_dir = OUTPUT_DIR,
                               Δx = 500.0, Δz = 1.0,                 # grid spacing (m), as in step 2
                               tracer_advection = scheme,             # the only thing that changes
                               momentum_advection = WENO(order=5),    # fixed, as in step 2
                               νh = 1e-2, νz = 1e-4,                  # viscosities (m^2/s), as in step 2
                               stop_time = 17hours,
                               cfl = 0.2,
                               save_interval = 10minutes,
                               consistency_tracer = true,
                               labels = Dict("tracer_advection" => label, "momentum_advection" => "WENO(order=5)"))
    println(result)
end

# The centred run stops at t = 59 min under the default limit (b more than one DELTA_B outside the starting range).
# To see whether it really blows up or just carries large overshoots, I repeat it with a much wider limit.
# NaN and Inf still stop the run.
result = run_lock_exchange(name = "step3_tracer_centered2_wide_limit",
                           output_dir = OUTPUT_DIR,
                           Δx = 500.0, Δz = 1.0,
                           tracer_advection = Centered(order=2),
                           momentum_advection = WENO(order=5),
                           νh = 1e-2, νz = 1e-4,
                           stop_time = 17hours,
                           cfl = 0.2,
                           save_interval = 10minutes,
                           consistency_tracer = true,
                           b_range_limit = 10.0,                 # allow b down to -10 DELTA_B and up to 11 DELTA_B
                           labels = Dict("tracer_advection" => "Centered(order=2)", "momentum_advection" => "WENO(order=5)"))
println(result)
