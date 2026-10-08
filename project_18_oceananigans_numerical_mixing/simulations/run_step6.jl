# run_step6.jl
# Step 6: time step and stability. The baseline setup (WENO tracer and momentum advection, dx = 500 m, dz = 1 m,
# nu_h = 1e-2 m^2/s) with FIXED time steps chosen for advective CFL numbers of about 0.2, 0.5, 1.0 and 1.5, plus one
# extra adaptive run with a CFL target of 0.5 (the adaptive run with target 0.2 is the step 2 baseline).
# The fixed steps come from the baseline's largest value of u/dx + w/dz (0.0073083 1/s, reached about 36 minutes in):
# dt = CFL / 0.0073083. These runs check for NaN or out-of-range values at EVERY time step.
# Run from the project folder with:   julia --project=. simulations/run_step6.jl

include("lock_exchange.jl")

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo
const MAX_RATE = 0.0073083088            # largest u/dx + w/dz in the step 2 baseline (1/s)

target_cfls = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]   # 2.0 to 4.0 added after the 1.5 run still finished,
                                                        # to find where the run actually becomes unstable

for cfl_target in target_cfls
    Δt = cfl_target / MAX_RATE                                # fixed time step (s) for this CFL target
    name = "step6_fixed_cfl" * replace(string(cfl_target), "." => "p")
    result = run_lock_exchange(name = name,
                               output_dir = OUTPUT_DIR,
                               Δx = 500.0, Δz = 1.0,
                               tracer_advection = WENO(order=5),
                               momentum_advection = WENO(order=5),
                               νh = 1e-2, νz = 1e-4,
                               stop_time = 17hours,
                               fixed_Δt = Δt,                       # no wizard: one fixed time step
                               save_interval = 10minutes,
                               check_every = 1,                     # check every single step
                               consistency_tracer = true,
                               labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => "WENO(order=5)"))
    println(result)
end

# adaptive time stepping with a higher target than the baseline
result = run_lock_exchange(name = "step6_adaptive_cfl0p5",
                           output_dir = OUTPUT_DIR,
                           Δx = 500.0, Δz = 1.0,
                           tracer_advection = WENO(order=5),
                           momentum_advection = WENO(order=5),
                           νh = 1e-2, νz = 1e-4,
                           stop_time = 17hours,
                           cfl = 0.5,                                # TimeStepWizard target
                           save_interval = 10minutes,
                           check_every = 1,
                           consistency_tracer = true,
                           labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => "WENO(order=5)"))
println(result)

# The CFL 2.0 and 2.5 runs stopped on my out-of-range rule. To see whether they really blow up (NaN or Inf) or
# only carry large errors, I repeat them with the range check effectively switched off; only NaN/Inf
# (or a speed above 10 m/s) stops these runs.
for cfl_target in [2.0, 2.5]
    Δt = cfl_target / MAX_RATE
    name = "step6_fixed_cfl" * replace(string(cfl_target), "." => "p") * "_nan_only"
    result = run_lock_exchange(name = name,
                               output_dir = OUTPUT_DIR,
                               Δx = 500.0, Δz = 1.0,
                               tracer_advection = WENO(order=5),
                               momentum_advection = WENO(order=5),
                               νh = 1e-2, νz = 1e-4,
                               stop_time = 17hours,
                               fixed_Δt = Δt,
                               save_interval = 10minutes,
                               check_every = 1,
                               consistency_tracer = true,
                               b_range_limit = 1e6,                 # effectively no range check
                               labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => "WENO(order=5)"))
    println(result)
end
