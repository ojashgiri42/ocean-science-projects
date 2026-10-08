# run_step5.jl
# Step 5: lateral viscosity and the grid Reynolds number. WENO (5th order) tracer advection on the baseline grid
# (dx = 500 m, dz = 1 m), lateral viscosity nu_h = 0.01 ... 200 m^2/s (the values of Ilicak et al. 2012, plus their
# base value 0.01), vertical viscosity 1e-4 m^2/s, 17 hours.
#   Set 1: centred (2nd order) momentum advection, which adds no damping of its own, so the explicit viscosity
#          alone controls the grid-scale noise in the velocity.
#   Set 2: WENO (5th order) momentum advection, which has built-in (implicit) damping. The nu_h = 0.01 case is the
#          step 2 baseline, so it is not repeated.
# Run from the project folder with:   julia --project=. simulations/run_step5.jl

include("lock_exchange.jl")

const OUTPUT_DIR = joinpath(homedir(), "oceananigans_runs", "numerical_mixing")   # raw output, outside the repo

viscosities = [0.01, 0.1, 1.0, 10.0, 100.0, 200.0]          # lateral viscosity nu_h (m^2/s)
momentum_sets = [("cmom", Centered(order=2), "Centered(order=2)"),   # set 1
                 ("wmom", WENO(order=5),     "WENO(order=5)")]       # set 2

for (short, momentum, label) in momentum_sets, νh in viscosities
    if short == "wmom" && νh == 0.01
        println("skipping WENO momentum, nu_h = 0.01: this is the step 2 baseline")
        continue
    end
    name = "step5_" * short * "_nu" * replace(string(νh), "." => "p")
    result = run_lock_exchange(name = name,
                               output_dir = OUTPUT_DIR,
                               Δx = 500.0, Δz = 1.0,                  # baseline grid (m)
                               tracer_advection = WENO(order=5),      # fixed tracer scheme
                               momentum_advection = momentum,         # set 1 or set 2
                               νh = νh, νz = 1e-4,                    # lateral viscosity varies (m^2/s)
                               stop_time = 17hours,
                               cfl = 0.2,
                               save_interval = 10minutes,
                               consistency_tracer = true,
                               labels = Dict("tracer_advection" => "WENO(order=5)", "momentum_advection" => label))
    println(result)
end
