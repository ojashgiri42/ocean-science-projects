# Tutorial: numerical mixing in the lock-exchange test (Oceananigans.jl)

These notes go through each step in order. I wrote them so I can explain every choice and every line.

## Step 1: Setup and a first run (`simulations/lock_exchange.jl`, `simulations/run_step1.jl`, `notebooks/01_setup_and_first_run.ipynb`)

### The question

Does the whole chain work: Julia runs an Oceananigans simulation, writes a NetCDF file outside the repository, and Python reads it and makes a plot? And which model settings do I need so that any mixing I measure later comes only from the numerics?

### Julia basics I needed

- Julia was installed with juliaup (version 1.13.1). The current Oceananigans needs Julia 1.10 or newer, so this is fine.
- A project environment is a folder with two files. `Project.toml` lists the packages I asked for; `Manifest.toml` records the exact version of every package, including the ones they depend on. Committing both means anyone can recreate exactly my set of versions. I created it inside the project folder with:

      julia --project=. -e 'using Pkg; Pkg.add(["Oceananigans", "NCDatasets"])'

  `--project=.` tells Julia to use the environment in the current folder instead of the global one. Installing took 1 minute 44 seconds and put about 1 GB in `~/.julia`.
- Installed versions: Oceananigans 0.113.5 and NCDatasets 0.14.15. NCDatasets is needed for NetCDF output: Oceananigans' `NetCDFWriter` is a package extension that only loads after `using NCDatasets` (my first run failed with exactly that error message).
- Julia compiles code the first time it runs it. So the first run of a script is slow and later runs in the same session are fast. My tiny test took about 32 seconds in total, of which the 444 time steps themselves are a small part; repeating an identical run inside the same Julia session took well under a second.
- A script is run from the project folder with `julia --project=. simulations/run_step1.jl`.

### The test case, and how it compares with Ilicak et al. (2012)

I checked the setup against the paper itself (Section 3.1):

- 2D channel, 64 km long and 20 m deep, no rotation, closed walls (no flow, free slip). Same in my setup.
- Dense water on the left half, light water on the right half; the paper uses temperatures of 5 and 35 °C with a linear equation of state and a density difference of 5 kg/m^3. I use buoyancy directly: b = -g (rho - rho0)/rho0 with rho0 = 1000 kg/m^3, so the dense water has b = 0 and the light water b = 9.81 x 5 / 1000 = 0.04905 m/s^2. This is the same density contrast.
- Zero tracer diffusivity, vertical viscosity 1e-4 m^2/s, base lateral viscosity 1e-2 m^2/s (larger values in their viscosity tests). Same in my setup.
- Grids: the paper uses dx = 500 m to 4 km and dz = 1 to 10 m, varied separately. I will vary dx and dz together (step 4).
- Time step: the paper uses a fixed 1 s step in all models. I use an adaptive step (step 6 tests fixed steps).
- Their models used free surfaces; MITgcm and MOM used a z* vertical coordinate. My choice is explained next.

### Model choices, and two problems I found on the way

- Hydrostatic model. Ilicak's models, MITgcm and MOM are all hydrostatic, and the lock exchange here is 64 km long but only 20 m deep, so I used `HydrostaticFreeSurfaceModel`.
- Fixed vertical grid instead of z*. Oceananigans' own lock-exchange example uses a z* grid (a vertical grid that stretches with the free surface). I tested both with an extra passive tracer that starts at exactly 1 everywhere: a consistent advection scheme must keep it at exactly 1. On the fixed grid it stayed at exactly 1; on the z* grid (with the split-explicit free surface) it drifted by 3e-4 within 3 hours. Such a drift would show up as false overshoots, so I used the fixed grid.
- A stiff free surface. With the real gravity in the free surface, total buoyancy was not conserved: it changed by 0.14% in 2 hours. The cause: on a fixed grid with a linear free surface, water moving through the top as the surface goes up and down carries buoyancy in and out. Adding the buoyancy in that surface layer (b at the top times the surface height) accounted for 95% of the change. That error would be about 3e-6 of the total mass, comparable to the RPE changes I want to measure (about 1e-5). So the free surface uses a gravity 100 000 times larger than the real one (`ImplicitFreeSurface(gravitational_acceleration = 9.81e5)`). The surface then moves by at most 3e-7 m, the model behaves like it has a rigid lid (the classic setup for this test), and total buoyancy changes by only 1.3e-8. The implicit free-surface solver stays stable with any gravity. The largest speed after 2 hours changed from 0.653 to 0.632 m/s, because the real free surface adds some fast surface sloshing.
- The time-step wizard. Oceananigans' `TimeStepWizard` adjusts the time step to keep the advective CFL number near a target. It can also respect a viscous limit, but in this version it computes that limit from the smallest horizontal spacing including the y-direction, which is `Flat` here and counted as 1 m. That capped the time step at 0.2 x (1 m)^2 / (0.01 m^2/s) = 20 s, even though the real cells are 2 km wide. I switched the wizard's viscous check off and gave it the real viscous limit, 0.2 x min(dx^2/nu_h, dz^2/nu_z), as its maximum time step.
- Output precision. `NetCDFWriter` saves 32-bit numbers by default (about 7 digits). The RPE changes are about 1e-5 of the RPE, so I save 64-bit numbers.

### The code, block by block (`lock_exchange.jl`)

- Constants: gravity, rho0 = 1000 kg/m^3, the density difference 5 kg/m^3, the buoyancy difference DELTA_B = 0.04905 m/s^2, the channel size, and the free-surface gravity.
- `run_lock_exchange(...)` takes the run name, output folder, grid spacing, the tracer and momentum advection schemes, the lateral and vertical viscosity, the stop time, the CFL target or a fixed time step, the snapshot interval and how often to check the run.
- Grid: `RectilinearGrid` with `topology = (Bounded, Flat, Bounded)`: walls at both ends and top and bottom, nothing in y.
- Closure: `HorizontalScalarDiffusivity(ν = νh, κ = 0)` and `VerticalScalarDiffusivity(ν = νz, κ = 0)`: viscosity for momentum, but zero diffusivity for buoyancy. Why zero: with no explicit diffusion, the only thing that can mix the two water masses is the error of the advection scheme, so any mixing I measure is numerical.
- Initial state: b = 0 for x < 32 km, b = DELTA_B elsewhere, water at rest.
- Output: `NetCDFWriter` saves b, u and w every 10 minutes (64-bit) to `fields.nc` in the run folder.
- A health check every 10 iterations: it writes the time, time step, smallest and largest b, total buoyancy, largest speeds and advective CFL number to `log.csv`, and stops the run cleanly (and records why) if anything is NaN, if b leaves [-DELTA_B, 2 DELTA_B], or if a speed exceeds 10 m/s. Checking often, not just at the end, is what my MITgcm run in project 15 was missing.
- At the end it writes `run_info.txt` with the settings, versions, status, number of iterations and wall-clock time.

### The first run

`run_step1.jl`: 32 x 5 cells (dx = 2 km, dz = 4 m), WENO (5th order) for buoyancy and momentum, 2 model hours, snapshots every 10 minutes. It completed in 444 iterations. The output folder holds `fields.nc` (105 kB), `log.csv` and `run_info.txt`.

The notebook reads the file with xarray. The dimension names show the staggered grid: b at cell centres (`x_caa`, `z_aac`), u on the x-faces (`x_faa`), w on the z-faces (`z_aaf`). There are 13 snapshots, stored as 64-bit numbers. The checks: total buoyancy changed by at most 2.09e-8 of its value; b stayed between -7.1e-7 and 0.04905071 m/s^2, so WENO created under- and overshoots of only 1.4e-5 of the buoyancy difference; the largest speed was 0.643 m/s.

### What to look for in the figure

- `figures/01_first_run_buoyancy.png`: at t = 0 the two water masses meet at x = 32 km. After 2 hours the dense water has moved about 4 km to the right along the bottom and the light water about 4 km to the left along the top. At 0.5 m/s, 2 hours gives 3.6 km. The interface between them is spread over a few 2 km cells.

### Three extra checks

- The free surface is implicit (`ImplicitFreeSurface`, recorded in every `run_info.txt`). An explicit free surface would need time steps shorter than dx / sqrt(g H) for its surface waves; with the stiff gravity and dx = 500 m that is 0.113 s. The baseline run (step 2) used steps of up to 43.9 s, so the stiff surface does not limit the time step.
- The y-direction is `Flat`: the grid's topology is `(Bounded, Flat, Bounded)`. Asking Oceananigans for the smallest y spacing still returns 1.0 m, which is why the time-step wizard's viscous limit came out as (1 m)^2 / nu. I saved a minimal example that reproduces this in `simulations/mwe_flat_y_diffusion_timescale.jl`, with its printed output next to it: the diffusion timescale is 100 s instead of dx^2 / nu = 4e8 s, and a 1000 s time step is cut to 20 s. I have not reported it to the Oceananigans developers.
- The uniform-tracer test is now part of the model function (`consistency_tracer = true` adds a passive tracer c = 1 and records its largest deviation from 1). In the baseline run the deviation was 2.2e-16, round-off.

### What I learned

The chain works, but getting a clean test needed three settings that are not the defaults: a fixed vertical grid, a stiff free surface, and a corrected time-step limit. Each one was found by checking a number that should have been constant or obvious (a uniform tracer, the total buoyancy, the time step), which is the habit this whole project depends on.

## Step 2: The baseline lock exchange (`simulations/run_step2.jl`, `notebooks/02_lock_exchange_baseline.ipynb`)

### The question

How much numerical mixing does Oceananigans' high-order WENO scheme create in the standard lock exchange, and does the flow itself behave as theory says?

### What I ran

Ilicak's base case: dx = 500 m, dz = 1 m (128 x 20 cells), lateral viscosity 1e-2 m^2/s, vertical viscosity 1e-4 m^2/s, zero tracer diffusivity, WENO (5th order) for buoyancy and momentum, 17 model hours, adaptive time step with advective CFL target 0.2, snapshots every 10 minutes, and the uniform check tracer. It completed in 1907 iterations; the whole script took 33 s, most of it compiling. The output file is 6.5 MB.

### How the analysis works (`notebooks/mixing_tools.py`)

- Front speed. Benjamin (1968) showed that a gravity current in a channel of depth H moves at u_f = 0.5 sqrt(g' H), where g' = g d(rho)/rho0 is the reduced gravity, here the buoyancy difference 0.04905 m/s^2. In plain words: the dense water's extra weight pushes it under the light water, and for an energy-conserving current filling half the depth, the speed comes out as half of sqrt(g' H), the speed of a long wave on the interface. Here u_f = 0.4952 m/s. I find each nose on every snapshot (the right-most bottom cell with b below half the difference, and the left-most top cell with b above it) and fit a straight line to its position between 3 and 15 hours.
- Reference potential energy (Winters et al. 1995). Imagine taking all the water and rearranging it, without mixing, into the most stable state possible: the densest water at the bottom, the lightest at the top. The potential energy of that sorted state is the RPE. Moving water around (gravity currents, waves) cannot change it, because the same set of densities sorts into the same state. Only mixing two densities together can raise it. In the code every cell has the same area, so I sort all the cell densities from heaviest to lightest and fill the grid row by row from the bottom, then add up g x density x height x area. Density is rho = rho0 (1 - b/g). Height is measured from the bottom. Like Ilicak et al. I report (RPE - RPE0)/RPE0; this number depends on where height zero is (they note this in their Appendix A), so I also report the rate dRPE/dt per square metre of channel floor, which does not.
- Buoyancy variance: the spread of b around its mean. Advection only moves water, so without diffusion the variance would stay constant; every loss is mixing.

### Results

- The flow: the dense water slides to the right along the bottom and the light water to the left along the top. At 17 h the noses are at 60.75 km and 3.25 km, just short of the walls. The two currents are exact mirror images (in this Boussinesq setup light water rising behaves like dense water sinking).
- Front speed: 0.4793 m/s for both currents, 3.2% slower than Benjamin's 0.4952 m/s. Ilicak et al. also found the models slightly slower than theory, with more mixing making them slower.
- RPE: up by 3.543e-5 of its starting value at 17 h, rising faster and faster (curving upward). dRPE/dt at 17 h is 1.890e-3 W/m^2.
- Buoyancy variance: down 18.56% at 17 h, falling at a nearly steady rate after the first half hour.
- The snapshots show where the mixing is: the interface between the two layers widens into a band of intermediate buoyancy about 10 m thick by 17 h. With zero diffusivity, all of that comes from the numerics.

### Sanity checks

- Total buoyancy: largest relative change 1.83e-8.
- RPE never decreases: the smallest change between snapshots is +7.62e-8 of RPE0.
- Uniform tracer: largest deviation from 1 is 2.2e-16.
- Overshoots: b went from -1.035e-5 to 0.04906035 m/s^2, so 2.1e-4 of the buoyancy difference beyond each end of the starting range. WENO is not strictly bounded, but nearly.
- Time step: up to 43.9 s; largest advective CFL number 0.257 (my log uses an upper bound that adds the largest u and w terms, even if they are in different cells).

### What to look for in the figures

- `figures/02_baseline_snapshots.png`: the two gravity currents and the growing band of intermediate (white) buoyancy between them.
- `figures/02_baseline_front_position.png`: both fronts move in steps of one 500 m cell, slightly slower than the dashed theory line.
- `figures/02_baseline_rpe_variance.png`: the RPE curving upward, the variance falling almost in a straight line.

### What I learned

Even a high-order scheme with tiny overshoots mixes the two water masses noticeably at this resolution and low viscosity: the interface becomes a 10 m thick band, the buoyancy variance drops by almost a fifth in 17 hours, and the RPE rises steadily. The physics still looks right (front speed within 3.2% of theory), so a model can get the flow right and the water masses wrong at the same time.

## Step 3: Tracer advection schemes (`simulations/run_step3.jl`, `notebooks/03_advection_schemes.ipynb`)

### The question

How much does the tracer advection scheme change the numerical mixing, and does any scheme create buoyancy values outside the starting range?

### What I changed

Only the scheme used for buoyancy: 2nd-order centred (`Centered(order=2)`), 3rd-order upwind-biased (`UpwindBiased(order=3)`), and 5th-order WENO (the step 2 baseline, not repeated). Momentum advection (WENO), viscosities, grid and run length are the same as in step 2.

What the schemes do, in plain words. To move b across a cell face, each scheme estimates b at the face from the cells around it.
- Centred: the average of the two neighbours. It has no built-in damping, so at a sharp jump it produces wiggles: values above the maximum and below the minimum.
- Upwind-biased: uses more cells from the side the flow comes from. This adds some damping, which reduces the wiggles but does not remove them.
- WENO (weighted essentially non-oscillatory): builds several estimates and weights them by how smooth the field is around each one, so near a jump it leans on the smooth side. This almost removes the wiggles, at the price of acting like a lower-order, more diffusive scheme right at the jump.

### A finding about my stopping rule

My health check stops a run if b goes more than one buoyancy difference outside its starting range. The centred run hit this after 59 minutes (b from -0.0494 to 0.0985). That limit was my own choice, so I made it a parameter (`b_range_limit`, default 1) and repeated the centred run with a limit of 10 (NaN still stops a run). It then ran all 17 hours without NaN, with b reaching -0.134 and 0.183 m/s^2. So the centred scheme did not blow up here: it carried huge, persistent overshoots. All results below use that run.

### Results at 17 hours

| scheme | RPE change | dRPE/dt (W/m^2) | variance change | overshoot each side | front speed |
| --- | ---: | ---: | ---: | ---: | ---: |
| centred, 2nd order | -9.485e-5 | -4.132e-3 | +0.75% | 2.74 x the difference | 0.5125 m/s (+3.5%) |
| upwind-biased, 3rd order | +1.577e-5 | +1.337e-3 | -15.24% | 0.379 x | 0.4881 m/s (-1.4%) |
| WENO, 5th order | +3.543e-5 | +1.890e-3 | -18.56% | 2.1e-4 x | 0.4793 m/s (-3.2%) |

Counting cells more than 1% of the density difference (0.05 kg/m^3) away from the two starting densities: the centred run has 42.6% of its cells mixed and 17.7% outside the starting range, the upwind run 31.3% mixed and 11.2% outside, and WENO 35.9% mixed and 0.0% outside. All three conserve total buoyancy to 1.8e-8, and the uniform check tracer stayed at 1 to within 3.3e-16 in all of them.

### Why the RPE can go down

The centred run's RPE falls and its variance rises, which mixing alone can never do, and my "RPE never decreases" check fails for it (smallest step between snapshots -4.0e-6 of RPE0) and, mildly, for the upwind run (-1.5e-6). The reason is the overshoots. The centred scheme creates water denser than the densest and lighter than the lightest starting water. In the sorted state that extra-dense water goes to the very bottom and the extra-light water to the top, which lowers the centre of mass: the opposite of mixing ("anti-diffusion"). It partly or completely hides the real mixing, which the density histogram shows is the largest of the three. Ilicak et al. (2012, Section 3.5) found the same with an unlimited Prather scheme: it "may appear to have less mixing from an integrated energetic point of view", but only through new extrema, which they call "generally unacceptable in realistic modelling". So the RPE is a fair measure of numerical mixing only for schemes that do not create new extremes; for the others it must be read together with the overshoots and the density distribution.

The upwind scheme shows a smaller RPE rise than WENO for both reasons: fewer mixed cells (31.3% against 35.9%), and some anti-diffusion from its 38% overshoots.

### Link to my MITgcm blow-up (project 15)

In project 15, one hosing run blew up in a single water column next to the 80°N wall, under an extremely fresh surface layer, with the tutorial's centred advection scheme, and stayed stable when rerun with a flux-limited scheme. This step shows the mechanism in a clean test: carrying a sharp front, the centred scheme puts values far outside the physical range (here up to 2.7 buoyancy differences beyond it) because nothing stops the wiggles at the jump. In a salinity field next to a very fresh layer, such wiggles mean unphysical salinities and densities, and those drive the flow. The lock exchange survived them; a run with stronger forcing that keeps sharpening a front may not.

### What to look for in the figures

- `figures/03_schemes_snapshots.png`: the centred run is speckled with blue and red cells inside the wrong layer; the upwind and WENO runs look similar, with a smooth band of mixed water.
- `figures/03_schemes_overshoots.png`: overshoots of order 1 (centred), 0.3 (upwind) and 1e-4 (WENO), reached within the first hour and then roughly steady. The "below 0" and "above the maximum" lines lie on top of each other because the two currents are mirror images.
- `figures/03_schemes_rpe_variance.png`: the centred run's RPE goes negative and its variance stays flat; the other two rise and fall steadily.
- `figures/03_schemes_density_pdf.png`: the centred and upwind distributions have tails outside 995-1000 kg/m^3; WENO's does not.

### What I learned

A scheme can look like it mixes less on a single number while it actually mixes more and also invents new water masses. WENO was the only scheme here that kept the water within its starting range, and it paid for that with a little more mixing than the upwind scheme. To judge a scheme I need the overshoots and the density distribution as well as the RPE.

## Step 4: Grid resolution (`simulations/run_step4.jl`, `notebooks/04_resolution.ipynb`)

### The question

Does the numerical mixing shrink when the grid gets finer, and how fast?

### What I changed

For the centred (2nd order) and WENO (5th order) tracer schemes I ran dx = 2000, 1000, 250 and 125 m, with dz = dx / 500 so every grid has the same cell shape (the 500 m runs come from steps 2 and 3). This differs from Ilicak et al. (2012), who varied dx and dz separately. The coarsest grid has only 5 vertical levels; the finest has 512 x 80 cells. Everything else is as in the baseline, including the lateral viscosity of 1e-2 m^2/s; the centred runs use the wide stopping limit from step 3.

Cost: the 125 m WENO run took 23 257 time steps and 123 s, and its output file is 102 MB (the 0.25 m layers make the time step short). All the step 4 runs together took about 4.5 minutes, and the raw output of the project is now 300 MB. The wall times in the tables include about 20 s of compiling for the first run of each script, so only the larger differences between them mean anything.

### Results at 17 hours

| dx (m) | grid | WENO RPE change | WENO dRPE/dt (W/m^2) | WENO mixed cells | centred overshoot | centred RPE change |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2000 | 32 x 5 | 5.126e-5 | 4.720e-3 | 61.3% | 0.81 | -2.014e-4 |
| 1000 | 64 x 10 | 4.392e-5 | 2.083e-3 | 40.9% | 1.71 | -1.219e-4 |
| 500 | 128 x 20 | 3.543e-5 | 1.890e-3 | 35.9% | 2.74 | -9.485e-5 |
| 250 | 256 x 40 | 3.324e-5 | 1.855e-3 | 34.0% | 3.63 | -9.080e-5 |
| 125 | 512 x 80 | 3.205e-5 | 1.827e-3 | 33.6% | 5.85 | -9.133e-5 |

(Overshoots are in units of the buoyancy difference; "mixed cells" are cells more than 0.05 kg/m^3 away from both starting densities.)

- WENO: the RPE change at 17 h falls only like dx^0.18 (dRPE/dt like dx^0.29). Making the grid 16 times finer reduced it by a factor of only 1.60, and between 250 and 125 m it hardly changed (local slope 0.05). The buoyancy variance loss stayed between 18.2% and 18.7% for every grid from 1000 m down. WENO stayed within its starting range on every grid (overshoots of 3.6e-5 to 3.2e-4 of the difference, growing slightly on finer grids), and the front speed stayed between 0.479 and 0.483 m/s.
- Centred: the RPE fell on every grid, so it cannot measure this scheme's mixing (step 3). Its overshoots grew on finer grids, from 0.81 to 5.85 buoyancy differences (roughly like dx^-0.68): a finer grid carries sharper fronts, and the centred scheme answers a sharper front with bigger wiggles. Between 17.7% and 23.8% of its cells were outside the starting range on every grid.

### Why the WENO mixing shrinks so slowly

In all these runs the lateral viscosity is 1e-2 m^2/s, so the grid Reynolds number U dx / nu (with U = 0.495 m/s) is between 6200 (125 m) and 99 000 (2000 m), always far above 100. Ilicak et al. (2012, Section 3.4 and Fig. 6) found that in this regime the spurious mixing is "high and steady (capped)", and that "changing horizontal grid spacing or viscosity does not reduce the spurious mixing in this regime", presumably because grid-scale noise in the velocity keeps stirring the interface. My WENO results fit that picture. The snapshots also show that the band of mixed water is almost as thick at 125 m as at 500 m, so the mixing is not simply confined to the coarse cells. From these runs alone I cannot tell how much of it is grid noise that does not go away and how much is set by the nearly inviscid flow itself; step 5 tests the viscosity directly.

### What to look for in the figures

- `figures/04_resolution_snapshots.png`: on the WENO side the interface gets smoother but not much thinner as the grid is refined; on the centred side the noise gets finer and denser, never smoother.
- `figures/04_resolution_convergence.png`: the WENO RPE change against dx on log-log axes, far flatter than a slope-1 line; and the overshoots, rising on finer grids for the centred scheme.
- `figures/04_resolution_rpe_time.png`: the WENO RPE curves for 500, 250 and 125 m lie close together; the 2000 and 1000 m curves rise faster.

### What I learned

A finer grid does not automatically remove numerical mixing. At this very high grid Reynolds number, a 16 times finer grid reduced the WENO mixing by only a factor of 1.6, and it made the centred scheme's overshoots worse. That points to the momentum settings (viscosity, momentum scheme) as the next thing to test, as Ilicak et al. argued.

## Step 5: Viscosity and the grid Reynolds number (`simulations/run_step5.jl`, `notebooks/05_viscosity_and_grid_reynolds_number.ipynb`)

### The question

Does the momentum side of the model, the viscosity and the momentum advection scheme, change how much the tracer is mixed? Ilicak et al. (2012) found that spurious mixing grows with the grid Reynolds number and is "capped" above about 100.

### The grid Reynolds number in plain words

Re = U dx / nu compares how fast the flow carries momentum across one grid cell (U dx) with how fast viscosity smooths it over that cell (nu). If Re is small, viscosity smooths every cell-to-cell wiggle in the velocity before it can grow. If Re is large, wiggles at the scale of the grid survive. A classic analysis of centred schemes (cited by Ilicak et al.) says Re must be below about 2 for the grid-scale noise to be damped. Noisy velocities stir the tracer at the grid scale, and that is where advection schemes make their errors, so more noise means more numerical mixing. I use Ilicak's definition, U = the theoretical front speed 0.495 m/s.

### What I changed

The tracer scheme stayed WENO (5th order) on the baseline grid. I varied the lateral viscosity: 0.01, 0.1, 1, 10, 100 and 200 m^2/s, giving Re from 24 761 down to 1.24.
- Set 1: centred (2nd order) momentum advection, which adds no damping of its own, so the explicit viscosity alone controls the velocity noise.
- Set 2: WENO momentum advection, which has built-in damping (the nu = 0.01 case is the step 2 baseline).

All 11 new runs completed (none of the low-viscosity centred runs blew up), in 42 s together. WENO kept the tracer within 3.2e-4 of its starting range in all of them, and the RPE never decreased, so the RPE is a fair measure here.

To test the explanation directly, I also measured the grid-scale noise in the velocity: the root-mean-square of (u[i-1] - 2 u[i] + u[i+1]) / 4 at 17 h, which is large for a cell-to-cell zigzag and small for a smooth field.

### Results

| nu_h (m^2/s) | Re | dRPE/dt, centred momentum (W/m^2) | dRPE/dt, WENO momentum (W/m^2) | velocity noise, centred / WENO (m/s) |
| ---: | ---: | ---: | ---: | ---: |
| 0.01 | 24761 | 4.068e-3 | 1.890e-3 | 0.0146 / 0.0104 |
| 0.1 | 2476 | 4.056e-3 | 1.882e-3 | 0.0145 / 0.0104 |
| 1 | 248 | 3.938e-3 | 1.861e-3 | 0.0137 / 0.0104 |
| 10 | 24.8 | 3.317e-3 | 1.673e-3 | 0.0122 / 0.0100 |
| 100 | 2.48 | 1.093e-3 | 7.070e-4 | 0.0055 / 0.0050 |
| 200 | 1.24 | 4.710e-4 | 4.381e-4 | 0.0029 / 0.0028 |

- My runs show what Ilicak et al. found. For Re above 100 the mixing is flat ("capped"): it changes by only 3% (centred momentum) and 2% (WENO momentum) over a factor of 100 in viscosity. Below Re of about 25 it falls quickly: at Re = 1.24 it is 8.6 times smaller than the capped value with centred momentum and 4.3 times smaller with WENO momentum. Ilicak et al. reported about an order of magnitude between Re above 100 and Re below 5.
- The momentum scheme matters as much as the viscosity. In the capped regime, WENO momentum gives 2.15 times less mixing than centred momentum, with the same tracer scheme. Its built-in damping acts like an extra viscosity at the grid scale: the velocity noise was 0.0104 m/s against 0.0146 m/s. At high viscosity the two sets meet (ratio 1.08 at nu = 200), because the explicit viscosity then does the damping. This matches Ilicak's explanation of why ROMS, with an upwind-biased momentum scheme, mixed less than MITgcm and MOM with centred momentum.
- The velocity noise follows the mixing in both sets: flat at high Re, falling at low Re. That supports the idea that grid-scale velocity noise drives the spurious mixing.
- More mixing also slows the currents: with centred momentum the front speed was 0.452 m/s at Re = 24 761 (8.7% below theory) and 0.485 m/s at Re = 1.24.
- The grid Reynolds number with the measured largest speed instead of the theoretical front speed is 1.4 to 1.8 times larger; it shifts the curves sideways but does not change the picture.
- The cost: a viscosity of 100 to 200 m^2/s is far larger than any real ocean viscosity at 500 m, and Ilicak et al. warn that it also smooths the real shear. Low numerical mixing bought with high viscosity changes the physics.

### What to look for in the figures

- `figures/05_mixing_vs_grid_reynolds.png`: left, dRPE/dt against Re on log-log axes, flat to the right of Re = 100 and falling to the left, with WENO momentum lower throughout; right, the velocity noise with the same shape.
- `figures/05_viscosity_snapshots.png`: at nu = 0.01 the interface is broad and disturbed, with bumps near 13 and 50 km; at nu = 100 and 200 it is thin and smooth.
- `figures/05_front_speed_vs_grid_reynolds.png`: the centred-momentum fronts are slowest at high Re and approach theory as Re falls.

### What I learned

The tracer scheme is only half of the story. With the same tracer scheme, the spurious mixing changed by a factor of 9 just through the viscosity, and by a factor of 2 through the momentum scheme, because what the tracer scheme has to handle is set by how noisy the velocity is at the grid scale. This explains step 4: refining the grid did little because the grid Reynolds number stayed far above 100.

## Step 6: Time step and stability (`simulations/run_step6.jl`, `notebooks/06_time_step_and_stability.ipynb`)

### The question

How large can the time step be before the results change, and where does the run become unstable?

### The CFL number

CFL = u dt / dx + w dt / dz: how many grid cells the water moves in one time step. An explicit scheme builds each new value from nearby old values; if the water moves further than that in one step, errors are amplified every step and the run blows up. The safe limit depends on the time-stepping method (here the default `QuasiAdamsBashforth2`) and the advection scheme. In this test the vertical part dominates: at the fastest moment (about 36 minutes in) w/dz was 0.0058 per second and u/dx only 0.0015.

### What I ran

The baseline setup (WENO tracer and momentum, dx = 500 m, dz = 1 m, nu_h = 1e-2 m^2/s) with fixed time steps dt = CFL / 0.0073083 for CFL 0.2, 0.5, 1.0 and 1.5. The 1.5 run still finished, so I added 2.0, 2.5, 3.0 and 4.0. Every run was checked at every step. CFL 2.0 and 2.5 were repeated with only the NaN and speed checks, to separate "badly wrong" from "blowing up". I also ran the adaptive TimeStepWizard with a target of 0.5 (target 0.2 is the step 2 baseline).

### Results

| fixed CFL target | dt (s) | largest measured CFL | outcome | overshoot | RPE change at 17 h |
| ---: | ---: | ---: | --- | ---: | ---: |
| 0.2 | 27.4 | 0.201 | completed | 2.1e-4 | 3.55e-5 |
| 0.5 | 68.4 | 0.505 | completed | 2.1e-4 | 3.53e-5 |
| 1.0 | 136.8 | 1.038 | completed | 2.3e-4 | 3.65e-5 |
| 1.5 | 205.2 | 1.666 | completed, degraded | 7.5e-2 | 3.41e-5 |
| 2.0 | 273.7 | 2.37 | stopped at 9.8 h (out of range); without the range check it ran 17 h with overshoots of 3.27 | | |
| 2.5 | 342.1 | 3.52 | stopped at 1.4 h; without the range check the largest speed grew to 38.7 m/s by 2.4 h: a real blow-up | | |
| 3.0 | 410.5 | 4.57 | stopped at 1.4 h | | |
| 4.0 | 547.3 | 7.35 | stopped at 1.0 h | | |

- Up to CFL 1 the results hardly change (RPE change 3.53e-5 to 3.65e-5, front speed 0.4793 m/s in every run).
- At CFL 1.5 WENO stops being clean: overshoots grow from 2e-4 to 7.5e-2 of the buoyancy difference.
- At CFL 2.0 the run survives but is badly wrong. From 2.5 on it is unstable: the speed grows from 4.7 to 7.2 to 38.7 m/s in three steps.
- Adaptive against fixed: with target 0.2 the wizard needed 1907 steps against 2244 for the fixed step, with the same result (RPE change 3.543e-5 against 3.553e-5), because it takes longer steps once the flow slows after the first hour. With target 0.5 it needed 1024 steps against 918, and its largest measured CFL was 0.769, above the target, because it only adjusts every 10 iterations and limits each change to 10%.

### Link to project 15

There the tracer time step was 1 day; in the column that blew up, the local vertical CFL rose from 0.03 to 0.52 and the centred tracer scheme created new extremes. Here WENO stayed clean up to about CFL 1, but step 3 showed the centred scheme overshooting badly even at CFL 0.2. The vertical velocity set the limit in both cases.

### What I learned

For this scheme the time step barely matters up to CFL 1, degrades the water masses before it breaks the run, and causes a runaway from about 2.5. Checking every step caught each failure within one step of it happening.

## Step 7: Summary (`notebooks/07_summary.ipynb`)

### The question

Which numerical choices create the most spurious mixing, taken together?

### What I did

No new runs. The notebook reads all 34 run folders from steps 2 to 6 (the step 1 test run is left out) and computes the same measures for each, so every run is treated the same way. The table is in `data/summary_all_runs.csv`, with one row per run: tracer and momentum scheme, dx, dz, viscosity, grid Reynolds number, largest measured CFL number, RPE change and dRPE/dt at 17 h, variance change, overshoot, front speed, buoyancy drift, the uniform-tracer check, iterations, wall time and status. 28 runs completed and 6 were stopped by the health check (their measures at 17 h are left empty).

One detail needed care. The table flags whether the RPE is a fair measure for each run (`rpe_valid`). My first rule, "RPE never falls between two snapshots", also flagged the 2000 m WENO run and the CFL 1.5 run, but their dips were about 1e-8 of the RPE, the size of rounding, against dips of about 4e-6 in the centred runs. So the rule became: RPE never falls by more than 1% of its total change over the run. Now only the runs with real new extremes fail it: all centred and upwind-biased tracer runs, and the CFL 2.0 run.

The summary figure compares the total RPE change over 17 h rather than dRPE/dt at 17 h, because the rate at one moment is noisy on the coarsest grid (dRPE/dt there was 2.5 times the baseline, the total change 1.45 times).

### Results

| choice, changed from the baseline | RPE change / baseline | largest overshoot |
| --- | ---: | ---: |
| baseline (WENO, 500 m, nu 0.01) | 1.00 | 2.1e-4 |
| tracer: upwind-biased 3rd order | not valid | 0.38 |
| tracer: centred 2nd order | not valid | 2.7 |
| grid: dx = 2000 m (5 levels) | 1.45 | 3.6e-5 |
| grid: dx = 125 m | 0.90 | 3.2e-4 |
| momentum: centred, nu 0.01 | 1.99 | 3.1e-4 |
| viscosity: nu = 200, WENO momentum | 0.28 | 2.2e-4 |
| viscosity: nu = 200, centred momentum | 0.30 | 2.2e-4 |
| time step: fixed CFL 1.0 | 1.03 | 2.3e-4 |
| time step: fixed CFL 1.5 | 0.96 | 7.5e-2 |

- The momentum settings changed the amount of mixing the most: a factor of 6.5 in the total RPE change (8.6 in the rate at 17 h) between nu = 0.01 and 200 with centred momentum, and a factor of 2 between centred and WENO momentum at low viscosity.
- The tracer scheme decided whether new water masses appear. That is a different kind of error from mixing, and a single mixing number hides it.
- Resolution helped little while the grid Reynolds number stayed far above 100 (a factor of 1.6 over 16 times refinement).
- The time step did not matter up to CFL 1, and then the water masses degraded before the run failed.
- Checks over all runs: total buoyancy changed by at most 2.2e-8, and the uniform tracer stayed at 1 to within 3e-15 in every run up to CFL 1.0.

### What to look for in the figure

- `figures/07_summary.png`: on the left, only the two viscosity bars are far below 1 and only the centred-momentum bar is far above it; the resolution and time-step bars stay near 1 except the coarse 2000 m grid. On the right (log scale), the upwind, centred and CFL 1.5 bars are about 350 to 13 000 times the baseline overshoot.

### Why it matters for climate models

A climate model runs for centuries, and spurious mixing acts at every time step. Mixing across the boundary between two water masses lifts dense water and lowers light water: the stratification weakens and the potential energy rises, which is exactly what the RPE measures. Over centuries this can erode the properties of deep water masses and change the density differences that drive the overturning circulation. A model's real mixing parameters are tuned on top of its numerical mixing, so if the numerical part is large, the tuning hides it. My test only shows the mechanism, in 2D and over 17 hours; it does not say how large the effect is in a climate model.

### What I learned

In this test, what decided the amount of mixing was mostly not the tracer scheme but the velocity field it had to carry: how much grid-scale noise the viscosity and the momentum scheme allowed. The tracer scheme decided something else, whether the model creates water that should not exist. Both need checking, and they need different measures.

## What I would do next

### Spurious mixing in my ice-shelf cavity runs (project 16)

My project 16 MITgcm runs used the default tracer advection (centred second order, since no scheme was set in `data`) with a horizontal diffusivity of 100 m^2/s, and the ISOMIP+ run used a flux-limited scheme (33) with 1 m^2/s. A cavity is a hard case for numerical mixing: the meltwater plume is a thin, buoyant layer flowing along a sloping ice base, cut across by the steps of a z-level grid, which is the situation Griffies et al. (2000) described. I would:

- add a passive tracer with no sources, and a tracer that starts at exactly 1 everywhere as a consistency check, as in this project;
- compute a tracer variance budget for the cavity: without any mixing the volume integral of the squared tracer only changes through the open boundary, explicit diffusion removes a known amount, and what is left over is the numerical mixing;
- compare the centred and the flux-limited schemes, and two resolutions, and see how much the melt rate and the meltwater plume change when the numerical part changes.

The RPE method from this project does not carry over directly, because the cavity has restoring, melting and an open boundary, so the RPE also changes for real reasons.

### Numerical mixing as its own term in water-mass budgets

An ocean model's water-mass budget in density (or temperature) classes says how much water is turned from one class into another, by surface fluxes and by mixing. The surface fluxes and the explicit (parameterized) mixing can be computed from the model's own output. The numerical mixing usually cannot, because it is hidden inside the advection scheme. It is often found as the residual: the change of volume in each density class, minus what the surface fluxes and the explicit mixing explain. I would like to learn how to compute this term in a realistic model and how its size compares with the explicit mixing in the deep ocean, where the physical mixing is small and the numerical part could matter most.
