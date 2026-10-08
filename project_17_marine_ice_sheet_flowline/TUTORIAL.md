# Tutorial: a 1D marine ice-sheet flowline model

These notes go through each notebook in order. I wrote them so I can explain every step and every line of the model.

The setup is MISMIP experiment 3 (Pattyn et al. 2012). I checked every parameter against the paper's Tables 4 and 5 and its Eq. 16 (the bed):
- ice density 900 kg/m^3, seawater 1000 kg/m^3, g = 9.8 m/s^2;
- Glen's exponent n = 3, sliding exponent m = 1/3, sliding coefficient C = 7.624e6 Pa m^(-1/3) s^(1/3);
- accumulation a = 0.3 m/yr, ice divide at x = 0;
- A from 3e-25 down to 2.5e-26 Pa^-3 s^-1 and back up, in 13 steps.

I also checked them against Schoof (2007). His Table 1 has the same densities, g, n, m, C and a, his Eq. 10 is the same bed, and his Eq. 16 is the grounding-line flux formula I use. He notes that with this C a basal stress of 80 kPa gives a sliding speed of about 35 m/yr; my conversion gives 36.46 m/yr. His own values of A (his Table 2: 4.227e-25, 2.478e-25, 1.370e-25, 7.433e-26 and 3.935e-26 Pa^-3 s^-1, chosen from ice temperatures of -12 to -32 °C) differ from the MISMIP list; I use the MISMIP list.

Neither paper gives the length of a year in seconds or a calving-front position for experiment 3, so I chose them: 31 556 926 s (365.2422 days), and the front at 1600 km (step 1).

Units: the code works in metres and years. Velocities are in m/yr, A is converted to Pa^-3 yr^-1 (multiply by 31 556 926), and C to Pa m^(-1/3) yr^(1/3) (divide by 31 556 926^(1/3)), which gives 2.4126e4.

The model functions are in `notebooks/flowline.py`, next to the notebooks, so each notebook can simply write `from flowline import ...`.

## Step 1: Schoof's steady states (`notebooks/01_schoof_steady_states.ipynb`)

### The question

Can the ice sheet on this bed have more than one steady grounding-line position for the same ice softness A, and are some of those positions unstable? I answer it without a time-dependent model, using Schoof's (2007) formula for the ice flux across the grounding line.

### The equations in plain words

- The bed: b(x) = 729 - 2184.8 (x/750 km)^2 + 1031.72 (x/750 km)^4 - 151.72 (x/750 km)^6 metres. It is above sea level near the divide and crosses sea level at 478.7 km. Between 973.7 km (b = -748.9 m) and 1265.7 km (b = -629.7 m) it gets deeper going inland, by 119.2 m in total. That is the overdeepening.
- Flotation thickness: H_f = -(rho_w/rho_i) b. A column of ice floats when it weighs as much as the seawater it pushes aside. At the grounding line the ice is exactly this thick.
- Schoof's flux: q(x_g) = K h(x_g)^((m+n+3)/(m+1)), with K = [A (rho_i g)^(n+1) (1 - rho_i/rho_w)^n / (4^n C)]^(1/(m+1)). For n = 3 and m = 1/3 the exponent is 4.75. Where it comes from: the floating shelf is pulled apart by the difference between the ice pressure and the water pressure, and that pull grows with the thickness h. Upstream, in a short "boundary layer", the ice has to slow from shelf-like stretching to sliding against the bed. Matching the two gives the flux. Softer ice (larger A) flows faster, and more basal friction (larger C) holds it back.
- Steady state: everything that falls between the divide and the grounding line must leave across it, a x_g = q(x_g).
- Stability: if the grounding line advances a little, the snowfall upstream grows by a per metre and the outflow by dq/dx_g per metre. If dq/dx_g > a, the sheet loses mass and the grounding line retreats back: stable. Otherwise it keeps moving: unstable.

Why a bed that deepens inland is dangerous: if the grounding line retreats there, the ice at the new grounding line is thicker (the water is deeper), so the outflow grows, more ice leaves, and the retreat speeds up. On a bed that deepens seaward, retreat brings thinner ice and less outflow, which slows the retreat.

### The code, block by block

1. Parameters. The constants live in `flowline.py`; the notebook prints them and sets its own parameters (the 13 MISMIP values of A, the example A = 1e-25, the search range and spacing).
2. A unit check. The flux in year units equals the flux in second units times the seconds in a year, to a relative difference of 1.1e-16. This tests my conversion of A and C.
3. The bed and the flotation thickness. `brentq` (a root finder that needs a bracket where the function changes sign) finds where b = 0 and where the bed slope is zero.
4. The flux formula: `schoof_prefactor`, `schoof_exponent`, `schoof_flux`, and `schoof_flux_slope` (dq/dx_g = 4.75 q/h dh/dx_g, from the chain rule).
5. Steady states for one A. `steady_positions` evaluates q(x_g) - a x_g every 100 m, looks for sign changes, and refines each one with `brentq`. Each is then classified by comparing dq/dx_g with a.
6. The sweep of A. Rather than repeating the search for many A, I solve the steady-state condition for A: q is proportional to A^(1/(m+1)), so A(x_g) = (a x_g / h^4.75)^(m+1) x 4^n C / ((rho_i g)^(n+1) (1 - rho_i/rho_w)^n). This traces every steady state at once. At a steady state, dq/dx_g > a is the same as A(x_g) falling with x_g, so the folds are where A(x_g) has a minimum or maximum. I find them with `brentq` on the derivative of ln A(x_g), which is (m+1)/x_g - (m+n+3) h'/h.
7. The jumps. At each fold I search for the other stable steady state at the same A.
8. The MISMIP values of A, the saved data, and the calving front.

### Results

For A = 1e-25 Pa^-3 s^-1 there are three steady states:

| x_g (km) | bed (m) | dq/dx_g (m/yr) | |
| ---: | ---: | ---: | --- |
| 799.8 | -644.4 | 2.147 | stable |
| 1124.3 | -692.3 | -1.421 | unstable |
| 1376.3 | -722.4 | 5.285 | stable |

The two folds:

| fold | x_g (km) | A (Pa^-3 s^-1) | jump |
| --- | ---: | ---: | --- |
| lower (A falling) | 948.4 | 4.9296e-26 | advance to 1413.0 km (464.6 km) |
| upper (A rising) | 1275.0 | 2.1449e-25 | retreat to 741.3 km (533.7 km) |

For 4.930e-26 < A < 2.145e-25 there are two stable states. Four of the MISMIP values of A (2e-25, 1.5e-25, 1e-25 and 5e-26) are inside this range. 3e-25, 2.5e-25 and 2.5e-26 have one stable state each (721.9, 732.1 and 1440.7 km).

The folds are a little outside the overdeepening (948 and 1275 km, against 974 and 1266 km). Just outside it, the bed still deepens seaward but so gently that the outflow grows by less than a per metre, which is still unstable.

The calving front: the farthest steady grounding line in the MISMIP range is 1440.7 km (A = 2.5e-26). I put the front at 1600 km, 159 km beyond it, where the bed is at -2147 m. In 1D the shelf length does not affect the grounding line, so I only need some room.

### What to look for in the figures

- `figures/01_bed_and_flotation.png`: the bed rising toward the ocean inside the grey band, and the flotation thickness dipping there.
- `figures/01_flux_crossings.png`: the straight snowfall line and the steep flux curve cross three times. At the filled circles the red curve crosses the blue line from below (outflow grows faster than snowfall); at the open circle it crosses from above.
- `figures/01_steady_states_vs_A.png`: an S-shaped curve. The upper and lower solid branches are stable, the dashed middle is unstable, and the red arrows show where the grounding line must jump at each fold. Going down in A it jumps forward at 4.93e-26; going back up it stays on the upper branch until 2.14e-25 and then jumps back. That difference is the hysteresis.

### What I learned

The steady states follow from one balance, snowfall upstream against flux out, exactly like the box-model equilibria in project 14. On a bed that deepens inland the flux falls as the grounding line advances, so the balance gives an unstable branch and two folds, and the ice sheet's position depends on its history over a wide range of A.

## Step 2: The velocity solver and its test (`notebooks/flowline.py`, `notebooks/02_velocity_solver_and_shelf_test.ipynb`)

### The question

Given the ice thickness, how fast does the ice move? And does my solver get the right answer where the exact answer is known?

### The equation in plain words

In fast-sliding ice and in floating ice shelves, the ice moves almost like a plug: the velocity hardly changes with depth. This is the shallow-shelf approximation (SSA). The forces on each column of ice must balance:

d/dx [ 4 eta H du/dx ]  -  beta u  =  rho_i g H ds/dx

- Stretching stress (first term): neighbouring ice pulls on the column when the ice stretches. eta is the viscosity, H the thickness, du/dx the stretching rate.
- Basal friction (second term): beta u = C |u|^(m-1) u, only where the ice is grounded.
- Driving stress (right-hand side): the weight of the ice pushes it down the surface slope.
- Glen's flow law: eta = 0.5 A^(-1/n) |du/dx|^((1-n)/n), so ice that deforms faster is softer. This makes the equation nonlinear.
- At the divide u = 0. At the calving front the stretching stress balances the push of the ice cliff minus the push of the seawater on its submerged part: 0.5 rho_i g H^2 - 0.5 rho_w g D^2, where D is the depth of the ice base below sea level. For floating ice this is 0.5 rho_i g (1 - rho_i/rho_w) H^2.

This is the form in MISMIP (Pattyn et al. 2012, Eq. A5), where 4 eta H du/dx is written 2 A^(-1/n) H |du/dx|^(1/n - 1) du/dx.

### The code, block by block

- `make_grid`: a staggered grid. Thickness at the N cell centres, velocity at the N + 1 cell edges; edge 0 is the divide and edge N the calving front. Why staggered: the ice flux u H across an edge is exactly what moves ice between cells (step 3); the surface slope between two centres is an edge value, where the driving stress acts; and du/dx between two edges is a centre value, where eta and H are. If everything sat at the same points I would need averages everywhere, and the solution could develop a zigzag that the equations do not see.
- `surface_elevation`: s = max(b + H, (1 - rho_i/rho_w) H). Grounded ice sits on the bed; floating ice has a fraction rho_i/rho_w of its thickness below sea level. Whichever is higher is the real surface.
- `effective_viscosity`: eta at the centres from du/dx between the two edges, with a floor: |du/dx| is replaced by sqrt((du/dx)^2 + 1e-6^2). Without it eta is infinite where the ice does not stretch (at the divide, and everywhere at the first iteration, when u = 0). 1e-6 per year is far below the real stretching rates (4.2e-3 per year in the test).
- `basal_drag_coefficient`: beta = C |u|^(m-1) at the edges, with a floor of 1e-3 m/yr on the speed for the same reason (m - 1 < 0, so it is infinite at u = 0), multiplied by the grounded fraction of each edge (0 under floating ice; step 3 will use values between 0 and 1).
- `calving_front_stress`: the cliff push at the front.
- `solve_ssa`: Picard iteration. With eta and beta taken from the current guess, each edge's equation involves only u there and at its two neighbours, so the system is tridiagonal. I store its three diagonals and solve it with `scipy.linalg.solve_banded`, which is fast and exact for this shape. Then I repeat with the new u until the largest change is below a tolerance times the largest speed.

### Results

- Picard iteration: starting from u = 0, each iteration reduces the change by a factor of 0.667 (about 1 - 1/n for n = 3), so it takes about 77 iterations to reach a tolerance of 1e-13. In the time-dependent model each solve can start from the previous velocity, which should need far fewer.
- Test 1, uniform shelf (H = 500 m, A = 1e-25 Pa^-3 s^-1, 200 km long). The exact solution is du/dx = A (rho_i g (1 - rho_i/rho_w) H / 4)^n = 4.2289e-3 per year, a front speed of 845.786 m/yr. The error was 5.592e-8 of the front speed at all six grid spacings (20 km to 0.625 km). The exact velocity is a straight line, which my centred differences reproduce exactly, so the only error left is the floor on the strain rate. It makes the ice stiffer by about (floor / strain rate)^2: with floors of 1e-4, 1e-5, 1e-6 and 1e-7 per year the errors were 5.585e-4, 5.592e-6, 5.592e-8 and 5.590e-10, matching (floor / strain rate)^2 = 5.592e-4, 5.592e-6, 5.592e-8 and 5.592e-10.
- Test 2, a shelf thinning linearly from 1000 m to 400 m. On a freely floating shelf the driving stress rho_i g H d/dx[(1 - rho_i/rho_w) H] is exactly the derivative of 0.5 rho_i g (1 - rho_i/rho_w) H^2, so the uniform-shelf formula holds at every point with the local thickness. Integrating with n = 3 gives u(x) = k (H(x)^4 - H1^4) / (4 H'), a front speed of 2747.11 m/yr. The error fell from 1.552e-3 (dx = 20 km) to 1.506e-6 (dx = 0.625 km), a factor of 4 for each halving: the convergence rate is 2.001 (between neighbouring spacings 2.000 to 2.007). Second order is what centred differences should give.
- A solve took 0.003 s with 10 cells and 0.021 s with 320 cells.

A note from Cuffey and Paterson (2010, Sect. 8.9.3): this formula is for a shelf that cannot spread sideways, which is the 1D flowline situation. For a shelf free to spread in both directions (Weertman 1957) the stretching is slower; for n = 3 the confined shelf spreads 9/8 times as fast, the "about 10%" they mention.

### What to look for in the figures

- `figures/02_picard_convergence.png`: the change per iteration falls along a straight line on the log axis: steady, linear convergence.
- `figures/02_shelf_velocity.png`: the solver's points sit on the exact curve, even at 20 km spacing.
- `figures/02_shelf_convergence.png`: the thinning-shelf errors lie on the slope-2 line; the uniform-shelf errors are flat at the regularization level.

### What I learned

A test needs an exact solution that the method cannot reproduce by accident: the uniform shelf only tested the regularization, while the thinning shelf showed the second-order convergence. The floors that keep the viscosity and friction finite have a measurable but tiny effect, and I can predict its size.

## Step 3: The time-dependent model and the grid-resolution problem (`notebooks/03_time_dependent_model_and_resolution.ipynb`)

### The question

When I run my model forward in time until it stops changing, does the grounding line end up where Schoof's theory says? And how does the answer depend on the grid spacing?

### The equations in plain words

- Mass conservation: dH/dt = a - d(uH)/dx. A column of ice gets thicker from snowfall and from ice flowing in, and thinner from ice flowing out.
- Grounded or floating: ice is grounded where it is thicker than the flotation thickness, H > -(rho_w/rho_i) b. Its surface is b + H when grounded and (1 - rho_i/rho_w) H when floating.
- At the calving front, the ice that flows across it leaves the domain (it calves).

### The code, block by block

- `edge_flux`: the flux u H at each edge, with H taken from the cell the ice comes from (upwind). Using the average of both neighbours instead would let small wiggles grow, because a cell's thickness would then depend partly on ice that has not reached it yet.
- `stable_time_step`: an explicit update is stable only if no ice moves more than one cell per step, dt < dx / max|u|. I use half of that (CFL = 0.5), capped at 10 years. The fastest ice, at the calving front (about 2000 m/yr in steady state), sets the step, so halving dx halves dt and doubles the number of cells: each halving costs about 4 times more.
- `run_to_steady_state`: at every step it decides which edges are grounded, solves the SSA for u (starting Picard from the previous u), takes a stable time step and updates H in every cell. Every 500 years it checks whether the grounding line moved less than 0.1 m/yr and the thickness changed less than 0.1 mm/yr on average; if both, the run is steady. It also keeps the ice volume, the total snowfall and the total calving for the mass check in step 5.
- `grounded_fraction_edges`, plain scheme: each edge is either grounded (full friction) or floating (none), from the average thickness above flotation of its two neighbouring cells.
- `grounded_fraction_edges`, sub-grid scheme (after Gladstone et al. 2010 and Leguy et al. 2014): the thickness above flotation, H - H_f, is positive on grounded and negative on floating ice. Between the last grounded and first floating centre, I place the grounding line where a straight line between the two values crosses zero, and give the edge between them friction times the grounded fraction of its stretch. Example from the notebook: with +25 m and -15 m above flotation at two centres 8 km apart, the grounding line is 25/40 = 0.625 of the way across, and that edge gets 0.625 of its friction.
- `grounding_line_position`: the same interpolation, used to report x_g for both schemes, so the reported position is not limited to grid points.
- Starting state, as in MISMIP: a 10 m layer of ice everywhere. It floats where the bed is below -9 m, so the starting grounding line is at 481.6 km.

### Results

I used A = 3e-25 Pa^-3 s^-1 (the first MISMIP value). Schoof's only steady state for it is at 721.90 km, where the bed is -530.2 m and deepens toward the ocean, well away from the overdeepening. Every run reached steady state after 18 500 to 22 500 years.

| scheme | dx (km) | x_g (km) | error (km) | time steps | run time (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| plain | 8 | 615.44 | -106.45 | 7943 | 3.1 |
| plain | 4 | 639.79 | -82.10 | 17180 | 7.0 |
| plain | 2 | 659.94 | -61.96 | 36042 | 16.8 |
| plain | 1 | 676.99 | -44.91 | 71529 | 45.2 |
| plain | 0.5 | 691.50 | -30.40 | 146230 | 136.9 |
| sub-grid | 8 | 691.84 | -30.06 | 8421 | 3.5 |
| sub-grid | 4 | 705.78 | -16.11 | 18723 | 8.0 |
| sub-grid | 2 | 710.98 | -10.92 | 36291 | 17.4 |
| sub-grid | 1 | 715.48 | -6.42 | 74438 | 43.7 |
| sub-grid | 0.5 | 717.73 | -4.16 | 152782 | 111.5 |

- Every run stops short of Schoof's position.
- With the plain scheme the error shrinks only like dx^0.45: even at 0.5 km the grounding line is 30.4 km short.
- The sub-grid scheme is much better: at 8 km it is as good as the plain scheme at 0.5 km, and its error shrinks like dx^0.70, to 4.16 km at 0.5 km. It is still not converged.
- In steady state the ice is up to 2948 m thick and moves at up to 2037 m/yr, at the calving front.

An extra check, advance against retreat. For the three coarsest grids I ran to steady state with A = 1e-25 first (Schoof: 799.8 km), then switched to A = 3e-25 so the grounding line had to retreat into place:

| scheme | dx (km) | start (km) | retreat ends (km) | advance ended (km) | difference (km) |
| --- | ---: | ---: | ---: | ---: | ---: |
| plain | 8 | 663.60 | 660.66 | 615.44 | +45.22 |
| plain | 4 | 691.95 | 690.16 | 639.79 | +50.37 |
| plain | 2 | 719.95 | 718.79 | 659.94 | +58.85 |
| sub-grid | 8 | 763.82 | 758.90 | 691.84 | +67.06 |
| sub-grid | 4 | 777.94 | 735.53 | 705.78 | +29.75 |
| sub-grid | 2 | 786.98 | 725.70 | 710.98 | +14.72 |

On the plain grid the retreat hardly happens: the grounding line moved only 1.2 to 2.9 km and the runs were steady after 1500 years. Notice also that the plain runs with A = 1e-25 stopped far short of 799.8 km. With the sub-grid scheme, advance and retreat end on opposite sides of Schoof's position (at 2 km: 710.98 and 725.70 km, around 721.90), and the gap shrinks with resolution (67.06, 29.75, 14.72 km).

### Why a plain fixed grid struggles

1. Friction switches off abruptly. On one side of an edge the ice has full friction, on the other none. In reality the change happens across Schoof's boundary layer, which is much shorter than a coarse grid cell.
2. The grounding line can only move a whole cell at a time. To advance one cell, the next cell must thicken until it grounds; until then the friction does not change. To retreat, a whole cell must thin to flotation. The grid can therefore hold the grounding line in many places, depending on where it came from. Pattyn et al. (2012) found the same in MISMIP: fixed-grid models were the most affected, and advance and retreat did not agree.

So on a fixed grid the computed steady state depends on the grid spacing and on the history, not only on the physics. That matters for step 4: some of the hysteresis a coarse model shows could come from the grid.

### What to look for in the figures

- `figures/03_resolution_error.png`: both lines fall with dx, the sub-grid one lower and steeper, but neither reaches slope 1.
- `figures/03_gl_vs_time.png`: all runs start together; the coarse plain runs advance in steps (one cell at a time) and stop early.
- `figures/03_steady_profiles.png`: the coarse plain profile is shorter and thinner; the close-up shows the floating shelf starting well before Schoof's grounding line.

### What I learned

A fixed grid can hold the grounding line in the wrong place. The error falls only slowly with resolution, and it depends on whether the grounding line advanced or retreated into position. Interpolating the grounding line and scaling the friction in that one cell cuts the error by a factor of 3.5 to 7.3 at the same grid spacing, but even at 0.5 km my model is about 4 km short of Schoof's position.

## Step 4: Hysteresis in the time-dependent model (`notebooks/04_hysteresis_experiment.ipynb`)

### The question

Does my time-dependent model show the hysteresis loop that Schoof's theory predicts in step 1? If A is lowered step by step and then raised again, does the grounding line jump across the overdeepening, and at different values of A going down and coming back?

### What I did

- The configuration: the sub-grid scheme at 1 km, the best I could run in reasonable time (step 3: 6.4 km short of Schoof's position at A = 3e-25).
- The sequence: the 13 values of A in MISMIP EXP 3a (Table 5): 3e-25 down to 2.5e-26, then back up to 3e-25 (Pa^-3 s^-1).
- Each value starts from the final ice sheet of the previous one (the first from a 10 m layer) and runs to steady state, with the same rule as step 3, checked every 100 years: the grounding line moves less than 0.1 m/yr and the thickness changes less than 0.1 mm/yr on average. Instead of MISMIP's fixed run lengths (15 000 or 30 000 years per step), each step runs until it is steady, up to 80 000 years.
- `run_to_steady_state` can now keep the thickness profile at every check (`keep_profiles=True`), so I can plot how a jump unfolds.
- A "jump" means the grounding line starts on one side of the overdeepening (973.7 to 1265.7 km) and ends on the other.

### Results

Every step reached steady state. The whole sequence covered 221 500 model years and took 5.5 minutes.

| step | A (Pa^-3 s^-1) | direction | x_g (km) | Schoof, same branch (km) | error (km) | years |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 3.0e-25 | start | 715.5 | 721.9 | -6.4 | 21 700 |
| 2 | 2.5e-25 | down | 725.5 | 732.1 | -6.6 | 7 500 |
| 3 | 2.0e-25 | down | 739.4 | 745.7 | -6.3 | 9 300 |
| 4 | 1.5e-25 | down | 758.5 | 765.5 | -7.0 | 9 500 |
| 5 | 1.0e-25 | down | 792.5 | 799.8 | -7.3 | 12 300 |
| 6 | 5.0e-26 | down | 908.5 | 926.1 | -17.6 | 31 600 |
| 7 | 2.5e-26 | down | 1433.5 | 1440.7 | -7.3 | 28 400 (jump) |
| 8 | 5.0e-26 | up | 1407.1 | 1412.4 | -5.3 | 13 400 |
| 9 | 1.0e-25 | up | 1372.1 | 1376.3 | -4.2 | 16 500 |
| 10 | 1.5e-25 | up | 1343.0 | 1346.1 | -3.1 | 17 000 |
| 11 | 2.0e-25 | up | 1307.9 | 1307.8 | +0.1 | 19 700 |
| 12 | 2.5e-25 | up | 730.9 | 732.1 | -1.2 | 26 700 (jump) |
| 13 | 3.0e-25 | up | 720.9 | 721.9 | -0.9 | 7 900 |

- The loop is there. Going down, the grounding line stays on the inner branch until A = 5e-26 and jumps forward across the overdeepening only at 2.5e-26. Going up, it stays on the outer branch through 2e-25 and jumps back only at 2.5e-25. For A = 5e-26, 1e-25, 1.5e-25 and 2e-25 the model has two different steady states depending on its history, as Schoof's theory predicts.
- The model's folds lie between the same MISMIP values as Schoof's: 4.93e-26 is between 5e-26 and 2.5e-26, and 2.145e-25 between 2e-25 and 2.5e-25.
- The errors are mostly -3 to -7 km, as in step 3, and a little smaller on the way up (-5.3 to +0.1 km) than on the way down (-6.3 to -7.3 km), the same advance/retreat difference as in step 3.
- The exception is step 6 (A = 5e-26): -17.6 km, and 31 600 years to settle. This steady state (926.1 km) is close to Schoof's fold (948.4 km at A = 4.93e-26). Near a fold the curve is almost vertical, so a small error in the flux moves the steady state a long way, and the ice sheet approaches it very slowly. This "critical slowing down" is the same effect behind the early warning signals before the Cessi fold in project 14.
- The two jumps:
  - Advance (step 7, A lowered to 2.5e-26): 908.5 to 1433.5 km. The grounding line crossed the overdeepening between year 1600 and 7400 (5800 years), on average 50 m/yr. Its fastest 100-year average, 80 m/yr, was right after A changed.
  - Retreat (step 12, A raised to 2.5e-25): 1307.9 to 730.9 km. It took 3800 years to reach the overdeepening, then crossed it between year 3800 and 11 300 (7500 years, 39 m/yr on average), and moved fastest, 167 m/yr, between year 11 100 and 11 200 near 982 km, at the landward end of the overdeepening.

### What to look for in the figures

- `figures/04_hysteresis.png`: the model's triangles sit on Schoof's stable branches; down-pointing triangles stay on the lower branch until the last step, up-pointing ones on the upper branch until A = 2.5e-25. No model state lies on the dashed unstable branch.
- `figures/04_jump_gl_vs_time.png`: S-shaped curves. The retreat is slow at first, then speeds up through the overdeepening and stops sharply on the normal bed beyond it: the instability feeds itself until the bed slope changes.
- `figures/04_jump_profiles.png`: during the advance the whole ice sheet thickens and the grounding line moves through the grey band without stopping.

### What I learned

The time-dependent model reproduces the hysteresis that the flux formula predicts: for the same A it can sit in two different places, and which one depends on its history. Crossing the overdeepening takes thousands of years, but once a retreat starts on the reverse slope it speeds up on its own. Near a fold the model is slow and least accurate, which matters because tipping points are exactly where we most want accurate answers.

## Step 5: Checks and summary (`notebooks/05_checks_and_summary.ipynb`)

### The question

Does my model conserve ice exactly? What do all the runs say together? And how does this project connect to projects 14 and 16?

### The mass check

The ice volume per metre of width is V = sum of H dx. Summing the update H_new = H + dt (a - (flux out - flux in)/dx) over all cells, every flux between two cells appears once with + and once with -, so only the two ends are left: the divide (flux 0) and the calving front. So the change in V must equal the total snowfall minus the total ice that crossed the front. `run_to_steady_state` adds up both as it goes.

For the run from a 10 m layer to steady state (sub-grid, 4 km, A = 3e-25), the volume grew by 1.826951e9 m^2 and snowfall minus calving was 1.826951e9 m^2: the largest difference over all checks is 1.87e-14 of the change, which is rounding error. In the last 500 years, the snowfall (480 000 m^2/yr) and the calving (479 841.55 m^2/yr) differed by 0.033%, so the steady state is a true balance. The runs with shelf melting (below) also pass, with the melted ice added to the budget (5.8e-15 and 8.0e-16).

### Summary of all runs

| what | result |
| --- | --- |
| SSA solver, thinning shelf (step 2) | error falls from 1.552e-3 (20 km) to 1.506e-6 (0.625 km), convergence rate 2.001 |
| plain fixed grid, A = 3e-25 (step 3) | 106.45 km short of Schoof at 8 km, 30.40 km at 0.5 km |
| sub-grid scheme, A = 3e-25 (step 3) | 30.06 km short at 8 km, 4.16 km at 0.5 km |
| hysteresis, sub-grid 1 km (step 4) | 13 steps, 221 500 model years, 5.5 minutes; error from -17.6 to +0.1 km, median 6.3 km |

The full table is in `data/summary_table.csv`, and the run times are in the step 3 table.

### The same method as project 14

In step 1 I found the steady states by plotting the snowfall upstream (a x_g) and the outflow (Schoof's q) against the grounding-line position and looking for crossings, and judged their stability by which curve rises faster. In project 14 I did the same for Stommel's box model. In both, the steady states against the control parameter (A here, the freshwater forcing there) form an S-shaped curve with two folds, and between the folds the history decides which of two stable states the system is in.

### Why, in 1D, the shelf cannot hold the grounding line back

In step 2 I showed that on a freely floating shelf the stretching stress at every point equals the push of an ice cliff of the local thickness. At the grounding line the thickness is fixed by the water depth, so the stress there does not depend on the shelf downstream. Two tests:

- Halving the thickness of the shelf (mean 257 to 128 m) and solving for the velocity again changed the speed of the grounded ice by at most 7.2e-10 m/yr (the grounded ice moves at up to 276 m/yr). Only the shelf itself slowed down: 2038 to 591 m/yr at the calving front, because thinner floating ice stretches more slowly.
- Melting the shelf from below at 2 and 10 m/yr thinned it from a mean of 258 m to 29.6 and 8.5 m, and moved the grounding line by -0.11 and 0.00 km.

A numerical lesson: my first version melted every floating cell, including the first one after the grounding line. With 10 m/yr the grounding line then retreated by 36.81 km in 5000 years and jumped back and forth. That came from the grid: on a 4 km grid the first floating cell's centre can be almost 4 km from the grounding line, so melting it removes ice that is still almost grounded. Not melting that cell, as recommended by Seroussi and Morlighem (2018) and Wang et al. (2024), removes the retreat.

### What this means for the real, 3D world (project 16)

Real ice shelves sit in embayments and rest on pinning points, so the side walls and bumps carry part of the shelf's push and the shelf holds the grounded ice back (buttressing). The stress at the grounding line is then smaller than the free-floating value. Schoof (2007, Eq. 29) writes this as a factor theta < 1 on the stress, which multiplies the flux by theta^(n/(m+1)) = theta^(9/4). If ocean melting thins a buttressing shelf, theta rises toward 1, the flux increases, and the grounding line can retreat, possibly onto a reverse slope where the instability of step 4 takes over. That is why the ocean melting I modelled in project 16 matters for the ice sheet, even though in my 1D model it has no effect.

### What to look for in the figures

- `figures/05_mass_check.png`: the volume change and snowfall minus calving lie on top of each other.
- `figures/05_shelf_melt_profiles.png`: the shelf almost disappears with melting, but the grounded ice and the grounding line do not change.
- `figures/05_melt_next_to_grounding_line.png`: melting the first floating cell makes the grounding line jump back and forth and retreat; not melting it keeps it still.

### What I learned

The model conserves ice to rounding error. In 1D the shelf has no influence on the grounding line, so any effect of ocean melting on the ice sheet needs the second horizontal dimension (buttressing). And where a model applies melt near the grounding line can create a retreat that is purely numerical.

## What I would do next

- Couple this model to an ocean box model, in which the ocean temperature in front of the ice shelf sets a melt rate near the grounding line. To make that melt matter in 1D, I would add Schoof's buttressing factor theta (optional step 6) and let it depend on the shelf thickness, so that melting reduces buttressing.
- Compare synchronous coupling (ocean and ice exchange information every time step) with asynchronous coupling (the ocean is updated only every few decades or centuries, which is cheaper). The question is whether infrequent coupling changes when the grounding line crosses a fold, which matters near a tipping point where the system responds slowly.
- Repeat the hysteresis experiment with slowly changing forcing instead of steps, to see how far the system lags behind its steady states near the folds.
