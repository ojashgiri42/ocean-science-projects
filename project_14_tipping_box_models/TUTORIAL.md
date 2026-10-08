# Tutorial: tipping points in simple box models

These notes walk through each notebook in order. I wrote them so I can explain every line of the project.

The notebooks are in `notebooks/`, and they all import the shared functions in `src/models.py`. Every function there takes its parameters as arguments, so each notebook can choose its own values. Each notebook is split into small cells, with a short note before each cell saying what it does.

## Step 1: Stommel's model and its equilibria (`notebooks/01_stommel_equilibria.ipynb`)

### The question

Can a simple two-box ocean have more than one steady circulation for the same forcing? If so, which starting states end up in which one?

### The equations in plain words

Stommel has two boxes of water, one at low latitude and one at high latitude, joined by a pipe. T is the temperature difference between them and S is the salinity difference. Both are dimensionless.

    dT/dt = 1 - T - |f| T
    dS/dt = delta (1 - S) - |f| S
    f = (-T + R S) / lambda

- `1 - T`: the atmosphere pulls the temperature difference back toward 1. This happens fast.
- `delta (1 - S)`: the salinity difference is pulled back toward 1. In Stommel's setup each box exchanges salt with a large reservoir through a porous wall, so this is a restoring condition, not a fixed freshwater flux. In ocean terms it stands for evaporation and rain. With delta = 1/6 it is 6 times slower than for temperature. That matters: Stommel showed that three equilibria need R delta < 1 (when R > 1), so salinity must be restored more weakly than temperature. Here R delta = 1/3.
- `|f| T` and `|f| S`: the flow mixes the two boxes and wears down both differences. Only the speed |f| matters, because mixing works the same in either direction.
- `f = (-T + R S) / lambda`: the flow is driven by the density difference. Warm water is light and salty water is heavy, so T and S push in opposite directions. R = 2 says how strong salinity is compared with temperature, and lambda is a resistance to the flow.
- f < 0 means the flow is temperature-driven (like today's overturning). f > 0 means it is salinity-driven and runs the other way.

### The code, block by block

1. Parameters at the top: R, delta and lambda, the f range to search, the time step and run length, and Stommel's published values for comparison.
2. Graphical method. At steady state the T and S equations give T = 1/(1+|f|) and S = delta/(delta+|f|). Putting these into the flow law gives one equation in f alone: lambda f = -1/(1+|f|) + R delta/(delta+|f|). I compute both sides on a fine grid of f. Stommel did the same thing by hand.
3. Finding crossings. Where the two sides cross, their difference changes sign. I look for each sign change between neighbouring grid points, then call `scipy.optimize.brentq` on that small interval. I use brentq because it is guaranteed to converge once you give it an interval with a sign change. Newton's method can jump away from the root.
4. Classifying. For each equilibrium I compute the 2x2 Jacobian (how dT/dt and dS/dt change when T or S is nudged) with centred finite differences, then its eigenvalues. Both negative and real: stable node. Complex with negative real part: stable spiral. One positive and one negative: saddle. I use finite differences because |f| makes the algebra messy.
5. Time integration. Forward Euler from 10 starting states: the corners and edge midpoints of the unit square, plus two points in between. Each path is coloured by the sign of f at the end.
6. Time step check. I rerun every path with dt/2. The end states agree to about 1e-13, but that alone proves little. A fixed point of the Euler update satisfies x + dt F(x) = x, which is just F(x) = 0, the true equilibrium, for any dt. So I also compare the whole path at the same times. The largest difference was 1.04e-2 with dt = 0.01 (from a test outside the notebook) and 1.91e-3 with dt = 0.002. That is about 5 times smaller for a 5 times smaller dt, which is the first-order behaviour expected from Euler. I kept dt = 0.002.

### Results

| | f (mine) | f (Stommel) | T | S | eigenvalues | type (mine) | type (Stommel) |
| :--- | ---: | ---: | ---: | ---: | :--- | :--- | :--- |
| 1 | -1.0679 | -1.1 | 0.4836 | 0.1350 | -3.610, -0.761 | stable node | stable node |
| 2 | -0.3070 | -0.30 | 0.7651 | 0.3518 | -2.849, +0.761 | saddle | saddle |
| 3 | +0.2191 | +0.23 | 0.8203 | 0.4321 | -0.912 ± 1.823i | stable spiral | stable spiral |

My values are close to Stommel's, which he read off a hand-drawn graph. All three types agree with his.

Eight of the 10 starting states ended in the temperature-driven state (f = -1.068). Two ended in the salinity-driven state (f = +0.219): the ones starting at (S, T) = (1.0, 1.0) and (0.5, 1.0).

### What to look for in the figures

- `figures/01_stommel_graphical.png`: the blue curve is the density difference at steady state, and the red line is lambda f. They cross three times. The middle crossing is the saddle.
- `figures/01_stommel_phase_plane.png`: blue paths go to the square (temperature-driven node) and orange paths go to the diamond (salinity-driven spiral). The orange paths loop around the diamond, which is what the complex eigenvalues mean. The saddle (X) sits on the boundary between the two groups. The dashed line is f = 0.

### What I learned

The same forcing allows two different stable circulations. Which one you get depends only on where you start. The saddle is the dividing point between the two basins, and the salinity-driven state is a spiral because its eigenvalues are complex.

## Step 2: Tipping and hysteresis in Stommel's model (`notebooks/02_stommel_hysteresis.ipynb`)

### The question

If I slowly make the flow harder to drive, does the temperature-driven circulation suddenly switch off? And if I undo the change, does it come back?

### The equations in plain words

The equations are the same as in step 1. The only change is that lambda, the resistance to the flow, is now a slowly changing control knob. A larger lambda means a weaker flow for the same density difference. That gives salinity more time to build up, and salinity works against the temperature-driven flow.

### The code, block by block

1. Shared root finder. I moved the step 1 method (scan for sign changes, then brentq) into `src/models.py` as `stommel_find_equilibria`, because this step calls it hundreds of times.
2. Equilibrium branches. For 701 values of lambda between 0.1 and 0.45 I find every equilibrium and classify it with the Jacobian. Each one goes into a branch by its sign of f and its stability. Gaps are filled with NaN so matplotlib does not draw lines between separate branches.
3. Finding the fold by bisection. At lambda = 0.1 a stable f < 0 state exists, and at 0.45 it does not. I keep halving the interval between the two and ask the root finder each time whether the state still exists. After 40 halvings the interval is far smaller than anything that matters.
4. An independent check of the fold. For f < 0 the equilibrium condition can be solved for lambda directly: lambda(f) = rhs(f) / f. The fold is the largest lambda on that curve, because past it no f < 0 works. I find this maximum with `scipy.optimize.minimize_scalar`. The two methods agree to 1e-6.
5. Stommel's own estimate. Stommel (1961, p. 229) read the fold off his hand-drawn graph as "slightly greater than 1/4". I count the f < 0 equilibria at lambda = 0.26 and 0.30, and evaluate lambda(f) at f = -0.5 as a quick hand estimate.
6. The sweep. I start on the temperature-driven equilibrium at lambda = 0.1, raise lambda in steps of 0.005 up to 0.45, then lower it back to 0.1. At each lambda I run forward Euler for 200 time units, starting from where the previous lambda ended. This is the key point. Starting from the previous end state is what a slowly changing climate does, and it is why the system remembers which branch it is on. That memory is the hysteresis. If I restarted from the same guess at every lambda, there would be no hysteresis at all.
7. Settling check. At the end of each lambda step I print how far from steady the state still is (the largest |dT/dt| or |dS/dt|). My first try used 100 time units per step and left a residual of 2.5e-2 (re-checked later outside the notebook: 2.47e-2). I traced it to lambda = 0.335, just past the fold. In a separate check run (not part of the notebook), f only drifted from -0.53 to -0.48 in the first 50 time units, was -0.25 at t = 100, and reached +0.20 by t = 150. The system lingers where the stable state used to be before it leaves. With 200 time units per step every step settles (largest residual 1e-11).
8. Time step check. I repeat both sweeps with dt/2. The recorded f values agree to 1e-12. As in step 1, this mainly shows that both runs ended on the same branch at every lambda, since settled Euler end states are exact equilibria for any dt.

### Results

- Fold of the temperature-driven state: lambda = 0.333800 by bisection and 0.333801 from the maximum of lambda(f), at f = -0.4837. My hand estimate, lambda(f) at f = -0.5, is exactly 1/3 = 0.3333, close to the maximum.
- Stommel's estimate, "slightly greater than 1/4", is too low: at lambda = 0.26 and 0.30 there are still two equilibria with f < 0 (the node and the saddle). He read it off a hand-drawn graph.
- Sweep up: f was -0.5331 at lambda = 0.330. At lambda = 0.335, the first step past the fold, it jumped to +0.2042.
- Sweep down: f stayed positive the whole way (smallest value +0.1937). At lambda = 0.1 it ended at +0.2328, even though the sweep started at -1.8542 at that same lambda.
- At lambda = 0.1 three equilibria exist: f = -1.8542, -0.2729 and +0.2328. Going up and back down moved the system from the first to the third.

### What to look for in the figure

- `figures/02_stommel_hysteresis.png`: the solid line along the bottom is the stable temperature-driven branch. The dashed line is the saddle. They meet at the fold (dotted line) and vanish together. The red circles follow the lower branch, then jump (arrow) to the upper branch. The blue crosses on the way back stay on the upper branch all the way to lambda = 0.1.

### What I learned

The temperature-driven state does not fade away. It disappears at a fold, where it meets the saddle, and the flow jumps to the salinity-driven state. Undoing the change does not bring it back, because the salinity-driven state is stable over the whole range. The system remembers its history. Stommel already made this point in 1961: after a slight increase in lambda the flow would jump to the salinity-driven state, and "would then stay there even when lambda was restored to its original value". Right past the fold the system moves very slowly before it jumps, which hints at the early warnings in step 5.

## Step 3: Cessi's model and its landscape (`notebooks/03_cessi_landscape.ipynb`)

### The question

Cessi reduced the two-box model to a single equation for the salinity difference. Does it show the same two states and the same hysteresis? Can I draw it as a ball rolling in a landscape?

### The equation in plain words

    dy/dt = p - y (1 + mu2 (1 - y)^2),   mu2 = 6.2

- y is the scaled salinity difference between the boxes. A small y means a strong overturning, and y near 1 means a weak one. The temperature difference is fixed at 1, so the density difference that drives the flow is proportional to 1 - y.
- p is the freshwater forcing (evaporation minus rain). It builds up the salinity difference.
- `y * 1` is slow diffusion, which wears the salinity difference down.
- `y * mu2 (1 - y)^2` is the salt carried away by the overturning flow. The flow is proportional to the density difference, 1 - y, and Cessi made the exchange rate grow with the square of the flow (Stommel used |f|). So a strong flow (small y) mixes salt away efficiently. That is the feedback: a weaker flow mixes less, so salinity builds up, so the flow gets weaker still. Squaring the flow keeps the equation smooth, which is why it has a potential (below).
- Time is measured in diffusion times, about 219 years each. mu2 = t_d / t_a is the diffusion time over the advection time (219 and 35 years in Cessi's parameter values, as listed by Yang et al. 2021). The temperature difference can be fixed at 1 because it relaxes in about 25 days, far faster than salinity.

Because there is only one variable, the equation can be written as dy/dt = -dV/dy with a potential V(y). The system behaves like a ball rolling downhill in V, with strong friction. Valleys are stable states and hills are unstable ones.

### The code, block by block

1. New functions in `src/models.py`: the tendency, p(y), its slope dp/dy, the potential V, the equilibria and the fold points.
2. Equilibrium diagram. Finding y for a given p means solving a cubic, but finding p for a given y is one formula: p = y (1 + mu2 (1 - y)^2). So I pick 2001 values of y and compute p, which draws the whole curve in one go, including the unstable middle part. A point is stable where dp/dy > 0. The tendency is p minus p(y), so if y is nudged up, p(y) rises above p and the tendency pushes y back down.
3. Fold points. The folds are where dp/dy = 0. dp/dy is a quadratic in y, so the quadratic formula gives both folds exactly. No root finder is needed.
4. Valleys and hill at p = 1.1. The equilibria solve the cubic mu2 y^3 - 2 mu2 y^2 + (1 + mu2) y - p = 0, so `np.roots` gives all three at once. The barrier height is V(hill) minus V(valley). Step 4 needs these numbers.
5. Hysteresis ramp. p goes from 0.8 to 1.5 over 500 diffusion times (109,500 years), then back down over the same time. Here p changes a little at every Euler step, so the system never fully rests. The ramp is slow compared with the relaxation time (about half a diffusion time), so y stays close to its branch until the branch ends. I locate each jump as the p where y changes fastest.
6. Time step check. This time y is always moving, so comparing dt with dt/2 tests a real transient, unlike the end-state checks in steps 1 and 2. The largest difference in y along the ramp was 3.42e-3, and the jump locations agree to 4 decimals.
7. Slower ramp. The jumps happened a bit past the folds. I reran with a ramp 4 times slower to see if this is a delay. Near a fold the pull back toward equilibrium is very weak, so y takes time to leave.

### Results

| | mine | published | source |
| :--- | ---: | ---: | :--- |
| fold where the strong branch ends | p = 1.2962 (y = 0.4272) | p about 1.3 | Cessi (1994), Fig. 4 |
| fold where the weak branch ends | p = 0.9556 (y = 0.9061) | none found | quadratic formula only |
| strong valley at p = 1.1 | y = 0.2402 | 0.2402 | Yang et al. (2021) |
| hill at p = 1.1 | y = 0.6911 | 0.6911 | Yang et al. (2021) |
| weak valley at p = 1.1 | y = 1.0687 | 1.0687 | Yang et al. (2021) |

Yang et al. (2021) use Cessi's parameters (mu^2 = 6.2, p = 1.1). Cessi's Fig. 4 shows the potential at p = 1.3 with mu^2 = 6.2, where the strong valley has just become a flat shoulder.

- At p = 1.1 the barrier from the strong valley to the hill is 0.05710. From the weak valley it is 0.03560. So the strong state is the deeper valley at this forcing.
- Ramp up: the jump came at p = 1.3144, 0.0182 past the fold. Ramp down: at p = 0.9375, 0.0182 below the other fold.
- With a ramp 4 times slower, the overshoot fell to 0.0072 in both directions. The overshoot is a delay caused by the finite ramp speed, not an error. 0.0182 / 0.0072 = 2.53, which is close to 4^(2/3) = 2.52, the scaling expected for a slow ramp through a fold. I only tried two ramp speeds, so this is a consistency check, not a proof.

### What to look for in the figures

- `figures/03_cessi_equilibria.png`: an S-shaped curve. Between p = 0.956 and 1.296 there are three equilibria, two stable (solid) and one unstable (dashed). Outside that range there is only one.
- `figures/03_cessi_potential.png`: at p = 1.0 the strong valley (small y) is the deep one. As p grows, the landscape tilts toward large y. At p = 1.35 the strong valley has disappeared, leaving only a flat shoulder, so a ball sitting there must roll to the weak state.
- `figures/03_cessi_hysteresis.png`: red (p increasing) follows the lower branch and jumps up just after the right fold. Blue (p decreasing) follows the upper branch and jumps down just after the left fold. The loop between them is the hysteresis.

### What I learned

Cessi's single equation has the same structure as Stommel's model: two stable states, an unstable one between them, and two folds. The potential makes this easy to picture. Changing p tilts the landscape until a valley disappears. A ramp always overshoots the fold a little, and the overshoot gets smaller as the ramp gets slower.

## Step 4: Adding weather noise (`notebooks/04_cessi_noise.ipynb`)

### The question

Real freshwater forcing is not smooth. It has random weather on top. Can noise alone push the circulation from one state to the other? How does the waiting time between flips depend on the noise strength?

### The equation in plain words

    dy = (p - y (1 + mu2 (1 - y)^2)) dt + sigma dW

The first part is Cessi's model from step 3, with p = 1.1, inside the range where two states exist. The new part, sigma dW, is a random kick (white noise) standing for weather. sigma sets its size. In the landscape picture, the ball is now being shaken. Usually it rattles around inside its valley, but now and then a run of kicks in the same direction pushes it over the hill.

### The code, block by block

1. Landscape at p = 1.1. I recompute the valleys and the hill as in step 3, and the two flip thresholds: low = 0.4656 (halfway from the strong valley to the hill) and high = 0.8799 (halfway from the hill to the weak valley).
2. Euler-Maruyama. Each step does `y_new = y + F(y) dt + sigma sqrt(dt) * N(0,1)`. The noise term scales with sqrt(dt), not dt. Independent random kicks add up like a random walk, so their spread grows like the square root of the number of steps. For the noise to have the same total effect over a fixed time whatever dt is, the variance of each kick must be proportional to dt, so its standard deviation is proportional to sqrt(dt). With dt instead, the noise would fade away as dt gets smaller. I generate all the random numbers at once with a fixed seed. That is faster than one call per step, and every run gives the same numbers.
3. Counting flips with two thresholds. A flip from strong to weak counts only when y climbs above the high threshold, and a flip back only when it drops below the low one. With a single threshold at the hill, every wiggle across the hilltop would count as a flip. The time spent in the strong state runs from entering it (crossing the low threshold downward) to leaving it (crossing the high threshold upward).
4. Histogram vs theory. For a one-variable model like this, the long-run distribution of y is known exactly: P(y) proportional to exp(-2 V(y) / sigma^2). Deep valleys are visited most. I normalize it to unit area with `np.trapezoid` so it can sit on top of a `density=True` histogram. For the uncertainty of the fraction of time below the hill I use batch means: I split the run into 20 equal parts (about 20 flips each) and use the spread of their 20 fractions.
5. Kramers test. At 6 noise levels from 0.18 to 0.30 I run 20,000 diffusion times (4.38 million years) each, with a different seed per level, and record the mean time spent in the strong state. Kramers' (1940) theory says the mean waiting time is tau = 2 pi / sqrt(V''(valley) |V''(hill)|) * exp(2 barrier / sigma^2). So ln(tau) against 1/sigma^2 should be a straight line with slope 2 x barrier. I fit a straight line with `scipy.stats.linregress`, which also gives the standard error of the slope from the scatter of the 6 points, and draw Kramers' full formula. V'' is the same dp/dy used for stability in step 3.
6. The exact answer. For one variable, the mean time to go from a to b has an exact formula, the mean first-passage time: T = (2/sigma^2) times the integral from a to b of exp(2V(u)/sigma^2) times [the integral from minus infinity to u of exp(-2V(z)/sigma^2)]. Kramers' formula is its weak-noise limit. I compute it from the low threshold, where my timer starts, and from the bottom of the valley.
7. Time step check. A long run gives only a few hundred waits, so its mean is uncertain by 3-9%, too much to see an effect of a few percent. So I start 20,000 walkers exactly at the low threshold at sigma = 0.30, time how long each takes to cross the high threshold, and compare dt = 0.01, dt = 0.0025 and the exact value. (My first version compared two long runs at sigma = 0.24: 28.17 ± 1.46 with dt = 0.01 and 26.79 ± 1.19 with dt = 0.005. They agreed, but they were too uncertain to show anything.)

### Results

- Example run (sigma = 0.2, 4.38 million years): 411 flips. Mean time in the strong state was 71.5 diffusion times (15.7 kyr), and in the weak state 25.4 diffusion times (5.6 kyr). y was below the hill 73.9 ± 2.4% of the time, against 72.4% from the theoretical distribution, so they agree within the uncertainty.
- Kramers test:

| sigma | 1/sigma^2 | flips | mean wait | std error | exact, from threshold | exact, from valley | Kramers formula |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.180 | 30.86 | 234 | 130.74 | 11.69 | 135.51 | 142.72 | 136.41 |
| 0.204 | 24.03 | 492 | 56.50 | 3.37 | 60.00 | 64.69 | 62.50 |
| 0.228 | 19.24 | 827 | 32.04 | 1.55 | 33.23 | 36.61 | 36.15 |
| 0.252 | 15.75 | 1212 | 21.00 | 0.86 | 21.19 | 23.79 | 24.27 |
| 0.276 | 13.13 | 1499 | 16.85 | 0.63 | 14.82 | 16.91 | 18.00 |
| 0.300 | 11.11 | 1997 | 11.91 | 0.39 | 11.06 | 12.80 | 14.29 |

  (Waiting times in diffusion times. "Exact" is the exact mean first-passage time to the high threshold.)

- Fitted slope: 0.1188 ± 0.0026 (standard error; 95% interval 0.1117 to 0.1260). 2 x barrier from step 3: 0.1142, inside the interval.
- Time step check at sigma = 0.30 (20,000 walkers): 11.518 ± 0.086 with dt = 0.01, 11.144 ± 0.083 with dt = 0.0025, exact 11.057. So dt = 0.01 makes the passage about 4% too slow at this noise level, and the error shrinks with dt.

### What to look for in the figures

- `figures/04_cessi_noise_timeseries.png`: the first 438 kyr of the run. y rattles around one valley, then suddenly jumps to the other. Stays in the strong state (near 0.24) are longer than stays in the weak state (near 1.07).
- `figures/04_cessi_noise_histogram.png`: two peaks at the two valleys and a dip at the hill. The red curve is exp(-2V/sigma^2) with no fitting at all. The simulated weak-state peak is slightly lower than the theory. The fraction of time below the hill agrees with theory within its uncertainty (0.739 ± 0.024 against 0.724), so sampling error is enough to explain it.
- `figures/04_cessi_kramers.png`: the points lie close to a straight line. The fitted line (blue) and Kramers' formula (red dashed) have nearly the same slope. The purple dotted line is the exact mean passage time between my two thresholds. Kramers' formula sits above the points for two reasons. First, it describes escape from the bottom of the valley, while my timer starts at the low threshold, partway up the barrier: the exact time from the valley bottom is 5-16% longer than from the threshold. Second, it assumes the barrier is large compared with the noise, which it is not here (2 x barrier / sigma^2 = 1.3 to 3.5). Compared with the exact time from the valley bottom, Kramers' formula is 4% too short at sigma = 0.18 and 12% too long at sigma = 0.30. The simulations agree with the exact line within 1.1 standard errors at the four weaker noise levels. At sigma = 0.276 and 0.30 they are 14% and 8% above it (3.2 and 2.2 standard errors); the time step check shows that about 4% of this comes from dt = 0.01.

### What I learned

Noise alone can flip the circulation between its two states, with no change in the average forcing. The waiting time grows exponentially with barrier / sigma^2. A slightly deeper valley or slightly weaker noise means a much longer stay, which is why the strong state (barrier 0.057) holds the system longer than the weak state (barrier 0.036). The measured slope, 0.119 ± 0.003, is consistent with 2 x barrier = 0.114. When the noise is not weak, the exact mean first-passage time is a better benchmark than Kramers' formula, and the place where the timer starts matters.

## Step 5: Early warning signals (`notebooks/05_early_warning.ipynb`)

### The question

If the forcing slowly pushes the circulation toward a tipping point, can we see it coming in the fluctuations, before the jump happens?

### The equation in plain words

This is the noisy Cessi model from step 4, but now p rises slowly from 1.0 to 1.4, past the fold at p = 1.296 where the strong state disappears. The noise is weak (sigma = 0.04). It is too small to cause random flips far from the fold. It just nudges y a little all the time, and those nudges let me measure how strongly the system is pulled back.

The key quantity is the recovery rate r: how fast a small nudge away from the strong state decays. It is minus the slope of the tendency at the strong equilibrium, which is dp/dy from step 3. At a fold dp/dy = 0, so r goes to zero. A system that recovers more slowly lets each kick last longer, so its fluctuations get bigger (higher variance) and each value looks more like the one before (higher autocorrelation). This is called critical slowing down.

### The code, block by block

1. Recovery rate. For p from 1.0 to the fold I find the strong equilibrium and compute r = dp/dy there.
2. Noisy ramp. Euler-Maruyama as in step 4, with p rising a little every step over 4000 diffusion times (876,000 years). It starts at the strong equilibrium for p = 1.0.
3. Tipping time. The tip is the first time y goes above 0.7. The strong branch never goes above 0.43 (its fold) and the weak branch is above about 0.9, so 0.7 is clearly between them.
4. Sampling like a record. I keep one value every 0.5 diffusion times (about 110 years). With every Euler step (dt = 0.01), neighbouring values would be almost identical, and the lag-1 autocorrelation would sit near 1 the whole time. I stop the analysis 10 diffusion times before the tip. On my first try the statistics shot up in the last few points, because the escape itself (y climbing out of the valley) was getting into the windows. That is the tip, not a warning of it.
5. Removing the trend. As p rises, the strong state itself moves to larger y, so the mean of y creeps upward. That slow drift is not a fluctuation. Left in, it would add to the variance and push the autocorrelation toward 1 even if the recovery rate did not change. So I subtract a centred rolling mean 50 diffusion times wide. This is much longer than the recovery time (under 3 diffusion times here) but short compared with the ramp.
6. Rolling statistics. In a trailing window of 250 diffusion times (500 samples) I compute the variance and the lag-1 autocorrelation (the correlation between each value and the next). A trailing window uses only data up to its own time. The test as a whole is still retrospective: the detrending is centred (it uses up to 25 diffusion times of later data), and I cut the record using the known tip time. Dakos et al. (2008) analysed past climate records in the same way, and Held and Kleinen (2004) proposed watching the recovery rate of the overturning circulation for the same reason.
7. Kendall's tau. A rank-based number for the trend: +1 if a statistic only ever goes up, 0 if there is no trend. It is the usual way to put a number on an early warning.
8. Linear theory. Close to a stable state, small fluctuations behave like a random walk pulled back at rate r. For that process, variance = sigma^2 / (2 r) and lag-1 autocorrelation over a gap of 0.5 is exp(-0.5 r). r changes within a window, so I evaluate the theory at every sample and average it over the same window as the data, and draw that as dotted lines. (My first version used r at each window's last point, which overestimates the variance near the fold.)
9. 100 runs and a control. One run can be lucky, so I repeat the whole analysis for 100 runs with different random numbers. As a control I make 100 records at a fixed p = 1.1, 2700 diffusion times long, where nothing approaches a fold, and analyse them in the same way. The function that does the analysis is first checked against the single run.

### Results

- The noisy system tipped at p = 1.2695 (t = 2695.3 diffusion times, 590 kyr). That is 0.0267 before the fold at p = 1.2962. Near the fold the barrier is so low that the noise pushed y over the hill early.
- Recovery rate at the strong state:

| p | r | recovery time 1/r |
| ---: | ---: | ---: |
| 1.00 | 2.947 | 74 years |
| 1.10 | 2.316 | 95 years |
| 1.20 | 1.540 | 142 years |
| 1.25 | 1.023 | 214 years |
| 1.28 | 0.580 | 378 years |
| 1.29 | 0.349 | 627 years |

- First window (ending at p = 1.025): variance 2.94e-4 (theory 2.79e-4), lag-1 autocorrelation 0.223 (theory 0.238).
- Last window (ending at p = 1.268): variance 9.54e-4 (theory 8.56e-4), lag-1 autocorrelation 0.655 (theory 0.625).
- Kendall tau: 0.880 for variance and 0.735 for lag-1 autocorrelation.
- 100 ramp runs: Kendall tau 0.775 to 0.912 for the variance and 0.664 to 0.855 for the autocorrelation (5-95% range; smallest values 0.712 and 0.629). 100 control runs at fixed p = 1.1: -0.301 to 0.337 and -0.322 to 0.350 (largest values 0.556 and 0.488). The runs tipped at p = 1.2769 (median), 5-95% range 1.2628 to 1.2863, so the run shown tipped a little early.

### What to look for in the figures

- `figures/05_early_warning_panels.png`: top, y creeps up along the strong branch, then jumps at the dashed line. Middle and bottom, the variance and autocorrelation rise toward the tip, roughly along the dotted theory curves. Near the tip the measured variance is about 10% above the window-averaged theory. Linear theory treats the valley as a parabola, which becomes less accurate as the valley flattens near the fold; I did not test this further. The autocorrelation is noisier than the variance, and has a dip around 150 kyr and a bump around 200 kyr. The theory curve is smooth there, so these are most likely sampling noise from a single run.
- `figures/05_early_warning_ensemble.png`: Kendall tau for the 100 ramp runs (blue) and the 100 control runs (grey). The two groups do not overlap: every ramp run has a larger tau than every control run, for both statistics. The dashed line is the single run above.
- `figures/05_recovery_rate.png`: r falls from about 3 to 0 at the fold. It drops steeply in the last part, so most of the warning comes late. The blue dotted line shows where the noisy run actually tipped, before r reached zero.

### What I learned

The warnings come from one cause: the recovery rate goes to zero at the fold, so the system responds more and more slowly to small kicks. Both variance and autocorrelation rose before the tip and roughly followed simple linear theory, and over 100 runs the rise was always clear. Three cautions. Removing the trend matters, because a drifting mean would fake a rising variance. The test is retrospective, because it uses the known tip time. And in the run shown, noise tipped the system before the fold was reached. When noise alone causes a jump far from any fold, these warnings need not appear at all (Ditlevsen and Johnsen 2010). I only tried one noise level and one ramp speed.

## Step 6 (optional): A toy Antarctic ice-shelf margin (`notebooks/06_antarctic_toy_model.ipynb`)

### The question

Near an Antarctic ice shelf, could warm water at depth and meltwater at the surface create a threshold in the far-field ocean temperature, the same kind of fold as in steps 2 and 3? This model is my own toy. It is not taken from a paper, and every parameter is an order-of-magnitude guess. I wanted to see whether a threshold appears with reasonable values, and to report honestly if it does not. Thresholds of this kind are known from regional ocean models: Hellmer et al. (2012) found a cold ice-shelf cavity warming abruptly in the 21st century, and Hazel and Stewart (2020) found that the Filchner-Ronne cavity can sit in a warm or a cold state under the same forcing. There the switch depends on which water mass fills the cavity. My toy model tests a different feedback, through surface freshening and vertical mixing.

### The equations in plain words

    dT/dt = (T_far - T)/tau_in  -  kappa(F) T  -  b m T^2
    dF/dt = c m T^2  -  F/tau_F

- T is the subsurface temperature above freezing (deg C). F is the surface freshening (psu).
- `(T_far - T)/tau_in`: warm water flowing in from the open ocean pulls T toward the far-field temperature T_far, over about tau_in = 1 year.
- `kappa(F) T`: vertical mixing carries heat up to the cold surface layer, which sits at the freezing point (T = 0). Freshening makes the surface lighter and weakens this mixing.
- `m T^2`: the basal melt rate, quadratic in temperature above freezing. Melting uses heat (the `b` term) and freshens the surface (the `c` term).
- `F/tau_F`: currents and sea-ice growth remove the freshening over about tau_F = 1 year.
- The feedback: warmer water melts more ice, which freshens the surface, which weakens mixing, which traps more heat below.

The coefficients b and c come from simple budgets, worked out in the parameter cell of the notebook. c = 0.0154 psu of freshening per m of ice melted, from the salinity, the area of ice-shelf base compared with the shelf sea, the depth of the surface layer, and the density of ice compared with seawater (melting 1 m of ice adds 917 kg of fresh water per m^2, not 1 m of water). b = 0.0187 deg C of cooling per m of melt, from the latent heat of ice and the heat capacity of the subsurface layer.

### The code, block by block

1. New functions in `src/models.py`: the mixing law, the melt rate, the tendencies, the Jacobian, and T_far as a function of T at steady state. The parameters travel in one dictionary because there are many of them.
2. Two mixing laws. kappa0 / (1 + F/F0) and kappa0 exp(-F/F0) both say "freshening weakens mixing". No data I have can tell them apart, so I run both.
3. Equilibrium curve. Setting dF/dt = 0 gives F as a function of T. Putting that into dT/dt = 0 and solving for T_far gives T_far(T) directly, the same "sideways" trick as Cessi's p(y). There is a fold only where the slope of T_far(T) changes sign. I classify stability with the Jacobian, as in step 1.
4. Sweep. T_far goes from 0 to 5 deg C and back in steps of 0.05, with 100 years at each step, each step starting from the previous end state, as in step 2. My first try used 30 years per step, and one point on the way down was caught halfway, just past a fold (re-checked later outside the notebook, exponential law: the largest residual was 9.3e-2 per year with 30 years per step, against 1.04e-4 with 100). That is the same slow escape I saw in step 2.
5. Smallest mixing rate for a fold. Bisection on kappa0, with F0 fixed at its baseline value.
6. Parameter scan. For 60 x 60 combinations of kappa0 (0.5 to 30 per year) and F0 (0.02 to 2 psu) I check for a fold, and whether it falls inside 0-4 deg C, roughly the range of far-field temperatures around Antarctica.
7. Why the two laws differ, by hand. Without the melt cooling term (b = 0) the steady state is T_far = T + tau_in kappa(F(T)) T, with F growing like T^2. A fold needs the heat loss kappa T to fall as T rises, faster than 1/tau_in. Writing u for F/F0: for the rational law kappa T falls at most at the rate kappa0/8 (at u = 3), and for the exponential law at most at 2 exp(-1.5) kappa0 = 0.446 kappa0 (at u = 1.5). So without melt cooling a fold needs kappa0 tau_in > 8 or > 2.24, whatever F0, c, m and tau_F are. I check both numbers with the bisection, with b set to 0.
8. A like-for-like comparison. At the same F0 the exponential law always mixes less than the rational one (exp(-u) < 1/(1+u)). So I also run the exponential law with F0 = 0.2 / ln 2 = 0.289 psu, which halves the mixing at F = 0.2 psu, just like the rational law with F0 = 0.2 psu.

### Results

- With the rational mixing law and my baseline guesses (kappa0 = 5 per year, F0 = 0.2 psu), there is no fold. The sweeps up and down lie exactly on top of each other. A fold needs kappa0 of at least 13.66 per year. With tau_in = 1 year, that means vertical mixing would have to remove heat about 14 times faster than the inflow brings it in.
- With the exponential mixing law and the same parameters, there is a fold. Thresholds sit at T_far = 3.505 deg C (going up) and 3.149 deg C (going down), and the sweeps differ by up to 1.40 deg C. Here a fold needs only kappa0 of at least 3.39 per year.
- Without melt cooling, the smallest kappa0 for a fold is 8.000 and 2.241 per year, the same as the hand values. Melt cooling raises them to 13.66 and 3.39.
- The exponential law with F0 = 0.289 psu still has a fold, at 3.923 and 4.258 deg C, so the difference between the laws is not only the choice of F0. But the threshold moves to warmer water.
- At T_far = 2.00 on the way up, both versions give T near 0.36 deg C, F near 0.02 psu and a melt rate of 1.3 m/yr. The shelf water is much colder than the far field, because strong mixing removes most of the heat.
- In the scan, thresholds inside 0-4 deg C appear only in a band of parameter space: small F0 for the rational law, and a wider band for the exponential law.

### What to look for in the figures

- `figures/06_antarctic_sweeps.png`: left (rational law), a single smooth curve with no jump. Right (exponential law), an S-shaped curve with a cold branch, a warm branch, an unstable middle, and a jump at each fold.
- `figures/06_antarctic_parameter_scan.png`: grey means no threshold. Light orange means there is a threshold, but at a far-field temperature above 4 deg C. Red means a threshold inside 0-4 deg C. The star is my baseline. It sits in grey on the left and in red on the right. F0 means something slightly different in the two panels: the rational law halves the mixing at F = F0, the exponential law reduces it by a factor e.

### What I learned

The feedback can create a threshold, but whether it does depends on things I cannot pin down: how strong vertical mixing is, and how sharply freshening weakens it. The hand calculation shows what decides it: how fast the vertical heat loss can fall as the shelf warms, compared with the inflow rate. With my baseline guesses, one reasonable mixing law gives a threshold near 3-3.5 deg C and the other gives none. So this toy model cannot say whether a real ice-shelf margin has a threshold. It does show which processes would decide it. There is also a tension. A threshold needs strong vertical heat loss, but strong heat loss keeps the shelf water much colder than the far field (0.36 deg C against 2 deg C here). So in this model, thresholds belong to cold shelves, not to shelves already flooded with warm water.

## What I would do next

- Run the same hysteresis experiment in a real ocean model. In MITgcm I would set up an idealized basin, add a slowly increasing freshwater flux at high latitudes, then reverse it, and see whether the overturning shows a fold and hysteresis like the box models. A full ocean model has winds, rotation and geography, and these could change where the fold is or whether there is one at all.
- Add an ice sheet. Meltwater from an ice sheet would make the freshwater forcing depend on the ocean state itself, and that extra feedback could move the tipping point or create new ones.
