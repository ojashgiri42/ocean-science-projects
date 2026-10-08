# Tutorial: melting under an ice shelf in MITgcm

These notes walk through each step of the project in order. I wrote them so I can explain every command and every choice.

The MITgcm source code is the clone in `~/models/MITgcm` from project 15 (release `checkpoint69q`, commit 853761d). I never change files inside it. The `isomip` experiment's `code/`, `input/` and `results/` folders are copied to `~/mitgcm_runs/isomip/source_copy/`, and every build and run lives in `~/mitgcm_runs/isomip/`, outside the git repo. The repo only holds scripts, the files I change, notebooks, figures and small results.

## Step 1: Build and verify (`scripts/build.sh`, `scripts/run_experiment.sh`, `notebooks/01_build_and_verify.ipynb`)

### The question

Does my build of MITgcm's `isomip` experiment give the same answer as the official reference run? And what exactly is this setup?

### What I did

Nothing in the model was changed in this step.

- I copied `build.sh` and `run_experiment.sh` from project 15, with all their fixes. The most important one is that `build.sh` adds `/opt/homebrew/lib` to `LIBRARY_PATH`, so the linker finds the netCDF library. `isomip` also compiles `pkg/mnc`, so it needs this too.
- I built the experiment with the same build-options file as in project 15, `darwin_arm64_gfortran`, because `uname -m` still says `arm64`. I built it twice: an "IEEE" build (strict arithmetic, no optimization) to compare with the reference, and a "fast" build (`-O3`) for long runs. Both built first time.
- I ran the default test (20 steps) exactly as provided and compared it with `results/output.txt`.
- I timed one model month.

### The setup in plain words

- Domain: from 80°S to 70°S and over 15° of longitude, on a latitude-longitude grid with 0.1° steps north-south and 0.3° east-west (100 x 50 points). That is about 1112 km north-south and 291-569 km east-west, because longitude lines get closer toward the pole. The cells are 11.1 km north-south and 5.8-11.4 km east-west.
- Vertical: 30 levels, each 30 m thick, 900 m in total. The seafloor is flat at 900 m.
- Walls: the first column (west) and first row (south) are land. MITgcm's grid wraps around by default (east connects to west, north to south), so these walls also close the eastern and northern sides. The domain is a sealed box.
- The ice shelf covers the whole domain. Its base (the ice draft) is about 694 m deep at the southern wall, rises steadily to 200 m by about 76°S, then stays flat at 200 m all the way to 70°S. The water under the ice is therefore about 219 m thick near the southern wall and 700 m thick in the flat part. This is ISOMIP "Experiment 1".
- Time step: `deltaT = 1800` s (30 minutes), the same for everything. There is no accelerated tracer step like in project 15. The default test runs 20 steps, 10 model hours.
- Initial state: the ocean is at rest, with uniform temperature -1.9 °C and salinity 34.4 g/kg everywhere (`Tref`, `Sref`).
- Forcing: none. There is no restoring of temperature or salinity anywhere, no open boundary and no surface forcing. The only thing that can drive a circulation is melting and freezing at the ice base. So the far field is not held fixed: anything that happens slowly uses up the water's heat.
- Other settings: equation of state `JMD95Z`, `rhoConst = 1030` kg/m^3, heat capacity 3974 J/kg/K, horizontal mixing 100 m^2/s for heat and salt, vertical mixing 5e-5 m^2/s.

### The settings in `data.shelfice`

I read the package's source (`pkg/shelfice/shelfice_readparms.F`) to see the default of every setting the file does not give.

- `SHELFICEboundaryLayer = .TRUE.`: the melt is computed from the temperature and salinity averaged over a layer one cell thick right under the ice, instead of just the top wet cell. Under a sloping ice base the top cell can be very thin, and its values jump around. Averaging over a fixed thickness makes the melt rate smoother (Losch 2008).
- `SHELFICEtopoFile = 'icetopo.exp1'`: the depth of the ice base at every point.
- `SHELFICEloadAnomalyFile = 'phi0surf.exp1.jmd95z'`: the pressure of the floating ice on the water below, as an anomaly. A floating ice shelf pushes down on the ocean with its weight. The model needs this pressure, otherwise the water would try to flow up into the "space" the ice takes.
- `useISOMIPTD = .TRUE.`: use the simple melt formula of the original ISOMIP protocol (more below). The default is `.FALSE.`, which gives the three-equation model.
- `SHELFICEuseGammaFrict = .FALSE.`: the heat and salt transfer coefficients are constants. With `.TRUE.` they would depend on the flow speed next to the ice (the friction velocity), following Holland and Jenkins (1999).
- `no_slip_shelfice = .false.`: the ice base is free-slip: water can slide along it without friction from the wall itself. A drag is still applied through the drag coefficient.
- `SHELFICEwriteState = .TRUE.`: write the package's fields, such as the melt (freshwater) flux, when the model writes snapshots.

Settings that are not in the file, so their defaults apply:

- `SHELFICEheatTransCoeff` (gamma_T) = 1.0e-4 m/s. This is the heat transfer coefficient: how fast heat moves from the water through the thin boundary layer to the ice.
- `SHELFICEsaltTransCoeff` (gamma_S) = `SHELFICEsaltToHeatRatio` x gamma_T = 5.05e-3 x 1.0e-4 = 5.05e-7 m/s. Salt crosses the boundary layer about 200 times more slowly than heat, because salt diffuses much more slowly than heat in water.
- Latent heat of ice 334000 J/kg, ice density 917 kg/m^3, ice heat capacity 2000 J/kg/K, heat conduction into the ice 1.54e-6 m^2/s, and an ice surface temperature of -20 °C (used for the heat conducted up into the ice).
- `shiCdrag = 0.0015` and other Holland-and-Jenkins constants, only used when `SHELFICEuseGammaFrict = .TRUE.`.

So in this default setup the melt does not depend on the flow speed next to the ice: gamma_T and gamma_S are constants.

### How a floating ice shelf fits into a z-level ocean model

MITgcm normally has fixed horizontal levels from the sea surface down (a "z-level" model). To put an ice shelf in, the model treats the ice base as a lid pushed down into the ocean:

- At each point, the cells above the ice base are "dry" (they are ice, not water). The first wet cell under the ice (`kTopC` in the code) plays the role the surface cell plays in the open ocean. Heat, salt and freshwater fluxes from melting enter there.
- The ice base usually falls between two levels. Then the top wet cell is only partly open, like the seafloor's partial cells (`hFacMin = 0.05` here, so a cell can be as little as 5% open). The ice base becomes a staircase of full and partial cells.
- The ice's weight is applied as a pressure on the water (the load anomaly file). Then a floating shelf is in balance with the water underneath it.
- Melting is applied as a flux of heat and fresh water at the ice base, and the ice is fixed in shape. The ice does not actually get thinner in the model. That is a key limitation.

### The "three-equation" melt model

At the ice base, three things have to be true at the same time. Unknowns: the melt rate m, the temperature at the ice base T_b and the salinity at the ice base S_b.

1. Heat balance. Heat carried from the ocean to the ice base equals heat used for melting plus heat conducted up into the cold ice:
   rho_w c_w gamma_T (T_ocean - T_b) = rho_i L m + (heat conducted into the ice)
   gamma_T sets how fast heat crosses the thin layer of water touching the ice.
2. Salt balance. The salt carried from the ocean to the ice base is used to salt the fresh meltwater up to S_b (the ice itself has no salt):
   rho_w gamma_S (S_ocean - S_b) = rho_i m S_b
3. Freezing point. The ice base is at the freezing temperature of the water there, which depends on its salinity and pressure:
   T_b = a S_b + b p + c
   In the code: a = -0.0575 °C per (g/kg), b = -7.61e-4 °C per dbar, c = 0.0901 °C.

Putting the three together gives a quadratic equation for S_b, which the code solves at every point under the ice every time step (`shelfice_thermodynamics.F`). Then T_b and m follow.

The simple ISOMIP formula that this experiment uses by default skips the salt balance. It computes the freezing point from the ocean's own salinity, and turns all the heat that reaches the ice into melt:

heat flux = rho c gamma_T (T_ocean - T_freezing(S_ocean, p))

I verify the default exactly as provided. From step 2 on, I switch to the three-equation model (`useISOMIPTD = .FALSE.`), because it is what most current ice-ocean models use.

### Why the freezing point drops with depth, and why -1.9 °C water can melt ice deep down

Pressure makes it harder for ice to form, so seawater freezes at a lower temperature at depth, about 0.75 °C lower per 1000 dbar (roughly per 1000 m). For salinity 34.4 g/kg (ISOMIP formula):

- 0 m: freezing point -1.888 °C;
- 200 m (202 dbar): -2.040 °C;
- 694 m (701 dbar): -2.416 °C.

The model's water starts at -1.9 °C. At the surface that is just below freezing, by 0.012 °C. Under the ice it is above the local freezing point: by 0.140 °C at the 200 m ice base and by 0.516 °C at the 694 m deep end. So water that would freeze at the surface can still melt ice deep in the cavity. That is the start of the "ice pump", which I look at in step 2.

### Comparing with the reference

The model starts at rest with perfectly uniform water, so in the first steps many monitored quantities are pure rounding noise. For example, the spread of temperature at step 0 is 8e-14, from a field that should be perfectly uniform. Comparing the digits of noise means nothing, so I compare monitor times 10 to 20, after the flow has started. I also leave out the mean sea-surface height, which the model keeps at zero; its values are noise of about 1e-21. Quantities printed only at the start (grid statistics) are left out too, which leaves 50 quantities.

- IEEE build: 13 of 50 identical in every digit, the worst agrees to 9.6 digits, the median is 13.4 digits.
- Fast build: 10 of 50 identical, the worst agrees to 7.5 digits, the median is 10.8 digits.

This is a bit less agreement than in project 15 (10-16 digits). At the end of the test the flow is still very weak (largest speed 1.79e-4 m/s), so tiny arithmetic differences show up more. Agreement to 7.5-9.6 digits still means my build solves the same equations the same way.

### Timing

One model month (1440 steps of 30 minutes, `experiments/timing_1month/`):

- fast build: 39.3 s, so about 8.0 minutes per 365-day model year on one core;
- IEEE build: 413.8 s, 10.5 times slower.

So the long runs use the fast build. This setup has no calendar, so I use a 365-day year for melt rates and planning, and say so wherever it matters.

### What to look for in the figures

- `figures/01_geometry_section.png`: the ice base sloping up from about 694 m at the southern wall to 200 m at 76°S, then flat, over a flat 900 m seafloor. It is the same at every longitude.
- `figures/01_freezing_point_depth.png`: the freezing point gets colder with depth. The two formulas almost overlap. The red line (-1.9 °C) is to the right of both curves everywhere below about 15 m, so the water starts above its local freezing point everywhere under the ice.
- `figures/01_reference_match.png`: matching digits for all 50 compared quantities, sorted. The IEEE build agrees better than the fast build, as expected.

### What I learned

The default `isomip` case is a sealed box completely covered by an ice shelf, with no forcing at all. Melting and freezing are the only things that can make the water move. The pressure-dependent freezing point is the key: the same -1.9 °C water is too cold to melt ice at the surface but warm enough to melt it at 700 m. My build reproduces the reference run, and the default uses a simplified melt formula, so I will switch to the three-equation model before the real experiments.

## Step 2: The baseline cavity (`experiments/baseline/`, `notebooks/02_baseline_cavity.ipynb`)

### The question

Where and how fast does the ice shelf melt in the baseline setup? How does the water under the ice move, and is there an "ice pump"?

### What I changed and why

Compared with the original `isomip` setup:

- `data.shelfice`: `useISOMIPTD = .FALSE.`, which switches to the three-equation melt model (step 1). The transfer coefficients stay at their constant defaults (gamma_T = 1e-4 m/s, gamma_S = 5.05e-7 m/s).
- `data`:
  - 5 model years of 365 days (87600 steps);
  - one output file per field (`globalFiles = .TRUE.`);
  - a yearly restart file;
  - the monitor and its bounds check every 30 days.
- `data.pkg` and `data.diagnostics`: monthly (30-day) means of the melt flux (`SHIfwFlx`), the heat flux into the ice, the transfer coefficient and the friction velocity; yearly means of 3D temperature, salinity and velocity; and a monthly text file of global-mean temperature and salinity.
- The build includes the NaN-safe bounds check from project 15 (`experiments/safety_code/mon_solution.F`), so a run that blows up stops instead of passing silently.

A 60-day test (`experiments/baseline_test/`) checked the output files and the sign of the melt flux. It also showed that the melt was falling fast: 0.92 m/yr in the first month, 0.35 m/yr in the second. The baseline run took 43 minutes and wrote 108 MB.

### The unit conversion

The model gives the melt as a fresh-water mass flux `SHIfwFlx` in kg/m^2/s, positive upward (out of the ocean), so melting is negative. To get metres of ice per year:

melt (m of ice per year) = -SHIfwFlx / 917 kg/m^3 x 31 536 000 s

917 kg/m^3 is the ice density, and 31 536 000 s is a 365-day year (this setup has no calendar). For example, -1e-5 kg/m^2/s is 0.344 m of ice per year.

### The notebook, block by block

1. Checks: the run ended normally, there is no NaN in the log, and none of the 65 diagnostic files contain non-finite values.
2. Geometry: maps of the ice draft and the seafloor, and a section showing how open each model cell is (`hFacC`). Every one of the 4851 water columns has exactly one partly open cell, the top one: the smooth ice base is a staircase in the model.
3. Melt over time, from the area-mean of every 30-day mean, saved to `data/melt_series_baseline.csv`.
4. A melt map (mean of the last year), and melt against the depth of the ice base.
5. Temperature, salinity and thermal driving (temperature minus the local freezing point, using the three-equation model's formula) along the cavity.
6. Circulation:
   - the barotropic streamfunction, from adding up the depth-integrated eastward transport from the southern wall northward;
   - the north-south overturning, from adding up the northward transport across each latitude from the top down.

### Results

The area-mean melt rate fell quickly, then slowly:

| after | area-mean melt (m of ice per year) |
| ---: | ---: |
| 1 month | 0.9198 |
| 0.5 years | 0.1205 |
| 1 year | 0.0843 |
| 2 years | 0.0591 |
| 3 years | 0.0472 |
| 4 years | 0.0404 |
| 4.9 years | 0.0369 |

The mean of the last 12 months, 0.0386 m/yr, was 10.5% lower than the 12 months before (0.0431 m/yr). So the melt rate is not steady after 5 years.

Why it cannot become steady in this setup: the box is sealed and nothing adds heat. The water started at -1.9 °C, above its local freezing point under the ice, and every bit of melting uses up some of that heat. The melt rate keeps falling as the heat runs out. The truly steady state of a sealed box is one where melting at depth and refreezing higher up exactly cancel, and the model approaches it slowly.

Where the ice melts and refreezes (last year):

- Melting happens under the sloping part of the ice, south of about 76.5°S. It is strongest at the deep southern end (ice base around 681 m: 0.87 m/yr, zonal mean) and in the south-east corner, where the largest local value is 1.652 m/yr. It gets weaker as the ice base gets shallower: 0.32 m/yr at 594 m, 0.23 at 494 m, 0.11 at 394 m, 0.02 at 294 m.
- Refreezing happens where the slope meets the flat ice, near 74.5-76.5°S, strongest on the western side (down to -0.668 m/yr). At the 200 m flat ice next to the slope, the zonal mean is -0.16 m/yr. Over most of the flat ice further north, the melt rate is close to zero, slightly negative in places. Counting every cell with a negative value, 68.0% of the area is refreezing, but most of that is very weak.
- The area-mean over the whole shelf, 0.0386 m/yr, is small because melting and refreezing nearly cancel.

Temperature and salinity (last year):

- A layer of cold, fresh water lies along the ice base and gets colder and fresher toward the deep southern end, down to -2.43 °C and 34.19 g/kg. This is meltwater mixed with seawater.
- Deep water, especially in the north, is still close to the starting -1.9 °C and 34.4 g/kg.
- Next to the ice, the thermal driving is close to zero: the water touching the ice is close to its freezing point. Deeper down it reaches 0.67 °C. Its lowest value is -0.10 °C, meaning a little water is supercooled, colder than its local freezing point.

Circulation (last year):

- The barotropic streamfunction shows one main circulation cell in the southern half, under the sloping ice, strongest in the south-west near 76°S (up to 544 mSv). Speeds in the yearly mean reach 0.047 m/s.
- The overturning streamfunction is positive everywhere, up to 160 mSv near 500 m depth and 75.8°S. Positive means northward transport above that depth. So water moves south toward the deep grounding line at depth, and north (back up the slope) higher up, along the ice base.

### The ice pump

Putting these together:

1. Water at depth is above its local freezing point, because the freezing point is lower at high pressure. Where it meets the deep ice near the southern wall, it melts the ice.
2. The meltwater is fresh, so the mixture is lighter than the water around it. It rises along the sloping ice base as a plume and flows north, up the slope. The overturning streamfunction shows that flow.
3. As it rises, the pressure drops and the local freezing point rises. The water is cooled toward the freezing point at depth, so higher up it can be colder than its new freezing point: it becomes supercooled.
4. Near the top of the slope and under the shallow flat ice, the supercooled water refreezes onto the ice base (marine ice). The melt map shows refreezing exactly there.

So ice is moved from deep to shallow, driven only by the pressure dependence of the freezing point: melting at depth, rising fresh meltwater, refreezing higher up. In this sealed box, that is the only thing driving the circulation.

### What to look for in the figures

- `figures/02_geometry_maps.png`: the ice draft changes only from south to north (bands), and the seafloor is flat. White is the land wall.
- `figures/02_geometry_section_cells.png`: the staircase. The partly open top cells follow the line of the true ice base.
- `figures/02_melt_timeseries.png`: the fast drop in the first months, then a slow decline that has not levelled off.
- `figures/02_melt_map.png`: red (melting) in the south under the deep ice, strongest in the south-east corner; blue (refreezing) near the top of the slope on the western side; nearly zero under the flat ice.
- `figures/02_melt_vs_draft.png`: melting is strongest at the deepest ice and changes to refreezing at the shallowest. This is the ice pump in one plot.
- `figures/02_ts_sections.png`: the cold, fresh meltwater layer along the ice base. Thermal driving near zero at the ice, positive below.
- `figures/02_barotropic_streamfunction.png`: one circulation cell under the sloping ice.
- `figures/02_overturning.png`: a single overturning cell, south at depth and north along the ice base.

### What I learned

Even with no forcing at all, the pressure dependence of the freezing point drives a circulation and an ice pump: deep ice melts, the meltwater rises along the ice base, and it refreezes higher up. In a sealed box the heat for melting runs out, so the melt rate keeps falling and never becomes truly steady. To study how melt depends on ocean temperature, the far field has to be kept warm, which is what step 3 needs.

## Step 3: Melt rate against ocean temperature (`experiments/rbcs_*`, `experiments/open_*`, `notebooks/03_melt_vs_ocean_temperature.ipynb`)

### The question

How does the melt rate change when the ocean outside the cavity is warmer? I compare four far-field temperatures: -1.9 °C (the baseline), and 0.5, 1 and 2 °C warmer (-1.4, -0.9 and +0.1 °C). Does the melt follow a power law in the thermal forcing, and is the exponent bigger than 1, as Holland, Jenkins and Holland (2008) found?

### Why the setup has to change

Step 2 showed that the default box is sealed and has no heat source, so its melt just keeps falling. To hold the ocean warm, something has to keep adding heat. MITgcm's `pkg/rbcs` does this: in the cells where a mask is non-zero it nudges temperature and salinity toward a target,

tendency = -mask x (T - T_target) / tau

with a time scale tau. I used it for both temperature and salinity. One trap: mask 1 is used for temperature and mask 2 for salinity, and if mask 2 is not set, salinity is quietly not restored (the model only prints a warning). So `data.rbcs` sets both.

To build with `pkg/rbcs`, I added it to `packages.conf` in a new code folder (`experiments/rbcs_code/`, which also has the NaN-safe bounds check) and built `build_rbcs_fast`.

### My first design, and why it failed

The first design (`experiments/rbcs_*`) left the ice where it was and restored temperature and salinity in the northernmost 1° (north of 71°S), under the flat 200 m ice, with tau = 10 days. The mask rose from 0 at 71°S to 1 at the northern wall. Each case started uniformly at its far-field temperature and ran 5 years.

It went wrong in three ways (Part B1 of the notebook):

1. The restoring zone itself melted the ice above it. Restoring kept warm water right against the ice base, so the melt there was 0.550, 2.587, 4.758 and 9.524 m/yr in the four cases, against 0.048, 0.270, 0.553 and 1.234 m/yr in the rest of the cavity. Most of the "extra" melt was made by my forcing, not by the cavity's own circulation.
2. The warm water never reached the deep cavity. In `rbcs_base`, 1189 of the 4851 columns, all south of 74.55°S, had the same melt as the sealed baseline to within 1e-6 m/yr. The restoring did nothing there.
3. The restoring zone did not hold its target: in the warmest case it sat around -0.07 °C against a target of +0.1 °C, because the melting above it kept cooling it. The melt rates were still falling after 5 years.

With these runs, a power-law fit gave an exponent of 1.92, but it mostly measured how my restoring zone melts its own ice lid. So I kept these runs only as a lesson and redesigned the setup.

### The final design: an open-ocean restoring zone

The fix is the same idea ISOMIP+ uses: put the restoring in open ocean, not under ice (`experiments/open_*`).

- Part A of the notebook removes the ice north of 71°S (rows 90-99): the ice draft is set to 0 there, and the ice load too. I checked that the load depends only on the ice draft (the largest difference between cells with the same draft is exactly 0), so a load of 0 is right where the ice is gone. The new files are `experiments/open_geometry/icetopo_open_north.bin` and `phi0surf_open_north.bin`. The ice shelf now ends in a 200 m ice front at 71°S, with 490 columns of open ocean to the north and 4361 ice-covered columns.
- The restoring mask rises linearly from 0 at the ice front to 1 at the northern wall, at every depth, only in the open ocean (`rbcs_mask_open.bin`). There is no ice above it, so it cannot melt anything directly.
- tau = 1 day, so the zone stays very close to its target.
- Each case runs 10 model years from the original cold start (-1.9 °C, 34.4 g/kg, at rest). Only the target temperature changes between cases; the target salinity is always 34.4 g/kg.
- `scripts/run_open_cases.sh` runs the four cases at once, one per core. They took 2 hours 54 minutes each.

Before the long runs, a 30-day test showed that the open zone has exactly zero melt, the sea-surface height stays within 4.5 cm, and the largest CFL number is 0.08.

### Thermal forcing

Thermal forcing (TF) is how far the far-field water is above its freezing point at the depth of the deepest ice base (693.75 m, 701 dbar). With the three-equation freezing formula and S = 34.4 g/kg the freezing point there is -2.421 °C, so TF = 0.521, 1.021, 1.521 and 2.521 °C. The restoring zone sat exactly on its targets in the last year (to 0.001 °C), so the actual thermal forcing equals the nominal one.

### Results

Checks: all runs ended normally, there is no NaN in any log, and every monthly mean is finite. The open zone has zero melt, and its temperature sat on its target to 0.001 °C.

Steadiness. I call a run steady when the mean melt of its last 12 months differs from the 12 months before by less than 2%. After 10 years the two warm cases passed; `open_base` (-6.25% per year) and `open_p05` (-2.17%) did not. So I continued those two for 10 more years from their year-10 pickup files (`experiments/open_base_ext`, `experiments/open_p05_ext`, run with `scripts/run_extension.sh`). The only change in their `data` file is `nIter0 = 175200`, the step to start from.

| case | far field | TF (°C) | years run | melt, last 12 months (m/yr) | change in last year |
| --- | ---: | ---: | ---: | ---: | ---: |
| open_base | -1.9 °C | 0.521 | 20 | 0.0252 | -3.68% |
| open_p05 | -1.4 °C | 1.021 | 20 | 0.1039 | +0.20% |
| open_p10 | -0.9 °C | 1.521 | 10 | 0.2074 | -1.98% |
| open_p20 | +0.1 °C | 2.521 | 10 | 0.4149 | +0.31% |

`open_p05` settled in its second decade. `open_base` is still falling after 20 years, more slowly than before (-6.25% per year at year 10, -3.68% at year 20), so its value is an upper limit for its steady melt. `open_p10` passes the test only just. `open_p20` has a bump between years 2 and 7 (up to about 0.5 m/yr) before it settles, which is when its circulation changes (below).

The power law. Fitting log(melt) against log(TF):

- all four cases: melt = 0.089 x TF^1.79;
- the three steady cases only: melt = 0.104 x TF^1.52.

The exponent between neighbouring cases is 2.11, 1.73 and 1.37: it falls as the ocean gets warmer. The first of these uses the case that is still falling, and will grow if that case keeps falling. So the exponent is between 1 and 2, about 1.5 for the steady cases, and it is not a single number over this range.

### Why the exponent is above 1: more melt pulls in more warm water

With constant transfer coefficients, the melt at each point is proportional to how warm the water touching the ice is. So the exponent can only be above 1 if the water next to the ice warms faster than the far field. What does that is the exchange with the open ocean: meltwater is fresh and buoyant, so more melt drives a stronger circulation across the ice front, which brings in more heat. The southward flow into the cavity across 71°S was 11.0, 79.2, 115.8 and 267.8 mSv, growing like TF^1.98. That is the feedback Holland, Jenkins and Holland (2008) describe, although in their model it also works through the transfer coefficient, which grows with the current speed (step 4 tests that).

### Where the extra melt happens, and a surprise

The melt maps show that the extra melt is not spread evenly:

- The single row of ice right behind the ice front melts fastest: 0.49, 2.31, 4.41 and 8.80 m/yr. It gives 28-32% of all the melt from only 49 of the 4361 ice columns, and its exponent is 1.83.
- The rest of the ice has an exponent of 1.77.
- The deep ice (draft deeper than 300 m, near the grounding line) hardly responds at all: 0.113, 0.139, 0.113 and 0.109 m/yr, an exponent of -0.04. In the two warm cases it melts less than at +0.5 °C.

The deep southern cavity never warms. Below 600 m and south of 76°S the water is -2.00 to -2.08 °C in every case, even when the far field is +0.1 °C. I think this is a density effect. My restored water always has salinity 34.4 g/kg, so warming it makes it lighter. At 700 dbar the restored water is 0.038 and 0.018 kg/m^3 denser than the deep cavity water in the two cold cases, but 0.020 and 0.080 kg/m^3 lighter in the +1 and +2 °C cases. Light water cannot sink under the dense, cold water in the deep cavity, so it flows in near the top, under the flat ice. The overturning shows the change: in the two cold cases there is one ice-pump cell (up to 103 and 110 mSv), while in the warmest case a new cell of the opposite sign (down to -163 mSv) fills the cavity from about 76°S to the ice front, and the ice-pump cell is pushed into the deep south and weakened (47 mSv). The -0.9 °C case is in between, with both cells (-77 and +69 mSv). This regime change also explains why the exponent falls with warming.

In the real Antarctic it is the other way round: the warm water that melts ice shelves (Circumpolar Deep Water) is also saltier, so it is denser and flows down to the grounding line. My uniform-salinity forcing removes that, so these numbers describe "warm, light water at the ice front", not warm deep water. ISOMIP+ (step 7) uses a warm profile in which salinity rises with depth.

### What to look for in the figures

- `figures/03_open_geometry_mask.png`: the ice front at 71°S and the green restoring zone north of it, with no ice above.
- `figures/03_melt_timeseries.png`: all cases start near 1 m/yr (the first-month spike from step 2), then separate. The -1.4 °C case levels off in its second decade; the -1.9 °C case keeps falling slowly; the +0.1 °C case has its bump between years 2 and 7.
- `figures/03_melt_vs_thermal_forcing.png`: on log-log axes the points lie close to a straight line, between the n = 1 and n = 2 lines; the open circle is the case that is still falling, and the two fits (all cases, steady cases) differ mainly because of it.
- `figures/03_overturning_cases.png`: the single red ice-pump cell in the cold case, and the blue cell near the ice front in the warm case.
- `figures/03_melt_maps.png`: in the cold case, melting at the deep south and at the ice front and refreezing on the western side (37% of the area refreezes, mostly weakly); in the warm case, very strong melting in the row behind the ice front (up to 14.6 m/yr in one cell) and much weaker melt elsewhere.

### What I learned

- Where the forcing is applied matters as much as how strong it is. Restoring under ice melts the ice directly and hides the cavity's own response.
- Melt grows faster than linearly with ocean temperature (exponent about 1.5 for the steady cases), because the meltwater drives an exchange flow that brings in more heat.
- A cold cavity can take decades to settle: 20 years were not enough for the coldest case.
- A single area-mean number hides a lot: the extra melt is concentrated at the ice front, while the deep grounding zone barely feels the warming because the warm water is too light to reach it. The density of the warm water, not just its temperature, decides where it goes.

## Step 4: Sensitivity to the melt settings (`experiments/open_g05`, `open_g20`, `open_frict`, `notebooks/04_melt_parameter_sensitivity.ipynb`)

### The question

The three-equation model needs transfer coefficients (how fast heat and salt cross the thin water layer touching the ice), and nobody knows their values well. How much does the melt change if I halve or double them, or if I let them depend on the current speed next to the ice?

### What I changed

All three cases are copies of `open_p10` from step 3 (open-ocean restoring zone, far field -0.9 °C, 10 years from the cold start). I chose `open_p10` as the reference because it was steady in step 3, while the colder cases were still drifting, and because the step 6 resolution test uses the same far field. (I first started these runs as copies of `open_base`, at -1.9 °C, then stopped them after 8 minutes and restarted them at -0.9 °C when step 3 showed that `open_base` was not steady.)

Only `data.shelfice` changes:

- `open_g05`: `SHELFICEheatTransCoeff = 0.5E-4` (gamma_T halved). The salt coefficient follows, because `SHELFICEsaltToHeatRatio` (5.05e-3) is unchanged.
- `open_g20`: `SHELFICEheatTransCoeff = 2.0E-4` (gamma_T doubled).
- `open_frict`: `SHELFICEuseGammaFrict = .TRUE.`, so gamma_T and gamma_S are computed from the friction velocity u* next to the ice with the Holland and Jenkins (1999) formula. Turning this on also changes the default drag on the ice base to `shiCdrag` = 1.5e-3, so I set `SHELFICEDragQuadratic = 2.5E-3` to keep the drag on the flow the same as in the other runs.

The three runs ran at once with `bash scripts/run_open_cases.sh ~/mitgcm_runs/isomip/inputs_rbcs g05 g20 frict` (2 hours 56 minutes).

### Results

All runs ended normally with 121 finite monthly means. Like the reference (-1.98% in the last year), each was still drifting down by about 2% per year (-2.01%, -2.13% and -2.95%), so I compare them at the same time (the last 12 months of year 10).

| case | setting | melt (m/yr) | / reference |
| --- | --- | ---: | ---: |
| open_g05 | gamma_T = 0.5e-4 m/s | 0.1797 | 0.866 |
| open_p10 | gamma_T = 1.0e-4 m/s (default) | 0.2074 | 1 |
| open_g20 | gamma_T = 2.0e-4 m/s | 0.2336 | 1.126 |
| open_frict | velocity-dependent | 0.0587 | 0.283 |

Halving or doubling gamma_T changes the melt by only about 13%. If melt = c x gamma^k, then k = 0.21 from halved to default and 0.17 from default to doubled (0.19 for a fit through all three). So the melt is far from proportional to gamma (k = 1).

Why: the water touching the ice adjusts. With a bigger coefficient the ice takes heat out of that water faster, so the water next to the ice gets colder and closer to its freezing point, and that limits the melt. The mean thermal driving in the top wet cell (temperature minus the freezing point at that cell's salinity) was +0.005, -0.017 and -0.030 °C for halved, default and doubled gamma_T. In the end the melt is limited more by how fast the circulation brings heat into the cavity than by how fast heat crosses the boundary layer.

The velocity-dependent coefficients cut the melt to 28% of the reference. The currents next to the ice are very weak in this cavity, and the model's mean u* was 6.9e-5 m/s, so gamma_T averaged 9.9e-7 m/s, about 100 times smaller than the constant default (range 5.9e-7 to 4.9e-6 m/s; nowhere above 1e-4). Here the melt does depend on gamma, because gamma is so small that the boundary layer really is the bottleneck: the water next to the ice stays 0.29 °C above freezing on average.

The code computes u* = sqrt(shiCdrag x max(1e-6, u^2 + v^2)), so currents are never counted as slower than 1 mm/s. With shiCdrag = 1.5e-3, the lowest possible u* is 3.87e-5 m/s, and 25.6% of the ice area was at that floor in every one of the last 12 months. So the frict result partly depends on this floor, which is a choice in the code, not physics. ISOMIP+ (step 7) uses a similar idea with a fixed "tidal" speed of 1 cm/s.

The circulation hardly changed with the constant coefficients (largest overturning 75.4, 77.1 and 77.1 mSv); with velocity-dependent coefficients it was 99.0 mSv.

### What to look for in the figures

- `figures/04_melt_timeseries.png`: the three constant-gamma curves start far apart (0.64 to 1.32 m/yr in the first month) and come close together within months; the velocity-dependent curve is much lower and starts low, because the water is still at rest at first.
- `figures/04_frict_gamma_ustar.png`: u* and gamma_T are largest along the western side, where the strongest currents run; they are smallest in large areas that sit at the floor value.
- `figures/04_melt_maps.png`: the same pattern for the three constant cases, with the ice-front strip getting stronger with gamma; the velocity-dependent case has almost no melt anywhere on this colour scale (largest 0.545 m/yr).

### What I learned

- With constant coefficients the answer is not very sensitive to their exact value (a factor of 4 in gamma changes the melt by 30%), because the ocean, not the boundary layer, limits the heat supply.
- Whether the coefficients depend on the current speed matters much more (a factor of 3.5 in melt here), and in a weakly moving cavity that answer depends on how the model treats very slow currents. This is why ISOMIP+ tunes its coefficient to a target melt rate (step 7).

## Step 5: Checks and summary (`notebooks/05_checks_and_summary.ipynb`)

### The question

Does the model's bookkeeping add up (is the fresh water from melting really what changes the salt)? And what do all the experiments say together, and about the melt rule I assumed in project 14?

### The salt check

I used the sealed step 2 baseline, because nothing but the ice can change its salt: no restoring and no open boundary. In this setup meltwater is not added as volume; the model takes salt out instead (a "virtual salt flux"). At the ice base the three-equation model balances salt, so

salt removed per second = SHIfwFlx x S_b / rhoConst

where S_b is the salinity right at the ice base. S_b is not an output, but the same balance gives it: gamma_S (S_b - S) = SHIfwFlx x S_b / rhoConst, so S_b = gamma_S S / (gamma_S - SHIfwFlx / rhoConst), with S the salinity of the top wet cell. I added this up for every 30-day mean and compared it with the model's own volume-mean salinity (`ts_global_stats`).

- After 1800 days, the model's volume-mean salinity had dropped by 1.8765e-2 g/kg, and the prediction from the melt flux was 1.8780e-2 g/kg: a difference of 0.08%.
- With the top-cell salinity S instead of S_b, the prediction would be 3.71% too big. My first version of this check made exactly that mistake. I had assumed S_b was within a few hundredths of S, but the meltwater freshens the ice base a lot: in the first month S - S_b averaged 1.80 g/kg, and 0.10 g/kg by the end.
- The melting over the 1800 days added 1.72e11 m^3 of fresh water, 0.37 m spread over the ice base.

### The summary table

`data/summary_table.csv`, every melt rate compared with the steady `open_p10` (far field -0.9 °C, default settings):

| experiment | what changed | TF (°C) | melt (m of ice/yr) | / open_p10 |
| --- | --- | ---: | ---: | ---: |
| baseline | sealed box, no restoring | - | 0.0386 | - |
| open_base | far field -1.9 °C (still falling) | 0.52 | 0.0252 | 0.12 |
| open_p05 | far field -1.4 °C | 1.02 | 0.1039 | 0.50 |
| open_p10 | far field -0.9 °C | 1.52 | 0.2074 | 1 |
| open_p20 | far field +0.1 °C | 2.52 | 0.4149 | 2.00 |
| open_g05 | gamma_T halved | 1.52 | 0.1797 | 0.87 |
| open_g20 | gamma_T doubled | 1.52 | 0.2336 | 1.13 |
| open_frict | velocity-dependent gamma | 1.52 | 0.0587 | 0.28 |
| fine_p10 | half the grid spacing (year 2) | 1.52 | 0.2594 | 1.09 (vs coarse, same years) |

ISOMIP+ Ocean0 (step 7) is a different geometry and unit, so it is listed separately: Gamma_T = 0.0483 gives 29.53 m/yr of water over the deep ice.

So in this setup, 1 °C of ocean warming (from -0.9 to +0.1 °C) doubles the melt, while doubling the transfer coefficient adds only 13%. The choice between constant and velocity-dependent coefficients matters more than either (a factor of 3.5).

### Relation to the project 14 melt rule

In project 14 my toy model used melt = 10 m/yr per °C^2 x T^2, where T is the shelf water's temperature above freezing.

- The exponent: 1.52 for the steady cases here (1.79 with the drifting cold case), not 2. From -1.4 to +0.1 °C the thermal forcing grows 2.47 times; the melt grew 3.99 times, where T^2 would give 6.09 times and T^1 2.47 times.
- The coefficient: melt / TF^2 is 0.065 to 0.100 m/yr per °C^2 in the isomip box, about 100 times smaller than my project 14 value. For ISOMIP+ Ocean0 it is 2.70 for the deep ice (32.2 m of ice per year at a thermal forcing of 3.45 °C) and 0.99 for the whole shelf. These are lower limits for a rule in terms of shelf-water temperature, because I divide by the far-field thermal forcing, which is larger.

So the project 14 rule had the right shape (faster than linear) but an exponent that is a little too large, and its coefficient cannot be taken from one model: it changes by a factor of 30 to 40 between these two setups, depending on the geometry and on whether the warm water can reach the deep ice. My value of 10 is higher than even the ISOMIP+ deep-ice value.

### What to look for in the figures

- `figures/05_salt_budget.png`: the model's salinity change (black) and the prediction with S_b (red dashed) lie on top of each other; the prediction with the top-cell salinity (blue dotted) drifts away.
- `figures/05_summary.png`: left, all the melt rates against thermal forcing, with the step 4 cases stacked at TF = 1.52 °C; right, the size of each change on one log scale.

### What I learned

The model conserves salt to 0.08% once the check uses the right salinity. In this setup the ocean temperature is the biggest control on melt, the form of the transfer coefficient is second, and its exact value and the grid spacing matter much less. A melt rule taken from one model's area-mean can be off by orders of magnitude in another geometry.

## Step 6 (optional): Resolution test (`experiments/fine_code/`, `experiments/fine_p10/`, `notebooks/06_resolution_test.ipynb`)

### The question

Does the melt rate depend on the grid spacing? The sloping ice base is a staircase of cells (step 2), and the strongest melting is in a narrow strip behind the ice front (step 3). Both should be better resolved on a finer grid.

### What I changed

I repeated `open_p10` (far field -0.9 °C) on a grid with half the spacing in both directions: 0.15° x 0.05°, so 100 x 200 points instead of 50 x 100. Everything else is unchanged: 30 levels of 30 m, the 1800 s time step, the open-ocean restoring zone, the three-equation melt model with constant transfer coefficients, and the cold start.

- `experiments/fine_code/SIZE.h`: the grid size is set when the model is compiled, so the fine grid needs its own build (`build_fine_fast`), with tiles of 50 x 50 points.
- `experiments/fine_p10/input/data`: `delX = 100*0.15`, `delY = 200*0.05`, the new bathymetry file, and 2 years (`nTimeSteps = 35040`) instead of 10, because the fine grid is 4 times bigger and so about 4 times slower.
- I kept the time step. The horizontal viscosity (600 m^2/s) is the term most likely to become unstable on smaller cells; its stability number Ah x dt x (1/dx^2 + 1/dy^2) is 0.163 on the smallest fine cells (2.90 km x 5.56 km near 80°S), below the limit of 0.25.
- Part A of the notebook makes the fine input files (written outside the repo, to `~/mitgcm_runs/isomip/inputs_fine/`):
  - the ice draft from the rule behind the original file: it rises 125 m per degree from -700 m at 80°S and stops at -200 m. My first version used `np.maximum` instead of `np.minimum`, and the check against the original file caught it with a 744 m mismatch. The corrected rule matches the original to 1.8e-12 m;
  - the ice load, interpolated from the original draft-to-load table. It matches exactly at the original drafts. The deepest fine row (-690.6 m) is slightly deeper than any coarse cell, so its load is extrapolated in a straight line from the last two table points;
  - the restoring mask and targets, with the ice removed north of 71°S, as in step 3.
- `scripts/run_fine_case.sh` runs it. A 2-day test was stable, and the first month of the real run had a largest CFL number of 0.08. The run took 2 hours 12 minutes.

### Results

The fine grid has 17721 ice-covered columns, not 4 x 4361 = 17444, and 1.4% more ice-covered area than the coarse grid. The reason is the walls: as in the original setup they are the first row and first column, which on the fine grid is half as wide.

| year | coarse grid (m/yr) | fine grid (m/yr) | fine / coarse |
| ---: | ---: | ---: | ---: |
| 1 | 0.3472 | 0.3684 | 1.061 |
| 2 | 0.2388 | 0.2594 | 1.086 |

So the fine grid melts 6-9% more, and the difference grew a little from year 1 to year 2. Where it comes from (year 2):

- the strip within 0.1° of the ice front: 4.637 m/yr on the coarse grid, 5.549 m/yr on the fine grid (20% more). It gives 28.2% and 31.0% of the total melt;
- the rest of the ice: 0.1739 against 0.1817 m/yr (4.5% more);
- the deep ice (draft deeper than 300 m): 0.2666 against 0.2708 m/yr (1.6% more).

So the cavity interior is close to converged at this resolution, but the melt right behind the ice front is not: the finer grid puts more of the incoming warm water against the ice there. The largest single-cell melt rose from 6.4 to 8.8 m/yr.

Both runs are still adjusting after 2 years (the coarse run took about 10 years to settle), so this compares the same stage of spin-up, not two steady states.

### What to look for in the figures

- `figures/06_resolution_melt.png`: the two curves have the same shape; the fine one sits a little higher after the first months.
- `figures/06_resolution_maps.png`: the same pattern on both grids; the fine grid's ice-front strip is thinner and more intense, and the zonal means only separate at the ice front.

### What I learned

At this resolution, the melt in the cavity interior hardly depends on the grid, but the melt at the ice front does. Since the front strip produces almost a third of the total melt in the warm cases, the area-mean melt rates of step 3 would probably be somewhat higher (by about 10% at this stage of the runs) on a finer grid, and the strongest dependence on the grid is exactly where the warm water first touches the ice.

## Step 7 (optional): ISOMIP+ Ocean0 at reduced resolution (`experiments/isomip_plus_code/`, `experiments/isoplus_*`, `notebooks/07_isomip_plus_ocean0.ipynb`)

### The question

Can I set up a standard ice-shelf test case from its published specification, follow it as closely as my laptop allows, and reach its target melt rate? ISOMIP+ Ocean0 (Asay-Davis et al. 2016) has a realistic ice-shelf shape from an ice-sheet model, a WARM far field in which the water gets warmer AND saltier with depth (so the warm water is dense), and a melt formula whose transfer coefficient Gamma_T is tuned until the mean melt over the deep part of the shelf (draft deeper than 300 m) is 30 +/- 2 m/yr over the last 6 months.

### What I used

- The paper (Sections 2.1, 3.1 and 3.2.1, Tables 1, 4 and 6). I downloaded the PDF and pulled the text out with a small Python script (zlib only), because the usual PDF tool (poppler) is not installed.
- The official geometry file `Ocean1_input_geom_v1.01.nc` (1 km grid, 480 x 80 points), downloaded from the GFZ archive (doi:10.5880/PIK.2016.002). Ocean0 and Ocean1 share it.

### What I set up

The setup is written from scratch, not as changes to `isomip`, because almost everything differs.

Compile time (`experiments/isomip_plus_code/`, build `build_isoplus_fast`):

- `SIZE.h`: 121 x 21 points in one tile, 36 levels. The ocean is 320-800 km by 0-80 km at 4 km spacing (120 x 20), plus one land column and one land row, because MITgcm's grid wraps around.
- `packages.conf`: the isomip list plus `rbcs`, and the NaN-safe `mon_solution.F`.
- A copy of `shelfice_thermodynamics.F` with every change commented:
  - the ISOMIP+ liquidus (-0.0573 °C per g/kg, +0.0832 °C, -7.53e-4 °C per dbar);
  - no conversion from potential to in-situ temperature (that conversion is for a nonlinear equation of state; ISOMIP+ uses a linear one);
  - gamma_T = Gamma_T x u* and gamma_S = Gamma_S x u*, with u*^2 = Cd (u^2 + u_tidal^2), Cd = 2.5e-3 and u_tidal = 0.01 m/s (Eqs. 24-27). MITgcm's own velocity-dependent option uses the Holland and Jenkins (1999) form instead.
  - I checked that the build really compiled my copy, not the original.

Run time (`experiments/isoplus_common/` and one `data.shelfice` per run):

- Cartesian grid, f-plane at 75°S (f0 = -1.4087e-4 1/s), 36 layers of 20 m to 720 m.
- Linear equation of state: rho = 1027.51 (1 - 3.733e-5 (T + 1) + 7.843e-4 (S - 34.2)); rhoConst = 1028.
- Mixing: viscAh 6, diffKh 1, viscAz 1e-3, diffKz 5e-5 m^2/s, convective mixing 0.1 m^2/s; no-slip side walls; quadratic bottom and top drag 2.5e-3.
- Initial state and restoring target: the WARM profile, -1.9 °C and 33.8 g/kg at the surface to 1.0 °C and 34.7 g/kg at 720 m. Restoring between 790 and 800 km, rising to 10 per day (`rbcs`, tau = 0.1 day).
- Melt: three-equation model with the boundary layer option, Gamma_S = Gamma_T / 35, perfectly insulating ice (`SHELFICEkappa = 0`).
- Part A of the notebook makes the input files from the geometry file: averages onto 4 km cells (the ice draft only over floating 1 km cells), land where more than half of a cell is grounded, ice thinner than 100 m calved (81 cells), and floating cells with less than 40 m of water made land (50 cells). This leaves 665 ice-covered cells, 157 of them with draft deeper than 300 m; the draft ranges from -565 to -90 m. The ice load is g times the integral of (rho - rhoConst) down to the ice base with the WARM profile: -2976 to -652 Pa.

### Problems on the way

- The analytic bedrock formula in the paper (MISMIP+ Eq. 1) does not match the bedrock in the official file everywhere: for x < 640 km the median difference is 0.23 m and 95% of points are within 8.1 m, but the largest difference is 90 m, and beyond 640 km 104 m, where the file sits at the -720 m floor. I used the file's bedrock and list this as a deviation.
- The model refused to start: "Cannot mix z, p and r". I had given the layer thicknesses as `delR` together with other z-coordinate settings; `delZ` fixed it.
- The first test blew up to NaN within 10 days. The vertical CFL number in thin partial cells under the ice (`advcfl_W_hf`) reached 1.38. I raised `hFacMin` from 0.1 to 0.2 (no cell less than 20% open) and halved the time step from 600 s to 300 s. After that the runs were stable.
- A comment I wrote in `data` contained the word "NaN". The model copies its input files into its log, so my NaN check would have failed a good run. I changed the comment and made `run_experiment.sh` ignore the copied input lines.

### Calibration

Four 1-year runs with Gamma_T = 0.025, 0.05, 0.1 and 0.2 (`scripts/run_isoplus_cases.sh`, about 65 minutes each, four at once). Mean melt where the draft is deeper than 300 m, last six 30-day means, in m/yr of water:

| Gamma_T | melt |
| ---: | ---: |
| 0.025 | 15.72 |
| 0.05 | 31.06 |
| 0.1 | 50.99 |
| 0.2 | 67.46 |

The curve bends: between neighbouring runs the exponent falls from 0.98 to 0.72 to 0.40. At small Gamma_T the thin boundary layer is the bottleneck, so the melt is almost proportional to Gamma_T; at large Gamma_T the ocean cannot bring in heat fast enough, as in step 4. So I interpolated (in log-log) between the two runs that bracket 30 m/yr, which gives Gamma_T = 0.0483. A single fit through all four runs would have given 0.0545.

### The final run

Gamma_T = 0.0483 (`experiments/isoplus_final/`), 2 years (210380 steps of 300 s), 1 hour 32 minutes.

- It ended normally with no NaN. The largest vertical CFL number in partial cells was 0.49, the largest speed 0.71 m/s, and the sea surface moved by at most 4.5 cm.
- Mean melt where the draft is deeper than 300 m, last 6 months: 29.53 m/yr, inside the 30 +/- 2 m/yr target, and the same as in the 6 months before (29.53), so the run is steady.
- Mean over the whole ice shelf: 10.79 m/yr. The largest local melt is 68.3 m/yr, near the deep grounding line (small x, where the ice is thickest); almost nowhere refreezes (lowest -0.15 m/yr).

### Comparison with the literature

ISOMIP+ fixes the deep melt rate, so the model-dependent number is the Gamma_T needed to reach it. The paper suggests Gamma_T = 2.2e-2 as a first guess, and reports that POP2x needs about 0.11 for Ocean0. My value, 0.0483, is 2.2 times the first guess (which would give about 13.9 m/yr in my model) and 0.44 times the POP2x value. So my model melts more easily than POP2x for the same Gamma_T. Possible reasons, which I have not tested: my 4 km grid and boundary-layer averaging, and the partial-cell settings near the ice, all change the velocity and temperature that the melt formula sees.

### Deviations from the specification

The full list is at the end of notebook 07 (13 items). The main ones: 4 km instead of 2 km resolution; the ISOMIP+ melt formula through a changed copy of the MITgcm code; the file's bedrock instead of Eq. 1; `hFacMin = 0.2` and a 300 s step for stability; a flux-limited advection scheme; a virtual salt flux instead of a volume flux; only two cells in the restoring zone (mask 0.4 and 0.8); 30-day averaging periods; and 1-year calibration runs.

### What to look for in the figures

- `figures/07_isoplus_geometry.png`: the ice shelf between about 460 and 640 km, deepest at its grounding line (small x), with land around it; the restoring zone at 790 km.
- `figures/07_isoplus_calibration.png`: the four calibration runs on log-log axes, flattening at large Gamma_T, the target band, and the final run on it.
- `figures/07_isoplus_final.png`: the deep melt starts at 33-35 m/yr, enters the target band in month 4, and from month 8 on stays between 29.3 and 29.8 m/yr; the melt map shows the strongest melting along the deep grounding line at small x.

### What I learned

- A published specification still leaves many choices to the modeller (minimum water column, partial cells, boundary layer, how to handle the restoring zone at coarse resolution), and each is a deviation to document.
- The melt formula is not universal: the same target melt needs a different Gamma_T in different models (0.0483 here, 0.11 in POP2x), which is exactly why ISOMIP+ calibrates.
- With warm water that is also salty, the warm water reaches the deep grounding zone, and that is where the melting is strongest: the opposite of my step 3 runs with uniform salinity.

## What I would do next

- Couple the ocean to an ice-sheet model. In every run here the ice shape is fixed, so melting never thins the ice, the cavity never grows, and the grounding line never moves. In a coupled model (as in the ISOMIP+ Ocean3/Ocean4 and MISOMIP IceOcean experiments) the melt thins the shelf, the ice upstream speeds up, and the new cavity shape changes the circulation and the melt again. MITgcm can do this with the `shelfice` and `streamice` packages, or through coupling to an external ice-sheet model.
- Look for tipping points in the coupled ice-ocean system. Project 14 had a toy threshold model of an ice-shelf margin. With a real cavity, the question becomes whether there is a far-field temperature at which warm water suddenly gets access to the deep grounding zone (for example when it becomes denser than the cavity water, which my step 3 runs show matters), and whether, once the ice retreats onto a bed that deepens inland, the retreat keeps going on its own (the marine ice-sheet instability). That needs the coupled model above and runs that slowly raise and then lower the ocean temperature, to look for hysteresis as in project 14.
- Repeat step 3 with warm water that is also salty (a two-layer far field like the ISOMIP+ WARM profile), so that the warm water is dense and can reach the grounding line, and compare the melt exponent with the one I found.
- Finish ISOMIP+ properly: run Ocean0 at the 2 km common resolution and Ocean1 (cold start, warm forcing, 20 years), and compare with the published model results.
