# Tutorial: freshwater hosing in MITgcm

These notes walk through each step of the project in order. I wrote them so I can explain every command and every choice.

The model source code is in `~/models/MITgcm` and every run is in `~/mitgcm_runs/hosing/`, both outside the git repo. The repo only holds the scripts, the changed parameter files, the notebooks, the figures and small processed results.

## Step 1: Build and verify (`scripts/build.sh`, `scripts/run_experiment.sh`, `notebooks/01_build_and_verify.ipynb`)

### The question

Does my build of MITgcm on this Mac give the same answer as the official reference run? I need to know this before I change anything. Otherwise I could not tell a real effect of my changes from a problem with my compiler or setup.

### What I did

I made no changes to the model in this step. I cloned MITgcm at the release tag `checkpoint69q` (commit 853761d8), built the `tutorial_global_oce_latlon` experiment, and ran its short test exactly as provided.

The setup is a global ocean on a 4 x 4 degree latitude-longitude grid, 90 x 40 points from 80°S to 80°N, with 15 levels. The levels get thicker with depth, from 50 m at the surface to 690 m at the bottom. It is forced by monthly winds, heat flux and freshwater flux, plus relaxation of surface temperature and salinity toward observations.

### The main commands

1. `uname -m` prints `arm64`, which means an Apple-silicon processor. MITgcm's build-options files are named by system and processor, so I use `darwin_arm64_gfortran` (darwin = macOS, gfortran = the compiler). `darwin_amd64_gfortran` would be for an Intel Mac.
2. `scripts/build.sh` runs three commands in an empty build folder:
   - `genmake2` reads the experiment's `code/` folder (grid size in `SIZE.h`, package list in `packages.conf`) and the build-options file, and writes a Makefile.
   - `make depend` works out which source files include which headers.
   - `make` compiles everything into one executable, `mitgcmuv`.
3. My first build failed at the very end with `ld: library 'netcdf' not found`. The tutorial compiles `pkg/mnc` (netCDF output), so the linker needs the netCDF C library. Homebrew's `nf-config` only reports the folder of the netCDF-Fortran library, but the C library lives in `/opt/homebrew/lib`. I added that folder to `LIBRARY_PATH` in `build.sh`, which gfortran searches when linking. I did not change any MITgcm file.
4. I built twice:
   - an "IEEE" build (`-ieee`, no optimization, `-O0`) for the closest possible comparison with the reference;
   - a "fast" build with the optimized flags from the build-options file (`-O3 -ftree-vectorize -funroll-loops`), for the long runs.
5. `scripts/run_experiment.sh` makes a run folder, links the tutorial's original input files into it, copies the executable in, and runs it. If I give it a folder of changed files, those are copied over the originals. It refuses to run in a folder that already has output, so I cannot overwrite a run by accident.

### Comparing with the reference

The model prints a block of `%MON` (monitor) lines at every monitor time: statistics such as the mean, maximum and standard deviation of temperature, salinity and velocity, and the CFL numbers. The experiment's `results/output.txt` holds the same lines from the reference run. In the notebook I read all 165 monitored quantities at all 21 monitor times. For each quantity I find the smallest number of significant digits that agree, the same idea as MITgcm's own `testreport`.

- IEEE build: 125 of 165 quantities are identical in every printed digit. The worst of the rest agrees to 11.3 digits.
- Fast build: 117 of 165 are identical. The worst of the rest agrees to 10.0 digits.
- One quantity seems not to match at all: `dynstat_wvel_mean`, the global mean vertical velocity. Its values are around 1e-22 m/s, and they should be exactly zero. They are rounding noise, so comparing their relative size means nothing, and I leave this one out.

Matching to 10-16 digits means my build solves the same equations in the same way as the reference. The small differences come from the order of floating-point operations, which depends on the compiler and the optimization level. In an ocean model, tiny differences like these can grow over time, so a long run will not stay digit-for-digit identical to another machine's run. But the climate of the run, its averages and circulation, should be the same.

### The main parameters in the `data` file

All times are in seconds.

- `deltaTmom = 1800` (30 minutes): the time step for the momentum equations (velocity).
- `deltaTtracer = 86400` (1 day): the time step for the tracers, temperature and salinity.
- `deltaTClock = 86400`: how far the model clock moves each step. It is set to the tracer step, so one step means one day for T and S, and a model year is 360 steps.
- `deltaTfreesurf = 86400`: the time step for the free surface (sea level).
- `nIter0 = 0`: the starting step number. 0 means start from the initial state (Levitus temperature and salinity, at rest). A nonzero value means start from a pickup (restart) file saved at that step.
- `nTimeSteps = 20`: how many steps to run. The test only runs 20 days. A comment in the file says about 100 years gives "a reasonable flow field".
- `pChkptFreq = 1728000` (20 days): how often to write a permanent pickup file, which I can restart from later.
- `dumpFreq = 864000` (10 days): how often to write snapshots of the model state.
- `monitorFreq = 1`: how often to print the monitor block. A value smaller than one step means every step.
- `tauThetaClimRelax = 5184000` (60 days) and `tauSaltClimRelax = 15552000` (180 days): the time scales for relaxing surface temperature and salinity toward the climatology files `lev_sst.bin` and `lev_sss.bin`. The model adds a flux to the top layer that pulls it back toward the observed value. With a shorter time scale the pull is stronger. Temperature is held about three times more tightly than salinity.
- `externForcingPeriod = 2592000` (30 days) and `externForcingCycle = 31104000` (360 days): the forcing files hold 12 monthly fields, and the model cycles through them every year.

Why is the tracer time step 48 times longer than the momentum step? The deep ocean takes thousands of years to adjust, mostly because temperature and salinity are slow to spread through it. The fastest things in the model are in the velocity field (waves and currents), and they limit the momentum time step. With "asynchronous" time stepping, both are advanced once per model step, but tracers advance a full day and velocity only 30 minutes. In effect the velocity field responds as if it had 48 times more inertia. A steady state does not depend on how fast the model gets there, so the final equilibrium is still correct, and it is reached far sooner. The cost is that anything time-dependent is distorted: the speed of adjustment, how long a transient lasts, the timing of a response. So in this project I take the equilibrium states literally, but not the time scales of change, such as how many years the AMOC takes to weaken or recover.

### CFL numbers

The monitor block includes CFL numbers: how far the flow moves in one time step, as a fraction of a grid cell (speed x time step / cell size). If a CFL number gets near 1, water moves more than one cell per step and the explicit scheme can become unstable.

- `advcfl_uvel_max`, `advcfl_vvel_max`, `advcfl_wvel_max` and `advcfl_W_hf_max` are computed from the current velocity. I checked the source (`pkg/monitor/monitor.F`): they use the larger of the momentum and tracer steps, here 1 day, so they are a cautious check.
- `trAdv_CFL_u_max`, `_v_max` and `_w_max` come from the latest tracer advection step.
- In the 20-step test the largest values are about 0.05 to 0.09 (`advcfl_W_hf_max` = 0.0911), far below 1. The test starts from rest, though, so the currents are still weak. I will check the CFL numbers again once the circulation is fully spun up.

### Timing

I timed a 2-year run (720 steps) with output reduced to roughly what the long runs will use: monitor once a model year, no snapshots, no pickups. The changed `data` file is in `experiments/timing_2yr/`.

- Fast build: 4.54 s per model year, including start-up. That is about 1.3 hours for 1000 model years and 2.5 hours for 2000.
- IEEE build: 49.45 s per model year, 10.9 times slower. So I will use the fast build for the long runs. Its agreement with the reference (at least 10 digits) is good enough for that.

### What to look for in the figure

- `figures/01_reference_match.png`: each point is one monitored quantity, sorted from worst to best agreement. Most sit at 16 (identical). The rest agree to between 10 and 14 digits. The IEEE build (blue) agrees slightly better than the fast build (orange), as expected.

### What I learned

My build reproduces the reference run to at least 10 significant digits. So any difference I see later comes from my changes, not from the build. The accelerated time stepping lets the deep ocean reach equilibrium in a few hours on a laptop, but it means I should only trust equilibrium states, not the speed of changes. The only build problem was a library search path, which I fixed in the build script without touching MITgcm.

## Step 2: Spin-up (`experiments/spinup/`, `notebooks/02_spinup.ipynb`)

### The question

Under the tutorial's forcing and surface restoring, does the Atlantic overturning circulation (AMOC) settle to a steady state? What does it look like? Every later experiment starts from this state, so I need to know it well, including how much it is still drifting.

### How the surface forcing is supplied

Before changing anything I read the input files and the source code. All forcing files hold 12 monthly 90 x 40 fields (single precision, big-endian). The model interpolates between months and repeats the cycle every 360-day year.

- Wind stress: `trenberth_taux.bin` and `trenberth_tauy.bin`, in N/m^2. Positive taux pushes eastward. Values range from -0.45 to +0.41 N/m^2.
- Heat flux: `ncep_qnet.bin`, in W/m^2. Positive means heat leaving the ocean (cooling). Values range from -225 to +442 W/m^2, with a global mean of 0.00.
- Freshwater flux: `ncep_emp.bin`, evaporation minus precipitation minus runoff. The file is in m/s (a volume flux), and the model multiplies it by `rhoConstFresh` = 1000 kg/m^3 as it reads it (`external_fields_load.F`), so inside the model EmPmR is in kg/m^2/s. Positive means net evaporation, which makes the ocean saltier. Values range from -1.8e-7 to +0.9e-7 m/s, and the global integral is 0.0000 Sv. So adding freshwater means a negative EmPmR anomaly.
- Restoring: the top layer (50 m) is pulled toward `lev_sst.bin` over 60 days and toward `lev_sss.bin` over 180 days. The salinity restoring flux is -(S_top - SSS_obs) x 50 m / 180 days (`forcing_surf_relax.F`).
- How freshwater changes salinity: the tutorial has a linear free surface and `useRealFreshWaterFlux = .TRUE.`, which sets `convertFW2Salt = -1`. EmPmR does not change the thickness of the tracer cells. Instead the model turns it into a "virtual salt flux", EmPmR times the local surface salinity.

### What I changed and why

- I added `pkg/diagnostics`, the standard way to save time averages. This is a compile-time change. `experiments/spinup/code/` holds only `packages.conf` (the tutorial's list plus `diagnostics`) and `DIAGNOSTICS_SIZE.h` (`numDiags` raised from 15 to 150, so all the output fits). `build.sh` now copies the tutorial's code folder, puts my changed files on top, and builds from that combined folder, so it is clear which file wins.
- `experiments/spinup/input/data`:
  - 720000 steps (2000 years);
  - one global output file instead of one per tile (`globalFiles`);
  - pickups every 100 years;
  - no 10-day snapshots;
  - monitor once a year.

  The physics is unchanged.
- `data.pkg` turns diagnostics on.
- `data.diagnostics` lists the output:
  - every year: annual-mean 3D meridional velocity (for the AMOC); annual-mean surface temperature and salinity; the surface fluxes `SRELAX`, `TRELAX`, `SFLUX` and `oceFWflx`;
  - every 10 years: 10-year means of 3D temperature, salinity and velocity;
  - every 100 years: a 12-month climatology of `SRELAX`, averaged over those 100 years (`averagingFreq` = 30 days, `repeatCycle` = 12), which is what step 5 needs;
  - every year: a small text file of global-mean temperature and salinity at every level.
- A 2-year test (`experiments/spinup_test/`) showed that every output list was written correctly. It also showed that `GM_PsiY`, the eddy-induced streamfunction from the GM scheme, was zero everywhere. MITgcm only fills it when GM uses its "advective form", and this tutorial uses the "skew-flux" form. I did not change the GM scheme, because that would change the model. So my AMOC is the resolved (Eulerian) overturning only, without the eddy-induced part.
- I installed MITgcmutils from the MITgcm copy I cloned (`pip install --no-deps`), to read the binary output.

The run took 2 h 36 min (4.7 s per model year) and wrote about 1.0 GB, all in `~/mitgcm_runs/hosing/spinup/`.

### The notebook, block by block

1. Grid. Temperature and salinity sit at cell centres (`XC`, `YC`). The meridional velocity v sits on the southern face of each cell, at latitude `YG`, with face width `DXG` and open fraction `hFacS`.
2. Atlantic mask, first try: every ocean point from 34°S northward, between 100°W (260°E) and 20°E. 34°S is about the latitude of the tip of Africa, and 20°E is the usual boundary with the Indian Ocean there.
3. Fixing the mask. A text map of the box showed it also caught Pacific points, because the Pacific coast of the Americas is east of 100°W. I removed everything west of 290°E from 34°S to 10°N (off Chile, Peru, Ecuador and Colombia), and 260-270°E at 14°N (the Pacific side of Central America). That removed 77 points and left 460. The map shows the Gulf of Mexico, the Caribbean, the Nordic Seas, Baffin Bay and Hudson Bay inside the mask. The Mediterranean does not exist at this resolution, and the 20°E cut keeps the Barents Sea out.
4. A v-point counts as Atlantic only if the cells on both sides of it are Atlantic (404 v-points).
5. Streamfunction. At each level and latitude I add up v x face width x level thickness x hFacS across the Atlantic, which is the northward transport in that level. Then I add those up from the surface downward. So psi at depth z is the total northward transport above z, positive for the AMOC (northward above, southward below). My first version subtracted instead of adding, which flipped the sign. The first plot showed the main cell as negative, and that is how I caught it.
6. AMOC index: the maximum of psi between 20°N and 60°N and below 500 m, for every model year. It is saved to `data/amoc_index_spinup.csv`.
7. Drift: from the yearly global statistics I fit straight lines over the last 500 years.
8. Salt budget for the final year (see below).

### Results

- Final state (mean of years 1901-2000): AMOC index 19.80 Sv, at 52°N and 1080 m depth.
- In the last 100 years the AMOC index barely changes from year to year (standard deviation 0.008 Sv). The forcing repeats every year and there is no atmosphere or weather, so there is almost no natural variability.
- AMOC history: it started near 17-18 Sv, dropped to 12.66 Sv by year 9 as the model adjusted from its initial state, then rose quickly and kept creeping up more slowly. Over the last 500 years the trend is +0.042 Sv per century.
- Global mean temperature rose from 3.6171 °C to 4.4727 °C. The trend over the last 500 years is +0.0072 °C per century, so it is still warming slowly.
- Global mean salinity rose to 34.73907 g/kg at year 488, then fell steadily to 34.72521 g/kg. The trend over the last 500 years is -0.00118 g/kg per century.
- CFL numbers in the last monitor block are all below 0.04 (largest: `advcfl_uvel_max` = 0.0396), even with the currents fully spun up. So the time steps are safely stable.
- There are 20 pickup files (every 100 years; the last is `pickup.0000720000`) and 20 SRELAX climatologies.

### The salt budget

Why does salinity fall in a straight line while everything else levels off? In the final year:

- `SFLUX` (the total surface salt flux) is +4013 t/s. That would make salinity rise by 0.0091 g/kg per century.
- `SRELAX` (restoring) is -550 t/s. That would make salinity fall by 0.00125 g/kg per century.
- The measured trend is a fall of 0.00118 g/kg per century, which matches the restoring flux alone.

So the large virtual salt flux from EmPmR (+4563 t/s, the difference between the two) does not change the total salt content. My reading is that with a linear free surface and a real freshwater flux, the model also moves salt through the fixed surface with the surface vertical velocity (w x S), and that cancels the virtual flux in the global total. `SFLUX` leaves that term out. I have not proven this from the code yet. I will test it directly with the hosing runs in step 7. The steady fall in salinity is the restoring slowly removing salt, a drift of about 0.003% per century.

### What to look for in the figures

- `figures/02_atlantic_mask.png`: blue is my Atlantic mask, red is the Pacific points I removed from the first box.
- `figures/02_amoc_index_spinup.png`: the quick drop and recovery in the first decades, then a slow rise that flattens out near 19.8 Sv.
- `figures/02_streamfunction_final.png`: one main clockwise cell. Water flows north in the upper ~1000 m, sinks between about 45°N and 60°N, and returns south between about 1000 and 3000 m. Below about 3000 m there is a weak cell of the opposite sign. The star marks the AMOC index.
- `figures/02_global_drift.png`: temperature still rising slowly. Salinity rising for about 500 years, then falling in a straight line.

### What I learned

After 2000 years the AMOC is close to steady at 19.8 Sv. Its remaining drift (+0.042 Sv per century) is tiny compared with the changes hosing should cause. The deep ocean is still adjusting slowly, which is normal for a model with only 2000 years of spin-up. The salt budget showed that a diagnostic named "total salt flux" does not always match the salt content. With this free-surface setup, only the restoring changes the total salt, and that will matter when I check the hosing runs.

## Step 3: Make the hosing forcing (`notebooks/03_make_hosing_forcing.ipynb`)

### The question

How do I add exactly 0.1, 0.3 or 0.5 Sv of extra freshwater to the subpolar North Atlantic in this model? I need to use the right units and sign, and a file the model reads exactly like its original freshwater forcing.

### What I changed and why

I made three new freshwater-flux files, `experiments/hosing_forcing/emp_hose010.bin`, `emp_hose030.bin` and `emp_hose050.bin`. Each is the original `ncep_emp.bin` plus a constant anomaly over the hosing region. The model can only read one EmPmR file, so each hosing file holds the normal freshwater flux and the extra water together. No model code changes. A run uses a hosing file by setting `EmPmRFile` in its `data` file.

The hosing is not balanced by removing freshwater anywhere else, so the ocean gains water. Some studies compensate it elsewhere. I kept it simple and uncompensated.

### The notebook, block by block

1. Grid and masks. I read the cell areas (`RAC`) and the surface ocean mask from the spin-up. The Atlantic mask comes from `data/atlantic_mask.csv`, which notebook 02 now saves, so both notebooks use exactly the same mask.
2. Region. Every Atlantic ocean cell with its centre between 50°N and 70°N and between 60°W and 20°E: 87 cells, 8.50 million km^2 (cell edges from 48°N to 72°N). The 60°W edge keeps Hudson Bay and Hudson Strait out, and keeps the Labrador Sea, the Irminger Sea, the Iceland Basin and the Nordic Seas in, where the model's deep water forms. The region mask is saved to `data/hosing_region_mask.csv` for later steps.
3. The original file has 43200 values, 12 months x 40 x 90, stored as single-precision big-endian numbers (`>f4`), which is what `readBinaryPrec=32` expects.
4. The anomaly. Volume flux / area = speed, so a flux of F Sv spread over area A gives -F x 1e6 / A m/s at every region cell. The minus sign is because positive EmPmR means water leaving the ocean. The anomaly is the same in all 12 months. I add it to every month of the original field and write the result in the same format.
   - 0.1 Sv: -1.1762e-08 m/s, or 0.366 m of extra water per year over the region.
   - 0.3 Sv: -3.5285e-08 m/s, or 1.097 m per year.
   - 0.5 Sv: -5.8808e-08 m/s, or 1.829 m per year.

   For comparison, the original EmPmR over the ocean ranges from -1.796e-07 to 8.936e-08 m/s.
5. The check. I read each file back from disk, subtract the original, multiply by the model's cell areas over ocean points, and add up. Every month of every file gives 0.10000000, 0.30000000 and 0.50000000 Sv, to 8 decimals. The largest change outside the region is exactly 0. All three files are 172800 bytes, the same as the original.

### What to look for in the figure

- `figures/03_hosing_region.png`: dark blue is the hosing region, light blue is the rest of the Atlantic mask. The model's 4° coastline is coarse, so some cells near Greenland and Iceland are ocean in the model, though not on the real map.

### What I learned

Getting a forcing perturbation right needs three things checked separately: the units the model expects (the file is m/s, and the model converts to kg/m^2/s), the sign convention (positive = evaporation), and the file format (record count, precision, byte order). Integrating the file back from disk tests all three at once.

## Step 4: Hosing with salinity restoring (`scripts/run_hosing_set.sh`, `notebooks/04_hosing_with_restoring.ipynb`)

### The question

How does the AMOC respond to 0.1, 0.3 and 0.5 Sv of extra freshwater in the subpolar North Atlantic, while surface salinity is still restored toward observations? Does it recover when the hosing stops? And how much of the added freshwater does the restoring cancel?

### What I ran

All runs start from the spin-up pickup at step 720000 (end of model year 2000) and use the spin-up's executable, `data.pkg` and `data.diagnostics`.

- `restore_control`: 500 years (2001-2500) with the original forcing.
- `restore_hoseXXX_on`: 200 years (2001-2200) with `EmPmRFile` set to the hosing file from step 3.
- `restore_hoseXXX_off`: 300 years (2201-2500) with the original EmPmR again, restarted from the "on" run's pickup at step 792000 (year 2200). Pickups are written every 36000 steps (100 years), and 792000 = 22 x 36000, so this pickup is always there. The script checks that it exists before starting an "off" run.

Each `experiments/restore_*/input/data` differs from the spin-up's `data` only in `nIter0` (start step), `nTimeSteps` (length) and, for the "on" runs, `EmPmRFile`.

`scripts/run_hosing_set.sh restore experiments/hosing_forcing ~/mitgcm_runs/hosing/spinup` runs the control and the three on-then-off chains at the same time: 4 model runs, one per performance core of the M1. With 4 runs at once, each ran at about 10 s per model year instead of 4.7 s, because they share memory bandwidth, not just cores. The whole set took 1 h 30 min.

Two script changes made this possible:

- `run_experiment.sh` now takes a list of change folders or files, copied in order so that later ones win, and an optional run folder to restart from. It links only the pickup for the run's start step (`nIter0`).
- A 1-year test of the 0.5 Sv run exposed a bug in the script. To read `nIter0` I kept every digit on that line, including the "0" in the name `nIter0`, so 720000 became "0720000", which `printf` reads as an octal number. The pickup check stopped the run before the model started. I fixed it by reading only the number after the "=".

The same 1-year test confirmed that the model really receives the hosing. Over the year, the freshwater flux into the ocean (`oceFWflx`) was 0.5000 Sv higher than in the last spin-up year, all of it inside the hosing region.

### The notebook, block by block

1. Shared functions. The AMOC calculation now lives in `src/amoc_tools.py`, so notebooks 04, 06 and 07 all use exactly the same code. I checked that it reproduces all 2000 AMOC values from notebook 02 exactly.
2. Checks: every run ended normally, has the right number of annual means, and the three year-2200 pickups exist.
3. AMOC index for every year of every run. Each hosing run's "on" and "off" parts are joined into one 500-year series, saved to `data/amoc_index_restore.csv`.
4. Summary table. "End of hosing" is the mean of years 2191-2200 and "end of recovery" the mean of years 2491-2500, each compared with the control over the same years.
5. Restoring cancellation. The extra restoring salt flux (hosing run minus control, `SRELAX`) is turned into the freshwater flux it is equivalent to. A salt flux F_s in g/m^2/s equals removing F_s / (1000 kg/m^3 x S) m/s of fresh water, where S is the local surface salinity in g/kg. That is the same conversion the model uses for EmPmR (`convertFW2Salt = -1`). Summed over an area, this gives Sv. I sum it over the hosing region, and over the whole ocean. The series are saved to `data/restoring_cancellation_restore.csv`.

### Results

The control AMOC is 19.87 Sv over 2191-2200 and 19.94 Sv over 2491-2500.

| hosing | end of hosing (2191-2200) | change vs control | lowest value (year) | end of recovery (2491-2500) | change vs control |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 Sv | 17.60 Sv | -2.27 Sv (-11.4%) | 17.21 Sv (2018) | 19.99 Sv | +0.04 Sv |
| 0.3 Sv | 14.27 Sv | -5.60 Sv (-28.2%) | 13.52 Sv (2020) | 20.06 Sv | +0.11 Sv |
| 0.5 Sv | 11.70 Sv | -8.17 Sv (-41.1%) | 10.77 Sv (2021) | 20.15 Sv | +0.20 Sv |

- The AMOC weakened in every run, more for stronger hosing, but less than in proportion: 2.27, 1.87 and 1.63 Sv of weakening per 0.1 Sv of hosing.
- It reached its lowest value about 20 years after the hosing started, then crept back up slightly while the hosing continued. Because of the accelerated time stepping (step 1), I do not take "20 years" literally.
- It recovered fully in every run once the hosing stopped. It even overshot: the peaks after 2200 were 20.21 Sv (2218), 20.64 Sv (2223) and 21.01 Sv (2244) for 0.1, 0.3 and 0.5 Sv, and at 2491-2500 it was still 0.20 Sv above the control.

How much of the hosing the restoring cancels (freshwater removed by the extra restoring, in Sv):

| hosing | where | year 2001 | mean 2001-2200 | mean 2191-2200 | fraction of hosing, 2191-2200 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 Sv | hosing region | 0.044 | 0.081 | 0.080 | 80% |
| 0.1 Sv | whole ocean | 0.045 | 0.070 | 0.064 | 64% |
| 0.3 Sv | hosing region | 0.132 | 0.250 | 0.250 | 83% |
| 0.3 Sv | whole ocean | 0.136 | 0.223 | 0.204 | 68% |
| 0.5 Sv | hosing region | 0.224 | 0.437 | 0.437 | 87% |
| 0.5 Sv | whole ocean | 0.229 | 0.406 | 0.373 | 75% |

So within a few years, the restoring in the hosing region removes 80-87% of the added freshwater. Over the whole ocean the cancellation is a little smaller late in the hosing period, which means the restoring outside the region is working slightly the other way there.

After the hosing stops, the whole-ocean restoring becomes negative for a few decades (lowest values -0.029, -0.067 and -0.102 Sv around 2217-2220). Restoring is then removing salt, which means the surface is saltier than observed somewhere outside the hosing region. A likely cause of the overshoot is salty water that built up in the subtropical Atlantic while the AMOC was weak and is then carried north. I have not tested this yet. I will look at salinity maps in step 6.

### A correction about salt

From the step 2 budget I first guessed that freshwater fluxes do not change the total salt in this setup, and that the AMOC overshoot came from extra salt added by restoring. The hosing runs show this is wrong. Global mean salinity is lower in every hosing run than in the control, by 0.00487, 0.01288 and 0.01647 g/kg at year 2200 for 0.1, 0.3 and 0.5 Sv.

A rough estimate fits the 0.5 Sv case:

- 0.5 Sv for 200 years with no restoring would lower the mean salinity by about 0.08 g/kg.
- The whole-ocean restoring cancelled 0.406 Sv on average, leaving a net 0.094 Sv.
- That net 0.094 Sv predicts about -0.015 g/kg, close to the measured -0.016.

The explanation that fits both the spin-up and these runs: the model's virtual salt flux (EmPmR x local salinity, +4527 t/s in the last spin-up year when computed directly from `oceFWflx` and surface salinity) is cancelled by salt carried through the fixed sea surface by the surface vertical velocity (w x S), except for the part of the freshwater that goes into changing global sea level. So the total salt only responds to the net freshwater the ocean gains or loses. In the spin-up the net is zero, so only restoring changed the salt. With hosing the net is not zero, so the ocean is diluted. I have not yet confirmed this line by line in the code. Step 7's salt check will test it properly.

### What to look for in the figures

- `figures/04_amoc_restoring.png`: a fast drop in the first years of hosing, a slight recovery while the hosing continues, a fast recovery after 2200, then an overshoot that slowly fades. The three runs keep the same order throughout.
- `figures/04_restoring_cancellation.png`: solid lines are the restoring inside the hosing region, dashed lines the whole ocean, dotted lines the hosing strength. The solid lines climb close to the dotted lines within a few years and stay there. After 2200 they drop to zero, and the dashed lines go negative for a few decades.

### What I learned

With salinity restoring, the model cancels 80-87% of the added freshwater in the hosing region, yet the AMOC still weakens by 11-41%. The restoring flux depends on how far the surface salinity is from the observed value, so the surface has to stay fresher than observed for the restoring to push back. That remaining freshening is what weakens the AMOC. But it can never run away, because a fresher surface means a stronger push back. So the AMOC weakens in proportion to the hosing, and recovers completely when the hosing stops. That is the behaviour of a system with only one stable state.

## Step 5: Switch to mixed boundary conditions (`experiments/mixed_forcing/`, `experiments/mixed_control/`, `notebooks/05_mixed_boundary_conditions.ipynb`)

### The question

Strong salinity restoring ties the surface salinity to observations. That works against the salt-advection feedback, which is the feedback that can give the AMOC two stable states. It is like a large "delta" in Stommel's model, where salinity is held so tightly that only one equilibrium exists. Can I remove the restoring without changing the state the model is in? And does the model stay close to that state on its own?

### What I changed and why

"Mixed boundary conditions" means keeping temperature restoring (a stand-in for the atmosphere damping temperature anomalies, which is real) while replacing salinity restoring with a fixed freshwater flux (because nothing in reality pushes surface salinity back). The fixed flux is chosen to do on average what the restoring did at the end of the spin-up. That way the switch should change as little as possible at first, but salinity anomalies are no longer damped.

- Salinity restoring off: `tauSaltClimRelax = 0`. I checked in the source (`set_parms.F`) that this sets `doSaltClimRelax` to false, so the restoring rate is exactly zero. The run's configuration printout confirms `doThetaClimRelax = T` and `doSaltClimRelax = F`.
- `EmPmRFile = 'emp_mixed.bin'`: the original EmPmR plus the freshwater flux equivalent to the spin-up's restoring.
- `mixed_control`: 500 years (2001-2500) from the spin-up pickup, the same start as the restoring runs. It also serves as the control for step 6.

### The notebook, block by block

Part A builds the forcing.

1. The spin-up saved a 12-month climatology of `SRELAX` over its last 100 years (1901-2000). The restoring has a large seasonal cycle, from +14481 t/s (month 4) to -34170 t/s (month 9), against an annual mean of only -519.6 t/s. That is why I keep all 12 months instead of using the annual mean.
2. The model converts EmPmR to salt with the local surface salinity, so I use the model's own surface salinity averaged over the same 100 years (29.82 to 37.44 g/kg). The annual output does not resolve the seasonal cycle of salinity, so this is an approximation.
3. The conversion: a salt flux F_s (g/m^2/s) has the same local effect as an extra EmPmR of F_s / (1000 x S) m/s. The sign stays the same: restoring that adds salt becomes extra evaporation. This is the same formula as in step 4.
   - The equivalent flux ranges from -2.457e-07 to 3.349e-07 m/s.
   - Its annual-mean global total is +0.0160 Sv (net evaporation). This is not a mistake. The restoring removed salt on balance (-519.6 t/s), but most of its salt-adding happened in fresher water, where each gram of salt is worth more fresh water, so after dividing by the local salinity the sum comes out as net evaporation.
   - Over the hosing region it is -0.0993 Sv. There the restoring was freshening the surface, because the model's subpolar surface was saltier than observed.

   I kept the net amount rather than forcing it to zero, so the new forcing is exactly what the restoring did.
4. `emp_mixed.bin` = original + equivalent flux, same format. Read back from disk, it matches the intended field to within 9.38e-15 m/s.
5. The three step-6 hosing files are `emp_mixed.bin` plus the step-3 anomaly. Each adds 0.10000000, 0.30000000 and 0.50000000 Sv relative to `emp_mixed.bin`, in every month.

I ran Part A on its own before the model run, to make the files. When the whole notebook ran afterwards, it wrote byte-for-byte identical files, which I checked with `md5`.

Part B analyses the mixed control.

6. Check: it ended normally, and `SRELAX` is 0.00e+00 in every one of its 500 annual means.
7. AMOC index for the mixed and restoring controls, saved to `data/amoc_index_mixed_control.csv`.
8. Surface salinity in the hosing region and the whole Atlantic, and the global volume-mean salinity, in both controls.
9. A salt prediction (below).

### Results

| | year 2001 | mean 2491-2500 | lowest | trend, last 300 years |
| :--- | ---: | ---: | ---: | ---: |
| restoring control | 19.82 Sv | 19.94 Sv | 19.82 Sv | +0.025 Sv per century |
| mixed control | 19.78 Sv | 19.74 Sv | 16.51 Sv | -0.023 Sv per century |

The mixed control did not stay close to the spun-up state at first, and I report that as a result, not a failure.

- With no hosing at all, the AMOC fell below 19 Sv by 2009 and reached 16.51 Sv in 2038, a drop of about 3.3 Sv.
- It recovered above 19 Sv by 2072 and peaked at 19.64 Sv in 2087.
- It dipped again to 19.20 Sv in 2127.
- It then settled at 19.74 Sv, 0.21 Sv below the restoring control at the end.

So switching the boundary condition is a shock in itself, even with a flux built to match the restoring. Possible reasons, which I have not tested:

- The fixed flux matches the restoring's 100-year average, but not its year-to-year response to the actual surface salinity.
- Once the restoring is gone, nothing damps the salinity anomalies the switch creates.

The model did not collapse, and it found a state close to the original. But the first ~150 years of a mixed-BC run start with this adjustment, which step 6 needs to keep in mind.

Surface salinity drifted too. At 2491-2500 the mixed control's surface is 0.0154 g/kg saltier in the hosing region and 0.0482 g/kg saltier over the Atlantic than the restoring control. Its global volume-mean salinity is 0.01247 g/kg higher.

The salt prediction. In step 4 I found that the total salt in this model only responds to the net freshwater the ocean gains or loses. If that is right, I can predict the 0.01247 g/kg difference:

- The restoring control loses salt through `SRELAX`: -545.9 t/s on average, which is -0.00620 g/kg over 500 years.
- The mixed control has a net evaporation of +0.0160 Sv, which concentrates the salt by +0.00631 g/kg over 500 years. This includes the factor 1000/1035 explained in step 7.
- The predicted difference is +0.01251 g/kg. The measured difference is +0.01247 g/kg, within 0.3%. Without the 1000/1035 factor the prediction was +0.01273 g/kg, 2% too large.

This is good support for the salt explanation from step 4.

### What to look for in the figures

- `figures/05_restoring_equivalent_flux.png`: the freshwater flux that replaces the restoring, in m/yr (red = extra evaporation, blue = extra freshwater). There is extra freshwater along the subpolar North Atlantic and the Kuroshio region, and extra evaporation in the western tropical Pacific. These are the places where the model's surface salinity disagreed with observations.
- `figures/05_amoc_mixed_control.png`: the restoring control (black) is flat. The mixed control (green) drops by about 3.3 Sv over the first 40 years, recovers, dips again, and settles slightly lower.
- `figures/05_sss_mixed_minus_restore.png`: by the end, most of the ocean surface is slightly saltier in the mixed control, but the Labrador Sea and the area near Hudson Strait are fresher.

### What I learned

Removing salinity restoring is not a neutral change, even when the replacement flux is built to match it. The model went through a 3 Sv dip and a slow adjustment before settling about 0.2 Sv weaker. Without restoring, nothing pulls the salinity back, so whatever the system does next is free to persist. The salt prediction matched to within 0.3%, which confirms that this model's total salt only follows the net freshwater input.

## Step 6: Hosing with mixed boundary conditions (`experiments/mixed_hose*/`, `notebooks/06_hosing_with_mixed_bc.ipynb`)

### The question

Does the AMOC respond differently to the same hosing when the surface salinity is no longer restored? Does it recover when the hosing stops?

### What I ran

The same design as step 4: 0.1, 0.3 and 0.5 Sv for 200 years (2001-2200), then 300 years without hosing (2201-2500).

- Each run's `data` file is a copy of `mixed_control`'s: salinity restoring off, temperature restoring on.
- The "on" runs use the mixed-BC hosing files from step 5 (`emp_mixed.bin` plus the anomaly). The "off" runs use `emp_mixed.bin`.
- `mixed_control` from step 5 is the control, so `run_hosing_set.sh` got a `skip_control` option, and only 3 runs ran at once.
- The command was `bash scripts/run_hosing_set.sh mixed experiments/mixed_forcing ~/mitgcm_runs/hosing/spinup experiments/mixed_forcing/emp_mixed.bin skip_control`. The set took 1 h 14 min.

### A numerical failure

All six runs printed "Execution ended Normally". But the 0.5 Sv run's fields became NaN (not a number) in model year 2090, and stayed NaN in the "off" part too, which restarted from a pickup that was already NaN.

- The free-surface solver (`cg2d`) shows what happened. Its right-hand side, normally about 2.5, grew slowly over about 100 time steps, then exploded over the last 6: 1.6e2, 1.2e5, 7.8e10, 4.6e27, 1.0e74, 4.1e204, then NaN.
- The monitor statistics just before gave no warning. The CFL numbers were about 0.04, and the lowest salinity was 23.5 g/kg.
- It was not just the rising sea level. The hosing is not compensated, so sea level rises: `eta` was 1.9-4.6 m when it failed, but the restoring 0.5 Sv run reached 9.7 m without trouble.
- In the 0.3 Sv run, the same solver value never went above 2.6.

To find the cause, I reran the run in two stages (`experiments/diagnose_hose050/`, appendix of notebook 06).

- Stage 1 reran it from year 2000 to step 752150, about 150 steps before the failure. It took 409 s.
- Stage 2 restarted from stage 1's final checkpoint. It was saved as a rolling checkpoint, `pickup.ckptA`, so stage 2 uses `pickupSuff = 'ckptA'`. It ran the last steps with a full snapshot and a monitor block at every step.

Stage 2 reproduced the original solver values exactly, so the blow-up is deterministic. With a monitor check at every step, MITgcm stopped it at step 752303, when temperatures reached -590 to 4776 °C. In the original run that bounds check (`MON_SOLUTION`) only ran once a model year. By the next check the fields were NaN, and every comparison with NaN is false, so the check passed and the run "ended normally".

What the snapshots show:

- Where: one column at 78°N, 2°E, in the northernmost row of the grid, next to the model's closed wall at 80°N (the Fram Strait / Greenland Sea area). Its cells are 92 km wide, the narrowest in the model, because the 4° longitude lines converge there.
- The column before it failed (step 752150): a very fresh, cold surface layer (27.7 g/kg, -1.8 °C) on top of warm Atlantic water (7.2 °C at 300-450 m), with salinity rising from 27.7 to 34.7 g/kg over the top 450 m. That fresh cap is what the collapsed AMOC and the hosing leave behind. Nothing mixes it away under mixed boundary conditions.
- What grew: the vertical velocity through the whole column, from about 5e-5 to 9e-4 m/s over about 140 steps. With the 1-day tracer time step, the local vertical CFL number (|w| x 1 day / layer thickness) rose from 0.03 to 0.52.
- What then went wrong: new values appeared that were nowhere in the column before: 12.5 °C at 290 m, salinity 36.4 g/kg at depth, and -2.5 °C at 290 m (below freezing). The tutorial uses centred second-order advection for tracers (`tempAdvScheme = saltAdvScheme = 2`). This scheme can create new maxima and minima when the flow crosses a sharp gradient at a large CFL number. Once that started, temperature and sea level in the column ran away within a few steps.

My reading: a numerical instability, not a physical result. It needed three things together: an extremely sharp fresh cap (only in the collapsed, strongly hosed state), the long 1-day tracer time step, and a non-monotonic advection scheme, in the narrowest grid cell next to an artificial wall. I don't know why the vertical velocity started to grow. The sensitivity test below shows that the advection scheme matters. The model blew up after the AMOC had already collapsed. I report the 0.5 Sv run only up to 2089. The two lessons:

- "Ended normally" does not prove a run is healthy, so notebook 06 now also checks that every value is a finite number.
- A run that becomes unstable must not be quietly dropped or patched.

#### Integrity check of all runs (appendix B of notebook 06)

After this failure I checked every run, not just that one. For each run I scanned:

- the model log for NaN;
- every output file (not counting pickups) for non-finite values;
- the yearly temperature and salinity range from the monitor;
- temperatures below the local freezing point. I used the UNESCO formula (Millero 1978) on the annual-mean surface fields and on the 10-year 3D means.

| run | NaN lines in log | output files with non-finite values | T range (°C) | S range (g/kg) | lowest surface T - Tf (°C) |
| :--- | ---: | ---: | :--- | :--- | ---: |
| spinup | 0 | 0 / 6263 | -1.97 to 30.44 | 29.74 to 37.53 | 0.36 |
| restore_control | 0 | 0 / 1598 | -1.97 to 30.42 | 29.74 to 37.51 | 0.49 |
| restore_hose010/030/050, on and off | 0 | 0 in each | -1.97 to 30.42 | 29.74 to 37.52 | 0.49 to 0.51 |
| mixed_control | 0 | 0 / 1598 | -1.97 to 30.42 | 29.52 to 37.55 | 0.48 |
| mixed_hose010_on / off | 0 / 0 | 0 / 665, 0 / 976 | -1.98 to 30.42 | 27.87 to 37.77 | 0.48 / 0.58 |
| mixed_hose030_on / off | 0 / 0 | 0 / 665, 0 / 976 | -1.98 to 30.43 | 25.78 to 37.56 | 0.48 / 0.52 |
| mixed_hose050_on | 44577 | 356 / 665 | -1.99 to 30.42 (finite part) | 23.49 to 37.53 | 0.48 |
| mixed_hose050_off | 121242 | 949 / 976 | all NaN | all NaN | - |

- Every run except the 0.5 Sv mixed chain is free of NaN.
- No annual-mean surface cell was more than 0.1 °C below its freezing point in any run.
- In the 10-year 3D means, the only value below freezing was -0.17 °C, early in the spin-up, while the model adjusted from its initial state.
- These are annual and 10-year means, so a short instantaneous dip below freezing would not show here.

#### Early warning signs? (appendix C)

There were none at yearly resolution. Up to 2089, the 0.5 Sv run's largest vertical velocity (2.09e-05 m/s) and its CFL number (0.024) were the same as in the healthy 0.3 Sv run. The only differences were how fresh it was:

- lowest salinity anywhere: 23.49 against 26.56 g/kg;
- surface salinity at 78°N, 2°E: 27.63 against 29.22 g/kg.

The instability itself grew within about 150 model days, too fast for yearly output to catch.

#### Sensitivity test of the advection scheme (appendix D)

I restarted from a clean state and changed only the advection scheme. These runs are kept out of the comparison tables.

- Stage 0 reran the 0.5 Sv mixed-BC run from year 2000 to the end of model year 2088 (step 751680), about 620 steps before the instability, and saved that state.
- `orig` restarted from it with the original scheme. All 630 solver values were identical to the original run, step by step, so the blow-up reproduced exactly. This time the run was stopped at step 752310 with "SOLUTION IS HEADING OUT OF BOUNDS: tMin,tMax = NaN NaN", by the new NaN-safe check described below.
- `fluxlim` restarted from the same state with a flux-limited, monotonic scheme (`tempAdvScheme = saltAdvScheme = 33`) and nothing else changed. It ran to year 2200 with no NaN: 1345 monthly checks, temperature -1.98 to 30.42 °C, lowest salinity 23.00 g/kg. The AMOC stayed collapsed, between -0.23 and -0.11 Sv. The surface salinity at 78°N, 2°E went from 27.66 to 26.95 g/kg.

So the original scheme was needed for the blow-up, and changing it removed the blow-up without changing the collapsed state. This supports my reading that the failure was numerical. The fresh cap and the long tracer time step may still have been needed too. I only changed one thing.

#### Checking for instability in future runs

Following this failure, every future run will:

- run the monitor and its bounds check every model month (`monitorFreq = 2592000`) instead of every year;
- use a build with a NaN-safe bounds check. In `experiments/safety_code/mon_solution.F`, the test `(tMax-tMin) .GT. limit` is rewritten as `.NOT.((tMax-tMin) .LE. limit)`, which is also true when the range is NaN. This only changes when a run stops, not what it computes: the `orig` rerun with this build reproduced the original solver values exactly;
- be marked as failed by `run_experiment.sh` if "NaN" appears anywhere in `output.txt`, even if the model says it ended normally.

A monthly check alone would not have been enough. In this run, the temperatures were out of range for only about 4 steps before turning into NaN, so a check every 30 steps could easily miss that window. With the NaN-safe check, it does not matter whether the check lands before or after the NaNs appear.

#### A limitation this exposed

The model domain ends at 80°N with a closed wall. There is no Arctic Ocean and no outflow of fresh water from the Nordic Seas into it. The hosing water and the fixed freshening from `emp_mixed.bin` can only leave the Nordic Seas southward, through a circulation that has collapsed. So the fresh cap in the Nordic Seas, where the instability started, was probably much stronger than it would be in the real ocean, or in a model with an Arctic.

Along the way I also caught a false alarm of my own. A quick check showed a surface salinity of -1.37 g/kg in the 0.3 Sv run, but that came from a land point, where output values mean nothing. Over ocean cells the lowest value is 25.80 g/kg.

### The notebook, block by block

1. Checks: every run ended normally, every annual-mean velocity value is finite (all except the 0.5 Sv run after 2089), and the year-2200 pickups exist.
2. AMOC index for the mixed-BC runs, saved to `data/amoc_index_mixed.csv`. The restoring runs come from step 4's CSV file.
3. A side-by-side figure and a summary table. Each run is compared with its own control. I call a run "recovered" if its last 10 years are within 1 Sv of its control.
4. The strongest valid case. The 0.5 Sv mixed-BC run is not valid after 2089, so for the maps and streamfunctions I use 0.3 Sv, the strongest run that is valid throughout, for both boundary conditions, so they are compared at the same strength.
5. The step-4 overshoot check: surface salinity in the subtropical Atlantic during the restoring run's overshoot.

### Results

| set | hosing | end of hosing (2191-2200) | change | lowest (year) | end of recovery (2491-2500) | change | recovered? |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| restoring | control | 19.87 Sv | | | 19.94 Sv | | |
| restoring | 0.1 Sv | 17.60 Sv | -2.27 Sv | 17.21 (2018) | 19.99 Sv | +0.04 Sv | yes |
| restoring | 0.3 Sv | 14.27 Sv | -5.60 Sv | 13.52 (2020) | 20.06 Sv | +0.11 Sv | yes |
| restoring | 0.5 Sv | 11.70 Sv | -8.17 Sv | 10.77 (2021) | 20.15 Sv | +0.20 Sv | yes |
| mixed | control | 19.66 Sv | | | 19.74 Sv | | |
| mixed | 0.1 Sv | -0.04 Sv | -19.70 Sv | -0.04 (2200) | 0.00 Sv | -19.73 Sv | no |
| mixed | 0.3 Sv | -0.14 Sv | -19.80 Sv | -0.14 (2199) | -0.02 Sv | -19.76 Sv | no |
| mixed | 0.5 Sv | numerical failure after 2089 | | -0.15 (2089) | | | |

- Under mixed boundary conditions the AMOC collapsed in every run, even with 0.1 Sv.
  - It fell below 5 Sv in 2058 (0.1 Sv), 2021 (0.3 Sv) and 2017 (0.5 Sv).
  - By the end of hosing the AMOC index was zero within 0.15 Sv. The streamfunction shows no deep overturning cell at all, only a shallow tropical cell.
- It did not recover. 300 years after the hosing stopped, the AMOC index was still 0.00 and -0.02 Sv.
- Surface salinity in the hosing region (0.3 Sv runs, hosing run minus its control):
  - With restoring: -0.304 g/kg at the end of hosing, and -0.001 g/kg at the end of recovery.
  - With mixed BC: -5.729 g/kg at the end of hosing, and still -2.794 g/kg at the end of recovery. Over the whole Atlantic the mixed-BC surface was 3.013 g/kg fresher at the end of hosing and 1.534 g/kg fresher at the end.
- The step-4 overshoot idea (salty water building up in the subtropical Atlantic, then moving north) is not supported at the surface. In the restoring 0.5 Sv run, the subtropical Atlantic (10-40°N) surface was only 0.013 g/kg saltier at the end of hosing, and 0.008 g/kg during the overshoot (2221-2230). Restoring keeps the surface close to observations, so if salt built up, it would be below the surface. I did not check that, so the overshoot is still unexplained.

### A caveat

The hosing runs start at the same moment as the switch to mixed boundary conditions. Step 5 showed that this switch alone made the AMOC dip by 3.3 Sv over the first 40 years. So the hosing hit a circulation that was already adjusting. The mixed control recovered from that dip, while the 0.1 Sv run kept falling all the way to zero, so the collapse is caused by the hosing. But a cleaner design would start the hosing from the mixed control after it has settled, for example at year 2300.

### What to look for in the figures

- `figures/06_amoc_restoring_vs_mixed.png`: on the left, every restoring run weakens and then recovers. On the right, every mixed-BC run falls to zero and stays there after 2200. The cross marks where the 0.5 Sv run stopped being valid. The mixed control (black) shows the step-5 dip and recovery.
- `figures/06_sss_change_030.png`: each row has its own colour scale, because the changes differ ten-fold. With restoring, the freshening (up to about 0.5 g/kg) stays in the hosing region and is gone by the end. With mixed BC, the whole North Atlantic surface becomes several g/kg fresher, spreading south to the subtropics, and it is still clearly fresher 300 years after the hosing stopped.
- `figures/06_streamfunction_mixed_030.png`: the control's deep cell (about 20 Sv) is gone at the end of hosing, and still gone at the end of recovery. Only a shallow cell near the tropics remains.
- `figures/06_early_warning_check.png`: the 0.5 and 0.3 Sv runs' largest vertical velocity, lowest salinity and surface salinity at 78°N, 2°E, year by year. The velocity curves lie on top of each other until the blow-up.
- `figures/06_sensitivity_fluxlim.png`: the 0.5 Sv run continued with the flux-limited scheme from 2089 stays collapsed, like the 0.3 Sv run.
- `figures/06_blowup_column.png`: the column at 78°N, 2°E in the last ~150 steps of the 0.5 Sv mixed-BC run. Temperature and salinity profiles develop values that were not in the column before, while the local vertical CFL number grows from 0.03 to 0.52.
- `figures/06_restoring_overshoot_sss.png`: the restoring 0.5 Sv run's surface salinity change at the end of hosing and during the overshoot, with almost no change in the subtropics.

### What I learned

The boundary condition decides the answer. With restoring, the AMOC weakens in proportion to the hosing and always recovers: one stable state. With mixed boundary conditions, even 0.1 Sv for 200 years switches the AMOC off, and it stays off after the hosing stops. That is a second stable state, as in Stommel's model. The fresh surface layer that a weak AMOC leaves behind keeps the AMOC weak, and nothing pushes the salinity back. I also learned to check runs for NaNs, not just for "ended normally".

## Step 7: Comparison and checks (`notebooks/07_compare_and_check.ipynb`)

### The question

Taken together, how much does the AMOC weaken and does it recover, for each hosing strength under each boundary condition? Does the salt budget of every hosing run add up? And what in this model plays the role of Stommel's delta?

### What I did

No new model runs. The notebook reads the AMOC series saved by notebooks 04-06, the yearly global salinity statistics, and the annual surface fields.

### Results: the AMOC

| boundary condition | hosing | change at end of hosing (2191-2200) | change at end of recovery (2491-2500) | result |
| :--- | ---: | ---: | ---: | :--- |
| restoring | 0.1 Sv | -2.27 Sv (-11.4%) | +0.04 Sv | recovered |
| restoring | 0.3 Sv | -5.60 Sv (-28.2%) | +0.11 Sv | recovered |
| restoring | 0.5 Sv | -8.17 Sv (-41.1%) | +0.20 Sv | recovered |
| mixed | 0.1 Sv | -19.70 Sv (-100.2%) | -19.73 Sv | not recovered |
| mixed | 0.3 Sv | -19.80 Sv (-100.7%) | -19.76 Sv | not recovered |
| mixed | 0.5 Sv | -19.78 Sv (-100.6%) in 2089 | n/a | numerical failure after 2089 |

Each run is compared with its own control (19.87 Sv for restoring, 19.66 Sv for mixed, at 2191-2200). The table is saved to `data/summary_table.csv`.

### Results: the salt check

If the total salt only follows the net freshwater (steps 4 and 5), the change in global mean salinity of a hosing run, relative to its control, should be:

- dilution by the hosing: -S x H x (1000/1035) x t / V, where S is the run's global area-mean surface salinity, H the hosing, t the time it has been on and V the ocean volume (1.3227e18 m^3);
- plus, in the restoring runs, the salt added by the extra restoring: the time integral of the global SRELAX difference, divided by (rhoConst x V).

The yearly statistics are annual averages, so I compare them with the prediction at the middle of each year.

My first version left out the factor 1000/1035. It predicted 3.3-3.4% too much salinity loss in every mixed run, at every time. In the restoring runs the leftover was exactly the same in absolute terms as in the mixed run with the same hosing (for example +0.00055 g/kg at 2200 for 0.1 Sv in both). That told me the restoring part was right and the whole error was in the dilution term. A constant 3.4% matched 1 - 1000/1035 = 0.034. The model turns a freshwater mass flux into volume with `mass2rUnit = 1/rhoConst` (`ini_parms.F`), so 0.5 Sv of fresh water (5e8 kg/s) enters the model as 0.483 Sv of seawater volume. In mass terms that is also the physically right dilution.

With the factor included, measured and predicted agree to within 0.1% in every run, at 2100, 2200 and 2500. Examples at 2200: restoring 0.5 Sv -0.01647 measured against -0.01646 predicted; mixed 0.3 Sv -0.04684 against -0.04684. The same factor also brought the step 5 prediction from 2% to 0.3%.

So the salt budget closes, and it tells a clear story:

- In the mixed-BC runs, the global salinity falls in a straight line while the hosing is on, by the full dilution, and stays flat afterwards.
- In the restoring runs, the restoring puts back most of the salt the hosing removes. At 2200, the 0.5 Sv run had lost only -0.01647 g/kg, against -0.0817 g/kg for dilution alone.

### What plays the role of Stommel's delta?

In Stommel's model, delta compares how fast salinity is pulled back to its forcing value with how fast temperature is. In this model:

- With restoring, temperature relaxes over 60 days and salinity over 180 days, so delta = 60/180 = 0.333. Over the 50 m top layer, the salinity restoring acts like a "piston velocity" of 100 m per year.
- With mixed boundary conditions, `tauSaltClimRelax = 0`, so delta = 0. Salinity is set by a fixed flux, and nothing pulls it back.

The ratio alone is not the whole story. In project 14, Stommel's model already had two equilibria with delta = 1/6, smaller than the 1/3 here. What matters is how strong the salinity restoring is compared with the circulation that carries salt into the sinking region. Step 4 measured that directly: the restoring cancelled 80-87% of the hosing in the region within a few years, so the salt-advection feedback could never take hold. With delta = 0 nothing opposes that feedback, and the model behaves like Stommel's model with a fixed salt flux: a weakened AMOC carries less salt north, the surface gets fresher, and the AMOC weakens further, down to the other stable state.

### What to look for in the figures

- `figures/07_summary.png`: on the left, the restoring runs weaken in a near-straight line with hosing strength, while every mixed-BC run sits at about -20 Sv, completely off. On the right, after 300 years without hosing, the restoring runs are back at zero change and the mixed-BC runs are still off.
- `figures/07_salt_check.png`: measured (solid) and predicted (dashed) salinity change. The dashed lines are hidden under the solid ones because they agree so closely. The restoring runs (left) bend as the restoring puts salt back, and the mixed runs (right) fall in straight lines and then stay flat.

### What I learned

The surface boundary condition decides whether this model has one stable AMOC state or two. Restoring salinity, which has no real counterpart in nature, removes the salt-advection feedback almost completely. The salt check closed to 0.1% only after I found a factor of 1000/1035 in how the model converts freshwater mass to volume. A constant percentage error across all runs pointed straight to it.

## What I would do next

- Run a similar experiment in CESM, a coupled climate model with an atmosphere, sea ice and an Arctic Ocean. There, the "boundary condition" is not a choice: the atmosphere responds to the ocean, sea ice forms and melts, and freshwater can leave through the Arctic. That would show which of the two behaviours here is closer to a full climate system. It would need many more computing resources than a laptop.
- Add an Antarctic ice shelf. Meltwater from Antarctica changes the density of the deep water that forms in the Southern Ocean. That deep water sits below the North Atlantic Deep Water and competes with it, so Antarctic melt could change how the AMOC responds to North Atlantic hosing. In MITgcm this could start with the `shelfice` package, in a configuration that includes the ice-shelf cavities.
- Within this setup:
  - Start the mixed-BC hosing from the mixed control after it has settled (for example year 2300), so the hosing and the switch of boundary condition are not mixed together.
  - Use weaker hosing (0.01-0.1 Sv) to find the threshold, because even 0.1 Sv collapsed the AMOC here.
  - Then look for hysteresis by stepping the hosing up and back down, like the optional step 8.
