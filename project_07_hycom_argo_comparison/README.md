# 07. HYCOM vs Argo Temperature Profile (Model–Observation Comparison)

**Question.** How does a HYCOM temperature profile at the nearest model grid point compare with an Argo float temperature profile?

**Data.** Argo float 6901814, first profile (2013-01-10, 61.50°S, 178.66°E), and HYCOM `GLBy0.08/expt_93.0` via OPeNDAP at the nearest grid point, first model time step (`isel(time=0)`). Upper 1000 m.

**Method.** Convert the Argo longitude to 0–360°, select the nearest HYCOM grid point, and plot both temperature profiles. The model time is **not** matched to the Argo date, because the HYCOM dataset starts years after 2013. Only temperature is compared, Argo pressure (dbar) is plotted against model depth (m) without conversion, and no bias statistic is computed.

**Finding.** Both profiles show a subsurface temperature minimum near 130–200 m and a warmer layer below. The model surface layer is about 1 °C colder than the Argo profile (about 2.2 vs 3.3 °C), and the model is warmer by up to ~0.5 °C between ~250 and 500 m. The profiles converge below ~750 m. Because the two profiles are years apart, these differences cannot be attributed to model error.

![07. HYCOM vs Argo Temperature Profile (Model–Observation Comparison)](figures/model_validation.png)