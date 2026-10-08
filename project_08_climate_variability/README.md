# 08. Niño 3.4 SST Anomalies from NOAA OISST

**Question.** How has monthly SST in the Niño 3.4 region varied since the early 1980s?

**Data.** NOAA OISST v2 monthly mean SST (1°), an observational analysis product (Optimum Interpolation SST), accessed via OPeNDAP from NOAA PSL. Niño 3.4 region: 5°S–5°N, 170–120°W.

**Method.** Take a cos(latitude)-weighted area mean, then subtract the monthly climatology computed over the full record. Months with anomaly ≥ +0.5 °C are shaded red and months ≤ −0.5 °C blue. This is a simple monthly threshold, not the NOAA ONI definition (no 3-month running mean, no persistence requirement, no fixed base period).

**Finding.** The largest positive anomalies (about +2.5 to +3 °C) occur in 1982–83, 1997–98 and 2015–16. The strongest negative anomalies (about −2 to −2.4 °C) occur in 1988–89, 1999–2000 and 2007–08/2010–11.

![08. Niño 3.4 SST Anomalies from NOAA OISST](figures/climate_variability_timeseries.png)
