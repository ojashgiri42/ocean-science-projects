"""Shared functions for reading the hosing runs and computing the AMOC.

The streamfunction and AMOC index are the same calculation that notebook 02
explains step by step. They are kept here so notebooks 04, 06 and 07 all use
exactly the same code.
"""

import os
import glob

import numpy as np
from MITgcmutils import rdmds

STEPS_PER_YEAR = 360          # deltaTClock = 1 day, 360-day model year
SV = 1.0e6                    # 1 Sverdrup = 1e6 m^3/s
RHO_FRESH = 1000.0            # rhoConstFresh in the data file (kg/m^3)
AMOC_LAT_MIN = 20.0           # AMOC index: max of psi between 20 and 60 N ...
AMOC_LAT_MAX = 60.0
AMOC_DEPTH_MIN = 500.0        # ... and below 500 m


def read_grid(run_dir, data_dir):
    """Grid fields and masks needed for the AMOC and for area integrals, in one dictionary."""
    grid = {}
    grid["xc"] = rdmds(os.path.join(run_dir, "XC"))
    grid["yc"] = rdmds(os.path.join(run_dir, "YC"))
    grid["lat_v"] = rdmds(os.path.join(run_dir, "YG"))[:, 0]
    grid["dxg"] = rdmds(os.path.join(run_dir, "DXG"))
    grid["drf"] = rdmds(os.path.join(run_dir, "DRF")).ravel()
    grid["hfacc"] = rdmds(os.path.join(run_dir, "hFacC"))
    grid["hfacs"] = rdmds(os.path.join(run_dir, "hFacS"))
    grid["rac"] = rdmds(os.path.join(run_dir, "RAC"))
    grid["ocean"] = grid["hfacc"][0] > 0
    grid["z_interface"] = np.concatenate([[0.0], np.cumsum(grid["drf"])])
    grid["atlantic"] = np.loadtxt(os.path.join(data_dir, "atlantic_mask.csv"), delimiter=",").astype(bool)
    grid["region"] = np.loadtxt(os.path.join(data_dir, "hosing_region_mask.csv"), delimiter=",").astype(bool)
    # A v-point is Atlantic only if the cells on both sides of it are Atlantic (as in notebook 02)
    atlantic_v = np.zeros_like(grid["atlantic"])
    for j in range(1, atlantic_v.shape[0]):
        atlantic_v[j, :] = grid["atlantic"][j, :] & grid["atlantic"][j - 1, :]
    grid["atlantic_v"] = atlantic_v
    return grid


def atlantic_streamfunction(v, grid):
    """Atlantic overturning streamfunction (Sv) at level boundaries: northward transport above each depth."""
    n_levels, n_lat, n_lon = v.shape
    transport = np.zeros((n_levels, n_lat))
    for k in range(n_levels):
        for j in range(n_lat):
            for i in range(n_lon):
                if grid["atlantic_v"][j, i]:
                    area = grid["dxg"][j, i] * grid["drf"][k] * grid["hfacs"][k, j, i]
                    transport[k, j] += v[k, j, i] * area
    psi = np.zeros((n_levels + 1, n_lat))
    for k in range(n_levels):
        psi[k + 1, :] = psi[k, :] + transport[k, :] / SV
    return psi


def amoc_index(psi, grid):
    """Maximum of psi between AMOC_LAT_MIN and AMOC_LAT_MAX and below AMOC_DEPTH_MIN (Sv)."""
    lat_ok = (grid["lat_v"] >= AMOC_LAT_MIN) & (grid["lat_v"] <= AMOC_LAT_MAX)
    depth_ok = grid["z_interface"] >= AMOC_DEPTH_MIN
    return psi[depth_ok, :][:, lat_ok].max()


def output_iterations(run_dir, prefix):
    """Step numbers of all output files with this prefix, sorted."""
    files = sorted(glob.glob(os.path.join(run_dir, prefix + ".*.data")))
    return [int(os.path.basename(f).split(".")[1]) for f in files]


def amoc_series(run_dirs, grid):
    """AMOC index for every annual mean in one or more runs (in order). Returns (years, amoc)."""
    years = []
    amoc = []
    for run_dir in run_dirs:
        for it in output_iterations(run_dir, "vvel_annual"):
            v = rdmds(os.path.join(run_dir, "vvel_annual"), it)
            years.append(it / STEPS_PER_YEAR)
            amoc.append(amoc_index(atlantic_streamfunction(v, grid), grid))
    return np.array(years), np.array(amoc)


def annual_field(run_dir, prefix, field, iteration):
    """One field from a multi-field annual output file, e.g. SRELAX from surf_flux_annual."""
    data, its, meta = rdmds(os.path.join(run_dir, prefix), iteration, returnmeta=True)
    return data[meta["fldlist"].index(field)]
