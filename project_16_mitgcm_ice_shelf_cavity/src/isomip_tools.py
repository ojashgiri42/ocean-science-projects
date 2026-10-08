"""Shared functions for reading the isomip runs (project 16).

Plain functions only, used by notebooks 02-05.
"""

import os
import glob

import numpy as np
from MITgcmutils import rdmds

NX, NY, NR = 50, 100, 30
RHO_ICE = 917.0                  # ice density used for melt in metres of ice (rhoShelfIce default, kg/m^3)
SECONDS_PER_YEAR = 365 * 86400.0 # this setup has no calendar; I use a 365-day year for melt rates
DELTA_T = 1800.0                 # model time step (s)
STEPS_PER_YEAR = 17520           # 365 days of 1800 s steps


def read_grid(run_dir, input_dir, ice_file="icetopo.exp1"):
    """Grid, masks and geometry in one dictionary.

    ice_file is the ice-draft file inside input_dir (step 3 uses a version
    with the ice removed in the northern restoring zone).
    """
    g = {}
    g["xc"] = rdmds(os.path.join(run_dir, "XC"))
    g["yc"] = rdmds(os.path.join(run_dir, "YC"))
    g["rac"] = rdmds(os.path.join(run_dir, "RAC"))
    g["dxg"] = rdmds(os.path.join(run_dir, "DXG"))
    g["dyg"] = rdmds(os.path.join(run_dir, "DYG"))
    g["drf"] = rdmds(os.path.join(run_dir, "DRF")).ravel()
    g["rc"] = rdmds(os.path.join(run_dir, "RC")).ravel()          # depth of cell centres (negative, m)
    g["rf"] = rdmds(os.path.join(run_dir, "RF")).ravel()          # depth of cell faces (negative, m)
    g["hfacc"] = rdmds(os.path.join(run_dir, "hFacC"))
    g["hfacw"] = rdmds(os.path.join(run_dir, "hFacW"))
    g["hfacs"] = rdmds(os.path.join(run_dir, "hFacS"))
    g["depth"] = rdmds(os.path.join(run_dir, "Depth"))           # seafloor depth (positive, m)
    ny, nx = g["depth"].shape                                     # 100 x 50 here, 200 x 100 on the step 6 fine grid
    g["ice_draft"] = np.fromfile(os.path.join(input_dir, ice_file), dtype=">f8").reshape(ny, nx)   # negative, m
    g["iced"] = (g["ice_draft"] < 0) & (g["depth"] > 0)          # ocean columns covered by ice
    g["wet"] = g["hfacc"].sum(axis=0) > 0                         # columns that contain water
    g["volume"] = np.sum(g["hfacc"] * g["drf"][:, None, None] * g["rac"][None, :, :])
    # index of the top wet cell under the ice in each column (-1 for land)
    ktop = np.full((ny, nx), -1)
    for j in range(ny):
        for i in range(nx):
            wet_levels = np.where(g["hfacc"][:, j, i] > 0)[0]
            if len(wet_levels) > 0:
                ktop[j, i] = wet_levels[0]
    g["ktop"] = ktop
    return g


def output_iterations(run_dir, prefix):
    files = sorted(glob.glob(os.path.join(run_dir, prefix + ".*.data")))
    return [int(os.path.basename(f).split(".")[1]) for f in files]


def read_field(run_dir, prefix, field, iteration):
    """One field from a multi-field output file."""
    data, its, meta = rdmds(os.path.join(run_dir, prefix), iteration, returnmeta=True)
    return data[meta["fldlist"].index(field)]


def melt_rate(fw_flux):
    """Melt rate in metres of ice per (365-day) year from SHIfwFlx.

    SHIfwFlx is a fresh-water mass flux in kg/m^2/s, positive UPWARD (out of the
    ocean), so melting is negative. Dividing by the ice density gives metres of ice
    per second; multiplying by the seconds in a year gives m/yr. The minus sign
    makes melting positive and refreezing negative.
    """
    return -fw_flux / RHO_ICE * SECONDS_PER_YEAR


def area_mean(field, grid, mask=None):
    """Area-weighted mean over the wet columns, or over `mask` if it is given."""
    if mask is None:
        mask = grid["wet"]
    w = grid["rac"] * mask
    return np.sum(field * w) / np.sum(w)


def melt_series(run_dirs, grid, mask=None):
    """Area-mean melt rate (m/yr) for every monthly mean, joined over run parts. Returns (days, melt).

    mask: columns to average over (default: all wet columns).
    """
    days = []
    melt = []
    for run_dir in run_dirs:
        for it in output_iterations(run_dir, "shelfice_monthly"):
            fw = read_field(run_dir, "shelfice_monthly", "SHIfwFlx", it)
            days.append(it * DELTA_T / 86400.0)
            melt.append(area_mean(melt_rate(fw), grid, mask))
    return np.array(days), np.array(melt)


def freezing_point(salt, pressure_dbar):
    """Freezing point used by the three-equation model in pkg/shelfice (linear in S and p)."""
    return -0.0575 * salt + 0.0901 - 7.61e-4 * pressure_dbar


def pressure_dbar(depth_m, rho=1030.0, g=9.81):
    """Hydrostatic pressure (dbar) at a depth (m, positive down), as the model's rhoConst gives."""
    return rho * g * np.abs(depth_m) / 1e4
