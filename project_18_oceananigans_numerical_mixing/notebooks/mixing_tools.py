"""Shared analysis functions for project 18 (numerical mixing in the lock exchange).

Plain functions, used by the notebooks in this folder. Units: metres, seconds, m/s^2 for buoyancy,
kg/m^3 for density, J per metre of channel width for energies (the runs are 2D).
"""

import os

import numpy as np
import pandas as pd
import xarray as xr

GRAVITY = 9.81            # m/s^2
RHO_0 = 1000.0            # reference density (kg/m^3)
DELTA_B = GRAVITY * 5.0 / RHO_0   # buoyancy difference between the two water masses (m/s^2)
LENGTH = 64e3             # channel length (m)
DEPTH = 20.0              # channel depth (m)
RUNS = os.path.expanduser("~/oceananigans_runs/numerical_mixing")


def load_run(name):
    """Open one run: returns (fields dataset, log table, run_info dictionary)."""
    run_dir = os.path.join(RUNS, name)
    ds = xr.open_dataset(os.path.join(run_dir, "fields.nc"))
    log = pd.read_csv(os.path.join(run_dir, "log.csv"))
    info = {}
    for line in open(os.path.join(run_dir, "run_info.txt")):
        key, value = line.split(" = ", 1)
        info[key.strip()] = value.strip()
    return ds, log, info


def density(b):
    """Density (kg/m^3) from buoyancy, b = -g (rho - rho0) / rho0, so rho = rho0 (1 - b / g)."""
    return RHO_0 * (1.0 - b / GRAVITY)


def reference_potential_energy(b, dx, dz):
    """Reference potential energy (J per metre of width) of one buoyancy snapshot b[z, x].

    Sorting (Winters et al. 1995): take every cell's density and rearrange the cells into the most
    stable possible state, densest at the bottom, lightest at the top, without mixing anything.
    All cells here have the same area dx * dz, so I sort the densities from heaviest to lightest and
    fill the grid row by row from the bottom. RPE = g * sum(rho_sorted * height * cell area), with the
    height measured upward from the bottom of the channel.
    """
    nz, nx = b.shape
    rho_sorted = np.sort(density(np.asarray(b)).ravel())[::-1]        # heaviest first
    rows = rho_sorted.reshape(nz, nx)                                 # row 0 = bottom row
    heights = (np.arange(nz) + 0.5) * dz                              # cell-centre heights above the bottom (m)
    return GRAVITY * np.sum(rows * heights[:, None]) * dx * dz


def rpe_series(ds):
    """RPE (J per metre of width) at every snapshot."""
    dx = float(ds.x_caa[1] - ds.x_caa[0])
    dz = float(ds.z_aac[1] - ds.z_aac[0])
    b = ds["b"].values
    return np.array([reference_potential_energy(b[n], dx, dz) for n in range(b.shape[0])])


def buoyancy_variance(ds):
    """Volume-mean variance of b about its mean, (m/s^2)^2, at every snapshot."""
    b = ds["b"].values
    mean = b.mean(axis=(1, 2), keepdims=True)
    return ((b - mean) ** 2).mean(axis=(1, 2))


def front_positions(ds):
    """Positions (m) of the two gravity-current noses at every snapshot.

    The dense current runs to the right along the bottom: its nose is the right-most bottom cell
    with b below half the difference. The light current runs to the left along the top: its nose is
    the left-most top cell with b above half the difference. Returns (right_front, left_front).
    """
    x = ds.x_caa.values
    bottom = ds["b"].isel(z_aac=0).values      # bottom row, [time, x]
    top = ds["b"].isel(z_aac=-1).values        # top row
    right, left = [], []
    for n in range(bottom.shape[0]):
        dense = np.where(bottom[n] < DELTA_B / 2)[0]
        light = np.where(top[n] > DELTA_B / 2)[0]
        right.append(x[dense.max()] if len(dense) else np.nan)
        left.append(x[light.min()] if len(light) else np.nan)
    return np.array(right), np.array(left)


def front_speed(times, positions, t_start, t_end):
    """Speed (m/s) from a straight-line fit of position against time between t_start and t_end (s)."""
    use = (times >= t_start) & (times <= t_end)
    slope, intercept = np.polyfit(times[use], positions[use], 1)
    return slope


def rpe_rate_at(times, rpe, t_target, window):
    """dRPE/dt (W per metre of width) from a straight-line fit over [t_target - window, t_target]."""
    use = (times >= t_target - window) & (times <= t_target + 1e-6)
    slope, intercept = np.polyfit(times[use], rpe[use], 1)
    return slope
