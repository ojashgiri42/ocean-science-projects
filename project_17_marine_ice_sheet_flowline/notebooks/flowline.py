"""A 1D marine ice-sheet flowline model (project 17), written as small plain functions.

Units used everywhere in this file:
    length    metres (m), with x measured from the ice divide (x = 0)
    time      years (yr); one year = SECONDS_PER_YEAR seconds
    velocity  m/yr
    stress    Pa
    A         Pa^-3 yr^-1  (Glen's law rate factor, converted from Pa^-3 s^-1)
    C         Pa m^(-1/3) yr^(1/3)  (sliding coefficient, converted from Pa m^(-1/3) s^(1/3))

Parameter values follow MISMIP experiment 3 (Pattyn et al. 2012, Tables 4 and 5).
"""

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

# ----------------------------------------------------------------------------------
# Constants (MISMIP, Pattyn et al. 2012, Table 4)
# ----------------------------------------------------------------------------------
SECONDS_PER_YEAR = 31556926.0   # 365.2422 days; my choice, the MISMIP text I could read does not state it
RHO_I = 900.0                   # ice density (kg/m^3)
RHO_W = 1000.0                  # seawater density (kg/m^3)
G = 9.8                         # gravitational acceleration (m/s^2)
N_GLEN = 3.0                    # exponent n in Glen's flow law
M_SLIDE = 1.0 / 3.0             # exponent m in the sliding law, basal stress = C |u|^(m-1) u
C_SI = 7.624e6                  # sliding coefficient (Pa m^(-1/3) s^(1/3))
ACCUMULATION = 0.3              # snowfall minus melt at the surface, a (m of ice per year)


def a_per_year(a_si):
    """Convert Glen's rate factor A from Pa^-3 s^-1 to Pa^-3 yr^-1."""
    return a_si * SECONDS_PER_YEAR


def c_per_year(c_si):
    """Convert the sliding coefficient C from Pa m^(-1/3) s^(1/3) to Pa m^(-1/3) yr^(1/3).

    Basal stress = C_si * u_si^m with u in m/s. A velocity in m/yr is u_yr = u_si * SECONDS_PER_YEAR,
    so C_si * u_si^m = C_si * (u_yr / SECONDS_PER_YEAR)^m = (C_si / SECONDS_PER_YEAR^m) * u_yr^m.
    """
    return c_si / SECONDS_PER_YEAR ** M_SLIDE


C_YR = c_per_year(C_SI)


# ----------------------------------------------------------------------------------
# Bed and flotation (step 1)
# ----------------------------------------------------------------------------------
def bed(x):
    """Bed elevation b(x) in m (positive up, sea level = 0), MISMIP experiment 3 (Pattyn et al. 2012, Eq. 16).

    x is in m. The polynomial uses x / 750 km.
    """
    s = x / 750e3
    return 729.0 - 2184.8 * s ** 2 + 1031.72 * s ** 4 - 151.72 * s ** 6


def bed_slope(x):
    """db/dx (dimensionless), the derivative of bed(x)."""
    s = x / 750e3
    return (-2.0 * 2184.8 * s + 4.0 * 1031.72 * s ** 3 - 6.0 * 151.72 * s ** 5) / 750e3


def flotation_thickness(x):
    """Ice thickness (m) at which ice would just float at x: H_f = -(rho_w / rho_i) * b.

    A column of ice floats when its weight, rho_i * g * H, equals the weight of the seawater it pushes
    aside, rho_w * g * (depth of the bed below sea level). Where the bed is above sea level, ice can
    never float there, so I return 0.
    """
    b = bed(x)
    return np.where(b < 0.0, -(RHO_W / RHO_I) * b, 0.0)


# ----------------------------------------------------------------------------------
# Schoof (2007) grounding-line flux (step 1)
# ----------------------------------------------------------------------------------
def schoof_prefactor(a_yr, c_yr):
    """The factor K in Schoof's flux formula q = K * h^((m+n+3)/(m+1)), in year units.

    K = [A (rho_i g)^(n+1) (1 - rho_i/rho_w)^n / (4^n C)]^(1/(m+1))   (Pattyn et al. 2012, Eq. A12)
    """
    inside = a_yr * (RHO_I * G) ** (N_GLEN + 1.0) * (1.0 - RHO_I / RHO_W) ** N_GLEN / (4.0 ** N_GLEN * c_yr)
    return inside ** (1.0 / (M_SLIDE + 1.0))


def schoof_exponent():
    """The exponent (m + n + 3) / (m + 1) on the grounding-line thickness (4.75 for n = 3, m = 1/3)."""
    return (M_SLIDE + N_GLEN + 3.0) / (M_SLIDE + 1.0)


def schoof_flux(x_g, a_yr, c_yr=C_YR):
    """Ice flux (m^2/yr, volume per year per metre of width) across a grounding line at x_g (m)."""
    h = flotation_thickness(x_g)
    return schoof_prefactor(a_yr, c_yr) * h ** schoof_exponent()


def schoof_flux_slope(x_g, a_yr, c_yr=C_YR):
    """dq/dx_g (m/yr): how fast the grounding-line flux changes when the grounding line moves.

    q = K h^p with h = -(rho_w/rho_i) b(x_g), so dq/dx_g = p * q / h * dh/dx_g
    and dh/dx_g = -(rho_w/rho_i) * db/dx.
    """
    h = flotation_thickness(x_g)
    q = schoof_flux(x_g, a_yr, c_yr)
    dh_dx = -(RHO_W / RHO_I) * bed_slope(x_g)
    return schoof_exponent() * q / h * dh_dx


def steady_a_for_position(x_g, c_yr=C_YR, accumulation=ACCUMULATION):
    """The rate factor A (Pa^-3 yr^-1) for which a grounding line at x_g (m) is a steady state.

    Steady state: everything that falls upstream leaves across the grounding line, a * x_g = q(x_g).
    q is proportional to A^(1/(m+1)), so I can solve for A directly:
    A = (a x_g / h^p)^(m+1) * 4^n C / ((rho_i g)^(n+1) (1 - rho_i/rho_w)^n).
    """
    h = flotation_thickness(x_g)
    ratio = accumulation * x_g / h ** schoof_exponent()
    return ratio ** (M_SLIDE + 1.0) * 4.0 ** N_GLEN * c_yr / ((RHO_I * G) ** (N_GLEN + 1.0) * (1.0 - RHO_I / RHO_W) ** N_GLEN)



def schoof_steady_positions(a_yr, x_min, x_max, dx_search=100.0, c_yr=C_YR, accumulation=ACCUMULATION):
    """All grounding-line positions (m) between x_min and x_max where snowfall upstream equals
    Schoof's flux, a x_g = q(x_g). I evaluate q - a x_g every dx_search metres, look for sign
    changes and refine each with brentq. Returns a list of (x_g, stable) pairs, where stable means
    dq/dx_g > a."""
    xs = np.arange(x_min, x_max, dx_search)
    mismatch = schoof_flux(xs, a_yr, c_yr) - accumulation * xs
    states = []
    for i in np.where(np.sign(mismatch[1:]) != np.sign(mismatch[:-1]))[0]:
        xg = brentq(lambda x: schoof_flux(x, a_yr, c_yr) - accumulation * x, xs[i], xs[i + 1])
        states.append((xg, schoof_flux_slope(xg, a_yr, c_yr) > accumulation))
    return states

# ----------------------------------------------------------------------------------
# Grid (step 2)
# ----------------------------------------------------------------------------------
def make_grid(length, dx):
    """A staggered grid from the divide (x = 0) to the calving front (x = length), both in m.

    Thickness lives at the N cell centres, velocity at the N + 1 cell edges:

        edge:    0     1     2    ...   N-1    N
                 |--*--|--*--|--  ...  --|--*--|
        centre:     0     1      ...       N-1

    Edge 0 is the ice divide, edge N the calving front. Returns (x_edges, x_centres).
    """
    n_cells = int(round(length / dx))
    x_edges = np.arange(n_cells + 1) * dx
    x_centres = (np.arange(n_cells) + 0.5) * dx
    return x_edges, x_centres


def surface_elevation(h, b):
    """Ice surface elevation (m). Grounded ice sits on the bed (s = b + H); floating ice floats
    with a fraction rho_i/rho_w of its thickness below sea level (s = (1 - rho_i/rho_w) H).
    Ice is grounded where b + H is the higher of the two, so the surface is the maximum."""
    return np.maximum(b + h, (1.0 - RHO_I / RHO_W) * h)


# ----------------------------------------------------------------------------------
# Shallow-shelf (SSA) velocity solver (step 2)
# ----------------------------------------------------------------------------------
def effective_viscosity(u, dx, a_yr, strain_rate_reg):
    """Effective viscosity eta (Pa yr) at the cell centres, from Glen's flow law.

    eta = 0.5 * A^(-1/n) * (strain rate)^((1 - n)/n), with the strain rate du/dx taken between
    the two edges of each cell. Where the ice barely stretches (du/dx near 0, e.g. at the divide)
    this would be infinite, so I use sqrt((du/dx)^2 + reg^2) instead of |du/dx|: a tiny floor on
    the strain rate (strain_rate_reg, 1/yr) that keeps eta finite and changes nothing where the
    ice really stretches.
    """
    dudx = (u[1:] - u[:-1]) / dx
    strain_rate = np.sqrt(dudx ** 2 + strain_rate_reg ** 2)
    return 0.5 * a_yr ** (-1.0 / N_GLEN) * strain_rate ** ((1.0 - N_GLEN) / N_GLEN)


def basal_drag_coefficient(u, c_yr, grounded_fraction, speed_reg):
    """beta (Pa yr/m) at the cell edges, so that basal stress = beta * u.

    The sliding law is basal stress = C |u|^(m-1) u. I write it as beta * u with
    beta = C |u|^(m-1), taken from the previous velocity (Picard). Since m - 1 < 0, beta would be
    infinite where u = 0 (the divide), so I use sqrt(u^2 + reg^2) instead of |u| (speed_reg in
    m/yr). grounded_fraction (0 to 1 at each edge) switches friction off under floating ice.
    """
    speed = np.sqrt(u ** 2 + speed_reg ** 2)
    return grounded_fraction * c_yr * speed ** (M_SLIDE - 1.0)


def calving_front_stress(h_front, b_front):
    """Depth-integrated stress the ice front must carry (Pa m): ice pressure minus water pressure.

    The ice column pushes outward with 0.5 rho_i g H^2; seawater pushes back on the submerged
    part (depth D below sea level) with 0.5 rho_w g D^2. For floating ice D = (rho_i/rho_w) H, which
    gives 0.5 rho_i g (1 - rho_i/rho_w) H^2.
    """
    s = surface_elevation(h_front, b_front)
    depth_below_sea = max(0.0, -(s - h_front))
    return 0.5 * RHO_I * G * h_front ** 2 - 0.5 * RHO_W * G * depth_below_sea ** 2


def solve_ssa(h, b_centres, dx, a_yr, grounded_fraction, u_guess=None, c_yr=C_YR,
              strain_rate_reg=1e-6, speed_reg=1e-3, tol=1e-10, max_iter=200):
    """Solve the 1D shallow-shelf momentum balance for the velocity at the cell edges (m/yr).

    At every interior edge j (1 to N-1):
        d/dx [ 4 eta H du/dx ]  -  beta u  =  rho_i g H ds/dx
        (stretching stress)       (friction)   (driving stress)
    with u = 0 at the divide (edge 0) and, at the calving front (edge N),
        4 eta H du/dx = calving_front_stress.

    eta and beta depend on u, so the equations are nonlinear. Picard iteration: compute eta and
    beta from the current u, solve the resulting linear (tridiagonal) system, repeat until the
    largest change in u is below tol times the largest speed.

    h, b_centres: thickness and bed at the N cell centres (m). grounded_fraction: N + 1 values.
    Returns (u, number of iterations, list of relative changes).
    """
    n_cells = len(h)
    s = surface_elevation(h, b_centres)
    if u_guess is None:
        u = np.zeros(n_cells + 1)
    else:
        u = u_guess.copy()

    # driving stress at the interior edges, from the surface slope between neighbouring centres
    h_edge = 0.5 * (h[:-1] + h[1:])
    driving = RHO_I * G * h_edge * (s[1:] - s[:-1]) / dx          # edges 1 .. N-1
    front = calving_front_stress(h[-1], b_centres[-1])

    changes = []
    for iteration in range(1, max_iter + 1):
        eta = effective_viscosity(u, dx, a_yr, strain_rate_reg)    # centres 0 .. N-1
        beta = basal_drag_coefficient(u, c_yr, grounded_fraction, speed_reg)   # edges 0 .. N
        k = 4.0 * eta * h / dx ** 2                                # stiffness of each cell

        # Unknowns: u_1 .. u_N (u_0 = 0), so edge j is row j - 1. solve_banded wants the three
        # diagonals stored as rows: ab[0] = upper diagonal, ab[1] = main diagonal, ab[2] = lower diagonal,
        # with ab[0, c] = matrix[c - 1, c] and ab[2, c] = matrix[c + 1, c].
        ab = np.zeros((3, n_cells))
        rhs = np.zeros(n_cells)
        # interior edges j = 1 .. N-1 (rows 0 .. N-2):
        #   k[j] (u_{j+1} - u_j) - k[j-1] (u_j - u_{j-1}) - beta_j u_j = driving_j
        ab[1, :n_cells - 1] = -(k[1:] + k[:-1]) - beta[1:n_cells]       # coefficient of u_j
        ab[0, 1:] = k[1:]                                               # coefficient of u_{j+1}, j = 1 .. N-1
        ab[2, :n_cells - 2] = k[1:n_cells - 1]                          # coefficient of u_{j-1}, j = 2 .. N-1
        rhs[:n_cells - 1] = driving
        # calving front, edge N (last row): 4 eta H (u_N - u_{N-1}) / dx = front
        ab[1, n_cells - 1] = 4.0 * eta[-1] * h[-1] / dx
        ab[2, n_cells - 2] = -4.0 * eta[-1] * h[-1] / dx
        rhs[n_cells - 1] = front

        u_new = np.zeros(n_cells + 1)
        u_new[1:] = solve_banded((1, 1), ab, rhs)
        change = np.max(np.abs(u_new - u)) / max(np.max(np.abs(u_new)), 1e-12)
        changes.append(change)
        u = u_new
        if change < tol:
            break
    return u, iteration, changes


# ----------------------------------------------------------------------------------
# Grounded or floating, and the grounding line (step 3)
# ----------------------------------------------------------------------------------
def thickness_above_flotation(h, x):
    """H - H_f (m): positive where the ice is grounded, negative where it floats."""
    return h - flotation_thickness(x)


def grounded_fraction_edges(h, x_centres, scheme):
    """Fraction of the friction to apply at each of the N + 1 edges (0 = floating, 1 = grounded).

    Each edge j controls the stretch of ice between the centres j - 1 and j.
    scheme = "plain": all or nothing. The edge counts as grounded if the thickness above flotation,
        averaged from its two neighbouring centres, is positive.
    scheme = "subgrid": the grounding line is placed between two centres where the thickness above
        flotation changes sign, by linear interpolation, and the edge gets the fraction of its
        stretch that is grounded (the idea of Gladstone et al. 2010 and Leguy et al. 2014).
    """
    p = thickness_above_flotation(h, x_centres)
    n_cells = len(h)
    fraction = np.zeros(n_cells + 1)
    fraction[0] = 1.0 if p[0] > 0 else 0.0                 # divide: same as the first centre
    fraction[n_cells] = 1.0 if p[-1] > 0 else 0.0          # calving front: same as the last centre
    left = p[:-1]                                          # centre j - 1 for edges j = 1 .. N-1
    right = p[1:]                                          # centre j
    if scheme == "plain":
        fraction[1:n_cells] = np.where(0.5 * (left + right) > 0, 1.0, 0.0)
    elif scheme == "subgrid":
        inner = np.where((left > 0) & (right > 0), 1.0, 0.0)
        # grounding line between the two centres: grounded part = distance from the grounded centre
        # to where p = 0, as a fraction of dx
        down = (left > 0) & (right <= 0)                   # grounded upstream, floating downstream
        up = (left <= 0) & (right > 0)                     # floating upstream, grounded downstream
        inner[down] = left[down] / (left[down] - right[down])
        inner[up] = right[up] / (right[up] - left[up])
        fraction[1:n_cells] = inner
    else:
        raise ValueError("scheme must be 'plain' or 'subgrid'")
    return fraction


def grounding_line_position(h, x_centres):
    """Grounding-line position (m): where the thickness above flotation changes from positive to
    negative for the first time going seaward, found by linear interpolation between two centres.
    Returns nan if the ice is grounded everywhere or floating everywhere."""
    p = thickness_above_flotation(h, x_centres)
    crossings = np.where((p[:-1] > 0) & (p[1:] <= 0))[0]
    if len(crossings) == 0:
        return np.nan
    i = crossings[0]
    return x_centres[i] + (x_centres[i + 1] - x_centres[i]) * p[i] / (p[i] - p[i + 1])


# ----------------------------------------------------------------------------------
# Mass conservation and time stepping (step 3)
# ----------------------------------------------------------------------------------
def edge_flux(u, h):
    """Ice flux u H (m^2/yr) at the N + 1 edges, taking H from the upstream cell (upwind).

    The velocity is positive (seaward) everywhere here, so the upstream cell of edge j is
    cell j - 1. At the divide u = 0, so the flux is 0. At the calving front (edge N) the ice that
    crosses it leaves the domain: that is the calving flux.
    """
    flux = np.zeros(len(u))
    flux[1:] = np.maximum(u[1:], 0.0) * h + np.minimum(u[1:], 0.0) * np.append(h[1:], h[-1])
    return flux


def stable_time_step(u, dx, cfl, dt_max):
    """Largest time step (yr) for which no ice moves more than cfl x dx in one step."""
    u_max = np.max(np.abs(u))
    if u_max == 0.0:
        return dt_max
    return min(cfl * dx / u_max, dt_max)


def run_to_steady_state(h0, x_centres, dx, a_yr, scheme, max_years, check_every=100.0, xg_rate_tol=1.0,
                        dhdt_tol=1e-3, cfl=0.5, dt_max=10.0, accumulation=ACCUMULATION, u_guess=None,
                        picard_tol=1e-6, keep_profiles=False, shelf_melt=0.0, h_min_shelf=1.0, melt_first_floating_cell=False,
                        verbose=False):
    """Step the ice sheet forward in time until it stops changing.

    Each time step: solve the SSA for u with the current thickness, take the largest stable time
    step, and update every cell's thickness with   H_new = H + dt (a - (flux_out - flux_in) / dx).

    Every check_every years I compare the grounding-line position and thickness with the last check.
    The run is steady when the grounding line moved by less than xg_rate_tol m/yr AND no thickness
    changed by more than dhdt_tol m/yr on average over that interval.

    Returns a dictionary with the final h and u, the time series (time, grounding line, ice volume,
    total accumulation, total calving) at every check, and whether it reached steady state.
    With keep_profiles=True it also keeps the thickness at every check ("profiles").

    shelf_melt (m/yr) melts floating ice from below (step 5). It can only remove ice down to
    h_min_shelf, so a cell never reaches zero thickness; the ice actually melted is counted in
    "melted" for the mass check. By default the first floating cell after the grounding line is not
    melted (melt_first_floating_cell=False); step 5 shows what happens if it is.
    """
    b = bed(x_centres)
    h = h0.copy()
    u = u_guess
    t = 0.0
    total_accumulated = 0.0          # m^2 of ice added by accumulation (per metre of width)
    total_calved = 0.0               # m^2 of ice that crossed the calving front
    total_melted = 0.0               # m^2 of ice melted under the floating shelf
    times, xgs, volumes, accumulated, calved = [0.0], [grounding_line_position(h, x_centres)], [np.sum(h) * dx], [0.0], [0.0]
    melted = [0.0]
    profiles = [h.copy()] if keep_profiles else []
    next_check = check_every
    h_last_check = h.copy()
    steps = 0
    steady = False
    while t < max_years:
        fraction = grounded_fraction_edges(h, x_centres, scheme)
        u, n_iter, changes = solve_ssa(h, b, dx, a_yr, fraction, u_guess=u, tol=picard_tol)
        dt = stable_time_step(u, dx, cfl, dt_max)
        if t + dt > next_check:
            dt = next_check - t                     # land exactly on the check time
        flux = edge_flux(u, h)
        h = h + dt * (accumulation - (flux[1:] - flux[:-1]) / dx)
        if shelf_melt > 0.0:
            p = thickness_above_flotation(h, x_centres)
            floating = p <= 0.0
            # no melt in the first floating cell after the grounding line: on a coarse grid its centre
            # can be almost a whole cell from the grounding line, and melting it would thin ice that
            # is still (partly) grounded (Wang et al. 2024, Seroussi and Morlighem 2018)
            crossing = np.where((p[:-1] > 0) & (p[1:] <= 0))[0]
            if len(crossing) > 0 and not melt_first_floating_cell:
                floating[crossing[0] + 1] = False
            melt = np.where(floating, np.minimum(dt * shelf_melt, np.maximum(h - h_min_shelf, 0.0)), 0.0)
            h = h - melt
            total_melted += np.sum(melt) * dx
        total_accumulated += dt * accumulation * dx * len(h)
        total_calved += dt * flux[-1]
        t += dt
        steps += 1
        if t >= next_check - 1e-9:
            xg = grounding_line_position(h, x_centres)
            times.append(t)
            xgs.append(xg)
            volumes.append(np.sum(h) * dx)
            accumulated.append(total_accumulated)
            calved.append(total_calved)
            melted.append(total_melted)
            if keep_profiles:
                profiles.append(h.copy())
            xg_rate = abs(xgs[-1] - xgs[-2]) / check_every
            dhdt = np.mean(np.abs(h - h_last_check)) / check_every
            if verbose:
                print(f"t = {t:8.0f} yr | x_g = {xg / 1e3:8.2f} km | |dx_g/dt| = {xg_rate:8.3f} m/yr | mean |dH/dt| = {dhdt:.2e} m/yr")
            h_last_check = h.copy()
            next_check += check_every
            if xg_rate < xg_rate_tol and dhdt < dhdt_tol:
                steady = True
                break
    return {"h": h, "u": u, "time": np.array(times), "x_g": np.array(xgs), "volume": np.array(volumes),
            "accumulated": np.array(accumulated), "calved": np.array(calved), "melted": np.array(melted), "steps": steps,
            "steady": steady, "years": t, "profiles": profiles}
