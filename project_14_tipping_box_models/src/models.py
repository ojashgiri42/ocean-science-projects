
import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------------
# Stommel (1961) two-box model, dimensionless form
#   T = temperature difference, S = salinity difference (scaled)
#   f = flow strength; f < 0 is temperature-driven, f > 0 salinity-driven
# ---------------------------------------------------------------------

def stommel_flow(T, S, R, lam):
    """Flow f = (-T + R S) / lambda. The density difference drives the flow."""
    return (-T + R * S) / lam


def stommel_tendency(T, S, R, delta, lam):
    """Return (dT/dt, dS/dt).

    1 - T and delta*(1 - S) pull each box toward the forcing value 1.
    |f|*T and |f|*S are mixing by the flow, which works the same in either
    direction, so only the speed |f| matters here.
    """
    f = stommel_flow(T, S, R, lam)
    dTdt = 1.0 - T - abs(f) * T
    dSdt = delta * (1.0 - S) - abs(f) * S
    return dTdt, dSdt


def stommel_equilibrium_mismatch(f, R, delta, lam):
    """Zero exactly at an equilibrium.

    At steady state T = 1/(1+|f|) and S = delta/(delta+|f|). Putting these
    into lambda*f = -T + R*S gives one equation in f alone.
    """
    rhs = -1.0 / (1.0 + abs(f)) + R * delta / (delta + abs(f))
    return rhs - lam * f


def stommel_find_equilibria(R, delta, lam, f_min, f_max, n_scan):
    """All equilibrium flows f in [f_min, f_max], found as in step 1.

    Scan for sign changes of the mismatch on a grid, then refine each one
    with brentq.
    """
    f_scan = np.linspace(f_min, f_max, n_scan)
    mismatch = np.zeros(n_scan)
    for i in range(n_scan):
        mismatch[i] = stommel_equilibrium_mismatch(f_scan[i], R, delta, lam)

    roots = []
    for i in range(n_scan - 1):
        if mismatch[i] * mismatch[i + 1] < 0:
            root = brentq(stommel_equilibrium_mismatch, f_scan[i], f_scan[i + 1],
                          args=(R, delta, lam), xtol=1e-12)
            roots.append(root)
    return roots


def stommel_state_from_f(f, delta):
    """Steady-state T and S for a given equilibrium flow f."""
    T = 1.0 / (1.0 + abs(f))
    S = delta / (delta + abs(f))
    return T, S


def stommel_jacobian(T, S, R, delta, lam, h=1e-6):
    """2x2 Jacobian by central finite differences.

    We do it numerically because |f| makes the algebra messy, and a small
    centred step is accurate enough away from f = 0.
    """
    J = np.zeros((2, 2))

    dT_plus, dS_plus = stommel_tendency(T + h, S, R, delta, lam)
    dT_minus, dS_minus = stommel_tendency(T - h, S, R, delta, lam)
    J[0, 0] = (dT_plus - dT_minus) / (2 * h)
    J[1, 0] = (dS_plus - dS_minus) / (2 * h)

    dT_plus, dS_plus = stommel_tendency(T, S + h, R, delta, lam)
    dT_minus, dS_minus = stommel_tendency(T, S - h, R, delta, lam)
    J[0, 1] = (dT_plus - dT_minus) / (2 * h)
    J[1, 1] = (dS_plus - dS_minus) / (2 * h)

    return J


def classify_eigenvalues(eigenvalues):
    """Name the equilibrium type from the two eigenvalues of the Jacobian."""
    re = np.real(eigenvalues)
    im = np.imag(eigenvalues)
    is_complex = abs(im[0]) > 1e-10

    if is_complex:
        if re[0] < 0:
            return "stable spiral"
        return "unstable spiral"
    if re[0] < 0 and re[1] < 0:
        return "stable node"
    if re[0] > 0 and re[1] > 0:
        return "unstable node"
    return "saddle"


def stommel_run_euler(T0, S0, R, delta, lam, dt, t_end):
    """Forward Euler from (T0, S0). Returns arrays of T and S over time."""
    n_steps = int(round(t_end / dt))
    T = np.zeros(n_steps + 1)
    S = np.zeros(n_steps + 1)
    T[0] = T0
    S[0] = S0
    for i in range(n_steps):
        dTdt, dSdt = stommel_tendency(T[i], S[i], R, delta, lam)
        T[i + 1] = T[i] + dt * dTdt
        S[i + 1] = S[i] + dt * dSdt
    return T, S


# ---------------------------------------------------------------------
# Cessi (1994) one-equation model
#   y = scaled salinity difference; small y = strong overturning
#   time is in units of the diffusion time (about 219 years)
# ---------------------------------------------------------------------

def cessi_tendency(y, p, mu2):
    """dy/dt = p - y (1 + mu2 (1 - y)^2).

    p is the freshwater forcing that builds up the salinity difference.
    y * 1 is slow diffusion. y * mu2 (1 - y)^2 is the salt carried by the
    overturning: the flow is proportional to the density difference, 1 - y
    (temperature is fixed at 1), and Cessi made the exchange rate grow with
    the square of the flow. Squaring, instead of Stommel's |f|, keeps the
    equation smooth, which is why it has a potential V(y).
    """
    return p - y * (1.0 + mu2 * (1.0 - y) ** 2)


def cessi_p_of_y(y, mu2):
    """Forcing p needed to hold the system steady at y (set dy/dt = 0)."""
    return y * (1.0 + mu2 * (1.0 - y) ** 2)


def cessi_dpdy(y, mu2):
    """Slope of p(y). An equilibrium is stable where dp/dy > 0.

    The tendency is p - p(y), so its slope at an equilibrium is -dp/dy.
    A positive dp/dy means a small nudge decays.
    """
    return 1.0 + mu2 * (1.0 - 4.0 * y + 3.0 * y ** 2)


def cessi_potential(y, p, mu2):
    """V(y) with dy/dt = -dV/dy. Valleys are stable states, hills unstable."""
    return mu2 * (y ** 4 / 4.0 - 2.0 * y ** 3 / 3.0 + y ** 2 / 2.0) + y ** 2 / 2.0 - p * y


def cessi_equilibria(p, mu2):
    """Real equilibria at forcing p, sorted from small to large y.

    Setting dy/dt = 0 gives a cubic: mu2 y^3 - 2 mu2 y^2 + (1 + mu2) y - p = 0,
    so np.roots finds all of them at once.
    """
    roots = np.roots([mu2, -2.0 * mu2, 1.0 + mu2, -p])
    real_roots = []
    for r in roots:
        if abs(r.imag) < 1e-9:
            real_roots.append(r.real)
    real_roots.sort()
    return real_roots


def cessi_fold_points(mu2):
    """The two y values where dp/dy = 0, and the p at each.

    dp/dy = 0 is the quadratic 3 mu2 y^2 - 4 mu2 y + (1 + mu2) = 0.
    """
    a = 3.0 * mu2
    b = -4.0 * mu2
    c = 1.0 + mu2
    disc = np.sqrt(b ** 2 - 4.0 * a * c)
    y_low = (-b - disc) / (2.0 * a)
    y_high = (-b + disc) / (2.0 * a)
    return y_low, cessi_p_of_y(y_low, mu2), y_high, cessi_p_of_y(y_high, mu2)


def cessi_run_ramp(y0, p_values, mu2, dt):
    """Forward Euler with p changing at every step. Returns y at each step."""
    n = len(p_values)
    y = np.zeros(n)
    y[0] = y0
    for i in range(n - 1):
        y[i + 1] = y[i] + dt * cessi_tendency(y[i], p_values[i], mu2)
    return y


# ---------------------------------------------------------------------
# Toy Antarctic ice-shelf margin (step 6). Not from a paper.
#   T = subsurface temperature above freezing (deg C)
#   F = surface freshening (psu); time in years
# The parameters are passed in one dictionary, because there are many.
# ---------------------------------------------------------------------

def antarctic_kappa(F, par):
    """Vertical mixing rate (1/yr). Freshening stabilizes the surface and weakens it.

    Two shapes are allowed, because the answer depends on which one we pick:
    "rational":    kappa0 / (1 + F/F0)
    "exponential": kappa0 * exp(-F/F0)
    """
    if par["form"] == "rational":
        return par["kappa0"] / (1.0 + F / par["F0"])
    if par["form"] == "exponential":
        return par["kappa0"] * np.exp(-F / par["F0"])
    raise ValueError("form must be 'rational' or 'exponential'")


def antarctic_melt(T, par):
    """Basal melt rate (m of ice per year), quadratic in temperature above freezing."""
    T_pos = max(T, 0.0)    # no melting (and no freezing) below the freezing point in this toy model
    return par["m"] * T_pos ** 2


def antarctic_tendency(T, F, T_far, par):
    """Return (dT/dt, dF/dt).

    dT/dt: inflow pulls T toward the far-field value, vertical mixing loses
    heat to the cold surface (T = 0 there), and melting uses heat.
    dF/dt: meltwater freshens the surface, and the freshening is removed
    over a time tau_F.
    """
    melt = antarctic_melt(T, par)
    dTdt = (T_far - T) / par["tau_in"] - antarctic_kappa(F, par) * T - par["b"] * melt
    dFdt = par["c"] * melt - F / par["tau_F"]
    return dTdt, dFdt


def antarctic_Tfar_of_T(T, par):
    """Far-field temperature needed for a steady state at subsurface temperature T.

    Steady F equation: F = c * tau_F * melt(T). Put that into the steady
    T equation and solve for T_far. Works on arrays, like cessi_p_of_y.
    """
    T_pos = np.maximum(T, 0.0)
    melt = par["m"] * T_pos ** 2
    F = par["c"] * par["tau_F"] * melt
    if par["form"] == "rational":
        kappa = par["kappa0"] / (1.0 + F / par["F0"])
    else:
        kappa = par["kappa0"] * np.exp(-F / par["F0"])
    return T + par["tau_in"] * (kappa * T + par["b"] * melt)


def antarctic_jacobian(T, F, T_far, par, h=1e-6):
    """2x2 Jacobian of (dT/dt, dF/dt) by central finite differences."""
    J = np.zeros((2, 2))
    a_plus, b_plus = antarctic_tendency(T + h, F, T_far, par)
    a_minus, b_minus = antarctic_tendency(T - h, F, T_far, par)
    J[0, 0] = (a_plus - a_minus) / (2 * h)
    J[1, 0] = (b_plus - b_minus) / (2 * h)
    a_plus, b_plus = antarctic_tendency(T, F + h, T_far, par)
    a_minus, b_minus = antarctic_tendency(T, F - h, T_far, par)
    J[0, 1] = (a_plus - a_minus) / (2 * h)
    J[1, 1] = (b_plus - b_minus) / (2 * h)
    return J
