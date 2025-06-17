from .physics import temperature_dry_adiabatic, pressure_profile, terminal_velocity, ventilation_factor, evaporation_rate_exact
from .constants import Q_, D_vap, eta_air

def simulate_raindrop_descent(r0, z0, planet, dz=Q_(0.1, "m")):
    """
    Simulates a raindrop's descent and evaporation through the atmosphere.

    Parameters
    ----------
    r0 : pint.Quantity
        Initial drop radius.
    z0 : pint.Quantity
        Initial height.
    planet : dict
        Planetary parameters (T0, Gamma_d, p0, g, rho_air, RH, T_LCL).
    dz : pint.Quantity, optional
        Step size in height (default: 0.1 m).

    Returns
    -------
    list of tuple
        List of (altitude [m], radius [m]) during descent.
    """
    r, z = r0, z0
    traj = []
    while r > Q_(1e-6, "m") and z > Q_(0, "m"):
        T = temperature_dry_adiabatic(z, planet["T0"], planet["Gamma_d"])
        p = pressure_profile(z, planet["p0"], planet["T0"], planet["g"])
        v = terminal_velocity(r)

        Re = (planet["rho_air"] * v * 2 * r) / eta_air
        Sc = eta_air / (planet["rho_air"] * D_vap)
        f_V = ventilation_factor(Re.magnitude, Sc.magnitude)

        drdt = evaporation_rate_exact(r, T, planet["RH"], p, f_V, planet["T_LCL"])
        drdz = drdt / v
        r = (r - drdz * dz).to("m")
        z = z - dz
        traj.append((z.magnitude, r.magnitude))
    return traj
