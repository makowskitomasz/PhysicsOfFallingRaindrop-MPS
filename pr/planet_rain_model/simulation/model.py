from .physics import *

def simulate_raindrop_descent(r0, z0, planet, dz=Q_(0.1, "m")):
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
