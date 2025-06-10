from .physics import *

def simulate_raindrop_descent(r0, z0, planet, dz=Q_(0.1, "m")):
    """
    Simulates the descent of a raindrop through the atmosphere.

    Parameters:
        r0 (Quantity): Initial radius of the raindrop (in meters).
        z0 (Quantity): Initial altitude of the raindrop (in meters).
        planet (dict): Dictionary containing planetary parameters:
            - T0 (Quantity): Surface temperature (in Kelvin).
            - Gamma_d (Quantity): Dry adiabatic lapse rate (in Kelvin/meter).
            - p0 (Quantity): Surface pressure (in Pascals).
            - g (Quantity): Gravitational acceleration (in m/s^2).
            - rho_air (Quantity): Density of air (in kg/m^3).
            - RH (float): Relative humidity (as a fraction).
            - T_LCL (Quantity): Temperature at the lifting condensation level (in Kelvin).
        dz (Quantity): Altitude step for the simulation (default is 0.1 meters).

    Returns:
        list: A trajectory list containing tuples of altitude (z) and radius (r) at each step.
    """
    r, z = r0, z0
    traj = []
    while r > Q_(1e-6, "m") and z > Q_(0, "m"):
        # Compute temperature and pressure at the current altitude
        T = temperature_dry_adiabatic(z, planet["T0"], planet["Gamma_d"])
        p = pressure_profile(z, planet["p0"], planet["T0"], planet["g"])
        
        # Calculate terminal velocity of the raindrop
        v = terminal_velocity(r)

        # Compute Reynolds and Schmidt numbers
        Re = (planet["rho_air"] * v * 2 * r) / eta_air
        Sc = eta_air / (planet["rho_air"] * D_vap)
        
        # Calculate ventilation factor
        f_V = ventilation_factor(Re.magnitude, Sc.magnitude)

        # Compute evaporation rate and update radius and altitude
        drdt = evaporation_rate_exact(r, T, planet["RH"], p, f_V, planet["T_LCL"])
        drdz = drdt / v
        r = (r - drdz * dz).to("m")
        z = z - dz
        
        # Append current state to trajectory
        traj.append((z.magnitude, r.magnitude))
    return traj
