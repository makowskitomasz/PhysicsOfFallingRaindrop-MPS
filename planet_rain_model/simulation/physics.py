from scipy.optimize import brentq, fixed_point
import numpy as np
from .constants import *

def saturation_vapor_pressure(T):
    """
    Calculates the saturation vapor pressure for a given temperature.

    Parameters:
        T (Quantity): Temperature in Kelvin.

    Returns:
        Quantity: Saturation vapor pressure in Pascals.
    """
    T_C = T.to("degC").magnitude
    p_sat = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))
    return Q_(p_sat, "Pa")

def temperature_dry_adiabatic(z, T0, Gamma_d):
    """
    Computes the temperature at a given altitude assuming a dry adiabatic lapse rate.

    Parameters:
        z (Quantity): Altitude in meters.
        T0 (Quantity): Surface temperature in Kelvin.
        Gamma_d (Quantity): Dry adiabatic lapse rate in Kelvin/meter.

    Returns:
        Quantity: Temperature at altitude z in Kelvin.
    """
    return T0 - Gamma_d * z

def pressure_profile(z, p0, T0, g):
    """
    Computes the pressure at a given altitude using the barometric formula.

    Parameters:
        z (Quantity): Altitude in meters.
        p0 (Quantity): Surface pressure in Pascals.
        T0 (Quantity): Surface temperature in Kelvin.
        g (Quantity): Gravitational acceleration in m/s^2.

    Returns:
        Quantity: Pressure at altitude z in Pascals.
    """
    H = (R * T0) / (M_v * g)
    return p0 * np.exp(-z / H)


def calculate_shape_ratio(r_eq):
    """
    Calculates the shape ratio of a raindrop based on its equivalent radius.

    Parameters:
        r_eq (Quantity): Equivalent radius of the raindrop in meters.

    Returns:
        Quantity: Shape ratio (dimensionless).
    """
    def F(x):
        # Unfortunately, scipy does not work well with pint units, forcing us to use
        # dimensionless quantities
        return np.sqrt(sigma_water_air.magnitude / (gravity.magnitude * (rho_water.magnitude - rho_air.magnitude))) * x**(-1/6) * \
               np.sqrt(x**(-2) - 2 * x**(-1/3) + 1) - r_eq.magnitude

    return Q_(brentq(F, 1e-9, 1.0 - 1e-9), "dimensionless")

def calculate_fSA(shape_ratio):
    """
    Computes the fSA parameter based on the shape ratio.

    Parameters:
        shape_ratio (float): Shape ratio (dimensionless).

    Returns:
        float: fSA parameter (dimensionless).
    """
    epsilon = np.sqrt(1 - shape_ratio**2)
    if shape_ratio < 1:
        return 0.5 * (shape_ratio ** (-2/3)) + shape_ratio ** (4/3) * \
            np.log((1 + epsilon) / (1 - epsilon)) / (4 * epsilon)
    elif shape_ratio == 1:
        return 1.0
    else:
        raise ValueError("Shape ratio must be in (0, 1]")

def calculate_C_shape(fSA):
    """
    Computes the C_shape parameter based on the fSA parameter.

    Parameters:
        fSA (float): fSA parameter (dimensionless).

    Returns:
        float: C_shape parameter (dimensionless).
    """
    return 1 + 1.5 * (fSA - 1)**0.5 + 6.7 * (fSA - 1)

def calculate_CD(C_shape, vT, r_eq):
    """
    Computes the drag coefficient (Cd) for a raindrop.

    Parameters:
        C_shape (float): C_shape parameter (dimensionless).
        vT (Quantity): Terminal velocity in m/s.
        r_eq (Quantity): Equivalent radius in meters.

    Returns:
        float: Drag coefficient (dimensionless).
    """
    Re = vT * 2 * r_eq * rho_air / air_viscosity
    return ((24 / Re) * (1 + 0.15 * Re**0.687) + 0.42 * (1 + 4.25 * 10**4 * Re**-1.16)**-1) * C_shape

def calculate_terminal_velocity(Cd, shape_ratio, r_eq):
    """
    Computes the terminal velocity of a raindrop.

    Parameters:
        Cd (float): Drag coefficient (dimensionless).
        shape_ratio (float): Shape ratio (dimensionless).
        r_eq (Quantity): Equivalent radius in meters.

    Returns:
        Quantity: Terminal velocity in m/s.
    """
    return np.sqrt(8/3 * ((rho_water - rho_air) / rho_air) * (gravity / Cd) * shape_ratio**(2/3) * r_eq)

def find_terminal_velocity(r_eq, C_shape, shape_ratio):
    """
    Finds the terminal velocity of a raindrop using a fixed-point method.

    Parameters:
        r_eq (Quantity): Equivalent radius in meters.
        C_shape (float): C_shape parameter (dimensionless).
        shape_ratio (float): Shape ratio (dimensionless).

    Returns:
        Quantity: Terminal velocity in m/s.
    """
    v0 = 0.001
    def f(vT):
        # Unfortunately, scipy does not work well with pint units, forcing us to use
        # dimensionless quantities and converting them back later on
        vT = np.abs(vT)
        Cd = calculate_CD(C_shape, Q_(vT, "m/s"), r_eq)
        return calculate_terminal_velocity(Cd, shape_ratio, r_eq).magnitude

    return Q_(fixed_point(f, v0), "m/s")

def terminal_velocity(r_eq):
    """
    Computes the terminal velocity of a raindrop based on its equivalent radius.

    Parameters:
        r_eq (Quantity): Equivalent radius in meters.

    Returns:
        Quantity: Terminal velocity in m/s.
    """
    shape_ratio = calculate_shape_ratio(r_eq)
    fSA = calculate_fSA(shape_ratio)
    C_shape = calculate_C_shape(fSA)
    v_terminal = find_terminal_velocity(r_eq, C_shape, shape_ratio)
    return v_terminal

def ventilation_factor(Re, Sc):
    """
    Computes the ventilation factor for a raindrop.

    Parameters:
        Re (float): Reynolds number (dimensionless).
        Sc (float): Schmidt number (dimensionless).

    Returns:
        float: Ventilation factor (dimensionless).
    """
    if Re < 1e-6:
        return 1.0
    elif Re < 2000:
        return 1 + 0.108 * Re**0.5 * Sc**(1/3)
    else:
        return 0.78 * Re**0.308

def evaporation_rate_exact(r_eq, T_air, RH, p, f_V, T_LCL):
    """
    Computes the evaporation rate of a raindrop.

    Parameters:
        r_eq (Quantity): Equivalent radius in meters.
        T_air (Quantity): Air temperature in Kelvin.
        RH (float): Relative humidity (dimensionless).
        p (Quantity): Pressure in Pascals.
        f_V (float): Ventilation factor (dimensionless).
        T_LCL (Quantity): Temperature at the lifting condensation level in Kelvin.

    Returns:
        Quantity: Evaporation rate in m/s.
    """
    T_drop = T_air - 0.5 * (T_air - T_LCL)
    p_sat_air = saturation_vapor_pressure(T_air)
    p_sat_drop = saturation_vapor_pressure(T_drop)

    delta = RH * (p_sat_air / T_air) - (p_sat_drop / T_drop)

    numerator = f_V * D_vap * M_v
    denominator = r_eq * rho_water * R

    term1 = numerator / denominator

    return (term1 * delta).to("m/s")
