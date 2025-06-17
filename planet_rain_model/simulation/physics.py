from scipy.optimize import brentq, fixed_point
import numpy as np
from .constants import Q_, R, M_v, rho_air, rho_water, sigma_water_air, gravity, air_viscosity, D_vap

def saturation_vapor_pressure(T):
    """
    Returns saturation vapor pressure over liquid water at temperature T.
    """
    T_C = T.to("degC").magnitude
    p_sat = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))
    return Q_(p_sat, "Pa")

def temperature_dry_adiabatic(z, T0, Gamma_d):
    """
    Returns temperature at altitude z using dry adiabatic lapse rate.
    """
    return T0 - Gamma_d * z

def pressure_profile(z, p0, T0, g):
    """
    Returns atmospheric pressure at altitude z assuming isothermal scale height.
    """
    H = (R * T0) / (M_v * g)
    return p0 * np.exp(-z / H)


def calculate_shape_ratio(r_eq):
    """
    Solves for drop shape ratio based on equilibrium radius.
    """
    def F(x):
        # Unfortunately, scipy does not work well with pint units, forcing us to use
        # dimensionless quantities
        return np.sqrt(sigma_water_air.magnitude / (gravity.magnitude * (rho_water.magnitude - rho_air.magnitude))) * x**(-1/6) * \
               np.sqrt(x**(-2) - 2 * x**(-1/3) + 1) - r_eq.magnitude

    return Q_(brentq(F, 1e-9, 1.0 - 1e-9), "dimensionless")

def calculate_fSA(shape_ratio):
    """
    Computes surface area correction factor fSA based on shape ratio.
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
    Computes drag shape correction factor from fSA.
    """
    return 1 + 1.5 * (fSA - 1)**0.5 + 6.7 * (fSA - 1)

def calculate_CD(C_shape, vT, r_eq):
    """
    Calculates drag coefficient (CD) based on velocity, shape, and radius.
    """
    Re = vT * 2 * r_eq * rho_air / air_viscosity
    return ((24 / Re) * (1 + 0.15 * Re**0.687) + 0.42 * (1 + 4.25 * 10**4 * Re**-1.16)**-1) * C_shape

def calculate_terminal_velocity(Cd, shape_ratio, r_eq):
    """
    Calculates terminal velocity from drag coefficient, shape and radius.
    """
    return np.sqrt(8/3 * ((rho_water - rho_air) / rho_air) * (gravity / Cd) * shape_ratio**(2/3) * r_eq)

def find_terminal_velocity(r_eq, C_shape, shape_ratio):
    """
    Iteratively solves for terminal velocity of a drop.
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
    Computes terminal velocity of a drop with given equivalent radius.
    """
    shape_ratio = calculate_shape_ratio(r_eq)
    fSA = calculate_fSA(shape_ratio)
    C_shape = calculate_C_shape(fSA)
    v_terminal = find_terminal_velocity(r_eq, C_shape, shape_ratio)
    return v_terminal

def ventilation_factor(Re, Sc):
    """
    Returns ventilation factor based on Reynolds (Re) and Schmidt (Sc) numbers.
    """
    if Re < 1e-6:
        return 1.0
    elif Re < 2000:
        return 1 + 0.108 * Re**0.5 * Sc**(1/3)
    else:
        return 0.78 * Re**0.308

def evaporation_rate_exact(r_eq, T_air, RH, p, f_V, T_LCL):
    """
    Computes dr/dt using Eq. from Figure 1:
    dr/dt = (f_V * D * M_v) / (r * rho_l * R) * (RH * p_sat(T_air)/T_air - p_sat(T_drop)/T_drop)
    """
    T_drop = T_air - 0.5 * (T_air - T_LCL)
    p_sat_air = saturation_vapor_pressure(T_air)
    p_sat_drop = saturation_vapor_pressure(T_drop)

    delta = RH * (p_sat_air / T_air) - (p_sat_drop / T_drop)

    numerator = f_V * D_vap * M_v
    denominator = r_eq * rho_water * R

    term1 = numerator / denominator

    return (term1 * delta).to("m/s")
