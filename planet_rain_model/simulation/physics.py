import numpy as np
from .constants import *

def saturation_vapor_pressure(T):
    T_C = T.to("degC").magnitude
    p_sat = 610.94 * np.exp((17.625 * T_C) / (T_C + 243.04))
    return Q_(p_sat, "Pa")

def temperature_dry_adiabatic(z, T0, Gamma_d):
    return T0 - Gamma_d * z

def pressure_profile(z, p0, T0, g):
    H = (R * T0) / (M_v * g)
    return p0 * np.exp(-z / H)

def terminal_velocity(r_eq):
    a = Q_(130, "m**0.5 / s")
    return (a * r_eq**0.5).to("m/s")

def ventilation_factor(Re, Sc):
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
