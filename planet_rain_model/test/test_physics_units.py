from simulation.physics import (
    saturation_vapor_pressure,
    temperature_dry_adiabatic,
    pressure_profile,
    terminal_velocity,
    evaporation_rate_exact,
    ventilation_factor
)
from simulation.constants import Q_, eta_air, D_vap
from planets.earth import earth

import pint


def test_terminal_velocity_unit():
    """
    Ensure terminal_velocity returns a pint Quantity with correct velocity units.
    """
    r = Q_(1e-3, "m")
    v = terminal_velocity(r)
    assert isinstance(v, pint.Quantity)
    assert v.check("[length] / [time]")


def test_saturation_vapor_pressure_units():
    """
    Ensure saturation_vapor_pressure returns pressure in proper units.
    """
    T = Q_(300, "K")
    p_sat = saturation_vapor_pressure(T)
    assert isinstance(p_sat, pint.Quantity)
    assert p_sat.check("[pressure]")


def test_temperature_profile_unit():
    """
    Check that dry adiabatic temperature profile has temperature units.
    """
    z = Q_(500, "m")
    T = temperature_dry_adiabatic(z, earth["T0"], earth["Gamma_d"])
    assert T.check("[temperature]")


def test_pressure_profile_unit():
    """
    Check that pressure_profile returns a quantity with pressure units.
    """
    z = Q_(1000, "m")
    p = pressure_profile(z, earth["p0"], earth["T0"], earth["g"])
    assert p.check("[pressure]")


def test_evaporation_rate_unit():
    """
    Ensure evaporation_rate_exact returns a rate with correct units (m/s).
    """
    r = Q_(0.0005, "m")
    z = Q_(100, "m")
    T = temperature_dry_adiabatic(z, earth["T0"], earth["Gamma_d"])
    p = pressure_profile(z, earth["p0"], earth["T0"], earth["g"])
    v = terminal_velocity(r)
    Re = (earth["rho_air"] * v * 2 * r) / eta_air
    Sc = eta_air / (earth["rho_air"] * D_vap)
    f_V = ventilation_factor(Re.magnitude, Sc.magnitude)

    drdt = evaporation_rate_exact(r, T, RH=earth["RH"], p=p, f_V=f_V, T_LCL=earth["T_LCL"])
    assert drdt.check("[length] / [time]")
