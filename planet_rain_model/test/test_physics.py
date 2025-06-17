import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from simulation.physics import saturation_vapor_pressure, terminal_velocity, pressure_profile
from simulation.constants import Q_

def test_saturation_pressure_reasonable():
    """
    Test that saturation vapor pressure at 300 K falls within a physically realistic range.
    """
    T = Q_(300, "K")
    p_sat = saturation_vapor_pressure(T)
    assert 3000 < p_sat.magnitude < 5000 

def test_terminal_velocity_scaling():
    """
    Ensure that larger drops fall faster than smaller ones and terminal velocity is positive.
    """
    r_small = Q_(1e-4, "m")
    r_big = Q_(1e-3, "m")
    v_small = terminal_velocity(r_small)
    v_big = terminal_velocity(r_big)
    assert v_big > v_small
    assert v_small.magnitude > 0

def test_pressure_decreases_with_height():
    """
    Confirm that atmospheric pressure decreases with increasing altitude.
    """
    z_low = Q_(0, "m")
    z_high = Q_(5000, "m")
    p0 = Q_(101325, "Pa")
    T0 = Q_(300, "K")
    g = Q_(9.82, "m/s^2")
    p1 = pressure_profile(z_low, p0, T0, g)
    p2 = pressure_profile(z_high, p0, T0, g)
    assert p2 < p1
