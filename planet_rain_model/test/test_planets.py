import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from planets.earth import earth
from planets.mars import mars
from planets.titan import titan
from planets.jupiter import jupiter
from planets.saturn import saturn
from planets.k2_18b import k2_18b

def test_earth_parameters():
    """
    Verify Earth's basic parameters: temperature, gravity, and relative humidity.
    """
    assert earth['T0'].magnitude == 300
    assert 0 < earth['RH'] <= 1
    assert earth['g'].magnitude == 9.82

def test_mars_gravity_and_pressure():
    """
    Check Mars' surface gravity and pressure are set to expected values.
    """
    assert mars['g'].magnitude == 3.71
    assert mars['p0'].magnitude == 2e5

def test_all_planets_have_required_fields():
    """
    Ensure all defined planets contain the required atmospheric parameters.
    """
    planets = [earth, mars, titan, jupiter, saturn, k2_18b]
    required_fields = ["T0", "RH", "p0", "g", "Gamma_d", "rho_air", "T_LCL"]
    
    for planet in planets:
        for field in required_fields:
            assert field in planet, f"Missing field {field} in planet {planet.get('name', '?')}"
