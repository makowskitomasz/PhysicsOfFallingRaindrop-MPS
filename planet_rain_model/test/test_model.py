import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from simulation.model import simulate_raindrop_descent
from simulation.constants import Q_
from planets.earth import earth

def test_evaporated_drop_stops_early():
    """
    Test that a very small droplet evaporates before reaching the ground.
    """
    r0 = Q_(0.01, "mm").to("m")
    z0 = Q_(600, "m")
    traj = simulate_raindrop_descent(r0, z0, earth)
    final_z = traj[-1][0]
    final_r = traj[-1][1]
    assert final_z > 0
    assert final_r < 1e-8

def test_large_drop_survives_to_ground():
    """
    Test that a large droplet reaches the ground with non-zero radius.
    """
    r0 = Q_(0.5, "mm").to("m")
    z0 = Q_(600, "m")
    traj = simulate_raindrop_descent(r0, z0, earth)
    final_z = traj[-1][0]
    final_r = traj[-1][1]
    assert final_z <= 1.0
    assert final_r > 1e-8
