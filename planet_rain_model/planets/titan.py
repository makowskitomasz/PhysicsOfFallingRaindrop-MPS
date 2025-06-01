from simulation.constants import Q_

titan = {
    "T0": Q_(90, "K"),
    "RH": 0.75,
    "p0": Q_(1.5e5, "Pa"),
    "g": Q_(1.35, "m/s^2"),
    "Gamma_d": Q_(9.8e-3, "K/m"),
    "rho_air": Q_(5.4, "kg/m^3"),
    "T_LCL": Q_(273 + 18, "K")
}
