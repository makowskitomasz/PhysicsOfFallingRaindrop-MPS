from simulation.constants import Q_

mars = {
    "T0": Q_(290, "K"),
    "RH": 0.75,
    "p0": Q_(2e5, "Pa"),
    "g": Q_(3.71, "m/s^2"),
    "Gamma_d": Q_(9.8e-3, "K/m"),
    "rho_air": Q_(0.02, "kg/m^3"),
    "T_LCL": Q_(273 + 14.5, "K")
}
