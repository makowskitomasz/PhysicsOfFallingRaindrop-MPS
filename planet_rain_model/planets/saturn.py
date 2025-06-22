from simulation.constants import Q_

saturn = {
    "T0": Q_(284, "K"),
    "RH": 1.0,
    "p0": Q_(1.04e6, "Pa"),
    "g": Q_(10.47, "m/s^2"),
    "Gamma_d": Q_(9.8e-3, "K/m"),
    "rho_air": Q_(0.18, "kg/m^3"),
    "T_LCL": Q_(273 + 99.2, "K")
}
