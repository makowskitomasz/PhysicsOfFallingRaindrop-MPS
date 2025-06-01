from simulation.constants import Q_

earth = {
    "T0": Q_(300.0, "K"),
    "RH": 0.75,
    "p0": Q_(101325, "Pa"),
    "g": Q_(9.82, "m/s^2"),
    "Gamma_d": Q_(9.8e-3, "K/m"),
    "rho_air": Q_(1.225, "kg/m^3"),
    "T_LCL": Q_(275, "K")
}
