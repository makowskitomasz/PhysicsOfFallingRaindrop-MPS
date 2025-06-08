from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

# Physical constants
R = Q_(8.314, "J/(mol*K)")
M_v = Q_(18e-3, "kg/mol")
rho_water = Q_(1000, "kg/m^3")
eta_air = Q_(1.8e-5, "kg/(m*s)")
D_vap = Q_(2.5e-5, "m^2/s")

sigma_water_air = Q_(0.073, "N/m")
gravity = Q_(9.81, "m/s^2")
rho_air = Q_(1.205, "kg/m^3")
air_viscosity = Q_(1.81e-5, "kg/(m*s)")