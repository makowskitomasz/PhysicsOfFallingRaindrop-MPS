from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

# Physical constants
R = Q_(8.314, "J/(mol*K)")  # Universal gas constant
M_v = Q_(18e-3, "kg/mol")  # Molar mass of water
rho_water = Q_(1000, "kg/m^3")  # Density of liquid water
eta_air = Q_(1.8e-5, "kg/(m*s)")  # Dynamic viscosity of air
D_vap = Q_(2.5e-5, "m^2/s")  # Diffusion coefficient of water vapor in air

sigma_water_air = Q_(0.073, "N/m")  # Surface tension between water and air
gravity = Q_(9.81, "m/s^2")  # Acceleration due to gravity
rho_air = Q_(1.205, "kg/m^3")  # Density of air
air_viscosity = Q_(1.81e-5, "kg/(m*s)")  # Viscosity of air