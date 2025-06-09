Raindrop Descent Simulator
==========================

This project simulates the descent and evaporation of raindrops through planetary atmospheres.
It combines physical models, numerical methods, and visualization tools to study how droplet radius
and altitude evolve under gravity, atmospheric drag, and evaporation processes.

🌧 Project Goals
----------------
- Model the descent of raindrops through a planetary atmosphere.
- Predict evaporation based on thermodynamic and diffusion-driven processes.
- Compute terminal velocity accounting for drop deformation.
- Compare numerical results with analytical and experimental data.
- Support multiple planets like Earth, Mars, Titan, Jupiter, Saturn, and K2-18b.

🧠 Physical Assumptions and Models
----------------------------------

1. Raindrop Terminal Velocity

Terminal velocity is not constant and depends on:
- Drop size (equivalent radius)
- Drag coefficient (C_D) that is shape-dependent
- Shape deformation, represented via the shape ratio (b/a)

To find terminal velocity:
- Estimate the shape ratio from r_eq using a non-linear equation.
- Use the shape ratio to compute the surface area correction f_SA
  and then the drag shape factor C_shape.
- Use an iterative fixed-point solver to converge to the correct
  velocity v_t, since C_D itself depends on velocity.

2. Evaporation Rate

Evaporation is modeled using a modified diffusion equation:

```math
\frac{dr}{dt} = \left( \frac{f_V \cdot D \cdot M_v}{r \cdot \rho_l \cdot R} \right) \cdot \left( RH \cdot \frac{p_{\text{sat}}(T_{\text{air}})}{T_{\text{air}}} - \frac{p_{\text{sat}}(T_{\text{drop}})}{T_{\text{drop}}} \right)
```

Where:
- f_V is the ventilation factor
- D is the diffusion coefficient
- T_drop is assumed to be cooler than ambient air (based on lifting condensation level)

3. Atmospheric Profiles

- Temperature is calculated with a dry adiabatic lapse rate:
      T(z) = T0 - Gamma_d * z

- Pressure is computed assuming an isothermal atmosphere:
      p(z) = p0 * exp(-z / H)

🔬 Validation
-------------

✅ Stokes’ Law (Analytical Benchmark)

For small droplets, terminal velocity is compared with the classical formula:

```math
v_t = \frac{2 \cdot r^2 \cdot g \cdot (\rho_{\text{water}} - \rho_{\text{air}})}{9 \cdot \mu}
```

Works only for very small droplets (≲ 100 µm).

✅ Experimental Data

Numerical results are validated against tabulated experimental measurements of terminal velocity
across a range of radii.

📊 Visualization
----------------

Function: plot_radius_vs_altitude()

- Plots the evolution of radius vs. altitude for a set of initial drop sizes.
- Differentiates between surviving and evaporated drops with colors and markers.

Example:

    from simulation.plotting import plot_radius_vs_altitude
    radii0 = Q_(np.geomspace(0.12, 0.9, 6), "mm")
    z0 = Q_(600.0, "m")
    plot_radius_vs_altitude(radii0, z0, earth, simulate_raindrop_descent)

🪐 Planetary Support
---------------------

You can simulate droplet behavior in atmospheres of:

- 🌍 Earth
- ♂ Mars
- 🪐 Saturn
- 🪐 Jupiter
- 🟤 Titan
- 🌊 Exoplanet K2-18b

Each planet is defined with:

- Surface temperature T0
- Relative humidity RH
- Pressure p0
- Gravitational acceleration g
- Lapse rate Gamma_d
- Air density rho_air
- Lifting condensation level temperature T_LCL

🧪 Tests
--------

Test suite verifies:
- Unit consistency using `pint` (e.g., velocity in m/s, pressure in Pa)
- Physical correctness, such as pressure decreasing with altitude
- Impact of droplet size on evaporation and terminal velocity
- Full evaporation of small droplets before ground impact
- Survival of large droplets to the ground
- Validity of evaporation rate equation units
- Completeness and plausibility of planetary definitions (e.g. temperature, gravity, humidity)

📂 Project Structure
--------------------
```
      .
      ├── simulation/
      │   ├── model.py           # Drop descent simulation logic
      │   ├── physics.py         # Physics utilities (velocity, pressure, evaporation)
      │   ├── plotting.py        # Visualization functions
      │   ├── constants.py       # Physical constants
      ├── planets/
      │   ├── earth.py           # Planet definitions
      │   └── ...
      ├── tests/
      │   └── test_*.py          # Unit tests
      ├── notebooks/
      │   └── exploration.ipynb  # Validation and analysis
      └── README.md              # You're here!
```
📦 Requirements
----------------

- Python 3.8+
- numpy
- matplotlib
- scipy
- pint
- pytest (for running tests)

Install dependencies:

    pip install -r requirements.txt

✅ Running the Simulation
-------------------------

    from simulation.model import simulate_raindrop_descent
    from planets.earth import earth
    from simulation.constants import Q_

    r0 = Q_(0.3, "mm")
    z0 = Q_(600, "m")
    trajectory = simulate_raindrop_descent(r0, z0, earth)

📈 Example Plot
----------------




<p align="center">
  <img src="https://github.com/user-attachments/assets/fd13a570-9a54-4e5e-838f-b9a611eb1c90" alt="Raindrop Simulation Plot" width="600"/><br/>
  <em>Figure: Droplet radius vs. altitude for various initial sizes.</em>
</p>


📘 License
----------

MIT License. See LICENSE file for details.

