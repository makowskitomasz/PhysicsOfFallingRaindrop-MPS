# Raindrop Descent Simulator


This project simulates the descent and evaporation of raindrops through planetary atmospheres.
It combines physical models, numerical methods, and visualization tools to study how droplet radius
and altitude evolve under gravity, atmospheric drag, and evaporation processes.

## Project Goals

Goal: Compare the results of the figure and the achieved values.

Steps:
* Reconstruct the initial values $r_0$ as shown in the figure.

* Check if the droplets exceed the threshold value $​r_{min}$.

* Compare the final values $r_{eq} (z = 0)$ with the results presented in the article.

## Physical Assumptions and Models

1. Raindrop Terminal Velocity

Terminal velocity is not constant and depends on:
- Drop size (equivalent radius)
- Drag coefficient ($C_D$) that is shape-dependent
- Shape deformation, represented via the shape ratio ($\frac{b}{a}$)

To find terminal velocity:
- Estimate the shape ratio from $r_{eq}$ using a non-linear equation.
- Use the shape ratio to compute the surface area correction $f_{SA}$
  and then the drag shape factor $C_{shape}$.
- Use an iterative fixed-point solver to converge to the correct
  velocity $v_t$, since $C_D$ itself depends on velocity.

2. Evaporation Rate

Evaporation is modeled using a modified diffusion equation:

```math
\frac{dr}{dt} = \left( \frac{f_V \cdot D \cdot M_v}{r \cdot \rho_l \cdot R} \right) \cdot \left( RH \cdot \frac{p_{\text{sat}}(T_{\text{air}})}{T_{\text{air}}} - \frac{p_{\text{sat}}(T_{\text{drop}})}{T_{\text{drop}}} \right)
```

Where:
- $f_V$ is the ventilation factor
- D is the diffusion coefficient
- $T_{drop}$ is assumed to be cooler than ambient air (based on lifting condensation level)


## Raindrop shape
"In isolation, a precipitating particle does two things: (1) it falls and (2) it evaporates. To calculate the rate at which a particle falls and evaporates requires knowledge of the relationship between particle mass and shape. Unlike solid precipitating particles, whose forms vary widely, raindrops have equilibrium shapes that can be uniquely calculated for a given mass of liquid condensible, air density, and surface gravity. A unique shape allows us, in a known external environment, to describe a raindrop with only a single size variable. Here, we use equivalent radius $r_{eq}$" - The Physics of Falling Raindrops in Diverse Planetary Atmospheres

In order to calculate terminal velocity we need to find $\frac{a}{b}$, which is uniquely determined from $r_{eq}$ with following equation:

$$ r_{eq} = \sqrt{\frac{\sigma_{c-air}}{g(\rho_{c,l} - \rho_{air})}}\left(\frac{b}{a}\right)^{-\frac{1}{6}}\sqrt{\left(\frac{b}{a}\right)^{-2} - 2\left(\frac{b}{a}\right)^{-\frac{1}{3}} + 1}$$

This can be done by transforming it into form $F(x) = 0$ where the variable $x = \frac{b}{a}$ and finding roots of $F$.

$$ F\left(\frac{b}{a}\right) = \sqrt{\frac{\sigma_{c-air}}{g(\rho_{c,l} - \rho_{air})}}\left(\frac{b}{a}\right)^{-\frac{1}{6}}\sqrt{\left(\frac{b}{a}\right)^{-2} - 2\left(\frac{b}{a}\right)^{-\frac{1}{3}} + 1} - r_{eq} = 0$$

Where:
- $r_{eq}$ – equivalent radius of the drop  
- $a$, $b$ – semi-major and semi-minor axes of the ellipsoidal drop  
- $\sigma_{c\text{-}air}$ – surface tension between cloud water and air  
- $g$ – gravitational acceleration  
- $\rho_{c,l}$ – density of cloud liquid water  
- $\rho_{air}$ – density of air


## Tests

To validate the predictions of the physical model presented by Loftus & Wordsworth (2021) by verifying whether real-world observations include raindrops exceeding the maximum survivable size as constrained by evaporation and stability. Open-access data from the DISDRODB project or NASA website provide high-temporal-resolution measurements from video disdrometers, used to analyze raindrop size and velocity.

Data includes:
* raindrop diameter

* fall speed

* time of observation

* type of precipitation

It allows comparison between theoretical models and actual raindrop size distributions


## Model Prediction

The model defines a narrow range of viable raindrop radii that are:
* Large enough to avoid complete evaporation during descent,

* Small enough to remain stable and avoid breakup due to aerodynamic forces.



## Visualization

Function: plot_radius_vs_altitude()

- Plots the evolution of radius vs. altitude for a set of initial drop sizes.
- Differentiates between surviving and evaporated drops with colors and markers.

Example:

    from simulation.plotting import plot_radius_vs_altitude
    radii0 = Q_(np.geomspace(0.12, 0.9, 6), "mm")
    z0 = Q_(600.0, "m")
    plot_radius_vs_altitude(radii0, z0, earth, simulate_raindrop_descent)

## Planetary Support


You can simulate droplet behavior in atmospheres of:

- Earth
- Mars
- Saturn
- Jupiter
- Titan
- Exoplanet K2-18b

Each planet is defined with:

- Surface temperature $T_0$
- Relative humidity RH
- Pressure $p_0$
- Gravitational acceleration g
- Lapse rate $Gamma_d$
- Air density $rho_{air}$
- Lifting condensation level temperature $T_{LCL}$

## Example Plot

<p align="center">
  <img src="https://github.com/user-attachments/assets/fd13a570-9a54-4e5e-838f-b9a611eb1c90" alt="Raindrop Simulation Plot" width="600"/><br/>
  <em>Figure: Droplet radius vs. altitude for various initial sizes.</em>
</p>



