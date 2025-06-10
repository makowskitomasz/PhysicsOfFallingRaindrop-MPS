import matplotlib.pyplot as plt
import numpy as np

def plot_radius_vs_altitude(radii0, z0, planet, model_fn, figsize=(5, 5)):
    """
    Plots the radius of raindrops as a function of altitude during their descent.

    Parameters:
        radii0 (list of Quantity): Initial radii of raindrops (in meters).
        z0 (Quantity): Initial altitude of raindrops (in meters).
        planet (dict): Dictionary containing planetary parameters used in the model.
        model_fn (function): Function to simulate the raindrop descent.
        figsize (tuple): Size of the figure (default is (5, 5)).

    Returns:
        None: Displays the plot.
    """
    plt.style.use("dark_background")
    plt.figure(figsize=figsize)

    for r0 in radii0:
        # Simulate the trajectory of the raindrop
        traj = model_fn(r0, z0, planet)
        traj = np.array(traj)
        if len(traj) == 0:
            continue
        
        # Extract altitude and radius values
        z_vals, r_vals = traj[:, 0], traj[:, 1]
        
        # Determine if the raindrop survives the descent
        survived = z_vals[-1] <= 1.0 and r_vals[-1] > 0.01
        color = "purple" if survived else "gray"

        # Plot the trajectory
        plt.plot(r_vals * 1e3, z_vals, color=color, lw=2.0)
        plt.text(r_vals[0] * 1e3, z_vals[0], "$r_0$", fontsize=8, color=color)
        if not survived:
            plt.scatter(r_vals[-1] * 1e3, z_vals[-1], color=color, marker="x", s=30)

    # Configure plot settings
    plt.xscale("log")
    plt.xlim(1e-2, 1e0)
    plt.ylim(0, 600)
    plt.xlabel(r"$r_{\mathrm{eq}}(z)$ [mm]")
    plt.ylabel(r"$z$ [m]")
    plt.xticks([1e-2, 1e-1, 1e0], ["$10^{-2}$", "$10^{-1}$", "$10^0$"])
    plt.yticks([100, 200, 300, 400, 500, 600])
    plt.grid(False)
    plt.tight_layout()
    plt.show()
