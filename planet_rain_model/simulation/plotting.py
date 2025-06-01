import matplotlib.pyplot as plt
import numpy as np

def plot_radius_vs_altitude(radii0, z0, planet, model_fn):
    plt.figure(figsize=(6, 4))

    for r0 in radii0:
        traj = model_fn(r0, z0, planet)
        traj = np.array(traj)
        if len(traj) == 0:
            continue
        z_vals, r_vals = traj[:, 0], traj[:, 1]
        survived = z_vals[-1] <= 1.0 and r_vals[-1] > 0.01
        color = "purple" if survived else "gray"

        plt.plot(r_vals * 1e3, z_vals, color=color, lw=2.0)
        plt.text(r_vals[0] * 1e3, z_vals[0], "$r_0$", fontsize=8, color=color)
        if not survived:
            plt.scatter(r_vals[-1] * 1e3, z_vals[-1], color=color, marker="x", s=30)

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
