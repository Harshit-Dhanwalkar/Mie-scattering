import math

import matplotlib
import numpy as np
import PyMieScatt as mie  # Importing miepython module for mie calculation
from scipy.integrate import trapezoid as trapz  # Importing scipy module for integration

matplotlib.use("Qt5Agg")  # Requires PyQt5 or PySide2
import matplotlib.pyplot as plt

m = 1.7 + 0.5j
w = 532
d = 5000

theta, SL, SR, SU = mie.ScatteringFunction(m, w, d)
qR, SLQ, SRQ, SUQ = mie.ScatteringFunction(m, w, d, space="qspace")

fig1 = plt.figure(figsize=(10, 6))
ax1 = fig1.add_subplot(1, 2, 1)
ax2 = fig1.add_subplot(1, 2, 2)

# Subplot 1 - Scattering Intensity
ax1.semilogy(theta, SL, "b", ls="dashdot", lw=1, label="Parallel Polarization")
ax1.semilogy(theta, SR, "r", ls="dashed", lw=1, label="Perpendicular Polarization")
ax1.semilogy(theta, SU, "k", lw=1, label="Unpolarized")

x_label = [
    "0",
    r"$\mathregular{\frac{\pi}{4}}$",
    r"$\mathregular{\frac{\pi}{2}}$",
    r"$\mathregular{\frac{3\pi}{4}}$",
    r"$\mathregular{\pi}$",
]
x_tick = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
ax1.set_xticks(x_tick)
ax1.set_xticklabels(x_label, fontsize=14)
ax1.tick_params(which="both", direction="in")
ax1.set_xlabel(r"$\theta$", fontsize=16)
ax1.set_ylabel(r"Intensity ($\mathregular{|S|^2}$)", fontsize=16, labelpad=10)

# Subplot 2 - qR Space
ax2.loglog(qR, SLQ, "b", ls="dashdot", lw=1, label="Parallel Polarization")
ax2.loglog(qR, SRQ, "r", ls="dashed", lw=1, label="Perpendicular Polarization")
ax2.loglog(qR, SUQ, "k", lw=1, label="Unpolarized")

ax2.tick_params(which="both", direction="in")
ax2.set_xlabel("qR", fontsize=14)

# Legend and Titles
handles, labels = ax1.get_legend_handles_labels()
fig1.legend(handles, labels, fontsize=14, ncol=3, loc="lower center")
fig1.suptitle("Scattering Intensity Functions", fontsize=18)

# Adjust layout and display
plt.tight_layout(rect=[0.01, 0.05, 0.915, 0.95])
plt.show()
