import math

import matplotlib
import numpy as np
import PyMieScatt as mie  # Importing miepython module for mie calculation
from scipy.integrate import trapezoid as trapz  # Importing scipy module for integration

matplotlib.use("Qt5Agg")  # Requires PyQt5 or PySide2
import matplotlib.pyplot as plt

m = 1.536  # Refractive index of NaCl
wavelength = 405  # Replace with the laser wavelength (nm)

dp_g = 85  # Geometric mean diameter (nm)
sigma_g = 1.5  # Geometric standard deviation (unitless)
N = 1e5  # Total number of particles (cm^-3)

# Calculate optical properties
B = mie.Mie_Lognormal(
    m, wavelength, sigma_g, dp_g, N, numberOfBins=1000, returnDistribution=True
)
S = mie.SF_SD(m, wavelength, B[7], B[8])

# Plot the scattering intensity
plt.close("all")

fig1 = plt.figure(figsize=(10.0, 6))
ax1 = fig1.add_subplot(1, 1, 1)

# Plot Polarizations
ax1.plot(S[0], S[1], "b", ls="dashdot", lw=1, label="Parallel Polarization")
ax1.plot(S[0], S[2], "r", ls="dashed", lw=1, label="Perpendicular Polarization")
ax1.plot(S[0], S[3], "k", lw=1, label="Unpolarized")

# X-axis labels and ticks
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
ax1.set_xlabel("Scattering Angle ϴ", fontsize=16)
ax1.set_ylabel(r"Intensity ($\mathregular{|S|^2}$)", fontsize=16, labelpad=10)

handles, labels = ax1.get_legend_handles_labels()
fig1.legend(handles, labels, fontsize=14, ncol=3, loc=8)
fig1.suptitle("Scattering Intensity Functions", fontsize=18)

# Highlight specific angles
sixty = [0.96 < x < 1.13 for x in S[0]]
ninety = [1.48 < x < 1.67 for x in S[0]]
ax1.fill_between(S[0], 1e2, S[3], where=sixty, color="g", alpha=0.15)
ax1.fill_between(S[0], 1e2, S[3], where=ninety, color="g", alpha=0.15)

# Use log scale
ax1.set_yscale("log")

# Compute integrals for specified ranges
int_sixty = trapz(S[3][110:130], S[0][110:130])
int_ninety = trapz(S[3][169:191], S[0][169:191])

# Annotate with integral results
ax1.annotate(
    f"Integrated value = {int_sixty:.3f}",
    xy=(np.pi / 3, S[3][120]),
    xycoords="data",
    xytext=(np.pi / 3, 12000),
    textcoords="data",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
)

ax1.annotate(
    f"Integrated value = {int_ninety:.3f}",
    xy=(np.pi / 2, S[3][180]),
    xycoords="data",
    xytext=(np.pi / 2, 8000),
    textcoords="data",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
)

plt.tight_layout(rect=[0.01, 0.05, 0.915, 0.95])
plt.show()
