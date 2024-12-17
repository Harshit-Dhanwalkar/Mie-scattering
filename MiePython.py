#!/usr/bin/env python
# coding: utf-8

# Python Mie code

# PyMieScatt is a python mie code based on based on Bohren and Huffman's Mie Theory derivations which calculates extinction efficiency (Qext),
# scattering efficiency (Qsca), backscattering efficiency (Qback) and asymmetry parameter (g). It requires the refractive index of the scattering particle and the size parameter.
#
# PyMieScatt is available at
# https://github.com/bsumlin/PyMieScatt
# and can be installed by running the following command in terminal.
#
# _pip install PyMieScatt_ or _conda install -c conda-forge pymiescatt_

import math

import matplotlib
import numpy as np
import PyMieScatt as mie  # Importing miepython module for mie calculation
from scipy.integrate import trapezoid as trapz  # Importing scipy module for integration

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

# Single particle
#
# Mie code calculates the mie efficiencies and asymmetry parameter from the particle refractive index, particle diameter and wavelength of interacting radiation

wavelength = 0.5  # Wavelength of incident radiation
radius = 0.1  # Radius of the scattering particle
ref_index = (
    1.45 + 0.45j
)  # Refractive index of the scattering particle in the format n+ik
mie_params = mie.MieQ(ref_index, wavelength, radius * 2, asDict=True)
print(mie_params)


# Polydisperse particles
#
# In atmosphere particles exist in different size and composition. In order to assess the scattering properties of atmospheric particles, we have to take account of their size distribution. Size distribution of atmospheric particles can be best assumed as lognormal which can be described by its
#
# - mod radius
# - standard deviation
# - upper and lower bounds of particle size
#  <br>${dn(r)/dr} = {(n/(\sqrt{2\pi}\ln{\sigma}r)\exp{[-\ln^{2}(r/r_{m})/2(ln\sigma)^{2}]}} $

n = 100  # Number of particles
rm = 0.0118  # Mod radius
sigma = 2  # Standard deviation
rmin = 0.005  # Minimum radius
rmax = 10  # Maximum radius
rad_array = np.arange(rmin, rmax, 0.001)
dnr_array = []
for r in rad_array:
    dnr = (
        math.sqrt(math.pi / 2.0)
        * (n / np.log(sigma))
        * r
        * math.exp(-0.5 * (((np.log(r / rm)) ** 2) / (np.log(sigma)) ** 2))
    )
    dnr_array.append(dnr)
fig = plt.figure(figsize=(5, 5))
plt.plot(rad_array, dnr_array)
plt.xscale("log")
plt.xlabel("radius (um)")
plt.ylabel("dn(r)/dr")
plt.title("Number-size lognormal distribution")
plt.grid()
plt.show()
# Or plt.savefig("lognormal_distribution.png") to save

wavelength = 0.5  # Wavelength of incident radiation
ref_index = 1.45 + 0.45j  # Refractive index of the scattering particle in
# the format n+ik
n = 1  # Number of particles
rm = 0.0118  # Mod radius
sigma = 2  # Standard deviation
rmin = 0.005  # Minimum radius
rmax = 10  # Maximum radius
rad_array = np.arange(rmin, rmax, 0.001)


def miecoeff(n, sigma, rad_array, lamb, rm, m):
    bext_array = []  # Arrays to store f(x)
    bsca_array = []  # of each radius value
    babs_array = []  # for each radius value
    g_array = []

    for r in rad_array:
        mie_params = mie.MieQ(
            ref_index, wavelength, radius * 2, asDict=True
        )  # Calling mie program to calculate
        # mie parameters
        bext = (
            math.sqrt(math.pi / 2.0)
            * (n / np.log(sigma))
            * r
            * mie_params["Qext"]
            * math.exp(-0.5 * (((np.log(r / rm)) ** 2) / (np.log(sigma)) ** 2))
        )  # log-normal distribution function
        bsca = (
            math.sqrt(math.pi / 2.0)
            * (n / np.log(sigma))
            * r
            * mie_params["Qsca"]
            * math.exp(-0.5 * (((np.log(r / rm)) ** 2) / (np.log(sigma)) ** 2))
        )
        babs = (
            math.sqrt(math.pi / 2.0)
            * (n / np.log(sigma))
            * r
            * mie_params["Qabs"]
            * math.exp(-0.5 * (((np.log(r / rm)) ** 2) / (np.log(sigma)) ** 2))
        )
        g = (
            math.sqrt(math.pi / 2.0)
            * (n / np.log(sigma))
            * r
            * mie_params["Qsca"]
            * mie_params["g"]
            * math.exp(-0.5 * (((np.log(r / rm)) ** 2) / (np.log(sigma)) ** 2))
        )

        bext_array.append(bext)  # Appending values to array
        bsca_array.append(bsca)
        babs_array.append(babs)
        g_array.append(g)
    return bext_array, bsca_array, babs_array, g_array


bext_array, bsca_array, babs_array, g_array = miecoeff(
    n, sigma, rad_array, wavelength, rm, ref_index
)
# Integrating vext,bsca,babs and g using trapezoidal rule
bext = trapz(bext_array, rad_array)
bsca = trapz(bsca_array, rad_array)
babs = trapz(babs_array, rad_array)
g = trapz(g_array, rad_array) / bsca

print(
    "bext = "
    + str(bext)
    + "  bsca = "
    + str(bsca)
    + "  babs = "
    + str(babs)
    + "  g = "
    + str(g)
)

# Visualization of mie parameters
#
# Variation of mie efficiencies with size parameter

wavelength = 0.5
ref_index = 1.45 + 0.05j
rmin = 0.005
rmax = 5
rad_array = np.arange(rmin, rmax, 0.001)
qext_array = []
qsca_array = []
qabs_array = []
size_param_array = []
for r in rad_array:
    mie_params = mie.MieQ(ref_index, wavelength, r * 2, asDict=True)
    size_param = (2 * math.pi * r) / wavelength
    qext_array.append(mie_params["Qext"])
    qabs_array.append(mie_params["Qabs"])
    qsca_array.append(mie_params["Qsca"])
    size_param_array.append(size_param)

fig = plt.Figure(figsize=(5, 5))
plt.grid()
plt.plot(size_param_array, qext_array, label="Qext")
plt.plot(size_param_array, qabs_array, label="Qabs")
plt.plot(size_param_array, qsca_array, label="Qsca")
plt.xlabel("Size parameter")
plt.ylabel("Qext/Qabs/Qsca")
plt.grid()
plt.legend()


# Variation of SSA with refractive index

wavelength = 0.5
ref_index = 1.33 + 0.0001j
rmin = 0.005
rmax = 10
rad_array = np.arange(rmin, rmax, 0.001)
qext_array = []
qsca_array = []
qabs_array = []
size_param_array = []
ssa_array = []
for r in rad_array:
    mie_params = mie.MieQ(ref_index, wavelength, r * 2, asDict=True)
    size_param = (2 * math.pi * r) / wavelength
    qext_array.append(mie_params["Qext"])
    qabs_array.append(mie_params["Qabs"])
    qsca_array.append(mie_params["Qsca"])
    size_param_array.append(size_param)
ssa_array = np.divide(qsca_array, qext_array)

fig = plt.Figure(figsize=(5, 5))
plt.plot(size_param_array, ssa_array, label=str(ref_index))
plt.plot(size_param_array, qext_array, label="Qext")
plt.plot(size_param_array, qabs_array, label="Qabs")
plt.plot(size_param_array, qsca_array, label="Qsca")
plt.xlabel("Size parameter")
plt.ylabel("Qext/Qabs/Qsca")
plt.grid()
plt.legend()
plt.show()
# plt.savefig("mie_efficiencies.png") to save

ref_index = 1.33 + 0.01j
qext_array = []
qsca_array = []
qabs_array = []
size_param_array = []
ssa_array = []
for r in rad_array:
    mie_params = mie.MieQ(ref_index, wavelength, r * 2, asDict=True)
    size_param = (2 * math.pi * r) / wavelength
    qext_array.append(mie_params["Qext"])
    qabs_array.append(mie_params["Qabs"])
    qsca_array.append(mie_params["Qsca"])
    size_param_array.append(size_param)
ssa_array = np.divide(qsca_array, qext_array)
plt.plot(size_param_array, ssa_array, label=str(ref_index))

ref_index = 1.33 + 0.1j
qext_array = []
qsca_array = []
qabs_array = []
size_param_array = []
ssa_array = []
for r in rad_array:
    mie_params = mie.MieQ(ref_index, wavelength, r * 2, asDict=True)
    size_param = (2 * math.pi * r) / wavelength
    qext_array.append(mie_params["Qext"])
    qabs_array.append(mie_params["Qabs"])
    qsca_array.append(mie_params["Qsca"])
    size_param_array.append(size_param)
ssa_array = np.divide(qsca_array, qext_array)
plt.plot(size_param_array, ssa_array, label=str(ref_index))

plt.xlabel("Size parameter")
plt.ylabel("SSA")
plt.legend()


# Useful links
#
# https://scattport.org/index.php/programs-menu/mie-type-codes-menu
# http://www.philiplaven.com/mieplot.htm
