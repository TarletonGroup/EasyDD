# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:55:48 2017

@author: daniel
"""
# %%
import matplotlib.pyplot as plt
import numpy as np
# Close all figures.
plt.close("all")
# Matplotlib parameters to use TeX font and customise text size.
plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=18)
plt.rcParams["lines.linewidth"] = 4
plt.rcParams["text.latex.preamble"] = r"\usepackage{bm}"
plt.rcParams.update({'figure.autolayout': True})

data = np.loadtxt(
    "D:\DPhil\OneDrive - Nexus365\EasyDD\src\wip\DanielCelisGarza\Accurate_evaluation_of_dislocation_tractions\paper\images\\test.txt")
theta = data[0:, 0]
fx = data[0:, 1]
fy = data[0:, 2]
fz = data[0:, 3]


fig = plt.figure(1, figsize=(10, 10/1.618))
ax = fig.add_subplot(111)
ax.plot(theta, fx, "-", label=r"$F_{x}(\theta)$")
ax.plot(theta, fy, "-", label=r"$F_{y}(\theta)$")
ax.plot(theta, fz, "-", label=r"$F_{z}(\theta)$")
ax.legend(bbox_to_anchor=(0.3, 0.4))
ax.set_ylabel(r"Force, A.U.")
ax.set_xlabel(r"$\theta$, deg")
ax.set_xlim([-180, 180])
ax.set_yticks(np.arange(-0.05, 0.06, 0.01))
ax.set_xticks(np.arange(-180, 210, 30))
ax.grid(True, which='both')
plt.show()
plt.savefig('ftot_rotation_lin_rect.pdf', format='pdf')

fig = plt.figure(2, figsize=(10, 10/1.618))
ax = fig.add_subplot(111)
ax.plot(theta, fx, "-", label=r"$F_{x}(\theta)$")
ax.plot(theta, fy, "-", label=r"$F_{y}(\theta)$")
ax.plot(theta, fz, "-", label=r"$F_{z}(\theta)$")
ax.legend(bbox_to_anchor=(0.3, 0.4))
ax.set_ylabel(r"Force, A.U.")
ax.set_xlabel(r"$\theta$, deg")
ax.set_xlim([-10, 10])
ax.set_yticks(np.arange(-0.05, 0.06, 0.01))
ax.set_xticks(np.arange(-10, 11, 2))
ax.grid(True, which='both')
plt.show()
plt.savefig('ftot_rotation_lin_rect_zoom.pdf', format='pdf')

# %%
