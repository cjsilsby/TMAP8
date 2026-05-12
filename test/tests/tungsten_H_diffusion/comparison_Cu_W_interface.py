import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
from scipy import special
from numpy import sin, cos, tan, sqrt, exp
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)


def get_lambdas_analytical(k, l, a):
    # Calculate lambda values for analytical solution
    lambda_range = np.arange(1e-12, 1e2, 1e-5)
    f = 1 / k * sin(lambda_range) * cos(lambda_range * l / a * k)
    g = cos(lambda_range) * sin(lambda_range * l / a * k)
    idx = np.where(np.diff(np.sign(f + g)))
    lambdas = np.expand_dims(lambda_range[idx][::1], axis=0)
    return lambdas


# ========= Comparison of concentration as a function of time in W side ===================

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

tmap_sol_tmap7 = pd.read_csv("./Cu_W_interface_csv.csv")

# tmap_time_tmap4 = tmap_sol_tmap4["time"]
# tmap_conc_tmap4 = tmap_sol_tmap4["concentration_at_x_SiC"]
tmap_time_tmap7 = tmap_sol_tmap7["time"]
tmap_conc_tmap7 = tmap_sol_tmap7["concentration_at_x_W"]
tmap_conc_tmap7_Cu = tmap_sol_tmap7["concentration_at_x_Cu"]

# ax.plot(tmap_time_tmap4, tmap_conc_tmap4, label=r"TMAP8-SiC (TMAP4 case)", c="tab:gray")
ax.plot(tmap_time_tmap7, tmap_conc_tmap7, label=r"TMAP8-W", c="tab:brown")

# Analytical parameters
t0 = 0.0001
c0 = 0.1  # concentration at the Cu free surface (moles/m^3)
a = 0.05e-3  # thickness of the Cu layer (m)
l = 0.025e-3  # thickness of the W layer (m)
temp = 1500  # K
Do_Cu = 1.74e-6  # diffusivity in Cu (m^2/s)
Do_W = 7.44e-8  # diffusivity in W (m^2/s)

D_Cu = Do_Cu * np.exp(-42000 / 8.31446261815324 / temp)  # diffusivity in Cu (m^2/s)
D_W = Do_W * np.exp(-0.13 / 8.61733e-5 / temp)  # diffusivity in W (m^2/s)
k = sqrt(D_Cu / D_W)


ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(0, 0.2)
ax.set_ylim(0, 0.001)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_inter_comparison_time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ========= Comparison of concentration as a function of time in Cu side ===================

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(tmap_time_tmap7, tmap_conc_tmap7_Cu, label=r"TMAP8-Cu", c="tab:brown")

ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(0, 0.2)
ax.set_ylim(0, 0.10)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_inter_comparison_Cu_time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Comparison of concentration as a function of distance ============
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

csv_folder_Cu = "./Cu_W_interface_vector_postproc_Cu_line_Cu_0083.csv"  # 0056.csv"
csv_folder_W = "./Cu_W_interface_vector_postproc_W_line_W_0083.csv"  # 0056.csv"
tmap_sol_Cu = pd.read_csv(csv_folder_Cu)
tmap_sol_W = pd.read_csv(csv_folder_W)
tmap_distance_Cu = tmap_sol_Cu["x"]
tmap_distance_W = tmap_sol_W["x"]
tmap_distance_Cu_microns = tmap_distance_Cu * 1e3
tmap_distance_W_microns = tmap_distance_W * 1e3
tmap_conc_Cu = tmap_sol_Cu["u_Cu"]
tmap_conc_W = tmap_sol_W["u_W"]
ax.plot(
    tmap_distance_Cu_microns,
    tmap_conc_Cu,
    label=r"TMAP8",
    c="tab:brown",
)

ax.plot(
    tmap_distance_W_microns,
    tmap_conc_W,
    # label=r"TMAP8",
    c="tab:brown",
)

ax.set_xlabel(r"Distance (mm)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend(loc="best")
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_inter_comparison.png", bbox_inches="tight", dpi=300)
plt.close(fig)
