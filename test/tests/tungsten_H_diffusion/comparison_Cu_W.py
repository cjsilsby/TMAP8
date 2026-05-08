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

tmap_sol_tmap7 = pd.read_csv("./Cu_W_Diff_csv.csv")

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

# Parameters for TMAP 7 Analytical solution
# l = 66e-6  # thickness of the SiC layer (m)
lambdas = get_lambdas_analytical(k, l, a)
t = np.expand_dims(tmap_time_tmap7, axis=0)
x = 0.02e-3  # depth into W layer from Cu/W interface (D_ver)
# where we compare analytical and numerical model concentration predictions (m)
x2 = x + a

summation = (
    (
        D_Cu * l * sin(lambdas) * sin(k * l / a * lambdas) * (cos(lambdas) - 1)
        + D_W
        * sin(lambdas)
        * (
            k * l * sin(lambdas) * cos(k * l / a * lambdas)
            - a * sin(k * l / a * lambdas)
        )
    )
    / (
        lambdas
        * (a * D_W + l * D_Cu)
        * (np.power(sin(k * l / a * lambdas), 2) + l / a * np.power(sin(lambdas), 2))
    )
    * sin(k * lambdas * (l + a - x2) / a)
    * exp(-D_Cu * np.power(lambdas / a, 2) * t.transpose())
)
sums = np.sum(summation, axis=1)

analytical_conc_tmap7 = c0 * (D_Cu * (l + a - x2) / (l * D_Cu + a * D_W) + 2 * sums)

idx = np.where(tmap_time_tmap7 >= t0)[0]
RMSE = np.sqrt(np.mean((tmap_conc_tmap7[idx] - analytical_conc_tmap7[idx]) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc_tmap7[idx])
ax.text(0.1, 0.07, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_time_tmap7,
    analytical_conc_tmap7,
    label=r"Analytical-W",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)


ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(0, 0.2)
ax.set_ylim(0, 0.10)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_comparison_time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ========= Comparison of concentration as a function of time in Cu side ===================

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(tmap_time_tmap7, tmap_conc_tmap7_Cu, label=r"TMAP8-Cu", c="tab:brown")

# Parameters for TMAP 4 Analytical solution
# l = 63e-6  # thickness of the SiC layer (m)
lambdas = get_lambdas_analytical(k, l, a)
t = np.expand_dims(tmap_time_tmap7, axis=0)
x = -0.005e-3  # depth into Cu layer from Cu / W interface (mm) (D_ver_Cu - T_Cu)
# where we compare analytical and numerical model concentration predictions (mm)
x1 = x + a  # (D_ver_Cu)

summation = (
    (
        D_Cu * l * np.power(sin(k * l / a * lambdas), 2) * (cos(lambdas) - 1)
        + D_W
        * sin(k * l / a * lambdas)
        * (
            k * l * sin(lambdas) * cos(k * l / a * lambdas)
            - a * sin(k * l / a * lambdas)
        )
    )
    / (
        lambdas
        * (a * D_W + l * D_Cu)
        * (np.power(sin(k * l / a * lambdas), 2) + l / a * np.power(sin(lambdas), 2))
    )
    * sin(lambdas * x1 / a)
    * exp(-D_Cu * np.power(lambdas / a, 2) * t.transpose())
)
sums = np.sum(summation, axis=1)

analytical_conc_tmap7 = c0 * (
    (D_Cu * l + (a - x1) * D_W) / (l * D_Cu + a * D_W) + 2 * sums
)

idx = np.where(tmap_time_tmap7 >= t0)[0]
RMSE = np.sqrt(np.mean((tmap_conc_tmap7_Cu[idx] - analytical_conc_tmap7[idx]) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc_tmap7[idx])
ax.text(0.1, 0.07, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_time_tmap7,
    analytical_conc_tmap7,
    label=r"Analytical-Cu",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)

ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(0, 0.2)
ax.set_ylim(0, 0.10)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_comparison_Cu_time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Comparison of concentration as a function of distance ============
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

csv_folder = "./Cu_W_Diff_vector_postproc_line_0083.csv"  # 0056.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_distance_tmap7 = tmap_sol["x"]
tmap_distance_tmap7_microns = tmap_distance_tmap7 * 1e3
tmap_conc_tmap7 = tmap_sol["u"]
ax.plot(
    tmap_distance_tmap7_microns,
    tmap_conc_tmap7,
    label=r"TMAP8",
    c="tab:brown",
)

# TMAP 7 Analytical solution
c0 = 0.1  # concentration at the Cu free surface (moles/mm^3)

x = tmap_distance_tmap7
# Cu_conc = c0 * (1 + (x / l) * ((a * D_Cu) / (a * D_Cu + l * D_W) - 1))
# W_conc = c0 * (((a + l - x) / l) * (a * D_Cu) / (a * D_Cu + l * D_W))
Cu_conc = c0 * (1 - x * D_W / (a * D_W + l * D_Cu))
W_conc = c0 * D_Cu * (a + l - x) / (a * D_W + l * D_Cu)
analytical_conc_tmap7 = (x < a) * Cu_conc + (x >= a) * W_conc
analytical_conc_tmap7 = (x < a) * Cu_conc + (x >= a) * W_conc

RMSE = np.sqrt(np.mean((tmap_conc_tmap7 - analytical_conc_tmap7) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc_tmap7)
ax.text(0.03, 0.02, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_distance_tmap7_microns,
    analytical_conc_tmap7,
    label=r"Analytical",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)

ax.set_xlabel(r"Distance (mm)")
ax.set_ylabel(r"Concentration (moles H/m$^3$)")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.legend(loc="best")
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

ax.minorticks_on()
plt.savefig("Cu_W_dist_comparison.png", bbox_inches="tight", dpi=300)
plt.close(fig)
