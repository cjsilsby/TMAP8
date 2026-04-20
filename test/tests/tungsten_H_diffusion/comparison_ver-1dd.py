import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

num_summation_terms = 2000
N_o = 3.1622e13  # (atom/mm^3) For convenience lattice density chosen as 3.1622e22 atom/m^3
c_o = 0.0001  # dissolved gas atom fraction (-)
D = 1.632  # diffusivity (mm^2/min) ---- 1.632 m^2/min, but input as if s. Corrected before graphing.
l = 0.025  # slab thickness (mm) --- 0.025 mm

# Extract data from 'gold' TMAP8 run
# if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
#     csv_folder = "../../../../test/tests/ver-1dd/gold/ver-1dd_out.csv"
# else:  # if in test folder
#     csv_folder = "./gold/ver-1dd_out.csv"
# csv_folder = "./ver-1dd_out.csv"
csv_folder = "./ver-1dd_out.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_time = np.array(tmap_sol["time"]) # minutes
tmap_prediction = np.array(tmap_sol["scaled_outflux"])
In_minus_Out = np.array(tmap_sol["time_integrated_flux"])
Accumulation = np.array(tmap_sol["mass_in_domain"])


idx = np.where(tmap_time >= 0.0000001)[0][0]

# Calculate the breakthrough time from numerical solution
tmap_slope = (tmap_prediction[idx + 1 :] - tmap_prediction[idx:-1]) / (
    tmap_time[idx + 1 :] - tmap_time[idx:-1]
)
tmap_intercept = tmap_time[np.argmax(tmap_slope) + idx] - (
    tmap_prediction[int(np.argmax(tmap_slope) + idx)]
) / np.max(tmap_slope)
output_line = f"The breakthrough time from numerical solution is {tmap_intercept:.2e} s"
print(output_line)

# analytical solution
analytical_time = np.array(tmap_time)
tau_be = l**2 / (2 * (np.pi) ** 2 * D)  # calculate analytical breakthrough time
output_line = f"The breakthrough time from analytical solution is {tau_be:.2e} s"
print(output_line)


def summation_term(num_terms, time):
    sum = 0.0
    for m in range(1, num_terms):
        sum += (-1) ** m * np.exp(-1 * m**2 * time / (2 * tau_be))
    return sum


# Calculate the analytical solution
analytical_flux = (
    N_o * (c_o * D / l) * (1 + 2 * summation_term(num_summation_terms, analytical_time))
)

# Plot figure for verification
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

# Correct time to minutes
tmap_time = tmap_time * 60
analytical_time = analytical_time * 60

ax.plot(tmap_time, tmap_prediction, label=r"TMP8", c="tab:gray")  # numerical solution
ax.plot(
    analytical_time, analytical_flux, label=r"Analytical", c="k", linestyle="--"
)  # analytical solution
# ax.plot(
#     [tmap_intercept, 0.221],
#     [0, (0.221 - tmap_intercept) * np.max(tmap_slope)],
#     label=r"Numerical breakthrough time",
#     c="tab:brown",
# )
ax.plot(
    [tau_be, tau_be],
    [0, 2.0e11],
    label=r"Analytical breakthrough time",
    c="tab:brown",
    linestyle="--",
)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Flux (atom/mm$^2$s)") # For tungsten in atom/mm^2, for TZM atom/mm^2
ax.legend(loc="best")
# ax.set_xlim(left=0, right=1.5)
ax.set_ylim(bottom=0)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
RMSE = np.sqrt(np.mean((tmap_prediction - analytical_flux)[idx:] ** 2))
RMSPE = RMSE * 100 / np.mean(analytical_flux[idx:])
ax.text(0.015, 1.5e11, "RMSPE = %.2f " % RMSPE + "%", fontweight="bold")
# ax.text(
#     0.25e-5,
#     0.5e14,
#     "Numerical breakthrough time = %.2e " % tmap_intercept + "s",
#     fontweight="bold",
# )
# ax.text(
#     0.25e-5,
#     0.18e14,
#     "Analytical breakthrough time = %.2e " % tau_be + "s",
#     fontweight="bold",
# )
ax.minorticks_on()
plt.savefig("ver-1dd_comparison_diffusion.png", bbox_inches="tight", dpi=300)
plt.close(fig)


################# Mass Balance Check ##########################

plt.figure(figsize=(6.5,5.5))
plt.plot(tmap_time, Accumulation, c="k", linestyle="--", label = 'Accumulation')
plt.plot(tmap_time, In_minus_Out, label = 'In - Out')
RMSE = np.sqrt(np.mean((Accumulation-In_minus_Out)**2) )
RMSPE = RMSE*100/np.mean(In_minus_Out)
# print(f'RMSPE = %.2f '%RMSPE+'%')
plt.text(0.01,0.006, 'RMSPE = %.2f '%RMSPE+'%',fontweight='bold')
plt.xlabel('Time (s)')
plt.ylabel(r'$\mu$mol H/mm')
plt.title(f'In - Out = Accumulation')
plt.xlim(0,tmap_time.max())
plt.ylim(0)
plt.legend()
plt.grid(True)
plt.tight_layout()
# plt.show()

plt.savefig("ver-1dd_MassBal.png", bbox_inches="tight", dpi=300)
plt.close(fig)