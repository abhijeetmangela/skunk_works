"""
===========================================================
UAV Weight Estimation Script

Outputs:
    • Final MTOW
    • Battery weight
    • Empty weight
    • Payload weight
    • Power requirements
    • Energy usage
    • Convergence plot (saved)
    • Empty weight regression plot (saved)

Author: Abhijeet Mangela (AE25M034)
===========================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


# =========================================================
# ------------------- CONSTANTS ----------------------------
# =========================================================

rho = 1.225
g = 9.81

b = 2
AR = 9.686999647
S = b**2 / AR

payload_mass_drop = 1
payload_mass = 2.1

V_climb = 20
V_cruise = 22
V_des = 20
V_loiter = 22

cruise_dist_1 = 40000
cruise_dist_2 = 40000

height_cruise = 200

gamma_climb = 8
gamma_des = 8

C_d_0 = 0.03
e = 0.8

energy_density = 3600 * 120

Initial_weight = 8.5


# =========================================================
# ---------------- POWER MODELS ---------------------------
# =========================================================

def Power_climb(V_cl, gamma_climb, Weight):
    """Climb power requirement"""
    C_l = Weight*np.cos(np.deg2rad(gamma_climb))/(0.5*rho*V_cl**2*S)
    C_d = C_d_0 + C_l**2/(np.pi*AR*e)
    D_cl = 0.5*rho*V_cl**2*S*C_d
    T_cl = D_cl + Weight*np.sin(np.deg2rad(gamma_climb))
    return T_cl * V_cl


def Power_cruise(V_cr, Weight):
    return Power_climb(V_cr, 0, Weight)


def Power_descent(V_des, gamma_des, Weight):
    return max(Power_climb(V_des, -gamma_des, Weight), 0)


def Power_loiter(V_loi, Weight):
    return Power_cruise(V_loi, Weight)


# =========================================================
# ----------------- TIME CALCULATIONS ---------------------
# =========================================================

t_cr_1 = cruise_dist_1 / V_cruise
t_cr_2 = cruise_dist_2 / V_cruise
t_loiter = 15 * 60

t_cl_1 = height_cruise / (V_climb*np.sin(np.deg2rad(gamma_climb)))
t_des_1 = height_cruise / (V_des*np.sin(np.deg2rad(gamma_des)))


# =========================================================
# ----------------- ENERGY MODEL --------------------------
# =========================================================

def Total_energy(Weight):
    """Total mission energy"""
    energy_climb = Power_climb(V_climb, gamma_climb, Weight) * t_cl_1
    energy_cruise_1 = Power_cruise(V_cruise, Weight) * t_cr_1
    energy_cruise_2 = Power_cruise(V_cruise, Weight-payload_mass_drop*g) * t_cr_2
    energy_descent = Power_descent(V_des, gamma_des, Weight-payload_mass_drop*g) * t_des_1
    energy_loiter = Power_loiter(V_loiter, Weight) * t_loiter

    return energy_climb + energy_cruise_1 + energy_cruise_2 + energy_descent + energy_loiter


# =========================================================
# ----------------- LOAD EXCEL DATA -----------------------
# =========================================================

print("\nLoading dataset...")
Data = pd.read_excel("Drone_dataset.xlsx")

# Same filtering as notebook
Data = Data[Data["Empty Weight (kg)"]/Data["MTOW (kg)"] <= 0.55]

Mtow_array = Data["MTOW (kg)"].to_numpy()
Empty_weight_array = Data["Empty Weight (kg)"].to_numpy()

Empty_weight_ratios = Empty_weight_array / Mtow_array


# =========================================================
# ------------- EMPTY WEIGHT REGRESSION -------------------
# =========================================================

logx = np.log10(Mtow_array)
logy = np.log10(Empty_weight_ratios)

A = np.column_stack((logx, np.ones_like(logx)))
Soln, _, _, _ = np.linalg.lstsq(A, logy, rcond=None)

L = Soln[0]
loga = Soln[1]
a = 10**loga

print("\nRegression results")
print("------------------")
print("a =", a)
print("L =", L)


# =========================================================
# ------------- ITERATIVE WEIGHT SOLVER -------------------
# =========================================================

W_0 = Initial_weight
Weight_array = []

for _ in range(100):
    Weight_array.append(W_0)

    W_b = Total_energy(W_0*g) / energy_density

    W_0 = (payload_mass / (1 - (1.2*W_b/W_0) - a*W_0**L)).item()


final_weight = Weight_array[-1]


# =========================================================
# ------------- COMPONENT BREAKDOWN -----------------------
# =========================================================

Battery_weight = Total_energy(final_weight*g) / energy_density
Empty_weight = final_weight*a * final_weight**L
Payload = payload_mass

print("\n================ FINAL RESULTS =================")
print(f"Final MTOW        : {final_weight:.3f} kg")
print(f"Battery weight    : {Battery_weight:.3f} kg")
print(f"Empty weight      : {Empty_weight:.3f} kg")
print(f"Payload           : {Payload:.3f} kg")


# Power reporting
print("\nPower Requirements (W)")
print("----------------------")
print("Climb   :", Power_climb(V_climb, gamma_climb, final_weight*g))
print("Cruise  :", Power_cruise(V_cruise, final_weight*g))
print("Descent :", Power_descent(V_des, gamma_des, final_weight*g))
print("Loiter  :", Power_loiter(V_loiter, final_weight*g))


# =========================================================
# ----------------- PLOTS --------------------------------
# =========================================================

os.makedirs("figures", exist_ok=True)

# ---------------------------------------------------------
# Convergence Plot
# ---------------------------------------------------------

iters = list(range(1, len(Weight_array) + 1))

plt.figure(figsize=(8, 5))
plt.plot(range(len(Weight_array)), Weight_array, 'o-')
plt.xlabel("Iteration")
plt.ylabel("Weight (kg)")
plt.title("Weight Convergence")
plt.annotate(f'{final_weight:.3f} kg', xy=(iters[-1], final_weight), xytext=(-60, 8), textcoords='offset points',
                 arrowprops=dict(arrowstyle='->', color='tab:blue'), fontsize=10, color='tab:blue')
plt.grid(True)
# plt.show()

plt.savefig("figures/weight_convergence.png", dpi=300)
plt.close()


# ---------------------------------------------------------
# Empty weight regression plot
# ---------------------------------------------------------
plot_weight = np.linspace(3, 20, 1000)
plot_empt_rat = a * plot_weight**L

plt.figure(figsize=(8, 5))
plt.scatter(Mtow_array, Empty_weight_ratios, label="Data")
plt.plot(plot_weight, plot_empt_rat, label="Fit")
plt.xlabel("MTOW (kg)")
plt.ylabel("Empty Weight Ratio")
plt.title("Empty Weight vs MTOW")
plt.grid(True)
plt.legend()
# plt.show()

plt.savefig("figures/empty_weight_fit.png", dpi=300)
plt.close()


print("\nPlots saved inside: ./figures/")
print("Done!")
