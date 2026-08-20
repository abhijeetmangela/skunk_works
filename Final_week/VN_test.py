import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# UAV V-n DIAGRAM
# ============================================================

# -----------------------------
# Aircraft parameters
# -----------------------------

W = 6.592 * 9.81       # Weight [N]
S = 0.5307             # Wing area [m^2]
rho = 1.220            # Air density [kg/m^3]

CL_max_pos = 1.3       # Maximum positive lift coefficient
CL_max_neg = -0.8      # Maximum negative lift coefficient

n_max_pos = 3.5        # Maximum positive load factor
n_max_neg = -1.5       # Maximum negative load factor

V_max = 24.0           # Maximum operating speed [m/s]
V_d = 1.25 * V_max     # Design dive speed [m/s]


# ============================================================
# CHARACTERISTIC VELOCITIES
# ============================================================

# Positive 1-g stall speed
V_s = np.sqrt(
    (2 * W) / (rho * S * CL_max_pos)
)

# Positive maneuvering speed
V_a = V_s * np.sqrt(n_max_pos)

# Negative-g stall speed
V_g = np.sqrt(
    (2 * abs(n_max_neg) * W) /
    (rho * S * abs(CL_max_neg))
)


# ============================================================
# VELOCITY ARRAY
# ============================================================

V = np.arange(0.0, V_d + 0.1, 0.1)


# ============================================================
# V-n CURVES
# ============================================================

# Positive stall boundary
n_curve_pos = (
    rho * V**2 * S * CL_max_pos
) / (2 * W)

# Negative stall boundary
n_curve_neg = (
    rho * V**2 * S * CL_max_neg
) / (2 * W)


# Only retain positive curve up to Va
n_curve_pos[V > V_a] = np.nan

# Only retain negative curve up to Vg
n_curve_neg[V > V_g] = np.nan


# ============================================================
# VELOCITY DATA FOR THE V-n DIAGRAM
# ============================================================

# Positive stall curve velocities
V_positive = V[~np.isnan(n_curve_pos)]
n_positive = n_curve_pos[~np.isnan(n_curve_pos)]

# Negative stall curve velocities
V_negative = V[~np.isnan(n_curve_neg)]
n_negative = n_curve_neg[~np.isnan(n_curve_neg)]


# ============================================================
# PRINT IMPORTANT VELOCITIES
# ============================================================

print("\n==============================")
print("      UAV V-n DIAGRAM")
print("==============================")

print(f"Vs   = {V_s:.2f} m/s")
print(f"Va   = {V_a:.2f} m/s")
print(f"Vg   = {V_g:.2f} m/s")
print(f"Vmax = {V_max:.2f} m/s")
print(f"Vd   = {V_d:.2f} m/s")

print("==============================\n")


# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(10, 7))

# Positive stall curve
plt.plot(
    V_positive,
    n_positive,
    linewidth=2,
    label="Positive stall boundary"
)

# Negative stall curve
plt.plot(
    V_negative,
    n_negative,
    linewidth=2,
    label="Negative stall boundary"
)


# Maximum positive load factor
plt.plot(
    [V_a, V_d],
    [n_max_pos, n_max_pos],
    linestyle="--",
    linewidth=2
)

# Maximum negative load factor
plt.plot(
    [V_g, V_d],
    [n_max_neg, n_max_neg],
    linestyle="--",
    linewidth=2
)


# Vmax line
plt.plot(
    [V_max, V_max],
    [n_max_neg, n_max_pos],
    linestyle="--",
    linewidth=1.5
)


# Vd line
plt.plot(
    [V_d, V_d],
    [n_max_neg, n_max_pos],
    linewidth=3
)


# Mark Va
plt.scatter(
    V_a,
    n_max_pos,
    marker="o",
    s=60,
    zorder=5
)

# Mark Vg
plt.scatter(
    V_g,
    n_max_neg,
    marker="o",
    s=60,
    zorder=5
)

# Mark Vd
plt.scatter(
    V_d,
    0,
    marker="o",
    s=60,
    zorder=5
)


# ============================================================
# LABELS
# ============================================================

plt.text(
    V_s,
    0.15,
    f"$V_s$ = {V_s:.2f} m/s"
)

plt.text(
    V_a,
    n_max_pos + 0.15,
    f"$V_a$ = {V_a:.2f} m/s"
)

plt.text(
    V_g,
    n_max_neg + 0.15,
    f"$V_g$ = {V_g:.2f} m/s"
)

plt.text(
    V_max,
    0.2,
    f"$V_{{max}}$ = {V_max:.2f} m/s"
)

plt.text(
    V_d,
    0.2,
    f"$V_d$ = {V_d:.2f} m/s"
)


# ============================================================
# AXES / GRID
# ============================================================

plt.xlabel("Velocity, V [m/s]", fontweight="bold")
plt.ylabel("Load Factor, n", fontweight="bold")

plt.title("UAV V-n Diagram", fontweight="bold")

plt.xlim(0, V_d + 5)
plt.ylim(n_max_neg - 0.5, n_max_pos + 0.5)

plt.axhline(0, linewidth=0.8)

plt.grid(True, alpha=0.3)

plt.legend()

plt.tight_layout()

plt.show()
