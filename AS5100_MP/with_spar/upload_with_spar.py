# Code for wing with spar
# Author: Abhijeet, Mahesh

# Shear flow code : Abhijeet 
# Buckling analysis: Mahesh


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('font', size=14)
np.set_printoptions(precision=5, suppress=True)

# ---------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------

airfoil_file = "s7055.dat"

# Wing / section data
chord = 0.265
skin_thickness = 0.5e-3

# Stringer geometry
stringer_thickness = 0.5e-3
stringer_length = 20e-3

# Number of stringers to add.
# Change this to match your section.
num_stringers = 2

# Stringer positions as normalized chordwise coordinates (x/c)
#
# Enter upper and lower stringer locations separately.
#
# Example:
# upper_stringer_x_norm = np.array([0.10, 0.30, 0.70, 0.90])
# lower_stringer_x_norm = np.array([0.10, 0.30, 0.70, 0.90])
#
# x/c = 0 -> Leading Edge
# x/c = 1 -> Trailing Edge

upper_stringer_x_norm = np.array([0.65])
lower_stringer_x_norm = np.array([0.65])

# Vertical shear force used only to obtain the shear-center moment.
# The resulting shear-center location is independent of its magnitude.
V = 28.1  # N

spar_location_norm = 0.3174
spar_location = spar_location_norm*chord
theo_spar_height = 0.028
flange_length = 0.03
flange_thickness = 1e-3
web_thickness = 1e-3
Xref = 0.25
Cmref = -0.073
CL = 0.72
Xcp = Xref - Cmref/CL
print(f"Chordwise Position of Center of Pressure (from LE) = {Xcp:.2g}")
# LOAD AIRFOIL
# ---------------------------------------------------------------------

airfoil = np.loadtxt(airfoil_file) * chord

# Selig coordinates contain a duplicate trailing-edge point.
# Remove the final duplicate if it is coincident with the first point.
if np.allclose(airfoil[0], airfoil[-1]):
    airfoil = airfoil[:-1]

x = airfoil[:, 0]
y = airfoil[:, 1]

n = len(x)

# Closed panels
x_next = np.roll(x, -1)
y_next = np.roll(y, -1)

dx = x_next - x
dy = y_next - y
ds = np.hypot(dx, dy)

# Panel midpoints
x_mid = 0.5 * (x + x_next)
y_mid = 0.5 * (y + y_next)

n_points = len(airfoil)
print(n_points)
# ---------------------------------------------------------------------
# 1. SKIN AREA IDEALISATION
# ---------------------------------------------------------------------

A_skin_panel = skin_thickness * ds
A_skin = np.sum(A_skin_panel)

x_centroid_skin = np.sum(A_skin_panel * x_mid) / A_skin
y_centroid_skin = np.sum(A_skin_panel * y_mid) / A_skin

print(f"Skin area       = {A_skin:.6e} m^2")
print(f"Skin centroid x = {x_centroid_skin:.6e} m")
print(f"Skin centroid y = {y_centroid_skin:.6e} m")

yc_skin = y - y_centroid_skin

# Check for points too close to the reference axis.
if np.any(np.abs(yc_skin) < 1e-10):
    raise ValueError(
        "A boom lies essentially on the skin centroidal axis. "
        "The standard boom-area formula becomes singular. "
        "Use a different boom/reference discretisation."
    )

l_prev = np.roll(ds, 1)
l_next = ds

y_prev = np.roll(yc_skin, 1)
y_next = np.roll(yc_skin, -1)

B_skin = (
    skin_thickness / 6.0
    * (
        l_prev * (2.0 + y_prev / yc_skin)
        + l_next * (2.0 + y_next / yc_skin)
    )
)


# ---------------------------------------------------------------------
# 3. ADD STRINGER BOOM AREAS
# ---------------------------------------------------------------------

B_stringer = np.zeros(n)

# Stringer boom area
A_stringer = stringer_thickness * stringer_length


def find_surface_index(x_norm, surface):
    """
    Find the S7055 coordinate closest to the requested x/c
    on either the upper or lower surface.
    """

    x_target = x_norm * chord

    if surface == "upper":
        candidates = np.where(y >= 0)[0]
    elif surface == "lower":
        candidates = np.where(y <= 0)[0]
    else:
        raise ValueError("surface must be 'upper' or 'lower'")

    if len(candidates) == 0:
        raise ValueError(f"No points found on {surface} surface.")

    idx = candidates[np.argmin(np.abs(x[candidates] - x_target))]

    return idx


# -------------------------
# Upper stringers
# -------------------------

upper_stringer_indices = []

for x_norm in upper_stringer_x_norm:

    if not 0 <= x_norm <= 1:
        raise ValueError(
            f"Upper stringer x/c = {x_norm} is outside 0 <= x/c <= 1."
        )

    idx = find_surface_index(x_norm, "upper")

    B_stringer[idx] += A_stringer

    upper_stringer_indices.append(idx)


# -------------------------
# Lower stringers
# -------------------------

lower_stringer_indices = []

for x_norm in lower_stringer_x_norm:

    if not 0 <= x_norm <= 1:
        raise ValueError(
            f"Lower stringer x/c = {x_norm} is outside 0 <= x/c <= 1."
        )

    idx = find_surface_index(x_norm, "lower")

    B_stringer[idx] += A_stringer

    lower_stringer_indices.append(idx)


# Total boom area
B_total = B_skin + B_stringer

upper_spar_index = find_surface_index(spar_location_norm, "upper")
lower_spar_index = find_surface_index(spar_location_norm, "lower")

upper_flange_location = spar_location + flange_length
lower_flange_location = spar_location + flange_length

upper_flange_index = find_surface_index(upper_flange_location / chord, "upper")
lower_flange_index = find_surface_index(lower_flange_location / chord, "lower")

print(f"Upper spar index: {upper_spar_index}")
print(f"Upper flange index: {upper_flange_index}")
print(f"Lower spar index: {lower_spar_index}")
print(f"Lower flange index: {lower_flange_index}")
spar_height = y[upper_spar_index] - y[lower_spar_index]
print(f"Spar height: {spar_height}" , 
      f"theoretical spar height: {theo_spar_height}" )
# Adding Bigger web booms

Boom_upper_spar = web_thickness*spar_height*(2+y[lower_spar_index]/y[upper_spar_index])/6 \
    + flange_thickness*ds[upper_spar_index-1]*(2+y[upper_spar_index-1]/y[upper_spar_index])/6

Boom_lower_spar = web_thickness*spar_height*(2+y[upper_spar_index]/y[lower_spar_index])/6 \
    + flange_thickness*ds[lower_spar_index-1]*(2+y[lower_spar_index-1]/y[lower_spar_index])/6

B_total[upper_spar_index] += Boom_upper_spar
B_total[lower_spar_index] += Boom_lower_spar
# adding upper flange booms

for i in range(upper_flange_index, upper_spar_index):
    if i == upper_flange_index:
        B_total[i] += flange_thickness*ds[i]*(2+y[i+1]/y[i])/6 
    else:
        B_total[i] += flange_thickness*ds[i-1]*(2+y[i-1]/y[i])/6 + flange_thickness*ds[i]*(2+y[i+1]/y[i])/6

# adding lower flange booms

for i in range(lower_spar_index + 1, lower_flange_index + 1):
    if i == lower_flange_index:
        B_total[i] += flange_thickness*ds[i-1]*(2+y[i-1]/y[i])/6 
    else:
        B_total[i] += flange_thickness*ds[i-1]*(2+y[i-1]/y[i])/6 + flange_thickness*ds[i]*(2+y[i+1]/y[i])/6
x_centroid_pt = np.sum(B_total * x_mid) / np.sum(B_total)
y_centroid_pt = np.sum(B_total * y_mid) / np.sum(B_total)

x_centroid = x - x_centroid_pt
y_centroid = y - y_centroid_pt

I_xx_centroid = np.sum(B_total*y_centroid**2)
I_yy_centroid = np.sum(B_total*x_centroid**2)
I_xy_centroid = np.sum(B_total*x_centroid*y_centroid)
D = I_xx_centroid*I_yy_centroid - I_xy_centroid**2
# we will make cut on the spar and the last panel

q_b_skin = np.zeros(n_points)

te_cut_panel = n_points-1
# ================================================================
# BASIC SHEAR FLOW WITH CUT AT TRAILING-EDGE PANEL
# ================================================================

q_b_panel = np.zeros(n)

# Panel 79 = trailing-edge cut
# Therefore reference basic shear flow at the cut is zero.
q_b_panel[te_cut_panel] = 0.0

q_running = 0.0

# Panels 0 -> 78
for i in range(n - 1):

    j = i + 1

    # Basic shear flow on panel i
    q_b_panel[i] = q_running

    # Cross the boom at point j
    dq = (
        V / D
        * B_total[j]
        * (
            I_xy_centroid * x_centroid[j]
            - I_xx_centroid * y_centroid[j]
        )
    )

    q_running += dq

# ================================================================
# TWO-CELL IDEALISATION
# ================================================================

# ------------------------------------------------
# 1. CELL PANEL INDICES
# ------------------------------------------------

# Front cell:
# upper spar -> leading edge -> lower spar
front_panels = np.arange(
    upper_spar_index,
    lower_spar_index
)

# Rear cell:
# upper spar -> upper TE
# lower TE -> lower spar
#
# The last panel (n-1 -> 0) is the cut panel.
rear_panels_upper = np.arange(
    0,
    upper_spar_index
)

rear_panels_lower = np.arange(
    lower_spar_index,
    n - 1
)

rear_panels = np.concatenate([
    rear_panels_upper,
    rear_panels_lower
])

# Last panel is the TE cut
te_cut_panel = n - 1


print("\nFront cell panels:")
print(front_panels)

print("\nRear cell panels:")
print(rear_panels)

print("\nTrailing-edge cut panel:")
print(te_cut_panel)
# ================================================================
# 2. CELL AREAS
# ================================================================

def polygon_area(xp, yp):
    """
    Shoelace formula.
    """
    return 0.5 * abs(
        np.sum(xp * np.roll(yp, -1))
        - np.sum(yp * np.roll(xp, -1))
    )


# -------------------------
# Front cell
# -------------------------

front_points = np.arange(
    upper_spar_index,
    lower_spar_index + 1
)

x_front = x[front_points]
y_front = y[front_points]

A1 = polygon_area(x_front, y_front)

# -------------------------
# Rear cell
# -------------------------

rear_points = np.concatenate([
    np.arange(lower_spar_index, n),
    np.arange(0, upper_spar_index + 1)
])

x_rear = x[rear_points]
y_rear = y[rear_points]

A2 = polygon_area(x_rear, y_rear)


print("\nCell areas")
print(f"A1 (front cell) = {A1:.8e} m^2")
print(f"A2 (rear cell)  = {A2:.8e} m^2")
print(f"A1 + A2         = {A1 + A2:.8e} m^2")
# ================================================================
# 3. BASIC TWIST INTEGRALS
# ================================================================

# ------------------------------------------------
# Front cell
# ------------------------------------------------

Iqb_1 = np.sum(
    q_b_panel[front_panels]
    * ds[front_panels]
    / skin_thickness
)


# ------------------------------------------------
# Rear cell
# ------------------------------------------------

# Important:
# The rear-cell positive circulation direction is opposite
# to the global airfoil point ordering.
#
# Therefore q_b for the rear cell has the opposite sign.

Iqb_2 = np.sum(
    (-q_b_panel[rear_panels])
    * ds[rear_panels]
    / skin_thickness
)


print("\nBasic twist integrals")
print(f"Iqb_1 = {Iqb_1:.8e}")
print(f"Iqb_2 = {Iqb_2:.8e}")
# ================================================================
# 4. CELL TWIST STIFFNESS TERMS
# ================================================================

R1 = np.sum(
    ds[front_panels] / skin_thickness
)

R2 = (
    np.sum(ds[rear_panels] / skin_thickness)
    + ds[te_cut_panel] / skin_thickness
)

Rw = spar_height / web_thickness


print("\nTwist stiffness terms")
print(f"R1 = {R1:.8e}")
print(f"R2 = {R2:.8e}")
print(f"Rw = {Rw:.8e}")
# ================================================================
# 5. TORQUE PRODUCED BY BASIC SHEAR FLOW
# ================================================================

T_b = 0.0

for i in range(n - 1):

    # TE cut panel is excluded
    if i == te_cut_panel:
        continue

    j = i + 1

    dx_panel = x[j] - x[i]
    dy_panel = y[j] - y[i]

    xmid_rel = 0.5 * (x_centroid[i] + x_centroid[j])
    ymid_rel = 0.5 * (y_centroid[i] + y_centroid[j])

    dT = q_b_panel[i] * (
        xmid_rel * dy_panel
        - ymid_rel * dx_panel
    )

    T_b += dT

print("\nBasic shear-flow torque")
print(f"T_b = {T_b:.8e} N m")
# ================================================================
# 6. SOLVE FOR q_0_1 AND q_0_2
# ================================================================

# ------------------------------------------------
# Torque equation
#
# 2*A1*q01 + 2*A2*q02 = -T_b
# ------------------------------------------------

a11 = 2.0 * A1
a12 = 2.0 * A2

b1 = -T_b


# ------------------------------------------------
# Compatibility equation
#
# theta_1 = theta_2
# ------------------------------------------------

a21 = (
    (R1 + Rw) / (2.0 * A1)
    + Rw / (2.0 * A2)
)

a22 = (
    -Rw / (2.0 * A1)
    - (R2 + Rw) / (2.0 * A2)
)

b2 = (
    Iqb_2 / (2.0 * A2)
    - Iqb_1 / (2.0 * A1)
)


# ------------------------------------------------
# Matrix
# ------------------------------------------------

M = np.array([
    [a11, a12],
    [a21, a22]
])

rhs = np.array([
    b1,
    b2
])


q_0_1, q_0_2 = np.linalg.solve(M, rhs)


print("\n========================================")
print("REDUNDANT CELL SHEAR FLOWS")
print("========================================")
print(f"q_0_1 (front cell) = {q_0_1:.8e} N/m")
print(f"q_0_2 (rear cell)  = {q_0_2:.8e} N/m")
# ================================================================
# 7. WEB SHEAR FLOW
# ================================================================

q_web = q_0_1 - q_0_2

print("\n========================================")
print("WEB SHEAR FLOW")
print("========================================")
print(f"q_web = {q_web:.8e} N/m")
# ================================================================
# 8. CHECK TWIST COMPATIBILITY
# ================================================================

q_web = q_0_1 - q_0_2

theta_1 = (
    Iqb_1
    + q_0_1 * R1
    + q_web * Rw
) / (2.0 * A1)


theta_2 = (
    Iqb_2
    + q_0_2 * R2
    - q_web * Rw
) / (2.0 * A2)


print("\n========================================")
print("TWIST COMPATIBILITY CHECK")
print("========================================")
print(f"theta'_1 = {theta_1:.8e}")
print(f"theta'_2 = {theta_2:.8e}")
print(f"difference = {theta_1 - theta_2:.8e}")
# ================================================================
# 9. ACTUAL SHEAR FLOWS
# ================================================================

q_front = q_b_panel[front_panels] + q_0_1

q_rear = -q_b_panel[rear_panels] + q_0_2

print("\n========================================")
print("ACTUAL CELL SHEAR FLOW")
print("========================================")
print(
    f"Front cell q range: "
    f"{np.min(q_front):.8e} to {np.max(q_front):.8e} N/m"
)

print(
    f"Rear cell q range:  "
    f"{np.min(q_rear):.8e} to {np.max(q_rear):.8e} N/m"
)

print(f"Web q = {q_web:.8e} N/m")
# ================================================================
# ACTUAL GLOBAL SHEAR FLOW ON AIRFOIL SKIN
# ================================================================

q_skin_global = np.zeros(n)

# ------------------------------------------------
# Front cell
# ------------------------------------------------

for i in front_panels:
    q_skin_global[i] = q_b_panel[i] + q_0_1


# ------------------------------------------------
# Rear cell
# ------------------------------------------------

for i in rear_panels:
    q_skin_global[i] = q_b_panel[i] - q_0_2


# ------------------------------------------------
# Trailing-edge cut panel
# ------------------------------------------------

# q_b = 0 at the cut
q_skin_global[te_cut_panel] = -q_0_2
# ================================================================
# MOMENT OF SKIN SHEAR FLOW ABOUT CENTROID
# ================================================================

M_skin = 0.0

for i in range(n):

    j = (i + 1) % n

    M_panel = q_skin_global[i] * (
        x_centroid[i] * y_centroid[j]
        - y_centroid[i] * x_centroid[j]
    )

    M_skin += M_panel


print("\n========================================")
print("SKIN SHEAR-FLOW MOMENT")
print("========================================")
print(f"M_skin = {M_skin:.8e} N m")
# ================================================================
# MOMENT OF WEB SHEAR FLOW
# ================================================================

x_web = (
    0.5 * (
        x[upper_spar_index]
        + x[lower_spar_index]
    )
    - x_centroid_pt
)

M_web = q_web * x_web * spar_height


print("\n========================================")
print("WEB SHEAR-FLOW MOMENT")
print("========================================")
print(f"x_web = {x_web:.8e} m")
print(f"M_web = {M_web:.8e} N m")
# ================================================================
# TOTAL INTERNAL MOMENT
# ================================================================

M_internal = M_skin + M_web

print("\n========================================")
print("TOTAL INTERNAL MOMENT")
print("========================================")
print(f"M_skin     = {M_skin:.8e} N m")
print(f"M_web      = {M_web:.8e} N m")
print(f"M_internal = {M_internal:.8e} N m")
# ================================================================
# SHEAR CENTER LOCATION
# ================================================================

eta = -M_internal / V

print("\n========================================")
print("SHEAR CENTER")
print("========================================")
print(f"eta = {eta:.8e} m")
print(f"eta = {eta*1000:.4f} mm")
x_SC = x_centroid_pt + eta

print(f"x_centroid = {x_centroid_pt:.8e} m")
print(f"x_SC       = {x_SC:.8e} m")
print(f"x_SC/c     = {x_SC/chord:.6f}")
x_CP = Xcp*chord
T_CP = V * (x_CP - x_SC)

print("\n========================================")
print("CENTER-OF-PRESSURE TORQUE")
print("========================================")
print(f"x_CP = {x_CP:.8e} m")
print(f"x_SC = {x_SC:.8e} m")
print(f"T_CP = {T_CP:.8e} N m")
# ================================================================
# ADDITIONAL CELL SHEAR FLOWS DUE TO T_CP
# ================================================================

# ------------------------------------------------
# Torque equilibrium
# ------------------------------------------------

a11_T = 2.0 * A1
a12_T = 2.0 * A2

b1_T = T_CP


# ------------------------------------------------
# Twist compatibility
# ------------------------------------------------

a21_T = (
    (R1 + Rw) / (2.0 * A1)
    + Rw / (2.0 * A2)
)

a22_T = (
    -Rw / (2.0 * A1)
    - (R2 + Rw) / (2.0 * A2)
)

b2_T = 0.0


# ------------------------------------------------
# Solve
# ------------------------------------------------

M_T = np.array([
    [a11_T, a12_T],
    [a21_T, a22_T]
])

rhs_T = np.array([
    b1_T,
    b2_T
])

q_01_T, q_02_T = np.linalg.solve(M_T, rhs_T)


print("\n========================================")
print("TORQUE-INDUCED CELL FLOWS")
print("========================================")
print(f"q_01_T = {q_01_T:.8e} N/m")
print(f"q_02_T = {q_02_T:.8e} N/m")
# ================================================================
# TORQUE-INDUCED WEB SHEAR FLOW
# ================================================================

q_web_T = q_01_T - q_02_T

print("\n========================================")
print("TORQUE-INDUCED WEB FLOW")
print("========================================")
print(f"q_web_T = {q_web_T:.8e} N/m")
# ================================================================
# SHEAR FLOW CAUSED ONLY BY T_CP
# ================================================================

q_T = np.zeros(n)

# ------------------------------------------------
# Front cell
# ------------------------------------------------

for i in front_panels:
    q_T[i] = q_01_T


# ------------------------------------------------
# Rear cell
# ------------------------------------------------

for i in rear_panels:
    q_T[i] = -q_02_T


# ------------------------------------------------
# Trailing-edge cut panel
# ------------------------------------------------

q_T[te_cut_panel] = -q_02_T
# ================================================================
# FINAL SHEAR FLOW
# ================================================================

q_final = q_skin_global + q_T

print("\n========================================")
print("FINAL SHEAR FLOW")
print("========================================")
print(
    f"q_final range = "
    f"{np.min(q_final):.8e} to {np.max(q_final):.8e} N/m"
)
# ================================================================
# FINAL WEB SHEAR FLOW
# ================================================================

q_web_final = q_web + q_web_T

print(f"q_web_final = {q_web_final:.8e} N/m")
# ================================================================
# COLOURMAP OF SHEAR-FLOW DISTRIBUTION
# ================================================================

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

# ------------------------------------------------
# 1. SHEAR FLOW TO PLOT
# ------------------------------------------------
#
# Use the final shear flow if it has already been
# calculated.
#
# q_final = q_skin_global + q_T
#
# If you have not added the CP-induced torque yet,
# use:
#
# q_final = q_skin_global
# ------------------------------------------------

# q_final = q_skin_global + q_T


# ------------------------------------------------
# 2. CREATE LINE SEGMENTS FOR EVERY AIRFOIL PANEL
# ------------------------------------------------

segments = []

for i in range(n):

    j = (i + 1) % n

    segments.append([
        [x[i], y[i]],
        [x[j], y[j]]
    ])


segments = np.asarray(segments)


# ------------------------------------------------
# 3. COLOUR NORMALISATION
# ------------------------------------------------

q_min = np.min(q_final)
q_max = np.max(q_final)

norm = Normalize(
    vmin=q_min,
    vmax=q_max
)


# ------------------------------------------------
# 4. AIRFOIL COLOURMAP
# ------------------------------------------------

fig, ax = plt.subplots(figsize=(12, 6))

lc = LineCollection(
    segments,
    cmap="turbo",
    norm=norm,
    linewidth=4
)

lc.set_array(q_final)

ax.add_collection(lc)


# ------------------------------------------------
# 5. DRAW THE SPAR
# ------------------------------------------------

ax.plot(
    [
        x[upper_spar_index],
        x[lower_spar_index]
    ],
    [
        y[upper_spar_index],
        y[lower_spar_index]
    ],
    color="black",
    linewidth=3,
    label="C-section spar",
    zorder=5
)


# ------------------------------------------------
# 6. DRAW STRINGERS
# ------------------------------------------------

# Upper stringers

ax.scatter(
    x[upper_stringer_indices],
    y[upper_stringer_indices],
    s=80,
    marker="o",
    facecolors="white",
    edgecolors="black",
    linewidths=2,
    label="Stringers",
    zorder=6
)


# Lower stringers

ax.scatter(
    x[lower_stringer_indices],
    y[lower_stringer_indices],
    s=80,
    marker="o",
    facecolors="white",
    edgecolors="black",
    linewidths=2,
    zorder=6
)


# ------------------------------------------------
# 7. MARK SHEAR CENTER
# ------------------------------------------------

ax.scatter(
    x_SC,
    y_centroid_pt,
    s=120,
    marker="x",
    color="black",
    linewidths=3,
    label="Shear center",
    zorder=7
)


# ------------------------------------------------
# 8. MARK CENTROID
# ------------------------------------------------

ax.scatter(
    x_centroid_pt,
    y_centroid_pt,
    s=100,
    marker="+",
    color="black",
    linewidths=3,
    label="Centroid",
    zorder=7
)


# ------------------------------------------------
# 9. COLOURBAR
# ------------------------------------------------

cbar = fig.colorbar(lc, ax=ax)

cbar.set_label(
    "Shear flow, q [N/m]",
    rotation=270,
    labelpad=20
)


# ------------------------------------------------
# 10. AXES / FORMAT
# ------------------------------------------------

ax.set_aspect("equal", adjustable="box")

ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")

ax.set_title(
    "Two-Cell Airfoil Shear-Flow Distribution"
)

ax.grid(True, alpha=0.25)

ax.legend(
    loc="best"
)

# plt.tight_layout()
plt.show()
# ================================================================
# COLOURMAP OF SHEAR-FLOW DISTRIBUTION
# ================================================================

font_size = 18
title_size = 18

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

# ------------------------------------------------
# 1. SHEAR FLOW
# ------------------------------------------------

# q_final = q_skin_global + q_T


# ------------------------------------------------
# 2. CREATE LINE SEGMENTS
# ------------------------------------------------

segments = []

for i in range(n):

    j = (i + 1) % n

    segments.append([
        [x[i], y[i]],
        [x[j], y[j]]
    ])

segments = np.asarray(segments)


# ------------------------------------------------
# 3. COLOUR NORMALISATION
# ------------------------------------------------

q_min = np.min(q_final)
q_max = np.max(q_final)

norm = Normalize(
    vmin=q_min,
    vmax=q_max
)


# ------------------------------------------------
# 4. CREATE FIGURE
# ------------------------------------------------

fig, ax = plt.subplots(
    figsize=(12, 7),
    constrained_layout=True
)


# ------------------------------------------------
# 5. SHEAR-FLOW COLOURMAP
# ------------------------------------------------

lc = LineCollection(
    segments,
    cmap="turbo",
    norm=norm,
    linewidth=5,
    capstyle="round"
)

lc.set_array(q_final)

ax.add_collection(lc)


# ------------------------------------------------
# 6. DRAW SPAR
# ------------------------------------------------

ax.plot(
    [
        x[upper_spar_index],
        x[lower_spar_index]
    ],
    [
        y[upper_spar_index],
        y[lower_spar_index]
    ],
    color="black",
    linewidth=3,
    label="C-section spar",
    zorder=5
)


# ------------------------------------------------
# 7. DRAW STRINGERS
# ------------------------------------------------

ax.scatter(
    x[upper_stringer_indices],
    y[upper_stringer_indices],
    s=65,
    marker="o",
    facecolors="white",
    edgecolors="black",
    linewidths=1.8,
    label="Stringers",
    zorder=6
)

ax.scatter(
    x[lower_stringer_indices],
    y[lower_stringer_indices],
    s=65,
    marker="o",
    facecolors="white",
    edgecolors="black",
    linewidths=1.8,
    zorder=6
)


# ------------------------------------------------
# 8. SHEAR CENTRE
# ------------------------------------------------

ax.scatter(
    x_SC,
    y_centroid_pt,
    s=130,
    marker="x",
    color="black",
    linewidths=3,
    label="Shear center",
    zorder=7
)


# ------------------------------------------------
# 9. CENTROID
# ------------------------------------------------

ax.scatter(
    x_centroid_pt,
    y_centroid_pt,
    s=110,
    marker="+",
    color="black",
    linewidths=3,
    label="Centroid",
    zorder=7
)


# ------------------------------------------------
# 10. COLOURBAR
# ------------------------------------------------

cbar = fig.colorbar(
    lc,
    ax=ax,
    orientation="horizontal",
    pad=0.02
)

cbar.set_label(
    "Shear flow, $q$ [N/m]",
    fontsize=font_size,
    labelpad=15
)

cbar.ax.tick_params(
    labelsize=font_size
)


# ------------------------------------------------
# 11. AXES
# ------------------------------------------------

ax.set_aspect(
    "equal",
    adjustable="box"
)

ax.set_xlabel(
    "x [m]",
    fontsize=font_size
)

ax.set_ylabel(
    "y [m]",
    fontsize=font_size
)

ax.tick_params(
    axis="both",
    labelsize=font_size
)


# ------------------------------------------------
# 12. TITLE
# ------------------------------------------------

ax.set_title(
    "Two-Cell Airfoil Shear-Flow Distribution",
    fontsize=title_size,
    pad=12
)


# ------------------------------------------------
# 13. GRID
# ------------------------------------------------

ax.grid(
    True,
    alpha=0.25,
    linewidth=0.8
)


# ------------------------------------------------
# 14. LEGEND
# ------------------------------------------------

ax.legend(
    fontsize=font_size,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=True,
    framealpha=0.9
)


# ------------------------------------------------
# 15. LIMITS
# ------------------------------------------------

margin_x = 0.02 * (np.max(x) - np.min(x))
margin_y = 0.05 * (np.max(y) - np.min(y))

ax.set_xlim(
    np.min(x) - margin_x,
    np.max(x) + margin_x
)

ax.set_ylim(
    np.min(y) - margin_y,
    np.max(y) + margin_y
)


plt.show()
# Position from Wing Root
wing_root_tip_len = 850e-3
aileron_inner_edge = 550e-3     # From Wing Root
# rib_spanwise_position = np.array([0, 100, 200, 300, 425, 550, 700, 850])*1e-3

# rib_spanwise_position = np.array([0, 125, 250, 400, 550, 700, 850])*1e-3

rib_spanwise_position = np.array([0, 250, 550, 850])*1e-3
# rib_spanwise_position = np.array([0, 850])*1e-3
# Finding the Panel Length between Stringers
# Note this is without spar case

# print(stringers_spacing)
stringers_idx = np.concatenate((np.flip(upper_stringer_indices), lower_stringer_indices))
print(stringers_idx)
# stringers_idx = np.insert(stringers_idx, [4], [42])
stringer_spar_idx = np.insert(stringers_idx, [1], [25, 42, 57])

num_spar_stringer = num_stringers + 2

# print(np.array(range(1, num_stringers)))
print(stringer_spar_idx)

# num_stringer_spar = num_stringers + 2
stringers_spacing = np.zeros(num_spar_stringer + 2)

for i in range(1, num_spar_stringer + 1):
    stringers_spacing[i] = np.sum(ds[:stringer_spar_idx[i] + 1]) - np.sum(ds[:stringer_spar_idx[i-1] + 1])

stringers_spacing[0] = np.sum(ds[:stringer_spar_idx[0] + 1])
stringers_spacing[num_spar_stringer + 1] = np.sum(ds) - np.sum(ds[:stringer_spar_idx[-1] + 1])
# print(stringers_idx)
print(stringers_spacing)

print(np.sum(ds) - np.sum(stringers_spacing))
# Wing Parameters
a = rib_spanwise_position[1:] - rib_spanwise_position[:-1]
b = stringers_spacing
t = 5e-4

# Material Aluminum
E = 69e9 
nu = 0.33

Ra = 1
Rb = 1

n_fos = 1.5
n_vn = 3.5
n_fatigue = 1.5
k_stress_conc_factor = 1.5


def Pcr(a_, b_):
    Pcr_array = np.zeros((len(a_), len(b_)))
    for i in range(len(a_)):
        for j in range(len(b_)):
            Kss = 5.34 + 4*(b_[j]/a_[i])**2
            Pcr_array[i, j] = Kss*(np.pi**2*E/(12*(1 - nu**2)))*(t/b_[j])**2*(Ra + (Rb - Ra)*(b_[j]/a_[i])**2/2)
    return Pcr_array/(n_fos*n_vn*n_fatigue*k_stress_conc_factor)
Pcr = Pcr(a, b)
qmax_panel = np.zeros(len(stringer_spar_idx) + 1)

print(np.array(range(1, len(qmax_panel) - 1)))

for i in range(1, len(qmax_panel) - 1):

    panel_start_idx = stringer_spar_idx[i - 1]
    panel_end_idx = stringer_spar_idx[i]
    print(i, panel_start_idx, panel_end_idx, "length =", panel_end_idx - panel_start_idx)

    qmax_panel[i] = np.max(np.abs(q_final[panel_start_idx:panel_end_idx]))

qmax_panel[0] = np.max(np.abs(q_final[:stringer_spar_idx[0]]))
qmax_panel[-1] = np.max(np.abs(q_final[stringer_spar_idx[-1]:]))

tau_max_panel = qmax_panel/skin_thickness
print(Pcr)
print(tau_max_panel)
Buck_output = np.zeros_like(Pcr)

for i in range(len(a)):
    for j in range(len(b)):
        Buck_output[i, j] = Pcr[i, j]/tau_max_panel[j]

Buck_status = np.where(Buck_output < 1, 0, 1)

# 0: Buckles, 1: Safe
Buck_status