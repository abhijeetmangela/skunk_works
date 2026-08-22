"""
S7055 Airfoil Idealised Structural Analysis
-------------------------------------------
Modified from the original FX63137 script so the airfoil discretisation is
determined automatically from the .dat file.

The original script assumed exactly 97 coordinate points / 96 panels and
therefore contained hard-coded indices such as 18, 19, 48, 49, 78 and 96.
S7055 in s7055(2).dat contains 81 coordinate points, so those indices cannot
be used directly. This version finds the leading edge and stringer locations
automatically, while placing the C-section spar exactly at x = 0.070 m.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# USER PARAMETERS
# =============================================================================

FILENAME = "s7055.dat"       # Selig-format airfoil file
c = 0.234                       # chord [m]
Vy = 1000.0                     # vertical shear force [N]
a_param = 0.20                  # rib spacing [m]

# Stringers are specified as fractions of the distance along each surface.
# 0.15 and 0.30 reproduce approximately the old MATLAB locations 8 and 15
# for the 49-point surface used by the FX63137 model.
STRINGER_FRACTIONS = [0.15, 0.30]

# Stringer cross-sectional area [m^2]
str_area = 10e-3 * 1e-3

# Skin thickness [m]
t_skin = 0.5e-3

# C-section spar thickness [m]
t_spar = 1e-3

# C-section spar geometry
# The web is vertical and the two flanges point in the +x direction.
SPAR_X = 0.070                 # fixed spar location [m]
C_FLANGE_LENGTH = 30e-3       # flange length [m]
C_SPAR_HEIGHT_EST = 28e-3     # estimated web height [m]

# Use the actual airfoil thickness at SPAR_X as the web length. This keeps the
# web between the upper and lower skins. Set False to use the 28 mm estimate.
USE_AIRFOIL_WEB_LENGTH = True

# =============================================================================
# LOAD AIRFOIL
# =============================================================================

def load_airfoil(filename):
    rows = []
    with open(filename, "r") as f:
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            try:
                rows.append((float(p[0]), float(p[1])))
            except ValueError:
                continue

    if len(rows) < 10:
        raise ValueError("Could not read enough numeric x-y points from the airfoil file.")

    return np.asarray(rows, dtype=float)


raw = load_airfoil(FILENAME)

# Selig files are normally TE -> upper surface -> LE -> lower surface -> TE.
# We insert the exact SPAR_X station into both surfaces so the spar is placed
# exactly at x = 0.070 m rather than at the nearest .dat coordinate.

if not (0.0 < SPAR_X < c):
    raise ValueError(f"SPAR_X must lie between 0 and chord ({c:.6f} m).")

SPAR_X_NORM = SPAR_X / c

le0_raw = int(np.argmin(raw[:, 0]))
upper_raw = raw[:le0_raw + 1].copy()
lower_raw = raw[le0_raw:].copy()

def insert_x_point(surface, x_target):
    """Insert an interpolated point at x_target into one airfoil surface."""
    x = surface[:, 0]

    # Already present.
    exact = np.where(np.isclose(x, x_target, atol=1e-12))[0]
    if exact.size:
        return surface, int(exact[0])

    # Find the segment crossing x_target.
    for j in range(len(surface) - 1):
        x1, x2 = x[j], x[j + 1]
        if (x1 - x_target) * (x2 - x_target) < 0.0:
            y_target = np.interp(x_target, [x1, x2],
                                 [surface[j, 1], surface[j + 1, 1]])
            point = np.array([[x_target, y_target]])
            out = np.vstack((surface[:j + 1], point, surface[j + 1:]))
            return out, j + 1

    raise ValueError(
        f"SPAR_X/c = {x_target:.6f} is outside the supplied airfoil surface."
    )

upper_raw, u_insert = insert_x_point(upper_raw, SPAR_X_NORM)
lower_raw, l_insert = insert_x_point(lower_raw, SPAR_X_NORM)

# Rebuild the complete Selig-order coordinate array, sharing the LE only once.
raw = np.vstack((upper_raw, lower_raw[1:]))
le0 = len(upper_raw) - 1
L = len(raw)

print(f"Airfoil points after spar insertion : {L}")
print(f"Leading-edge point                : {le0 + 1} / {L}")
print(f"Fixed spar location               : x = {SPAR_X:.6f} m (x/c = {SPAR_X_NORM:.4f})")

# Scale coordinates by chord.
data = np.zeros((L + 1, 2))
data[1:, :] = raw * c

# =============================================================================
# PANEL LENGTHS AND MIDPOINTS
# =============================================================================

dis = np.zeros(L + 1)
mid_data = np.zeros((L + 1, 2))

for i in range(1, L):
    dis[i] = np.hypot(data[i + 1, 0] - data[i, 0],
                      data[i + 1, 1] - data[i, 1])
    mid_data[i] = 0.5 * (data[i] + data[i + 1])

# =============================================================================
# UPPER / LOWER SURFACES AND FIXED SPAR
# =============================================================================

# Upper surface: TE -> LE
upper_idx = np.arange(1, le0 + 2)
upper_xy = data[upper_idx]

# Lower surface: LE -> TE
lower_idx = np.arange(le0 + 1, L + 1)
lower_xy = data[lower_idx]

# Because the exact SPAR_X point was inserted above, the spar indices are exact.
u_spar = int(upper_idx[u_insert])
l_spar = int(lower_idx[l_insert])

# Extract the exact upper/lower y-coordinates at SPAR_X.
y_spar_top_airfoil = float(data[u_spar, 1])
y_spar_bottom_airfoil = float(data[l_spar, 1])
web_length_airfoil = abs(y_spar_top_airfoil - y_spar_bottom_airfoil)

t_max = web_length_airfoil
x_spar = SPAR_X

print(f"Spar upper point                   : {u_spar}")
print(f"Spar lower point                   : {l_spar}")
print(f"Airfoil web length at spar         : {web_length_airfoil*1000:.3f} mm")

# =============================================================================
# C-SECTION SPAR GEOMETRY
# =============================================================================

# The requested C-section is defined by:
#   - vertical web at the automatically selected spar x-location
#   - upper flange of length 30 mm extending in +x
#   - lower flange of length 30 mm extending in +x
#   - nominal overall web height of 28 mm
#
# A C-section line model is used for the section-property calculation.
# Its top and bottom flange centrelines are placed symmetrically about the
# midpoint of the local airfoil thickness.

spar_mid_y = 0.5 * (y_spar_top_airfoil + y_spar_bottom_airfoil)

# The requested 28 mm value is an estimate; by default the actual web length
# between the S7055 skins at x = 0.070 m is used for the C-section web.
if USE_AIRFOIL_WEB_LENGTH:
    c_spar_height = web_length_airfoil
else:
    c_spar_height = C_SPAR_HEIGHT_EST

spar_x = SPAR_X
spar_y_top = spar_mid_y + 0.5 * c_spar_height
spar_y_bottom = spar_mid_y - 0.5 * c_spar_height

# C-section corner/end points.
C_web_bottom = np.array([spar_x, spar_y_bottom])
C_web_top = np.array([spar_x, spar_y_top])
C_flange_top_end = np.array([spar_x + C_FLANGE_LENGTH, spar_y_top])
C_flange_bottom_end = np.array([spar_x + C_FLANGE_LENGTH, spar_y_bottom])

print(f"C-section flange length           : {C_FLANGE_LENGTH*1000:.1f} mm")
print(f"Estimated web height              : {C_SPAR_HEIGHT_EST*1000:.1f} mm")
print(f"Actual airfoil web length         : {web_length_airfoil*1000:.3f} mm")
print(f"C-section web length used         : {c_spar_height*1000:.3f} mm")
print(f"Flanges extend in +x by           : {C_FLANGE_LENGTH*1000:.1f} mm")
if abs(c_spar_height - C_SPAR_HEIGHT_EST) > 0.5e-3:
    print(
        f"NOTE: actual web length differs from the 28 mm estimate by "
        f"{abs(c_spar_height-C_SPAR_HEIGHT_EST)*1000:.3f} mm."
    )

# Thin-walled C-section area and centroid.
A_web = t_spar * c_spar_height
A_flange = t_spar * C_FLANGE_LENGTH
A_C = A_web + 2.0 * A_flange

# Centroid of the C-section measured from the spar web centreline.
x_C_local = (A_web * 0.0 + 2.0 * A_flange * (0.5 * C_FLANGE_LENGTH)) / A_C

y_C_local = spar_mid_y

x_C = spar_x + x_C_local

print(f"C-section area          : {A_C*1e6:.3f} mm^2")
print(f"C-section centroid      : x = {x_C:.6f} m, y = {y_C_local:.6f} m")

# -----------------------------------------------------------------------------
# C-SPAR BOOM IDEALISATION
# -----------------------------------------------------------------------------
# The C-section is represented by FOUR concentrated booms:
#
#   B1 o====================o B2
#      |                    |
#      |       web          |
#   B3 o====================o B4
#
# B1/B3 = web-side corners
# B2/B4 = flange ends (30 mm in +x)
#
# To avoid double-counting material, the web and each flange are split between
# their two end booms. This preserves the total C-section area and its first
# moments, while providing explicit flange-end booms for the idealised model.
#
# The exact thin-walled C-section geometry above is still retained for reporting
# its geometric area/centroid. The boom representation below is used in the
# shear-flow idealisation.
A_web_half = 0.5 * A_web
A_flange_half = 0.5 * A_flange

# A_boom_web_top = A_web_half + A_flange_half
# A_boom_flange_top = A_flange_half
# A_boom_web_bottom = A_web_half + A_flange_half
# A_boom_flange_bottom = A_flange_half

A_boom_flange_top = t_spar * C_FLANGE_LENGTH* (2+ 1)/6
A_boom_web_top = t_spar * c_spar_height* (2+ 1)/6 + t_spar*c_spar_height* (2+(spar_y_bottom/spar_y_top))/6

A_boom_flange_bottom = t_spar * C_FLANGE_LENGTH* (2+ 1)/6
A_boom_web_bottom = t_spar * c_spar_height* (2+ 1)/6 + t_spar*c_spar_height* (2+(spar_y_top/spar_y_bottom))/6
spar_boom_points = np.array([
    [spar_x, spar_y_top],
    [spar_x + C_FLANGE_LENGTH, spar_y_top],
    [spar_x, spar_y_bottom],
    [spar_x + C_FLANGE_LENGTH, spar_y_bottom],
], dtype=float)

spar_boom_areas = np.array([
    A_boom_web_top,
    A_boom_flange_top,
    A_boom_web_bottom,
    A_boom_flange_bottom,
], dtype=float)

print("\nC-SPAR BOOMS")
for name, pt, area in zip(
    ["Upper web corner", "Upper flange end",
     "Lower web corner", "Lower flange end"],
    spar_boom_points,
    spar_boom_areas,
):
    print(
        f"{name:20s}: x = {pt[0]:.6f} m, y = {pt[1]:.6f} m, "
        f"A = {area*1e6:.3f} mm^2"
    )

print(f"Total C-spar boom area : {np.sum(spar_boom_areas)*1e6:.3f} mm^2")
print(f"Geometric C-section area: {A_C*1e6:.3f} mm^2")

# =============================================================================
# CENTROID AND SECOND MOMENTS
# =============================================================================

sum_dis = np.sum(dis[1:L])

cg1 = (
    np.sum(mid_data[1:L, 0] * dis[1:L] * t_skin)
    + A_C * x_C
) / (
    np.sum(dis[1:L] * t_skin) + A_C
)

cg2 = (
    np.sum(mid_data[1:L, 1] * dis[1:L] * t_skin)
    + A_C * y_C_local
) / (
    np.sum(dis[1:L] * t_skin) + A_C
)

print(f"Centroid: x = {cg1:.6f} m, y = {cg2:.6f} m")

# Exact thin-walled C-section contribution to the boom-level section properties.
# I_x/I_y here follow the same convention as the rest of the original script:
#   I_x = integral(y^2 dA), I_y = integral(x^2 dA), I_xy = integral(x y dA).

x_web = spar_x - cg1

y_top = spar_y_top - cg2

y_bottom = spar_y_bottom - cg2
x_flange_c = spar_x + 0.5 * C_FLANGE_LENGTH - cg1

I_xy = (
    A_web * x_web * (0.5 * (spar_y_top + spar_y_bottom) - cg2)
    + A_flange * x_flange_c * y_top
    + A_flange * x_flange_c * y_bottom
)

# Web about the vertical y-axis plus flange contributions.
I_y = A_web * x_web**2 + 2.0 * A_flange * x_flange_c**2
I_x = (
    A_web * ((c_spar_height**2) / 12.0 +
             (0.5 * (spar_y_top + spar_y_bottom) - cg2)**2)
    + A_flange * y_top**2
    + A_flange * y_bottom**2
)

# Re-centre the contour coordinates.
data[:, 0] -= cg1
data[:, 1] -= cg2

for i in range(1, L):
    x1, y1 = data[i]
    x2, y2 = data[i + 1]

    I_xy += 0.25 * (x1 + x2) * (y1 + y2) * dis[i] * t_skin
    I_y += 0.25 * (x1 + x2)**2 * dis[i] * t_skin
    I_x += 0.25 * (y1 + y2)**2 * dis[i] * t_skin

DR = I_y * I_x - I_xy**2
Kxy = I_xy / DR
K_y = I_y / DR
K_x = I_x / DR
tan_alpha = -I_xy / I_y

# =============================================================================
# BOOM AREAS
# =============================================================================

pd = np.zeros(L + 1)
p_d = np.zeros(L + 1)

for i in range(1, L):
    px = 0.5 * (data[i, 0] + data[i + 1, 0])
    py = 0.5 * (data[i, 1] + data[i + 1, 1])

    pd[i] = abs(py - tan_alpha * px) / np.sqrt(1 + tan_alpha**2)

for i in range(1, L + 1):
    p_d[i] = abs(data[i, 1] - tan_alpha * data[i, 0]) / np.sqrt(1 + tan_alpha**2)

a_boom = np.zeros(L + 1)

for i in range(1, L):
    # Protect against an exactly-zero neutral-axis distance.
    if pd[i] < 1e-14:
        continue

    if i < L - 1 and pd[i + 1] > 1e-14:
        a_boom[i] += (1/6) * dis[i] * t_skin * (2 + p_d[i + 1]/pd[i])
        a_boom[i + 1] += (1/6) * dis[i] * t_skin * (2 + p_d[i]/pd[i + 1])

    elif i == L - 1 and pd[1] > 1e-14:
        a_boom[i] += (1/6) * dis[i] * t_skin * (2 + p_d[1]/pd[i])
        a_boom[1] += (1/6) * dis[i] * t_skin * (2 + p_d[i]/pd[1])

# Add the two web-side C-spar booms to the airfoil contour.
# The two flange-end booms are internal spar booms and are kept separately in
# spar_boom_points / spar_boom_areas so they are not incorrectly placed on the skin.
a_boom[u_spar] += A_boom_web_top
a_boom[l_spar] += A_boom_web_bottom

# =============================================================================
# OPEN SHEAR FLOW
# =============================================================================

def find_shear_spar(
    data2,
    a_boom2,
    spar_boom_points2,
    spar_boom_areas2,
    Vy,
    L,
    u_spar,
    l_spar,
    le_point,
    dis,
):
    """
    Generic version of the original two-cell shear-flow calculation.

    Cell 1:
        upper TE -> upper spar -> spar -> lower spar -> lower TE

    Cell 2:
        upper spar -> upper LE -> lower LE -> lower spar -> spar
    """

    # Combine the skin booms with the four explicit C-spar booms.
    boom_x = np.concatenate((data2[1, 1:L], spar_boom_points2[:, 0]))
    boom_y = np.concatenate((data2[2, 1:L], spar_boom_points2[:, 1]))
    boom_A = np.concatenate((a_boom2[1:L], spar_boom_areas2))

    total_A = np.sum(boom_A)

    cg_new_x = np.sum(boom_x * boom_A) / total_A
    cg_new_y = np.sum(boom_y * boom_A) / total_A

    x = data2[1] - cg_new_x
    y = data2[2] - cg_new_y

    spar_x_rel = spar_boom_points2[:, 0] - cg_new_x
    spar_y_rel = spar_boom_points2[:, 1] - cg_new_y

    Ixy2 = (
        np.sum(x[1:L] * y[1:L] * a_boom2[1:L])
        + np.sum(spar_x_rel * spar_y_rel * spar_boom_areas2)
    )
    Ix2 = (
        np.sum(y[1:L]**2 * a_boom2[1:L])
        + np.sum(spar_y_rel**2 * spar_boom_areas2)
    )
    Iy2 = (
        np.sum(x[1:L]**2 * a_boom2[1:L])
        + np.sum(spar_x_rel**2 * spar_boom_areas2)
    )

    den2 = Ix2 * Iy2 - Ixy2**2
    Kxy2 = Ixy2 / den2
    Ky2 = Iy2 / den2

    Qy = np.zeros(L + 1)
    Qx = np.zeros(L + 1)
    qs = np.zeros(L + 1)

    # Upper surface: TE -> LE.
    for i in range(1, le_point):
        Qy[i + 1] = Qy[i] + a_boom2[i + 1] * x[i + 1]
        Qx[i + 1] = Qx[i] + a_boom2[i + 1] * y[i + 1]
        qs[i + 1] = Vy * (Kxy2 * Qy[i + 1] - Ky2 * Qx[i + 1])

    # Lower surface: LE -> TE.
    for i in range(le_point, L):
        Qy[i + 1] = Qy[i] + a_boom2[i + 1] * x[i + 1]
        Qx[i + 1] = Qx[i] + a_boom2[i + 1] * y[i + 1]
        qs[i + 1] = Vy * (Kxy2 * Qy[i + 1] - Ky2 * Qx[i + 1])

    # The spar joins two points on the contour. The open-section flow
    # discontinuity at the upper spar is used as the spar flow.
    q_spar = qs[u_spar - 1] - qs[u_spar]

    # Cell 1 skin lengths:
    cell1_upper = np.sum(dis[1:u_spar])
    cell1_lower = np.sum(dis[l_spar:L])

    # Cell 2 skin lengths:
    cell2_upper = np.sum(dis[u_spar:le_point])
    cell2_lower = np.sum(dis[le_point:l_spar])

    A_mat = np.array([
        [cell1_upper + cell1_lower + c_spar_height, -c_spar_height],
        [-c_spar_height, cell2_upper + cell2_lower + c_spar_height]
    ])

    B_vec = np.array([
        np.sum(qs[1:u_spar] * dis[1:u_spar])
        + np.sum(qs[l_spar:L] * dis[l_spar:L])
        + q_spar * c_spar_height,

        np.sum(qs[u_spar:le_point] * dis[u_spar:le_point])
        + np.sum(qs[le_point:l_spar] * dis[le_point:l_spar])
        - q_spar * c_spar_height
    ])

    q0 = np.linalg.solve(A_mat, B_vec)

    q_spar = q_spar + q0[0] - q0[1]

    return qs, q_spar, q0


data2 = np.zeros((3, L + 1))
data2[1] = data[:, 0]
data2[2] = data[:, 1]

qs, q_spar, q0 = find_shear_spar(
    data2, a_boom, spar_boom_points, spar_boom_areas,
    Vy, L, u_spar, l_spar, le0 + 1, dis
)

# Apply redundant cell flows.
shear = qs.copy()
shear[1:u_spar] += q0[0]
shear[u_spar:le0 + 1] += q0[1]
shear[le0 + 1:l_spar] += q0[1]
shear[l_spar:L] += q0[0]

# =============================================================================
# BUCKLING CHECK
# =============================================================================

dis2 = np.array([
    np.sum(dis[1:u_spar]),
    np.sum(dis[u_spar:le0 + 1]),
    np.sum(dis[le0 + 1:l_spar]),
    np.sum(dis[l_spar:L]),
    c_spar_height
])

D = (
    np.pi**2 * 70e9 * (t_skin**3) * 1.62
    / (12 * 2.5 * 1.5 * 1.5 * 1.5 * (1 - 0.33**2))
)

Kss = 5.34 + 4.0 * dis2**2 / a_param**2
N = Kss * D / dis2**2

q_max = np.array([
    np.max(np.abs(shear[1:u_spar])),
    np.max(np.abs(shear[u_spar:le0 + 1])),
    np.max(np.abs(shear[le0 + 1:l_spar])),
    np.max(np.abs(shear[l_spar:L])),
    abs(q_spar)
])

print("\n--- WITHOUT STRINGERS ---")
print("q_max [N/m] =", q_max)
print("N [N/m]     =", N)
print("N/q_max     =", N/q_max)

# =============================================================================
# PLOT: AIRFOIL WITH SPAR
# =============================================================================

x_plot = data[1:, 0] + cg1
y_plot = data[1:, 1] + cg2

fig, ax = plt.subplots(figsize=(14, 7))
ax.plot(x_plot, y_plot, linewidth=2)
# C-section spar: vertical web + two flanges extending in +x.
ax.plot(
    [C_web_bottom[0], C_web_top[0]],
    [C_web_bottom[1], C_web_top[1]],
    linewidth=3, label="C-section web"
)
ax.plot(
    [C_web_top[0], C_flange_top_end[0]],
    [C_web_top[1], C_flange_top_end[1]],
    linewidth=3
)
ax.plot(
    [C_web_bottom[0], C_flange_bottom_end[0]],
    [C_web_bottom[1], C_flange_bottom_end[1]],
    linewidth=3
)
ax.scatter(
    [data[u_spar, 0] + cg1, data[l_spar, 0] + cg1],
    [data[u_spar, 1] + cg2, data[l_spar, 1] + cg2],
    s=100, label="Web-side booms"
)

# Explicit flange-end booms.
ax.scatter(
    spar_boom_points[[1, 3], 0],
    spar_boom_points[[1, 3], 1],
    s=120, marker="s", label="Flange-end booms"
)
# Show the requested +x flange direction and dimensions.
ax.annotate(
    "+x flange",
    xy=(C_flange_top_end[0], C_flange_top_end[1]),
    xytext=(C_flange_top_end[0] + 0.008, C_flange_top_end[1] + 0.006),
    arrowprops=dict(arrowstyle="->")
)
ax.text(
    spar_x + 0.5 * C_FLANGE_LENGTH,
    spar_y_top + 0.002,
    f"30 mm flange",
    ha="center", fontsize=10
)
ax.text(
    spar_x + 0.002,
    spar_mid_y,
    f"C height = {c_spar_height*1000:.1f} mm",
    rotation=90, va="center", fontsize=10
)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("S7055 Airfoil with C-Section Spar")
ax.axis("equal")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.savefig("s7055_x070_Cspar_geometry.png", dpi=150)
plt.close(fig)

# =============================================================================
# PLOT: SHEAR FLOW WITHOUT STRINGERS
# =============================================================================

fig, ax = plt.subplots(figsize=(14, 7))
ax.plot(
    data[1:le0 + 1, 0] + cg1,
    -shear[1:le0 + 1],
    linewidth=2,
    label="Upper Surface"
)
ax.plot(
    data[le0 + 1:L, 0] + cg1,
    -shear[le0 + 1:L],
    linewidth=2,
    label="Lower Surface"
)
ax.set_xlabel("x [m]")
ax.set_ylabel("q_s [N/m]")
ax.set_title("S7055 Shear Flow Distribution with Spar")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.savefig("s7055_x070_shear_flow_no_stringers.png", dpi=150)
plt.close(fig)

# =============================================================================
# ADD STRINGERS
# =============================================================================

a_boom2 = a_boom.copy()

# Convert fractions to valid 1-based upper-surface point indices.
# The first and last points are always the TE and LE boundaries.
u_stringers = [
    int(round(1 + f * (le0 - 0)))
    for f in STRINGER_FRACTIONS
]
u_stringers = sorted(set(
    max(2, min(le0, i)) for i in u_stringers
))

# Lower stringers use the same fractional position measured LE -> TE.
lower_count = L - le0
l_stringers = [
    int(round((le0 + 1) + f * (L - le0 - 1)))
    for f in STRINGER_FRACTIONS
]
l_stringers = sorted(set(
    max(le0 + 1, min(L - 1, i)) for i in l_stringers
))

for i in u_stringers:
    a_boom2[i] += str_area

for i in l_stringers:
    a_boom2[i] += str_area

print("\n--- STRINGER LOCATIONS ---")
print("Upper stringers:")
for i in u_stringers:
    print(f"  point {i:3d}, x/c = {(data[i,0]+cg1)/c:.4f}")

print("Lower stringers:")
for i in l_stringers:
    print(f"  point {i:3d}, x/c = {(data[i,0]+cg1)/c:.4f}")

qs2, q_spar2, q02 = find_shear_spar(
    data2, a_boom2, spar_boom_points, spar_boom_areas,
    Vy, L, u_spar, l_spar, le0 + 1, dis
)

shear2 = qs2.copy()
shear2[1:u_spar] += q02[0]
shear2[u_spar:le0 + 1] += q02[1]
shear2[le0 + 1:l_spar] += q02[1]
shear2[l_spar:L] += q02[0]

# =============================================================================
# PANEL / STRINGER BUCKLING SEGMENTS
# =============================================================================

# Upper boundary points: TE -> stringers -> spar -> LE
upper_bounds = sorted(set([1] + u_stringers + [u_spar, le0 + 1]))

# Lower boundary points: LE -> spar -> stringers -> TE
lower_bounds = sorted(set([le0 + 1] + l_stringers + [l_spar, L]))

segments = []
segment_names = []

# Upper segments
for j in range(len(upper_bounds) - 1):
    lo, hi = upper_bounds[j], upper_bounds[j + 1]
    segments.append(np.sum(dis[lo:hi]))
    segment_names.append(f"Upper {j + 1}")

# Lower segments
for j in range(len(lower_bounds) - 1):
    lo, hi = lower_bounds[j], lower_bounds[j + 1]
    segments.append(np.sum(dis[lo:hi]))
    segment_names.append(f"Lower {j + 1}")

segments.append(c_spar_height)
segment_names.append("Spar")

q_segment = []

for j in range(len(upper_bounds) - 1):
    lo, hi = upper_bounds[j], upper_bounds[j + 1]
    q_segment.append(np.max(np.abs(shear2[lo:hi])))

for j in range(len(lower_bounds) - 1):
    lo, hi = lower_bounds[j], lower_bounds[j + 1]
    q_segment.append(np.max(np.abs(shear2[lo:hi])))

q_segment.append(abs(q_spar2))

segments = np.asarray(segments)
q_segment = np.asarray(q_segment)

Kss2 = 5.34 + 4.0 * segments**2 / a_param**2
N2 = Kss2 * D / segments**2

print("\n--- WITH STRINGERS ---")
for name, length, q, n in zip(segment_names, segments, q_segment, N2):
    print(f"{name:12s}: length = {length:.5f} m, "
          f"q_max = {q:.3f} N/m, N = {n:.3f} N/m, N/q = {n/q:.3f}")

# =============================================================================
# PLOT: SECTION WITH STRINGERS
# =============================================================================

fig, ax = plt.subplots(figsize=(14, 7))
ax.plot(x_plot, y_plot, linewidth=2)

# Spar
# C-section spar: vertical web + two flanges extending in +x.
ax.plot(
    [C_web_bottom[0], C_web_top[0]],
    [C_web_bottom[1], C_web_top[1]],
    linewidth=3, label="C-section web"
)
ax.plot(
    [C_web_top[0], C_flange_top_end[0]],
    [C_web_top[1], C_flange_top_end[1]],
    linewidth=3
)
ax.plot(
    [C_web_bottom[0], C_flange_bottom_end[0]],
    [C_web_bottom[1], C_flange_bottom_end[1]],
    linewidth=3
)

# Upper stringers
for k, idx in enumerate(u_stringers, 1):
    ax.scatter(
        data[idx, 0] + cg1,
        data[idx, 1] + cg2,
        s=80
    )
    ax.text(
        data[idx, 0] + cg1,
        data[idx, 1] + cg2 + 0.004,
        f"U{k}",
        fontsize=11
    )

# Lower stringers
for k, idx in enumerate(l_stringers, 1):
    ax.scatter(
        data[idx, 0] + cg1,
        data[idx, 1] + cg2,
        s=80
    )
    ax.text(
        data[idx, 0] + cg1,
        data[idx, 1] + cg2 - 0.006,
        f"L{k}",
        fontsize=11
    )

ax.scatter(
    [data[u_spar, 0] + cg1, data[l_spar, 0] + cg1],
    [data[u_spar, 1] + cg2, data[l_spar, 1] + cg2],
    s=100
)
ax.scatter(
    spar_boom_points[[1, 3], 0],
    spar_boom_points[[1, 3], 1],
    s=120, marker="s", label="Flange-end booms"
)

ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("S7055 Section with C-Section Spar and Stringers")
ax.axis("equal")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.savefig("s7055_x070_section_with_stringers.png", dpi=150)
plt.close(fig)

# =============================================================================
# PLOT: SHEAR FLOW WITH STRINGERS
# =============================================================================

fig, ax = plt.subplots(figsize=(14, 7))
ax.plot(
    data[1:le0 + 1, 0] + cg1,
    -shear2[1:le0 + 1],
    linewidth=2,
    label="Upper Surface"
)
ax.plot(
    data[le0 + 1:L, 0] + cg1,
    -shear2[le0 + 1:L],
    linewidth=2,
    label="Lower Surface"
)
ax.set_xlabel("x [m]")
ax.set_ylabel("q_s [N/m]")
ax.set_title("S7055 Shear Flow Distribution with Spar and Stringers")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.savefig("s7055_x070_shear_flow_with_stringers.png", dpi=150)
plt.close(fig)

print("\nDone.")
print("Generated:")
print("  s7055_x070_Cspar_geometry.png")
print("  s7055_x070_shear_flow_no_stringers.png")
print("  s7055_x070_section_with_stringers.png")
print("  s7055_x070_shear_flow_with_stringers.png")
print("  Four C-spar booms included: 2 web-side corners + 2 flange ends.")
