import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# INPUT / INITIALISATION
# ============================================================

# Airfoil coordinate file
filename = "s7055.dat"

# Airfoil chord [m]
c = 0.265

# Skin thickness [m]
skin_thickness = 0.0005

# Stringer area [m^2]
str_area = 10.0 * 0.001 * 0.001

# Material / plate constants
E = 69e9
nu = 0.33

# Shear force [N]
Vy = 100.0

# Rib spacing [m]
a = 0.20

# Number of stringers on upper surface
iu = 4

# Upper stringer indices
# MATLAB indexing starts at 1, Python starts at 0.
# These values are therefore kept as Python indices.
# ============================================================
# STRINGER / PANEL INPUTS
# ============================================================

# Number of stringers
iu = 4
il = 4

# # MATLAB-style panel boundary indices:
# # 1 = first point
# # 49 = last point
# iu_ind_matlab = np.array([1, 13, 25, 49])
# il_ind_matlab = np.array([1, 13, 25, 49])

# # Convert MATLAB indices to Python indices
# iu_ind = iu_ind_matlab - 1
# il_ind = il_ind_matlab - 1

# Plot settings
font_size = 18
title_size = 24


# ============================================================
# FUNCTION: SHEAR FLOW
# ============================================================

def find_shear(data2, a_boom2, Vy, L):
    """
    Calculate the basic shear-flow distribution.

    Parameters
    ----------
    data2 : ndarray, shape (2, L)
        x and y coordinates.
    a_boom2 : ndarray
        Boom areas.
    Vy : float
        Applied vertical shear force [N].
    L : int
        Number of airfoil coordinate points.

    Returns
    -------
    shear : ndarray
        Shear flow distribution [N/m].
    """

    # Centroid of boom idealisation
    cg_new_x = np.sum(data2[0, :L-1] * a_boom2) / np.sum(a_boom2)
    cg_new_y = np.sum(data2[1, :L-1] * a_boom2) / np.sum(a_boom2)

    # Coordinates relative to boom centroid
    data3 = data2.copy()
    data3[0, :] -= cg_new_x
    data3[1, :] -= cg_new_y

    # Section properties
    Ixy2 = np.sum(
        data3[0, :L-1] *
        data3[1, :L-1] *
        a_boom2
    )

    Ix2 = np.sum(
        data3[1, :L-1] ** 2 *
        a_boom2
    )

    Iy2 = np.sum(
        data3[0, :L-1] ** 2 *
        a_boom2
    )

    den2 = Ix2 * Iy2 - Ixy2 ** 2

    Kxy2 = Ixy2 / den2
    Ky2 = Ix2 / den2

    # These are retained from the MATLAB formulation
    Kx2 = Iy2 / den2

    # First moments
    Qy2 = np.zeros(L - 1)
    Qx2 = np.zeros(L - 1)

    # Shear flow
    qs2 = np.zeros(L - 1)

    for i in range(L - 2):
        Qy2[i + 1] = (
            Qy2[i] +
            data3[0, i + 1] * a_boom2[i + 1]
        )

        Qx2[i + 1] = (
            Qx2[i] +
            data3[1, i + 1] * a_boom2[i + 1]
        )

        qs2[i + 1] = Vy * (
            Kxy2 * Qy2[i + 1] -
            Ky2 * Qx2[i + 1]
        )

    return qs2


# ============================================================
# READ AIRFOIL DATA
# ============================================================

data1 = np.loadtxt(filename)

# Convert nondimensional DAT coordinates to physical coordinates
data = data1 * c

L = len(data)
# mid = (L - 1) //2

# iu_ind = np.linspace(0,   mid,   iu + 2).round().astype(int)
# il_ind = np.linspace(mid, L - 1, il + 2).round().astype(int)

iu_ind = np.array([0, 12, 16, 20, 30, 42])   # pick whichever point indices you want
il_ind = np.array([42, 52, 61, 65, 69, 80])

stringer_idx = np.concatenate([iu_ind[1:-1], il_ind[1:-1]])
plt.figure(figsize=(10, 4))
plt.plot(data[:, 0], data[:, 1], 'k-', linewidth=1.5, label="Airfoil")
plt.scatter(data[stringer_idx, 0], data[stringer_idx, 1],
            s=80, c='red', zorder=5, label="Stringers")
 
for i in stringer_idx:
    plt.annotate(f"{data[i,0]:.3f}", (data[i,0], data[i,1]),
                 textcoords="offset points", xytext=(0, 8), fontsize=8, ha='center')
 
plt.axis("equal")
plt.grid(True, alpha=0.3)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Stringer Positions on Airfoil")
plt.legend()
plt.tight_layout()
plt.show()


# Allocate arrays
dis = np.zeros(L - 1)
mid_data = np.zeros((L - 1, 2))


# ============================================================
# PANEL LENGTHS AND MIDPOINTS
# ============================================================

for i in range(L - 1):

    dx = data[i + 1, 0] - data[i, 0]
    dy = data[i + 1, 1] - data[i, 1]

    dis[i] = np.sqrt(dx**2 + dy**2)

    mid_data[i, :] = (
        data[i, :] +
        data[i + 1, :]
    ) * 0.5


# ============================================================
# GEOMETRIC CENTROID
# ============================================================

cg = np.zeros(2)

cg[0] = (
    np.sum(mid_data[:, 0] * dis) /
    np.sum(dis)
)

cg[1] = (
    np.sum(mid_data[:, 1] * dis) /
    np.sum(dis)
)


# ============================================================
# SHIFT AIRFOIL TO GEOMETRIC CENTROID
# ============================================================

I_num = 0.0
I_den = 0.0

data[:, 0] -= cg[0]
data[:, 1] -= cg[1]


# ============================================================
# NEUTRAL AXIS / TAN(ALPHA)
# ============================================================

for i in range(L - 1):

    I_num += (
        (data[i, 0] + data[i + 1, 0]) *
        (data[i, 1] + data[i + 1, 1]) *
        0.25 *
        dis[i]
    )

    I_den += (
        (data[i, 1] + data[i + 1, 1]) *
        (data[i, 1] + data[i + 1, 1]) *
        0.25 *
        dis[i]
    )

tan_alpha = I_num / I_den


# ============================================================
# DISTANCE FROM NEUTRAL AXIS
# ============================================================

p_d = np.zeros(L - 1)

for i in range(L - 1):

    p_d[i] = abs(
        data[i, 1] -
        tan_alpha * data[i, 0]
    ) / np.sqrt(1 + tan_alpha**2)


# ============================================================
# SKIN BOOM AREAS
# ============================================================

a_boom = np.zeros(L - 1)

for i in range(L - 1):

    if i < L - 2:

        a_boom[i] += (
            (1 / 6) *
            dis[i] *
            skin_thickness *
            (2 + p_d[i + 1] / p_d[i])
        )

        a_boom[i + 1] += (
            (1 / 6) *
            dis[i] *
            skin_thickness *
            (2 + p_d[i] / p_d[i + 1])
        )

    else:

        a_boom[i] += (
            (1 / 6) *
            dis[i] *
            skin_thickness *
            (2 + p_d[0] / p_d[i])
        )

        a_boom[0] += (
            (1 / 6) *
            dis[i] *
            skin_thickness *
            (2 + p_d[i] / p_d[0])
        )


# ============================================================
# SHEAR FLOW WITHOUT STRINGERS
# ============================================================

data2 = data.T

qs = find_shear(
    data2,
    a_boom,
    Vy,
    L
)


# ============================================================
# PLOT: SHEAR FLOW WITHOUT STRINGERS
# ============================================================

mid = (L - 1) // 2

plt.figure(figsize=(12, 8))

plt.plot(
    data2[0, :mid] + cg[0],
    -qs[:mid],
    linewidth=2,
    label="Upper Surface"
)

plt.plot(
    data2[0, mid:L-1] + cg[0],
    -qs[mid:L-1],
    linewidth=2,
    label="Lower Surface"
)

plt.grid(True)
plt.xlabel("x [m]")
plt.ylabel(r"$q_s$ [N/m]")
plt.legend(fontsize=font_size)
plt.title(
    "Shear Flow Distribution without spar",
    fontsize=title_size,
    fontweight="bold"
)

plt.tick_params(labelsize=font_size)
plt.tight_layout()


# ============================================================
# CLOSED-SECTION CONSTANT SHEAR FLOW
# ============================================================

q0 = -np.sum(qs * dis) / np.sum(dis)

qs = qs + q0

shear = qs


# ============================================================
# SKIN BUCKLING CHECK
# ============================================================

D = (
    np.pi**2 *
    E *
    skin_thickness**3 *
    1.62
) / (
    12 *
    2.5 *
    1.5 *
    1.5 *
    1.5 *
    (1 - nu**2)
)

Kss = 5.34 + (4 * c**2) / a**2

N = (Kss * D) / c**2

max_shear = np.max(np.abs(shear))

print("\n========================================")
print("SKIN BUCKLING CHECK")
print("========================================")
print(f"Maximum shear flow = {max_shear:.6e} N/m")
print(f"Critical buckling load = {N:.6e} N/m")

if N > max_shear:
    print("No stringer needed without spar")
else:
    print("Stringer needed without spar")


# ============================================================
# ADD STRINGERS
# ============================================================

a_boom2 = a_boom.copy()

# Upper surface stringers
for i in range(iu):
    idx = iu_ind[i + 1]
    a_boom2[idx] += str_area

# Lower surface stringers
for i in range(il):
    idx = il_ind[i + 1]
    a_boom2[idx] += str_area


# ============================================================
# SHEAR FLOW WITH STRINGERS
# ============================================================

qs2 = find_shear(
    data2,
    a_boom2,
    Vy,
    L
)

q02 = -np.sum(qs2 * dis) / np.sum(dis)


# ============================================================
# PANEL SHEAR-FLOW CHECK
# ============================================================

q_max2 = []
dis2 = []

print("\n========================================")
print("UPPER STRINGER / PANEL RESULTS")
print("========================================")

for i in range(iu + 1):

    start = iu_ind[i]
    end = iu_ind[i + 1]

    q_panel = qs2[start:end] + q02

    q_max2.append(
        np.max(np.abs(q_panel))
    )

    dis2.append(
        np.sum(dis[start:end])
    )

    print(
        f"Panel {i + 1}: "
        f"x = {data[start, 0] + cg[0]:.6f} m"
    )


print("\n========================================")
print("LOWER STRINGER / PANEL RESULTS")
print("========================================")

for i in range(il + 1):

    start = il_ind[i]
    end = il_ind[i + 1] 

    q_panel = qs2[start:end] + q02

    q_max2.append(
        np.max(np.abs(q_panel))
    )

    dis2.append(
        np.sum(dis[start:end])
    )

    print(
        f"Panel {iu + 2 + i}: "
        f"x = {data[start, 0] + cg[0]:.6f} m"
    )


q_max2 = np.asarray(q_max2)
dis2 = np.asarray(dis2)


# ============================================================
# STRINGER PANEL BUCKLING CHECK
# ============================================================

Kss2 = 5.34 + (4 * dis2**2) / a**2

N2 = (
    Kss2 * D /
    dis2**2
)

buckling_ratio = N2 / q_max2

print("\n========================================")
print("STRINGER PANEL BUCKLING CHECK")
print("========================================")

for i in range(len(N2)):
    print(
        f"Panel {i + 1}: "
        f"Ncr/qs,max = {buckling_ratio[i]:.6f}"
    )


# ============================================================
# PLOT: SHEAR FLOW WITH STRINGERS
# ============================================================

plt.figure(figsize=(12, 8))

plt.plot(
    data2[0, :mid] + cg[0],
    -qs2[:mid],
    linewidth=2,
    label="Upper Surface"
)

plt.plot(
    data2[0, mid:L-1] + cg[0],
    -qs2[mid:L-1],
    linewidth=2,
    label="Lower Surface"
)

plt.grid(True)
plt.xlabel("x [m]")
plt.ylabel(r"$q_s$ [N/m]")
plt.legend(fontsize=font_size)
plt.title(
    "Shear Flow Distribution without spar with stringers",
    fontsize=title_size,
    fontweight="bold"
)

plt.tick_params(labelsize=font_size)
plt.tight_layout()


# ============================================================
# PLOT: AIRFOIL WITH PANELS
# ============================================================

plt.figure(figsize=(12, 8))

for i in range(iu + 1):

    start = iu_ind[i]
    end = iu_ind[i + 1]

    plt.plot(
        data[start:end + 1, 0] + cg[0],
        data[start:end + 1, 1] + cg[1],
        linewidth=2
    )

    centre = int(round(0.5 * (start + end)))

    plt.text(
        data[centre, 0] + cg[0],
        data[centre, 1] + cg[1] + 0.00234,
        f"Panel {i + 1}",
        fontsize=font_size,
        fontweight="bold"
    )


for i in range(il + 1):

    start = il_ind[i] 
    end = il_ind[i + 1] 

    plt.plot(
        data[start:end + 1, 0] + cg[0],
        data[start:end + 1, 1] + cg[1],
        linewidth=2
    )

    centre = int(round(0.5 * (start + end)))

    plt.text(
        data[centre, 0] + cg[0],
        data[centre, 1] + cg[1] + 0.00234,
        f"Panel {iu + 2 + i}",
        fontsize=font_size,
        fontweight="bold"
    )


plt.axis("equal")
plt.grid(True)

plt.xlim(
    np.min(data[:, 0] + cg[0]) - 0.01,
    np.max(data[:, 0] + cg[0]) + 0.01
)

plt.ylim(
    np.min(data[:, 1] + cg[1]) - 0.011,
    np.max(data[:, 1] + cg[1]) + 0.011
)

plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title(
    "Section 3 with Stringers",
    fontsize=title_size,
    fontweight="bold"
)

plt.tick_params(labelsize=font_size)
plt.tight_layout()


# ============================================================
# PLOT: IDEALISED AIRFOIL / BOOMS
# ============================================================

plt.figure(figsize=(12, 8))

plt.plot(
    data[:, 0] + cg[0],
    data[:, 1] + cg[1],
    linewidth=2,
    label="Airfoil"
)

plt.scatter(
    data[:, 0] + cg[0],
    data[:, 1] + cg[1],
    s=50,
    marker="o",
    label="Boom points"
)

plt.axis("equal")
plt.grid(True)

plt.xlim(
    np.min(data[:, 0] + cg[0]) - 0.01,
    np.max(data[:, 0] + cg[0]) + 0.01
)

plt.ylim(
    np.min(data[:, 1] + cg[1]) - 0.011,
    np.max(data[:, 1] + cg[1]) + 0.011
)

plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title(
    "FX63137 Airfoil Idealised",
    fontsize=title_size,
    fontweight="bold"
)

plt.tick_params(labelsize=font_size)
plt.legend()

plt.tight_layout()

plt.savefig(
    "airfoil_booms.png",
    dpi=300,
    bbox_inches="tight"
)


# ============================================================
# DISPLAY ALL FIGURES
# ============================================================

plt.show()
