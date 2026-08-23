# Code for wing without spar
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
num_stringers = 4

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

# upper_stringer_x_norm = np.array([0.16, 0.5, 0.65, 0.78])
# lower_stringer_x_norm = np.array([0.16, 0.5, 0.65, 0.78])

upper_stringer_x_norm = np.array([0.16, 0.65])
lower_stringer_x_norm = np.array([0.16, 0.65])

# Vertical shear force used only to obtain the shear-center moment.
# The resulting shear-center location is independent of its magnitude.
V = 13.3 # N

x_c_cp = 0.3174 # chord ratio not in m 


Xref = 0.3174
Cmref = -0.073
CL = 0.72

Xcp = Xref - Cmref / CL
x_cp = Xcp * chord

print()
print(f"Xcp = {Xcp:.6f} c")
print(f"x_cp = {x_cp:.6e} m")

# ---------------------------------------------------------------------
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

print(ds)

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

# ---------------------------------------------------------------------
# 2. BOOM AREAS
#
# Use the standard idealisation:
#
# B_i = t/6 * [
#       l_(i-1) (2 + y_(i-1)/y_i)
#     + l_i     (2 + y_(i+1)/y_i)
#     ]
#
# Here the boom formula is evaluated using the skin centroidal
# y-coordinate. The final centroid is then calculated from the
# complete idealised section.
#
# Points extremely close to the reference axis make this formula
# singular. For a robust shear-center calculation, use panel-area
# discretisation instead when such points occur.
# ---------------------------------------------------------------------

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

print(upper_stringer_x_norm)
print(lower_stringer_x_norm)

print(upper_stringer_indices)
print(lower_stringer_indices)

# ---------------------------------------------------------------------
# 4. FINAL CENTROID
# ---------------------------------------------------------------------

A_total = np.sum(B_total)

x_bar = np.sum(B_total * x) / A_total
y_bar = np.sum(B_total * y) / A_total

x_c = x - x_bar
y_c = y - y_bar

print()
print(f"Total idealised area = {A_total:.6e} m^2")
print(f"Final centroid x     = {x_bar:.6e} m")
print(f"Final centroid y     = {y_bar:.6e} m")

# ---------------------------------------------------------------------
# 5. SECOND MOMENTS OF AREA
# ---------------------------------------------------------------------

Ixx = np.sum(B_total * y_c**2)
Iyy = np.sum(B_total * x_c**2)
Ixy = np.sum(B_total * x_c * y_c)

D = Ixx * Iyy - Ixy**2

print()
print(f"Ixx = {Ixx:.6e} m^4")
print(f"Iyy = {Iyy:.6e} m^4")
print(f"Ixy = {Ixy:.6e} m^4")
print(f"D   = {D:.6e} m^8")


# ---------------------------------------------------------------------
# 6. BASIC SHEAR FLOW
#
# For a vertical shear force V:
#
# q_b = V/D * (Ixy*x - Ixx*y) integrated as a cumulative
# boom contribution.
#
# The first boom is taken as the cut/reference point.
#
# The panel flow is reconstructed by averaging adjacent nodal flows.
# ---------------------------------------------------------------------

# ================================================================
# 6. BASIC SHEAR FLOW
#    CUT = LAST PANEL (n-1 -> 0)
# ================================================================

te_cut_panel = n - 1

q_b_panel = np.zeros(n)

q_running = 0.0

for i in range(n - 1):

    q_b_panel[i] = q_running

    j = i + 1

    dq = (
        V / D
        * B_total[j]
        * (
            Ixy * x_c[j]
            - Ixx * y_c[j]
        )
    )

    q_running += dq

q_b_panel[te_cut_panel] = 0.0

# ---------------------------------------------------------------------
# 7. REDUNDANT CONSTANT SHEAR FLOW
#
# For a single-cell closed section:
#
# integral(q ds) = 0
#
# Therefore:
#
# q0 = - integral(q_b ds) / integral(ds)
# ---------------------------------------------------------------------



q0 = -np.sum(q_b_panel * ds) / np.sum(ds)

q_panel = q_b_panel + q0

print()
print(f"Basic-flow closure correction q0 = {q0:.6e} N/m")


# ---------------------------------------------------------------------
# 8. TORQUE ABOUT THE FINAL CENTROID
#
# For each straight panel:
#
# dT = q * (x dy - y dx)
#
# For a straight panel:
#
# integral(x dy - y dx)
#     = x_i*y_(i+1) - y_i*x_(i+1)
# ---------------------------------------------------------------------

# ================================================================
# TORQUE ABOUT CENTROID
# ================================================================

x_next = np.roll(x_c, -1)
y_next = np.roll(y_c, -1)

panel_moment_arm = (
    x_c * y_next
    - y_c * x_next
)

T_shear_flow = np.sum(
    q_panel * panel_moment_arm
)

# ================================================================
# SHEAR CENTER
# ================================================================

e = -T_shear_flow / V

x_sc_LE = x_bar + e

print(f"Torque about centroid = {T_shear_flow:.6e} N m")
print(f"Shear-center offset = {e:.6e} m")
print(f"Shear center from LE = {x_sc_LE:.6e} m")
print(f"Shear center x/c = {x_sc_LE/chord:.6f}")


# ================================================================
# TORQUE ABOUT SHEAR CENTER
# ================================================================

T_CP = V * (x_cp - x_sc_LE)

print(f"T_CP = {T_CP:.6e} N m")


# ================================================================
# MEDIAN-LINE AREA
# ================================================================

x_mid = 0.5 * (x + np.roll(x, -1))
y_mid = 0.5 * (y + np.roll(y, -1))

x_mid_next = np.roll(x_mid, -1)
y_mid_next = np.roll(y_mid, -1)

A_m = 0.5 * abs(
    np.sum(
        x_mid * y_mid_next
        - y_mid * x_mid_next
    )
)

print(f"A_m = {A_m:.6e} m^2")


# ================================================================
# TORQUE-INDUCED SHEAR FLOW
# ================================================================

q_CP = T_CP / (2.0 * A_m)

print(f"q_CP = {q_CP:.6e} N/m")


# ================================================================
# FINAL SHEAR FLOW
# ================================================================

q_final = q_panel + q_CP

print()
print(f"Minimum q_final = {np.min(q_final):.6e} N/m")
print(f"Maximum q_final = {np.max(q_final):.6e} N/m")
print(f"Maximum |q_final| = {np.max(np.abs(q_final)):.6e} N/m")
plt.plot(x_mid, q_final)
plt.xlabel('x (m)')
plt.ylabel('Shear flow q (N/m)')
plt.title('Closed-Section Shear Flow')
plt.grid()
plt.show()
# Author: Mahesh

# Check Critical Buckling Loads
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
# stringers_idx = np.insert(stringers_idx, [4], [42])
stringers_idx = np.insert(stringers_idx, [2], [42])
print(np.array(range(1, num_stringers)))
print(stringers_idx)

# num_stringer_spar = num_stringers + 2
stringers_spacing = np.zeros(num_stringers + 2)

for i in range(1, num_stringers + 1):
    stringers_spacing[i] = np.sum(ds[:stringers_idx[i] + 1]) - np.sum(ds[:stringers_idx[i-1] + 1])

stringers_spacing[0] = np.sum(ds[:stringers_idx[0] + 1])
stringers_spacing[num_stringers + 1] = np.sum(ds) - np.sum(ds[:stringers_idx[-1] + 1])
# print(stringers_idx)
print(stringers_spacing)

# print(np.sum(ds) - np.sum(stringers_spacing))
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
qmax_panel = np.zeros(len(stringers_idx) + 1)

print(np.array(range(1, len(qmax_panel) - 1)))

for i in range(1, len(qmax_panel) - 1):

    panel_start_idx = stringers_idx[i - 1]
    panel_end_idx = stringers_idx[i]
    # print(i, panel_start_idx, panel_end_idx, "length =", panel_end_idx - panel_start_idx)

    qmax_panel[i] = np.max(np.abs(q_final[panel_start_idx:panel_end_idx]))

qmax_panel[0] = np.max(np.abs(q_final[:stringers_idx[0]]))
qmax_panel[-1] = np.max(np.abs(q_final[stringers_idx[-1]:]))

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
print(len(a))   # Ribs Panel Spacing
print(len(b))   # Stringer Panel Spacing