### Author: Mahesh ###

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
plt.rc('font', size=16)
np.set_printoptions(precision=3)

from week1_sd import chord, s7055, skin_thickness


# Choosing the Position and Area Contribution of Stringers
stringers_xloc = np.array([0.15, 0.7])*chord
spar_xloc = 0.3174*chord
booms_xloc = np.concatenate((stringers_xloc, [spar_xloc]))
booms_skin_xloc = np.array([[0, 0.2],
                            [0.5, 1],
                            [0.2, 0.5]])    # Last Row corresponds to Spar
print(f"Booms Chordwise X location: {booms_xloc}")

spar_thickness = 1e-3
flange_length = 40e-3
stringer_thickness = 5e-4
l_stringer_length = 1.5e-2    # Total length of L shaped Stringer
stringer_area = l_stringer_length*stringer_thickness   # Area Contribution of Each Stringer (Identical Stringers are used)

# Splitting the Airfoil into 2 Sections
idx_TE = 0
idx_LE = 42
s7055_u = s7055[:idx_LE + 1, :]
s7055_l = s7055[idx_LE:, :]

# Length of Upper and Lower Airfoil Surfaces
s7055_u_dl = np.linalg.norm(s7055_u[:-1, :] - s7055_u[1:, :], axis=1)
s7055_u_len = np.sum(s7055_u_dl)

s7055_l_dl = np.linalg.norm(s7055_l[:-1, :] - s7055_l[1:, :], axis=1)
s7055_l_len = np.sum(s7055_l_dl)

print(f"Length of Lower Surface of Airfoil = {s7055_l_len*1e3:.3g} mm")
print(f"Length of Upper Surface of Airfoil = {s7055_u_len*1e3:.3g} mm")


plt.scatter(s7055_u[:, 0], s7055_u[:, 1])
plt.scatter(s7055_l[:, 0], s7055_l[:, 1])
plt.scatter(spar_xloc, 0)
plt.axis('equal')


# Cubic Spline Interpolation (Can't be used near the Leading Edge due to High Curvature)
s7055_u_fit = CubicSpline(np.flip(s7055_u[:, 0]), np.flip(s7055_u[:, 1]))
s7055_l_fit = CubicSpline(s7055_l[:, 0], s7055_l[:, 1])
t = np.linspace(0, chord, 100)

# # Cubic Spline Output Plot
# plt.plot(t, s7055_l_fit(t))
# plt.plot(t, s7055_u_fit(t))

# Find the Position of Spar/Stringer Boom on the Wing Skin
def boom_yloc(boom_xloc):
    yloc_u = s7055_u_fit(boom_xloc)
    yloc_l = s7055_l_fit(boom_xloc)
    return np.array([yloc_u, yloc_l])

booms_yloc = boom_yloc(booms_xloc)

# plt.scatter(booms_xloc, booms_yloc[0, :], label="Upper")
# plt.scatter(booms_xloc, booms_yloc[1, :], label="Lower")
plt.legend()
# plt.show()

spar_skin_xloc = booms_skin_xloc[-1]*chord

booms_skin_yloc = boom_yloc(spar_skin_xloc)

for i in [12, 16, 20, 30, 52, 61, 65, 69]:
    plt.scatter(s7055[i, 0], s7055[i, 1], marker='^', c='r')

# plt.scatter(spar_skin_xloc, booms_skin_yloc[0, :], marker='x')
# plt.scatter(spar_skin_xloc, booms_skin_yloc[1, :], marker='x')
plt.show()

# s1: stringer 1 (near LE)
# s2: stringer 2 (near TE)

# Skin Surfaces in LE to TE order
s7055_u_x = np.flip(s7055_u[:, 0])
s7055_u_y = np.flip(s7055_u[:, 1])
s7055_l_x = s7055_l[:, 0]
s7055_l_y = s7055_l[:, 1]

# Boom Boundary x location (chordwise)
# [LE, Spar Boundary (Left), Spar Boundary (Right), TE]
x_bounds = np.array([0, spar_skin_xloc[0], spar_skin_xloc[1], chord])

# y location at the boundaries (from data and spline fit)
y_bounds_u = np.array([s7055_u_y[0], booms_skin_yloc[0, 0], booms_skin_yloc[0, 1], s7055_u_y[-1]])
y_bounds_l = np.array([s7055_l_y[0], booms_skin_yloc[1, 0], booms_skin_yloc[1, 1], s7055_l_y[-1]])


def skin_panel_length(x_arr, y_arr, x0, x1, y0, y1):
    """Arc length of the skin between chordwise locations x0 and x1,
    using actual coordinate points in between and exact boundary points at each end."""
    idx0 = np.searchsorted(x_arr, x0)
    idx1 = np.searchsorted(x_arr, x1)

    x_full = np.concatenate(([x0], x_arr[idx0:idx1], [x1]))
    y_full = np.concatenate(([y0], y_arr[idx0:idx1], [y1]))

    dx = np.diff(x_full)
    dy = np.diff(y_full)
    return np.sum(np.sqrt(dx**2 + dy**2))


# Stringer 1 boom (LE): panels from LE to Spar boundary (Left) 
s1_len_u = skin_panel_length(s7055_u_x, s7055_u_y, x_bounds[0], x_bounds[1], y_bounds_u[0], y_bounds_u[1])
s1_len_l = skin_panel_length(s7055_l_x, s7055_l_y, x_bounds[0], x_bounds[1], y_bounds_l[0], y_bounds_l[1])

# Spar boom: panels between Spar boundary (Left) and Spar boundary (RIght)
spar_len_u = skin_panel_length(s7055_u_x, s7055_u_y, x_bounds[1], x_bounds[2], y_bounds_u[1], y_bounds_u[2])
spar_len_l = skin_panel_length(s7055_l_x, s7055_l_y, x_bounds[1], x_bounds[2], y_bounds_l[1], y_bounds_l[2])

# Stringer 2 boom (TE): panels from Spar boundary (Right) to TE 
s2_len_u = skin_panel_length(s7055_u_x, s7055_u_y, x_bounds[2], x_bounds[3], y_bounds_u[2], y_bounds_u[3])
s2_len_l = skin_panel_length(s7055_l_x, s7055_l_y, x_bounds[2], x_bounds[3], y_bounds_l[2], y_bounds_l[3])

# Convert lengths to lumped skin areas
s1_skin_area_u = s1_len_u*skin_thickness
s1_skin_area_l = s1_len_l*skin_thickness
spar_skin_area_u = spar_len_u*skin_thickness
spar_skin_area_l = spar_len_l*skin_thickness
s2_skin_area_u = s2_len_u*skin_thickness
s2_skin_area_l = s2_len_l*skin_thickness

print(f"Stringer 1 boom - Upper skin length: {s1_len_u*1e3:.3g} mm, area: {s1_skin_area_u:.3g} m^2")
print(f"Stringer 1 boom - Lower skin length: {s1_len_l*1e3:.3g} mm, area: {s1_skin_area_l:.3g} m^2")
print(f"Spar boom - Upper skin length: {spar_len_u*1e3:.3g} mm, area: {spar_skin_area_u:.3g} m^2")
print(f"Spar boom - Lower skin length: {spar_len_l*1e3:.3g} mm, area: {spar_skin_area_l:.3g} m^2")
print(f"Stringer 2 boom - Upper skin length: {s2_len_u*1e3:.3g} mm, area: {s2_skin_area_u:.3g} m^2")
print(f"Stringer 2 boom - Lower skin length: {s2_len_l*1e3:.3g} mm, area: {s2_skin_area_l:.3g} m^2")


booms_skin_area_u = np.array([s1_skin_area_u, spar_skin_area_u, s2_skin_area_u])
booms_skin_area_l = np.array([s1_skin_area_l, spar_skin_area_l, s2_skin_area_l])

span_members_area_u = np.array([l_stringer_length*stringer_thickness, flange_length*spar_thickness, l_stringer_length*stringer_thickness])
span_members_area_l = np.array([l_stringer_length*stringer_thickness, flange_length*spar_thickness, l_stringer_length*stringer_thickness])

booms_total_area_u = booms_skin_area_u + span_members_area_u
booms_total_area_l = booms_skin_area_l + span_members_area_l
print(booms_total_area_u)   # in m^4
print(booms_total_area_l)   # in m^4


# No of Idealised Booms are 6 (Includes 2 Spar Booms)

# x_com = 


# Cell 1


