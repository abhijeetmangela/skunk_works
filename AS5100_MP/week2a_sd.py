import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
plt.rc('font', size=16)
np.set_printoptions(precision=3)

from week1_sd import chord, s7055, skin_thickness, s7055_dl, y_centroid

n_pt = s7055.shape[0] - 1

# Find Aerodynamic Center

# Calculating Center of Pressure (Quantities are in Chordwise Sense)
Xref = 0.25
Cmref = -0.073
CL = 0.72
Xcp = Xref - Cmref/CL
print(f"Chordwise Position of Center of Pressure (from LE) = {Xcp:.2g}")

# # u: Upper, l: Lower, Each Surface includes LE & TE
# idx_TE = 0
# idx_LE = 42
# s7055_u = s7055[:idx_LE + 1, :]
# s7055_l = s7055[idx_LE:, :]

# # Skin Surfaces in LE to TE order
# s7055_u_x = np.flip(s7055_u[:, 0])
# s7055_u_y = np.flip(s7055_u[:, 1])
# s7055_l_x = s7055_l[:, 0]
# s7055_l_y = s7055_l[:, 1]

tD = skin_thickness

# Boom Coordinates (Y coordinates are wrt Centroid)
s7055_x = s7055[:, 0]
s7055_y = s7055[:, 1]
s7055_y_c = s7055[:, 1] - y_centroid
s7055_dl = np.sqrt((np.diff(s7055_x)**2) + (np.diff(s7055_y_c)**2))
s7055_dl = np.append(s7055_dl, np.sqrt((s7055_x[-1] - s7055_x[0])**2 + (s7055_y_c[-1] - s7055_y_c[0])**2))  # Append the last element to make it the same length as s7055_x

# we will iignore the duplicate TE point in s7055 for the purpose of calculating boom areas

# Note that s7055 has a duplicate TE, but B_Area doesn't have duplicate TE

B_Area = np.zeros(n_pt)

for i in range(1, s7055.shape[0] - 1):
    Bi = (tD/6)*(s7055_dl[i - 1]*(2 + s7055_y_c[i - 1]/s7055_y_c[i]) + s7055_dl[i]*(2 + s7055_y_c[i + 1]/s7055_y_c[i]))
    B_Area[i] = Bi

s7055_x  = s7055_x[:-1] # delete duplicate TE point
s7055_y  = s7055_y[:-1] # delete duplicate TE point 
s7055_y_c = s7055_y_c[:-1] # delete duplicate TE point

B_Area[0] = (tD/6)*(s7055_dl[-1]*(2 + s7055_y_c[-2]/s7055_y_c[-1]) + s7055_dl[0]*(2 + s7055_y_c[1]/s7055_y_c[0]))
