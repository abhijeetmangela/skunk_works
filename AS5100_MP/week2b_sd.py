import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size=16)
np.set_printoptions(precision=3)

from week2a_sd import B_Area, s7055_x, s7055_y

stringer_thickness = 0.5e-3
stringer_length = 20e-3     # Total Length of the L shape

Skin_boom_array  = B_Area
Skin_boom_x = s7055_x
Skin_boom_y = s7055_y
num_top_stringers = 4
Stringer_boom_array = np.ones(num_top_stringers)*stringer_thickness
Stringer_boom_x = np.array([s7055_x[i] for i in [12, 16, 20, 30, 52, 61, 65, 69]])
Stringer_boom_y = np.array([s7055_y[i] for i in [12, 16, 20, 30, 52, 61, 65, 69]])

ds_array = np.sqrt((np.diff(Skin_boom_x)**2) + (np.diff(Skin_boom_y)**2))

Sheer_flow_y = 5

# compute the centroid of the cross-section

Centroid_y = (np.sum(Skin_boom_array[:, 0]*Skin_boom_y) + np.sum(Stringer_boom_array[:, 0]*Stringer_boom_y)) / (np.sum(Skin_boom_array[:, 0]) + np.sum(Stringer_boom_array[:, 0]))
Centroid_x = (np.sum(Skin_boom_array[:, 0]*Skin_boom_x) + np.sum(Stringer_boom_array[:, 0]*Stringer_boom_x)) / (np.sum(Skin_boom_array[:, 0]) + np.sum(Stringer_boom_array[:, 0]))

# Area moment of inertia about the centroidal axis

Skin_boom_y_centroid = Skin_boom_y - Centroid_y
Skin_boom_x_centroid = Skin_boom_x - Centroid_x

Stringer_boom_y_centroid = Stringer_boom_y - Centroid_y
Stringer_boom_x_centroid = Stringer_boom_x - Centroid_x

Centroidal_Ixx = np.sum(Skin_boom_array[:, 0]*Skin_boom_y_centroid**2) + np.sum(Stringer_boom_array[:, 0]*Stringer_boom_y**2)
Centroidal_Iyy = np.sum(Skin_boom_array[:, 0]*Skin_boom_x_centroid**2) + np.sum(Stringer_boom_array[:, 0]*Stringer_boom_x**2)
Centroidal_Ixy = np.sum(Skin_boom_array[:, 0]*Skin_boom_x_centroid*Skin_boom_y_centroid) + np.sum(Stringer_boom_array[:, 0]*Stringer_boom_x*Stringer_boom_y)

DR = Centroidal_Ixx*Centroidal_Iyy - Centroidal_Ixy**2

q_b_0 = 0

q_b_array = np.array([q_b_0])

for i in range(len(Skin_boom_array)):
    for j in range(len(Stringer_boom_array)):
        if j < num_top_stringers:
            while Skin_boom_x_centroid[i] > Stringer_boom_x_centroid[j]:
                q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i, 0]*(Centroidal_Ixy*Skin_boom_x_centroid[i] - Centroidal_Ixx*Skin_boom_y_centroid[i] )/DR)

            q_b_0 = q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j, 0]*(Centroidal_Ixy*Stringer_boom_x_centroid[j] - Centroidal_Ixx*Stringer_boom_y_centroid[j] )/DR
        else:
            while Skin_boom_x_centroid[i] < Stringer_boom_x_centroid[j]:
                q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i, 0]*(Centroidal_Ixy*Skin_boom_x_centroid[i] - Centroidal_Ixx*Skin_boom_y_centroid[i] )/DR)

            q_b_0 = q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j, 0]*(Centroidal_Ixy*Stringer_boom_x_centroid[j] - Centroidal_Ixx*Stringer_boom_y_centroid[j] )/DR


# Wing_Ixx = np.sum(Skin_boom_array*Skin_boom_y**2) + np.sum(Stringer_boom_array*Stringer_boom_y**2)
# Wing_Iyy = np.sum(Skin_boom_array*Skin_boom_x**2) + np.sum(Stringer_boom_array*Stringer_boom_x**2)
# Wing_Ixy = np.sum(Skin_boom_array*Skin_boom_x*Skin_boom_y) + np.sum(Stringer_boom_array*Stringer_boom_x*Stringer_boom_y)

# DR = Wing_Ixx*Wing_Iyy - Wing_Ixy**2

# q_b_0 = 0

# q_b_array = np.array([q_b_0])

# for i in range(len(Skin_boom_array)):
#     for j in range(len(Stringer_boom_array)):
#         if j < num_top_stringers:
#             while Skin_boom_x[i] > Stringer_boom_x[j]:
#                 q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i, 0]*(Wing_Ixy*Skin_boom_x[i] - Wing_Ixx*Skin_boom_y[i] )/DR)

#             q_b_0 = q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j]*(Wing_Ixy*Stringer_boom_x[j] - Wing_Ixx*Stringer_boom_y[j] )/DR
#         else:
#             while Skin_boom_x[i] < Stringer_boom_x[j]:
#                 q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i, 0]*(Wing_Ixy*Skin_boom_x[i] - Wing_Ixx*Skin_boom_y[i] )/DR)

#             q_b_0 = q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j]*(Wing_Ixy*Stringer_boom_x[j] - Wing_Ixx*Stringer_boom_y[j] )/DR

q_s_0 = np.sum(q_b_array*ds_array)/np.sum(ds_array)


            
    