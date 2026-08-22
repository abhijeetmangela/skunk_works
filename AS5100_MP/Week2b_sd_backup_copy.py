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
num_stringers = 8
num_top_stringers = 4
stringer_indices = np.array([12, 16, 20, 30, 52, 61, 65, 69])
Stringer_boom_array = np.ones(num_stringers)*stringer_thickness*stringer_length
Stringer_boom_x = np.array([s7055_x[i] for i in stringer_indices])
Stringer_boom_y = np.array([s7055_y[i] for i in stringer_indices])

Total_boom_array = Skin_boom_array.copy()

for j, idx in enumerate(stringer_indices):
    Total_boom_array[idx] += Stringer_boom_array[j]



ds_array = np.sqrt((np.diff(Skin_boom_x)**2) + (np.diff(Skin_boom_y)**2))
ds_array = np.append(ds_array, np.sqrt((Skin_boom_x[-1] - Skin_boom_x[0])**2 + (Skin_boom_y[-1] - Skin_boom_y[0])**2))  # Append the last element to make it the same length as Skin_boom_array

Shear_flow_y = 5

# compute the centroid of the cross-section

Centroid_y = (np.sum(Skin_boom_array*Skin_boom_y) + np.sum(Stringer_boom_array*Stringer_boom_y)) / (np.sum(Skin_boom_array) + np.sum(Stringer_boom_array))
Centroid_x = (np.sum(Skin_boom_array*Skin_boom_x) + np.sum(Stringer_boom_array*Stringer_boom_x)) / (np.sum(Skin_boom_array) + np.sum(Stringer_boom_array))

print(f"Centroid of the cross-section: ({Centroid_x:.3f}, {Centroid_y:.3f})")

# Area moment of inertia about the centroidal axis

Skin_boom_y_centroid = Skin_boom_y - Centroid_y
Skin_boom_x_centroid = Skin_boom_x - Centroid_x

Stringer_boom_y_centroid = Stringer_boom_y - Centroid_y
Stringer_boom_x_centroid = Stringer_boom_x - Centroid_x

Centroidal_Ixx = np.sum(Skin_boom_array*Skin_boom_y_centroid**2) + np.sum(Stringer_boom_array*Stringer_boom_y_centroid**2)
Centroidal_Iyy = np.sum(Skin_boom_array*Skin_boom_x_centroid**2) + np.sum(Stringer_boom_array*Stringer_boom_x_centroid**2)
Centroidal_Ixy = np.sum(Skin_boom_array*Skin_boom_x_centroid*Skin_boom_y_centroid) + np.sum(Stringer_boom_array*Stringer_boom_x_centroid*Stringer_boom_y_centroid)

DR = Centroidal_Ixx*Centroidal_Iyy - Centroidal_Ixy**2

print(f"Centroidal Ixx: {Centroidal_Ixx:.3f}, Centroidal Iyy: {Centroidal_Iyy:.3f}, Centroidal Ixy: {Centroidal_Ixy:.3f}, DR: {DR:.3f}")

q_b_0 = 0

q_b_array = np.array([q_b_0])
print("initialised q_b_array: ", q_b_array)

# for i in range(len(Skin_boom_array)):
#     for j in range(len(Stringer_boom_array)):
#         if j < num_top_stringers:
#             while Skin_boom_x_centroid[i] > Stringer_boom_x_centroid[j]:
#                 q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i]*(Centroidal_Ixy*Skin_boom_x_centroid[i] - Centroidal_Ixx*Skin_boom_y_centroid[i] )/DR)

#             q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j]*(Centroidal_Ixy*Stringer_boom_x_centroid[j] - Centroidal_Ixx*Stringer_boom_y_centroid[j] )/DR)
#         else:
#             while Skin_boom_x_centroid[i] < Stringer_boom_x_centroid[j]:
#                 q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Skin_boom_array[i]*(Centroidal_Ixy*Skin_boom_x_centroid[i] - Centroidal_Ixx*Skin_boom_y_centroid[i] )/DR)

#             q_b_array = np.append(q_b_array, q_b_array[-1] + Sheer_flow_y*Stringer_boom_array[j]*(Centroidal_Ixy*Stringer_boom_x_centroid[j] - Centroidal_Ixx*Stringer_boom_y_centroid[j] )/DR)

q_b_array = np.zeros(len(Total_boom_array))

for i in range(1, len(Total_boom_array)):

    q_increment = (
        Shear_flow_y
        * Total_boom_array[i]
        * (
            Centroidal_Ixy * Skin_boom_x_centroid[i]
            - Centroidal_Ixx * Skin_boom_y_centroid[i]
        )
        / DR
    )

    q_b_array[i] = q_b_array[i-1] + q_increment

print("Basic shear flow:")
print(q_b_array)

# print(f"Shear Flow in the booms: {q_b_array}")
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

q_s_0 = np.sum(q_b_array*ds_array)/np.sum(ds_array) # formula for q_s_0 is derived from the fact that the integral of q_s over the length of the skin must equal the integral of q_b over the length of the skin. This is because the total shear flow in the skin must be equal to the total shear flow in the booms, as they are connected and must balance each other out. Therefore, we can calculate q_s_0 by taking the average of q_b over the length of the skin, which is given by the formula above.



            
    