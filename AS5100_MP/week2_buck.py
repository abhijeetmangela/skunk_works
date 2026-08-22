import numpy as np
import matplotlib.pyplot as plt

plt.rc('font', size=16)
np.set_printoptions(precision=8)

from week1_sd import chord, s7055_dl
from week2b_sd import stringer_indices

print("-------------------------------")


print(stringer_indices.shape[0])

print(f"Chord Length = {chord*1e3:.3g} mm")
print(f"No of Stringers = {len(stringer_indices)}")
print(np.array(range(len(stringer_indices))))
# Position from Wing Root
wing_root_tip_len = 850e-3
aileron_inner_edge = 550e-3     # From Wing Root
rib_spanwise_position = np.array([0, 75, 150, 250, 350, 450, 550, 700, 850])*1e-3

def len_from_TE(indices):
    spacing = np.zeros(indices.shape[0] + 1)
    for i in range(len(indices)):
        if i != 0 or i != len(indices) - 1:
            spacing[i] = np.sum(s7055_dl[:indices[i] + 1]) - np.sum(s7055_dl[:indices[i-1] + 1])
        if i == 0:
            spacing[i] = np.sum(s7055_dl[:indices[i] + 1])
        if 

# Spacing from TE through Upper, then Lower Surface back to TE
stringers_spacing = []

# a: Spanwise Length (bw adjacent Ribs) 
# b: Chordwise Length (bw adjacent Stringers) 

# Wing Parameters
a = rib_spanwise_position[1:] - rib_spanwise_position[:-1]
b = 0.265
t = 5e-4

# Material Aluminum
E = 69e9 
nu = 0.33

Ra = 1
Rb = 1
Kss = 5.34 + 4*(b/a)**2

Pcr = Kss*(np.pi**2*E/(12*(1 - nu)))*(t/b)**2*(Ra + (Rb - Ra)*(b/a)**2/2)

