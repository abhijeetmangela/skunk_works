import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt


def sfd_bmd_fuse(L, p_pos, p_mag, npts=5000):
    
    # Convert inputs to 1D NumPy arrays
    p_pos = np.asarray(p_pos, dtype=float).flatten()
    p_mag = np.asarray(p_mag, dtype=float).flatten()


    # Discretized x-coordinate
    x_grid = np.linspace(0, L, npts)

    # Shear force
    # True wherever load position <= x
    ind = p_pos[:, np.newaxis] <= x_grid[np.newaxis, :]

    # Sum all loads that have been encountered
    V = -np.sum(p_mag[:, np.newaxis] * ind, axis=0)


    # Bending moment


    dx = x_grid[np.newaxis, :] - p_pos[:, np.newaxis]
    dx[dx < 0] = 0

    M = -np.sum(p_mag[:, np.newaxis] * dx, axis=0)

    # Plot SFD and BMD


    fig, axes = plt.subplots(2, 1, figsize=(10, 7))

    # Shear Force Diagram
    axes[0].plot(x_grid, V, linewidth=1.8)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("Shear V(x)")
    axes[0].set_title("Shear Force Diagram (No Supports)")
    axes[0].grid(True)

    # Bending Moment Diagram
    axes[1].plot(x_grid, M, linewidth=1.8)
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("Moment M(x)")
    axes[1].set_title("Bending Moment Diagram (No Supports)")
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()

    return x_grid, V, M

L = 10

p_pos = [2, 5, 6, 8]      # load locations in meters
p_mag = [100, 200,-450, 150]  # downward loads in Newtons

x, V, M = sfd_bmd_fuse(L, p_pos, p_mag)

V_max = np.max(np.abs(V))       # N
M_max = np.max(np.abs(M))       # N.m

h_f = 0.18 # maximum height. Change this 
t = 0.0005 # thickness of the plate

A = 1 # no longerons for now

I_skin = (h_f**4 / 12) - ((h_f - t)**4 / 12)
I_boom = 4 * A * h_f * h_f * 0.25

sig_yield = 95000000
FOS = 8.5 

sig_all = sig_yield / FOS

I_req = M_max * h_f / (2 * sig_all)

if I_req > I_skin + I_boom:
    print("Longeron is not enough to withstand bending moment of fuselage")
else:
    print("Longeron is enough to withstand bending moment")

# Shear Flow & F_cr Calculation of Fuselage

I_xx = I_skin + I_boom
# Maximum shear flow
q_max = (3 * V_max * t * h_f**2) / (8 * I_xx) # at the midpoint of the side webbs 

# Maximum shear stress
tau_max = q_max / t # Pa

a = 0.1 # m, distance between bulkheads, check this value
b = 0.36 # m, distance between longerons, 2* height of fuselage, check this value

E = 70e9
nu = 0.33

Ra = 1.62 # fixed support


Kss = 5.34 + 4 * (b / a)**2

F_cr = (
    Kss * np.pi**2 * E
    / (12 * (1 - nu**2)) 
) * (t / b)**2 * Ra # Pa

F_cr = F_cr / FOS 

# Buckling Check

if tau_max > F_cr:
    print("Buckled")
else:
    print("Not Buckled")

print("\n------ Results ------")
print(f"V_max   = {V_max:.4f} N")
print(f"M_max   = {M_max:.4f} N.m")
print(f"I_skin  = {I_skin:.6e} m^4")
print(f"I_boom  = {I_boom:.6e} m^4")
print(f"I_req   = {I_req:.6e} m^4")
print(f"I_xx    = {I_xx:.6e} m^4")
print(f"q_max   = {q_max:.6f} N/m")
print(f"tau_max = {tau_max:.6e} Pa")
print(f"Kss     = {Kss:.6f}")
print(f"F_cr    = {F_cr:.6e} Pa")

