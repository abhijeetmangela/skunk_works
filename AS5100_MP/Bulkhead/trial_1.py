import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt


def sfd_bmd_fuse(L, p_pos, p_mag, npts=5000):
    
    # Convert inputs to 1D NumPy arrays
    p_pos = np.asarray(p_pos, dtype=float).flatten()
    p_mag = np.asarray(p_mag, dtype=float).flatten()

    # # Check that positions and magnitudes have the same length
    # if len(p_pos) != len(p_mag):
    #     raise ValueError("p_pos and p_mag must have the same length.")

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

p_pos = [2, 5, 8]      # load locations
p_mag = [100, 200, 150]  # downward loads

x_grid, V, M = sfd_bmd_fuse(L, p_pos, p_mag)

