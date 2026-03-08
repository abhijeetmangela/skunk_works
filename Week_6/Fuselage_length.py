## Fuselage Sizing

## AUTHOR = ABHIJEET MANGELA (AE25M034)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Data = pd.read_excel("Similar UAV.xlsx")

Weights = Data["MTOW (kg)"].to_numpy()
Fuselage_lengths = Data["Fuselage Length (m)"].to_numpy()

# Power law fit
a, b = np.polyfit(np.log(Weights), np.log(Fuselage_lengths), 1)

# Smooth curve
x_array = np.linspace(min(Weights)*0.8, max(Weights)*1.2, 200)
y_array = np.exp(b) * x_array**a

# Our UAV
Our_weight = 7.021
Our_length = np.exp(b) * Our_weight**a

plt.figure(figsize=(7,5))

# Data points
plt.scatter(Weights, Fuselage_lengths,
            s=60, alpha=0.8, label="Existing UAVs")

# Fit curve
plt.plot(x_array, y_array,
         linewidth=2.5,
         label=f"Fit: L = {np.exp(b):.2f} W^{a:.2f}")

# Our UAV point
plt.scatter(Our_weight, Our_length,
            color='red', s=120, marker='*',
            label="Our UAV")

# Vertical line at our weight
plt.axvline(Our_weight, linestyle="--", color="red", alpha=0.7)

# Horizontal line at predicted length
plt.axhline(Our_length, linestyle="--", color="gray", alpha=0.5)

# Annotate value
plt.text(Our_weight*1.02, Our_length,
         f"L = {Our_length:.2f} m",
         verticalalignment='bottom')

plt.xlabel("Maximum Takeoff Weight (kg)")
plt.ylabel("Fuselage Length (m)")
plt.title("Power Law Relationship Between UAV MTOW and Fuselage Length")

plt.grid(True, linestyle="--", alpha=0.4)
plt.legend()
plt.tight_layout()

plt.show()