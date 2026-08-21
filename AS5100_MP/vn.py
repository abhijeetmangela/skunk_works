# -*- coding: utf-8 -*-
"""UAV_Drag_Calculations

AUTHOR = G. Vishwas (AE25M024), Divyam Pandey (AE25M021)

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#AUTHOR: DIVYAM PANDEY(AE25M021)
class UAVPerformanceCalculator:
    def __init__(self):
        # ==========================================
        # 1. ENVIRONMENT & AIRCRAFT PARAMETERS
        # ==========================================
        self.rho = 1.220  # Air density at cruise altitude (kg/m^3)
        self.mu = 1.626e-5  # Dynamic viscosity (Ns/m^2)
        self.V_cruise = 18.8  # Cruise velocity (m/s)
        self.Mach = self.V_cruise / 349.02  # Cruise Mach number

        self.W_o = 6.592 * 9.81  # Aircraft total weight in Newtons (~10kg)
        self.S_ref = 0.5307  # Wing reference area (m^2)
        self.AR = 7.53  # Aspect Ratio
        self.e = 0.7648  # Oswald Efficiency factor

        # Battery & Propulsion (from report)
        self.P_av = 490.6  # Max Power Available (W)
        self.Energy_Bat = 1.75e6  # Battery Energy (Joules, 1.75 MJ)
        self.eta_p = 0.7  # Propeller efficiency
        self.eta_m = 0.7  # Motor efficiency
        self.eta_b = 0.95  # Battery efficiency

        # Structural Limits
        self.n_max = 3.0
        self.n_min = -1.5
        self.CL_max = 1.4
        self.CL_min = -0.4961

        # ==========================================
        # 2. COMPONENT GEOMETRY (FROM CAD)
        # ==========================================
        # Dictionary format: [characteristic_length (m), S_wet (m^2), t/c, (x/c)_m, sweep (rad), Q_factor]
        self.wing = {
            "l": 0.2653,
            "Swet": 0.93,
            "tc": 0.105,
            "xcm": 0.318,
            "sweep": 0.0,
            "Q": 1.0,
        }
        self.htail = {
            "l": 0.154,
            "Swet": 0.22,
            "tc": 0.12,
            "xcm": 0.3,
            "sweep": 0.0,
            "Q": 1.0,
        }
        self.vtail = {
            "l": 0.178,
            "Swet": 0.09,
            "tc": 0.12,
            "xcm": 0.3,
            "sweep": 0.0,
            "Q": 1.05,
        }

        # Fuselage: [length (m), S_wet (m^2), A_max (m^2), Q_factor, square_sided_penalty]
        self.fuselage = {
            "l": 0.8,
            "Swet": 0.75,
            "tc": 0.25,
            "xcm": 0.3,
            "sweep": 0.0,
            "Q": 1.0,
            # "Amax": 0.06,
            "sq_penalty": 0.2,
        }

    def calc_skin_friction(self, length):
        """Calculates skin friction coefficient Cf based on Reynolds number."""
        Re = (self.rho * self.V_cruise * length) / self.mu
        return 0.455 / ((np.log10(Re)) ** 2.58)

    def calc_form_factor_lifting(self, comp):
        """Calculates Form Factor for Wing and Tails."""
        term1 = 1 + (0.6 / comp["xcm"]) * comp["tc"] + 100 * (comp["tc"] ** 4)
        term2 = 1 #1.34 * (self.Mach**0.18) * (np.cos(comp["sweep"]) ** 0.28)
        return term1 * term2

    def calc_form_factor_fuselage(self, comp):

        """Calculates Form Factor for Wing and Tails."""
        term1 = 1 + (0.6 / comp["xcm"]) * comp["tc"] + 100 * (comp["tc"] ** 4)
        term2 = 1 #1.34 * (self.Mach**0.18) * (np.cos(comp["sweep"]) ** 0.28)
        # return term1 * term2

        # """Calculates Form Factor for Fuselage."""
        # d = np.sqrt((4 / np.pi) * self.fuselage["Amax"])
        # f = self.fuselage["l"] / d
        # FF = 0.9 + (5 / (f**1.5)) + (f / 400)
        # Apply square-sided penalty if applicable
        FF = term1*term2
        penalty = self.fuselage["sq_penalty"]
        FF = FF + (FF - 1) * penalty
        return FF

    def execute_parasitic_drag_table(self):
        """Calculates CD0 and outputs the component breakdown table."""
        components = ["Wing", "Horizontal Tail", "Vertical Tail"]
        data = [self.wing, self.htail, self.vtail]

        results = []
        total_cd0_area = 0

        # Lifting surfaces
        for name, comp in zip(components, data):
            Cf = self.calc_skin_friction(comp["l"])
            FF = self.calc_form_factor_lifting(comp)
            area_val = Cf * FF * comp["Q"] * comp["Swet"]
            total_cd0_area += area_val
            results.append(
                [
                    name,
                    round(Cf, 5),
                    round(FF, 3),
                    comp["Q"],
                    comp["Swet"],
                    round(area_val, 5),
                ]
            )

        # Fuselage
        Cf_fuse = self.calc_skin_friction(self.fuselage["l"])
        FF_fuse = self.calc_form_factor_fuselage(self.fuselage)
        area_val_fuse = Cf_fuse * FF_fuse * self.fuselage["Q"] * self.fuselage["Swet"]
        total_cd0_area += area_val_fuse
        results.append(
            [
                "Fuselage",
                round(Cf_fuse, 5),
                round(FF_fuse, 3),
                self.fuselage["Q"],
                self.fuselage["Swet"],
                round(area_val_fuse, 5),
            ]
        )

        self.CD0 = total_cd0_area / self.S_ref

        print("\n--- COMPONENT PARASITIC DRAG BREAKDOWN ---")
        df = pd.DataFrame(
            results,
            columns=[
                "Component",
                "Cf",
                "Form Factor (FF)",
                "Interference (Q)",
                "S_wet (m^2)",
                "Cf*FF*Q*Swet",
            ],
        )
        print(df.to_markdown(index=False))
        print(f"\nTotal Aircraft CD0 = {self.CD0:.5f}")

    def calculate_performance(self):
        """Calculates Range, Endurance, and max L/D."""
        self.K = 1 / (np.pi * self.AR * self.e)
        self.LD_max = 0.5 * np.sqrt(1 / (self.K * self.CD0))

        # Max Range Conditions (CD0 = CDi)
        V_Rmax = np.sqrt((2 * self.W_o) / (self.rho * self.S_ref)) * (
            (self.K / self.CD0) ** 0.25
        )
        P_Rmax = self.rho * self.S_ref * self.CD0 * (V_Rmax**3)
        R_max = (
            self.eta_p * self.eta_m * self.eta_b * self.Energy_Bat * (V_Rmax / P_Rmax)
        )

        # Max Endurance Conditions (3*CD0 = CDi)
        V_Emax = np.sqrt((2 * self.W_o) / (self.rho * self.S_ref)) * (
            (self.K / (3 * self.CD0)) ** 0.25
        )
        P_Emax = 2 * self.rho * self.S_ref * self.CD0 * (V_Emax**3)
        E_max = (self.eta_p * self.eta_m * self.eta_b * self.Energy_Bat) / P_Emax

        print("\n--- PERFORMANCE METRICS ---")
        print(f"(L/D)max: {self.LD_max:.2f}")
        print(f"Max Range Velocity: {V_Rmax:.2f} m/s")
        print(f"Max Range: {R_max/1000:.2f} km")
        print(f"Max Endurance Velocity: {V_Emax:.2f} m/s")
        print(f"Max Endurance: {E_max/3600:.2f} hrs")

    def plot_graphs(self):
        """Generates all the required plots."""
        v_arr = np.linspace(1, 35, 100)

        # 1. Power vs Velocity
        P_req = 0.5 * self.rho * (v_arr**3) * self.S_ref * self.CD0 + (
            2 * self.K * (self.W_o**2)
        ) / (self.rho * v_arr * self.S_ref)

        plt.figure(figsize=(10, 6))
        plt.plot(v_arr, P_req, "b-", linewidth=2, label="Power Required")
        plt.axhline(self.P_av, color="k", linestyle="--", label="Power Available")
        plt.title("Power vs Velocity")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Power Required (W)")
        plt.ylim(0, 1000)
        plt.grid(True)
        plt.legend()
        plt.show()

        # 2. ROC vs Velocity
        ROC = (self.P_av - P_req) / self.W_o

        plt.figure(figsize=(10, 6))
        plt.plot(v_arr, ROC, "b-", linewidth=2)
        plt.title("ROC vs Velocity")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("ROC (m/s)")
        plt.ylim(0, 5)
        plt.grid(True)
        plt.show()

        # 3. Drag Polar
        CL_arr = np.linspace(0, 1.4, 100)
        CD_arr = self.CD0 + self.K * (CL_arr**2)

        plt.figure(figsize=(10, 6))
        plt.plot(CD_arr, CL_arr, "b-", linewidth=2)
        plt.title("Drag Polar")
        plt.xlabel("$C_D$")
        plt.ylabel("$C_L$")
        plt.xlim(0.02, 0.12)
        plt.ylim(0, 1.6)
        plt.grid(True)
        plt.show()




# Run the program
if __name__ == "__main__":
    uav = UAVPerformanceCalculator()
    uav.execute_parasitic_drag_table()
    uav.calculate_performance()
    uav.plot_graphs()