import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ============================================================
# DEFAULT VALUES
# Taken from the original VN_test.py
# ============================================================

DEFAULTS = {
    "mass": 6.592,          # kg
    "S": 0.5307,            # m^2
    "rho": 1.220,           # kg/m^3

    "CL_max_pos": 1.3,
    "CL_max_neg": -0.8,

    "n_max_pos": 3.5,
    "n_max_neg": -1.5,

    "V_max": 24.0,          # m/s
    "dive_factor": 1.25
}

G = 9.81


# ============================================================
# MAIN GUI CLASS
# ============================================================

class VNDiagramGUI:

        # ========================================================
    # EXPORT PLOT
    # ========================================================

    def export_plot(self):

        # ----------------------------------------------------
        # File type selection
        # ----------------------------------------------------

        file_path = filedialog.asksaveasfilename(
            title="Export V-n Diagram",
            defaultextension=".png",
            filetypes=[
                ("PNG Image", "*.png"),
                ("JPEG Image", "*.jpg"),
                ("TIFF Image", "*.tiff"),
                ("PDF Document", "*.pdf"),
                ("SVG Vector Image", "*.svg"),
                ("All Files", "*.*")
            ]
        )

        if not file_path:
            return

        # ----------------------------------------------------
        # Create resolution dialog
        # ----------------------------------------------------

        resolution_window = tk.Toplevel(self.root)

        resolution_window.title("Export Settings")
        resolution_window.geometry("350x230")
        resolution_window.resizable(False, False)

        ttk.Label(
            resolution_window,
            text="Export Settings",
            font=("Arial", 14, "bold")
        ).pack(pady=15)

        # ----------------------------------------------------
        # DPI
        # ----------------------------------------------------

        dpi_frame = ttk.Frame(resolution_window)

        dpi_frame.pack(
            pady=5
        )

        ttk.Label(
            dpi_frame,
            text="Resolution (DPI):"
        ).grid(
            row=0,
            column=0,
            padx=5
        )

        dpi_entry = ttk.Entry(
            dpi_frame,
            width=10
        )

        dpi_entry.insert(
            0,
            "300"
        )

        dpi_entry.grid(
            row=0,
            column=1,
            padx=5
        )

        # ----------------------------------------------------
        # Transparent background
        # ----------------------------------------------------

        transparent_var = tk.BooleanVar(
            value=False
        )

        transparent_check = ttk.Checkbutton(
            resolution_window,
            text="Transparent background",
            variable=transparent_var
        )

        transparent_check.pack(
            pady=10
        )

        # ----------------------------------------------------
        # Export function
        # ----------------------------------------------------

        def save_file():

            try:

                dpi = int(
                    dpi_entry.get()
                )

                if dpi <= 0:

                    raise ValueError

            except ValueError:

                messagebox.showerror(
                    "Invalid DPI",
                    "Please enter a positive integer for DPI."
                )

                return

            try:

                self.figure.savefig(
                    file_path,
                    dpi=dpi,
                    bbox_inches="tight",
                    transparent=transparent_var.get()
                )

                resolution_window.destroy()

                messagebox.showinfo(
                    "Export Successful",
                    f"Plot exported successfully:\n\n"
                    f"{file_path}\n\n"
                    f"Resolution: {dpi} DPI"
                )

            except Exception as e:

                messagebox.showerror(
                    "Export Error",
                    f"Could not export the plot:\n\n{e}"
                )

        # ----------------------------------------------------
        # Buttons
        # ----------------------------------------------------

        button_frame = ttk.Frame(
            resolution_window
        )

        button_frame.pack(
            pady=15
        )

        ttk.Button(
            button_frame,
            text="Export",
            command=save_file
        ).grid(
            row=0,
            column=0,
            padx=5
        )

        ttk.Button(
            button_frame,
            text="Cancel",
            command=resolution_window.destroy
        ).grid(
            row=0,
            column=1,
            padx=5
        )

    def __init__(self, root):

        self.root = root

        self.root.title("UAV V-n Diagram")
        self.root.geometry("1400x800")
        self.root.minsize(1100, 700)

        # ----------------------------------------------------
        # Main layout
        # ----------------------------------------------------

        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        # ====================================================
        # LEFT CONTROL PANEL
        # ====================================================

        control_frame = ttk.Frame(
            root,
            padding=15
        )

        control_frame.grid(
            row=0,
            column=0,
            sticky="ns"
        )

        # ====================================================
        # TITLE
        # ====================================================

        title = ttk.Label(
            control_frame,
            text="V-n Diagram Parameters",
            font=("Arial", 16, "bold")
        )

        title.grid(
            row=0,
            column=0,
            columnspan=2,
            pady=(0, 20)
        )

        # ====================================================
        # INPUT VARIABLES
        # ====================================================

        self.entries = {}

        parameters = [
            ("Mass", "mass", "kg"),
            ("Wing Area", "S", "m²"),
            ("Air Density", "rho", "kg/m³"),
            ("CL max (+)", "CL_max_pos", ""),
            ("CL max (-)", "CL_max_neg", ""),
            ("n max (+)", "n_max_pos", ""),
            ("n max (-)", "n_max_neg", ""),
            ("V max", "V_max", "m/s"),
            ("Dive Speed Factor", "dive_factor", "")
        ]

        # ====================================================
        # CREATE INPUT FIELDS
        # ====================================================

        for i, (label, key, unit) in enumerate(parameters, start=1):

            ttk.Label(
                control_frame,
                text=label
            ).grid(
                row=i,
                column=0,
                sticky="w",
                pady=5
            )

            entry = ttk.Entry(
                control_frame,
                width=15
            )

            entry.insert(
                0,
                str(DEFAULTS[key])
            )

            entry.grid(
                row=i,
                column=1,
                padx=(15, 5),
                pady=5
            )

            self.entries[key] = entry

            if unit:
                ttk.Label(
                    control_frame,
                    text=unit
                ).grid(
                    row=i,
                    column=2,
                    sticky="w"
                )

        # ====================================================
        # BUTTONS
        # ====================================================

        button_frame = ttk.Frame(control_frame)

        button_frame.grid(
            row=11,
            column=0,
            columnspan=3,
            pady=20
        )

        update_button = ttk.Button(
            button_frame,
            text="Update Diagram",
            command=self.update_plot
        )

        update_button.grid(
            row=0,
            column=0,
            padx=5
        )

        reset_button = ttk.Button(
            button_frame,
            text="Reset Defaults",
            command=self.reset_defaults
        )

        reset_button.grid(
            row=0,
            column=1,
            padx=5
        )

        export_button = ttk.Button(
            button_frame,
            text="Export Plot",
            command=self.export_plot
        )

        export_button.grid(
            row=0,
            column=2,
            padx=5
        )

        # ====================================================
        # CHARACTERISTIC VELOCITIES
        # ====================================================

        velocity_label = ttk.Label(
            control_frame,
            text="Characteristic Velocities",
            font=("Arial", 13, "bold")
        )

        velocity_label.grid(
            row=12,
            column=0,
            columnspan=3,
            pady=(20, 10)
        )

        # Velocity output variables

        self.velocity_vars = {
            "Vs": tk.StringVar(value="--"),
            "Va": tk.StringVar(value="--"),
            "Vg": tk.StringVar(value="--"),
            "Vmax": tk.StringVar(value="--"),
            "Vd": tk.StringVar(value="--")
        }

        velocity_names = [
            ("Vs", "Stall Speed"),
            ("Va", "Maneuvering Speed"),
            ("Vg", "Negative-g Stall Speed"),
            ("Vmax", "Maximum Speed"),
            ("Vd", "Dive Speed")
        ]

        for i, (key, name) in enumerate(
            velocity_names,
            start=13
        ):

            ttk.Label(
                control_frame,
                text=name
            ).grid(
                row=i,
                column=0,
                sticky="w",
                pady=3
            )

            ttk.Label(
                control_frame,
                textvariable=self.velocity_vars[key],
                font=("Arial", 10, "bold")
            ).grid(
                row=i,
                column=1,
                sticky="w",
                padx=15
            )

        # ====================================================
        # RIGHT SIDE: MATPLOTLIB FIGURE
        # ====================================================

        self.figure, self.ax = plt.subplots(
            figsize=(9, 7)
        )

        self.canvas = FigureCanvasTkAgg(
            self.figure,
            master=root
        )

        self.canvas.get_tk_widget().grid(
            row=0,
            column=1,
            sticky="nsew",
            padx=10,
            pady=10
        )

        # ====================================================
        # INITIAL PLOT
        # ====================================================

        self.update_plot()


    # ========================================================
    # READ INPUTS
    # ========================================================

    def get_parameters(self):

        try:

            params = {}

            for key, entry in self.entries.items():

                params[key] = float(
                    entry.get()
                )

            return params

        except ValueError:

            messagebox.showerror(
                "Invalid Input",
                "Please enter valid numerical values."
            )

            return None


    # ========================================================
    # CALCULATE V-n DIAGRAM
    # ========================================================

    def calculate(self, p):

        # ----------------------------------------------------
        # Aircraft weight
        # ----------------------------------------------------

        W = p["mass"] * G

        S = p["S"]
        rho = p["rho"]

        CL_pos = p["CL_max_pos"]
        CL_neg = p["CL_max_neg"]

        n_pos = p["n_max_pos"]
        n_neg = p["n_max_neg"]

        V_max = p["V_max"]

        V_d = p["dive_factor"] * V_max

        # ----------------------------------------------------
        # Positive 1-g stall speed
        # ----------------------------------------------------

        V_s = np.sqrt(
            (2 * W) /
            (rho * S * CL_pos)
        )

        # ----------------------------------------------------
        # Positive maneuvering speed
        # ----------------------------------------------------

        V_a = V_s * np.sqrt(
            n_pos
        )

        # ----------------------------------------------------
        # Negative-g stall speed
        # ----------------------------------------------------

        V_g = np.sqrt(
            (2 * abs(n_neg) * W) /
            (
                rho *
                S *
                abs(CL_neg)
            )
        )

        # ----------------------------------------------------
        # Velocity array
        # ----------------------------------------------------

        V = np.arange(
            0,
            V_d + 0.1,
            0.1
        )

        # ----------------------------------------------------
        # Positive stall boundary
        # ----------------------------------------------------

        n_curve_pos = (
            rho *
            V**2 *
            S *
            CL_pos
        ) / (2 * W)

        # ----------------------------------------------------
        # Negative stall boundary
        # ----------------------------------------------------

        n_curve_neg = (
            rho *
            V**2 *
            S *
            CL_neg
        ) / (2 * W)

        # ----------------------------------------------------
        # Only show positive stall curve up to Va
        # ----------------------------------------------------

        n_curve_pos[
            V > V_a
        ] = np.nan

        # ----------------------------------------------------
        # Only show negative stall curve up to Vg
        # ----------------------------------------------------

        n_curve_neg[
            V > V_g
        ] = np.nan

        return {
            "W": W,
            "V": V,
            "n_pos_curve": n_curve_pos,
            "n_neg_curve": n_curve_neg,
            "Vs": V_s,
            "Va": V_a,
            "Vg": V_g,
            "Vmax": V_max,
            "Vd": V_d,
            "n_pos": n_pos,
            "n_neg": n_neg
        }


    # ========================================================
    # UPDATE PLOT
    # ========================================================

    def update_plot(self):

        p = self.get_parameters()

        if p is None:
            return

        # ----------------------------------------------------
        # Check physical inputs
        # ----------------------------------------------------

        if p["mass"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Mass must be greater than zero."
            )
            return

        if p["S"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Wing area must be greater than zero."
            )
            return

        if p["rho"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Air density must be greater than zero."
            )
            return

        if p["CL_max_pos"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Positive CLmax must be greater than zero."
            )
            return

        if p["CL_max_neg"] >= 0:
            messagebox.showerror(
                "Invalid Input",
                "Negative CLmax must be less than zero."
            )
            return

        if p["n_max_pos"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Positive maximum load factor must be greater than zero."
            )
            return

        if p["n_max_neg"] >= 0:
            messagebox.showerror(
                "Invalid Input",
                "Negative maximum load factor must be less than zero."
            )
            return

        if p["V_max"] <= 0:
            messagebox.showerror(
                "Invalid Input",
                "Vmax must be greater than zero."
            )
            return

        if p["dive_factor"] <= 1:
            messagebox.showerror(
                "Invalid Input",
                "Dive speed factor should be greater than 1."
            )
            return

        # ----------------------------------------------------
        # Calculate
        # ----------------------------------------------------

        data = self.calculate(p)

        V = data["V"]
        n_pos_curve = data["n_pos_curve"]
        n_neg_curve = data["n_neg_curve"]

        V_s = data["Vs"]
        V_a = data["Va"]
        V_g = data["Vg"]
        V_max = data["Vmax"]
        V_d = data["Vd"]

        n_pos = data["n_pos"]
        n_neg = data["n_neg"]

        # ====================================================
        # CLEAR PREVIOUS PLOT
        # ====================================================

        self.ax.clear()

        # ====================================================
        # PLOT POSITIVE STALL CURVE
        # ====================================================

        self.ax.plot(
            V,
            n_pos_curve,
            linewidth=2.5,
            label="Positive stall boundary"
        )

        # ====================================================
        # PLOT NEGATIVE STALL CURVE
        # ====================================================

        self.ax.plot(
            V,
            n_neg_curve,
            linewidth=2.5,
            label="Negative stall boundary"
        )

        # ====================================================
        # MAXIMUM POSITIVE LOAD FACTOR
        # ====================================================

        self.ax.plot(
            [V_a, V_d],
            [n_pos, n_pos],
            linestyle="--",
            linewidth=2,
            label=r"$n_{max,+}$"
        )

        # ====================================================
        # MAXIMUM NEGATIVE LOAD FACTOR
        # ====================================================

        self.ax.plot(
            [V_g, V_d],
            [n_neg, n_neg],
            linestyle="--",
            linewidth=2,
            label=r"$n_{max,-}$"
        )

        # ====================================================
        # Vmax LINE
        # ====================================================

        self.ax.plot(
            [V_max, V_max],
            [n_neg, n_pos],
            linestyle="--",
            linewidth=1.5,
            label=r"$V_{max}$"
        )

        # ====================================================
        # Vd LINE
        # ====================================================

        self.ax.plot(
            [V_d, V_d],
            [n_neg, n_pos],
            linewidth=3,
            label=r"$V_d$"
        )

        # ====================================================
        # CHARACTERISTIC POINTS
        # ====================================================

        self.ax.scatter(
            V_s,
            1.0,
            s=70,
            zorder=5
        )

        self.ax.scatter(
            V_a,
            n_pos,
            s=70,
            zorder=5
        )

        self.ax.scatter(
            V_g,
            n_neg,
            s=70,
            zorder=5
        )

        self.ax.scatter(
            V_max,
            0,
            s=70,
            zorder=5
        )

        self.ax.scatter(
            V_d,
            0,
            s=70,
            zorder=5
        )

        # ====================================================
        # LABELS
        # ====================================================

        self.ax.annotate(
            f"$V_s$ = {V_s:.2f} m/s",
            xy=(V_s, 1.0),
            xytext=(V_s + 0.5, 1.2),
            arrowprops=dict(
                arrowstyle="->"
            )
        )

        self.ax.annotate(
            f"$V_a$ = {V_a:.2f} m/s",
            xy=(V_a, n_pos),
            xytext=(V_a - 5, n_pos + 0.3),
            arrowprops=dict(
                arrowstyle="->"
            )
        )

        self.ax.annotate(
            f"$V_g$ = {V_g:.2f} m/s",
            xy=(V_g, n_neg),
            xytext=(V_g - 5, n_neg - 0.5),
            arrowprops=dict(
                arrowstyle="->"
            )
        )

        self.ax.annotate(
            f"$V_{{max}}$ = {V_max:.2f} m/s",
            xy=(V_max, 0),
            xytext=(V_max + 0.5, 0.4),
            arrowprops=dict(
                arrowstyle="->"
            )
        )

        self.ax.annotate(
            f"$V_d$ = {V_d:.2f} m/s",
            xy=(V_d, 0),
            xytext=(V_d + 0.5, 0.4),
            arrowprops=dict(
                arrowstyle="->"
            )
        )

        # ====================================================
        # AXIS SETTINGS
        # ====================================================

        self.ax.set_xlabel(
            "Velocity, V [m/s]",
            fontsize=12,
            fontweight="bold"
        )

        self.ax.set_ylabel(
            "Load Factor, n",
            fontsize=12,
            fontweight="bold"
        )

        self.ax.set_title(
            "UAV V-n Diagram",
            fontsize=15,
            fontweight="bold"
        )

        self.ax.set_xlim(
            0,
            V_d + 5
        )

        self.ax.set_ylim(
            n_neg - 0.5,
            n_pos + 0.5
        )

        self.ax.axhline(
            0,
            linewidth=0.8
        )

        self.ax.grid(
            True,
            alpha=0.3
        )

        self.ax.legend(
            loc="best"
        )

        self.figure.tight_layout()

        # ====================================================
        # UPDATE VELOCITY OUTPUTS
        # ====================================================

        self.velocity_vars["Vs"].set(
            f"{V_s:.2f} m/s"
        )

        self.velocity_vars["Va"].set(
            f"{V_a:.2f} m/s"
        )

        self.velocity_vars["Vg"].set(
            f"{V_g:.2f} m/s"
        )

        self.velocity_vars["Vmax"].set(
            f"{V_max:.2f} m/s"
        )

        self.velocity_vars["Vd"].set(
            f"{V_d:.2f} m/s"
        )

        # ====================================================
        # REDRAW
        # ====================================================

        self.canvas.draw()


    # ========================================================
    # RESET TO DEFAULT VALUES
    # ========================================================

    def reset_defaults(self):

        for key, entry in self.entries.items():

            entry.delete(
                0,
                tk.END
            )

            entry.insert(
                0,
                str(DEFAULTS[key])
            )

        self.update_plot()


# ============================================================
# RUN APPLICATION
# ============================================================

if __name__ == "__main__":

    root = tk.Tk()

    app = VNDiagramGUI(root)

    root.mainloop()
