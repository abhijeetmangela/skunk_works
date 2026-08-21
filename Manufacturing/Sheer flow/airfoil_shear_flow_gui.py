import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ============================================================
# Airfoil parsing
# ============================================================

def read_dat_file(filename):
    """
    Read an airfoil .dat file.

    The parser accepts:
      - normal whitespace-separated x y coordinates
      - a text/header line at the beginning
      - blank lines
      - extra columns (only first two numeric columns are used)
    """
    points = []

    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().replace(",", " ").split()

            if len(parts) < 2:
                continue

            try:
                x = float(parts[0])
                y = float(parts[1])
                points.append([x, y])
            except ValueError:
                # Ignore airfoil name/header lines
                continue

    data = np.asarray(points, dtype=float)

    if len(data) < 3:
        raise ValueError(
            "The selected file does not contain enough numeric x-y coordinates."
        )

    # Remove consecutive duplicate points
    keep = np.ones(len(data), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(data, axis=0), axis=1) > 1e-14
    data = data[keep]

    if len(data) < 3:
        raise ValueError("Not enough unique airfoil points.")

    return data


# ============================================================
# Shear-flow calculation
# ============================================================

def find_shear(data2, a_boom2, Vy):
    """
    Equivalent of the MATLAB find_shear() function.

    data2      : 2 x L coordinate array
    a_boom2    : L-1 boom-area array
    Vy         : applied vertical shear force [N]

    Returns:
        qs : shear flow at the boom stations [N/m]
    """
    L = data2.shape[1]

    if len(a_boom2) != L - 1:
        raise ValueError("Boom-area array must have length L-1.")

    if np.any(a_boom2 <= 0):
        raise ValueError(
            "All boom areas must be positive. Check the airfoil geometry."
        )

    cg_new = np.array([
        np.sum(data2[0, :L - 1] * a_boom2) / np.sum(a_boom2),
        np.sum(data2[1, :L - 1] * a_boom2) / np.sum(a_boom2)
    ])

    data3 = data2 - cg_new[:, None]

    Ixy2 = np.sum(
        data3[0, :L - 1] *
        data3[1, :L - 1] *
        a_boom2
    )

    Ix2 = np.sum(
        data3[1, :L - 1] ** 2 *
        a_boom2
    )

    Iy2 = np.sum(
        data3[0, :L - 1] ** 2 *
        a_boom2
    )

    den2 = Ix2 * Iy2 - Ixy2 ** 2

    if abs(den2) < 1e-30:
        raise ValueError(
            "Section inertia matrix is singular. "
            "Check the airfoil coordinates and boom areas."
        )

    Kxy2 = Ixy2 / den2
    Ky2 = Iy2 / den2

    Qy2 = np.zeros(L - 1)
    Qx2 = np.zeros(L - 1)
    qs2 = np.zeros(L - 1)

    for i in range(L - 2):
        Qy2[i + 1] = (
            Qy2[i] +
            data3[0, i + 1] * a_boom2[i + 1]
        )

        Qx2[i + 1] = (
            Qx2[i] +
            data3[1, i + 1] * a_boom2[i + 1]
        )

        qs2[i + 1] = Vy * (
            Kxy2 * Qy2[i + 1] -
            Ky2 * Qx2[i + 1]
        )

    return qs2


def calculate_model(
    raw_data,
    chord,
    shear_force,
    rib_spacing,
    upper_indices,
    lower_indices,
    stringer_area=10e-6,
    skin_thickness=0.0005,
):
    """
    Main structural calculation corresponding to the MATLAB script.
    """

    # Scale normalized airfoil coordinates
    data = raw_data * chord
    L = len(data)

    if L < 4:
        raise ValueError("At least four airfoil coordinate points are required.")

    # Panel lengths
    dis = np.linalg.norm(np.diff(data, axis=0), axis=1)

    if np.any(dis <= 0):
        raise ValueError("Airfoil contains zero-length panels.")

    # Panel midpoints
    mid_data = 0.5 * (data[:-1] + data[1:])

    # Geometric centroid used in the original MATLAB code
    cg = np.sum(mid_data * dis[:, None], axis=0) / np.sum(dis)

    # Shift coordinates
    shifted = data - cg

    # Neutral-axis estimate
    I_num = np.sum(
        (shifted[:-1, 0] + shifted[1:, 0]) *
        (shifted[:-1, 1] + shifted[1:, 1]) *
        0.25 * dis
    )

    I_den = np.sum(
        (shifted[:-1, 1] + shifted[1:, 1]) ** 2 *
        0.25 * dis
    )

    if abs(I_den) < 1e-30:
        raise ValueError("Unable to calculate the neutral-axis slope.")

    tan_alpha = I_num / I_den

    # Distance of each point from neutral axis
    p_d = np.abs(
        shifted[:, 1] -
        tan_alpha * shifted[:, 0]
    ) / np.sqrt(1 + tan_alpha ** 2)

    p_d = np.maximum(p_d, 1e-15)

    # Boom areas
    a_boom = np.zeros(L - 1)

    for i in range(L - 1):
        if i < L - 2:
            a_boom[i] += (
                (1 / 6) *
                dis[i] *
                skin_thickness *
                (2 + p_d[i + 1] / p_d[i])
            )

            a_boom[i + 1] += (
                (1 / 6) *
                dis[i] *
                skin_thickness *
                (2 + p_d[i] / p_d[i + 1])
            )
        else:
            a_boom[i] += (
                (1 / 6) *
                dis[i] *
                skin_thickness *
                (2 + p_d[0] / p_d[i])
            )

            a_boom[0] += (
                (1 / 6) *
                dis[i] *
                skin_thickness *
                (2 + p_d[i] / p_d[0])
            )

    data2 = data.T

    # Shear flow without stringers
    qs = find_shear(data2, a_boom, shear_force)

    # Closed-section correction
    q0 = -np.sum(qs * dis) / np.sum(dis)
    shear = qs + q0

    # Stringer requirement
    D = (
        np.pi ** 2 *
        70e9 *
        skin_thickness ** 3 *
        1.62
    ) / (
        12 *
        2.5 *
        1.5 *
        1.5 *
        1.5 *
        (1 - 0.33 ** 2)
    )

    Kss = 5.34 + (4 * chord ** 2) / (rib_spacing ** 2)
    N = Kss * D / chord ** 2

    # Add stringers
    a_boom2 = a_boom.copy()

    half = int(round(0.5 * (L - 1)))

    for idx in upper_indices:
        if idx < 1 or idx > half:
            raise ValueError(
                f"Upper stringer index {idx} is outside 1-{half}."
            )
        a_boom2[idx - 1] += stringer_area

    for idx in lower_indices:
        # User index is 1-based within the lower surface
        if idx < 1 or idx > L - half - 1:
            raise ValueError(
                f"Lower stringer index {idx} is outside the lower-surface range."
            )
        actual = half + idx - 1
        if actual >= L - 1:
            raise ValueError(
                f"Lower stringer index {idx} maps outside the boom array."
            )
        a_boom2[actual] += stringer_area

    qs2 = find_shear(data2, a_boom2, shear_force)
    q02 = -np.sum(qs2 * dis) / np.sum(dis)
    shear2 = qs2 + q02

    # Panel data between stringers
    all_boundaries = [1] + sorted(set(upper_indices)) + [half]
    upper_boundaries = sorted(set(all_boundaries))

    lower_boundaries = [1] + sorted(set(lower_indices)) + [L - half - 1]
    lower_boundaries = sorted(set(lower_boundaries))

    q_max2 = []
    dis2 = []

    # Upper panels
    for j in range(len(upper_boundaries) - 1):
        s = upper_boundaries[j] - 1
        e = upper_boundaries[j + 1] - 1

        if e > s:
            q_max2.append(np.max(np.abs(shear2[s:e])))
            dis2.append(np.sum(dis[s:e]))

    # Lower panels
    for j in range(len(lower_boundaries) - 1):
        s = half + lower_boundaries[j] - 1
        e = half + lower_boundaries[j + 1] - 1

        if e > s and e <= len(dis):
            q_max2.append(np.max(np.abs(shear2[s:e])))
            dis2.append(np.sum(dis[s:e]))

    q_max2 = np.asarray(q_max2)
    dis2 = np.asarray(dis2)

    if len(dis2):
        Kss2 = 5.34 + (4 * dis2 ** 2) / rib_spacing ** 2
        N2 = Kss2 * D / dis2 ** 2
        stringer_margin = N2 / np.maximum(q_max2, 1e-30)
    else:
        N2 = np.array([])
        stringer_margin = np.array([])

    return {
        "data": data,
        "cg": cg,
        "dis": dis,
        "a_boom": a_boom,
        "a_boom2": a_boom2,
        "tan_alpha": tan_alpha,
        "qs": qs,
        "shear": shear,
        "qs2": qs2,
        "shear2": shear2,
        "N": N,
        "N2": N2,
        "q_max2": q_max2,
        "stringer_margin": stringer_margin,
        "half": half,
    }


# ============================================================
# GUI
# ============================================================

class ShearFlowGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Airfoil Shear Flow & Stringer Analysis")
        self.root.geometry("1250x820")
        self.root.minsize(1050, 700)

        self.data = None
        self.result = None

        self.build_gui()

    def build_gui(self):
        # Main layout
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        control = ttk.Frame(self.root, padding=12)
        control.grid(row=0, column=0, sticky="ns")

        plot_frame = ttk.Frame(self.root, padding=(0, 12, 12, 12))
        plot_frame.grid(row=0, column=1, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        # ----------------------------------------------------
        # File selection
        # ----------------------------------------------------

        file_box = ttk.LabelFrame(
            control,
            text="Airfoil DAT File",
            padding=10
        )
        file_box.pack(fill="x", pady=(0, 10))

        self.file_var = tk.StringVar()

        ttk.Entry(
            file_box,
            textvariable=self.file_var,
            width=34
        ).pack(fill="x", pady=(0, 7))

        ttk.Button(
            file_box,
            text="Browse .dat File",
            command=self.browse_file
        ).pack(fill="x")

        self.file_status = ttk.Label(
            file_box,
            text="No file selected",
            foreground="gray"
        )
        self.file_status.pack(anchor="w", pady=(7, 0))

        # ----------------------------------------------------
        # General inputs
        # ----------------------------------------------------

        input_box = ttk.LabelFrame(
            control,
            text="Analysis Inputs",
            padding=10
        )
        input_box.pack(fill="x", pady=(0, 10))

        self.chord_var = tk.StringVar(value="0.234")
        self.shear_var = tk.StringVar(value="100")
        self.rib_var = tk.StringVar(value="0.2")
        self.skin_var = tk.StringVar(value="0.0005")
        self.stringer_var = tk.StringVar(value="1.0e-5")

        self.add_entry(input_box, "Chord c [m]", self.chord_var)
        self.add_entry(input_box, "Shear force Vy [N]", self.shear_var)
        self.add_entry(input_box, "Rib spacing a [m]", self.rib_var)
        self.add_entry(input_box, "Skin thickness [m]", self.skin_var)
        self.add_entry(input_box, "Stringer area [m²]", self.stringer_var)

        # ----------------------------------------------------
        # Stringer inputs
        # ----------------------------------------------------

        stringer_box = ttk.LabelFrame(
            control,
            text="Stringer Locations",
            padding=10
        )
        stringer_box.pack(fill="x", pady=(0, 10))

        ttk.Label(
            stringer_box,
            text="Enter 1-based point indices separated by spaces."
        ).pack(anchor="w", pady=(0, 6))

        self.upper_var = tk.StringVar()
        self.lower_var = tk.StringVar()

        ttk.Label(stringer_box, text="Upper surface:").pack(anchor="w")
        ttk.Entry(
            stringer_box,
            textvariable=self.upper_var
        ).pack(fill="x", pady=(2, 7))

        ttk.Label(stringer_box, text="Lower surface:").pack(anchor="w")
        ttk.Entry(
            stringer_box,
            textvariable=self.lower_var
        ).pack(fill="x", pady=(2, 7))

        # ----------------------------------------------------
        # Buttons
        # ----------------------------------------------------

        ttk.Button(
            control,
            text="Run Analysis",
            command=self.run_analysis
        ).pack(fill="x", pady=(3, 5))

        ttk.Button(
            control,
            text="Clear Results",
            command=self.clear_results
        ).pack(fill="x", pady=5)

        ttk.Button(
            control,
            text="Save Current Plot",
            command=self.save_plot
        ).pack(fill="x", pady=5)

        # ----------------------------------------------------
        # Results
        # ----------------------------------------------------

        result_box = ttk.LabelFrame(
            control,
            text="Results",
            padding=10
        )
        result_box.pack(fill="both", expand=True, pady=(10, 0))

        self.result_text = tk.Text(
            result_box,
            width=42,
            height=18,
            wrap="word",
            font=("DejaVu Sans Mono", 9)
        )
        self.result_text.pack(fill="both", expand=True)

        # ----------------------------------------------------
        # Matplotlib figure
        # ----------------------------------------------------

        self.fig, self.ax = plt.subplots(
            2,
            2,
            figsize=(9, 7)
        )

        self.fig.tight_layout(pad=3)

        self.canvas = FigureCanvasTkAgg(
            self.fig,
            master=plot_frame
        )

        self.canvas.get_tk_widget().grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        self.show_empty_plot()

    def add_entry(self, parent, label, variable):
        ttk.Label(parent, text=label).pack(anchor="w")
        ttk.Entry(
            parent,
            textvariable=variable
        ).pack(fill="x", pady=(2, 7))

    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="Select Airfoil DAT File",
            filetypes=[
                ("DAT files", "*.dat"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )

        if not filename:
            return

        try:
            data = read_dat_file(filename)

            self.data = data
            self.file_var.set(filename)

            self.file_status.config(
                text=f"{len(data)} coordinate points loaded",
                foreground="green"
            )

            self.show_airfoil_preview()

        except Exception as e:
            messagebox.showerror(
                "File Error",
                str(e)
            )

    def parse_indices(self, text):
        text = text.strip()

        if not text:
            return []

        try:
            values = [int(x) for x in text.replace(",", " ").split()]
        except ValueError:
            raise ValueError(
                "Stringer indices must be integers separated by spaces."
            )

        if any(x < 1 for x in values):
            raise ValueError(
                "Stringer indices must be 1 or greater."
            )

        return values

    def run_analysis(self):
        try:
            if self.data is None:
                raise ValueError("Please select an airfoil .dat file first.")

            chord = float(self.chord_var.get())
            shear_force = float(self.shear_var.get())
            rib_spacing = float(self.rib_var.get())
            skin_thickness = float(self.skin_var.get())
            stringer_area = float(self.stringer_var.get())

            if chord <= 0:
                raise ValueError("Chord must be greater than zero.")

            if rib_spacing <= 0:
                raise ValueError("Rib spacing must be greater than zero.")

            if skin_thickness <= 0:
                raise ValueError("Skin thickness must be greater than zero.")

            if stringer_area < 0:
                raise ValueError("Stringer area cannot be negative.")

            upper_indices = self.parse_indices(
                self.upper_var.get()
            )

            lower_indices = self.parse_indices(
                self.lower_var.get()
            )

            self.result = calculate_model(
                self.data,
                chord,
                shear_force,
                rib_spacing,
                upper_indices,
                lower_indices,
                stringer_area,
                skin_thickness
            )

            self.update_results(
                chord,
                shear_force,
                rib_spacing,
                upper_indices,
                lower_indices
            )

            self.update_plots(
                upper_indices,
                lower_indices
            )

        except Exception as e:
            messagebox.showerror(
                "Analysis Error",
                str(e)
            )

    def update_results(
        self,
        chord,
        shear_force,
        rib_spacing,
        upper_indices,
        lower_indices
    ):
        r = self.result

        max_q = np.max(np.abs(r["shear"]))
        max_q_stringers = np.max(np.abs(r["shear2"]))

        no_stringer_required = r["N"] > max_q

        self.result_text.delete("1.0", tk.END)

        lines = [
            "AIRFOIL SHEAR FLOW ANALYSIS",
            "=" * 34,
            "",
            f"Number of points : {len(r['data'])}",
            f"Chord            : {chord:.6g} m",
            f"Shear force Vy   : {shear_force:.6g} N",
            f"Rib spacing      : {rib_spacing:.6g} m",
            "",
            "SECTION PROPERTIES",
            "-" * 34,
            f"Centroid x       : {r['cg'][0]:.8e} m",
            f"Centroid y       : {r['cg'][1]:.8e} m",
            f"tan(alpha)       : {r['tan_alpha']:.8e}",
            "",
            "WITHOUT STRINGERS",
            "-" * 34,
            f"Maximum |q|      : {max_q:.8e} N/m",
            f"Critical N       : {r['N']:.8e} N/m",
            "",
            (
                "STATUS: No stringer needed"
                if no_stringer_required
                else "STATUS: Stringer required"
            ),
            "",
            "STRINGERS",
            "-" * 34,
            f"Upper indices    : {upper_indices}",
            f"Lower indices    : {lower_indices}",
            f"Maximum |q|      : {max_q_stringers:.8e} N/m",
        ]

        if len(r["stringer_margin"]) > 0:
            lines.extend([
                "",
                "Panel N/q_max:",
                np.array2string(
                    r["stringer_margin"],
                    precision=4,
                    separator=", "
                )
            ])

        self.result_text.insert(
            "1.0",
            "\n".join(lines)
        )

    def show_empty_plot(self):
        for a in self.ax.flat:
            a.clear()
            a.grid(True, alpha=0.25)

        self.ax[0, 0].set_title("Airfoil")
        self.ax[0, 1].set_title("Shear Flow - Skin Only")
        self.ax[1, 0].set_title("Shear Flow - With Stringers")
        self.ax[1, 1].set_title("Airfoil & Stringer Locations")

        self.canvas.draw()

    def show_airfoil_preview(self):
        data = self.data

        for a in self.ax.flat:
            a.clear()
            a.grid(True, alpha=0.25)

        self.ax[0, 0].plot(
            data[:, 0],
            data[:, 1],
            linewidth=1.5
        )
        self.ax[0, 0].set_title("Loaded Airfoil")
        self.ax[0, 0].set_xlabel("x")
        self.ax[0, 0].set_ylabel("y")
        self.ax[0, 0].axis("equal")

        self.canvas.draw()

    def update_plots(self, upper_indices, lower_indices):
        r = self.result
        data = r["data"]
        cg = r["cg"]
        half = r["half"]

        x = data[:, 0] + cg[0]
        y = data[:, 1] + cg[1]

        for a in self.ax.flat:
            a.clear()
            a.grid(True, alpha=0.25)

        # ----------------------------------------------------
        # 1. Airfoil
        # ----------------------------------------------------

        self.ax[0, 0].plot(
            x,
            y,
            linewidth=1.7
        )

        self.ax[0, 0].set_title("Airfoil")
        self.ax[0, 0].set_xlabel("x [m]")
        self.ax[0, 0].set_ylabel("y [m]")
        self.ax[0, 0].axis("equal")

        # ----------------------------------------------------
        # 2. Skin-only shear flow
        # ----------------------------------------------------

        self.ax[0, 1].plot(
            x[:half],
            -r["shear"][:half],
            linewidth=1.8,
            label="Upper"
        )

        self.ax[0, 1].plot(
            x[half:-1],
            -r["shear"][half:],
            linewidth=1.8,
            label="Lower"
        )

        self.ax[0, 1].set_title(
            "Shear Flow - Skin Only"
        )
        self.ax[0, 1].set_xlabel("x [m]")
        self.ax[0, 1].set_ylabel("q [N/m]")
        self.ax[0, 1].legend()

        # ----------------------------------------------------
        # 3. Shear flow with stringers
        # ----------------------------------------------------

        self.ax[1, 0].plot(
            x[:half],
            -r["shear2"][:half],
            linewidth=1.8,
            label="Upper"
        )

        self.ax[1, 0].plot(
            x[half:-1],
            -r["shear2"][half:],
            linewidth=1.8,
            label="Lower"
        )

        self.ax[1, 0].set_title(
            "Shear Flow - With Stringers"
        )
        self.ax[1, 0].set_xlabel("x [m]")
        self.ax[1, 0].set_ylabel("q [N/m]")
        self.ax[1, 0].legend()

        # ----------------------------------------------------
        # 4. Stringer locations
        # ----------------------------------------------------

        self.ax[1, 1].plot(
            x,
            y,
            linewidth=1.5
        )

        # Upper stringers
        for idx in upper_indices:
            if 1 <= idx <= half:
                j = idx - 1
                self.ax[1, 1].scatter(
                    x[j],
                    y[j],
                    s=50,
                    marker="o",
                    label="Upper stringer"
                    if idx == upper_indices[0]
                    else None
                )

        # Lower stringers
        for idx in lower_indices:
            j = half + idx - 1
            if 0 <= j < len(x):
                self.ax[1, 1].scatter(
                    x[j],
                    y[j],
                    s=50,
                    marker="s",
                    label="Lower stringer"
                    if idx == lower_indices[0]
                    else None
                )

        self.ax[1, 1].set_title(
            "Airfoil & Stringer Locations"
        )
        self.ax[1, 1].set_xlabel("x [m]")
        self.ax[1, 1].set_ylabel("y [m]")
        self.ax[1, 1].axis("equal")

        handles, labels = self.ax[1, 1].get_legend_handles_labels()
        if handles:
            self.ax[1, 1].legend()

        self.fig.tight_layout(pad=3)
        self.canvas.draw()

    def save_plot(self):
        if self.result is None:
            messagebox.showwarning(
                "No Results",
                "Run the analysis before saving a plot."
            )
            return

        filename = filedialog.asksaveasfilename(
            title="Save Plot",
            defaultextension=".png",
            filetypes=[
                ("PNG", "*.png"),
                ("PDF", "*.pdf"),
                ("SVG", "*.svg")
            ]
        )

        if filename:
            self.fig.savefig(
                filename,
                dpi=300,
                bbox_inches="tight"
            )

    def clear_results(self):
        self.result = None

        self.result_text.delete(
            "1.0",
            tk.END
        )

        self.show_empty_plot()


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    root = tk.Tk()

    # Use a modern ttk theme when available
    try:
        style = ttk.Style()
        if "clam" in style.theme_names():
            style.theme_use("clam")
    except Exception:
        pass

    app = ShearFlowGUI(root)

    root.mainloop()
