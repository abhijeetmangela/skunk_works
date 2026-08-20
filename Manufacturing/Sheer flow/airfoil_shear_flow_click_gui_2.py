import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ============================================================
# DAT FILE READER
# ============================================================

def read_dat_file(filename):
    """
    Reads an airfoil .dat file.

    The first two numeric columns are interpreted as x and y.
    Header/text lines and blank lines are ignored.
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
                continue

    data = np.asarray(points, dtype=float)

    if len(data) < 4:
        raise ValueError(
            "The DAT file does not contain enough x-y coordinate points."
        )

    # Remove consecutive duplicate points
    keep = np.ones(len(data), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(data, axis=0), axis=1) > 1e-14
    data = data[keep]

    return data


# ============================================================
# SHEAR FLOW FUNCTION
# ============================================================

def find_shear(data2, a_boom2, Vy):
    """
    Equivalent to the MATLAB find_shear() function.
    """

    L = data2.shape[1]

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
            "Section inertia matrix is singular."
        )

    Kxy2 = Ixy2 / den2
    Ky2 = Iy2 / den2

    Qy2 = np.zeros(L - 1)
    Qx2 = np.zeros(L - 1)
    qs2 = np.zeros(L - 1)

    for i in range(L - 2):

        Qy2[i + 1] = (
            Qy2[i] +
            data3[0, i + 1] *
            a_boom2[i + 1]
        )

        Qx2[i + 1] = (
            Qx2[i] +
            data3[1, i + 1] *
            a_boom2[i + 1]
        )

        qs2[i + 1] = Vy * (
            Kxy2 * Qy2[i + 1]
            -
            Ky2 * Qx2[i + 1]
        )

    return qs2


# ============================================================
# MAIN STRUCTURAL CALCULATION
# ============================================================

def calculate_model(
    raw_data,
    chord,
    shear_force,
    rib_spacing,
    upper_indices,
    lower_indices,
    stringer_area,
    skin_thickness
):

    data = raw_data * chord
    L = len(data)

    dis = np.linalg.norm(
        np.diff(data, axis=0),
        axis=1
    )

    if np.any(dis <= 0):
        raise ValueError(
            "The DAT file contains zero-length panels."
        )

    # Panel midpoints
    mid_data = 0.5 * (
        data[:-1] + data[1:]
    )

    # Geometric centroid
    cg = (
        np.sum(
            mid_data * dis[:, None],
            axis=0
        )
        / np.sum(dis)
    )

    # Shift coordinates
    shifted = data - cg

    # Neutral axis
    I_num = np.sum(
        (shifted[:-1, 0] + shifted[1:, 0]) *
        (shifted[:-1, 1] + shifted[1:, 1]) *
        0.25 * dis
    )

    I_den = np.sum(
        (shifted[:-1, 1] + shifted[1:, 1]) ** 2 *
        0.25 * dis
    )

    tan_alpha = I_num / I_den

    # Distance from neutral axis
    p_d = np.abs(
        shifted[:, 1] -
        tan_alpha * shifted[:, 0]
    ) / np.sqrt(1 + tan_alpha ** 2)

    p_d = np.maximum(p_d, 1e-15)

    # --------------------------------------------------------
    # Boom areas
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Skin-only shear flow
    # --------------------------------------------------------

    qs = find_shear(
        data2,
        a_boom,
        shear_force
    )

    q0 = (
        -np.sum(qs * dis)
        / np.sum(dis)
    )

    shear = qs + q0

    # --------------------------------------------------------
    # Stringer requirement
    # --------------------------------------------------------

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

    Kss = (
        5.34 +
        4 * chord ** 2 /
        rib_spacing ** 2
    )

    N = (
        Kss * D /
        chord ** 2
    )

    # --------------------------------------------------------
    # Add stringer boom areas
    # --------------------------------------------------------

    a_boom2 = a_boom.copy()

    half = int(round(0.5 * (L - 1)))

    # Upper indices are supplied as 1-based point numbers
    for idx in upper_indices:

        j = idx - 1

        if j < 0 or j >= half:
            continue

        a_boom2[j] += stringer_area

    # Lower indices are relative to lower surface
    for idx in lower_indices:

        j = half + idx - 1

        if j < half or j >= L - 1:
            continue

        a_boom2[j] += stringer_area

    # --------------------------------------------------------
    # Shear flow with stringers
    # --------------------------------------------------------

    qs2 = find_shear(
        data2,
        a_boom2,
        shear_force
    )

    q02 = (
        -np.sum(qs2 * dis)
        / np.sum(dis)
    )

    shear2 = qs2 + q02

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
        "half": half
    }


# ============================================================
# GUI
# ============================================================

class ShearFlowGUI:

    def __init__(self, root):

        self.root = root

        self.root.title(
            "Airfoil Shear Flow & Stringer Analysis"
        )

        self.root.geometry(
            "1350x850"
        )

        self.root.minsize(
            1100,
            700
        )

        self.data = None
        self.result = None

        self.upper_indices = set()
        self.lower_indices = set()

        self.click_cid = None

        self.build_gui()


    # ========================================================
    # GUI LAYOUT
    # ========================================================

    def build_gui(self):

        self.root.columnconfigure(
            1,
            weight=1
        )

        self.root.rowconfigure(
            0,
            weight=1
        )

        # ----------------------------------------------------
        # Left panel
        # ----------------------------------------------------

        control = ttk.Frame(
            self.root,
            padding=12
        )

        control.grid(
            row=0,
            column=0,
            sticky="ns"
        )

        # ----------------------------------------------------
        # Right panel
        # ----------------------------------------------------

        plot_frame = ttk.Frame(
            self.root,
            padding=(0, 12, 12, 12)
        )

        plot_frame.grid(
            row=0,
            column=1,
            sticky="nsew"
        )

        plot_frame.rowconfigure(
            0,
            weight=1
        )

        plot_frame.columnconfigure(
            0,
            weight=1
        )

        # ====================================================
        # FILE
        # ====================================================

        file_box = ttk.LabelFrame(
            control,
            text="Airfoil DAT File",
            padding=10
        )

        file_box.pack(
            fill="x",
            pady=(0, 10)
        )

        self.file_var = tk.StringVar()

        ttk.Entry(
            file_box,
            textvariable=self.file_var,
            width=38
        ).pack(
            fill="x",
            pady=(0, 7)
        )

        ttk.Button(
            file_box,
            text="Browse .dat File",
            command=self.browse_file
        ).pack(
            fill="x"
        )

        self.file_status = ttk.Label(
            file_box,
            text="No file selected"
        )

        self.file_status.pack(
            anchor="w",
            pady=(6, 0)
        )

        # ====================================================
        # INPUTS
        # ====================================================

        input_box = ttk.LabelFrame(
            control,
            text="Analysis Inputs",
            padding=10
        )

        input_box.pack(
            fill="x",
            pady=(0, 10)
        )

        self.chord_var = tk.StringVar(
            value="0.234"
        )

        self.shear_var = tk.StringVar(
            value="100"
        )

        self.rib_var = tk.StringVar(
            value="0.2"
        )

        self.skin_var = tk.StringVar(
            value="0.0005"
        )

        self.stringer_var = tk.StringVar(
            value="1e-5"
        )

        self.add_entry(
            input_box,
            "Chord c [m]",
            self.chord_var
        )

        self.add_entry(
            input_box,
            "Shear force Vy [N]",
            self.shear_var
        )

        self.add_entry(
            input_box,
            "Rib spacing a [m]",
            self.rib_var
        )

        self.add_entry(
            input_box,
            "Skin thickness t [m]",
            self.skin_var
        )

        self.add_entry(
            input_box,
            "Stringer area [m²]",
            self.stringer_var
        )

        # ====================================================
        # STRINGER MODE
        # ====================================================

        stringer_box = ttk.LabelFrame(
            control,
            text="Stringer Placement",
            padding=10
        )

        stringer_box.pack(
            fill="x",
            pady=(0, 10)
        )

        ttk.Label(
            stringer_box,
            text=(
                "Click directly on the airfoil points.\n"
                "Left-click: add/remove stringer\n"
                "Right-click: clear selected point"
            ),
            wraplength=330
        ).pack(
            anchor="w",
            pady=(0, 8)
        )

        self.mode_var = tk.StringVar(
            value="Upper"
        )

        ttk.Label(
            stringer_box,
            text="Current surface:"
        ).pack(
            anchor="w"
        )

        mode_frame = ttk.Frame(
            stringer_box
        )

        mode_frame.pack(
            fill="x",
            pady=5
        )

        ttk.Radiobutton(
            mode_frame,
            text="Upper",
            variable=self.mode_var,
            value="Upper"
        ).pack(
            side="left",
            padx=(0, 15)
        )

        ttk.Radiobutton(
            mode_frame,
            text="Lower",
            variable=self.mode_var,
            value="Lower"
        ).pack(
            side="left"
        )

        self.upper_label = ttk.Label(
            stringer_box,
            text="Upper: None"
        )

        self.upper_label.pack(
            anchor="w",
            pady=(5, 0)
        )

        self.lower_label = ttk.Label(
            stringer_box,
            text="Lower: None"
        )

        self.lower_label.pack(
            anchor="w"
        )

        ttk.Button(
            stringer_box,
            text="Clear All Stringers",
            command=self.clear_stringers
        ).pack(
            fill="x",
            pady=(8, 0)
        )

        # ====================================================
        # BUTTONS
        # ====================================================

        ttk.Button(
            control,
            text="Run Analysis",
            command=self.run_analysis
        ).pack(
            fill="x",
            pady=(2, 5)
        )

        ttk.Button(
            control,
            text="Save Plot",
            command=self.save_plot
        ).pack(
            fill="x",
            pady=5
        )

        # ====================================================
        # RESULTS
        # ====================================================

        result_box = ttk.LabelFrame(
            control,
            text="Results",
            padding=8
        )

        result_box.pack(
            fill="both",
            expand=True,
            pady=(10, 0)
        )

        self.result_text = tk.Text(
            result_box,
            width=43,
            height=17,
            wrap="word",
            font=("DejaVu Sans Mono", 9)
        )

        self.result_text.pack(
            fill="both",
            expand=True
        )

        # ====================================================
        # MATPLOTLIB
        # ====================================================

        self.fig, self.ax = plt.subplots(
            2,
            2,
            figsize=(9, 7)
        )

        self.fig.tight_layout(
            pad=3
        )

        self.canvas = FigureCanvasTkAgg(
            self.fig,
            master=plot_frame
        )

        self.canvas.get_tk_widget().grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        # Connect mouse events
        self.click_cid = self.canvas.mpl_connect(
            "button_press_event",
            self.on_plot_click
        )

        self.show_empty_plot()


    # ========================================================
    # ENTRY HELPER
    # ========================================================

    def add_entry(
        self,
        parent,
        label,
        variable
    ):

        ttk.Label(
            parent,
            text=label
        ).pack(
            anchor="w"
        )

        ttk.Entry(
            parent,
            textvariable=variable
        ).pack(
            fill="x",
            pady=(2, 7)
        )


    # ========================================================
    # LOAD DAT
    # ========================================================

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

            data = read_dat_file(
                filename
            )

            self.data = data

            self.file_var.set(
                filename
            )

            self.file_status.config(
                text=f"{len(data)} points loaded"
            )

            self.upper_indices.clear()
            self.lower_indices.clear()

            self.update_stringer_labels()

            self.show_airfoil()

        except Exception as e:

            messagebox.showerror(
                "DAT File Error",
                str(e)
            )


    # ========================================================
    # AIRFOIL PLOT
    # ========================================================

    def show_airfoil(self):

        if self.data is None:
            return

        self.ax[0, 0].clear()

        x = self.data[:, 0]
        y = self.data[:, 1]

        self.ax[0, 0].plot(
            x,
            y,
            linewidth=1.5
        )

        self.ax[0, 0].scatter(
            x,
            y,
            s=14
        )

        self.ax[0, 0].set_title(
            "Click Points to Place Stringers"
        )

        self.ax[0, 0].set_xlabel(
            "x"
        )

        self.ax[0, 0].set_ylabel(
            "y"
        )

        self.ax[0, 0].axis(
            "equal"
        )

        self.ax[0, 0].grid(
            True,
            alpha=0.25
        )

        self.fig.tight_layout(
            pad=3
        )

        self.canvas.draw()


    # ========================================================
    # MOUSE CLICK
    # ========================================================

    def on_plot_click(self, event):

        if self.data is None:
            return

        # Only allow clicks in the first plot
        if event.inaxes != self.ax[0, 0]:
            return

        if event.xdata is None or event.ydata is None:
            return

        # ----------------------------------------------------
        # Right click: clear nearest selected stringer
        # ----------------------------------------------------

        if event.button == 3:

            self.remove_nearest_stringer(
                event.xdata,
                event.ydata
            )

            return

        # ----------------------------------------------------
        # Left click: add/remove nearest point
        # ----------------------------------------------------

        if event.button != 1:
            return

        points = self.data

        distances = np.sqrt(
            (points[:, 0] - event.xdata) ** 2 +
            (points[:, 1] - event.ydata) ** 2
        )

        nearest = int(
            np.argmin(distances)
        )

        self.toggle_stringer(
            nearest
        )


    # ========================================================
    # TOGGLE STRINGER
    # ========================================================

    def toggle_stringer(
        self,
        index
    ):

        L = len(self.data)

        half = int(
            round(0.5 * (L - 1))
        )

        # The stringer must be placed on a boom station
        if index >= L - 1:
            return

        # ----------------------------------------------------
        # Upper surface
        # ----------------------------------------------------

        if self.mode_var.get() == "Upper":

            if index >= half:
                messagebox.showwarning(
                    "Wrong Surface",
                    "The selected point is on the lower surface."
                )
                return

            point_number = index + 1

            if point_number in self.upper_indices:
                self.upper_indices.remove(
                    point_number
                )
            else:
                self.upper_indices.add(
                    point_number
                )

        # ----------------------------------------------------
        # Lower surface
        # ----------------------------------------------------

        else:

            if index < half:
                messagebox.showwarning(
                    "Wrong Surface",
                    "The selected point is on the upper surface."
                )
                return

            # Lower-surface index is relative to lower surface
            point_number = (
                index - half + 1
            )

            if point_number in self.lower_indices:
                self.lower_indices.remove(
                    point_number
                )
            else:
                self.lower_indices.add(
                    point_number
                )

        self.update_stringer_labels()
        self.show_airfoil()


    # ========================================================
    # REMOVE NEAREST STRINGER
    # ========================================================

    def remove_nearest_stringer(
        self,
        x,
        y
    ):

        candidates = []

        for idx in self.upper_indices:

            j = idx - 1

            if j < len(self.data):
                candidates.append(
                    ("upper", idx, j)
                )

        half = int(
            round(0.5 * (len(self.data) - 1))
        )

        for idx in self.lower_indices:

            j = half + idx - 1

            if j < len(self.data):
                candidates.append(
                    ("lower", idx, j)
                )

        if not candidates:
            return

        distances = []

        for surface, idx, j in candidates:

            d = (
                self.data[j, 0] - x
            ) ** 2 + (
                self.data[j, 1] - y
            ) ** 2

            distances.append(d)

        k = int(
            np.argmin(distances)
        )

        surface, idx, j = candidates[k]

        # Only remove if reasonably close
        if distances[k] > 0.01 ** 2:
            return

        if surface == "upper":
            self.upper_indices.remove(idx)
        else:
            self.lower_indices.remove(idx)

        self.update_stringer_labels()
        self.show_airfoil()


    # ========================================================
    # STRINGER LABELS
    # ========================================================

    def update_stringer_labels(self):

        upper = sorted(
            self.upper_indices
        )

        lower = sorted(
            self.lower_indices
        )

        self.upper_label.config(
            text=(
                "Upper: "
                + (
                    "None"
                    if not upper
                    else " ".join(
                        map(str, upper)
                    )
                )
            )
        )

        self.lower_label.config(
            text=(
                "Lower: "
                + (
                    "None"
                    if not lower
                    else " ".join(
                        map(str, lower)
                    )
                )
            )
        )


    # ========================================================
    # CLEAR STRINGERS
    # ========================================================

    def clear_stringers(self):

        self.upper_indices.clear()
        self.lower_indices.clear()

        self.update_stringer_labels()

        if self.data is not None:
            self.show_airfoil()


    # ========================================================
    # RUN ANALYSIS
    # ========================================================

    def run_analysis(self):

        try:

            if self.data is None:
                raise ValueError(
                    "Please load a .dat airfoil file first."
                )

            chord = float(
                self.chord_var.get()
            )

            shear_force = float(
                self.shear_var.get()
            )

            rib_spacing = float(
                self.rib_var.get()
            )

            skin_thickness = float(
                self.skin_var.get()
            )

            stringer_area = float(
                self.stringer_var.get()
            )

            if chord <= 0:
                raise ValueError(
                    "Chord must be greater than zero."
                )

            if rib_spacing <= 0:
                raise ValueError(
                    "Rib spacing must be greater than zero."
                )

            if skin_thickness <= 0:
                raise ValueError(
                    "Skin thickness must be greater than zero."
                )

            if stringer_area <= 0:
                raise ValueError(
                    "Stringer area must be greater than zero."
                )

            self.result = calculate_model(
                self.data,
                chord,
                shear_force,
                rib_spacing,
                sorted(self.upper_indices),
                sorted(self.lower_indices),
                stringer_area,
                skin_thickness
            )

            self.update_results(
                chord,
                shear_force,
                rib_spacing
            )

            self.update_plots()

        except Exception as e:

            messagebox.showerror(
                "Analysis Error",
                str(e)
            )


    # ========================================================
    # RESULTS
    # ========================================================

    def update_results(
        self,
        chord,
        shear_force,
        rib_spacing
    ):

        r = self.result

        max_q = np.max(
            np.abs(r["shear"])
        )

        max_q2 = np.max(
            np.abs(r["shear2"])
        )

        self.result_text.delete(
            "1.0",
            tk.END
        )

        text = f"""
AIRFOIL SHEAR FLOW ANALYSIS
================================

Number of points : {len(r["data"])}
Chord            : {chord:.6g} m
Shear force Vy   : {shear_force:.6g} N
Rib spacing      : {rib_spacing:.6g} m

SECTION PROPERTIES
--------------------------------
Centroid x       : {r["cg"][0]:.8e} m
Centroid y       : {r["cg"][1]:.8e} m
tan(alpha)       : {r["tan_alpha"]:.8e}

SKIN ONLY
--------------------------------
Maximum |q|      : {max_q:.8e} N/m
Critical N       : {r["N"]:.8e} N/m

{"NO STRINGER REQUIRED" if r["N"] > max_q else "STRINGER REQUIRED"}

STRINGER LOCATIONS
--------------------------------
Upper: {sorted(self.upper_indices) or "None"}
Lower: {sorted(self.lower_indices) or "None"}

WITH STRINGERS
--------------------------------
Maximum |q|      : {max_q2:.8e} N/m
"""

        self.result_text.insert(
            "1.0",
            text
        )


    # ========================================================
    # UPDATE ALL PLOTS
    # ========================================================

    def update_plots(self):

        r = self.result

        data = r["data"]
        cg = r["cg"]
        half = r["half"]

        x = data[:, 0] + cg[0]
        y = data[:, 1] + cg[1]

        for a in self.ax.flat:
            a.clear()
            a.grid(
                True,
                alpha=0.25
            )

        # ----------------------------------------------------
        # Plot 1: airfoil + stringers
        # ----------------------------------------------------

        self.ax[0, 0].plot(
            x,
            y,
            linewidth=1.5
        )

        self.ax[0, 0].scatter(
            x,
            y,
            s=12
        )

        # Upper stringers
        for idx in sorted(
            self.upper_indices
        ):

            j = idx - 1

            if j < len(x):

                self.ax[0, 0].scatter(
                    x[j],
                    y[j],
                    s=90,
                    marker="o",
                    facecolors="none",
                    linewidths=2
                )

                self.ax[0, 0].annotate(
                    str(idx),
                    (x[j], y[j]),
                    xytext=(5, 5),
                    textcoords="offset points"
                )

        # Lower stringers
        for idx in sorted(
            self.lower_indices
        ):

            j = half + idx - 1

            if j < len(x):

                self.ax[0, 0].scatter(
                    x[j],
                    y[j],
                    s=90,
                    marker="s",
                    facecolors="none",
                    linewidths=2
                )

                self.ax[0, 0].annotate(
                    str(idx),
                    (x[j], y[j]),
                    xytext=(5, 5),
                    textcoords="offset points"
                )

        self.ax[0, 0].set_title(
            "Airfoil & Stringer Locations"
        )

        self.ax[0, 0].set_xlabel(
            "x [m]"
        )

        self.ax[0, 0].set_ylabel(
            "y [m]"
        )

        self.ax[0, 0].axis(
            "equal"
        )

        # ----------------------------------------------------
        # Plot 2: skin only
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

        self.ax[0, 1].set_xlabel(
            "x [m]"
        )

        self.ax[0, 1].set_ylabel(
            "q [N/m]"
        )

        self.ax[0, 1].legend()

        # ----------------------------------------------------
        # Plot 3: with stringers
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

        self.ax[1, 0].set_xlabel(
            "x [m]"
        )

        self.ax[1, 0].set_ylabel(
            "q [N/m]"
        )

        self.ax[1, 0].legend()

        # ----------------------------------------------------
        # Plot 4: boom areas
        # ----------------------------------------------------

        self.ax[1, 1].plot(
            x[:-1],
            r["a_boom"],
            linewidth=1.5,
            label="Skin booms"
        )

        self.ax[1, 1].plot(
            x[:-1],
            r["a_boom2"],
            linewidth=1.5,
            label="With stringers"
        )

        self.ax[1, 1].set_title(
            "Boom Area Distribution"
        )

        self.ax[1, 1].set_xlabel(
            "x [m]"
        )

        self.ax[1, 1].set_ylabel(
            "Boom area [m²]"
        )

        self.ax[1, 1].legend()

        self.fig.tight_layout(
            pad=3
        )

        self.canvas.draw()


    # ========================================================
    # SAVE PLOT
    # ========================================================

    def save_plot(self):

        if self.result is None:

            messagebox.showwarning(
                "No Results",
                "Run the analysis first."
            )

            return

        filename = filedialog.asksaveasfilename(
            title="Save Analysis Plot",
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


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    root = tk.Tk()

    try:

        style = ttk.Style()

        if "clam" in style.theme_names():
            style.theme_use("clam")

    except Exception:
        pass

    app = ShearFlowGUI(root)

    root.mainloop()
