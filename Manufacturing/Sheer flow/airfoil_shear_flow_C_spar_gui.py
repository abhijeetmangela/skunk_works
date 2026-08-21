import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ============================================================
# DAT FILE READER
# ============================================================

def read_dat_file(filename):
    """Read the first two numeric columns from an airfoil DAT file."""

    points = []

    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().replace(",", " ").split()

            if len(parts) < 2:
                continue

            try:
                points.append([float(parts[0]), float(parts[1])])
            except ValueError:
                # Header/name lines are ignored.
                continue

    data = np.asarray(points, dtype=float)

    if len(data) < 5:
        raise ValueError(
            "The DAT file must contain at least five numeric x-y points."
        )

    # Remove consecutive duplicate points.
    keep = np.ones(len(data), dtype=bool)
    keep[1:] = np.linalg.norm(np.diff(data, axis=0), axis=1) > 1e-14
    data = data[keep]

    return data


# ============================================================
# AIRFOIL GEOMETRY
# ============================================================

def prepare_airfoil(raw_data):
    """
    Prepare an airfoil for the same ordering assumed by the MATLAB code:

        upper surface: first half
        lower surface: second half

    The code requires a closed coordinate sequence, i.e. first and last
    points should represent the trailing edge.
    """

    data = np.asarray(raw_data, dtype=float)

    # If the airfoil is not closed, close it.
    if np.linalg.norm(data[0] - data[-1]) > 1e-10:
        data = np.vstack([data, data[0]])

    L = len(data)

    # For the source formulation, L must be odd:
    # e.g. 97 points -> 96 panels.
    if (L - 1) % 2 != 0:
        raise ValueError(
            "The DAT file must have an odd number of coordinate points "
            "(including the repeated trailing-edge point), so that the "
            "upper and lower surfaces contain equal numbers of panels."
        )

    return data


def find_max_thickness(data):
    """
    Find maximum thickness using corresponding upper/lower points.

    Assumes the DAT file is ordered from TE -> upper surface -> LE ->
    lower surface -> TE, as in the supplied MATLAB program.
    """

    L = len(data)
    n = (L + 1) // 2

    upper = data[:n]
    lower = np.flip(data[n - 1:], axis=0)

    if len(upper) != len(lower):
        raise ValueError("Could not pair upper and lower surfaces.")

    thickness = np.abs(upper[:, 1] - lower[:, 1])

    index = int(np.argmax(thickness))

    x_spar = upper[index, 0]
    y_upper = upper[index, 1]
    y_lower = lower[index, 1]
    t_max = thickness[index]

    return {
        "upper": upper,
        "lower": lower,
        "index": index,
        "x_spar": x_spar,
        "y_upper": y_upper,
        "y_lower": y_lower,
        "t_max": t_max
    }


# ============================================================
# THIN-WALLED SECTION ANALYSIS
# ============================================================

def calculate_section(
    raw_data,
    chord,
    Vy,
    rib_spacing,
    skin_t,
    web_t,
    flange_t,
    flange_width
):
    """
    Calculate shear flow for an airfoil with a C-section spar.

    The spar is automatically positioned at the maximum-thickness station,
    following the supplied MATLAB model.

    The C-section is represented as:
        - vertical web of thickness web_t
        - upper and lower flanges of thickness flange_t
        - flange width flange_width

    The flanges are represented as concentrated boom areas at their
    skin/spar intersections, while the web is represented as an internal
    thin-walled spar segment.

    The external airfoil skin is the closed outer contour.
    """

    data = prepare_airfoil(raw_data)

    # Scale DAT coordinates by chord.
    data = data * chord

    L = len(data)
    n = (L + 1) // 2

    geo = find_max_thickness(data)

    # Recalculate paired surfaces after scaling.
    upper = geo["upper"] * chord
    lower = geo["lower"] * chord

    # Maximum thickness station.
    index = geo["index"]

    x_spar = upper[index, 0]
    y_u = upper[index, 1]
    y_l = lower[index, 1]

    web_height = abs(y_u - y_l)

    if web_height <= 0:
        raise ValueError("Maximum-thickness spar height is zero.")

    # --------------------------------------------------------
    # Panel lengths
    # --------------------------------------------------------

    dis = np.linalg.norm(
        np.diff(data, axis=0),
        axis=1
    )

    # --------------------------------------------------------
    # C-section areas
    # --------------------------------------------------------

    web_area = web_t * web_height

    upper_flange_area = flange_t * flange_width
    lower_flange_area = flange_t * flange_width

    # Effective spar boom areas at top and bottom.
    # Half of the web area is assigned to each flange attachment.
    B_upper = upper_flange_area + 0.5 * web_area
    B_lower = lower_flange_area + 0.5 * web_area

    # --------------------------------------------------------
    # Centroid of thin-walled skin + C spar
    # --------------------------------------------------------

    mid = 0.5 * (data[:-1] + data[1:])

    skin_area = dis * skin_t

    total_area = (
        np.sum(skin_area)
        + upper_flange_area
        + lower_flange_area
        + web_area
    )

    cg = (
        np.sum(mid * skin_area[:, None], axis=0)
        + upper_flange_area * np.array([x_spar, y_u])
        + lower_flange_area * np.array([x_spar, y_l])
        + web_area * np.array([x_spar, 0.5 * (y_u + y_l)])
    ) / total_area

    # --------------------------------------------------------
    # Section inertias
    # --------------------------------------------------------

    shifted = data - cg

    Ixy = np.sum(
        shifted[:-1, 0] *
        shifted[1:, 0] * 0.0
    )

    # Use thin-wall midpoint integration for the skin.
    mid_shifted = mid - cg

    I_x = np.sum(
        mid_shifted[:, 1] ** 2 *
        skin_area
    )

    I_y = np.sum(
        mid_shifted[:, 0] ** 2 *
        skin_area
    )

    I_xy = np.sum(
        mid_shifted[:, 0] *
        mid_shifted[:, 1] *
        skin_area
    )

    # Add C-section web using its centroid.
    web_yc = 0.5 * (y_u + y_l)
    dy_web = web_yc - cg[1]
    dx_web = x_spar - cg[0]

    # Web local I about x/y axes.
    I_x_web_local = web_t * web_height**3 / 12
    I_y_web_local = web_height * web_t**3 / 12

    I_x += I_x_web_local + web_area * dy_web**2
    I_y += I_y_web_local + web_area * dx_web**2
    I_xy += web_area * dx_web * dy_web

    # Add flanges as rectangular thin strips.
    for y_f, area in [
        (y_u, upper_flange_area),
        (y_l, lower_flange_area)
    ]:
        dx = x_spar - cg[0]
        dy = y_f - cg[1]

        I_x += area * dy**2
        I_y += area * dx**2
        I_xy += area * dx * dy

    DR = I_x * I_y - I_xy**2

    if abs(DR) < 1e-30:
        raise ValueError("Section inertia matrix is singular.")

    Kxy = I_xy / DR
    Ky = I_y / DR

    # Neutral-axis slope.
    tan_alpha = -I_xy / I_y

    # --------------------------------------------------------
    # Skin boom idealisation
    # --------------------------------------------------------

    shifted = data - cg

    p_d = np.abs(
        shifted[:, 1] -
        tan_alpha * shifted[:, 0]
    ) / np.sqrt(1 + tan_alpha**2)

    p_d = np.maximum(p_d, 1e-14)

    a_boom = np.zeros(L - 1)

    for i in range(L - 1):

        if i < L - 2:

            a_boom[i] += (
                (1 / 6)
                * dis[i]
                * skin_t
                * (2 + p_d[i + 1] / p_d[i])
            )

            a_boom[i + 1] += (
                (1 / 6)
                * dis[i]
                * skin_t
                * (2 + p_d[i] / p_d[i + 1])
            )

        else:

            a_boom[i] += (
                (1 / 6)
                * dis[i]
                * skin_t
                * (2 + p_d[0] / p_d[i])
            )

            a_boom[0] += (
                (1 / 6)
                * dis[i]
                * skin_t
                * (2 + p_d[i] / p_d[0])
            )

    # Add C-section flange/web boom contributions at the spar stations.
    upper_boom_index = index
    lower_boom_index = L - 1 - index

    a_boom[upper_boom_index] += B_upper
    a_boom[lower_boom_index] += B_lower

    # --------------------------------------------------------
    # Basic shear-flow calculation around the outer skin
    # --------------------------------------------------------

    data2 = data.T

    cg_boom = np.array([
        np.sum(data2[0, :L - 1] * a_boom) / np.sum(a_boom),
        np.sum(data2[1, :L - 1] * a_boom) / np.sum(a_boom)
    ])

    data3 = data2 - cg_boom[:, None]

    Ixy_b = np.sum(
        data3[0, :L - 1] *
        data3[1, :L - 1] *
        a_boom
    )

    Ix_b = np.sum(
        data3[1, :L - 1]**2 *
        a_boom
    )

    Iy_b = np.sum(
        data3[0, :L - 1]**2 *
        a_boom
    )

    den_b = Ix_b * Iy_b - Ixy_b**2

    if abs(den_b) < 1e-30:
        raise ValueError("Boom inertia matrix is singular.")

    Kxy_b = Ixy_b / den_b
    Ky_b = Iy_b / den_b

    Qy = np.zeros(L - 1)
    Qx = np.zeros(L - 1)
    q_basic = np.zeros(L - 1)

    # Traverse from TE over upper surface toward LE.
    for i in range(0, upper_boom_index - 1):
        Qy[i + 1] = (
            Qy[i]
            + data3[0, i + 1] * a_boom[i + 1]
        )
        Qx[i + 1] = (
            Qx[i]
            + data3[1, i + 1] * a_boom[i + 1]
        )
        q_basic[i + 1] = Vy * (
            Kxy_b * Qy[i + 1]
            - Ky_b * Qx[i + 1]
        )

    # Continue through the remaining skin stations.
    for i in range(max(upper_boom_index, 0), L - 2):
        Qy[i + 1] = (
            Qy[i]
            + data3[0, i + 1] * a_boom[i + 1]
        )
        Qx[i + 1] = (
            Qx[i]
            + data3[1, i + 1] * a_boom[i + 1]
        )
        q_basic[i + 1] = Vy * (
            Kxy_b * Qy[i + 1]
            - Ky_b * Qx[i + 1]
        )

    # --------------------------------------------------------
    # Two-cell closed-section compatibility
    #
    # Cell 1 = forward/leading-edge side of spar
    # Cell 2 = aft/trailing-edge side of spar
    #
    # The spar web is the common wall.
    # --------------------------------------------------------

    # Determine panel sets on each side of the spar.
    # The DAT ordering is TE -> upper -> LE -> lower -> TE.
    # The spar creates two closed cells by connecting upper/lower skin.
    #
    # Outer-skin panel indices:
    #   cell 1: upper path from TE to spar + lower path spar to TE
    #   cell 2: upper spar->LE + lower LE->spar

    # Upper station index is 'index'; corresponding lower point:
    lower_point = L - 1 - index

    # Panels on upper side:
    # 0 ... index-1 are TE->spar
    # index ... n-2 are spar->LE
    #
    # Lower side:
    # lower_point ... L-2 are spar->TE
    # n-1 ... lower_point-1 are LE->spar

    cell1_panels = list(range(0, index)) + list(
        range(lower_point, L - 1)
    )

    cell2_panels = list(
        range(index, n - 1)
    ) + list(
        range(n - 1, lower_point)
    )

    # Remove invalid/duplicate indices.
    cell1_panels = sorted(set(
        i for i in cell1_panels
        if 0 <= i < L - 1
    ))

    cell2_panels = sorted(set(
        i for i in cell2_panels
        if 0 <= i < L - 1
    ))

    if not cell1_panels or not cell2_panels:
        raise ValueError(
            "Could not form two cells. Check the DAT point ordering."
        )

    # Spar web length is actual maximum thickness.
    spar_length = web_height

    # Basic shear flow through the spar web is obtained from the jump
    # between the two skin-side basic flows at the spar.
    q_web_basic = (
        q_basic[lower_point]
        - q_basic[index]
    )

    # Cell areas using polygon area of each cell.
    # Approximate with line integration of the outer skin and spar.
    A1 = calculate_cell_area(
        data,
        index,
        lower_point,
        leading=False
    )

    A2 = calculate_cell_area(
        data,
        index,
        lower_point,
        leading=True
    )

    # Circumference integrals.
    J1 = np.sum(
        dis[cell1_panels] / skin_t
    )

    J2 = np.sum(
        dis[cell2_panels] / skin_t
    )

    # Add spar web to both cell compatibility equations.
    J1 += spar_length / web_t
    J2 += spar_length / web_t

    # Coupling through common spar.
    A = np.array([
        [J1, -spar_length / web_t],
        [-spar_length / web_t, J2]
    ])

    # Basic-flow compatibility RHS.
    Bc = np.array([
        -np.sum(
            q_basic[cell1_panels] *
            dis[cell1_panels] /
            skin_t
        ) - q_web_basic * spar_length / web_t,

        -np.sum(
            q_basic[cell2_panels] *
            dis[cell2_panels] /
            skin_t
        ) + q_web_basic * spar_length / web_t
    ])

    try:
        q0 = np.linalg.solve(A, Bc)
    except np.linalg.LinAlgError:
        q0 = np.zeros(2)

    # Correct outer skin flows.
    q_skin = q_basic.copy()

    q_skin[cell1_panels] += q0[0]
    q_skin[cell2_panels] += q0[1]

    # Common spar-web flow.
    q_web = q_web_basic + q0[0] - q0[1]

    # --------------------------------------------------------
    # Stringer locations are intentionally handled separately
    # in the GUI. The C spar itself remains fixed at max thickness.
    # --------------------------------------------------------

    return {
        "data": data,
        "cg": cg,
        "upper": upper,
        "lower": lower,
        "index": index,
        "x_spar": x_spar,
        "y_upper": y_u,
        "y_lower": y_l,
        "t_max": web_height,
        "web_height": web_height,
        "B_upper": B_upper,
        "B_lower": B_lower,
        "web_area": web_area,
        "upper_flange_area": upper_flange_area,
        "lower_flange_area": lower_flange_area,
        "dis": dis,
        "a_boom": a_boom,
        "q_basic": q_basic,
        "q_skin": q_skin,
        "q_web": q_web,
        "q0": q0,
        "I_x": I_x,
        "I_y": I_y,
        "I_xy": I_xy,
        "tan_alpha": tan_alpha,
        "cell1_panels": cell1_panels,
        "cell2_panels": cell2_panels,
        "A1": A1,
        "A2": A2
    }


def calculate_cell_area(data, index, lower_point, leading):
    """
    Approximate cell area bounded by the outer airfoil skin and the spar.

    leading=False -> trailing/aft cell
    leading=True  -> leading/front cell
    """

    n = (len(data) + 1) // 2

    if leading:
        # Upper spar -> LE -> lower spar
        pts = np.vstack([
            data[index],
            data[index:n],
            data[n - 1:lower_point + 1]
        ])
    else:
        # TE -> upper spar -> lower spar -> TE
        pts = np.vstack([
            data[0:index + 1],
            data[lower_point:],
            data[0:1]
        ])

    x = pts[:, 0]
    y = pts[:, 1]

    return 0.5 * abs(
        np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    )


# ============================================================
# GUI
# ============================================================

class CSparGUI:

    def __init__(self, root):

        self.root = root
        self.root.title(
            "Airfoil Shear Flow Analysis — C-Section Spar"
        )
        self.root.geometry("1400x880")
        self.root.minsize(1150, 720)

        self.data = None
        self.result = None

        self.upper_indices = set()
        self.lower_indices = set()

        self.build_gui()


    # --------------------------------------------------------
    # GUI
    # --------------------------------------------------------

    def build_gui(self):

        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        left = ttk.Frame(
            self.root,
            padding=12
        )
        left.grid(
            row=0,
            column=0,
            sticky="ns"
        )

        right = ttk.Frame(
            self.root,
            padding=(0, 12, 12, 12)
        )
        right.grid(
            row=0,
            column=1,
            sticky="nsew"
        )

        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        # ====================================================
        # FILE
        # ====================================================

        box = ttk.LabelFrame(
            left,
            text="Airfoil DAT File",
            padding=10
        )
        box.pack(fill="x", pady=(0, 10))

        self.file_var = tk.StringVar()

        ttk.Entry(
            box,
            textvariable=self.file_var,
            width=40
        ).pack(fill="x", pady=(0, 6))

        ttk.Button(
            box,
            text="Browse .dat File",
            command=self.browse_file
        ).pack(fill="x")

        self.status = ttk.Label(
            box,
            text="No file loaded"
        )
        self.status.pack(anchor="w", pady=(6, 0))

        # ====================================================
        # LOAD / AUTO SPAR
        # ====================================================

        box = ttk.LabelFrame(
            left,
            text="Spar Definition",
            padding=10
        )
        box.pack(fill="x", pady=(0, 10))

        ttk.Label(
            box,
            text=(
                "The spar is automatically placed at the "
                "maximum-thickness station."
            ),
            wraplength=340
        ).pack(anchor="w", pady=(0, 8))

        self.spar_position_var = tk.StringVar(
            value="Maximum thickness"
        )

        ttk.Label(
            box,
            textvariable=self.spar_position_var
        ).pack(anchor="w")

        # ====================================================
        # INPUTS
        # ====================================================

        box = ttk.LabelFrame(
            left,
            text="Analysis Inputs",
            padding=10
        )
        box.pack(fill="x", pady=(0, 10))

        self.chord_var = tk.StringVar(value="0.234")
        self.shear_var = tk.StringVar(value="100")
        self.rib_var = tk.StringVar(value="0.2")
        self.skin_var = tk.StringVar(value="0.0005")

        self.web_t_var = tk.StringVar(value="0.0005")
        self.flange_t_var = tk.StringVar(value="0.0005")
        self.flange_w_var = tk.StringVar(value="0.015")

        self.add_entry(box, "Chord c [m]", self.chord_var)
        self.add_entry(box, "Shear force Vy [N]", self.shear_var)
        self.add_entry(box, "Rib spacing a [m]", self.rib_var)
        self.add_entry(box, "Skin thickness [m]", self.skin_var)

        ttk.Separator(box).pack(
            fill="x",
            pady=5
        )

        ttk.Label(
            box,
            text="C-Section Spar"
        ).pack(anchor="w", pady=(0, 3))

        self.add_entry(
            box,
            "Web thickness [m]",
            self.web_t_var
        )

        self.add_entry(
            box,
            "Flange thickness [m]",
            self.flange_t_var
        )

        self.add_entry(
            box,
            "Flange width [m]",
            self.flange_w_var
        )

        # ====================================================
        # STRINGERS
        # ====================================================

        box = ttk.LabelFrame(
            left,
            text="Stringers",
            padding=10
        )
        box.pack(fill="x", pady=(0, 10))

        ttk.Label(
            box,
            text=(
                "Select Upper or Lower, then click an airfoil "
                "point. Click again to remove it."
            ),
            wraplength=340
        ).pack(anchor="w", pady=(0, 6))

        self.mode_var = tk.StringVar(
            value="Upper"
        )

        row = ttk.Frame(box)
        row.pack(fill="x")

        ttk.Radiobutton(
            row,
            text="Upper",
            variable=self.mode_var,
            value="Upper"
        ).pack(side="left", padx=(0, 15))

        ttk.Radiobutton(
            row,
            text="Lower",
            variable=self.mode_var,
            value="Lower"
        ).pack(side="left")

        self.upper_label = ttk.Label(
            box,
            text="Upper: None"
        )
        self.upper_label.pack(anchor="w", pady=(6, 0))

        self.lower_label = ttk.Label(
            box,
            text="Lower: None"
        )
        self.lower_label.pack(anchor="w")

        ttk.Button(
            box,
            text="Clear Stringers",
            command=self.clear_stringers
        ).pack(fill="x", pady=(7, 0))

        # ====================================================
        # ANALYSIS BUTTON
        # ====================================================

        ttk.Button(
            left,
            text="RUN ANALYSIS",
            command=self.run_analysis
        ).pack(fill="x", pady=(2, 6))

        ttk.Button(
            left,
            text="Save Figure",
            command=self.save_figure
        ).pack(fill="x")

        # ====================================================
        # RESULTS
        # ====================================================

        box = ttk.LabelFrame(
            left,
            text="Results",
            padding=7
        )
        box.pack(
            fill="both",
            expand=True,
            pady=(10, 0)
        )

        self.result_text = tk.Text(
            box,
            width=44,
            height=18,
            font=("DejaVu Sans Mono", 9),
            wrap="word"
        )
        self.result_text.pack(
            fill="both",
            expand=True
        )

        # ====================================================
        # FIGURE
        # ====================================================

        self.fig, self.ax = plt.subplots(
            2,
            2,
            figsize=(10, 7)
        )

        self.canvas = FigureCanvasTkAgg(
            self.fig,
            master=right
        )

        self.canvas.get_tk_widget().grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        self.canvas.mpl_connect(
            "button_press_event",
            self.on_plot_click
        )

        self.show_empty_plot()


    def add_entry(
        self,
        parent,
        label,
        variable
    ):

        ttk.Label(
            parent,
            text=label
        ).pack(anchor="w")

        ttk.Entry(
            parent,
            textvariable=variable
        ).pack(
            fill="x",
            pady=(2, 6)
        )


    # ========================================================
    # EMPTY PLOT
    # ========================================================

    def show_empty_plot(self):

        for a in self.ax.flat:
            a.clear()
            a.grid(True, alpha=0.25)

        self.ax[0, 0].set_title(
            "Load an Airfoil DAT File"
        )
        self.ax[0, 1].set_title(
            "Skin Shear Flow"
        )
        self.ax[1, 0].set_title(
            "C-Spar Shear Flow"
        )
        self.ax[1, 1].set_title(
            "Section / Spar"
        )

        self.fig.tight_layout(
            pad=3
        )
        self.canvas.draw()


    # ========================================================
    # LOAD FILE
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

            self.data = read_dat_file(
                filename
            )

            self.data = prepare_airfoil(
                self.data
            )

            self.file_var.set(
                filename
            )

            self.status.config(
                text=f"{len(self.data)} points loaded"
            )

            self.upper_indices.clear()
            self.lower_indices.clear()
            self.update_stringer_labels()

            geo = find_max_thickness(
                self.data
            )

            self.spar_position_var.set(
                f"Maximum thickness: "
                f"point {geo['index'] + 1}, "
                f"t/c = {geo['t_max']:.5f}"
            )

            self.plot_geometry()

        except Exception as e:

            messagebox.showerror(
                "DAT File Error",
                str(e)
            )


    # ========================================================
    # GEOMETRY PLOT
    # ========================================================

    def plot_geometry(self):

        if self.data is None:
            return

        geo = find_max_thickness(
            self.data
        )

        data = self.data

        ax = self.ax[0, 0]
        ax.clear()
        ax.grid(True, alpha=0.25)

        ax.plot(
            data[:, 0],
            data[:, 1],
            linewidth=1.6
        )

        ax.scatter(
            data[:, 0],
            data[:, 1],
            s=12
        )

        # Spar web
        ax.plot(
            [geo["x_spar"], geo["x_spar"]],
            [geo["y_lower"], geo["y_upper"]],
            linewidth=3
        )

        # Approximate C-section flanges.
        flange_width = float(
            self.flange_w_var.get()
        )

        ax.plot(
            [
                geo["x_spar"] - flange_width,
                geo["x_spar"]
            ],
            [
                geo["y_upper"],
                geo["y_upper"]
            ],
            linewidth=3
        )

        ax.plot(
            [
                geo["x_spar"] - flange_width,
                geo["x_spar"]
            ],
            [
                geo["y_lower"],
                geo["y_lower"]
            ],
            linewidth=3
        )

        ax.scatter(
            [geo["x_spar"], geo["x_spar"]],
            [geo["y_upper"], geo["y_lower"]],
            s=70,
            marker="o"
        )

        ax.set_title(
            "Airfoil with C-Section Spar"
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.axis("equal")

        self.plot_stringers(
            ax
        )

        self.fig.tight_layout(
            pad=3
        )
        self.canvas.draw()


    # ========================================================
    # STRINGER INTERACTION
    # ========================================================

    def on_plot_click(self, event):

        if self.data is None:
            return

        if event.inaxes != self.ax[0, 0]:
            return

        if event.xdata is None or event.ydata is None:
            return

        points = self.data

        distances = np.sqrt(
            (points[:, 0] - event.xdata)**2
            + (points[:, 1] - event.ydata)**2
        )

        index = int(
            np.argmin(distances)
        )

        L = len(points)
        half = int(
            round(0.5 * (L - 1))
        )

        if index >= L - 1:
            return

        if self.mode_var.get() == "Upper":

            if index >= half:
                return

            p = index + 1

            if p in self.upper_indices:
                self.upper_indices.remove(p)
            else:
                self.upper_indices.add(p)

        else:

            if index < half:
                return

            p = index - half + 1

            if p in self.lower_indices:
                self.lower_indices.remove(p)
            else:
                self.lower_indices.add(p)

        self.update_stringer_labels()
        self.plot_geometry()


    def plot_stringers(self, ax):

        if self.data is None:
            return

        geo = find_max_thickness(
            self.data
        )

        half = int(
            round(0.5 * (len(self.data) - 1))
        )

        for p in sorted(
            self.upper_indices
        ):

            j = p - 1

            if 0 <= j < half:

                ax.scatter(
                    self.data[j, 0],
                    self.data[j, 1],
                    s=80,
                    marker="o"
                )

                ax.annotate(
                    str(p),
                    (
                        self.data[j, 0],
                        self.data[j, 1]
                    ),
                    xytext=(5, 5),
                    textcoords="offset points"
                )

        for p in sorted(
            self.lower_indices
        ):

            j = half + p - 1

            if 0 <= j < len(self.data) - 1:

                ax.scatter(
                    self.data[j, 0],
                    self.data[j, 1],
                    s=80,
                    marker="s"
                )

                ax.annotate(
                    str(p),
                    (
                        self.data[j, 0],
                        self.data[j, 1]
                    ),
                    xytext=(5, 5),
                    textcoords="offset points"
                )


    def update_stringer_labels(self):

        self.upper_label.config(
            text=(
                "Upper: "
                + (
                    "None"
                    if not self.upper_indices
                    else " ".join(
                        str(i)
                        for i in sorted(
                            self.upper_indices
                        )
                    )
                )
            )
        )

        self.lower_label.config(
            text=(
                "Lower: "
                + (
                    "None"
                    if not self.lower_indices
                    else " ".join(
                        str(i)
                        for i in sorted(
                            self.lower_indices
                        )
                    )
                )
            )
        )


    def clear_stringers(self):

        self.upper_indices.clear()
        self.lower_indices.clear()

        self.update_stringer_labels()

        if self.data is not None:
            self.plot_geometry()


    # ========================================================
    # ANALYSIS
    # ========================================================

    def run_analysis(self):

        try:

            if self.data is None:
                raise ValueError(
                    "Load an airfoil DAT file first."
                )

            chord = float(
                self.chord_var.get()
            )

            Vy = float(
                self.shear_var.get()
            )

            rib = float(
                self.rib_var.get()
            )

            skin_t = float(
                self.skin_var.get()
            )

            web_t = float(
                self.web_t_var.get()
            )

            flange_t = float(
                self.flange_t_var.get()
            )

            flange_w = float(
                self.flange_w_var.get()
            )

            values = [
                chord,
                rib,
                skin_t,
                web_t,
                flange_t,
                flange_w
            ]

            if any(v <= 0 for v in values):
                raise ValueError(
                    "All geometry and spacing inputs must be positive."
                )

            self.result = calculate_section(
                self.data,
                chord,
                Vy,
                rib,
                skin_t,
                web_t,
                flange_t,
                flange_w
            )

            self.display_results(
                chord,
                Vy,
                rib,
                skin_t,
                web_t,
                flange_t,
                flange_w
            )

            self.display_analysis_plots()

        except Exception as e:

            messagebox.showerror(
                "Analysis Error",
                str(e)
            )


    # ========================================================
    # RESULTS
    # ========================================================

    def display_results(
        self,
        chord,
        Vy,
        rib,
        skin_t,
        web_t,
        flange_t,
        flange_w
    ):

        r = self.result

        q_skin_max = np.max(
            np.abs(r["q_skin"])
        )

        self.result_text.delete(
            "1.0",
            tk.END
        )

        text = f"""
C-SECTION SPAR SHEAR-FLOW ANALYSIS
============================================

INPUTS
--------------------------------------------
Chord                 : {chord:.6g} m
Shear force Vy        : {Vy:.6g} N
Rib spacing           : {rib:.6g} m
Skin thickness        : {skin_t:.6g} m

C-SPAR
--------------------------------------------
Web thickness         : {web_t:.6g} m
Flange thickness      : {flange_t:.6g} m
Flange width          : {flange_w:.6g} m
Web height            : {r["web_height"]:.6g} m

Spar x-location       : {r["x_spar"]:.6g} m
Maximum thickness     : {r["t_max"]:.6g} m

SECTION PROPERTIES
--------------------------------------------
Centroid x            : {r["cg"][0]:.8e} m
Centroid y            : {r["cg"][1]:.8e} m
I_x                   : {r["I_x"]:.8e} m^4
I_y                   : {r["I_y"]:.8e} m^4
I_xy                  : {r["I_xy"]:.8e} m^4
tan(alpha)            : {r["tan_alpha"]:.8e}

SPAR AREAS
--------------------------------------------
Web area              : {r["web_area"]:.8e} m²
Upper flange area     : {r["upper_flange_area"]:.8e} m²
Lower flange area     : {r["lower_flange_area"]:.8e} m²

SHEAR FLOW
--------------------------------------------
Maximum skin |q|      : {q_skin_max:.8e} N/m
C-spar web q          : {r["q_web"]:.8e} N/m

Cell 1 q0             : {r["q0"][0]:.8e} N/m
Cell 2 q0             : {r["q0"][1]:.8e} N/m

STRINGERS
--------------------------------------------
Upper indices         : {sorted(self.upper_indices) or "None"}
Lower indices         : {sorted(self.lower_indices) or "None"}

NOTE:
The C-section spar is fixed at the maximum-thickness
station, following the supplied MATLAB model.
"""

        self.result_text.insert(
            "1.0",
            text
        )


    # ========================================================
    # ANALYSIS PLOTS
    # ========================================================

    def display_analysis_plots(self):

        r = self.result

        data = r["data"]
        q = r["q_skin"]

        # Geometry
        self.plot_geometry()

        # Skin shear flow
        ax = self.ax[0, 1]
        ax.clear()
        ax.grid(True, alpha=0.25)

        half = int(
            round(0.5 * (len(data) - 1))
        )

        ax.plot(
            data[:half, 0],
            -q[:half],
            linewidth=1.8,
            label="Upper"
        )

        ax.plot(
            data[half:-1, 0],
            -q[half:],
            linewidth=1.8,
            label="Lower"
        )

        ax.set_title(
            "Shear Flow in Airfoil Skin"
        )
        ax.set_xlabel(
            "x [m]"
        )
        ax.set_ylabel(
            "q [N/m]"
        )
        ax.legend()

        # C-spar shear flow
        ax = self.ax[1, 0]
        ax.clear()
        ax.grid(True, alpha=0.25)

        geo = find_max_thickness(
            data
        )

        ax.plot(
            [
                r["x_spar"],
                r["x_spar"]
            ],
            [
                r["y_lower"],
                r["y_upper"]
            ],
            linewidth=5
        )

        ax.scatter(
            [r["x_spar"]],
            [0.5 * (
                r["y_lower"] +
                r["y_upper"]
            )],
            s=80
        )

        ax.set_title(
            f"C-Spar Web Shear Flow = "
            f"{r['q_web']:.3e} N/m"
        )
        ax.set_xlabel(
            "x [m]"
        )
        ax.set_ylabel(
            "y [m]"
        )
        ax.axis("equal")

        # Section overview
        ax = self.ax[1, 1]
        ax.clear()
        ax.grid(True, alpha=0.25)

        ax.plot(
            data[:, 0],
            data[:, 1],
            linewidth=1.5
        )

        ax.plot(
            [r["x_spar"], r["x_spar"]],
            [r["y_lower"], r["y_upper"]],
            linewidth=4
        )

        flange_w = float(
            self.flange_w_var.get()
        )

        ax.plot(
            [
                r["x_spar"] - flange_w,
                r["x_spar"]
            ],
            [
                r["y_upper"],
                r["y_upper"]
            ],
            linewidth=4
        )

        ax.plot(
            [
                r["x_spar"] - flange_w,
                r["x_spar"]
            ],
            [
                r["y_lower"],
                r["y_lower"]
            ],
            linewidth=4
        )

        ax.set_title(
            "Airfoil + C-Section Spar"
        )
        ax.set_xlabel(
            "x [m]"
        )
        ax.set_ylabel(
            "y [m]"
        )
        ax.axis("equal")

        self.fig.tight_layout(
            pad=3
        )

        self.canvas.draw()


    # ========================================================
    # SAVE
    # ========================================================

    def save_figure(self):

        if self.result is None:

            messagebox.showwarning(
                "No Results",
                "Run the analysis first."
            )

            return

        filename = filedialog.asksaveasfilename(
            title="Save Figure",
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

    app = CSparGUI(root)

    root.mainloop()
