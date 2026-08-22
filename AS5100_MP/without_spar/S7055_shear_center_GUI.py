import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

plt.rc("font", size=12)


# ================================================================
# GLOBAL AIRFOIL DATA
# ================================================================

airfoil = None
airfoil_file = None


# ================================================================
# LOAD AIRFOIL
# ================================================================

def load_airfoil():
    global airfoil, airfoil_file

    filename = filedialog.askopenfilename(
        title="Select airfoil coordinate file",
        filetypes=[
            ("DAT / TXT files", "*.dat *.txt"),
            ("All files", "*.*")
        ]
    )

    if not filename:
        return

    try:
        data = np.loadtxt(filename)

        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(
                "The selected file must contain at least two columns: x and y."
            )

        airfoil = data[:, :2]
        airfoil_file = filename

        file_label.config(
            text=f"Loaded: {filename}"
        )

        status_label.config(
            text=f"Loaded {len(airfoil)} coordinate points."
        )

    except Exception as e:
        messagebox.showerror(
            "File error",
            f"Could not load the airfoil file.\n\n{e}"
        )


# ================================================================
# PARSE INPUT ARRAY
# ================================================================

def parse_array(text):

    text = text.strip()

    if not text:
        return np.array([])

    try:
        return np.array([
            float(value.strip())
            for value in text.split(",")
            if value.strip()
        ])
    except ValueError:
        raise ValueError(
            "Stringer positions must be comma-separated numbers.\n"
            "Example: 0.10, 0.30, 0.50, 0.70, 0.90"
        )


# ================================================================
# CALCULATION
# ================================================================

def calculate():

    global airfoil

    if airfoil is None:
        messagebox.showwarning(
            "No airfoil loaded",
            "Load an airfoil coordinate file first."
        )
        return

    try:

        chord = float(chord_var.get())
        skin_t = float(skin_t_var.get())
        stringer_t = float(stringer_t_var.get())
        stringer_L = float(stringer_L_var.get())
        spar_xc = float(spar_xc_var.get())
        spar_t = float(spar_t_var.get())
        spar_flange_L = float(spar_flange_L_var.get())

        V = float(V_var.get())

        Xref = float(Xref_var.get())
        Cmref = float(Cmref_var.get())
        CL = float(CL_var.get())

        upper_xc = parse_array(upper_var.get())
        lower_xc = parse_array(lower_var.get())

        if not 0.0 <= spar_xc <= 1.0:
            raise ValueError("Spar x/c must be within 0 <= x/c <= 1.")

        if spar_t < 0.0 or spar_flange_L < 0.0:
            raise ValueError("Spar thickness and flange length must be non-negative.")

        # --------------------------------------------------------
        # AIRFOIL GEOMETRY
        # --------------------------------------------------------

        geom = airfoil * chord

        # Remove duplicate trailing-edge point
        if np.allclose(geom[0], geom[-1]):
            geom = geom[:-1]

        x = geom[:, 0]
        y = geom[:, 1]

        n = len(x)

        x_next = np.roll(x, -1)
        y_next = np.roll(y, -1)

        dx = x_next - x
        dy = y_next - y

        ds = np.hypot(dx, dy)

        x_mid = 0.5 * (x + x_next)
        y_mid = 0.5 * (y + y_next)

        # --------------------------------------------------------
        # SKIN BOOM IDEALISATION
        # --------------------------------------------------------

        A_panel = skin_t * ds
        A_skin = np.sum(A_panel)

        x_skin_centroid = (
            np.sum(A_panel * x_mid) / A_skin
        )

        y_skin_centroid = (
            np.sum(A_panel * y_mid) / A_skin
        )

        y_skin_c = y - y_skin_centroid

        if np.any(np.abs(y_skin_c) < 1e-10):

            raise ValueError(
                "A skin boom is essentially on the skin centroidal "
                "axis. The standard boom-area formula becomes "
                "singular."
            )

        B_skin = (
            skin_t / 6.0 *
            (
                np.roll(ds, 1) *
                (
                    2.0 +
                    np.roll(y_skin_c, 1) /
                    y_skin_c
                )
                +
                ds *
                (
                    2.0 +
                    np.roll(y_skin_c, -1) /
                    y_skin_c
                )
            )
        )

        # --------------------------------------------------------
        # STRINGER BOOMS
        # --------------------------------------------------------

        B_stringer = np.zeros(n)

        A_stringer = (
            stringer_t *
            stringer_L
        )

        upper_indices = []
        lower_indices = []

        def nearest_surface_index(xc, surface):

            if not 0.0 <= xc <= 1.0:

                raise ValueError(
                    f"{surface.capitalize()} stringer x/c = "
                    f"{xc} is outside 0 <= x/c <= 1."
                )

            target_x = xc * chord

            if surface == "upper":
                candidates = np.where(y >= 0)[0]
            else:
                candidates = np.where(y <= 0)[0]

            if len(candidates) == 0:

                raise ValueError(
                    f"Could not identify the {surface} surface."
                )

            return candidates[
                np.argmin(
                    np.abs(
                        x[candidates] -
                        target_x
                    )
                )
            ]

        for xc in upper_xc:

            i = nearest_surface_index(
                xc,
                "upper"
            )

            B_stringer[i] += A_stringer

            upper_indices.append(i)

        for xc in lower_xc:

            i = nearest_surface_index(
                xc,
                "lower"
            )

            B_stringer[i] += A_stringer

            lower_indices.append(i)

        # --------------------------------------------------------
        # SPAR BOOMS
        #
        # The upper and lower spar flanges are represented by two
        # discrete booms at the spar/skin intersections. The spar web
        # carries shear and is therefore not included as a bending boom.
        # --------------------------------------------------------

        B_spar = np.zeros(n)
        A_spar_boom = spar_t * spar_flange_L

        upper_spar_index = nearest_surface_index(spar_xc, "upper")
        lower_spar_index = nearest_surface_index(spar_xc, "lower")

        B_spar[upper_spar_index] += A_spar_boom
        B_spar[lower_spar_index] += A_spar_boom

        # --------------------------------------------------------
        # TOTAL BOOM AREA
        # --------------------------------------------------------

        B_total = (
            B_skin +
            B_stringer +
            B_spar
        )

        # --------------------------------------------------------
        # FINAL CENTROID
        # --------------------------------------------------------

        A_total = np.sum(B_total)

        x_bar = (
            np.sum(B_total * x) /
            A_total
        )

        y_bar = (
            np.sum(B_total * y) /
            A_total
        )

        x_c = x - x_bar
        y_c = y - y_bar

        # --------------------------------------------------------
        # SECOND MOMENTS
        # --------------------------------------------------------

        Ixx = np.sum(
            B_total * y_c**2
        )

        Iyy = np.sum(
            B_total * x_c**2
        )

        Ixy = np.sum(
            B_total * x_c * y_c
        )

        D = (
            Ixx * Iyy -
            Ixy**2
        )

        if D <= 0:

            raise ValueError(
                "Invalid section properties: "
                "Ixx*Iyy - Ixy^2 <= 0."
            )

        # --------------------------------------------------------
        # CENTRE OF PRESSURE
        # --------------------------------------------------------

        Xcp = (
            Xref -
            Cmref / CL
        )

        x_cp = Xcp * chord

        # --------------------------------------------------------
        # BASIC SHEAR FLOW
        # --------------------------------------------------------

        dq_b = (
            V / D *
            (
                Ixy * x_c -
                Ixx * y_c
            ) *
            B_total
        )

        q_b = np.zeros(n)

        for i in range(1, n):

            q_b[i] = (
                q_b[i - 1] +
                dq_b[i - 1]
            )

        # Panel values
        q_b_panel = 0.5 * (
            q_b +
            np.roll(q_b, -1)
        )

        # --------------------------------------------------------
        # REDUNDANT SHEAR FLOW
        # --------------------------------------------------------

        q_0 = (
            -np.sum(
                q_b_panel * ds
            ) /
            np.sum(ds)
        )

        q_shear = (
            q_b_panel +
            q_0
        )

        # --------------------------------------------------------
        # MEDIAN-LINE AREA
        # --------------------------------------------------------

        x_m = 0.5 * (
            x_c +
            np.roll(x_c, -1)
        )

        y_m = 0.5 * (
            y_c +
            np.roll(y_c, -1)
        )

        A_m = 0.5 * abs(
            np.sum(
                x_m *
                np.roll(y_m, -1)
                -
                y_m *
                np.roll(x_m, -1)
            )
        )

        # --------------------------------------------------------
        # MOMENT CAUSED BY LIFT AT CP
        # --------------------------------------------------------

        M_cp = (
            V *
            (x_cp - x_bar)
        )

        # Moment-induced circulating flow
        q_moment = (
            M_cp /
            (2.0 * A_m)
        )

        # --------------------------------------------------------
        # TOTAL FLOW
        # --------------------------------------------------------

        q_total = (
            q_shear +
            q_moment
        )

        # --------------------------------------------------------
        # SHEAR CENTER
        #
        # Use the shear-force component only.
        # --------------------------------------------------------

        panel_arm = (
            x_c *
            np.roll(y_c, -1)
            -
            y_c *
            np.roll(x_c, -1)
        )

        T_shear = np.sum(
            q_shear *
            panel_arm
        )

        if abs(V) > 1e-15:

            x_sc_offset = (
                T_shear /
                V
            )

        else:

            x_sc_offset = np.nan

        x_sc = (
            x_bar +
            x_sc_offset
        )

        # Actual twisting moment because
        # CP does not lie at the shear center
        M_twist = (
            V *
            (x_cp - x_sc)
        )

        # --------------------------------------------------------
        # SAVE RESULTS
        # --------------------------------------------------------

        result = {
            "x": x,
            "y": y,
            "x_mid": x_mid,
            "y_mid": y_mid,
            "q_shear": q_shear,
            "q_moment": q_moment,
            "q_total": q_total,
            "x_bar": x_bar,
            "y_bar": y_bar,
            "Ixx": Ixx,
            "Iyy": Iyy,
            "Ixy": Ixy,
            "A_m": A_m,
            "q_0": q_0,
            "Xcp": Xcp,
            "x_cp": x_cp,
            "M_cp": M_cp,
            "x_sc": x_sc,
            "x_sc_offset": x_sc_offset,
            "M_twist": M_twist,
            "T_shear": T_shear,
            "chord": chord,
            "upper_indices": upper_indices,
            "lower_indices": lower_indices,
            "upper_spar_index": upper_spar_index,
            "lower_spar_index": lower_spar_index,
            "A_spar_boom": A_spar_boom
        }

        root.result = result

        # --------------------------------------------------------
        # DISPLAY RESULTS
        # --------------------------------------------------------

        output.delete(
            "1.0",
            tk.END
        )

        output_text = f"""
S7055 SHEAR CENTER / SHEAR FLOW
============================================================

LOAD
------------------------------------------------------------
Actual vertical shear V       = {V:.6g} N

CENTRE OF PRESSURE
------------------------------------------------------------
Xcp                           = {Xcp:.6f} c
x_cp                          = {x_cp:.6e} m

FINAL CENTROID
------------------------------------------------------------
x_bar                         = {x_bar:.6e} m
y_bar                         = {y_bar:.6e} m

SECTION PROPERTIES
------------------------------------------------------------
Ixx                           = {Ixx:.6e} m^4
Iyy                           = {Iyy:.6e} m^4
Ixy                           = {Ixy:.6e} m^4

SPAR BOOMS
------------------------------------------------------------
Spar x/c                      = {spar_xc:.6f}
Area per spar boom            = {A_spar_boom:.6e} m^2

MEDIAN-LINE CELL AREA
------------------------------------------------------------
A_m                           = {A_m:.6e} m^2

SHEAR-FORCE SHEAR FLOW
------------------------------------------------------------
q_0                           = {q_0:.6e} N/m
max |q_shear|                 = {np.max(np.abs(q_shear)):.6e} N/m

CP-INDUCED MOMENT
------------------------------------------------------------
M_CP                          = {M_cp:.6e} N m
q_moment                      = {q_moment:.6e} N/m

TOTAL SHEAR FLOW
------------------------------------------------------------
max |q_total|                 = {np.max(np.abs(q_total)):.6e} N/m
min q_total                   = {np.min(q_total):.6e} N/m
max q_total                   = {np.max(q_total):.6e} N/m

SHEAR CENTER
------------------------------------------------------------
offset from centroid          = {x_sc_offset:.6e} m
x_SC from LE                  = {x_sc:.6e} m
x_SC/c                        = {x_sc/chord:.6f}

TORSIONAL EFFECT
------------------------------------------------------------
M_twist = V(x_CP - x_SC)      = {M_twist:.6e} N m
"""

        output.insert(
            tk.END,
            output_text
        )

        status_label.config(
            text="Calculation completed."
        )

    except Exception as e:

        messagebox.showerror(
            "Calculation error",
            str(e)
        )

        status_label.config(
            text="Calculation failed."


        )


# ================================================================
# PLOT
# ================================================================

def plot_results():

    if not hasattr(root, "result"):

        messagebox.showwarning(
            "No results",
            "Run Calculate first."
        )

        return

    r = root.result

    chord = r["chord"]

    # ------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(10, 5)
    )

    ax.plot(
        r["x"] / chord,
        r["y"] / chord,
        linewidth=1.5,
        label="S7055"
    )

    ax.scatter(
        r["x_bar"] / chord,
        r["y_bar"] / chord,
        s=70,
        label="Centroid"
    )

    ax.scatter(
        r["x_sc"] / chord,
        r["y_bar"] / chord,
        marker="x",
        s=100,
        label="Shear center"
    )

    ax.scatter(
        r["Xcp"],
        r["y_bar"] / chord,
        marker="^",
        s=80,
        label="Center of pressure"
    )

    # Stringers
    upper = r["upper_indices"]
    lower = r["lower_indices"]

    if len(upper):

        ax.scatter(
            r["x"][upper] / chord,
            r["y"][upper] / chord,
            marker="s",
            s=45,
            label="Upper stringers"
        )

    if len(lower):

        ax.scatter(
            r["x"][lower] / chord,
            r["y"][lower] / chord,
            marker="s",
            s=45,
            label="Lower stringers"
        )

    ax.set_aspect(
        "equal"
    )

    ax.set_xlabel(
        "x/c"
    )

    ax.set_ylabel(
        "y/c"
    )

    ax.set_title(
        "S7055 — CP, Centroid, Stringers and Shear Center"
    )

    ax.grid()
    ax.legend()

    # ------------------------------------------------------------
    # Shear flow
    # ------------------------------------------------------------

    fig2, ax2 = plt.subplots(
        figsize=(10, 5)
    )

    ax2.plot(
        r["x_mid"] / chord,
        r["q_total"],
        label="Total q"
    )

    ax2.plot(
        r["x_mid"] / chord,
        r["q_shear"],
        linestyle="--",
        label="Shear-force q"
    )

    ax2.axhline(
        r["q_moment"],
        linestyle=":",
        label="Moment q"
    )

    ax2.set_xlabel(
        "Panel x/c"
    )

    ax2.set_ylabel(
        "Shear flow q [N/m]"
    )

    ax2.set_title(
        "S7055 Shear Flow Distribution"
    )

    ax2.grid()
    ax2.legend()

    plt.show()


# ================================================================
# GUI
# ================================================================

root = tk.Tk()

root.title(
    "S7055 Shear Center / Shear Flow"
)

root.geometry(
    "1100x800"
)

main = ttk.Frame(
    root,
    padding=12
)

main.pack(
    fill="both",
    expand=True
)

# ---------------------------------------------------------------
# File loading
# ---------------------------------------------------------------

file_frame = ttk.LabelFrame(
    main,
    text="Airfoil File",
    padding=10
)

file_frame.pack(
    fill="x"
)

ttk.Button(
    file_frame,
    text="Load Airfoil .dat / .txt",
    command=load_airfoil
).pack(
    side="left",
    padx=5
)

file_label = ttk.Label(
    file_frame,
    text="No airfoil loaded."
)

file_label.pack(
    side="left",
    padx=10
)

# ---------------------------------------------------------------
# Variables
# ---------------------------------------------------------------

chord_var = tk.StringVar(
    value="0.265"
)

skin_t_var = tk.StringVar(
    value="0.0005"
)

stringer_t_var = tk.StringVar(
    value="0.0005"
)

stringer_L_var = tk.StringVar(
    value="0.020"
)

upper_var = tk.StringVar(
    value=""
)

lower_var = tk.StringVar(
    value=""
)

V_var = tk.StringVar(
    value="-20"
)

Xref_var = tk.StringVar(
    value="0.25"
)

Cmref_var = tk.StringVar(
    value="-0.073"
)

CL_var = tk.StringVar(
    value="0.72"
)

# ---------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------

input_frame = ttk.LabelFrame(
    main,
    text="Section and Load Inputs",
    padding=10
)

input_frame.pack(
    fill="x",
    pady=10
)

fields = [
    ("Chord [m]", chord_var),
    ("Skin thickness [m]", skin_t_var),
    ("Stringer thickness [m]", stringer_t_var),
    ("Stringer length [m]", stringer_L_var),
    ("Actual shear V [N]", V_var),
    ("Xref [x/c]", Xref_var),
    ("Cmref", Cmref_var),
    ("CL", CL_var)
]

for i, (label, variable) in enumerate(fields):

    row = i // 2
    col = (i % 2) * 2

    ttk.Label(
        input_frame,
        text=label
    ).grid(
        row=row,
        column=col,
        sticky="w",
        padx=5,
        pady=4
    )

    ttk.Entry(
        input_frame,
        textvariable=variable,
        width=20
    ).grid(
        row=row,
        column=col + 1,
        padx=5,
        pady=4
    )

# ---------------------------------------------------------------
# Stringers
# ---------------------------------------------------------------

stringer_frame = ttk.LabelFrame(
    main,
    text="Stringer Positions — normalized x/c",
    padding=10
)

stringer_frame.pack(
    fill="x",
    pady=5
)

ttk.Label(
    stringer_frame,
    text="Upper stringers x/c:"
).grid(
    row=0,
    column=0,
    sticky="w",
    padx=5,
    pady=4
)

ttk.Entry(
    stringer_frame,
    textvariable=upper_var,
    width=65
).grid(
    row=0,
    column=1,
    padx=5,
    pady=4
)

ttk.Label(
    stringer_frame,
    text="Lower stringers x/c:"
).grid(
    row=1,
    column=0,
    sticky="w",
    padx=5,
    pady=4
)

ttk.Entry(
    stringer_frame,
    textvariable=lower_var,
    width=65
).grid(
    row=1,
    column=1,
    padx=5,
    pady=4
)

ttk.Label(
    stringer_frame,
    text="Example: 0.10, 0.30, 0.50, 0.70, 0.90"
).grid(
    row=2,
    column=1,
    sticky="w",
    padx=5
)

# ---------------------------------------------------------------
# Results
# ---------------------------------------------------------------

result_frame = ttk.LabelFrame(
    main,
    text="Results",
    padding=8
)

result_frame.pack(
    fill="both",
    expand=True,
    pady=10
)

output = tk.Text(
    result_frame,
    font=("DejaVu Sans Mono", 10),
    height=22
)

output.pack(
    fill="both",
    expand=True
)

# ---------------------------------------------------------------
# Buttons
# ---------------------------------------------------------------

button_frame = ttk.Frame(
    main
)

button_frame.pack(
    fill="x",
    pady=5
)

ttk.Button(
    button_frame,
    text="Calculate",
    command=calculate
).pack(
    side="left",
    padx=5
)

ttk.Button(
    button_frame,
    text="Plot Results",
    command=plot_results
).pack(
    side="left",
    padx=5
)

ttk.Button(
    button_frame,
    text="Clear",
    command=lambda: output.delete(
        "1.0",
        tk.END
    )
).pack(
    side="left",
    padx=5
)

status_label = ttk.Label(
    main,
    text="Load an airfoil coordinate file."
)

status_label.pack(
    anchor="w"
)


# ================================================================
# START GUI
# ================================================================

root.mainloop()
