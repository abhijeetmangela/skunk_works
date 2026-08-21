#!/usr/bin/env python3
"""
Spanwise Load Distribution Generator
====================================

Based on the user's original wing-load calculation.

Outputs:
    span_load.csv
    span_load.txt
    span_load_plot.png

The CSV is directly usable by the wing sheet-flow GUI.

CSV columns:
    y_m, V_N, T_Nm

Here:
    y_m  = spanwise location from root
    V_N  = internal shear force at that location
    T_Nm = torsional moment about the shear centre

IMPORTANT:
The original code calculates lift distribution, shear force and bending moment,
but does not provide the aerodynamic centre / center-of-pressure offset needed
to calculate torsion. Therefore T_Nm is not invented here.

Set LOAD_OFFSET_M below if you know the distance between the lift resultant
and the shear centre:
    T(y) = V(y) * LOAD_OFFSET_M

Set LOAD_OFFSET_M = 0.0 if you want the output to contain zero torsion.

The output shear force includes the specified load factor.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# USER INPUTS
# ---------------------------------------------------------------------

# Wing parameters
b = 2.0                       # full span [m]
S = 0.53                      # wing area [m^2]
W = 6.9 * 9.81                # aircraft weight [N]

# Discretisation
N = 1000

# Load factor
LOAD_FACTOR = 1.0

# Torsional load:
# distance between aerodynamic load line and shear centre [m]
# Example: 0.02 means T(y) = V(y)*0.02
LOAD_OFFSET_M = 0.0

# Input/output files
AIRFOIL_FILE = "s7055.dat"
OUTPUT_CSV = "span_load.csv"
OUTPUT_TXT = "span_load.txt"
OUTPUT_PLOT = "span_load_distribution.png"


# ---------------------------------------------------------------------
# LOAD AIRFOIL
# ---------------------------------------------------------------------

def load_airfoil(path: str, chord: float) -> np.ndarray:
    """
    Load a two-column airfoil .dat file and scale coordinates by chord.

    Assumes coordinates are already normalized approximately to x/c, y/c.
    """
    data = np.loadtxt(path)

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"{path} must contain at least two columns: x/c and y/c."
        )

    return data[:, :2] * chord


# ---------------------------------------------------------------------
# SPANWISE LOAD CALCULATION
# ---------------------------------------------------------------------

def calculate_spanwise_loads():
    # Spanwise coordinate, half wing
    y = np.linspace(0.0, b / 2.0, N)

    # Rectangular wing chord
    c_rectangular = S / b

    # Elliptical wing chord distribution
    inside = np.clip(1.0 - (2.0 * y / b) ** 2, 0.0, None)
    c_elliptical = (
        4.0 * S * np.sqrt(inside) / (np.pi * b)
    )

    # Schrenk equivalent chord
    c_schrenk = 0.5 * (c_rectangular + c_elliptical)

    # Scale so that the lift over the half wing equals W/2
    k = (W / 2.0) / np.trapezoid(c_schrenk, y)

    # Distributed lift [N/m]
    q_schrenk = k * c_schrenk
    q_elliptical = k * c_elliptical

    # -----------------------------------------------------------------
    # Shear force and bending moment
    #
    # V(y) = - integral_y^(b/2) l(eta) d eta
    #
    # M(y) =   integral_y^(b/2) l(eta) [eta-y] d eta
    # -----------------------------------------------------------------

    V = np.zeros_like(y)
    M = np.zeros_like(y)

    for i in range(N):
        yi = y[i]

        y_tail = y[i:]
        l_tail = q_schrenk[i:]

        V[i] = -np.trapezoid(l_tail, y_tail)

        M[i] = np.trapezoid(
            l_tail * (y_tail - yi),
            y_tail
        )

    # Apply load factor to the design loads
    V_design = LOAD_FACTOR * V
    M_design = LOAD_FACTOR * M

    # Torsional moment about the shear centre, if the load offset is known
    T_design = V_design * LOAD_OFFSET_M

    return (
        y,
        c_schrenk,
        q_schrenk,
        q_elliptical,
        V_design,
        M_design,
        T_design,
    )


# ---------------------------------------------------------------------
# WRITE OUTPUT
# ---------------------------------------------------------------------

def save_span_load_file(
    y: np.ndarray,
    V: np.ndarray,
    T: np.ndarray,
    output_csv: str,
    output_txt: str,
):
    # CSV: directly readable by the sheet-flow GUI
    data = np.column_stack((y, V, T))

    header = "y_m,V_N,T_Nm"

    np.savetxt(
        output_csv,
        data,
        delimiter=",",
        header=header,
        comments="",
        fmt="%.10e",
    )

    # Human-readable text version
    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("# Spanwise wing load distribution\n")
        f.write("# Columns: y[m], V[N], T[Nm]\n")
        f.write("y[m]        V[N]        T[Nm]\n")

        for yi, Vi, Ti in zip(y, V, T):
            f.write(
                f"{yi:10.6f}  "
                f"{Vi:12.6f}  "
                f"{Ti:12.6f}\n"
            )


# ---------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------

def make_plot(
    y: np.ndarray,
    q_schrenk: np.ndarray,
    q_elliptical: np.ndarray,
    V: np.ndarray,
    M: np.ndarray,
    T: np.ndarray,
    output_path: str,
):
    fig, axes = plt.subplots(3, 1, figsize=(10, 14))

    # Lift distribution
    axes[0].plot(
        y,
        q_schrenk,
        label="Schrenk approximated",
        linewidth=2,
    )
    axes[0].plot(
        y,
        q_elliptical,
        label="Elliptical",
        linewidth=2,
    )
    axes[0].set_title("Half-Span Lift Distribution")
    axes[0].set_xlabel("Span y [m]")
    axes[0].set_ylabel("Distributed Load [N/m]")
    axes[0].grid(True)
    axes[0].legend()

    # Shear force
    axes[1].plot(
        y,
        V,
        linewidth=2,
        label=f"Max |V| = {np.max(np.abs(V)):.4g} N",
    )
    axes[1].set_title("Half-Span Shear Force Distribution")
    axes[1].set_xlabel("Span y [m]")
    axes[1].set_ylabel("Shear Force V [N]")
    axes[1].grid(True)
    axes[1].legend()

    # Bending moment
    axes[2].plot(
        y,
        M,
        linewidth=2,
        label=f"Max |M| = {np.max(np.abs(M)):.4g} N-m",
    )
    axes[2].set_title("Half-Span Bending Moment Distribution")
    axes[2].set_xlabel("Span y [m]")
    axes[2].set_ylabel("Bending Moment M [N-m]")
    axes[2].grid(True)
    axes[2].legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.show()


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():
    # Original wing chord for the airfoil model.
    chord = S / b

    # Load airfoil, although the airfoil coordinates are not needed for the
    # spanwise load calculation itself. This preserves the original workflow
    # and checks that the supplied .dat file is valid.
    airfoil_path = Path(AIRFOIL_FILE)

    if not airfoil_path.exists():
        print(
            f"WARNING: {AIRFOIL_FILE} was not found.\n"
            "The span-load calculation can still be performed because it "
            "does not mathematically require the airfoil coordinates."
        )
    else:
        airfoil = load_airfoil(AIRFOIL_FILE, chord)
        print(
            f"Loaded airfoil: {AIRFOIL_FILE} "
            f"({len(airfoil)} coordinate points)"
        )

    (
        y,
        c_schrenk,
        q_schrenk,
        q_elliptical,
        V,
        M,
        T,
    ) = calculate_spanwise_loads()

    # Save the file for the sheet-flow GUI
    save_span_load_file(
        y=y,
        V=V,
        T=T,
        output_csv=OUTPUT_CSV,
        output_txt=OUTPUT_TXT,
    )

    make_plot(
        y=y,
        q_schrenk=q_schrenk,
        q_elliptical=q_elliptical,
        V=V,
        M=M,
        T=T,
        output_path=OUTPUT_PLOT,
    )

    # Console summary
    print("\n" + "=" * 70)
    print("SPANWISE LOAD CALCULATION")
    print("=" * 70)
    print(f"Full span                  = {b:.4f} m")
    print(f"Half span                  = {b/2:.4f} m")
    print(f"Wing area                  = {S:.4f} m^2")
    print(f"Aircraft weight            = {W:.4f} N")
    print(f"Schrenk root chord         = {c_schrenk[0]:.6f} m")
    print(f"Load factor                = {LOAD_FACTOR:.4f}")
    print(f"Torsion offset             = {LOAD_OFFSET_M:.6f} m")
    print(f"Maximum |V|                = {np.max(np.abs(V)):.6f} N")
    print(f"Maximum |M|                = {np.max(np.abs(M)):.6f} N-m")
    print(f"Maximum |T|                = {np.max(np.abs(T)):.6f} N-m")
    print()
    print(f"Created: {OUTPUT_CSV}")
    print(f"Created: {OUTPUT_TXT}")
    print(f"Created: {OUTPUT_PLOT}")
    print()
    print(
        "The CSV can be selected in the sheet-flow GUI as the span-load file."
    )


if __name__ == "__main__":
    main()
