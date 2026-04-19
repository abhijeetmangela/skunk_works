"""
Flight Dynamics & Stability Analysis for the IITM fixed-wing UAV.

This version is cleaned up to use the UAV report data rather than the
unrelated GA-aircraft example values.

It computes:
  • Longitudinal static stability and trim
  • Longitudinal 1st-order A-matrix and eigenvalues
  • Lateral-directional 1st-order A-matrix and eigenvalues
  • Eigenvalue plots for both subsystems

Outputs:
  • Console summary
  • PNG plots saved next to the script
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar


D2R = np.pi / 180.0
R2D = 180.0 / np.pi


@dataclass
class UAVInputs:
    # Geometry / mass from the report
    S: float = 0.5307
    b: float = 2.0
    c_bar: float = 0.2653
    AR: float = 7.53
    m: float = 6.592
    g: float = 9.81
    rho: float = 1.142
    x_cg: float = 0.264

    # Wing / fuselage stability data from Chapter 6 of the report
    x_ac_w: float = 0.1985
    x_ac_f: float = 0.09
    Cmac_w: float = -0.0683
    Cmac_f: float = 0.0
    CL_w: float = 0.378058
    CL_f: float = 0.061015
    CLalpha_w: float = 3.6107
    CLalpha_f: float = 0.7952
    Cm_alpha_f: float = 0.130214

    # Tail sizing from the report
    S_t: float = 0.1327
    AR_t: float = 4.5
    l_t: float = 0.7
    CLalpha_t: float = 3.8310
    eta_t: float = 0.99

    S_v: float = 0.0478
    AR_v: float = 1.5
    l_v: float = 0.7
    CLalpha_v: float = 3.83
    eta_v: float = 0.95  # not explicitly stated in the report; standard assumption

    # Aerodynamic / propulsion data from the report
    e: float = 0.7648
    CD0: float = 0.02873
    alpha_trim_deg: float = 3.0
    i_w_deg: float = 3.0
    i_t_deg_guess: float = 0.0

    # Mass properties from the SolidWorks mass-properties screenshot
    # in Fig. 6.1 of the report.
    Ixx: float = 0.730
    Iyy: float = 1.161
    Izz: float = 0.493

    # Geometry assumptions needed only for the lateral-directional model
    # (these are not explicitly tabulated in the report).
    D_f: float = 0.30   # central body width proxy
    z_w: float = 0.00   # wing vertical offset proxy
    Gamma_deg: float = 0.0
    Lambda_c4_deg: float = 0.0


def print_section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def print_subsection(title: str) -> None:
    print(f"\n--- {title} ---")


def solve_tail_incidence(inp: UAVInputs, de_da: float, Cm_f: float) -> float:
    """Solve the tail incidence angle from the trim moment condition."""
    alpha = inp.alpha_trim_deg * D2R

    def residual(it_rad: float) -> float:
        return (
            inp.Cmac_w
            + inp.CL_w * (inp.x_cg - inp.x_ac_w) / inp.c_bar
            - inp.eta_t
            * (inp.S_t / inp.S)
            * (inp.CLalpha_t * (alpha + it_rad - de_da * alpha))
            * (inp.l_t / inp.c_bar)
            + Cm_f
        )

    # The report gives a value close to -0.31 deg.
    sol = root_scalar(
        residual,
        bracket=(np.deg2rad(-10.0), np.deg2rad(10.0)),
        method="brentq",
    )
    return float(sol.root)


def eigenvalue_plot(evals: np.ndarray, title: str, filename: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)
    ax.scatter(evals.real, evals.imag, s=70)

    for i, ev in enumerate(evals, start=1):
        ax.annotate(f"λ{i}", (ev.real, ev.imag), textcoords="offset points", xytext=(6, 6))

    ax.set_xlabel("Real part (1/s)")
    ax.set_ylabel("Imaginary part (1/s)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='datalim')
    fig.tight_layout()
    fig.savefig(filename, dpi=200)
    plt.close(fig)


def compute_longitudinal(inp: UAVInputs) -> Dict[str, float | np.ndarray]:
    print_section("LONGITUDINAL STABILITY")

    # Static stability derivatives
    de_da = (2 * inp.CLalpha_w) / (np.pi * inp.e * inp.AR)
    Cm_f = inp.Cmac_f + inp.CL_f * ((inp.x_cg - inp.x_ac_f) / inp.c_bar)
    V_H = (inp.S_t * inp.l_t) / (inp.S * inp.c_bar)
    X_np = inp.x_ac_w + inp.c_bar * (
        inp.eta_t * V_H * (inp.CLalpha_t / inp.CLalpha_w) * (1 - de_da)
        - inp.Cm_alpha_f / inp.CLalpha_w
    )
    Kn = (X_np - inp.x_cg) / inp.c_bar

    print(f"dε/dα = {de_da:.4f} /rad")
    print(f"Fuselage moment term Cm_f = {Cm_f:.5f}")
    print(f"Horizontal tail volume coefficient V_H = {V_H:.4f}")
    print(f"Neutral point X_np = {X_np:.5f} m")
    print(f"Static margin K_n = {Kn:.5f} ({Kn * 100:.2f} %)")

    it_rad = solve_tail_incidence(inp, de_da, Cm_f)
    alpha = inp.alpha_trim_deg * D2R
    alpha_t = alpha + it_rad - de_da * alpha
    CL_t = inp.CLalpha_t * alpha_t

    Cm_alpha = (
        ((inp.x_cg - inp.x_ac_w) / inp.c_bar) * inp.CLalpha_w
        - inp.eta_t * V_H * inp.CLalpha_t * (1 - de_da)
        + inp.Cm_alpha_f
    )
    Cm0 = -Cm_alpha * alpha
    Cm_trim = Cm0 + Cm_alpha * alpha

    # A first trim CL estimate using the report's wing + tail contributions
    CL_trim = inp.CL_w + inp.eta_t * (inp.S_t / inp.S) * CL_t
    V_trim = math.sqrt(2 * inp.m * inp.g / (inp.rho * inp.S * CL_trim))
    q_trim = 0.5 * inp.rho * V_trim**2
    Ma_trim = V_trim / 340.0

    K = 1.0 / (np.pi * inp.AR * inp.e)
    CD_trim = inp.CD0 + K * CL_trim**2
    T_req = q_trim * inp.S * CD_trim
    P_req = T_req * V_trim

    # Stability derivatives used in the first-order form
    qSc_Iyy = q_trim * inp.S * inp.c_bar / inp.Iyy
    qS_W = q_trim * inp.S / (inp.m * inp.g)
    g_V = inp.g / V_trim
    c_2V = inp.c_bar / (2 * V_trim)

    Cm_q1 = -2 * inp.eta_t * V_H * inp.CLalpha_t * (inp.l_t / inp.c_bar)
    CL_q1 = 2 * inp.eta_t * V_H * inp.CLalpha_t
    Cm_adot = 0.0
    CL_adot = 0.0
    CD_alpha = 2 * K * CL_trim * inp.CLalpha_w
    CD_q1 = 2 * K * CL_trim * CL_q1
    CL_Ma = CL_trim * Ma_trim / (1 - Ma_trim**2)
    CD_Ma = 0.0

    # Second-order approximate modes
    wn_SP = math.sqrt(max(1e-12, -qSc_Iyy * Cm_alpha))
    two_zeta_wn_SP = -qSc_Iyy * c_2V * (Cm_q1 + Cm_adot)
    zeta_SP = two_zeta_wn_SP / (2 * wn_SP)

    wn_Ph = math.sqrt(max(1e-12, (g_V**2) * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)))
    two_zeta_wn_Ph = g_V * qS_W * (Ma_trim * 0.0 + 2 * CD_trim)
    zeta_Ph = two_zeta_wn_Ph / (2 * wn_Ph)

    # Longitudinal A-matrix in the same first-order structure as the supplied code
    a11 = -g_V * qS_W * (Ma_trim * CD_Ma + 2 * CD_trim)
    a12 = -g_V
    a13 = -g_V * qS_W * CD_alpha
    a14 = -g_V * qS_W * CD_q1 * c_2V

    a21 = g_V * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    a22 = 0.0
    a23 = g_V * qS_W * (inp.CLalpha_w + inp.eta_t * (inp.S_t / inp.S) * inp.CLalpha_t * (1 - de_da))
    a24 = g_V * qS_W * CL_q1 * c_2V

    k_alpha_dot = 1 + g_V * qS_W * CL_q1 * c_2V
    k_gamma = (g_V**2) * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    k_alpha = ((g_V**2) * (qS_W**2) * (Ma_trim * CL_Ma + 2 * CL_trim) * CD_alpha + qSc_Iyy * Cm_alpha)
    k_alphadot_r = (
        qSc_Iyy * c_2V * Cm_q1
        - g_V * qS_W * (inp.CLalpha_w + inp.eta_t * (inp.S_t / inp.S) * inp.CLalpha_t * (1 - de_da))
        + (g_V**2) * (qS_W**2) * CD_q1 * c_2V * (Ma_trim * CL_Ma + 2 * CL_trim)
    )
    k_V_r = (g_V**2) * (qS_W**2) * (Ma_trim * CL_Ma + 2 * CL_trim) * (Ma_trim * CD_Ma + 2 * CD_trim)

    a41 = k_V_r / k_alpha_dot
    a42 = k_gamma / k_alpha_dot
    a43 = k_alpha / k_alpha_dot
    a44 = k_alphadot_r / k_alpha_dot

    A_long = np.array([
        [a11, a12, a13, a14],
        [a21, a22, a23, a24],
        [0.0, 0.0, 0.0, 1.0],
        [a41, a42, a43, a44],
    ])

    evals_long = np.linalg.eigvals(A_long)

    print(f"Tail incidence it = {it_rad * R2D:.4f} deg")
    print(f"Tail angle of attack α_t = {alpha_t:.5f} rad")
    print(f"Tail lift coefficient C_L,t = {CL_t:.5f}")
    print(f"Cm_alpha = {Cm_alpha:.5f} /rad")
    print(f"Cm0 = {Cm0:.5f}")
    print(f"Trim CL ≈ {CL_trim:.5f}")
    print(f"Trim speed V* ≈ {V_trim:.3f} m/s")
    print(f"Trim drag coefficient C_D ≈ {CD_trim:.5f}")
    print(f"Thrust required ≈ {T_req:.3f} N")
    print(f"Power required ≈ {P_req:.3f} W")

    print("\nLongitudinal A-matrix:")
    print(np.array2string(A_long, precision=5, suppress_small=True))

    print("\nLongitudinal eigenvalues:")
    for ev in evals_long:
        print(f"  {ev.real:+.6f} {ev.imag:+.6f}j")

    # Sort complex pair by frequency for reporting
    complex_evals = [ev for ev in evals_long if abs(ev.imag) > 1e-6]
    complex_evals.sort(key=lambda ev: abs(ev.imag), reverse=True)
    if len(complex_evals) >= 2:
        sp = complex_evals[0]
        ph = complex_evals[-1]
        print(f"\nShort-period pair: λ = {sp.real:+.5f} ± j{abs(sp.imag):.5f}")
        print(f"Phugoid pair:      λ = {ph.real:+.5f} ± j{abs(ph.imag):.5f}")

    return {
        "de_da": de_da,
        "Cm_f": Cm_f,
        "V_H": V_H,
        "X_np": X_np,
        "Kn": Kn,
        "it_rad": it_rad,
        "alpha_t": alpha_t,
        "CL_t": CL_t,
        "Cm_alpha": Cm_alpha,
        "Cm0": Cm0,
        "CL_trim": CL_trim,
        "V_trim": V_trim,
        "CD_trim": CD_trim,
        "T_req": T_req,
        "P_req": P_req,
        "wn_SP": wn_SP,
        "zeta_SP": zeta_SP,
        "wn_Ph": wn_Ph,
        "zeta_Ph": zeta_Ph,
        "A_long": A_long,
        "evals_long": evals_long,
        "q_trim": q_trim,
        "Ma_trim": Ma_trim,
        "qSc_Iyy": qSc_Iyy,
        "qS_W": qS_W,
        "g_V": g_V,
        "c_2V": c_2V,
        "Cm_q1": Cm_q1,
        "CL_q1": CL_q1,
        "CD_alpha": CD_alpha,
        "CD_q1": CD_q1,
        "CL_Ma": CL_Ma,
        "CD_Ma": CD_Ma,
        "CLalpha": inp.CLalpha_w,
        "K": K,
    }


def compute_lateral(inp: UAVInputs, long: Dict[str, float | np.ndarray]) -> Dict[str, float | np.ndarray]:
    print_section("LATERAL-DIRECTIONAL STABILITY")

    V_trim = float(long["V_trim"])
    q_trim = float(long["q_trim"])
    g_V = float(long["g_V"])
    CL_trim = float(long["CL_trim"])
    CD_trim = float(long["CD_trim"])
    CD_alpha = float(long["CD_alpha"])
    CLalpha = float(long["CLalpha"])
    alpha_t = float(long["alpha_t"])

    Gamma_rad = inp.Gamma_deg * D2R
    Lambda_c4_rad = inp.Lambda_c4_deg * D2R

    # Simplified sidewash relation from the book-style derivation.
    # D_f and z_w are report-unsupported proxies and are kept explicit here.
    sidewash = (
        0.724
        + 3.06 * (inp.S_v / inp.S) / (1 + np.cos(Lambda_c4_rad))
        + 0.4 * inp.z_w / inp.D_f
        + 0.009 * inp.AR
    )

    CY_beta_w = -0.0001 * inp.Gamma_deg * D2R
    CY_beta_VT = -(inp.eta_v * inp.CLalpha_v * (inp.S_v / inp.S) * sidewash)
    CY_beta = CY_beta_w + CY_beta_VT

    Cn_beta_w = -0.075 * inp.Gamma_deg * D2R * CL_trim
    Cn_beta_VT = -CY_beta_VT * (inp.l_v / inp.b)
    Cn_beta = Cn_beta_w + Cn_beta_VT

    z_VT = (4.0 / 9.0) * np.sqrt(inp.AR_v * inp.S_v)
    z_eff = z_VT * np.cos(alpha_t) - inp.l_v * np.sin(alpha_t)

    Cl_beta_w = -(Gamma_rad * inp.CLalpha_w * inp.c_bar * inp.b) / (6 * inp.S)
    Cl_beta_VT = (z_eff / inp.b) * CY_beta_VT
    Cl_beta = Cl_beta_w + Cl_beta_VT

    CY_p2_VT = (2 / inp.b) * (z_eff - z_VT) * CY_beta_VT
    Cl_p2_VT = 2 * (z_eff / inp.b) * ((z_eff - z_VT) / inp.b) * CY_beta_VT
    Cl_p2_w = -(1 / 6) * (CLalpha + CD_trim)
    Cl_p2 = Cl_p2_VT + Cl_p2_w

    lVT_eff = inp.l_v * np.cos(alpha_t) + z_VT * np.sin(alpha_t)
    Cn_p2_w = -(1 / 6) * (CL_trim - CD_alpha)
    Cn_p2_VT = (-(2 / inp.b) * lVT_eff * ((z_eff - z_VT) / inp.b) * CY_beta_VT)
    Cn_p2 = Cn_p2_w + Cn_p2_VT

    CY_r1_VT = -(2 / inp.b) * lVT_eff * CY_beta_VT
    Cl_r1_w = CL_trim / 3
    Cl_r1_VT = (-(2 / inp.b**2) * lVT_eff * z_eff * CY_beta_VT)
    Cl_r1 = Cl_r1_w + Cl_r1_VT
    Cn_r1_VT = (2 / inp.b**2) * lVT_eff**2 * CY_beta_VT
    Cn_r1 = Cn_r1_VT

    Cn_r2 = -CD_trim / 3
    Cl_r2 = CL_trim / 3

    qSb_Ixx = q_trim * inp.S * inp.b / inp.Ixx
    qSb_Izz = q_trim * inp.S * inp.b / inp.Izz
    qS_m = q_trim * inp.S / inp.m
    b_2V = inp.b / (2 * V_trim)

    N_beta = qSb_Izz * Cn_beta
    N_r2 = qSb_Izz * Cn_r2 * b_2V
    Y_beta = qS_m * CY_beta
    L_beta = qSb_Ixx * Cl_beta
    L_p2 = qSb_Ixx * Cl_p2 * b_2V
    N_r1 = qSb_Izz * Cn_r1 * b_2V
    L_r1 = qSb_Ixx * Cl_r1 * b_2V
    L_r2 = qSb_Ixx * Cl_r2 * b_2V

    lambda_roll = L_p2
    wn_DR_sq = N_beta + g_V * (Y_beta * N_r2 + L_beta / lambda_roll)
    wn_DR = math.sqrt(abs(wn_DR_sq))
    two_zeta_wn_DR = -N_r1 - g_V * (Y_beta + L_r1 / lambda_roll)
    zeta_DR = two_zeta_wn_DR / (2 * wn_DR)
    lambda_spiral = g_V * (L_beta * N_r2 - N_beta * L_r2) / (lambda_roll * wn_DR_sq)

    A_lat = np.array([
        [0.0, 1.0, 0.0, 0.0],
        [g_V * L_r2, L_p2, L_beta + g_V * Y_beta * L_r2, -L_r1],
        [0.0, 0.0, 0.0, 1.0],
        [-g_V * N_r2, g_V, -(N_beta + g_V * Y_beta * N_r2), N_r1 + g_V * Y_beta],
    ])

    evals_lat = np.linalg.eigvals(A_lat)

    print(f"Sidewash factor = {sidewash:.4f}")
    print(f"CY_beta = {CY_beta:.5f} /rad")
    print(f"Cn_beta = {Cn_beta:.5f} /rad")
    print(f"Cl_beta = {Cl_beta:.5f} /rad")
    print(f"Roll mode λ = {lambda_roll:.5f} 1/s")
    print(f"Dutch roll ω_n ≈ {wn_DR:.5f} rad/s")
    print(f"Dutch roll ζ ≈ {zeta_DR:.5f}")
    print(f"Spiral mode λ ≈ {lambda_spiral:.5f} 1/s")

    print("\nLateral-directional A-matrix:")
    print(np.array2string(A_lat, precision=5, suppress_small=True))

    print("\nLateral-directional eigenvalues:")
    for ev in evals_lat:
        print(f"  {ev.real:+.6f} {ev.imag:+.6f}j")

    return {
        "sidewash": sidewash,
        "CY_beta": CY_beta,
        "Cn_beta": Cn_beta,
        "Cl_beta": Cl_beta,
        "lambda_roll": lambda_roll,
        "wn_DR": wn_DR,
        "zeta_DR": zeta_DR,
        "lambda_spiral": lambda_spiral,
        "A_lat": A_lat,
        "evals_lat": evals_lat,
        "CY_p2_VT": CY_p2_VT,
        "Cl_p2": Cl_p2,
        "Cn_p2": Cn_p2,
        "CY_r1_VT": CY_r1_VT,
        "Cl_r1": Cl_r1,
        "Cn_r1": Cn_r1,
        "Cn_r2": Cn_r2,
        "Cl_r2": Cl_r2,
    }


def main() -> None:
    inp = UAVInputs()
    print_section("UAV STABILITY ANALYSIS")
    print("Using the fixed-wing UAV design data from the report:")
    print(f"  Wing area S = {inp.S} m²")
    print(f"  Wingspan b = {inp.b} m")
    print(f"  MTOW = {inp.m} kg")
    print(f"  Horizontal tail area = {inp.S_t} m²")
    print(f"  Static-margin target ≈ 14%\n")

    long = compute_longitudinal(inp)
    lat = compute_lateral(inp, long)

    # Save eigenvalue plots
    out_dir = Path(__file__).resolve().parent
    eigenvalue_plot(long["evals_long"], "Longitudinal eigenvalues", out_dir / "longitudinal_eigenvalues.png")
    eigenvalue_plot(lat["evals_lat"], "Lateral-directional eigenvalues", out_dir / "lateral_eigenvalues.png")

    print_section("SUMMARY")
    print(f"Static margin: {long['Kn'] * 100:.2f} %")
    print(f"Tail incidence: {long['it_rad'] * R2D:.3f} deg")
    print(f"Trim speed: {long['V_trim']:.3f} m/s")
    print(f"Trim CL: {long['CL_trim']:.4f}")
    print(f"Longitudinal short-period ω_n ≈ {long['wn_SP']:.4f} rad/s")
    print(f"Longitudinal phugoid ω_n ≈ {long['wn_Ph']:.4f} rad/s")
    print(f"Dutch-roll ω_n ≈ {lat['wn_DR']:.4f} rad/s")
    print(f"Dutch-roll ζ ≈ {lat['zeta_DR']:.4f}")
    print(f"Spiral λ ≈ {lat['lambda_spiral']:.5f} 1/s")
    print("\nSaved plots:")
    print(f"  {out_dir / 'longitudinal_eigenvalues.png'}")
    print(f"  {out_dir / 'lateral_eigenvalues.png'}")


if __name__ == "__main__":
    main()
