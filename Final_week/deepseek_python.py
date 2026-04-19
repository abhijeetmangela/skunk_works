"""
============================================================
  Flight Dynamics & Stability Analysis – UAV Design
  Based on Chapter 9 Appendix (Case Studies) and
  Group_11_Design_Report (IIT Madras)
============================================================
All formulas are taken directly from the supplied PDF.
Input parameters updated from the design report.
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
#  UTILITY
# ------------------------------------------------------------
D2R = np.pi / 180.0

def section_header(title):
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)

def sub_header(title):
    print(f"\n--- {title} ---")

# ------------------------------------------------------------
#  CASE 1 – UAV (from design report)
# ------------------------------------------------------------
def uav_aircraft():
    section_header("UAV – Stability Analysis (Group 11 Design)")

    # --------------------------------------------------------
    # (A) INPUT PARAMETERS (from report, pages 54,68,76,92,106)
    # --------------------------------------------------------
    sub_header("A) Input Parameters (Design Report)")

    # ── Geometry ──────────────────────────────────────────
    g        = 9.81          # m/s²
    S_w      = 0.5307        # m²   – wing planform area
    b_w      = 2.0           # m    – wing span
    c_r      = 0.2653        # m    – root chord (rectangular)
    c_t      = 0.2653        # m    – tip chord (rectangular)
    lam      = c_t / c_r     # taper ratio = 1.0
    Lambda_c4_deg = 0.0      # deg  – quarter‑chord sweep
    Gamma_deg     = 0.0      # deg  – dihedral
    i_w      = 3.0           # deg  – wing incidence

    S_HT     = 0.1327        # m²   – horizontal tail area
    AR_HT    = 4.5           #      – HT aspect ratio
    i_t      = -0.3126       # deg  – HT incidence (from page 89)
    X_AC_HT  = 0.964         # m    – HT aerodynamic centre from nose
                             #       (CG=0.264 + l_t=0.7)

    S_VT     = 0.0478        # m²   – vertical tail area
    AR_VT    = 1.5           #      – VT aspect ratio
    X_AC_VT  = 0.964         # m    – VT AC from nose

    l_f      = 1.2           # m    – fuselage length
    S_f      = 1.11          # m²   – fuselage side area (estimated)
    D_f      = 0.338         # m    – max fuselage diameter (from NACA0025)

    S_da     = 0.05          # m²   – aileron area (one side)
    l_da     = 0.76          # m    – aileron span
    S_de     = 0.04          # m²   – elevator area
    S_dr     = 0.01          # m²   – rudder area

    # ── Mass & Inertia ────────────────────────────────────
    m        = 6.298         # kg   – aircraft mass (page 85)
    I_xx     = 0.751         # kg·m² (CAD estimate)
    I_zz     = 0.384         # kg·m²
    I_yy     = 1.045         # kg·m²
    X_CG     = 0.264         # m    – CG from nose (page 102)

    # ── Aerodynamic parameters ───────────────────────────
    alpha_L0_deg = -3.5      # deg  – wing zero‑lift AoA (Selig S7055 approx)
    Cm_ac_w      = -0.005  #      – wing pitching moment at AC (page 76)
    SM           = 0.13      #      – static margin (13%, page 103)
    e_oswald     = 0.7648    #      – Oswald efficiency (page 30)
    C_D0         = 0.02873   #      – zero‑lift drag coefficient (page 106)

    # Lift‑curve slopes (2D airfoil, /rad)
    cl_alpha_w   = 3.6107    # /rad – Selig S7055 (page 76)
    cl_alpha_HT  = 6.24      # /rad – NACA0012 (0.109/deg)
    cl_alpha_VT  = 6.24      # /rad – NACA0012

    # ── Flight conditions ─────────────────────────────────
    rho          = 1.22      # kg/m³ (200 m altitude)
    a_sound      = 340.0     # m/s
    eta          = 0.99      # tail dynamic pressure ratio

    # ── Engine/propeller (not critical for stability) ────
    P_max        = 450       # N·m/s (max power)
    N_blades     = 4
    D_prop       = 0.3048    # m
    l_prop       = 0.12957   # m

    z_w          = 0.06      # m – vertical offset wing root

    print("  Inputs loaded (UAV design).")

    # --------------------------------------------------------
    # (B) DERIVED GEOMETRY
    # --------------------------------------------------------
    sub_header("B) Derived Geometry")

    # Mean aerodynamic chord (rectangular)
    c_bar = c_r   # since taper = 1
    # Wing aspect ratio
    AR_w = b_w**2 / S_w

    # Tail spans
    b_HT = np.sqrt(AR_HT * S_HT)
    b_VT = np.sqrt(AR_VT * S_VT)

    # VT AC height above fuselage reference line
    z_VT = (4/9) * b_VT

    # Moment arms
    l_HT = X_AC_HT - X_CG
    l_VT = X_AC_VT - X_CG

    # Tail volume ratios
    V_HT = (S_HT * l_HT) / (S_w * c_bar)
    V_VT = (S_VT * l_VT) / (S_w * b_w)

    # Non‑dimensional CG location
    h_CG = X_CG / c_bar

    # Induced drag factor
    K = 1.0 / (np.pi * AR_w * e_oswald)

    print(f"  c_bar   = {c_bar:.4f} m")
    print(f"  AR_w    = {AR_w:.4f}")
    print(f"  lam     = {lam:.2f}")
    print(f"  b_HT    = {b_HT:.3f} m,   b_VT = {b_VT:.3f} m")
    print(f"  z_VT    = {z_VT:.3f} m")
    print(f"  l_HT    = {l_HT:.3f} m,   l_VT = {l_VT:.3f} m")
    print(f"  V_HT    = {V_HT:.4f},      V_VT = {V_VT:.4f}")
    print(f"  h_CG    = {h_CG:.3f}")
    print(f"  K       = {K:.4f}")

    # --------------------------------------------------------
    # (C) AERODYNAMIC COEFFICIENTS
    # --------------------------------------------------------
    sub_header("C) Aerodynamic Coefficients")

    # Finite‑wing lift‑curve slopes (subsonic)
    CL_alpha_w  = cl_alpha_w  / (1 + cl_alpha_w  / (np.pi * AR_w))
    CL_alpha_HT = cl_alpha_HT / (1 + cl_alpha_HT / (np.pi * AR_HT))
    CL_alpha_VT = cl_alpha_VT / (1 + cl_alpha_VT / (np.pi * AR_VT))

    print(f"  CL_alpha_w   = {CL_alpha_w:.3f} /rad")
    print(f"  CL_alpha_HT  = {CL_alpha_HT:.3f} /rad")
    print(f"  CL_alpha_VT  = {CL_alpha_VT:.3f} /rad")

    # Wing lift coefficient at zero AoA
    CL0_w = CL_alpha_w * abs(alpha_L0_deg) * D2R
    print(f"  CL0_w        = {CL0_w:.4f}")

    # Downwash gradient
    eps_alpha = 2 * CL_alpha_w / (np.pi * AR_w)
    eps_0     = 2 * CL0_w / (np.pi * AR_w)
    print(f"  de/da        = {eps_alpha:.4f}")
    print(f"  eps_0        = {eps_0:.4f}")

    # --------------------------------------------------------
    # (D) TRIM CALCULATIONS
    # --------------------------------------------------------
    sub_header("D) Trim Calculations")

    # Neutral point (non‑dimensional)
    h_NP = h_CG + SM
    X_NP = h_NP * c_bar

    # Wing‑body AC location (from NP and tail contribution)
    h_AC_wb = h_NP - V_HT * (CL_alpha_HT / CL_alpha_w) * (1 - eps_alpha)
    X_AC_wb = h_AC_wb * c_bar

    print(f"  h_NP         = {h_NP:.3f},   X_NP   = {X_NP:.3f} m")
    print(f"  h_AC_wb      = {h_AC_wb:.3f},  X_AC_wb = {X_AC_wb:.3f} m")

    # Pitching moment derivatives
    Cm_alpha_w  = ((X_CG - X_AC_wb) / c_bar) * CL_alpha_w
    Cm_alpha_HT = -eta * V_HT * CL_alpha_HT * (1 - eps_alpha)
    Cm_alpha    = Cm_alpha_w + Cm_alpha_HT

    # Zero‑angle pitching moment (elevator fixed)
    Cm0 = Cm_ac_w + eta * V_HT * CL_alpha_HT * eps_0

    # Trim AoA (zero elevator)
    alpha_trim_rad = -Cm0 / Cm_alpha
    alpha_trim_deg = alpha_trim_rad / D2R

    # Whole‑aircraft CL_alpha and CL0
    CL0 = (CL_alpha_w * abs(alpha_L0_deg) * D2R
           + eta * (S_HT / S_w) * CL_alpha_HT * (i_t * D2R - eps_0))
    CL_alpha = (CL_alpha_w
                + eta * (S_HT / S_w) * CL_alpha_HT * (1 - eps_alpha))

    # Trim lift coefficient
    CL_trim = CL0 + CL_alpha * alpha_trim_rad

    # Trim speed from L = W
    V_trim = np.sqrt(2 / rho * (m * g / S_w) * (1 / CL_trim))

    # Dynamic pressure
    q_trim = 0.5 * rho * V_trim**2

    # Mach number
    Ma_trim = V_trim / a_sound

    # Trim drag
    CD_trim = C_D0 + K * CL_trim**2

    # Thrust and power required
    T_req = q_trim * S_w * CD_trim
    P_req = T_req * V_trim

    print(f"  Cm_alpha_w   = {Cm_alpha_w:.4f} /rad")
    print(f"  Cm_alpha_HT  = {Cm_alpha_HT:.4f} /rad")
    print(f"  Cm_alpha     = {Cm_alpha:.4f} /rad")
    print(f"  Cm0          = {Cm0:.4f}")
    print(f"  CL0          = {CL0:.4f}")
    print(f"  CL_alpha     = {CL_alpha:.4f} /rad")
    print(f"  Alpha_trim   = {alpha_trim_rad:.4f} rad  ({alpha_trim_deg:.2f} deg)")
    print(f"  CL_trim      = {CL_trim:.4f}")
    print(f"  V_trim       = {V_trim:.3f} m/s")
    print(f"  Ma_trim      = {Ma_trim:.4f}")
    print(f"  q_trim       = {q_trim:.3f} N/m²")
    print(f"  CD_trim      = {CD_trim:.4f}")
    print(f"  Thrust req.  = {T_req:.2f} N")
    print(f"  Power req.   = {P_req:.2f} W  ({P_req/745.7:.1f} hp)")

    # --------------------------------------------------------
    # (E) LONGITUDINAL STABILITY DERIVATIVES
    # --------------------------------------------------------
    sub_header("E) Longitudinal Stability Derivatives")

    Cm_q1 = -2 * eta * V_HT * CL_alpha_HT * (l_HT / c_bar)
    CL_q1 = 2 * eta * V_HT * CL_alpha_HT
    Cm_adot = 0.0
    CL_adot = 0.0

    CD_alpha = 2 * K * CL_trim * CL_alpha
    CD_q1    = 2 * K * CL_trim * CL_q1

    CL_Ma = CL_trim * Ma_trim / (1 - Ma_trim**2)
    CD_Ma = 0.0

    qSc_Iyy  = q_trim * S_w * c_bar / I_yy
    qS_W     = q_trim * S_w / (m * g)
    g_V      = g / V_trim
    c_2V     = c_bar / (2 * V_trim)

    print(f"  Cm_q1        = {Cm_q1:.4f} /rad")
    print(f"  CL_q1        = {CL_q1:.4f} /rad")
    print(f"  CD_alpha     = {CD_alpha:.4f} /rad")
    print(f"  CD_q1        = {CD_q1:.4f} /rad")
    print(f"  CL_Ma        = {CL_Ma:.4f}")
    print(f"  q̄*Sc/Iyy     = {qSc_Iyy:.4f} /s²")
    print(f"  q̄*S/W        = {qS_W:.4f}")
    print(f"  g/V*         = {g_V:.4f}")
    print(f"  c/(2V*)      = {c_2V:.5f} s")

    # --------------------------------------------------------
    # (F) LONGITUDINAL DYNAMIC MODES (2nd-order approx.)
    # --------------------------------------------------------
    sub_header("F) Longitudinal Dynamic Modes (Approximate)")

    # Short period
    wn_SP_sq = -qSc_Iyy * Cm_alpha
    wn_SP    = np.sqrt(wn_SP_sq)
    T_SP     = 2 * np.pi / wn_SP
    two_zeta_wn_SP = -qSc_Iyy * c_2V * (Cm_q1 + Cm_adot)
    zeta_SP  = two_zeta_wn_SP / (2 * wn_SP)

    # Phugoid
    wn_Ph_sq = (g_V**2) * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    wn_Ph    = np.sqrt(wn_Ph_sq)
    T_Ph     = 2 * np.pi / wn_Ph
    two_zeta_wn_Ph = g_V * qS_W * (Ma_trim * CD_Ma + 2 * CD_trim)
    zeta_Ph  = two_zeta_wn_Ph / (2 * wn_Ph)

    print(f"\n  SHORT PERIOD MODE")
    print(f"    omega_n   = {wn_SP:.4f} rad/s")
    print(f"    zeta      = {zeta_SP:.4f}")
    print(f"    Period    = {T_SP:.4f} s")

    print(f"\n  PHUGOID MODE")
    print(f"    omega_n   = {wn_Ph:.4f} rad/s")
    print(f"    zeta      = {zeta_Ph:.4f}")
    print(f"    Period    = {T_Ph:.4f} s")

    # --------------------------------------------------------
    # (G) LONGITUDINAL STATE-SPACE MATRIX (1st order)
    # --------------------------------------------------------
    sub_header("G) Longitudinal State-Space Matrix A_long")

    a11 = -g_V * qS_W * (Ma_trim * CD_Ma + 2 * CD_trim)
    a12 = -g_V
    a13 = -g_V * qS_W * CD_alpha
    a14 = -g_V * qS_W * CD_q1 * c_2V

    a21 =  g_V * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    a22 =  0.0
    a23 =  g_V * qS_W * CL_alpha
    a24 =  g_V * qS_W * CL_q1 * c_2V

    a31 = 0.0
    a32 = 0.0
    a33 = 0.0
    a34 = 1.0

    # k-factors for α̈ row
    k_alpha_dot  = 1 + g_V * qS_W * CL_q1 * c_2V
    k_gamma      = (g_V**2) * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    k_alpha      = ((g_V**2) * (qS_W**2)
                    * (Ma_trim * CL_Ma + 2 * CL_trim) * CD_alpha
                    + qSc_Iyy * Cm_alpha)
    k_alphadot_r = (qSc_Iyy * c_2V * Cm_q1
                    - g_V * qS_W * CL_alpha
                    + (g_V**2) * (qS_W**2) * CD_q1 * c_2V
                    * (Ma_trim * CL_Ma + 2 * CL_trim))
    k_V_r        = ((g_V**2) * (qS_W**2)
                    * (Ma_trim * CL_Ma + 2 * CL_trim)
                    * (Ma_trim * CD_Ma + 2 * CD_trim)
                    + qSc_Iyy * 0.0)   # Cm_Ma = 0

    a41 = k_V_r      / k_alpha_dot
    a42 = k_gamma    / k_alpha_dot
    a43 = k_alpha    / k_alpha_dot
    a44 = k_alphadot_r / k_alpha_dot

    A_long = np.array([
        [a11, a12, a13, a14],
        [a21, a22, a23, a24],
        [a31, a32, a33, a34],
        [a41, a42, a43, a44],
    ])

    print("\n  A_long =")
    print(np.array2string(A_long, precision=4, suppress_small=True,
                          formatter={'float_kind': lambda x: f"{x:10.4f}"}))

    evals_long = np.linalg.eigvals(A_long)
    print("\n  Eigenvalues of A_long:")
    for ev in evals_long:
        print(f"    λ = {ev.real:+.5f}  {'+' if ev.imag >= 0 else '-'}  j{abs(ev.imag):.5f}")

    # --------------------------------------------------------
    # (H) LATERAL-DIRECTIONAL AERODYNAMICS
    # --------------------------------------------------------
    sub_header("H) Lateral-Directional Aerodynamic Coefficients")

    Gamma_rad    = Gamma_deg * D2R
    alpha_t_rad  = alpha_trim_rad

    # Sidewash factor
    Lambda_c4_rad = Lambda_c4_deg * D2R
    sidewash = (0.724
                + 3.06 * (S_VT / S_w) / (1 + np.cos(Lambda_c4_rad))
                + 0.4 * z_w / D_f
                + 0.009 * AR_w)

    # Side-force coefficients
    CY_beta_w  = -0.0001 * Gamma_deg * D2R
    CY_beta_VT = (-eta * CL_alpha_VT * (S_VT / S_w) * sidewash)
    CY_beta    = CY_beta_w + CY_beta_VT

    # Yawing-moment coefficients
    Cn_beta_w  = -0.075 * Gamma_deg * D2R * CL_trim
    Cn_beta_VT = -CY_beta_VT * (l_VT / b_w)
    Cn_beta    = Cn_beta_w + Cn_beta_VT

    # Rolling-moment coefficients
    z_eff = z_VT * np.cos(alpha_t_rad) - l_VT * np.sin(alpha_t_rad)
    Cl_beta_w  = -(Gamma_rad * CL_alpha_w * c_bar * b_w) / (6 * S_w)
    Cl_beta_VT = (z_eff / b_w) * CY_beta_VT
    Cl_beta    = Cl_beta_w + Cl_beta_VT

    print(f"  Sidewash factor = {sidewash:.4f}")
    print(f"  CY_beta  = {CY_beta:.4f} /rad")
    print(f"  Cn_beta  = {Cn_beta:.4f} /rad")
    print(f"  Cl_beta  = {Cl_beta:.4f} /rad")

    # Roll-rate derivatives
    CY_p2_VT = (2 / b_w) * (z_eff - z_VT) * CY_beta_VT
    CY_p2    = CY_p2_VT

    Cl_p2_VT = 2 * (z_eff / b_w) * ((z_eff - z_VT) / b_w) * CY_beta_VT
    Cl_p2_w  = -(1/6) * (CL_alpha + CD_trim)
    Cl_p2    = Cl_p2_VT + Cl_p2_w

    lVT_eff = l_VT * np.cos(alpha_t_rad) + z_VT * np.sin(alpha_t_rad)
    Cn_p2_w  = -(1/6) * (CL_trim - CD_alpha)
    Cn_p2_VT = (-(2 / b_w) * lVT_eff * ((z_eff - z_VT) / b_w) * CY_beta_VT)
    Cn_p2    = Cn_p2_w + Cn_p2_VT

    print(f"\n  CY_p2    = {CY_p2:.4f} /rad")
    print(f"  Cl_p2    = {Cl_p2:.4f} /rad")
    print(f"  Cn_p2    = {Cn_p2:.4f} /rad")

    # Yaw-rate derivatives
    CY_r1_VT = -(2 / b_w) * lVT_eff * CY_beta_VT
    Cl_r1_w  = CL_trim / 3
    Cl_r1_VT = (-(2 / b_w**2) * lVT_eff * z_eff * CY_beta_VT)
    Cl_r1    = Cl_r1_w + Cl_r1_VT
    Cn_r1_VT = (2 / b_w**2) * lVT_eff**2 * CY_beta_VT
    Cn_r1    = Cn_r1_VT

    Cn_r2    = -CD_trim / 3
    Cl_r2    =  CL_trim / 3

    print(f"  CY_r1    = {CY_r1_VT:.4f} /rad")
    print(f"  Cl_r1    = {Cl_r1:.4f} /rad")
    print(f"  Cn_r1    = {Cn_r1:.4f} /rad")
    print(f"  Cn_r2    = {Cn_r2:.4f} /rad")
    print(f"  Cl_r2    = {Cl_r2:.4f} /rad")

    # --------------------------------------------------------
    # (I) LATERAL-DIRECTIONAL DYNAMIC MODES (2nd-order)
    # --------------------------------------------------------
    sub_header("I) Lateral-Directional Dynamic Modes (Approximate)")

    qSb_Ixx = q_trim * S_w * b_w / I_xx
    qSb_Izz = q_trim * S_w * b_w / I_zz
    qS_m    = q_trim * S_w / m
    b_2V    = b_w / (2 * V_trim)

    N_beta  = qSb_Izz * Cn_beta
    N_r2    = qSb_Izz * Cn_r2 * b_2V
    Y_beta  = qS_m    * CY_beta
    L_beta  = qSb_Ixx * Cl_beta
    L_p2    = qSb_Ixx * Cl_p2  * b_2V
    N_r1    = qSb_Izz * Cn_r1  * b_2V
    L_r1    = qSb_Ixx * Cl_r1  * b_2V
    L_r2    = qSb_Ixx * Cl_r2  * b_2V

    print(f"  N_beta  = {N_beta:.4f} /s²")
    print(f"  Y_beta  = {Y_beta:.4f} m/s²")
    print(f"  L_beta  = {L_beta:.4f} /s²")
    print(f"  L_p2    = {L_p2:.4f} /s²")
    print(f"  N_r1    = {N_r1:.4f} /s²")
    print(f"  L_r1    = {L_r1:.4f} /s²")

    # Roll mode
    lambda_roll = L_p2

    # Dutch roll
    wn_DR_sq = (N_beta
                + g_V * (Y_beta * N_r2 + L_beta / L_p2))
    wn_DR    = np.sqrt(abs(wn_DR_sq))
    two_zeta_wn_DR = (-N_r1
                      - g_V * (Y_beta + L_r1 / L_p2))
    zeta_DR  = two_zeta_wn_DR / (2 * wn_DR)

    # Spiral mode
    lambda_spiral = (g_V * (L_beta * N_r2 - N_beta * L_r2)
                     / (lambda_roll * wn_DR_sq))

    print(f"\n  ROLL MODE")
    print(f"    lambda_r  = {lambda_roll:.4f} /s")

    print(f"\n  DUTCH ROLL MODE")
    print(f"    omega_n   = {wn_DR:.4f} rad/s")
    print(f"    zeta      = {zeta_DR:.4f}")

    print(f"\n  SPIRAL MODE")
    print(f"    lambda_s  = {lambda_spiral:.5f} /s  "
          f"({'STABLE' if lambda_spiral < 0 else 'UNSTABLE'})")

    # --------------------------------------------------------
    # (J) LATERAL-DIRECTIONAL A-MATRIX
    # --------------------------------------------------------
    sub_header("J) Lateral-Directional State-Space Matrix A_lat-dir")

    Alat = np.array([
        [0.0,         1.0,         0.0,         0.0],
        [g_V * L_r2,  L_p2,        L_beta + g_V * Y_beta * L_r2,  -L_r1],
        [0.0,         0.0,         0.0,         1.0],
        [-g_V * N_r2, g_V,        -(N_beta + g_V * Y_beta * N_r2),  N_r1 + g_V * Y_beta],
    ])

    print("\n  A_lat-dir =")
    print(np.array2string(Alat, precision=4, suppress_small=True,
                          formatter={'float_kind': lambda x: f"{x:10.4f}"}))

    evals_lat = np.linalg.eigvals(Alat)
    print("\n  Eigenvalues of A_lat-dir:")
    for ev in evals_lat:
        print(f"    λ = {ev.real:+.5f}  {'+' if ev.imag >= 0 else '-'}  j{abs(ev.imag):.5f}")

    # --------------------------------------------------------
    # (K) CONTROL DERIVATIVES
    # --------------------------------------------------------
    sub_header("K) Control Effectiveness & Derivatives")

    tau_de = 0.5
    tau_dr = 0.5
    tau_da = 0.22

    CL_de   = eta * (S_HT / S_w) * CL_alpha_HT * tau_de
    Cm_de   = -eta * V_HT * CL_alpha_HT * tau_de

    y1, y2  = 0.5 * b_w, 0.9 * b_w
    integ   = ((y2**2 - y1**2) / 2
               + 2 * (y2**3 - y1**3) / (3 * b_w) * (lam - 1))
    Cl_da   = -2 * cl_alpha_w * tau_da / (S_w * b_w) * c_bar * integ

    Cn_da   = -(CD_alpha / cl_alpha_w) * Cl_da
    CY_da   = 0.0

    Cn_dr   = -V_VT * tau_dr * CL_alpha_VT
    CY_dr   = (S_VT / S_w) * tau_dr * CL_alpha_VT
    Cl_dr   = (S_VT / S_w) * tau_dr * CL_alpha_VT * (z_eff / b_w)

    print(f"  tau_de = {tau_de},  tau_dr = {tau_dr},  tau_da = {tau_da}")
    print(f"  CL_de  = {CL_de:.4f} /rad")
    print(f"  Cm_de  = {Cm_de:.4f} /rad")
    print(f"  Cl_da  = {Cl_da:.4f} /rad")
    print(f"  Cn_da  = {Cn_da:.4f} /rad")
    print(f"  Cn_dr  = {Cn_dr:.4f} /rad")
    print(f"  CY_dr  = {CY_dr:.4f} /rad")
    print(f"  Cl_dr  = {Cl_dr:.4f} /rad")

    # --------------------------------------------------------
    #  PLOT EIGENVALUES
    # --------------------------------------------------------
    def plot_eigenvalues(evals, title, filename):
        plt.figure(figsize=(8,6))
        real = [e.real for e in evals]
        imag = [e.imag for e in evals]
        plt.scatter(real, imag, color='red', s=60, zorder=5)
        # Add axes
        plt.axhline(0, color='black', linewidth=0.8)
        plt.axvline(0, color='black', linewidth=0.8)
        # Label each eigenvalue
        for i, ev in enumerate(evals):
            plt.annotate(f"{i+1}", (ev.real, ev.imag), textcoords="offset points",
                         xytext=(5,5), fontsize=9)
        plt.xlabel("Real part (1/s)")
        plt.ylabel("Imaginary part (1/s)")
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(filename, dpi=150)
        plt.show()

    plot_eigenvalues(evals_long, "Longitudinal Modes – Eigenvalues",
                     "longitudinal_eigenvalues.png")
    plot_eigenvalues(evals_lat, "Lateral-Directional Modes – Eigenvalues",
                     "lateral_eigenvalues.png")

    # --------------------------------------------------------
    #  SUMMARY TABLE
    # --------------------------------------------------------
    section_header("SUMMARY OF RESULTS")
    results = {
        "Trim AoA (deg)"            : alpha_trim_deg,
        "Trim Speed V* (m/s)"       : V_trim,
        "Trim CL"                   : CL_trim,
        "Trim CD"                   : CD_trim,
        "Thrust Required (N)"       : T_req,
        "Power Required (hp)"       : P_req / 745.7,
        "SP omega_n (rad/s) [approx]": wn_SP,
        "SP zeta        [approx]"   : zeta_SP,
        "Phugoid omega_n (rad/s)"   : wn_Ph,
        "Phugoid zeta"              : zeta_Ph,
        "Roll mode λ (1/s)"         : lambda_roll,
        "Dutch-roll omega_n (rad/s)": wn_DR,
        "Dutch-roll zeta"           : zeta_DR,
        "Spiral mode λ (1/s)"       : lambda_spiral,
    }
    for k, v in results.items():
        print(f"  {k:<40s}: {v:.5f}")

# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 60)
    print("  FLIGHT DYNAMICS & STABILITY ANALYSIS – UAV DESIGN")
    print("  Based on Group_11_Design_Report and Chapter 9 Appendix")
    print("=" * 60)

    uav_aircraft()

    print("\n" + "=" * 60)
    print("  ANALYSIS COMPLETE")
    print("  Eigenvalue plots saved as:")
    print("    - longitudinal_eigenvalues.png")
    print("    - lateral_eigenvalues.png")
    print("=" * 60)