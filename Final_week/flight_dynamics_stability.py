"""
============================================================
  Flight Dynamics & Stability Analysis
  Based on: Chapter 9 Appendix – Case Studies
  Case 1 : UAV aircraft.
============================================================
All formulas are taken directly from the supplied PDF.
No external formulas have been introduced.
"""

import numpy as np

# ─────────────────────────────────────────────────────────
#  UTILITY
# ─────────────────────────────────────────────────────────
D2R = np.pi / 180.0   # degrees → radians


def section_header(title):
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def sub_header(title):
    print(f"\n--- {title} ---")


# ╔══════════════════════════════════════════════════════════╗
# ║          CASE 1 – GA AIRCRAFT                           ║
# ╚══════════════════════════════════════════════════════════╝

def ga_aircraft():
    section_header("CASE 1 – UAV")

    # ─────────────────────────────────────────────────────
    # (A) INPUT PARAMETERS
    # ─────────────────────────────────────────────────────
    sub_header("A) Input Parameters")

    # ── Geometry ──────────────────────────────────────────
    g        = 9.81          # m/s^2  – gravitational acceleration
    S_w      = 0.5307         # m^2   – wing planform (reference) area
    b_w      = 2         # m     – wing span
    c_r      = 0.2653          # m     – wing root chord (NACA 23018)
    c_t      = 0.2650          # m     – wing tip chord  (NACA 23012)
    lam      = c_t / c_r    # taper ratio  (= 0.5)
    Lambda_c4_deg = 0     # deg   – wing quarter-chord sweep angle
    Gamma_deg     = 0     # deg   – wing dihedral angle
    i_w      = 3           # deg   – wing incidence setting angle

    S_HT     = 0.1327          # m^2   – horizontal tail planform area
    AR_HT    = 4.5           # –     – HT aspect ratio
    i_t      = -0.423           # deg   – HT incidence setting angle
    X_AC_HT  = 0.964          # m     – HT aerodynamic-centre from nose

    S_VT     = 0.0478          # m^2   – vertical tail planform area
    AR_VT    = 1.5           # –     – VT aspect ratio
    X_AC_VT  = 0.964          # m     – VT AC from nose

    l_f      = 1.12          # m     – fuselage length
    S_f      = 1.11          # m^2   – fuselage side projected area
    D_f      = 2*(np.sqrt(0.3*0.25*1.2/np.pi))          # m     – max fuselage diameter

    S_da     = 0.05          # m^2   – aileron surface area
    l_da     = 0.76          # m     – aileron span
    S_de     = 0.04          # m^2   – elevator surface area
    S_dr     = 0.01          # m^2   – rudder surface area

    # ── Mass & Inertia ────────────────────────────────────
    m        = 6.298       # kg    – aircraft mass
    I_xx     = 0.751       # kg·m² – roll inertia
    I_yy     = 0.384       # kg·m² – pitch inertia
    I_zz     = 1.045       # kg·m² – yaw inertia

    X_CG     = 0.264          # m     – CG location from nose

    # ── Aerodynamic parameters ───────────────────────────
    alpha_L0_deg = -3.5      # deg   – wing zero-lift AoA
    Cm_ac_w      = -0.022    # –     – residual pitching moment at AC
    SM           = 0.14      # –     – assumed static margin (15 %)
    e_oswald     = 1       # –     – Oswald efficiency factor
    C_D0         = 0.02873     # –     – zero-lift drag coefficient
    # Airfoil lift-curve slopes (rad⁻¹)
    cl_alpha_w   = 3.6107      # /rad  – NACA 23018 airfoil (0.1/deg)
    cl_alpha_HT  = 6.24      # /rad  – NACA 0012 (0.109/deg)
    cl_alpha_VT  = 6.24      # /rad  – NACA 0012

    # ── Flight conditions ─────────────────────────────────
    rho          = 1.22      # kg/m³ – air density at flight altitude
    a_sound      = 340.0     # m/s   – speed of sound
    eta          = 0.99      # –     – dynamic pressure ratio at tail

    # ── Engine/propeller data ─────────────────────────────
    P_max        = 450  # N·m/s – max engine power (360 hp)
    N_blades     = 4         # –     – number of propeller blades
    D_prop       = 0.3048      # m     – propeller diameter
    l_prop       = 0.12957      # m     – prop disc to CG distance

    # vertical location of wing root q-chord above CG line (for lat-dir)
    z_w          = 0.06     # m

    print("  Inputs loaded.")

    # ─────────────────────────────────────────────────────
    # (B) DERIVED GEOMETRY PARAMETERS
    # ─────────────────────────────────────────────────────
    sub_header("B) Derived Geometry Parameters")

    # Mean aerodynamic chord  (Eq. from PDF)
    c_bar = (2/3) * c_r * (1 + lam + lam**2) / (1 + lam)

    # Wing aspect ratio
    AR_w = b_w**2 / S_w

    # Tail spans
    b_HT = np.sqrt(AR_HT * S_HT)
    b_VT = np.sqrt(AR_VT * S_VT)

    # VT aerodynamic centre height above fuselage reference line
    z_VT = (4/9) * b_VT

    # Propeller disc area
    S_prop = np.pi * D_prop**2 / 4

    # Moment arm lengths
    l_HT = X_AC_HT - X_CG
    l_VT = X_AC_VT - X_CG

    # Tail volume ratios
    V_HT = (S_HT * l_HT) / (S_w * c_bar)
    V_VT = (S_VT * l_VT) / (S_w * b_w)

    # Non-dimensionalised CG location
    h_CG = X_CG / c_bar

    # Induced-drag factor
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

    # ─────────────────────────────────────────────────────
    # (C) AERODYNAMIC COEFFICIENT CALCULATIONS
    # ─────────────────────────────────────────────────────
    sub_header("C) Aerodynamic Coefficients")

    # Finite-wing lift-curve slopes  (low-subsonic simplified form, PDF Eq.)
    CL_alpha_w  = cl_alpha_w  / (1 + cl_alpha_w  / (np.pi * AR_w))
    CL_alpha_HT = cl_alpha_HT / (1 + cl_alpha_HT / (np.pi * AR_HT))
    CL_alpha_VT = cl_alpha_VT / (1 + cl_alpha_VT / (np.pi * AR_VT))

    print(f"  CL_alpha_w   = {CL_alpha_w:.3f} /rad")
    print(f"  CL_alpha_HT  = {CL_alpha_HT:.3f} /rad")
    print(f"  CL_alpha_VT  = {CL_alpha_VT:.3f} /rad")

    # Wing lift coefficient at zero AoA
    CL0_w = CL_alpha_w * (alpha_L0_deg * D2R) * (-1)   # note: alpha_L0 is negative
    # PDF formula: CL0_w = CL_alpha_w × |alpha_L=0| × pi/180
    CL0_w = CL_alpha_w * abs(alpha_L0_deg) * D2R
    print(f"  CL0_w        = {CL0_w:.4f}")

    # Downwash gradient  de/da
    eps_alpha = 2 * CL_alpha_w / (np.pi * AR_w)
    eps_0     = 2 * CL0_w / (np.pi * AR_w)
    print(f"  de/da        = {eps_alpha:.4f}")
    print(f"  eps_0        = {eps_0:.4f}")

    # ─────────────────────────────────────────────────────
    # (D) TRIM CALCULATIONS
    # ─────────────────────────────────────────────────────
    sub_header("D) Trim Calculations")

    # Neutral-point location (non-dim)
    h_NP = h_CG + SM
    X_NP = h_NP * c_bar

    # Wing-body AC location (from NP and tail contribution)
    h_AC_wb = h_NP - V_HT * (CL_alpha_HT / CL_alpha_w) * (1 - eps_alpha)
    X_AC_wb = h_AC_wb * c_bar

    print(f"  h_NP         = {h_NP:.3f},   X_NP   = {X_NP:.3f} m")
    print(f"  h_AC_wb      = {h_AC_wb:.3f},  X_AC_wb = {X_AC_wb:.3f} m")

    # Pitching moment derivatives
    Cm_alpha_w  = ((X_CG - X_AC_wb) / c_bar) * CL_alpha_w
    Cm_alpha_HT = -eta * V_HT * CL_alpha_HT * (1 - eps_alpha)
    Cm_alpha    = Cm_alpha_w + Cm_alpha_HT

    Cm0_HT = eta * V_HT * CL_alpha_HT * (i_t * D2R * (-1) + eps_0)
    Cm0    = Cm_ac_w + eta * V_HT * CL_alpha_HT * (eps_0)

    # Trim AoA (zero elevator deflection)
    alpha_trim_rad = -Cm0 / Cm_alpha
    alpha_trim_deg = alpha_trim_rad / D2R

    # Whole-aircraft CL_alpha and CL0
    CL0 = (CL_alpha_w * abs(alpha_L0_deg) * D2R
           + eta * (S_HT / S_w) * CL_alpha_HT * (i_t * D2R - eps_0))
    CL_alpha = (CL_alpha_w
                + eta * (S_HT / S_w) * CL_alpha_HT * (1 - eps_alpha))

    # Trim lift coefficient
    CL_trim = (CL_alpha_w * (alpha_trim_deg - alpha_L0_deg) * D2R
               + eta * (S_HT / S_w) * CL_alpha_HT
               * (alpha_trim_rad * (1 - eps_alpha) - eps_0 + i_t * D2R))

    # Trim speed from L = W
    V_trim = np.sqrt(2 / rho * (m * g / S_w) * (1 / CL_trim))

    # Dynamic pressure at trim
    q_trim = 0.5 * rho * V_trim**2

    # Mach number
    Ma_trim = V_trim / a_sound

    # Trim drag coefficient
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

    # ─────────────────────────────────────────────────────
    # (E) LONGITUDINAL STABILITY DERIVATIVES
    # ─────────────────────────────────────────────────────
    sub_header("E) Longitudinal Stability Derivatives")

    # Pitch-rate & alpha-dot derivatives
    Cm_q1 = -2 * eta * V_HT * CL_alpha_HT * (l_HT / c_bar)
    CL_q1 = 2 * eta * V_HT * CL_alpha_HT
    # alpha-dot term (PDF sets Cm_alpha_dot & CL_alpha_dot = 0 in this example)
    Cm_adot = 0.0
    CL_adot = 0.0

    # Drag derivatives
    CD_alpha = 2 * K * CL_trim * CL_alpha
    CD_q1    = 2 * K * CL_trim * CL_q1

    # Mach-number derivatives (compressibility)
    CL_Ma = CL_trim * Ma_trim / (1 - Ma_trim**2)
    CD_Ma = 0.0     # PDF sets this to zero for this aircraft

    # Useful ratios
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

    # ─────────────────────────────────────────────────────
    # (F) LONGITUDINAL DYNAMIC MODES (2nd-order approx.)
    # ─────────────────────────────────────────────────────
    sub_header("F) Longitudinal Dynamic Modes (Approximate)")

    # ── Short Period ──────────────────────────────────────
    # (omega_n)^2_SP = -(q̄*Sc/Iyy) * Cm_alpha
    wn_SP_sq = -qSc_Iyy * Cm_alpha
    wn_SP    = np.sqrt(wn_SP_sq)
    T_SP     = 2 * np.pi / wn_SP

    # 2*zeta*wn_SP = -(q̄*Sc/Iyy) * (c/2V*) * (Cm_q1 + Cm_adot)
    two_zeta_wn_SP = -qSc_Iyy * c_2V * (Cm_q1 + Cm_adot)
    zeta_SP  = two_zeta_wn_SP / (2 * wn_SP)

    # ── Phugoid ───────────────────────────────────────────
    # (omega_n)^2_Ph = (g/V*)^2 * [(q̄*S/W) * (Ma* * CL_Ma + 2*CL*)]
    wn_Ph_sq = (g_V**2) * qS_W * (Ma_trim * CL_Ma + 2 * CL_trim)
    wn_Ph    = np.sqrt(wn_Ph_sq)
    T_Ph     = 2 * np.pi / wn_Ph

    # 2*zeta*wn_Ph = (g/V*) * [(q̄*S/W) * (Ma* * CD_Ma + 2*CD*)]
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

    # ─────────────────────────────────────────────────────
    # (G) LONGITUDINAL A-MATRIX (1st-order form, Sec. 9.1.2)
    # ─────────────────────────────────────────────────────
    sub_header("G) Longitudinal State-Space Matrix A_long")

    # State vector: [ΔV/V*, Δγ, Δα, Δα̇]
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

    # k-factors for the α̈ row
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

    # Identify modes from eigenvalues
    for ev in evals_long:
        wn = abs(ev)
        if wn > 0:
            zeta_ev = -ev.real / wn if ev.imag != 0 else None
            if abs(ev.imag) > 0.05:
                label = "Short Period" if wn > 1.0 else "Phugoid"
                print(f"    → {label}: omega_n = {wn:.4f} rad/s, "
                      f"zeta = {(-ev.real/wn):.4f}")

    # ─────────────────────────────────────────────────────
    #  LATERAL-DIRECTIONAL AERODYNAMICS (Sec. 9.1.3)
    # ─────────────────────────────────────────────────────
    sub_header("H) Lateral-Directional Aerodynamic Coefficients")

    Gamma_rad    = Gamma_deg * D2R
    alpha_t_rad  = alpha_trim_rad

    # ── Sidewash factor ────────────────────────────────────
    # (1 + dsig/dbeta) * eta_VT = formula from PDF
    Lambda_c4_rad = Lambda_c4_deg * D2R
    sidewash = (0.724
                + 3.06 * (S_VT / S_w) / (2.0)
                + 0.4 * z_w / D_f
                + 0.009 * AR_w)
    # Note: PDF uses (S_VT/S_w)/(1+cos(Λ_c/4)) with cos≈1 for small sweep
    # Reproduced exactly as in PDF:
    sidewash = (0.724
                + 3.06 * (S_VT / S_w) / (1 + np.cos(Lambda_c4_rad))
                + 0.4 * z_w / D_f
                + 0.009 * AR_w)

    # ── Side-force coefficients ──────────────────────────
    CY_beta_w  = -0.0001 * Gamma_deg * D2R   # ≈ 0
    CY_beta_VT = (-eta * CL_alpha_VT * (S_VT / S_w) * sidewash)
    CY_beta    = CY_beta_w + CY_beta_VT

    # ── Yawing-moment coefficients ────────────────────────
    Cn_beta_w  = -0.075 * Gamma_deg * D2R * CL_trim
    Cn_beta_VT = -CY_beta_VT * (l_VT / b_w)
    Cn_beta    = Cn_beta_w + Cn_beta_VT

    # ── Rolling-moment coefficients ───────────────────────
    # z = z_VT*cos(α) – l_VT*sin(α)
    z_eff = z_VT * np.cos(alpha_t_rad) - l_VT * np.sin(alpha_t_rad)

    Cl_beta_w  = -(Gamma_rad * CL_alpha_w * c_r * b_w) / (6 * S_w)
    Cl_beta_VT = (z_eff / b_w) * CY_beta_VT
    Cl_beta    = Cl_beta_w + Cl_beta_VT

    print(f"  Sidewash factor (1+dsig/db)*eta_VT = {sidewash:.4f}")
    print(f"  CY_beta  = {CY_beta:.4f} /rad")
    print(f"  Cn_beta  = {Cn_beta:.4f} /rad")
    print(f"  Cl_beta  = {Cl_beta:.4f} /rad")

    # ── Roll-rate derivatives ─────────────────────────────
    CY_p2_VT = (2 / b_w) * (z_eff - z_VT) * CY_beta_VT
    CY_p2    = CY_p2_VT

    Cl_p2_VT = 2 * (z_eff / b_w) * ((z_eff - z_VT) / b_w) * CY_beta_VT
    Cl_p2_w  = -(1/6) * (CL_alpha + CD_trim)    # CL_alpha here is whole-ac
    Cl_p2    = Cl_p2_VT + Cl_p2_w

    lVT_eff = l_VT * np.cos(alpha_t_rad) + z_VT * np.sin(alpha_t_rad)
    Cn_p2_w  = -(1/6) * (CL_trim - CD_alpha)
    Cn_p2_VT = (-(2 / b_w) * lVT_eff
                * ((z_eff - z_VT) / b_w) * CY_beta_VT)
    Cn_p2    = Cn_p2_w + Cn_p2_VT

    print(f"\n  CY_p2    = {CY_p2:.4f} /rad")
    print(f"  Cl_p2    = {Cl_p2:.4f} /rad")
    print(f"  Cn_p2    = {Cn_p2:.4f} /rad")

    # ── Yaw-rate derivatives ──────────────────────────────
    CY_r1_VT = -(2 / b_w) * lVT_eff * CY_beta_VT
    Cl_r1_w  = CL_trim / 3
    Cl_r1_VT = (-(2 / b_w**2) * lVT_eff * z_eff * CY_beta_VT)
    Cl_r1    = Cl_r1_w + Cl_r1_VT
    Cn_r1_VT = (2 / b_w**2) * lVT_eff**2 * CY_beta_VT
    Cn_r1    = Cn_r1_VT

    # Secondary yaw/roll-rate terms
    Cn_r2    = -CD_trim / 3
    Cl_r2    =  CL_trim / 3

    print(f"  CY_r1    = {CY_r1_VT:.4f} /rad")
    print(f"  Cl_r1    = {Cl_r1:.4f} /rad")
    print(f"  Cn_r1    = {Cn_r1:.4f} /rad")
    print(f"  Cn_r2    = {Cn_r2:.4f} /rad")
    print(f"  Cl_r2    = {Cl_r2:.4f} /rad")

    # ─────────────────────────────────────────────────────
    # (I) LATERAL-DIRECTIONAL DYNAMIC MODES (2nd-order)
    # ─────────────────────────────────────────────────────
    sub_header("I) Lateral-Directional Dynamic Modes (Approximate)")

    qSb_Ixx = q_trim * S_w * b_w / I_xx
    qSb_Izz = q_trim * S_w * b_w / I_zz
    qS_m    = q_trim * S_w / m
    b_2V    = b_w / (2 * V_trim)

    # Dimensional derivatives
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

    # Roll mode eigenvalue  (λ_r ≈ L_p2)
    lambda_roll = L_p2

    # Dutch roll mode
    wn_DR_sq = (N_beta
                + g_V * (Y_beta * N_r2 + L_beta / L_p2))
    wn_DR    = np.sqrt(abs(wn_DR_sq))
    two_zeta_wn_DR = (-N_r1
                      - g_V * (Y_beta + L_r1 / L_p2))
    zeta_DR  = two_zeta_wn_DR / (2 * wn_DR)

    # Spiral mode eigenvalue
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

    # ─────────────────────────────────────────────────────
    # (J) LATERAL-DIRECTIONAL A-MATRIX (Sec. 9.1.4)
    # ─────────────────────────────────────────────────────
    sub_header("J) Lateral-Directional State-Space Matrix A_lat-dir")

    # State vector: [Δμ, Δμ̇, Δβ, Δβ̇]
    # from PDF Eq. 9.11
    Alat = np.array([
        [0.0,
         1.0,
         0.0,
         0.0],
        [g_V * L_r2,
         L_p2,
         L_beta + g_V * Y_beta * L_r2,
         -L_r1],
        [0.0,
         0.0,
         0.0,
         1.0],
        [-g_V * N_r2,
         g_V,
         -(N_beta + g_V * Y_beta * N_r2),
         N_r1 + g_V * Y_beta],
    ])

    print("\n  A_lat-dir =")
    print(np.array2string(Alat, precision=4, suppress_small=True,
                          formatter={'float_kind': lambda x: f"{x:10.4f}"}))

    evals_lat = np.linalg.eigvals(Alat)
    print("\n  Eigenvalues of A_lat-dir:")
    for ev in evals_lat:
        print(f"    λ = {ev.real:+.5f}  {'+' if ev.imag >= 0 else '-'}  j{abs(ev.imag):.5f}")

    # ─────────────────────────────────────────────────────
    # (K) CONTROL DERIVATIVES (Sec. 9.1.4 end)
    # ─────────────────────────────────────────────────────
    sub_header("K) Control Effectiveness & Derivatives")

    # Control effectiveness parameters (from PDF chart Fig 7.9)
    tau_de = 0.5    # elevator effectiveness
    tau_dr = 0.5    # rudder effectiveness
    tau_da = 0.22   # aileron effectiveness

    CL_de   = eta * (S_HT / S_w) * CL_alpha_HT * tau_de
    Cm_de   = -eta * V_HT * CL_alpha_HT * tau_de

    # Aileron (running 50-90% span)
    # Integral: ∫c(y)*y dy from 0.5b to 0.9b with linear taper
    y1, y2  = 0.5 * b_w, 0.9 * b_w
    integ   = ((y2**2 - y1**2) / 2
               + 2 * (y2**3 - y1**3) / (3 * b_w) * (lam - 1))
    Cl_da   = -2 * cl_alpha_w * tau_da / (S_w * b_w) * c_r * integ

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

    # ─────────────────────────────────────────────────────
    # SUMMARY TABLE
    # ─────────────────────────────────────────────────────
    section_header("AIRCRAFT – SUMMARY OF RESULTS")
    results = {
        "Trim AoA (deg)"            : alpha_trim_deg,
        "Trim Speed V* (m/s)"       : V_trim,
        "Trim CL"                   : CL_trim,
        "Trim CD"                   : CD_trim,
        "Thrust Required (N)"       : T_req,
        "Power Required (hp)"       : P_req / 745.7,
        "SP omega_n (rad/s) [approx]": wn_SP,
        "SP zeta        [approx]"   : zeta_SP,
        "SP omega_n (rad/s) [eigen]": abs(evals_long[np.argmax(
                                        [abs(e.imag) for e in evals_long])]),
        "Phugoid omega_n (rad/s)"   : wn_Ph,
        "Phugoid zeta"              : zeta_Ph,
        "Roll mode λ (1/s)"         : lambda_roll,
        "Dutch-roll omega_n (rad/s)": wn_DR,
        "Dutch-roll zeta"           : zeta_DR,
        "Spiral mode λ (1/s)"       : lambda_spiral,
    }
    for k, v in results.items():
        print(f"  {k:<40s}: {v:.5f}")




# ──────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  FLIGHT DYNAMICS & STABILITY ANALYSIS")
    print("  Chapter 9 Appendix – Case Studies")
    print("=" * 60)

    ga_aircraft()
  
    print("\n" + "=" * 60)
    print("  ANALYSIS COMPLETE")
    print("=" * 60)
