import math
import numpy as np


# =========================================================
# Characteristic parameters (paper equations)
# =========================================================

def thrust_characteristics(beta):
    J0T = 0.4559*beta + 0.1574
    dCT = 0.0205*beta - 0.0073
    JmT = 0.2275*beta + 0.0792
    CT0 = -0.034*beta**2 + 0.13306*beta - 0.00287
    return CT0, dCT, JmT, J0T


def power_characteristics(beta):
    J0P = 0.3839*beta + 0.3553
    dCP = 0.0249*beta - 0.0117
    JmP = 0.2670*beta + 0.1192
    CP0 = 0.017716*beta**2 + 0.0064*beta + 0.014084
    return CP0, dCP, JmP, J0P


# =========================================================
# Polynomial coefficients
# =========================================================

def thrust_poly_coeffs(beta):

    CT0, dCT, JmT, J0T = thrust_characteristics(beta)

    bT2 = -dCT/(JmT**2)
    bT1 = dCT/JmT - CT0/J0T
    bT0 = CT0

    return bT2, bT1, bT0


def power_poly_coeffs(beta):

    CP0, dCP, JmP, J0P = power_characteristics(beta)

    A = np.array([
        [J0P**3, J0P**2, J0P, 1],
        [JmP**3, JmP**2, JmP, 0],
        [3*JmP**2, 2*JmP, 1, 0],
        [0, 0, 0, 1]
    ])

    B = np.array([
        0,
        dCP - CP0*JmP/J0P,
        -CP0/J0P,
        CP0
    ])

    return np.linalg.solve(A, B)


# =========================================================
# Main sizing tool
# =========================================================

def size_prop_motor(mtow_kg, T_W, cruise_speed, battery_voltage, beta, rho=1.225):

    g = 9.81
    W = mtow_kg * g
    Treq = T_W * W
    V = cruise_speed

    # Optimal advance ratio (paper Eq. 22)
    Jopt = 0.2684*beta + 0.0879

    # Polynomial coefficients
    bT2, bT1, bT0 = thrust_poly_coeffs(beta)
    bP3, bP2, bP1, bP0 = power_poly_coeffs(beta)

    # Evaluate coefficients at Jopt
    CT = bT2*Jopt**2 + bT1*Jopt + bT0
    CP = bP3*Jopt**3 + bP2*Jopt**2 + bP1*Jopt + bP0

    # Diameter (paper Eq. 25)
    D = math.sqrt(Treq * Jopt**2 / (CT * rho * V**2))

    # RPM
    n = V/(Jopt*D)
    rpm = 60*n

    # Power (paper Eq. 23)
    P = CP * rho * n**3 * D**5

    Kv = rpm/(0.8*battery_voltage)

    print("\n===== PAPER-CORRECT RESULTS =====")
    print(f"CT(Jopt) = {CT:.4f}")
    print(f"CP(Jopt) = {CP:.4f}")
    print(f"Diameter = {D*39.37:.1f} in")
    print(f"RPM      = {rpm:.0f}")
    print(f"Power    = {P:.0f} W")
    print(f"Kv       = {Kv:.0f}")
    print("================================")

    return D, rpm, P

if __name__ == "__main__":

    size_prop_motor(
        mtow_kg=6.896,          # kg
        T_W=0.08468237455864157,            # thrust-to-weight
        cruise_speed=22,    # m/s
        battery_voltage= 14.8, # 14.8  # 4S LiPo
        beta=1.6
    )
