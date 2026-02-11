import math

# =========================================================
# UAV Propeller + Motor Preliminary Sizing Tool
# Based on:
#  - Raymer conceptual design
#  - Advance ratio method
#  - Empirical propeller relations
# =========================================================

def size_prop_motor(
    mtow_kg,
    T_W,
    cruise_speed,
    battery_voltage,
    beta=1.6,          # pitch ratio (1.5–1.8 typical)
    rho=1.225,
    Ct=0.10,           # thrust coefficient at Jopt (0.08–0.12 typical)
    eta_prop=0.75
):

    g = 9.81

    # -----------------------------------------------------
    # 1. Required thrust
    # -----------------------------------------------------
    W = mtow_kg * g
    T_req = T_W * W

    # -----------------------------------------------------
    # 2. Optimal advance ratio (from your paper)
    #    Jopt = 0.2684β + 0.0879
    # -----------------------------------------------------
    Jopt = 0.2684 * beta + 0.0879

    # -----------------------------------------------------
    # 3. Prop diameter from:
    #    D = sqrt( T J² / (Ct ρ V²) )
    # -----------------------------------------------------
    D = math.sqrt((T_req * Jopt**2) / (Ct * rho * cruise_speed**2))

    # -----------------------------------------------------
    # 4. RPM
    # -----------------------------------------------------
    n = cruise_speed / (Jopt * D)   # rev/s
    rpm = 60 * n

    # -----------------------------------------------------
    # 5. Pitch from β = p/R
    # -----------------------------------------------------
    R = D / 2
    pitch = beta * R

    # -----------------------------------------------------
    # 6. Motor Kv
    #    RPM ≈ 0.8 Kv V
    # -----------------------------------------------------
    Kv = rpm / (0.8 * battery_voltage)

    # -----------------------------------------------------
    # 7. Power required
    # -----------------------------------------------------
    P_req = T_req * cruise_speed / eta_prop

    # -----------------------------------------------------
    # Unit conversions
    # -----------------------------------------------------
    D_in = D * 39.37
    pitch_in = pitch * 39.37

    # -----------------------------------------------------
    # Print results
    # -----------------------------------------------------
    print("\n================= RESULTS =================")
    print(f"Required thrust       : {T_req:.2f} N")
    print(f"Optimal advance ratio : {Jopt:.3f}")
    print(f"Diameter              : {D:.3f} m  ({D_in:.1f} in)")
    print(f"Pitch                 : {pitch:.3f} m  ({pitch_in:.1f} in)")
    print(f"Target RPM            : {rpm:.0f}")
    print(f"Motor Kv              : {Kv:.0f} RPM/V")
    print(f"Power required        : {P_req:.0f} W")
    print("===========================================\n")

    return {
        "Thrust_N": T_req,
        "Jopt": Jopt,
        "Diameter_m": D,
        "Diameter_in": D_in,
        "Pitch_m": pitch,
        "Pitch_in": pitch_in,
        "RPM": rpm,
        "Kv": Kv,
        "Power_W": P_req
    }


# =========================================================
# Example usage
# =========================================================

if __name__ == "__main__":

    size_prop_motor(
        mtow_kg=6.896,          # kg
        T_W=0.08468237455864157,            # thrust-to-weight
        cruise_speed=22,    # m/s
        battery_voltage= 10, # 14.8  # 4S LiPo
        beta=1.6
    )
