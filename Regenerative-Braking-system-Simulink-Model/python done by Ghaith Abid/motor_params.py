
import math
import numpy as np
from scipy.optimize import brentq

# PARAMETER SOLVER FUNCTIONS


def _id_mtpa_at_Imax(L_d, L_q, lambda_m, I_max):
    
    delta_L = L_q - L_d    # positive for IPMSM
    return (lambda_m - math.sqrt(lambda_m**2 + 8*delta_L**2*I_max**2)) \
           / (4 * delta_L)


def solve_consistent_params(lambda_m, xi, target_gap, I_max=580.0,Ld_min=100e-6, Ld_max=2000e-6):
    """
    Find Ld such that the FW arc width equals target_gap exactly.

    FW arc width = id_mtpa(I_MAX) - I_CH
                 = id_mtpa(I_MAX) - (-λm/Ld)
                 must be > 0 (id_mtpa less negative than I_CH)
    """
    def residual(Ld):
        Lq    = xi * Ld
        id_m  = _id_mtpa_at_Imax(Ld, Lq, lambda_m, I_max)
        I_ch  = -lambda_m / Ld
        return (id_m - I_ch) - target_gap

    Ld_sol = brentq(residual, Ld_min, Ld_max)
    return Ld_sol, xi * Ld_sol


def verify_params(L_d, L_q, lambda_m, I_max, T_max, P_max, omega_base,p_poles, verbose=True):
    K_T      = 1.5 * p_poles
    I_ch     = -lambda_m / L_d
    xi       = L_q / L_d
    id_m     = _id_mtpa_at_Imax(L_d, L_q, lambda_m, I_max)
    iq_m     = math.sqrt(max(I_max**2 - id_m**2, 0.0))
    T_ach    = K_T * (lambda_m * iq_m + (L_d - L_q) * id_m * iq_m)
    gap      = id_m - I_ch
    ratio    = abs(I_ch) / I_max
    P_at_wb  = T_max * omega_base / 1e3

    if verbose:
        print("=" * 52)
        print("  PARAMETER VERIFICATION REPORT")
        print("=" * 52)
        print(f"  Ld              = {L_d*1e6:.1f} µH")
        print(f"  Lq              = {L_q*1e6:.1f} µH")
        print(f"  λm              = {lambda_m*1e3:.1f} mWb")
        print(f"  ξ  (Lq/Ld)     = {xi:.2f}   (target 2.0–4.0)")
        print(f"  I_CH            = {I_ch:.1f} A")
        print(f"  id_mtpa(I_MAX)  = {id_m:.1f} A")
        print(f"  FW arc width    = {gap:.1f} A   (must be > 0)")
        print(f"  |I_CH|/I_MAX    = {ratio:.3f}  (target 0.45–0.75)")
        print(f"  T_achieved      = {T_ach:.1f} Nm  (need ≥ {T_max:.0f} Nm)")
        print(f"  P @ ω_base      = {P_at_wb:.1f} kW  (target ≈ 200 kW)")

    passed = True
    checks = {
        "FW arc > 0"          : gap > 0,
        "|I_CH|/I_MAX in range": 0.45 <= ratio <= 0.75,
        "T_achieved >= T_max" : T_ach >= T_max,
        "ξ in range"          : 2.0 <= xi <= 4.0,
    }
    for name, ok in checks.items():
        status = "✓" if ok else "✗  FAIL"
        if verbose:
            print(f"  [{status}]  {name}")
        if not ok:
            passed = False

    if verbose:
        print("=" * 52)
        print(f"  Overall: {'PASSED' if passed else 'FAILED'}")
        print("=" * 52)

    return passed


# 1. Vehicle & Drivetrain
MOTOR_NAME       = "Mercedes CLA 250+ IPMSM"
P_max            = 200e3      # [W]   Peak mechanical power
T_max            = 335.0      # [Nm]  Peak torque
GEAR_RATIO_1     = 11.0       # [-]   First gear
GEAR_RATIO_2     = 5.0        # [-]   Second gear
WHEEL_RADIUS     = 0.330      # [m]   Tyre rolling radius

# 2. DC Bus
V_DC             = 706.0                   # [V]  Nominal HV battery voltage
V_MAX            = V_DC / math.sqrt(3)     # [V]  ≈ 407.6 V  (space-vector limit)

# 3. Pole Pairs
P_POLES          = 4          # [-]   Estimated 

# 4. Electrical Parameters — solved for consistency solve_consistent_params(λm=0.095, ξ=2.4, target_gap=50 A)
R_S              = 0.012      # [Ω]   Stator resistance (unchanged)

_LAMBDA_M_DESIGN = 0.095      # [Wb]  design choice
_XI_DESIGN       = 2.4        # [-]   saliency design choice
_TARGET_GAP      = 50.0       # [A]   desired FW arc width

L_D, L_Q = solve_consistent_params(_LAMBDA_M_DESIGN, _XI_DESIGN, _TARGET_GAP)
LAMBDA_M         = _LAMBDA_M_DESIGN

# 5. Current Limit
I_MAX            = 580.0      # [A]   Peak phase current

# 6. Speed Limits

OMEGA_BASE_RPM   = 5500.0     # [RPM]  ω_base ≈ P_max/T_max = 597 rad/s
OMEGA_MAX_RPM    = 14000.0    # [RPM]  back-calculated from top speed + gear ratio


# 7. Mechanical

J_ROTOR          = 0.025      # [kg·m²]  Rotor inertia (estimated)


# 8. Derived Constants 
K_T              = 1.5 * P_POLES                   # [-]     = 6.0
OMEGA_BASE       = OMEGA_BASE_RPM * 2*math.pi/60   # [rad/s] mechanical base speed
OMEGA_MAX        = OMEGA_MAX_RPM  * 2*math.pi/60   # [rad/s] mechanical max speed
I_CH             = -LAMBDA_M / L_D                 # [A]     characteristic current
SALIENCY         = L_Q / L_D                       # [-]     = ξ

if __name__ == "__main__":
    verify_params(L_D, L_Q, LAMBDA_M, I_MAX, T_max, P_max, OMEGA_BASE,
                  P_POLES, verbose=True)
    print(f"\n  L_D  = {L_D*1e6:.2f} µH  (use this in all modules)")
    print(f"  L_Q  = {L_Q*1e6:.2f} µH  (use this in all modules)")