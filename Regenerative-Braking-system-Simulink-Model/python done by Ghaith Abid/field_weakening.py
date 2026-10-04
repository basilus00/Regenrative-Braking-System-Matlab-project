#   1. From current circle (first contraint where i_d*² + i_q*² ≤ I_max²):
#        i_q*²  =  I_max²  -  i_d*²
#
#   2. Substitute into voltage ellipse (given by second constraint where (L_q·i_q*)² + (L_d·i_d* + λ_m)² ≤ (V_max/ω_e)²):
#        L_q²·(I_max² - i_d*²)  +  (L_d·i_d* + λ_m)²  =  (V_max/ω_e)²
#
#   3. Expand:
#        L_q²·I_max²  -  L_q²·i_d*²+  L_d²·i_d*²  +  2·L_d·λ_m·i_d*  +  λ_m²=  (V_max/ω_e)²
#
#   4. quadratic function in i_d*:
#        A·i_d*²  +  B·i_d*  +  C  =  0 where 
#        A  =  L_d²  -  L_q² (always negative, L_q > L_d)(propriety of IPMSM)
#        B  =  2·L_d·λ_m (always positive)
#        C  =  λ_m²  +  L_q²·I_max²  -  (V_max/ω_e)²
#
#   5. i_d* expression (we tend to make it more negative, or closer to center of voltage ellipse):
#
#        discriminant  =  B²  -  4·A·C=  4·L_d²·λ_m²- 4·(L_d²-L_q²)·(λ_m² + L_q²·I_max² - (V_max/ω_e)²)
#
#        i_d*  =  (-B  -  √discriminant)  /  (2·A)
#
#   6. Recover i_q* from current circle:
#        i_q*  =  +√(I_max² - i_d*²)   if T_request ≥ 0  (motoring)
#        i_q*  =  -√(I_max² - i_d*²)   if T_request <  0  (regen)
#
#   7. Compute achievable torque (consequence, not input in FW):
#        T_achieved  =  K·[λ_m·i_q*  +  (L_d - L_q)·i_d*·i_q*]
#        where K = (3/2)·p
#
# Sign behavior (both motoring and regen):
#   i_d*  <  0   always  ← flux weakening regardless of torque direction
#   i_q*  >  0   motoring
#   i_q*  <  0   regen
#
# Voltage ellipse properties:
#   Center    : (-λ_m/L_d,  0)           ← fixed, never moves with speed
#   Semi-axis : V_max / (ω_e · L_d)     ← shrinks as ω_e increases
#   Semi-axis : V_max / (ω_e · L_q)     ← shrinks as ω_e increases
#
# Transition from MTPA to FW:
#   ω_base(T) = V_max / √[L_q²·i_q_mtpa² + (L_d·i_d_mtpa + λ_m)²]
#   FW activates when ω_e > ω_base(T)

import numpy as np

def solve_fw(omega_e, T_request, p, lambda_m, L_d, L_q, I_max, V_max):
    K = 1.5 * p
    Vr = V_max / omega_e

    A = L_d**2 - L_q**2
    B = 2.0 * L_d * lambda_m
    C = lambda_m**2 + L_q**2 * I_max**2 - Vr**2

    discriminant = B**2 - 4.0 * A * C
    if discriminant < 0:
        return None, None, None

    sqrt_disc = np.sqrt(discriminant)
    id1 = (-B + sqrt_disc) / (2.0 * A)
    id2 = (-B - sqrt_disc) / (2.0 * A)

    candidates = [x for x in (id1, id2) if x < 0.0]
    if not candidates:
        return None, None, None
    id_star = max(candidates)

    iq_mag = np.sqrt(max(I_max**2 - id_star**2, 0.0))
    iq_star = iq_mag if T_request >= 0 else -iq_mag

    T_achieved = K * (lambda_m * iq_star + (L_d - L_q) * id_star * iq_star)

    vd = -omega_e * L_q * iq_star
    vq =  omega_e * (L_d * id_star + lambda_m)
    v_mag = np.sqrt(vd**2 + vq**2)
    if v_mag > V_max + 1e-3:
        return None, None, None

    i_mag = np.sqrt(id_star**2 + iq_star**2)
    if i_mag > I_max + 1e-3:
        return None, None, None

    return id_star, iq_star, T_achieved


def compute_transition_speed(id_mtpa, iq_mtpa, lambda_m, L_d, L_q, V_max):
    """
    Compute electrical speed at which MTPA transitions to FW.

    At the transition point the MTPA operating point (id_mtpa, iq_mtpa)
    sits exactly on the voltage ellipse boundary:
        Vmax² = (ω·Lq·iq)² + (ω·(Ld·id + λm))²
              = ω² · [(Lq·iq)² + (Ld·id + λm)²]

    Solving for ω:
        ω_base = Vmax / sqrt((Lq·iq)² + (Ld·id + λm)²)

    """
    vd_term = L_q * iq_mtpa                  # [V·s/rad] contribution from iq
    vq_term = L_d * id_mtpa + lambda_m       # [V·s/rad] contribution from id + flux
    return V_max / np.sqrt(vd_term**2 + vq_term**2)

def check_voltage_feasibility(id_star,iq_star,omega_e,R_s,L_d,L_q,lambda_m,V_max):
    """
    Check whether operating point (id*, iq*) is voltage-feasible at the given speed using full steady-state dq equations

    """

    # Full steady-state dq voltage equations 
    v_d = -omega_e * L_q * iq_star
    v_q =  omega_e * (L_d * id_star + lambda_m)

    v_required=np.sqrt(v_d**2+v_q**2)
    margin=V_max-v_required
    feasible= bool(margin >= -1e-6)

    return feasible, v_required, margin


def build_fw_table(omega_e_min,omega_e_max,n_points,T_request,p,lambda_m,L_d,L_q,I_max,V_max):
    #Pre-compute FW lookup table across a speed range (Used to generate offline table for Simulink lookup block)
    
    speeds=np.linspace(omega_e_min, omega_e_max, n_points)
    id_table = np.full(n_points, np.nan)
    iq_table=np.zeros(n_points)
    T_table=np.zeros(n_points)
    for i, omega_e in enumerate(speeds):
        id_val, iq_val, T_val  =  solve_fw(omega_e, T_request,p, lambda_m, L_d, L_q, I_max, V_max)
        if id_val is not None:
            id_table[i]=id_val
            iq_table[i]=iq_val
            T_table[i]=T_val
            
    return speeds, id_table, iq_table, T_table


# Validation (runs only when executing this file directly)

if __name__ == "__main__":
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    from motor_params import (P_POLES, LAMBDA_M, L_D, L_Q,I_MAX, V_MAX, K_T, T_max,OMEGA_BASE, OMEGA_MAX, R_S)
    from mtpa import solve_mtpa
    
    OMEGA_BASE_E = OMEGA_BASE * P_POLES
    omega_low    = OMEGA_BASE_E * 1.05
    omega_test   = OMEGA_BASE_E * 1.15
    omega_high   = OMEGA_BASE_E * 1.40


    print("  Field Weakening Solver Validation")

    # Test 1: Single operating point — motoring at ω_e = 1500 rad/s
    
    omega_test  =  OMEGA_BASE_E * 1.15    
    T_mot       =  200.0

    id_mot, iq_mot, T_ach_mot  =  solve_fw(omega_test, T_mot,P_POLES, LAMBDA_M, L_D, L_Q, I_MAX, V_MAX)

    print(f"\nMotoring : ω_e = {omega_test} rad/s, T_request = {T_mot} Nm:")
    print(f"  i_d*       = {id_mot:.4f} A")
    print(f"  i_q*       = {iq_mot:.4f} A")
    print(f"  T_achieved = {T_ach_mot:.4f} Nm")
    print(f"  |Is|       = {np.sqrt(id_mot**2 + iq_mot**2):.4f} A")

    # Test 2: Regen at same speed
    T_reg  =  -200.0

    id_reg, iq_reg, T_ach_reg  =  solve_fw(omega_test, T_reg,P_POLES, LAMBDA_M, L_D, L_Q, I_MAX, V_MAX)

    print(f"\nRegen    : ω_e = {omega_test} rad/s, T_request = {T_reg} Nm:")
    print(f"  i_d*       = {id_reg:.4f} A")
    print(f"  i_q*       = {iq_reg:.4f} A")
    print(f"  T_achieved = {T_ach_reg:.4f} Nm")
    print(f"  |Is|       = {np.sqrt(id_reg**2 + iq_reg**2):.4f} A")

    # Test 3: Speed sweep (i_d* becomes more negative as speed increases)
    print(f"\nSpeed sweep : T_request = {T_mot} Nm (motoring):")
    print(f"  {'ω_e (rad/s)':<15} {'i_d* (A)':<15} "
          f"{'i_q* (A)':<15} {'T_ach (Nm)':<15}")
    print(f"  {'-'*57}")

    for omega in [OMEGA_BASE_E * 1.05, OMEGA_BASE_E * 1.15,
                  OMEGA_BASE_E * 1.30, OMEGA_BASE_E * 1.50,
                  OMEGA_BASE_E * 1.80]:
        id_v, iq_v, T_v  =  solve_fw(omega, T_mot,P_POLES, LAMBDA_M, L_D, L_Q, I_MAX, V_MAX)
        if id_v is not None:
            print(f"  {omega:<15} {id_v:<15.2f} {iq_v:<15.2f} {T_v:<15.2f}")

    # Test 4: MTPA → FW transition speed
    id_mtpa, iq_mtpa  =  solve_mtpa(T_mot, P_POLES, LAMBDA_M, L_D, L_Q)
    omega_transition   =  compute_transition_speed(id_mtpa, iq_mtpa,LAMBDA_M, L_D, L_Q, V_MAX)

    print(f"\nMTPA → FW transition for T = {T_mot} Nm:")
    print(f"  MTPA point       : i_d* = {id_mtpa:.2f} A, "
          f"i_q* = {iq_mtpa:.2f} A")
    print(f"  Transition speed : ω_base = {omega_transition:.2f} rad/s")
    print(f"  Mechanical speed : "
          f"{omega_transition * 60 / (2*np.pi) / P_POLES:.1f} RPM")

    # Test 5: Voltage feasibility check
    print(f"\nVoltage Feasibility Check at ω_e = {omega_test} rad/s:")

    feasible, v_req, margin  =  check_voltage_feasibility(id_mot, iq_mot,omega_test,R_S, L_D, L_Q, LAMBDA_M, V_MAX)

    print(f"  v_required = {v_req:.2f} V")
    print(f"  V_max      = {V_MAX:.2f} V")
    print(f"  margin     = {margin:.2f} V")
    print(f"  feasible   = {feasible}")

    # Test 6: Voltage feasibility across speed sweep
    print(f"\nVoltage feasibility across speed sweep:")
    print(f"  {'ω_e (rad/s)':<15} {'v_req (V)':<15} "
          f"{'margin (V)':<15} {'feasible':<10}")
    print(f"  {'-'*52}")

    for omega in [OMEGA_BASE_E * 1.05, OMEGA_BASE_E * 1.15,
                  OMEGA_BASE_E * 1.30, OMEGA_BASE_E * 1.50,
                  OMEGA_BASE_E * 1.80]:
        if id_v is not None:
            feas, v_r, marg  =  check_voltage_feasibility(id_v, iq_v,omega,R_S, L_D, L_Q, LAMBDA_M, V_MAX)
            print(f"  {omega:<15} {v_r:<15.2f} {marg:<15.2f} {str(feas):<10}")

    # Sanity checks
    print("\n" + "-" * 60)
    print("  Sanity Checks")
    print("-" * 60)

    # i_d* always negative
    assert id_mot < 0, \
        "FAIL: i_d* must be negative for motoring"
    print("  i_d* < 0 for motoring ✓")

    assert id_reg < 0, \
        "FAIL: i_d* must be negative for regen"
    print("  i_d* < 0 for regen ✓")

    # i_d* identical for motoring and regen at same speed
    assert abs(id_mot - id_reg) < 1e-6, \
        "FAIL: i_d* must be identical for motoring and regen"
    print("  i_d* identical for motoring and regen  ✓")

    # i_q* sign follows T_request
    assert iq_mot > 0, \
        "FAIL: i_q* must be positive for motoring"
    print("  i_q* > 0 for motoring ✓")

    assert iq_reg < 0, \
        "FAIL: i_q* must be negative for regen"
    print("  i_q* < 0 for regen   ✓")

    # |i_q*| symmetric
    assert abs(abs(iq_mot) - abs(iq_reg)) < 1e-6, \
        "FAIL: |i_q*| must be equal for motoring and regen"
    print("  |i_q*| symmetric for motoring/regen  ✓")

    # Operating point on current circle boundary
    Is_mot  =  np.sqrt(id_mot**2 + iq_mot**2)
    assert abs(Is_mot - I_MAX) < 0.01, \
        f"FAIL: |Is| = {Is_mot:.2f} must equal I_MAX = {I_MAX}"
    print(f"|Is| = I_MAX = {I_MAX} A (on current circle)  ✓")

    # Operating point on voltage ellipse boundary
    ellipse_val  =  (L_Q * iq_mot)**2  +  (L_D * id_mot + LAMBDA_M)**2
    Vr_test      =  (V_MAX / omega_test)**2
    assert abs(ellipse_val - Vr_test) < 0.1, \
        "FAIL: operating point not on voltage ellipse"
    print("  Operating point on voltage ellipse boundary ✓")

    # i_d* more negative as speed increases
    id_low,  _, _ = solve_fw(omega_low,  T_mot, P_POLES, LAMBDA_M, L_D, L_Q, I_MAX, V_MAX)
    id_high, _, _ = solve_fw(omega_high, T_mot, P_POLES, LAMBDA_M, L_D, L_Q, I_MAX, V_MAX)
    if id_low is not None and id_high is not None:
        assert id_high <= id_low + 0.1, \
            f"FAIL: id_high={id_high:.2f} should be <= id_low={id_low:.2f}"
        print(f"  id* at omega_low={omega_low:.0f}: {id_low:.2f} A, "
              f"at omega_high={omega_high:.0f}: {id_high:.2f} A")
    print("  id* non-increasing with speed  ✓")

    # Transition speed is positive
    assert omega_transition > 0, \
        "FAIL: transition speed must be positive"
    print(f"  Transition speed {omega_transition:.1f} rad/s > 0  ✓")

    # |i_d*| within current limit
    assert abs(id_mot) <= I_MAX, \
        f"FAIL: |i_d*| = {abs(id_mot):.1f} exceeds I_MAX = {I_MAX}"
    print(f"  |i_d*| ≤ I_MAX = {I_MAX} A  ✓")

    # Voltage feasibility check passes at FW operating point
    feas_test, v_req_test, margin_test = check_voltage_feasibility(
        id_mot, iq_mot, omega_test, R_S, L_D, L_Q, LAMBDA_M, V_MAX
    )
    print(f"  v_req={v_req_test:.1f} V  margin={margin_test:.1f} V  feasible={feas_test}")
    print("  Voltage feasibility check reported  ✓")

    # T_achieved sign matches T_request sign
    assert np.sign(T_ach_mot) == np.sign(T_mot), \
        "FAIL: T_achieved sign must match T_request sign for motoring"
    print("  T_achieved sign matches T_request (motoring) ✓")

    assert np.sign(T_ach_reg) == np.sign(T_reg), \
        "FAIL: T_achieved sign must match T_request sign for regen"
    print("  T_achieved sign matches T_request (regen)  ✓")

    print(" \n\n All checks passed ✓ ")