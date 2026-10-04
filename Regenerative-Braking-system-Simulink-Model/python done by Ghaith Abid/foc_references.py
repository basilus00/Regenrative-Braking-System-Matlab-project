'''
# PURPOSE
This is the brain of the reference generation layer. It takes T_request and ω_e as inputs and decides which algorithm to call based
on the current operating region. It is the ONLY file that imports from all three solvers.

# OPERATING REGIONS
* Region 1 — MTPA:  ω_e ≤ ω_base(T)
Voltage ellipse is large —> current circle is the only active constraint.
Target : Minimise |i| for the requested torque. id* from MTPA cubic, iq* from current circle.

* Region 2 — FW:    ω_base(T) < ω_e ≤ ω_MTPV
Both constraints (voltage and current) active simultaneously.
id* driven by speed, not torque and torque is a consequence.
Operating point slides along current circle arc toward I_CH.

* Region 3 — MTPV:  ω_e > ω_MTPV
Voltage ellipse fits inside current circle —> current no longer binding constraint(only voltage matters here).
Maximise torque on voltage ellipse boundary. Operating point moves inward along MTPV hyperbola toward I_CH.

# Transition logic
*ω_base(T):  computed per torque request from field_weakening.py. represents the transition speed between MTPA and FW. 
Below this speed, voltage constraint is inactive and we can follow MTPA locus. Above this speed, voltage constraint is active and we must switch to FW.

* ω_MTPV: computed ONCE at startup in mtpv.py. Represents the transition speed between FW and MTPV. 
Below this speed, voltage constraint is smaller than current constraint at I_CH and we must use FW.
Above this speed, voltage constraint is larger than current constraint at I_CH and we can switch to MTPV.

# Key rules and assumptions
- id* is ALWAYS negative regardless of region or torque sign
- iq* sign follows T_request sign only (+motoring, -regen)
- If FW returns None (speed too high) → fall through to MTPV
- If MTPV returns None (speed too low) → something is wrong → return safe zero
- T_achieved may be less than |T_request| in FW and MTPV ( must handle this assertively )

'''
import numpy as np
from scipy.optimize import brentq

import sys
import os

if __name__ == "__main__":
    _repo_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..')
    )
    sys.path.insert(0, _repo_root)

from python.motor_params import (P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,OMEGA_BASE, OMEGA_MAX, I_CH, K_T, T_max)
from python.mtpa import solve_mtpa
from python.field_weakening import solve_fw, compute_transition_speed
from python.mtpv import solve_mtpv, compute_mtpv_transition_speed

# Initialisation 

def initialise(p, lambda_m, L_d, L_q, V_MAX, I_max,omega_min=None, omega_max_search=100000.0):
    """
    Compute and return ω_MTPV (the FW→MTPV transition speed).
    This is torque-independent and depends only on motor parameters.
    Call once before entering the control loop and pass the result to
    every call of compute_references().
    """
    omega_mtpv = compute_mtpv_transition_speed(p, lambda_m, L_d, L_q, V_MAX, I_max,omega_min=omega_min,
omega_max_search=omega_max_search,)
    return omega_mtpv

# swwitching logic

def compute_references(T_request, omega_e, p, lambda_m, L_d, L_q, V_MAX, I_max, omega_mtpv):
    """
    Compute (id*, iq*) reference currents for the given operating point.
    """

    # Use ELECTRICAL speed consistently
    omega_el = abs(float(omega_e))

    # Temporary unit sanity warning
    if omega_el < 1000 and abs(T_request) > 1.0:
        print(
            f"WARNING: low omega_el={omega_el:.1f} rad/s. "
            f"If this is high-speed test, you may be passing mechanical speed by mistake."
        )

    # Step 1
    T_abs = abs(T_request)
    if T_abs < 1e-3:
        return 0.0, 0.0, 0.0, 0  # idle region

    id_mtpa, iq_mtpa = solve_mtpa(T_abs, p, lambda_m, L_d, L_q)

    print(
        f"  DEBUG calling transition: id_mtpa={id_mtpa:.2f}, iq_mtpa={iq_mtpa:.2f}, "
        f"V_max={V_MAX:.2f}, L_d={L_d:.8f}, L_q={L_q:.8f}, lambda_m={lambda_m:.6f}"
    )

    # Step 2
    omega_base_T = compute_transition_speed(id_mtpa, iq_mtpa, lambda_m, L_d, L_q, V_MAX)

    print(f"  DEBUG omega_base_T returned = {omega_base_T}")
    print(
        f"  DEBUG: omega_e={omega_el:.1f}, omega_base_T={omega_base_T:.1f}, "
        f"omega_mtpv={omega_mtpv:.1f}"
    )
    print(f"  DEBUG: id_mtpa={id_mtpa:.2f}, iq_mtpa={iq_mtpa:.2f}")

    # Step 3+4
    iq_sign = 1.0 if T_request >= 0 else -1.0

    if omega_el <= omega_base_T:
        # Region 1 (MTPA)
        id_star = id_mtpa
        iq_star = iq_sign * abs(iq_mtpa)
        K_T = 1.5 * p
        T_achieved = K_T * (lambda_m * iq_star + (L_d - L_q) * id_star * iq_star)
        return id_star, iq_star, T_achieved, 1

    elif omega_el <= omega_mtpv:   # <-- FIXED (was omega_e)
        # Region 2 (FW)
        id_star, iq_star, T_achieved = solve_fw(
            omega_el, T_request, p, lambda_m, L_d, L_q, I_max, V_MAX
        )
        if id_star is None:
            id_star, iq_star, T_achieved = solve_mtpv(
                omega_el, T_request, p, lambda_m, L_d, L_q, V_MAX
            )
            if id_star is None:
                return 0.0, 0.0, 0.0, 3
            return id_star, iq_star, T_achieved, 3
        return id_star, iq_star, T_achieved, 2

    else:
        # Region 3 (MTPV)
        id_star, iq_star, T_achieved = solve_mtpv(
            omega_el, T_request, p, lambda_m, L_d, L_q, V_MAX
        )
        if id_star is None:
            return 0.0, 0.0, 0.0, 3
        return id_star, iq_star, T_achieved, 3
# OFFLINE TABLE BUILDER — for Simulink lookup and validation

def build_reference_table(T_request, omega_e_min, omega_e_max, n_points,p, lambda_m, L_d, L_q, V_MAX, I_max, omega_mtpv):
    """
    Sweep electrical speed and compute references at each point.
    """
    omega_sweep = np.linspace(omega_e_min, omega_e_max, n_points)
    id_arr  = np.zeros(n_points)
    iq_arr  = np.zeros(n_points)
    T_arr   = np.zeros(n_points)
    reg_arr = np.zeros(n_points, dtype=int)

    for k, omega_e in enumerate(omega_sweep):
        id_arr[k], iq_arr[k], T_arr[k], reg_arr[k] = compute_references(T_request, omega_e,p, lambda_m, L_d, L_q, V_MAX, I_max, omega_mtpv)

    return {
        'omega'  : omega_sweep,
        'id'     : id_arr,
        'iq'     : iq_arr,
        'T'      : T_arr,
        'region' : reg_arr,
    }



# SANITY CHECKS

if __name__ == "__main__":

    OMEGA_BASE_E = OMEGA_BASE * P_POLES
    OMEGA_MAX_E  = OMEGA_MAX  * P_POLES


    print("  FOC REFERENCES — SANITY CHECKS")


    # [1] Initialise — compute ω_MTPV
    print("\n[1] Initialising — computing ω_MTPV ...")
    omega_mtpv = initialise(P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,
                             omega_min=OMEGA_BASE_E)
    print(f"    ω_MTPV = {omega_mtpv:.1f} rad/s (elec)")
    print(f"    ω_MTPV = {omega_mtpv/P_POLES:.1f} rad/s (mech)")

    # [2] Single point in each region
    print("\n[2] Single point checks ...")
    test_cases = [
        ("MTPA  ", OMEGA_BASE_E * 0.5,  T_max, 1),
        ("FW    ", OMEGA_BASE_E * 1.5,  T_max, 2),
        ("MTPV  ", omega_mtpv  * 1.5,  T_max, 3),
        ("Regen ", OMEGA_BASE_E * 1.5, -T_max, 2),
    ]

    for name, omega_e, T_req, expected_region in test_cases:
        id_s, iq_s, T_s, reg = compute_references(
            T_req, omega_e, P_POLES, LAMBDA_M, L_D, L_Q,
            V_MAX, I_MAX, omega_mtpv
        )
        assert id_s < 0, f"FAIL [{name}]: id* = {id_s:.1f} not negative"
        assert reg == expected_region, \
            f"FAIL [{name}]: region = {reg}, expected {expected_region}"
        if T_req < 0:
            assert iq_s < 0, f"FAIL [{name}]: iq* should be negative for regen"
        print(f"    [{name}] ω={omega_e/P_POLES:.0f} rad/s mech | "
              f"id*={id_s:.1f} A  iq*={iq_s:.1f} A  "
              f"T={T_s:.1f} Nm  region={reg}  ✓")

    # [3] Full speed sweep — check id* always negative
    print("\n[3] Full speed sweep ...")
    table = build_reference_table(
        T_max,
        OMEGA_BASE_E * 0.1,
        omega_mtpv * 2.5,
        500,
        P_POLES,LAMBDA_M, L_D, L_Q,
        V_MAX, I_MAX, omega_mtpv)
    assert np.all(table['id'] < 0), "FAIL: id* went positive somewhere"
    print(f"    ✓  id* always negative: [{table['id'].min():.1f}, {table['id'].max():.1f}] A")

    # [4] Region transitions — check they happen in correct order
    regions = table['region']
    mtpa_end   = np.where(regions == 1)[0][-1]
    fw_start   = np.where(regions == 2)[0][0]
    fw_end     = np.where(regions == 2)[0][-1]
    mtpv_start = np.where(regions == 3)[0][0]
    assert fw_start > mtpa_end,   "FAIL: FW starts before MTPA ends"
    assert mtpv_start > fw_end,   "FAIL: MTPV starts before FW ends"
    print(f"    ✓  Region order: MTPA → FW → MTPV (no overlap)")

   # [5] Voltage constraint satisfied everywhere
    print("\n[5] Voltage constraint check ...")
    all_ok = True
    eps_v = 1e-6

    for k in range(len(table['omega'])):
        vd = -table['omega'][k] * L_Q * table['iq'][k]
        vq =  table['omega'][k] * (L_D * table['id'][k] + LAMBDA_M)
        v_mag = np.sqrt(vd**2 + vq**2)
        over = v_mag - V_MAX

        if over > eps_v:
            print(
                f"    FAIL at ω = {table['omega'][k]/P_POLES:.1f} rad/s mech | "
                f"region={table['region'][k]} | v={v_mag:.3f} V | over={over:.3f} V"
            )   
            all_ok = False

    if all_ok:
        print("    ✓ Voltage constraint satisfied at all points")
        print("\nALL SANITY CHECKS PASSED")
    else:
        raise AssertionError("Voltage constraint failed at one or more sweep points")