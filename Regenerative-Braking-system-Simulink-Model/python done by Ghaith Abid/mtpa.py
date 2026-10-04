"""
Problem:
    Minimize: f(i_d, i_q) =i_d²+i_q²
    Subject to: g(i_d, i_q) = K·i_q·[λ_m + (L_d-L_q)·i_d] - T_request = 0
    where K = (3/2)·p

Derivation:
    1. Form Lagrangian:
          L(i_d, i_q, μ) = (i_d² + i_q²) + μ·g(i_d, i_q)

    2. First-order conditions (∂L/∂i_d = 0, ∂L/∂i_q = 0, ∂L/∂μ = 0):
          (I)  : 2·i_d + μ·K·i_q·(L_d-L_q)          = 0
          (II) : 2·i_q + μ·K·[λ_m + (L_d-L_q)·i_d]  = 0
          (III): K·i_q·[λ_m + (L_d-L_q)·i_d]         = T_request

    3. Eliminate μ by dividing (I) by (II) → MTPA condition:
          (L_d-L_q)·i_d² + λ_m·i_d - (L_d-L_q)·i_q² = 0     

    4. Two-equation system:
          MTPA eq:(L_d-L_q)·i_d² + λ_m·i_d - (L_d-L_q)·i_q² = 0
          Torque eq: K·i_q·[λ_m + (L_d-L_q)·i_d]= T_request

    5. Let A=λ_m + (L_d-L_q)·i_d  :
          MTPA eq:   i_d·A = (L_d-L_q)·i_q²
          Torque eq: i_q·A = T_request/K

    6. Divide MTPA by Torque :
          i_d = [K·(L_d-L_q) / T_request] · i_q³ → MTPA curve

    7. Substitute into Torque equation → deg-4 polynomial in i_q:
          a·i_q⁴ + b·i_q + c = 0
          where:
             a = K·(L_d-L_q)²
             b = λ_m·T_request
             c = -T_request²/K

    8. Solve with Newton-Raphson → i_q*
        Recover i_d* directly from MTPA curve

Sign behavior (for later classification):
    Motoring (T_request > 0): i_q* > 0, i_d* < 0
    Regen    (T_request < 0): i_q* < 0, i_d* < 0  (always negative for IPMSM)
"""

import numpy as np


def solve_mtpa(
    T_request : float,
    p         : int,
    lambda_m  : float,
    L_d       : float,
    L_q       : float,
    tol       : float = 1e-6,
    max_iter  : int   = 50
) -> tuple:
    """
    Parameters:
    T_request (float)  Torque demand (Nm). Positive=motoring, Negative=regen.
    p         (int)    Number of pole pairs.
    lambda_m  (float)  PM flux linkage (Wb).
    L_d       (float)  d-axis inductance (H).
    L_q       (float)  q-axis inductance (H).
    tol       (float)  Convergence tolerance for Newton-Raphson algorithm (A).
    max_iter  (int)    Maximum Newton-Raphson iterations.

    function returns
    id_star : float  Optimal d-axis current (A). Always negative for IPMSM.
    iq_star : float  Optimal q-axis current (A). Sign follows T_request.
    """

    if T_request == 0.0: # Edge case: zero torque → no current needed
        return 0.0, 0.0
    
    # Step 1 — Constants
    K  = 1.5 * p
    dL = L_d - L_q
    
    # Step 2 — Polynomial coefficients
    # From substituting MTPA curve into torque equation (step 7 in documentation): a·i_q⁴ + b·i_q + c = 0
    a =  K * dL**2
    b =  lambda_m * T_request
    c = -(T_request**2) / K
    
    # Step 3 — Newton-Raphson to solve f(i_q) = 0 where f(i_q) = a·i_q⁴ + b·i_q + c and f'(i_q) = 4·a·i_q³ + b
    # Initial guess: pure PM torque assumption (i_d = 0) and so T_request ≈ K·λ_m·i_q  →  i_q⁰ = T_request / (K·λ_m)

    iq = T_request / (K * lambda_m)

    for i in range(max_iter):
        f=a*iq**4+b*iq+c  
        df=4*a*iq**3+b          
        corr=f/df       # Newton-Raphson correction
        iq-= corr
        if abs(corr)<tol:  # convergence check
            break

    # Step 4 — Recover i_d* directly from MTPA curve (step 6 in documentation) where i_d*=[K·(L_d-L_q)/T_request]·i_q*³
    id = (K * dL / T_request) * iq**3
    return id, iq

def build_mtpa_table(  
    T_max    : float,
    n_points : int,
    p        : int,
    lambda_m : float,
    L_d      : float,
    L_q      : float) -> tuple:
    """
    Pre-compute MTPA lookup table for T in [-T_max, +T_max]. This table is used to generate offline table for Simulink lookup block since 
    the FOC loop in Simulink would actually runs 10000 to 20000 times per second making running Newton-Raphson algorithm impractical 
    inside a real-time loop. 
    And so we solves the optimization once for all,offline and for every torque value from -T_max to +T_max,then, stores all the results in a table
    that will be exported and loaded in Simulink as 1D lookup table block .

    Parameters:
    T_max    (float)  Maximum torque magnitude (Nm).
    n_points (int)    Number of table points (odd number recommended so that T=0 is included too).
    p        (int)    Number of pole pairs.
    lambda_m (float)  PM flux linkage (Wb).
    L_d      (float)  d-axis inductance (H).
    L_q      (float)  q-axis inductance (H).

    this function returns:
    torques  (array)  Torque values       (Nm),shape(n_points,)
    id_table (array)  Optimal i_d* values (A),shape(n_points,)
    iq_table (array)  Optimal i_q* values (A),shape(n_points,)
    """

    # Symmetric range: negative (regen) → zero → positive (motoring)
    torques  = np.linspace(-T_max, T_max, n_points)
    id_table = np.zeros(n_points)
    iq_table = np.zeros(n_points)

    for i, T in enumerate(torques):
        id_table[i], iq_table[i] = solve_mtpa(T, p, lambda_m, L_d, L_q)
    return torques, id_table, iq_table

# Validation — runs only when executing this file directly

if __name__ == "__main__":

    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    from motor_params import (P_POLES, LAMBDA_M, L_D, L_Q, T_max, I_MAX, K_T)

    print("=" * 55)
    print("MTPA Solver Validation")
    print("=" * 55)

    # Test 1: Motoring (T > 0)
    T_mot = 200.0
    id_mot, iq_mot = solve_mtpa(T_mot, P_POLES, LAMBDA_M, L_D, L_Q)
    Is_mot = np.sqrt(id_mot**2 + iq_mot**2)
    print(f"\nMotoring:T_request = +{T_mot} Nm:")
    print(f"i_d*={id_mot:.4f}A")
    print(f"i_q*={iq_mot:.4f}A")
    print(f"|Is|={Is_mot:.4f}A")

    # Test 2: Regen (T < 0)
    T_reg = -200.0
    id_reg, iq_reg = solve_mtpa(T_reg, P_POLES, LAMBDA_M, L_D, L_Q)
    Is_reg = np.sqrt(id_reg**2 + iq_reg**2)
    print(f"\nRegen: T_request={T_reg} Nm:")
    print(f"i_d*={id_reg:.4f}A")
    print(f"i_q*={iq_reg:.4f}A")
    print(f"|Is|={Is_reg:.4f}A")
    
    # Test 3: Zero torque
    id_zero, iq_zero = solve_mtpa(0.0, P_POLES, LAMBDA_M, L_D, L_Q)
    print(f"\nZero:T_request = 0 Nm:")
    print(f"i_d*={id_zero:.4f} A")
    print(f"i_q*={iq_zero:.4f} A")

    # Sanity checks
    print("\n" + "-" * 55)
    print("  Sanity Checks")
    print("-" * 55)

    # i_d* always negative for IPMSM
    assert id_mot < 0, \
        "FAIL: i_d* must be negative for motoring"
    print("  i_d* < 0 for motoring                    ✓")

    assert id_reg < 0, \
        "FAIL: i_d* must be negative for regen"
    print("  i_d* < 0 for regen                       ✓")

    # i_q* sign follows torque sign
    assert iq_mot > 0, \
        "FAIL: i_q* must be positive for motoring"
    print("  i_q* > 0 for motoring                    ✓")

    assert iq_reg < 0, \
        "FAIL: i_q* must be negative for regen"
    print("  i_q* < 0 for regen                       ✓")

    # Motoring/regen symmetry
    assert abs(id_mot - id_reg) < 1e-6, \
        "FAIL: i_d* must be equal for +T and -T"
    print("  i_d* symmetric for motoring/regen        ✓")

    assert abs(abs(iq_mot) - abs(iq_reg)) < 1e-6, \
        "FAIL: |i_q*| must be equal for +T and -T"
    print("  |i_q*| symmetric for motoring/regen      ✓")

    # Current within inverter limit
    assert Is_mot <= I_MAX, \
        f"FAIL: |Is|={Is_mot:.1f} A exceeds I_MAX={I_MAX} A"
    print(f"  |Is| = {Is_mot:.1f} A ≤ I_MAX = {I_MAX} A        ✓")

    # Torque constraint satisfied — plug (i_d*, i_q*) back into torque eq
    # T = K·[λ_m·i_q + (L_d-L_q)·i_d·i_q]
    T_check = K_T * (LAMBDA_M * iq_mot + (L_D - L_Q) * id_mot * iq_mot)
    assert abs(T_check - T_mot) < 0.01, \
        f"FAIL: torque constraint violated. Got {T_check:.3f} Nm"
    print(f"  Torque constraint: {T_check:.3f} ≈ {T_mot} Nm      ✓")

    # MTPA curve relationship: i_d* = C · i_q*³
    C       = K_T * (L_D - L_Q) / T_mot
    id_from_curve = C * iq_mot**3
    assert abs(id_from_curve - id_mot) < 1e-6, \
        "FAIL: i_d* does not satisfy MTPA curve i_d = C·i_q³"
    print("  MTPA curve satisfied: i_d* = C·i_q*³    ✓")

    print("\n" + "=" * 55)
    print("  All checks passed ✓")
    print("=" * 55)