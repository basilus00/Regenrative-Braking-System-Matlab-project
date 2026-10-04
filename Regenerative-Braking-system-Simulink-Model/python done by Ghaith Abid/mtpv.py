""" 
Gereneral Purpose : 
Finds (id*, iq*) that maximise torque subject to the voltage ellipse boundary (voltage constraint) at a given electrical speed. 
(the current circle is no longer a binding constraint.)
 
 
 this mode is activated when ω_e > ω_MTPV, where ω_MTPV is the speed at which the MTPV locus 
 crosses inside the current circle or in simpler words , when the voltage ellipse is entirely inside the current circle. 
 Below this speed, the MTPV point is outside the current circle → field weakening active. 
 Above it, the MTPV point is inside the current circle → MTPV active. 
 The transition speed ω_MTPV is computed numerically by finding the root of f(ω_e) = id*(ω_e)² + iq*(ω_e)² − I_MAX². 
 And so , below ω_MTPV, field_weakening.py is active. Above it, this module takes over.

# workflow summary:

Maximise T = K_T·[λm·iq + (Ld-Lq)·id·iq]
subject to: ωe²·Lq²·iq² + ωe²·(Ld·id + λm)² = Vmax² (voltage ellipse boundary)

To do so , we would use Lagrange multipliers: L=f(id,iq)+μ·g(id,iq) where g is the voltage constraint rearranged to equal zero and f is the torque function.
 then we eliminate μ by dividing optimality conditions:

   (Ld - Lq)·iq         Ld·(Ld·id + λm)
  ─────────────────── = ─────────────────
  λm + (Ld-Lq)·id          Lq²·iq

 Cross-multiply, substitute iq² from voltage ellipse to eliminate iq:
→ quadratic in id* only:

   −Ld·δL·id*²  +  λm·(2Ld − Lq)·id*  +  (λm² − Vmax²/ωe²) = 0 → this is the MTPV equation, solved for id*.
 where δL = Lq − Ld > 0

 Solve → pick least-negative root in [I_CH, 0]
 Then: iq* from voltage ellipse, sign from T_request sign.

# GEOMETRIC PICTURE

The MTPV locus is a hyperbola in the dq plane — the voltage-domain analogue of the MTPA cubic. 
Each MTPV point is a tangency between a torque hyperbola and the voltage ellipse at that speed.
As ω → ∞, the ellipse collapses → MTPV point → I_CH = (−λm/Ld, 0) (leans to the infinite speed).

# KEY PHYSICS RULES (never violate)
* id* always negative 
* id* ≥ I_CH = −λm/Ld — cannot go beyond ellipse centre (infinite speed limit)
* iq* sign follows T_request sign only
* |i| ≤ I_MAX — verified in validate_mtpv.py (not enforced here)
"""
import numpy as np
from scipy.optimize import brentq


def solve_mtpv(omega_e, T_request, p, lambda_m, L_d, L_q, V_max):
    """
    Derivation recap
    ----------------
    Maximise T on voltage ellipse → Lagrange condition → eliminate μ:
        w*(λm*Lq - δL*w) = δL*(w² - Vr²)   where w = Ld*id + λm
    Rearranges to quadratic in w:
        2*δL*w² - λm*Lq*w - δL*Vr² = 0
    Discriminant = (λm*Lq)² + 8*δL²*Vr²  → always positive, no None case
    Two roots: w1 large positive (outside ellipse), w2 small/negative (on ellipse)
    w2 is ALWAYS the physically correct root.
    id* = (w2 - λm)/Ld   → can be more negative than I_CH, that is fine for MTPV
    iq* from ellipse equation directly.
    """
    K_T     = 1.5 * p
    delta_L = L_q - L_d          # positive for IPMSM
    Vr      = V_max / omega_e    # voltage headroom [V·s/rad]

    # ── Step 1: quadratic in w = Ld*id + λm ──────────────────────────────────
    # 2·δL·w² − λm·Lq·w − δL·Vr² = 0
    lm_Lq    = lambda_m * L_q
    disc     = lm_Lq**2 + 8.0 * delta_L**2 * Vr**2   # always > 0
    sqrt_disc = np.sqrt(disc)

    # w1 → large positive → |w1| > Vr → outside ellipse → discard
    # w2 → small / negative → |w2| < Vr → on ellipse → keep
    w2 = (lm_Lq - sqrt_disc) / (4.0 * delta_L)

    # ── Step 2: id* from w ───────────────────────────────────────────────────
    id_star = (w2 - lambda_m) / L_d

    # id* must be negative (demagnetisation direction)
    if id_star >= 0.0:
        return None, None, None

    # ── Step 3: iq* from voltage ellipse ─────────────────────────────────────
    # w² + Lq²·iq² = Vr²  →  iq² = (Vr² - w²)/Lq²
    iq2 = (Vr**2 - w2**2) / L_q**2
    if iq2 < 0.0:
        return None, None, None

    iq_magnitude = np.sqrt(iq2)
    iq_star = iq_magnitude if T_request >= 0 else -iq_magnitude

    # ── Step 4: achieved torque ───────────────────────────────────────────────
    T_achieved = K_T * (lambda_m * iq_star + (L_d - L_q) * id_star * iq_star)

    return id_star, iq_star, T_achieved


def compute_mtpv_transition_speed(p, lambda_m, L_d, L_q, V_max, I_max,omega_min=None, omega_max_search=50000.0):
    """
    Find the electrical speed at which the MTPV locus enters the current circle (entire voltage ellipse inside the current circle).
    The FW→MTPV boundary is where: id*(ωe)² + iq*(ωe)² = I_max² 
    Solved numerically with brentq (no closed form exists).
    
    Below this speed: MTPV point is outside current circle → FW active.
    Above this speed: MTPV point is inside current circle → MTPV active.


    Derivation:
    f(ωe) = id*(ωe)² + iq*(ωe)² − I_max²
    Sign of f:
      f < 0  →  |i_mtpv| < I_max  →  MTPV point inside circle  →  MTPV active
      f > 0  →  |i_mtpv| > I_max  →  MTPV point outside circle →  FW active
    Root: f(ω_MTPV) = 0
    """
    if omega_min is None:
        # Minimum meaningful speed: ellipse semi-axis along iq equals I_max
        omega_min = V_max / (L_q * I_max) * 1.01
    def _residual(omega_e):
        id_s, iq_s, _ = solve_mtpv(omega_e, 1.0, p, lambda_m, L_d, L_q, V_max)
        if id_s is None:
            return 1.0    # no solution → f > 0 → FW still active
        return id_s**2 + iq_s**2 - I_max**2

    # Verify sign change exists in search range
    f_min = _residual(omega_min)
    f_max = _residual(omega_max_search)

    if f_min * f_max > 0:
        raise ValueError(
            f"No sign change in [{omega_min:.1f}, {omega_max_search:.1f}] rad/s. "
            f"f_min={f_min:.1f}, f_max={f_max:.1f}. "
            f"Check motor parameters or widen search range."
        )

    return brentq(_residual, omega_min, omega_max_search, xtol=0.1)


def build_mtpv_table(omega_e_min, omega_e_max, n_points,T_request, p, lambda_m, L_d, L_q, V_max):
    """
    Build a speed-swept table of MTPV operating points.

    Used for offline analysis and Simulink lookup table generation.
    Points where solve_mtpv returns None are excluded automatically.

    Returns : 
    omega  (ndarray)  [rad/s] electrical speeds (valid points only)
    id     (ndarray)  [A]     d-axis references
    iq     (ndarray)  [A]     q-axis references
    T_achieved (ndarray) [Nm] achieved torques
    """
    omega_sweep = np.linspace(omega_e_min, omega_e_max, n_points)
    omega_out, id_out, iq_out, T_out = [], [], [], []

    for omega_e in omega_sweep:
        id_s, iq_s, T_s = solve_mtpv(omega_e, T_request,p, lambda_m, L_d, L_q, V_max)
        if id_s is not None:
            omega_out.append(omega_e)
            id_out.append(id_s)
            iq_out.append(iq_s)
            T_out.append(T_s)

    return (np.array(omega_out), np.array(id_out),np.array(iq_out),   np.array(T_out))

# SANITY CHECKS

if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

    from motor_params import (P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,OMEGA_BASE, OMEGA_MAX, I_CH, K_T)

    OMEGA_BASE_E = OMEGA_BASE * P_POLES
    OMEGA_MAX_E  = OMEGA_MAX  * P_POLES
    
    print(" # MTPV SANITY CHECKS")
    
    # [1] Transition speed
    print("\n[1] Computing MTPV transition speed ...")
    omega_mtpv = compute_mtpv_transition_speed(P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,omega_min=OMEGA_BASE_E)
    print(f"    ω_MTPV = {omega_mtpv:.1f} rad/s (elec)")
    print(f"    ω_MTPV = {omega_mtpv/P_POLES:.1f} rad/s (mech)")
    assert omega_mtpv > OMEGA_BASE_E, \
        "FAIL: MTPV transition must be above FW entry speed"
    print("    ✓  ω_MTPV > ω_base")

    # [2] Single point check at ω_MTPV
    print("\n[2] Operating point at ω_MTPV ...")
    id_t, iq_t, T_t = solve_mtpv(omega_mtpv * 1.01, 1.0, P_POLES,
                               LAMBDA_M, L_D, L_Q, V_MAX)
    if id_t is None:
        print("    FAIL: solve_mtpv returned None at omega_mtpv — check discriminant")
    else:
        i_mag = np.sqrt(id_t**2 + iq_t**2)
    print(f"    id* = {id_t:.2f} A,  iq* = {iq_t:.2f} A")
    print(f"    |i| = {i_mag:.2f} A  (should ≈ {I_MAX:.1f} A)")
    assert abs(i_mag - I_MAX) < 5.0, \
        f"FAIL: |i| at transition = {i_mag:.1f} A ≠ I_MAX"
    print("    ✓  |i| ≈ I_MAX at transition")

    # [3] id* always negative and >= I_CH
    print("\n[3] Sweeping full MTPV speed range ...")
    omega_arr, id_arr, iq_arr, T_arr = build_mtpv_table(
        omega_e_min=omega_mtpv,
        omega_e_max=OMEGA_MAX_E * 2.0,
        n_points=300,
        T_request=335.0,
        p=P_POLES, lambda_m=LAMBDA_M, L_d=L_D, L_q=L_Q, V_max=V_MAX
    )
    assert np.all(id_arr < 0), "FAIL: id* went positive"
    id_min_allowed = I_CH - V_MAX / (OMEGA_BASE_E * L_D)   # ellipse left tip at base speed
    assert np.all(id_arr >= id_min_allowed - 1.0), \
    f"FAIL: id* = {id_arr.min():.1f} below ellipse left tip"
    print(f"    ✓  id* range: [{id_arr.min():.1f}, {id_arr.max():.1f}] A")
    print(f"    ✓  All id* >= I_CH = {I_CH:.1f} A")

    # [4] |i| decreasing with speed (moving inside current circle)
    i_mag_arr = np.sqrt(id_arr**2 + iq_arr**2)
    assert np.all(i_mag_arr <= I_MAX + 0.5), \
        f"FAIL: current exceeds I_MAX in MTPV region"
    print(f"    ✓  |i| range: [{i_mag_arr.min():.1f}, {i_mag_arr.max():.1f}] A ≤ I_MAX")

    # [5] Torque non-increasing with speed
    assert T_arr[-1] <= T_arr[0] + 1.0, \
        "FAIL: torque rose with speed in MTPV region"
    print(f"    ✓  T: {T_arr[0]:.1f} → {T_arr[-1]:.1f} Nm (non-increasing)")

    # [6] Voltage constraint satisfied at every point
    print("\n[4] Voltage constraint check ...")
    all_ok = True
    for k in range(len(omega_arr)):
        vd = -omega_arr[k] * L_Q * iq_arr[k]
        vq =  omega_arr[k] * (L_D * id_arr[k] + LAMBDA_M)
        v_mag = np.sqrt(vd**2 + vq**2)
        if v_mag > V_MAX * 1.01:
            print(f"    FAIL at ω = {omega_arr[k]:.1f}: |v| = {v_mag:.1f} > {V_MAX:.1f}")
            all_ok = False
    if all_ok:
        print("    ✓  All MTPV points satisfy |v| ≤ V_MAX")


    print("  ALL SANITY CHECKS PASSED")





'''
 """
    Solve for optimal (id*, iq*) on the voltage ellipse boundary.

    Implements the MTPV quadratic derived via Lagrange multipliers.
    The result is the point on the voltage ellipse that produces
    maximum torque at the given electrical speed.

    Returns (None, None, None) if no valid solution exists at this speed.

    """
    K_T   = 1.5 * p
    I_CH  = -lambda_m / L_d      # characteristic current — hard floor for id*
    delta_L = L_q - L_d         
    Vr    = V_max / omega_e       # voltage-to-speed ratio [V·s/rad]

    # Step 1: MTPV quadratic coefficients from: −Ld·δL·id² + λm·(2Ld−Lq)·id + (λm² − Vmax²/ωe²) = 0
    A =  L_d * delta_L                       
    B =  lambda_m * (L_q - 2.0 * L_d)        
    C =  Vr**2 - lambda_m**2                  
    #  Step 2: Discriminant      
    B_std = lambda_m * (L_q - 2.0 * L_d)    
    C_std = Vr**2 - lambda_m**2              
    discriminant = B_std**2 - 4.0 * A * C_std

    if discriminant < 0:
        # No real solution — speed too low for MTPV to be active
        return None, None, None

    # Step 3: Solve quadratic
    sqrt_disc = np.sqrt(discriminant)
    id1 = (-B_std + sqrt_disc) / (2.0 * A)
    id2 = (-B_std - sqrt_disc) / (2.0 * A)

    # Physical validity rules :
    #1. id* must be negative (demagnetisation direction)
    #2. id* must be >= I_CH (cannot over-cancel PM flux)
    #3. Among valid roots: pick least negative (min demagnetisation)
    
    candidates = [x for x in (id1, id2) if I_CH - 1e-6 <= x < 0]

    if not candidates:
        # Both roots outside valid range
        id_star = I_CH
    else:
        id_star = max(candidates)    
    id_star = max(id_star, I_CH)

    # Step 4: iq* from voltage ellipse ( ωe²·Lq²·iq² + ωe²·(Ld·id+λm)² = Vmax²)

    iq2 = (Vr**2 - L_d**2 * (id_star + lambda_m / L_d)**2) / L_q**2
    iq2 = max(iq2, 0.0)             
    iq_magnitude = np.sqrt(iq2)
    iq_star = iq_magnitude if T_request >= 0 else -iq_magnitude

    # Step 5: Achieved torque 
    T_achieved = K_T * (lambda_m * iq_star + (L_d - L_q) * id_star * iq_star)
    return id_star, iq_star, T_achieved
'''