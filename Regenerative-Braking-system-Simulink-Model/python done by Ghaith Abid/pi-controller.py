"""
Purpose :
Convert current references (id_ref, iq_ref) and measured currents (id_meas, iq_meas) into voltage commands (vd_cmd, vq_cmd) with:
- PI regulation on each axis
- cross-coupling + back-EMF feedforward
- vector saturation to inverter limit
- anti-windup (back-calculation)

Notice: 
The implemented control architecture is intentionally split into two layers:
# The 1st layer : (MTPA/FW/MTPV reference generation) uses steady-state equations to compute ((i_d^, i_q^)) efficiently 
and robustly across speed regions. This keeps the optimization tractable and ensures deterministic region transitions.
Yet, this assumptions is but a simplification. It does not fully capture fast transients, discretization effects
, inverter saturation, parameter mismatch, and disturbance dynamics. Therefore, by itself, it is not sufficient for precise real-time current tracking.

The 2nd layer : (inner PI current controller) handles exactly those dynamic effects. 
It compares measured currents to references each sampling cycle and corrects residual errors through feedback. 
Feedforward decoupling terms improve axis decoupling, while saturation + anti-windup preserve stability when voltage limits are reached.

So in short , the design choice is deliberate:
steady-state model for optimal references, then dynamic feedback for real-world correction.

Conventionnaly, the voltage equations in the rotating d-q frame are: 
vd = Rs*id + Ld*did/dt - omega_e*Lq*iq
vq = Rs*iq + Lq*diq/dt + omega_e*(Ld*id + lambda_m)

So feedforward terms used here:
vd_ff = - omega_e * Lq * iq_meas
vq_ff = + omega_e * (Ld * id_meas + lambda_m)

"""
import numpy as np


def compute_pi_gains(R_s, L_d, L_q, omega_c):
    """
    Simple initial tuning based on bandwidth omega_c:
      Kp_d = Ld * omega_c
      Ki_d = Rs * omega_c
      Kp_q = Lq * omega_c
      Ki_q = Rs * omega_c
    """
    Kp_d = L_d * omega_c
    Ki_d = R_s * omega_c
    Kp_q = L_q * omega_c
    Ki_q = R_s * omega_c
    return Kp_d, Ki_d, Kp_q, Ki_q


def saturate_voltage_vector(v_d, v_q, V_max):
    """
    Vector saturation: keep voltage direction, limit magnitude to V_max.
    """
    v_mag = np.sqrt(v_d**2 + v_q**2)
    if v_mag <= V_max or v_mag < 1e-12:
        return v_d, v_q
    scale = V_max / v_mag
    return v_d * scale, v_q * scale


def pi_current_step(
    id_ref, iq_ref, id_meas, iq_meas, omega_e,
    int_d, int_q,
    R_s, L_d, L_q, lambda_m,
    V_max, Ts, omega_c, kaw=1.0,
    gains=None,
    omega_ff_min=0.0
):
    """
    One control step of dq current PI with feedforward + anti-windup.

    Inputs:
      gains: optional tuple (Kp_d, Ki_d, Kp_q, Ki_q) to avoid recomputing each call
      omega_ff_min: feedforward disabled when |omega_e| < omega_ff_min
    """

    # 1) PI gains
    if gains is None:
        Kp_d, Ki_d, Kp_q, Ki_q = compute_pi_gains(R_s, L_d, L_q, omega_c)
    else:
        Kp_d, Ki_d, Kp_q, Ki_q = gains

    # 2) Current errors
    e_d = id_ref - id_meas
    e_q = iq_ref - iq_meas

    # 3) Feedforward (optional low-speed gating)
    if abs(omega_e) >= omega_ff_min:
        v_d_ff = -omega_e * L_q * iq_meas
        v_q_ff =  omega_e * (L_d * id_meas + lambda_m)
    else:
        v_d_ff = 0.0
        v_q_ff = 0.0

    # 4) Unsaturated PI outputs
    v_d_unsat = Kp_d * e_d + int_d + v_d_ff
    v_q_unsat = Kp_q * e_q + int_q + v_q_ff

    # 5) Saturation
    v_d_cmd, v_q_cmd = saturate_voltage_vector(v_d_unsat, v_q_unsat, V_max)

    # 6) Anti-windup (back-calculation)
    aw_d = kaw * (v_d_cmd - v_d_unsat)
    aw_q = kaw * (v_q_cmd - v_q_unsat)

    # 7) Integrator update
    int_d_new = int_d + Ts * (Ki_d * e_d + aw_d)
    int_q_new = int_q + Ts * (Ki_q * e_q + aw_q)

    debug = {
        "e_d": e_d, "e_q": e_q,
        "v_d_ff": v_d_ff, "v_q_ff": v_q_ff,
        "v_d_unsat": v_d_unsat, "v_q_unsat": v_q_unsat,
        "v_mag_cmd": np.sqrt(v_d_cmd**2 + v_q_cmd**2),
        "aw_d": aw_d, "aw_q": aw_q
    }

    return v_d_cmd, v_q_cmd, int_d_new, int_q_new, debug


# VALIDATION (same style as your modules)
if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from motor_params import R_S, L_D, L_Q, LAMBDA_M, V_MAX

    print(" PI CONTROLLER — SANITY CHECKS")

    Ts = 100e-6
    omega_c = 2 * np.pi * 400
    kaw = 1.0
    omega_ff_min = 20.0   # rad/s, optional gating near standstill

    # Precompute gains once
    gains = compute_pi_gains(R_S, L_D, L_Q, omega_c)

    int_d = 0.0
    int_q = 0.0

    # Fixed test point (can be replaced by plant model later)
    id_ref = -220.0
    iq_ref = 330.0
    id_meas = -200.0
    iq_meas = 300.0
    omega_e = 3000.0

    n_steps = 100
    vmag_hist = []
    intd_hist = []
    intq_hist = []

    for _ in range(n_steps):
        v_d, v_q, int_d, int_q, dbg = pi_current_step(
            id_ref, iq_ref, id_meas, iq_meas, omega_e,
            int_d, int_q,
            R_S, L_D, L_Q, LAMBDA_M,
            V_MAX, Ts, omega_c, kaw,
            gains=gains,
            omega_ff_min=omega_ff_min
        )
        vmag_hist.append(dbg["v_mag_cmd"])
        intd_hist.append(int_d)
        intq_hist.append(int_q)

    print(f"  Final v_d* = {v_d:.2f} V")
    print(f"  Final v_q* = {v_q:.2f} V")
    print(f"  Final |v|  = {vmag_hist[-1]:.2f} V  (<= V_MAX={V_MAX:.2f} V)")
    print(f"  Final int_d={int_d:.4f}, int_q={int_q:.4f}")

    assert np.max(vmag_hist) <= V_MAX + 1e-6, "FAIL: voltage saturation failed"
    assert np.isfinite(int_d) and np.isfinite(int_q), "FAIL: integrator became non-finite"

    print("  ✓ Voltage vector saturation works across loop")
    print("  ✓ Integrator states evolve and remain finite")
    print("  ✓ Multi-step sanity check passed")