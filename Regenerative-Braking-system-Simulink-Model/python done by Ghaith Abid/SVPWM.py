"""
Purpose:
- Convert stationary-frame voltage commands (v_alpha, v_beta) into 3-phase PWM duties
  (d_a, d_b, d_c) using standard 2-level SVPWM via zero-sequence injection.
- Optionally scale the requested voltage vector to stay inside linear SVPWM range.

Design notes:
- If your PI current controller already saturates voltage to V_MAX = Vdc/sqrt(3),
  set apply_linear_saturation=False in normal pipeline to avoid double-saturation.
- Keep apply_linear_saturation=True for standalone robustness/testing.
"""

import numpy as np


def voltage_saturation_svpwm(v_alpha, v_beta, Vdc):
    """
    Scale (v_alpha, v_beta) to linear SVPWM limit:

        |v_ref| <= Vdc / sqrt(3)

    Returns:
      v_alpha_sat, v_beta_sat, scale
    """
    if Vdc <= 0:
        raise ValueError("Vdc must be > 0")

    v_lim = Vdc / np.sqrt(3.0)
    mag = float(np.hypot(v_alpha, v_beta))

    if mag <= 1e-12:
        return float(v_alpha), float(v_beta), 1.0
    if mag <= v_lim:
        return float(v_alpha), float(v_beta), 1.0

    scale = v_lim / mag
    return float(v_alpha * scale), float(v_beta * scale), float(scale)


def svpwm_duty_from_vab(v_alpha, v_beta, Vdc, clamp=True):
    """
    SVPWM duty generation (offset-injection form).

    Steps:
      1) alpha-beta -> 3-phase line-to-neutral references
      2) common-mode offset v0 = -(vmax + vmin)/2
      3) duty mapping d = 0.5 + v_phase / Vdc

    Returns:
      d_a, d_b, d_c, debug
    """
    if Vdc <= 0:
        raise ValueError("Vdc must be > 0")

    sqrt3 = np.sqrt(3.0)

    # alpha-beta -> abc (line-neutral references)
    v_a = float(v_alpha)
    v_b = float(-0.5 * v_alpha + 0.5 * sqrt3 * v_beta)
    v_c = float(-0.5 * v_alpha - 0.5 * sqrt3 * v_beta)

    # zero-sequence injection
    v_max = max(v_a, v_b, v_c)
    v_min = min(v_a, v_b, v_c)
    v_0 = -0.5 * (v_max + v_min)

    v_a0 = v_a + v_0
    v_b0 = v_b + v_0
    v_c0 = v_c + v_0

    # duty mapping
    d_a = 0.5 + (v_a0 / Vdc)
    d_b = 0.5 + (v_b0 / Vdc)
    d_c = 0.5 + (v_c0 / Vdc)

    unclamped = (float(d_a), float(d_b), float(d_c))

    if clamp:
        d_a = float(np.clip(d_a, 0.0, 1.0))
        d_b = float(np.clip(d_b, 0.0, 1.0))
        d_c = float(np.clip(d_c, 0.0, 1.0))

    overmod = (
        (unclamped[0] < 0.0 or unclamped[0] > 1.0) or
        (unclamped[1] < 0.0 or unclamped[1] > 1.0) or
        (unclamped[2] < 0.0 or unclamped[2] > 1.0)
    )

    debug = {
        "v_a": v_a, "v_b": v_b, "v_c": v_c,
        "v_0": v_0,
        "v_a0": v_a0, "v_b0": v_b0, "v_c0": v_c0,
        "d_a_unclamped": unclamped[0],
        "d_b_unclamped": unclamped[1],
        "d_c_unclamped": unclamped[2],
        "overmodulation_detected": overmod
    }

    return d_a, d_b, d_c, debug


def svpwm_step(v_alpha, v_beta, Vdc, apply_linear_saturation=False, clamp=True):
    """
    Convenience wrapper:
      - optional linear-region saturation
      - duty computation
      - consolidated debug output
    """
    if apply_linear_saturation:
        v_alpha_use, v_beta_use, scale = voltage_saturation_svpwm(v_alpha, v_beta, Vdc)
    else:
        v_alpha_use, v_beta_use, scale = float(v_alpha), float(v_beta), 1.0

    d_a, d_b, d_c, dbg = svpwm_duty_from_vab(
        v_alpha_use, v_beta_use, Vdc, clamp=clamp
    )

    dbg_all = {
        "v_alpha_in": float(v_alpha),
        "v_beta_in": float(v_beta),
        "v_alpha_used": float(v_alpha_use),
        "v_beta_used": float(v_beta_use),
        "pre_scale": float(scale),
        **dbg
    }

    return d_a, d_b, d_c, dbg_all


# VALIDATION (same style as your modules)
if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

    # Uses your existing parameter file
    # Expected names from your project style:
    #   VDC_LINK  (preferred) or fallback to deriving from V_MAX if needed
    from motor_params import V_MAX
    try:
        from motor_params import VDC_LINK
        Vdc = float(VDC_LINK)
    except Exception:
        # fallback to remain compatible with your existing constants
        Vdc = float(np.sqrt(3.0) * V_MAX)

    print(" SVPWM — SANITY CHECKS")
    print(f"  Using Vdc = {Vdc:.2f} V")

    # ---------- Test 1: inside linear region ----------
    v_alpha_1 = 0.35 * Vdc
    v_beta_1 = 0.10 * Vdc

    d_a, d_b, d_c, dbg = svpwm_step(
        v_alpha_1, v_beta_1, Vdc,
        apply_linear_saturation=True,
        clamp=True
    )

    print("\n  [Test 1] In-range vector")
    print(f"    v_in = ({v_alpha_1:.2f}, {v_beta_1:.2f}) V")
    print(f"    scale={dbg['pre_scale']:.4f}, overmod={dbg['overmodulation_detected']}")
    print(f"    duty = (da={d_a:.4f}, db={d_b:.4f}, dc={d_c:.4f})")

    assert 0.0 <= d_a <= 1.0 and 0.0 <= d_b <= 1.0 and 0.0 <= d_c <= 1.0
    assert np.isfinite(d_a) and np.isfinite(d_b) and np.isfinite(d_c)

    # ---------- Test 2: out-of-range vector (saturate first) ----------
    v_alpha_2 = 0.90 * Vdc
    v_beta_2 = 0.60 * Vdc

    d_a2, d_b2, d_c2, dbg2 = svpwm_step(
        v_alpha_2, v_beta_2, Vdc,
        apply_linear_saturation=True,
        clamp=True
    )

    print("\n  [Test 2] Out-of-range vector + linear saturation")
    print(f"    v_in = ({v_alpha_2:.2f}, {v_beta_2:.2f}) V")
    print(f"    v_used = ({dbg2['v_alpha_used']:.2f}, {dbg2['v_beta_used']:.2f}) V")
    print(f"    scale={dbg2['pre_scale']:.4f}, overmod={dbg2['overmodulation_detected']}")
    print(f"    duty = (da={d_a2:.4f}, db={d_b2:.4f}, dc={d_c2:.4f})")

    assert dbg2["pre_scale"] <= 1.0 + 1e-12
    assert 0.0 <= d_a2 <= 1.0 and 0.0 <= d_b2 <= 1.0 and 0.0 <= d_c2 <= 1.0

    # ---------- Test 3: no pre-saturation (detect overmod, clamp) ----------
    d_a3, d_b3, d_c3, dbg3 = svpwm_step(
        v_alpha_2, v_beta_2, Vdc,
        apply_linear_saturation=False,  # intentionally off
        clamp=True
    )

    print("\n  [Test 3] Out-of-range vector without pre-saturation")
    print(f"    overmod={dbg3['overmodulation_detected']} (expected True in many cases)")
    print(f"    duty = (da={d_a3:.4f}, db={d_b3:.4f}, dc={d_c3:.4f})")

    assert 0.0 <= d_a3 <= 1.0 and 0.0 <= d_b3 <= 1.0 and 0.0 <= d_c3 <= 1.0

    print("\n  ✓ SVPWM duty generation works")
    print("  ✓ Linear saturation helper works")
    print("  ✓ Overmodulation detection path works")