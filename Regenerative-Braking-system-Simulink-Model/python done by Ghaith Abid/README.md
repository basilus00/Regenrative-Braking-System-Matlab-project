# Regenerative Braking System — Mercedes CLA 250+

## Overview
Full simulation of a Regenerative Braking System (RBS) for the Mercedes-Benz CLA 250+
electric vehicle platform. The project implements a complete Field-Oriented Control (FOC)
pipeline for an Interior Permanent Magnet Synchronous Motor (IPMSM), developed in Python
for algorithm validation and translated to MATLAB/Simulink for final simulation.

The control architecture covers three operating regions:
- **MTPA** (Maximum Torque Per Ampere) — low speed, current-limited
- **Field Weakening (FW)** — intermediate speed, both constraints active
- **MTPV** (Maximum Torque Per Volt) — high speed, voltage-limited

**Motor:** Interior Permanent Magnet Synchronous Motor (IPMSM)  
**Platform:** Mercedes-Benz CLA 250+ (706V DC bus, 200kW, 335Nm)  
**Group:** Group 2

---

## Project Structure
```
RBS/
├── python/
│   ├── motor_params.py          # Single source of truth for ALL motor parameters
│   │                            # Contains two embedded solver functions:
│   │                            #   solve_consistent_params(λm, ξ, target_gap) → Ld, Lq
│   │                            #   verify_params(Ld, Lq, λm, ...) → consistency report
│   │                            # Never hardcode parameters anywhere else
│   │
│   ├── clarke_park.py           # Amplitude-invariant Clarke transform
│   │                            # Park and Inverse Park transforms
│   │                            # Only ia, ib measured — ic reconstructed via KCL
│   │
│   ├── mtpa.py                  # MTPA reference generator (low speed)
│   │                            # Newton-Raphson on degree-4 polynomial
│   │                            # Returns (id*, iq*) minimising |i| for given T
│   │
│   ├── field_weakening.py       # FW reference generator (intermediate speed)
│   │                            # Closed-form quadratic — intersection of
│   │                            # current circle AND voltage ellipse
│   │                            # Returns (id*, iq*, T_achieved)
│   │                            # Includes compute_transition_speed() for MTPA→FW boundary
│   │
│   ├── mtpv.py                  # MTPV reference generator (high speed)
│   │                            # Lagrange multipliers + w = Ld·id + λm substitution
│   │                            # Quadratic always has positive discriminant
│   │                            # Returns (id*, iq*, T_achieved)
│   │                            # Includes compute_mtpv_transition_speed() for FW→MTPV boundary
│   │
│   ├── foc_references.py        # Master switcher — brain of the reference layer
│   │                            # initialise() called once at startup → computes ω_MTPV
│   │                            # compute_references() called every control cycle
│   │                            # Only file that imports from all three solvers
│   │                            # Returns (id*, iq*, T_achieved, active_region)
│   │                            #   active_region: 1=MTPA, 2=FW, 3=MTPV
│   │
│   ├── pi_controller.py         # Inner current control loop
│   │                            # Full dq decoupling feedforward
│   │                            # Anti-windup + voltage saturation on (vd*,vq*) circle
│   │                            # Gains: Kp = L·ωc, Ki = Rs·ωc (pole-zero cancellation)
│   │
│   ├── svpwm.py                 # Space Vector PWM
│   │                            # (vα*, vβ*) → duty cycles (Da, Db, Dc) ∈ [0,1]
│   │                            # Works identically for motoring AND regen
│   │                            # Power direction set by iq* sign, not SVPWM logic
│   │
│   └── validation/
│       ├── validate_mtpa.py     # 4 plots: MTPA curve, current circle,
│       │                        # torque hyperbola, NR convergence
│       ├── validate_fw.py       # 3 plots: dq geometry (shrinking ellipses),
│       │                        # ellipse semi-axes vs speed,
│       │                        # id*/iq*/|i| vs speed
│       ├── validate_mtpv.py     # 1 plot: full three-region dq geometry
│       │                        # MTPA locus → FW arc → MTPV locus → I_CH
│       └── figures/             # Auto-generated .png files (not tracked by Git)
│
├── matlab/
│   ├── functions/               # MATLAB Function blocks
│   │   ├── solve_mtpa.m         ← translate after mtpa.py validated
│   │   ├── solve_fw.m           ← translate after field_weakening.py validated
│   │   ├── solve_mtpv.m         ← translate after mtpv.py validated
│   │   └── pi_controller.m      ← translate after pi_controller.py validated
│   │
│   ├── params/
│   │   └── motor_params.m       # Same values as motor_params.py, MATLAB syntax
│   │                            # Must be kept in sync with Python version
│   │
│   ├── models/
│   │   └── RBS_main.slx         # Main Simulink model
│   │                            # Full FOC pipeline:
│   │                            # sensors → Clarke/Park → foc_references
│   │                            # → PI → Inverse Park → SVPWM → inverter
│   │
│   └── scripts/
│       └── run_simulation.m     # Entry point: loads params, runs model, plots
│
└── docs/
    └── weekly_reports/          # Weekly progress reports (Word documents)
```

---

## Algorithm Status

| Algorithm | File | Status | Validation |
|-----------|------|--------|------------|
| Clarke + Park transforms | `clarke_park.py` | ✅ Done | ✅ Validated |
| MTPA solver | `mtpa.py` | ✅ Done | ✅ Validated |
| Field Weakening | `field_weakening.py` | ✅ Done | ✅ Validated |
| MTPV solver | `mtpv.py` | ✅ Done | ✅ Validated |
| FOC Reference switcher | `foc_references.py` | ✅ Done | ✅ Validated |
| PI Controller + Decoupling | `pi_controller.py` | ✅ Done | ⏳ Pending |
| SVPWM | `svpwm.py` | ✅ Done | ⏳ Pending |
| MATLAB translation | `matlab/functions/` | ⏳ Pending | ⏳ Pending |
| Simulink model | `RBS_main.slx` | ⏳ Pending | ⏳ Pending |
| Battery model | — | ⏳ Pending | ⏳ Pending |

---

## Control Architecture

```
T_request, ω_e
      │
      ▼
┌──────────────────────────────────┐
│       foc_references.py          │  ← master switcher
│                                  │
│  ω ≤ ω_base(T)    → MTPA        │  minimise |i| for given T
│  ω_base < ω ≤ ω_MTPV → FW      │  current circle ∩ voltage ellipse
│  ω > ω_MTPV       → MTPV        │  maximise T on voltage ellipse
│                                  │
│  returns: id*, iq*, T_ach, region│
└────────────┬─────────────────────┘
             │ id*, iq*
             ▼
┌─────────────────────┐
│   pi_controller.py  │  ← inner current loop
│                     │    vd* = Kp·ed + Ki·∫ed − ωe·Lq·iq
│   + dq decoupling   │    vq* = Kp·eq + Ki·∫eq + ωe·Ld·id + ωe·λm
│   + anti-windup     │    Gains: Kp=L·ωc, Ki=Rs·ωc
│   + V saturation    │
└────────┬────────────┘
         │ vd*, vq*
         ▼
┌─────────────────────┐
│  inverse_park.py    │  ← back to stationary frame
│                     │    vα* = vd·cos θ − vq·sin θ
│                     │    vβ* = vd·sin θ + vq·cos θ
└────────┬────────────┘
         │ vα*, vβ*
         ▼
┌─────────────────────┐
│     svpwm.py        │  → Da, Db, Dc ∈ [0,1]
│                     │    sector identification → dwell times → duty cycles
│                     │    identical for motoring and regen
└─────────────────────┘
         │
         ▼
    3-phase inverter → IPMSM
         │
         ▼ (feedback)
    ia, ib measured → ic = -ia-ib → Clarke → Park → id_meas, iq_meas
```

---

## How to Run

### Requirements
```bash
pip install numpy matplotlib scipy
```

### Verify motor parameters
```bash
python python/motor_params.py
```

### Run individual algorithm sanity checks
```bash
python python/mtpa.py
python python/field_weakening.py
python python/mtpv.py
python python/foc_references.py
python python/pi_controller.py
python python/svpwm.py
```

### Run validation plots
```bash
python python/validation/validate_mtpa.py
python python/validation/validate_fw.py
python python/validation/validate_mtpv.py
```

---

## Motor Parameters (CLA 250+)

| Parameter | Symbol | Value | Source |
|-----------|--------|-------|--------|
| Peak power | P_max | 200 kW | Official spec |
| Peak torque | T_max | 335 Nm | Official spec |
| DC bus voltage | V_dc | 706 V | Official spec |
| Max phase voltage | V_max | 407.6 V | Derived (V_dc/√3) |
| Gear ratio 1 | — | 11:1 | Official spec |
| Gear ratio 2 | — | 5:1 | Official spec |
| Pole pairs | p | 4 | Estimated |
| Stator resistance | R_s | 0.012 Ω | Estimated |
| d-axis inductance | L_d | 250 µH | Solver-derived* |
| q-axis inductance | L_q | 600 µH | Solver-derived* |
| PM flux linkage | λ_m | 0.095 Wb | Estimated |
| Max phase current | I_max | 580 A | Estimated |
| Base speed | ω_base | 576 rad/s (5500 RPM) | Derived (P/T) |
| Max speed | ω_max | 1466 rad/s (14000 RPM) | Estimated |
| Characteristic current | I_ch | −380 A | Derived (−λm/Ld) |
| Saliency ratio | ξ | 2.4 | Derived (Lq/Ld) |
| MTPV transition speed | ω_MTPV | 1047 rad/s mech | Computed numerically |

*L_d and L_q computed via `solve_consistent_params()` in `motor_params.py`
to guarantee a physically valid FW arc: `id_mtpa(I_max) > I_ch`.
True values proprietary to Mercedes-Benz.

---

## Key Physics Rules (never violate)
- `id*` is **always negative** — in MTPA, FW, MTPV, motoring and regen
- `iq*` sign follows `T_request` sign only (+motoring, −regen)
- `L_q > L_d` always — saliency condition for IPMSM
- `|I_ch| < I_max` always — guarantees high-speed capability
- Clarke transform uses **amplitude-invariant** form
- Only `ia` and `ib` measured — `ic = −ia − ib` always reconstructed
- `id* ≥ I_ch` in FW — cannot over-cancel PM flux beyond characteristic current
- `id*` **can** go below `I_ch` in MTPV — operating on ellipse boundary, not current circle
- All solvers work in **electrical** rad/s — divide by `P_POLES` for mechanical

## Speed Unit Convention
```
motor_params.py exports OMEGA_BASE, OMEGA_MAX in MECHANICAL rad/s
All solvers expect ELECTRICAL rad/s

OMEGA_BASE_E = OMEGA_BASE * P_POLES   # use this when calling solvers
omega_mech   = omega_e   / P_POLES   # use this when plotting
```

---

## Known Design Notes
- FW arc is narrow (~50 A) due to estimated parameters. Motor transitions
  MTPA → FW → MTPV quickly above ω_base. Documented as expected behaviour
  given proprietary parameters are unavailable.
- Voltage feasibility check returns `feasible=False` at some FW test points
  because `id*` is clamped to `I_ch`. This is physically correct.
- SVPWM works identically for motoring and regen. Power direction is determined
  entirely by the sign of `iq*`, not by any change in switching logic.

---

## Team
- Ghaith Abid
- Maissa Lajmi
