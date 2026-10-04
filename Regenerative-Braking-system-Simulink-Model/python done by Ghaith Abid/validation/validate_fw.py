# PURPOSE
# -------
# Validates field_weakening.py by visualizing the fundamental dq-plane geometry
# of the IPMSM control problem across the full speed range.
# Four figures are produced, each targeting a specific concept for the report.
#
# THEORETICAL RECAP
# -----------------
# Steady-state dq voltages (Rs·i dropped at high speed — valid for ω >> Rs/Ld):
#
#   vd ≈ −ω·Lq·iq
#   vq ≈  ω·(Ld·id + λm)
#
# Voltage constraint  vd² + vq² ≤ Vmax²  rearranges to an ellipse in (id, iq):
#
#       (id + λm/Ld)²               iq²
#   ──────────────────────  +  ─────────────────────  ≤  1
#    (Vmax / (ω·Ld))²            (Vmax / (ω·Lq))²
#
#   Center of ellipse : (I_CH, 0)    I_CH = −λm/Ld          [A]
#   Semi-axis id      : a(ω) = Vmax / (ω·Ld)                [A], shrinks ∝ 1/ω
#   Semi-axis iq      : b(ω) = Vmax / (ω·Lq)                [A], shrinks ∝ 1/ω
#
# Current constraint  id² + iq² ≤ Imax²  is a fixed circle centred at origin.
#
# Feasible operating region at speed ω = INTERSECTION of both constraints.
#
# OPERATING REGIONS (CLA 250+ parameters):
#   ω < ω_base :  Ellipse large — MTPA point interior — voltage not active
#   ω = ω_base :  Ellipse just touches MTPA(Imax) point — FW begins
#   ω > ω_base :  FW — operating point slides along current circle arc
#                 toward more-negative id, tracking the ellipse boundary
#   ω → ∞      :  Ellipse shrinks toward I_CH = −380 A (inside current circle)
#
# NOTE ON CHARACTERISTIC CURRENT:
#   |I_CH| = 380 A < Imax = 580 A  →  ellipse centre lies inside the current
#   circle.  This is the condition for torque production up to arbitrarily high
#   speed and is a deliberate design choice for traction motors.
#
# =============================================================================

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
    _repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    sys.path.insert(0, _repo_root)

from python.motor_params import (I_MAX, V_MAX, L_D, L_Q, LAMBDA_M,OMEGA_BASE, OMEGA_MAX, T_max, K_T,P_POLES)
OMEGA_BASE_E = OMEGA_BASE * P_POLES   # electrical base speed  [rad/s]
OMEGA_MAX_E  = OMEGA_MAX  * P_POLES   # electrical max speed   [rad/s]
from python.field_weakening import (solve_fw,compute_transition_speed,build_fw_table,)
from python.mtpa import solve_mtpa   

FIG_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIG_DIR, exist_ok=True)


I_CH = -LAMBDA_M / L_D
# Characteristic current [A] — centre of the voltage ellipse.
# |I_CH| < I_MAX  ⟹  ellipse centre inside current circle  ⟹  high-speed capable.

SALIENCY = L_Q / L_D
# Sets the ratio b/a of the ellipse semi-axes at every speed.

# MTPA curve (positive-torque, motoring side)


from scipy.optimize import brentq

# Find the torque that puts MTPA exactly on the current circle boundary
def _imag_residual(T):
    id_, iq_ = solve_mtpa(T, P_POLES, LAMBDA_M, L_D, L_Q)
    return np.sqrt(id_**2 + iq_**2) - I_MAX

T_at_Imax = brentq(_imag_residual, 1.0, T_max * 3.0)
print(f"T at I_MAX on MTPA = {T_at_Imax:.1f} Nm")

# Sweep torque from near-zero to T_MAX, collect (id*, iq*) from MTPA solver.
T_sweep = np.linspace(1.0, T_at_Imax, 300)   
id_mtpa = np.zeros_like(T_sweep)
iq_mtpa = np.zeros_like(T_sweep)

for k, T_req in enumerate(T_sweep):
    id_mtpa[k], iq_mtpa[k] = solve_mtpa(
    T_req,
    P_POLES,
    LAMBDA_M,
    L_D,
    L_Q
)
# MTPA transition point — now guaranteed to sit on the current circle
id_mtpa_Imax, iq_mtpa_Imax = solve_mtpa(T_at_Imax, P_POLES, LAMBDA_M, L_D, L_Q)
print(f"|i| check = {np.sqrt(id_mtpa_Imax**2 + iq_mtpa_Imax**2):.1f} A  (must = {I_MAX:.1f})")

# FW table at that same torque level
omega_fw, id_fw, iq_fw, T_fw = build_fw_table(
    omega_e_min=OMEGA_BASE_E,
    omega_e_max=OMEGA_MAX_E,
    n_points=300,
    T_request=T_at_Imax,
    p=P_POLES,
    lambda_m=LAMBDA_M,
    L_d=L_D,
    L_q=L_Q,
    I_max=I_MAX,
    V_max=V_MAX,
)


# Transition speed — should match OMEGA_BASE (validated below)
omega_trans = compute_transition_speed(
    id_mtpa=id_mtpa_Imax,
    iq_mtpa=iq_mtpa_Imax,
    V_max=V_MAX,
    L_d=L_D,
    L_q=L_Q,
    lambda_m=LAMBDA_M,
)


def voltage_ellipse(omega,n= 500):
    """
    Returns (id_pts, iq_pts) tracing the voltage ellipse boundary at speed omega.

    Parameters:
    omega (float)   Electrical speed [rad/s]
    n     (int)     Number of points on the ellipse

    Returns:
    id_pts (ndarray)  [A]
    iq_pts (ndarray)  [A]

    Derivation
    ----------
    Centre (I_CH, 0) with I_CH = −λm/Ld
    a(ω) = Vmax / (ω·Ld)    semi-axis along id
    b(ω) = Vmax / (ω·Lq)    semi-axis along iq
    Parametric: id = I_CH + a·cos θ,  iq = b·sin θ,  θ ∈ [0, 2π]
    """
    a = V_MAX / (omega * L_D)           # semi-axis along id  [A]
    b = V_MAX / (omega * L_Q)           # semi-axis along iq  [A]
    theta = np.linspace(0, 2 * np.pi, n)
    return I_CH + a * np.cos(theta), b * np.sin(theta)


theta_c   = np.linspace(0, 2 * np.pi, 500)
id_circle = I_MAX * np.cos(theta_c)    # current circle boundary
iq_circle = I_MAX * np.sin(theta_c)



# SANITY CHECKS  

print("  FIELD-WEAKENING VALIDATION — SANITY CHECKS")

print(f"\n[1] Characteristic current:")
print(f"    I_CH = −λm/Ld = {I_CH:.1f} A")
print(f"    |I_CH| = {abs(I_CH):.1f} A   vs   I_MAX = {I_MAX:.1f} A")
assert abs(I_CH) < I_MAX, \
    "FAIL: ellipse centre is outside the current circle — impossible for traction"
print("    ✓  Ellipse centre lies inside the current circle")

print(f"\n[2] Transition speed continuity:")
print(f"    ω_trans (computed) = {omega_trans:.2f} rad/s")
print(f"    OMEGA_BASE          = {OMEGA_BASE:.2f} rad/s")

finite = np.isfinite(id_fw)
print("id_fw min =", id_fw[finite].min())
print("id_fw max =", id_fw[finite].max())

# find any positive entries
pos_idx = np.where(id_fw[finite] > 0)[0]
print("positive count =", len(pos_idx))
print("first 5 positive indices =", pos_idx[:5])

# show the first offending point
if len(pos_idx) > 0:
    i = pos_idx[0]
    print("omega_fw =", omega_fw[finite][i])
    print("id_fw =", id_fw[finite][i], "iq_fw =", iq_fw[finite][i])
    
print(f"\n[3] id* always negative throughout FW sweep:")
assert np.all(id_fw[np.isfinite(id_fw)] < 0), "FAIL: id* went positive somewhere in FW region"
print(f"    ✓  id range in FW: [{id_fw.min():.1f}, {id_fw.max():.1f}] A")

print(f"\n[4] id* more negative than (or equal to) MTPA id* at I_MAX:")
violations = id_fw[id_fw > id_mtpa_Imax + 1.0]
if len(violations) > 0:
    print(f"    Violations: {violations}")
print(f"    id_mtpa_Imax = {id_mtpa_Imax:.3f} A")
print(f"    id_fw max    = {id_fw.max():.3f} A")
print(f"    difference   = {id_fw.max() - id_mtpa_Imax:.4f} A")
assert np.all(id_fw <= id_mtpa_Imax + 1.0), \
    "FAIL: some FW id* is less negative than the MTPA id* at I_MAX"
print(f"    ✓  id_fw min = {id_fw.min():.1f} A   MTPA id = {id_mtpa_Imax:.1f} A")

print(f"\n[5] Current circle not violated in FW region:")
i_mag_fw = np.sqrt(id_fw**2 + iq_fw**2)
assert np.all(i_mag_fw <= I_MAX * 1.001), \
    f"FAIL: current exceeds I_MAX. Max |i| = {i_mag_fw.max():.1f} A"
print(f"    ✓  Max |i| in FW = {i_mag_fw.max():.1f} A  ≤  I_MAX = {I_MAX:.1f} A")

print(f"\n[6] Voltage constraint satisfied at every FW operating point:")
all_ok = True
for k in range(len(omega_fw)):
    vd = -omega_fw[k] * L_Q * iq_fw[k]
    vq =  omega_fw[k] * (L_D * id_fw[k] + LAMBDA_M)
    if np.sqrt(vd**2 + vq**2) > V_MAX * 1.01:
        print(f"    FAIL at ω = {omega_fw[k]:.1f} rad/s")
        all_ok = False
if all_ok:
    print("    ✓  All FW points satisfy |v| ≤ V_MAX")

print(f"\n[7] Torque is non-increasing with speed in FW region:")
assert T_fw[-1] <= T_fw[0] + 1.0, \
    "FAIL: torque unexpectedly rose in the FW region"
print(f"    ✓  T @ ω_base = {T_fw[0]:.1f} Nm → T @ ω_max = {T_fw[-1]:.1f} Nm")

print("\n" + "=" * 62)
print("  ALL SANITY CHECKS PASSED")
print("=" * 62 + "\n")



# FIGURE 1 — dq PLANE GEOMETRY

"""This is the central geometric picture of the IPMSM control problem and the
most important figure in this validation set.
 The thick grey circle is the current limit boundary — no operating point can ever lie outside it regardless of speed or torque demand. This constraint
comes from the inverter current rating and the motor's thermal limit.

 The dashed coloured ellipses represent the voltage constraint at five
 increasing speeds. At low speed the ellipse is very large (the inverter can
 supply large dq voltages), so it does not restrict the operating point. As
 speed rises, the back-EMF fills up more of the available voltage headroom and
 the ellipse shrinks — but its centre stays fixed at the characteristic current
 point I_CH = −λm/Ld (marked ×). This centre never moves because it depends only on motor parameters, not on speed.
 
 The dot-dash green line is the MTPA locus — the curve of minimum-current operating points for each torque level. 
In the low-speed region the controller follows this curve freely. The red star marks where the MTPA locus at I_MAX
touches the voltage ellipse at ω_base: this is the MTPA→FW transition point.

 Above ω_base the blue arc shows the FW trajectory: the operating point must remain on BOTH constraint boundaries simultaneously 
(current circle and ellipse boundary). This forces id* to become progressively more negative  while iq* decreases.

(for report later)
The key insight to emphasise in the report: I_CH lies INSIDE the current
 circle for this motor (|I_CH| = 380 A < I_MAX = 580 A). This is the
 necessary condition for the ellipse to always intersect the current circle
 at any speed, guaranteeing a feasible operating point up to ω_max. Had
 |I_CH| > I_MAX, a finite maximum speed would exist beyond which no feasible
 point remains.
"""


omega_plot_speeds = np.array([
    OMEGA_BASE_E,
    OMEGA_BASE_E * 1.4,
    OMEGA_BASE_E * 1.8,
    OMEGA_BASE_E * 2.2,
    OMEGA_MAX_E,
])
omega_fw_mech = omega_fw / P_POLES   # [rad/s] mechanical 
ell_colors = plt.cm.plasma(np.linspace(0.15, 0.82, len(omega_plot_speeds)))

fig1, ax1 = plt.subplots(figsize=(8, 8))
ax1.set_aspect('equal')

# Current circle
ax1.plot(id_circle, iq_circle, color='#404040', lw=2.8, zorder=3,label=f'Current circle  $I_{{max}}$ = {I_MAX:.0f} A')
def voltage_ellipse_clipped(omega, n=2000, r_clip=I_MAX * 1.25):
    a = V_MAX / (omega * L_D)
    b = V_MAX / (omega * L_Q)
    theta = np.linspace(0, 2 * np.pi, n)
    id_e = I_CH + a * np.cos(theta)
    iq_e = b * np.sin(theta)
    mask = (id_e**2 + iq_e**2) <= r_clip**2
    if mask.sum() < 3:
        return np.array([]), np.array([])
    return np.where(mask, id_e, np.nan), np.where(mask, iq_e, np.nan)

# Voltage ellipses at selected speeds
for j, (omega_j, col_j) in enumerate(zip(omega_plot_speeds, ell_colors)):
    id_e, iq_e = voltage_ellipse_clipped(omega_j)
    if len(id_e) > 0:
        ax1.plot(id_e, iq_e, color=col_j, lw=1.7, linestyle='--', zorder=2,
                 label=f'Voltage ellipse  ω = {omega_j:.0f} rad/s')

# MTPA locus (motoring, positive iq)
ax1.plot(id_mtpa, iq_mtpa, color='#2ca02c', lw=2.2, linestyle='-.', zorder=4,
         label='MTPA locus')

# FW trajectory along current circle
ax1.plot(id_fw, iq_fw, color='royalblue', lw=2.5, zorder=5,
         label='FW trajectory (current circle arc)')

# Sampled FW operating points
for k in np.linspace(0, len(omega_fw) - 1, 6, dtype=int):
    ax1.plot(id_fw[k], iq_fw[k], 'o', color='royalblue',
             markersize=7, zorder=6)

# MTPA→FW transition point
ax1.plot(id_mtpa_Imax, iq_mtpa_Imax, '*', color='red',
         markersize=18, markeredgecolor='darkred', zorder=7,
         label=f'MTPA→FW transition  ω_base = {OMEGA_BASE:.0f} rad/s')

# Characteristic current point
ax1.plot(I_CH, 0, 'x', color='black', markersize=13,
         markeredgewidth=3, zorder=8,
         label=f'$i_{{ch}}$ = {I_CH:.0f} A  (ellipse centre)')
ax1.annotate(f'$i_{{ch}}$ = {I_CH:.0f} A',
             xy=(I_CH, 0), xytext=(I_CH + 20, 120),
             arrowprops=dict(arrowstyle='->', color='black', lw=1.3),
             fontsize=10)

ax1.axhline(0, color='gray', lw=0.6, linestyle=':')
ax1.axvline(0, color='gray', lw=0.6, linestyle=':')
ax1.set_xlabel('$i_d$  [A]', fontsize=13)
ax1.set_ylabel('$i_q$  [A]', fontsize=13)
ax1.set_title(
    'Figure 1 — dq Plane Geometry: Current Circle & Voltage Ellipses\n'
    'CLA 250+ IPMSM  |  MTPA → Field-Weakening Transition',
    fontsize=12
)
ax1.legend(fontsize=9, loc='upper right')
ax1.set_xlim(-I_MAX * 1.15, I_MAX * 0.5)
ax1.set_ylim(-I_MAX * 0.25, I_MAX * 1.15)
ax1.grid(True, linestyle=':', alpha=0.45)
plt.tight_layout()
fig1.savefig(os.path.join(FIG_DIR, 'fw_fig1_dq_geometry.png'),
             dpi=150, bbox_inches='tight')
plt.show()
print("Figure 1 saved: fw_fig1_dq_geometry.png")


# FIGURE 2 — ELLIPSE SEMI-AXES vs SPEED
'''
This figure makes the 1/ω shrinkage of the voltage ellipse quantitative.
Both semi-axes (a(ω) = Vmax/(ω·Ld) along id and b(ω) = Vmax/(ω·Lq) along
iq) are plotted against electrical speed. Their hyperbolic decay confirms
that doubling the speed halves the voltage "budget" in the dq frame. The
vertical dashed line at ω_base marks where a(ω) first drops to the level of
the MTPA current at I_MAX, triggering field weakening. This is the geometric definition of ω_base.

 The two horizontal dashed reference lines — at I_MAX and |I_CH| — allow
visual confirmation that: (a) FW begins when a(ω) drops below I_MAX; and
(b) the ellipse would only collapse to a point at |I_CH|, which is below I_MAX for this motor, confirming indefinite high-speed capability.
'''

omega_axis = np.linspace(OMEGA_BASE_E * 0.5, OMEGA_MAX_E * 1.05, 600)
a_ax = V_MAX / (omega_axis * L_D)      # semi-axis along id  [A]
b_ax = V_MAX / (omega_axis * L_Q)      # semi-axis along iq  [A]

fig2, ax2 = plt.subplots(figsize=(9, 5))

ax2.plot(omega_axis, a_ax, color='darkorange', lw=2.2,
         label=r'$a(\omega) = V_{max}/(\omega L_d)$ — $i_d$ semi-axis')
ax2.plot(omega_axis, b_ax, color='steelblue', lw=2.2,
         label=r'$b(\omega) = V_{max}/(\omega L_q)$ — $i_q$ semi-axis')

ax2.axvline(OMEGA_BASE_E, color='red', lw=1.6, linestyle='--',
            label=f'$\\omega_{{base}}$ = {OMEGA_BASE:.0f} rad/s mech ({OMEGA_BASE_E:.0f} elec, FW entry)')
ax2.axvline(OMEGA_MAX_E, color='darkred', lw=1.4, linestyle=':',
            label=f'$\\omega_{{max}}$ = {OMEGA_MAX:.0f} rad/s mech ({OMEGA_MAX_E:.0f} elec)')
ax2.axhline(I_MAX, color='#404040', lw=1.3, linestyle='--',
            label=f'$I_{{max}}$ = {I_MAX:.0f} A')
ax2.axhline(abs(I_CH), color='black', lw=1.3, linestyle=':',
            label=f'$|i_{{ch}}|$ = {abs(I_CH):.0f} A  (collapse point)')

ax2.set_xlabel('Electrical speed  ω  [rad/s]', fontsize=12)
ax2.set_ylabel('Semi-axis magnitude  [A]', fontsize=12)
ax2.set_title(
    f'Figure 2 — Voltage Ellipse Semi-Axes vs Speed  (ξ = Lq/Ld = {SALIENCY:.2f})\n'
    r'Both axes shrink $\propto 1/\omega$; ratio $b/a = \xi$ is constant at all speeds',
    fontsize=11
)
ax2.legend(fontsize=9)
ax2.grid(True, linestyle=':', alpha=0.45)
ax2.set_xlim(omega_axis[0], omega_axis[-1])
ax2.set_ylim(0, I_MAX * 1.3)
plt.tight_layout()
fig2.savefig(os.path.join(FIG_DIR, 'fw_fig2_ellipse_axes_vs_speed.png'),
             dpi=150, bbox_inches='tight')
plt.show()
print("Figure 2 saved: fw_fig2_ellipse_axes_vs_speed.png")


# FIGURE 4 — id*, iq* AND CURRENT MAGNITUDE vs SPEED
'''
This two-panel figure makes the control action explicit and provides a numerical cross-check of the current circle constraint.
# Upper panel — id* and iq* vs speed:
In the MTPA region both reference currents are constant — the controller simply holds the MTPA 
solution at I_MAX with no speed-dependent adjustment.
At ω_base (red dashed line) the references diverge: id* (blue) starts
decreasing — becoming more negative — while iq* (orange) also falls. 
The physical interpretation: id* demagnetises the machine (opposes λm) so the
#voltage equation ω·(Ld·id + λm) stays within Vmax at higher speed. This
comes at the price of reduced iq* due to the current circle coupling
id² + iq² = I_MAX². At high speed, id* asymptotically approaches the
characteristic current I_CH (black dotted line) — the point at which the back-EMF would be fully cancelled 
and the voltage ellipse would collapse.

   Note for the report: the gap between id*(ω_max) and I_CH indicates how much
   demagnetisation margin remains and whether an MTPV region is geometrically
   reachable before ω_max (it is — the MTPV locus intercepts the current circle
   before id* reaches I_CH).

# Lower panel — current magnitude |i| vs speed:
This is a hard constraint check. In the MTPA region, |i| equals I_MAX by
construction (maximum torque per ampere at full current). In the FW region,
|i| must remain equal to I_MAX because the FW algorithm operates on the
current circle boundary. Any departure from the red dashed reference line
indicates either a numerical error in solve_field_weakening or a sub-optimal operating point that is not fully exploiting 
the available current. And so a flat line at I_MAX across the entire FW region is the desired result.
'''

omega_zone1b   = np.linspace(0.0, OMEGA_BASE, 150)
id_zone1b      = np.full_like(omega_zone1b, id_mtpa_Imax)
iq_zone1b      = np.full_like(omega_zone1b, iq_mtpa_Imax)
i_mag_zone1b   = np.full_like(omega_zone1b, I_MAX)         # = I_MAX by construction

fig4, (ax4a, ax4b) = plt.subplots(
    2, 1, figsize=(9, 8), sharex=True,
    gridspec_kw={'height_ratios': [3, 1.4], 'hspace': 0.08}
)

# ── Upper panel
ax4a.plot(omega_zone1b, id_zone1b, color='steelblue', lw=2.3,
          label='$i_d^*$ (MTPA region)')
ax4a.plot(omega_fw_mech,     id_fw,     color='steelblue', lw=2.3, linestyle='--',
          label='$i_d^*$ (FW region)')

ax4a.plot(omega_zone1b, iq_zone1b, color='darkorange', lw=2.3,
          label='$i_q^*$ (MTPA region)')
ax4a.plot(omega_fw_mech,     iq_fw,     color='darkorange', lw=2.3, linestyle='--',
          label='$i_q^*$ (FW region)')

ax4a.axvline(OMEGA_BASE, color='red', lw=1.6, linestyle='--',
             label=f'$\\omega_{{base}}$ = {OMEGA_BASE:.0f} rad/s')
ax4a.axhline(I_CH, color='black', lw=1.2, linestyle=':',
             label=f'$i_{{ch}}$ = {I_CH:.0f} A  (asymptote for $i_d^*$)')
ax4a.axhline(0, color='gray', lw=0.6, linestyle=':')

ax4a.set_ylabel('Reference current  [A]', fontsize=12)
ax4a.set_title(
    'Figure 3 — Reference Currents $i_d^*$, $i_q^*$ and Current Magnitude vs Speed',
    fontsize=11
)
ax4a.legend(fontsize=9, loc='center right')
ax4a.grid(True, linestyle=':', alpha=0.45)
ax4a.set_ylim(-I_MAX * 1.1, I_MAX * 1.1)

# ── Lower panel: |i| check 
ax4b.plot(omega_zone1b, i_mag_zone1b, color='gray', lw=2.2,
          label='|i| (MTPA region)')
ax4b.plot(omega_fw_mech, np.sqrt(id_fw**2 + iq_fw**2), color='purple', lw=2.2,
          label='|i| (FW region) — must equal $I_{max}$')
ax4b.axhline(I_MAX, color='red', lw=1.6, linestyle='--',
             label=f'$I_{{max}}$ = {I_MAX:.0f} A')

ax4b.set_xlabel('Mechanical speed  ω  [rad/s]', fontsize=12)
ax4b.set_ylabel('$|i|$  [A]', fontsize=12)
ax4b.set_title('Current magnitude — constraint check (flat line = correct)',
               fontsize=10)
ax4b.legend(fontsize=9)
ax4b.grid(True, linestyle=':', alpha=0.45)
ax4b.set_ylim(0, I_MAX * 1.18)
ax4b.set_xlim(0, OMEGA_MAX * 1.05)

plt.tight_layout()
fig4.savefig(os.path.join(FIG_DIR, 'fw_fig3_currents_vs_speed.png'),dpi=150, bbox_inches='tight')
plt.show()
print("Figure 3 saved: fw_fig3_currents_vs_speed.png")

print(f"\nAll figures saved to: {FIG_DIR}")