'''
# MTPV Validation 
PURPOSE:
Validates mtpv.py by visualising the full three-region dq geometry in one comprehensive figure.
This is the central picture of the complete FOC reference generation strategy.

# WHAT TO LOOK FOR
This figure shows all three operating regions simultaneously:

1) MTPA locus (green dash-dot):
Starts near the origin, curves toward the current circle. 
Every point on this curve is the minimum-current solution for a given torque.
The locus ends at the MTPA→FW transition point (red star) where the voltage ellipse first touches the current circle at ω_base.

2)FW arc (blue, on current circle boundary):
The operating point is pinned to the current circle and slides toward more negative id* as speed rises.
Both the current circle AND voltage ellipse are active simultaneously.
The arc ends at the FW→MTPV transition point (orange star) where the MTPV locus enters the circle.

3) MTPV locus (red, inside current circle):
Above ω_MTPV the operating point detaches from the current circle and moves inward along the MTPV hyperbola.
Only the voltage ellipse is active.As ω → ∞ the ellipse collapses and the operating point approaches I_CH = (−λm/Ld, 0).

# KEY INSIGHTS 
1. The FW→MTPV transition (orange star) marks where the MTPV locus crosses inside the current circle.
Below this speed the MTPV point is unreachable (outside the circle). Above it, staying on the circle is sub-optimal.

2. The MTPV operating points move INSIDE the current circle — current is being intentionally under-used 
because the voltage ellipse is the only binding constraint. Pushing more current would violate V_max.

3. The voltage ellipses in the MTPV region are entirely contained within the current 
circle , confirming that the current constraint is inactive.


5. As speed increases, the MTPV operating point approaches I_CH = (-λm/Ld, 0) .
At that limit iq* → 0 and T → 0. This is the theoretical infinite-speed limit of the drive.
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq as _brentq

if __name__ == "__main__":
    _repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    sys.path.insert(0, _repo_root)

from python.motor_params import (P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,OMEGA_BASE, OMEGA_MAX, I_CH, K_T, T_max)
from python.mtpa import solve_mtpa
from python.field_weakening import build_fw_table, compute_transition_speed
from python.mtpv import solve_mtpv, compute_mtpv_transition_speed, build_mtpv_table

FIG_DIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# Electrical speeds
OMEGA_BASE_E = OMEGA_BASE * P_POLES
OMEGA_MAX_E  = OMEGA_MAX  * P_POLES

#  Transition speeds between regions
omega_mtpv = compute_mtpv_transition_speed(P_POLES, LAMBDA_M, L_D, L_Q, V_MAX, I_MAX,omega_min=OMEGA_BASE_E)
omega_mtpv_mech = omega_mtpv / P_POLES

#  MTPA locus 
def _imag_residual(T):
    id_, iq_ = solve_mtpa(T, P_POLES, LAMBDA_M, L_D, L_Q)
    return np.sqrt(id_**2 + iq_**2) - I_MAX

T_at_Imax = _brentq(_imag_residual, 1.0, T_max * 3.0)
id_mtpa_Imax, iq_mtpa_Imax = solve_mtpa(T_at_Imax, P_POLES, LAMBDA_M, L_D, L_Q)
T_sweep = np.linspace(1.0, T_at_Imax, 300)
id_mtpa_arr = np.zeros(300)
iq_mtpa_arr = np.zeros(300)
for k, T in enumerate(T_sweep):
    id_mtpa_arr[k], iq_mtpa_arr[k] = solve_mtpa(T, P_POLES, LAMBDA_M, L_D, L_Q)

#FW trajectory 
omega_fw, id_fw, iq_fw, T_fw = build_fw_table(
    omega_e_min=OMEGA_BASE_E,
    omega_e_max=omega_mtpv,
    n_points=200,
    T_request=T_at_Imax,
    p=P_POLES, lambda_m=LAMBDA_M, L_d=L_D, L_q=L_Q,
    I_max=I_MAX, V_max=V_MAX,
)

# FW→MTPV transition point (last FW point = first MTPV point)
id_fw_mtpv = id_fw[-1]
iq_fw_mtpv = iq_fw[-1]

# MTPV trajectory 
omega_mtpv_arr, id_mtpv_arr, iq_mtpv_arr, T_mtpv_arr = build_mtpv_table(
    omega_e_min=omega_mtpv,
    omega_e_max=OMEGA_MAX_E * 2.5,
    n_points=300,
    T_request=T_at_Imax,
    p=P_POLES, lambda_m=LAMBDA_M, L_d=L_D, L_q=L_Q, V_max=V_MAX,
)

# assisit functions
def voltage_ellipse_clipped(omega_e, n=3000, r_clip=I_MAX * 1.05):
    a = V_MAX / (omega_e * L_D)
    b = V_MAX / (omega_e * L_Q)
    theta = np.linspace(0, 2 * np.pi, n)
    id_e = I_CH + a * np.cos(theta)
    iq_e = b * np.sin(theta)
    mask = (id_e**2 + iq_e**2) <= r_clip**2
    if mask.sum() < 3:
        return np.array([]), np.array([])
    return np.where(mask, id_e, np.nan), np.where(mask, iq_e, np.nan)

theta_c   = np.linspace(0, 2*np.pi, 500)
id_circle = I_MAX * np.cos(theta_c)
iq_circle = I_MAX * np.sin(theta_c)

# MTPV ellipse speeds (pick 4 speeds  inside MTPV region) 
omega_mtpv_speeds = np.array([omega_mtpv * 1.1,omega_mtpv * 1.5,omega_mtpv * 2.0,OMEGA_MAX_E * 2.0,])
ell_colors = plt.cm.Reds(np.linspace(0.35, 0.85, len(omega_mtpv_speeds)))

# FIGURE ( 3 region dq Plane Geometry)
fig, ax = plt.subplots(figsize=(9, 9))
ax.set_aspect('equal')

# Current circle
ax.plot(id_circle, iq_circle, color='#303030', lw=2.8, zorder=3,label=f'Current circle  $I_{{max}}$ = {I_MAX:.0f} A')

# Voltage ellipses (MTPV region — entirely inside current circle)
for j, (omega_j, col_j) in enumerate(zip(omega_mtpv_speeds, ell_colors)):
    id_e, iq_e = voltage_ellipse_clipped(omega_j, r_clip=I_MAX * 1.5)
    if len(id_e) > 0:
        ax.plot(id_e, iq_e, color=col_j, lw=1.5, linestyle='--', zorder=2,
                label=f'Voltage ellipse  ω = {omega_j/P_POLES:.0f} rad/s (mech)')
# MTPA locus
ax.plot(id_mtpa_arr, iq_mtpa_arr, color='#2ca02c', lw=2.2,linestyle='-.', zorder=4, label='MTPA locus')
# FW arc
ax.plot(id_fw, iq_fw, color='royalblue', lw=2.5, zorder=5,label='FW trajectory (on current circle)')
# MTPV locus
ax.plot(id_mtpv_arr, iq_mtpv_arr, color='crimson', lw=2.5, zorder=6,label='MTPV locus (inside current circle)')

# Sampled MTPV operating points
for k in np.linspace(0, len(omega_mtpv_arr)-1, 5, dtype=int):
    ax.plot(id_mtpv_arr[k], iq_mtpv_arr[k], 'o',
            color='crimson', markersize=8, zorder=7)

# MTPA→FW transition point
ax.plot(id_mtpa_Imax, iq_mtpa_Imax, '*', color='red',markersize=16, markeredgecolor='darkred', zorder=8,
        label=f'MTPA→FW  ω_base = {OMEGA_BASE:.0f} rad/s')
# FW→MTPV transition point
ax.plot(id_fw_mtpv, iq_fw_mtpv, '*', color='darkorange',
        markersize=16, markeredgecolor='saddlebrown', zorder=8,
        label=f'FW→MTPV  ω = {omega_mtpv_mech:.0f} rad/s')
# Characteristic current point
ax.plot(I_CH, 0, 'x', color='black', markersize=14,
        markeredgewidth=3, zorder=9,
        label=f'$i_{{ch}}$ = {I_CH:.0f} A  (MTPV asymptote)')
ax.annotate(f'$i_{{ch}}$ = {I_CH:.0f} A',
            xy=(I_CH, 0), xytext=(I_CH - 20, 120),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.3),
            fontsize=10)

ax.axhline(0, color='gray', lw=0.6, linestyle=':')
ax.axvline(0, color='gray', lw=0.6, linestyle=':')

ax.set_xlabel('$i_d$  [A]', fontsize=13)
ax.set_ylabel('$i_q$  [A]', fontsize=13)
ax.set_title(
    'Full Three-Region dq Plane Geometry \n'
    'MTPA  →  Field Weakening  →  MTPV',
    fontsize=12
)
ax.legend(fontsize=8.5, loc='upper right')
ax.set_xlim(-I_MAX * 1.3, I_MAX * 0.5)
ax.set_ylim(-I_MAX * 0.25, I_MAX * 1.2)
ax.grid(True, linestyle=':', alpha=0.4)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, 'mtpv_fig1_dq_geometry.png'),dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved: mtpv_fig1_dq_geometry.png")