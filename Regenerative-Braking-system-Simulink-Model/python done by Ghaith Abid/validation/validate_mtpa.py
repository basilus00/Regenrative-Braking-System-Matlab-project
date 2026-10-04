
# validate_mtpa.py
# Produces 4 plots:
#   1. MTPA curve (i_d* vs i_q*)
#   2. MTPA curve overlaid on current limit circle
#   3. Torque hyperbola with MTPA operating point
#   4. Newton-Raphson convergence


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Path setup
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from motor_params import (P_POLES, LAMBDA_M, L_D, L_Q,T_max, I_MAX, K_T)
from mtpa import solve_mtpa, build_mtpa_table


# 1. Build MTPA table across full torque range

N_POINTS = 401                        
torques, id_table, iq_table = build_mtpa_table(T_max, N_POINTS, P_POLES, LAMBDA_M, L_D, L_Q)

# 2. Prepare figure — 2×2 grid

fig = plt.figure(figsize=(14, 11))
fig.suptitle("MTPA Validation",fontsize=14, fontweight='bold', y=0.98)
gs=gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35)

ax1 = fig.add_subplot(gs[0, 0])   # MTPA curve
ax2 = fig.add_subplot(gs[0, 1])   # current circle
ax3 = fig.add_subplot(gs[1, 0])   # torque hyperbola
ax4 = fig.add_subplot(gs[1, 1])   # convergence

# Plot 1 — MTPA Curve: i_d* vs i_q* with an expected shape: cubic  i_d = C·i_q³

ax1.plot(iq_table, id_table, 'b-', linewidth=2, label='MTPA locus')
ax1.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax1.axvline(0, color='k', linewidth=0.5, linestyle='--')
id_200, iq_200 = solve_mtpa(200.0, P_POLES, LAMBDA_M, L_D, L_Q)# Mark a specific operating point (T=200 Nm)
ax1.plot(iq_200, id_200, 'ro', markersize=8, label='T=+200 Nm (motoring)')
id_m200, iq_m200 = solve_mtpa(-200.0, P_POLES, LAMBDA_M, L_D, L_Q)
ax1.plot(iq_m200, id_m200, 'gs', markersize=8, label='T=−200 Nm (regen)')
ax1.set_xlabel('i_q* (A)')
ax1.set_ylabel('i_d* (A)')
ax1.set_title('Plot 1 — MTPA Curve')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)
ax1.set_facecolor('#f8f8f8')

# Plot 2 — MTPA curve inside current limit circle

theta = np.linspace(0, 2*np.pi, 500)
circle_id = I_MAX * np.sin(theta)
circle_iq = I_MAX * np.cos(theta)
ax2.plot(circle_iq, circle_id, 'r--',linewidth=1.5, label=f'Current limit |Is|={I_MAX}A')
ax2.plot(iq_table, id_table, 'b-',linewidth=2, label='MTPA locus')
ax2.plot(iq_200, id_200, 'ro', markersize=8)
ax2.plot(iq_m200, id_m200, 'gs', markersize=8)
ax2.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax2.axvline(0, color='k', linewidth=0.5, linestyle='--')
ax2.set_xlabel('i_q* (A)')
ax2.set_ylabel('i_d* (A)')
ax2.set_title('Plot 2 — MTPA Inside Current Circle')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_facecolor('#f8f8f8')
ax2.set_aspect('equal')

# Plot 3 — Torque hyperbola at T=200 Nm with MTPA operating point . PS: The MTPA point should be where the smallest current circle is tangent to the torque hyperbola

# Generate torque hyperbola: for each i_d, compute i_q from torque equation
id_range= np.linspace(-I_MAX * 0.8, 0, 300)
denom =LAMBDA_M + (L_D - L_Q) * id_range
valid = np.abs(denom) > 1e-6
iq_hyp = np.zeros_like(id_range)
iq_hyp[valid] = 200.0 / (K_T * denom[valid])
Is_hyp = np.sqrt(id_range**2 + iq_hyp**2)
mask=Is_hyp <= I_MAX * 1.1
ax3.plot(iq_hyp[mask], id_range[mask], 'm-',linewidth=2, label='Torque hyperbola T=200 Nm')

for Is_val in [100, 200, 270, 350, 450]:
    theta_c = np.linspace(0, 2*np.pi, 300)
    ax3.plot(Is_val * np.cos(theta_c),Is_val * np.sin(theta_c),'b-', linewidth=0.6, alpha=0.3)

ax3.plot(iq_200, id_200, 'ro', markersize=10, zorder=5,label=f'MTPA point\n(i_d={id_200:.1f}, i_q={iq_200:.1f})')
ax3.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax3.axvline(0, color='k', linewidth=0.5, linestyle='--')
ax3.set_xlim(-50, I_MAX * 0.8)
ax3.set_ylim(-I_MAX * 0.6, 0)
ax3.set_xlabel('i_q (A)')
ax3.set_ylabel('i_d (A)')
ax3.set_title('Plot 3 — Torque Hyperbola + MTPA Point')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)
ax3.set_facecolor('#f8f8f8')

# Plot 4 — Newton-Raphson convergence (shows how fast the solver reaches the solution)

def solve_mtpa_tracked(T_request, p, lambda_m, L_d, L_q):#Same as solve_mtpa but records error at each iteration.
    K=1.5*p
    dL=L_d-L_q
    a=K*dL**2
    b=lambda_m*T_request
    c=-(T_request**2)/K
    iq=T_request/(K * lambda_m)
    errors = []
    for _ in range(15):
        f=a*iq**4+b*iq+ c
        df=4*a*iq**3+b
        step=f/df
        iq-=step
        errors.append(abs(step))
        if abs(step) < 1e-6:
            break
    return errors
errors_mot = solve_mtpa_tracked(200.0,  P_POLES, LAMBDA_M, L_D, L_Q)
errors_reg = solve_mtpa_tracked(-200.0, P_POLES, LAMBDA_M, L_D, L_Q)

ax4.semilogy(range(1, len(errors_mot)+1), errors_mot,'b-o', linewidth=2, markersize=5, label='T=+200 Nm (motoring)')
ax4.semilogy(range(1, len(errors_reg)+1), errors_reg,'g-s', linewidth=2, markersize=5, label='T=−200 Nm (regen)')
ax4.axhline(1e-6, color='r', linestyle='--',linewidth=1, label='Tolerance = 1e-6')
ax4.set_xlabel('Iteration')
ax4.set_ylabel('Step size |Δi_q| (A) — log scale')
ax4.set_title('Plot 4 — Newton-Raphson Convergence')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3, which='both')
ax4.set_facecolor('#f8f8f8')

# Save and show
output_dir = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'mtpa_validation.png')
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"Figure saved → {output_path}")
plt.show()
''' 
## What to look for exactly in each plot

*Plot 1 — MTPA Curve:
-Shape should be cubic (not linear, not parabolic)
-Both motoring (red dot) and regen (green square) visible
-Symmetric about i_q = 0 axis
-i_d* always negative 

*Plot 2 — Current Circle:
-MTPA locus must stay entirely INSIDE the red circle
-If it crosses outside → current limit violated → problem with I_MAX

**Plot 3 — Torque Hyperbola:
-Red dot should sit exactly where the smallest blue circle is tangent to the purple hyperbola
-This is the geometric proof that MTPA is correct

**Plot 4 — Convergence:
-Should drop below tolerance (red dashed line) in 5-7 iterations
-Both motoring and regen should converge at same rate
-Log scale — each step should be roughly 2x smaller than previous'''