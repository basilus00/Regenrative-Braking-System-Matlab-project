function [d_a, d_b, d_c] = svpwm_step(v_alpha, v_beta, Vdc, apply_saturation)
% SVPWM_STEP - Space Vector PWM duty cycle generation
%
% Inputs:
%   v_alpha          : Alpha-axis voltage reference [V]
%   v_beta           : Beta-axis voltage reference [V]
%   Vdc              : DC bus voltage [V]
%   apply_saturation : 1 = scale vector to linear region, 0 = skip
%
% Outputs:
%   d_a    : Phase A duty cycle [0 to 1]
%   d_b    : Phase B duty cycle [0 to 1]
%   d_c    : Phase C duty cycle [0 to 1]
%   overmod: Overmodulation flag (1 = exceeded linear region)
%   scale  : Scaling factor applied (1.0 = no scaling)

    %% Input protection
    if Vdc <= 0
        error('Vdc must be greater than 0');
    end

    %% Step 1 — Voltage Saturation (optional)
    v_lim = Vdc / sqrt(3);
    mag   = sqrt(v_alpha^2 + v_beta^2);

    if apply_saturation && mag > v_lim && mag > 1e-12
        scale       = v_lim / mag;
        v_alpha_use = v_alpha * scale;
        v_beta_use  = v_beta  * scale;
    else
        scale       = 1.0;
        v_alpha_use = v_alpha;
        v_beta_use  = v_beta;
    end

    %% Step 2 — Alpha-Beta to ABC (line-neutral references)
    v_a =  v_alpha_use;
    v_b = -0.5 * v_alpha_use + (sqrt(3)/2) * v_beta_use;
    v_c = -0.5 * v_alpha_use - (sqrt(3)/2) * v_beta_use;

    %% Step 3 — Zero-sequence injection (SVPWM offset)
    v_max = max([v_a, v_b, v_c]);
    v_min = min([v_a, v_b, v_c]);
    v_0   = -0.5 * (v_max + v_min);

    v_a0 = v_a + v_0;
    v_b0 = v_b + v_0;
    v_c0 = v_c + v_0;

    %% Step 4 — Duty cycle mapping
    d_a_raw = 0.5 + (v_a0 / Vdc);
    d_b_raw = 0.5 + (v_b0 / Vdc);
    d_c_raw = 0.5 + (v_c0 / Vdc);

    %% Step 6 — Clamp duties to [0, 1]
    d_a = min(max(d_a_raw, 0), 1);
    d_b = min(max(d_b_raw, 0), 1);
    d_c = min(max(d_c_raw, 0), 1);

end