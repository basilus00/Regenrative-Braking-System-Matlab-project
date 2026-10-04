

%% Battery Update Function (standalone)
function [P_batt, P_loss] = ...
          update_battery(Vd, Vq, Id, Iq, SOC, V_batt, energy)

    % All parameters come from base workspace automatically
    % because this runs in the same script context
    
    global Rs V_batt_nom R_batt_int C_batt I_charge_max SOC_max dt

    %% 6.10
    P_motor  = 1.5 * (Vd*Id + Vq*Iq);
    P_switch = 0.02 * abs(P_motor);
    P_cond   = 1.5 * Rs * (Id^2 + Iq^2);
    P_loss   = P_switch + P_cond;
    P_batt   = P_motor - P_loss;
    I_batt   = P_batt / V_batt;

    if I_batt < I_charge_max
        I_batt = I_charge_max;
        P_batt = I_batt * V_batt;
    end

    V_batt = V_batt_nom - I_batt * R_batt_int;

    %% 6.11
    SOC = SOC + (-I_batt * dt / C_batt);
    SOC = min(max(SOC, 0), SOC_max);

    %% 6.12
    if P_batt < 0
        energy = energy + abs(P_batt) * dt;
    end
end