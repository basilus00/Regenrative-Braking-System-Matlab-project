function [Id_ref, Iq_ref] = calculate_mptv_currents(T_ref, Vs_max, omega_e, Rs, Ld, Lq, lambda_pm)
    % Maximum Torque Per Voltage - voltage limited operation
    
    delta_L = Ld - Lq;
    T_norm = T_ref / 1.5;
    
    % MPTV occurs at the intersection of constant torque curve and voltage ellipse
    % with maximum Id (most negative) for that torque
    
    % Start from characteristic current
    Id_char = -lambda_pm / Lq;
    Id_guess = Id_char;
    Iq_guess = T_norm / (lambda_pm + delta_L * Id_char);
    
    for iter = 1:100
        % Voltage
        Vd = Rs*Id_guess - omega_e*Lq*Iq_guess;
        Vq = Rs*Iq_guess + omega_e*Ld*Id_guess + omega_e*lambda_pm;
        Vs = sqrt(Vd^2 + Vq^2);
        
        % Torque
        T_calc = lambda_pm * Iq_guess + delta_L * Id_guess * Iq_guess;
        
        err_V = Vs - Vs_max;
        err_T = T_calc - T_norm;
        
        % Adjust
        dId = 0.5;
        Vd_plus = Rs*(Id_guess+dId) - omega_e*Lq*Iq_guess;
        Vq_plus = Rs*Iq_guess + omega_e*Ld*(Id_guess+dId) + omega_e*lambda_pm;
        Vs_plus = sqrt(Vd_plus^2 + Vq_plus^2);
        dVs_dId = (Vs_plus - Vs) / dId;
        
        if abs(dVs_dId) > 1e-10
            Id_guess = Id_guess - err_V / dVs_dId;
        end
        
        if abs(lambda_pm + delta_L * Id_guess) > 1e-10
            Iq_guess = T_norm / (lambda_pm + delta_L * Id_guess);
        end
        
        if abs(err_V) < 0.1 && abs(err_T) < 0.01
            break;
        end
    end
    
    Id_ref = Id_guess;
    Iq_ref = Iq_guess;
end
