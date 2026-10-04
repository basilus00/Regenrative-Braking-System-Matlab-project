function [Id_ref_mtpa, Iq_ref_mtpa] = calculate_mtpa_currents(T_ref, Ld, Lq, lambda_pm)
    % Calculate MTPA current references for given torque
    % For negative torque (regen), Iq is negative
    T_norm = T_ref / 1.5;
    
    % MTPA trajectory: Id = f(Iq)
    % Solve: T_norm = λ_pm*Iq + (Ld-Lq)*Id*Iq for Id given T and Iq
    % Using MTPA condition
    
    if abs(T_ref) < 0.01
        Id_ref_mtpa = 0;
        Iq_ref_mtpa = 0;
        return;
    end
    
    % Analytical MTPA solution
    delta_L = Ld - Lq;  % Negative for IPMSM with Lq > Ld
    
    % For IPMSM, MTPA gives:
    % Id = -(λ_pm)/(2*delta_L) - sign(Iq)*sqrt((λ_pm/(2*delta_L))^2 + Iq^2)
    
    % Solve iteratively for Iq given T
    Iq_guess = T_norm / lambda_pm;  % Initial guess (ignoring reluctance)
    
    for iter = 1:50
        Id_guess = -(lambda_pm)/(2*delta_L) - sign(Iq_guess)*sqrt((lambda_pm/(2*delta_L))^2 + Iq_guess^2);
        T_calc = lambda_pm * Iq_guess + delta_L * Id_guess * Iq_guess;
        error = T_calc - T_norm;
        
        % Numerical derivative
        dIq = 0.1 * sign(Iq_guess);
        Id_plus = -(lambda_pm)/(2*delta_L) - sign(Iq_guess+dIq)*sqrt((lambda_pm/(2*delta_L))^2 + (Iq_guess+dIq)^2);
        T_plus = lambda_pm * (Iq_guess+dIq) + delta_L * Id_plus * (Iq_guess+dIq);
        dT_dIq = (T_plus - T_calc) / dIq;
        
        if abs(dT_dIq) > 1e-10
            Iq_guess = Iq_guess - error / dT_dIq;
        else
            break;
        end
        
        if abs(error) < 1e-6
            break;
        end
    end
    
    Id_ref_mtpa = -(lambda_pm)/(2*delta_L) - sign(Iq_guess)*sqrt((lambda_pm/(2*delta_L))^2 + Iq_guess^2);
    Iq_ref_mtpa = Iq_guess;
end
