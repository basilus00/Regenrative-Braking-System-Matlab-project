function [Vd,Vq] = Current_PI_control_VoltageCalc(Id_ref,Id,Iq_ref,Iq)

    % ----- 6.5 FOC CURRENT CONTROL (PI Controllers) -----
    
    % Anti-windup limits
    Vd_max = Vs_max;
    Vq_max = Vs_max;
    
    % d-axis current error and PI
    err_Id = Id_ref - Id;
    int_err_Id = int_err_Id + err_Id * dt;
    
    % Anti-windup
    Vd_cmd = Kp_d * err_Id + Ki_d * int_err_Id;
    if abs(Vd_cmd) > Vd_max
        int_err_Id = int_err_Id - err_Id * dt;  % Revert integration
        Vd_cmd = sign(Vd_cmd) * Vd_max;
    end
    
    % q-axis current error and PI
    err_Iq = Iq_ref - Iq;
    int_err_Iq = int_err_Iq + err_Iq * dt;
    
    % Anti-windup
    Vq_cmd = Kp_q * err_Iq + Ki_q * int_err_Iq;
    if abs(Vq_cmd) > Vq_max
        int_err_Iq = int_err_Iq - err_Iq * dt;  % Revert integration
        Vq_cmd = sign(Vq_cmd) * Vq_max;
    end
    
    % Voltage limiting (circle constraint)
    Vs_cmd = sqrt(Vd_cmd^2 + Vq_cmd^2);
    if Vs_cmd > Vs_max
        Vd_cmd = Vd_cmd * Vs_max / Vs_cmd;
        Vq_cmd = Vq_cmd * Vs_max / Vs_cmd;
    end
    
    Vd = Vd_cmd;
    Vq = Vq_cmd;
    
end
