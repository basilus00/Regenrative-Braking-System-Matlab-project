%% =========================================================================
%  1. SYSTEM PARAMETERS (INPUTS)
%  =========================================================================

% ----- Vehicle Parameters -----
m_vehicle = 2135;           % Vehicle mass [kg]
r_wheel = 0.43;             % Wheel radius [m]
gear_ratio = 11;          % Gear ratio //the second is 5:1
C_drag = 0.21;               % Aerodynamic drag coefficient
A_frontal = 2.2;            % Frontal area [m²]
Den_air = 1.225;            % Air density [kg/m³] //standard according to ISA
C_rr = 0.01;                % Rolling resistance coefficient= force to move a rolling object devided by it's load
g = 9.81;                   % Gravitational acceleration [m/s²]

% ----- Motor IPMSM Parameters -----
Tmax=330;
Lq=0.0008;
Ld=0.0004;
lambda=0.095;
p=4;

% ----- Battery Parameters -----
Vmax=407;
V_batt_nom=700;
Soc_init=60;
Soc_max=95;
I_max=580;
I_charge_max=-100;
V_batt_nom = 700;           % Nominal battery voltage [V]
Q_batt = 130;                % Battery capacity [Ah]
C_batt = Q_batt * 3600;     % Battery capacity [As]
R_batt_int = 0.05;
energy=0;

% ----- Inverter Parameters -----
Imax=580;
Vdc_nom=750;

% ----- Simulation Parameters -----
dt = 0.0001;                % Time step [s]
t_end = 10;                 % Simulation time [s]
t = 0:dt:t_end;             % Time vector
N = length(t);


% ----- Parameters initialisation-----
%vehicle
mode=0;
Vspeed=100;

%motor
Rs=0.012;
Id = 0;                    % d-axis current [A]
Iq = 0;                    % q-axis current [A]
theta_e = 0;               % Electrical angle [rad]

% Control states
Vd = 0;                    % d-axis voltage [V]
Vq = 0;                    % q-axis voltage [V]
T_motor = 0;               % Motor torque [Nm]
T_ref = 0;                 % Reference torque [Nm]
f_pwm = 10000;

%batterry
Soc=Soc_init;
V_batt=V_batt_nom;
I_batt=0;
Vs_max=V_batt/sqrt(3);
Id_char = -lambda / Lq;
Iq_char = -sqrt(I_max*I_max - Id_char*Id_char)
% Current controller bandwidth
wc = 2000;  % [rad/s]

% d-axis PI gains
Kp_d = wc * Ld;
Ki_d = wc * Rs;

% q-axis PI gains
Kp_q = wc * Lq;
Ki_q = wc * Rs;





