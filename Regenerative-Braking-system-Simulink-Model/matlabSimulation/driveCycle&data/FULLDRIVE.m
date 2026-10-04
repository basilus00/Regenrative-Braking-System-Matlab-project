% ===================================================================
% EXTENDED PMSM DRIVE CYCLE (60 SECONDS)
% Generates Speed, AccelPedal, and BrakePedal signals for State Machine
% ===================================================================

clear; clc;

% 1. Configuration
Ts = 1e-3;                  % Sample time (1ms)
T_total = 60;               % Total duration (s)
t = 0:Ts:T_total;           % Time vector

% Initialize Signal Arrays
w_ref = zeros(size(t));     % Speed Reference (rad/s)
accel_pedal = zeros(size(t)); % Accelerator Pedal (0 to 1)
brake_pedal = zeros(size(t));  % Brake Pedal (0 to 1)

% Max Speed for this simulation
W_MAX = 2500; 

% =========================================================================
% DEFINE CYCLE SEGMENTS
% =========================================================================

% --- SEGMENT 1: Soft Start (0-5s) ---
% Gentle acceleration to 1000 rad/s
idx1 = t <= 5;
w_ref(idx1) = (1000/5) * t(idx1);
accel_pedal(idx1) = 0.4; % 40% throttle
brake_pedal(idx1) = 0;

% --- SEGMENT 2: Cruise (5-10s) ---
% Hold steady speed
idx2 = (t > 5) & (t <= 10);
w_ref(idx2) = 1000;
accel_pedal(idx2) = 0.1; % Maintain speed against friction
brake_pedal(idx2) = 0;

% --- SEGMENT 3: Coasting / Mild Regen (10-18s) ---
% Driver lifts foot completely. Motor coasts to 500 rad/s.
% This tests "Passive Regen" (if enabled) or just friction.
idx3 = (t > 10) & (t <= 18);
% Linear ramp down
w_ref(idx3) = 1000 - (500/8) * (t(idx3) - 10);
accel_pedal(idx3) = 0;
brake_pedal(idx3) = 0.0; % No active brake pedal, just coasting

% --- SEGMENT 4: Hard Acceleration (18-25s) ---
% Ramp up to max speed (2500)
idx4 = (t > 18) & (t <= 25);
w_ref(idx4) = 1000 + (1500/7) * (t(idx4) - 18);
accel_pedal(idx4) = 0.9; % 90% throttle
brake_pedal(idx4) = 0;

% --- SEGMENT 5: High Speed Cruise (25-30s) ---
idx5 = (t > 25) & (t <= 30);
w_ref(idx5) = 2500;
accel_pedal(idx5) = 0.15; % Higher load at high speed
brake_pedal(idx5) = 0;

% --- SEGMENT 6: Aggressive Braking (30-35s) ---
% Active braking to 0 rad/s. This is the main REGEN test.
idx6 = (t > 30) & (t <= 35);
w_ref(idx6) = 2500 - (2500/5) * (t(idx6) - 30);
accel_pedal(idx6) = 0;
brake_pedal(idx6) = 0.7; % 70% Brake Pedal

% --- SEGMENT 7: Full Stop (35-40s) ---
idx7 = (t > 35) & (t <= 40);
w_ref(idx7) = 0;
accel_pedal(idx7) = 0;
brake_pedal(idx7) = 0.1; % Holding brake

% --- SEGMENT 8: Dynamic Drive / Traffic (40-50s) ---
% Speed oscillates (Sine wave) to test tracking
idx8 = (t > 40) & (t <= 50);
base_speed = 1500;
w_ref(idx8) = base_speed + 500 * sin(2*pi*0.2*(t(idx8)-40)); 
% Simple logic: if speed is increasing, accel; if decreasing, brake
% (Simplified driver model)
slope = gradient(w_ref(idx8), Ts);
accel_pedal(idx8) = max(0, min(1, slope / 1000)); 
brake_pedal(idx8) = max(0, min(1, -slope / 1000));

% --- SEGMENT 9: Final Coast to Stop (50-60s) ---
idx9 = t > 50;
w_ref(idx9) = w_ref(find(t>=50,1)) * exp(-0.5*(t(idx9)-50)); % Exponential decay
accel_pedal(idx9) = 0;
brake_pedal(idx9) = 0.0; % Coast to stop

% =========================================================================
% PREPARE DATA FOR SIMULINK
% =========================================================================

% Create structure for 'From Workspace' blocks
DriveCycle.Time = t;
DriveCycle.Signals.Values = [w_ref', accel_pedal', brake_pedal']; % 3 columns: Speed, Accel, Brake
DriveCycle.Signals.Dimensions = 3;

% Save to workspace
assignin('base', 'DriveCycle', DriveCycle);
assignin('base', 'DriveCycle_Time', t);
assignin('base', 'DriveCycle_Speed', w_ref);
assignin('base', 'DriveCycle_Accel', accel_pedal);
assignin('base', 'DriveCycle_Brake', brake_pedal);

% =========================================================================
% VISUALIZATION
% =========================================================================
figure('Name', 'Extended Drive Cycle Profile', 'Color', 'w', 'Position', [100, 100, 1000, 600]);

% Plot 1: Speed
subplot(3,1,1);
plot(t, w_ref, 'b', 'LineWidth', 1.5);
grid on; ylabel('Speed (rad/s)'); title('Extended Drive Cycle Reference');
ylim([-500 3000]);

% Plot 2: Accel Pedal
subplot(3,1,2);
area(t, accel_pedal, 'FaceColor', [0 0.8 0]); % Green for Gas
grid on; ylabel('Accel Pedal (0-1)'); ylim([0 1.2]);

% Plot 3: Brake Pedal
subplot(3,1,3);
area(t, brake_pedal, 'FaceColor', [0.8 0 0]); % Red for Brake
grid on; ylabel('Brake Pedal (0-1)'); xlabel('Time (s)'); ylim([0 1.2]);

% Add Annotations for critical Regen points
hold on;
text(32.5, 0.8, 'AGGRESSIVE REGEN TEST', 'Color', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(14, 0.1, 'COASTING', 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

disp('Extended 60s Drive Cycle generated.');
disp('Workspace Variables created: DriveCycle (Struct), DriveCycle_Speed, DriveCycle_Accel, DriveCycle_Brake');