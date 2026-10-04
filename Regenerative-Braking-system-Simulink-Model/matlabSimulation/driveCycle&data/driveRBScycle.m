% ===================================================================
% EV VEHICLE SPEED DRIVE CYCLE GENERATOR (m/s)
% ===================================================================

% 1. Time Settings
T_stop = 10;            % Total simulation time (seconds)
Ts = 1e-4;              % Sample time (seconds)
t = 0:Ts:T_stop;

% 2. Initialize Speed Vector (Vehicle Speed in m/s)
Ref_Vehicle_Speed = zeros(size(t));

% 3. Define Cycle Phases
% -----------------------------------------------------------
% PHASE 1: HARD ACCELERATION (Building Kinetic Energy)
% We ramp up to 40 m/s (~144 km/h) quickly.
% -----------------------------------------------------------
t1 = 2.0;               % Duration of acceleration
v_target = 200;          % Target speed in m/s (Realistic highway speed)

idx_accel = t <= t1;
Ref_Vehicle_Speed(idx_accel) = (v_target / t1) * t(idx_accel);

% -----------------------------------------------------------
% PHASE 2: COASTING / STEADY STATE
% Brief period to verify stability before braking.
% -----------------------------------------------------------
t_coast_start = t1;
t_coast_end = 3.0;      % Start braking at 3 seconds
idx_coast = (t > t_coast_start) & (t <= t_coast_end);
Ref_Vehicle_Speed(idx_coast) = v_target;

% -----------------------------------------------------------
% PHASE 3: AGGRESSIVE BRAKING (Regen Test)
% Ramp down from 40 m/s to 0 m/s.
% The "Slope" here represents the braking deceleration required.
% -----------------------------------------------------------
t_brake_end = 6.0;      % Stop by 6 seconds
idx_brake = (t > t_coast_end) & (t <= t_brake_end);

% Calculate linear ramp down: v = v_start - slope * (t - t_start)
brake_slope = v_target / (t_brake_end - t_coast_end);
Ref_Vehicle_Speed(idx_brake) = v_target - brake_slope * (t(idx_brake) - t_coast_end);

% -----------------------------------------------------------
% PHASE 4: STOPPED
% -----------------------------------------------------------
idx_stop = t > t_brake_end;
Ref_Vehicle_Speed(idx_stop) = 0;

% 4. Visualization
figure('Name', 'Vehicle Braking Drive Cycle', 'Color', 'w');
plot(t, Ref_Vehicle_Speed, 'b', 'LineWidth', 2.5);
grid on;
hold on;

% Annotations
xline(t_coast_end, '--r', 'Brake Start', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
yline(v_target, '--k', 'Max Speed', 'LineWidth', 1.5);

title('EV Drive Cycle: Acceleration & Regenerative Braking');
xlabel('Time (s)'); ylabel('Vehicle Speed (m/s)');
legend('Reference Speed', 'Location', 'NorthEast');

% 5. Prepare Data for Simulink
% We format this as a structure compatible with Simulink 'From Workspace' blocks
VehicleCycle.Time = t;
VehicleCycle.Signals.Values = Ref_Vehicle_Speed;
VehicleCycle.Signals.Dimensions = 1;

% Save to Base Workspace
assignin('base', 'VehicleCycle', VehicleCycle);

disp('==============================================');
disp('Vehicle Drive Cycle Generated Successfully.');
disp(['Variable Name: ''VehicleCycle'' saved to workspace.']);
disp(['Max Speed: ' num2str(v_target) ' m/s (' num2str(v_target*3.6) ' km/h)']);
disp(['Braking Duration: ' num2str(t_brake_end - t_coast_end) ' s']);
disp('==============================================');

