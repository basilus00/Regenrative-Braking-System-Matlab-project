%% =========================================================================
%  SIMULINK MODEL CREATION SCRIPT
%  =========================================================================

% Create new Simulink model
new_system('RegenBrakingSystem');
open_system('RegenBrakingSystem');

% Set solver configuration
set_param('RegenBrakingSystem', 'SolverType', 'Fixed-step', ...
    'Solver', 'ode4', 'FixedStep', '0.0001');

fprintf('Simulink model created. Build the following subsystems:\n');
fprintf('1. Vehicle Dynamics\n');
fprintf('2. PMSM Motor Model\n');
fprintf('3. FOC Controller (with MTPA/FW/MPTV switching)\n');
fprintf('4. Battery Model\n');
fprintf('5. Inverter Model\n');