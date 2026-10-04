%% 1. Open and Prepare the Base Model
example_name = 'PMSMDriveExample';
openExample(sprintf('simscapeelectrical/%s', example_name));

% Save as a new model so we don't overwrite the example
newModelName = 'PMSM_Regen_StateMachine';
save_system(newModelName);
close_system(example_name); % Close the original example
open_system(newModelName);

% Get the controller subsystem path (Standard PMSM Drive structure)
ctrlSubsys = [newModelName '/Control'];

% Open the controller to modify it
open_system(ctrlSubsys);

%% 2. Switch Control Mode to Torque Control
% The standard block usually has a 'Control Mode' parameter.
% For regen, we need to command Torque directly, not Speed.
try
    % Attempt to find the PMSM Controller block
    ctrlBlock = find_system(ctrlSubsys, 'BlockType', 'SubSystem', 'Name', '*/PMSM Controller*');
    if ~isempty(ctrlBlock)
        % Set to Torque Control (Typically parameter value 'Torque Control')
        % Note: Exact parameter name depends on version, often 'ControlMode'
        set_param(ctrlBlock{1}, 'ControlMode', 'Torque Control');
    end
catch
    warning('Could not auto-set control mode. Please ensure PMSM Controller is set to Torque Control manually.');
end

%% 3. Create the State Machine (Stateflow)
% We will build a Stateflow chart programmatically

% Remove existing Speed Reference logic if it exists (Simplification)
refBlock = [ctrlSubsys '/Speed Reference'];
if exist(refBlock, 'block')
    delete_block(refBlock);
end

% Add Stateflow Chart
chartName = 'DriverLogic_StateMachine';
sfBlock = add_block('stateflow/Chart', [ctrlSubsys '/' chartName], ...
    'Position', [50 200 300 200]);

% Define Stateflow API
machine = sfroot;
chart = machine.find('-isa', 'Stateflow.Chart', 'Name', chartName);

% Define Events and Data
% Inputs
data_accel = Stateflow.Data(chart);
data_accel.Name = 'AccelPedal';
data_accel.Scope = 'Input';
data_accel.DataType = 'double';

data_brake = Stateflow.Data(chart);
data_brake.Name = 'BrakePedal';
data_brake.Scope = 'Input';
data_brake.DataType = 'double';

data_regen_en = Stateflow.Data(chart);
data_regen_en.Name = 'RegenEnable';
data_regen_en.Scope = 'Input';
data_regen_en.DataType = 'boolean';

% Outputs
data_tref = Stateflow.Data(chart);
data_tref.Name = 'TorqueRef';
data_tref.Scope = 'Output';
data_tref.DataType = 'double';

% Internal Constants
data_max_t = Stateflow.Data(chart);
data_max_t.Name = 'MaxTorque';
data_max_t.Scope = 'Constant';
data_max_t.DataType = 'double'; % Set to appropriate max torque for your motor

% Define States
% We create a textual representation of the logic for simplicity
chart.Text = sprintf([...
    'state machine:\n' ...
    '  state IDLE:\n' ...
    '    during: TorqueRef = 0;\n' ...
    '    transition: AccelPedal > 0 -> MOTORING;\n' ...
    '\n' ...
    '  state MOTORING:\n' ...
    '    during: TorqueRef = AccelPedal * MaxTorque;\n' ...
    '    transition: AccelPedal == 0 && BrakePedal > 0 -> BRAKING_DECISION;\n' ...
    '\n' ...
    '  state BRAKING_DECISION:\n' ...
    '    state REGEN:\n' ...
    '      during: TorqueRef = -BrakePedal * MaxTorque;\n' ...
    '      transition: ~RegenEnable -> COASTING;\n' ...
    '    state COASTING:\n' ...
    '      during: TorqueRef = 0;\n' ...
    '      transition: RegenEnable -> REGEN;\n' ...
    '    transition: BrakePedal == 0 -> IDLE;\n' ...
    '    default: REGEN\n' ... % Default to Regen if enabled at start of brake
    ]);

%% 4. Create Simulation Inputs (Driver & Cycle)
% Add signals to the top level
sys = get_param(newModelName, 'Handle');

% 1. Braking Cycle Signal (For Driver Model)
add_block('simulink/Sources/From Workspace', [newModelName '/Ref_Signal'], ...
    'VariableName', 'BrakingCycle_Data', 'Position', [50 50 150 70]);

% 2. Regen Enable Switch
add_block('simulink/Sources/Constant', [newModelName '/Regen_Switch'], ...
    'Value', '1', 'Position', [50 120 150 140]); % 1 = Enabled, 0 = Disabled

% 3. Driver Model (Simplified: Compares Ref Speed to Actual to generate Pedals)
% Ideally, we use a Driver block. For this demo, we will approximate the logic
% that the "Driver" block would produce using a MATLAB Function block or logic.
% Here, we will route the Reference Signal directly to a modified Controller path
% and let the User control the "Brake" manually for testing, or use a logic block.

% NOTE: The PMSMDriveExample uses a "Driver" block that outputs Throttle/Brake.
% We will bypass the complex Driver block for this specific Regen test 
% and inject signals manually to prove the State Machine works.

% Delete existing Driver block if it exists
driverBlk = find_system(newModelName, 'Name', 'Driver');
if ~isempty(driverBlk)
    delete_block(driverBlk{1});
end

% Create a Signal Builder or Manual Switch for user interaction testing
% For automation, let's create a Signal that simulates Acceleration then Braking
t = 0:0.01:8;
accel_sig = [zeros(1,300) ones(1,100) zeros(1,400)]; % Accel then coast/brake
brake_sig = [zeros(1,400) ones(1,300) zeros(1,100)]; % Brake later
assignin('base', 'test_accel', accel_sig);
assignin('base', 'test_brake', brake_sig);

blk_accel = add_block('simulink/Sources/From Workspace', [newModelName '/Accel_In'], ...
    'VariableName', 'test_accel', 'Position', [300 50 420 70]);
blk_brake = add_block('simulink/Sources/From Workspace', [newModelName '/Brake_In'], ...
    'VariableName', 'test_brake', 'Position', [300 100 420 120]);

%% 5. Connect the State Machine
% The State Machine is inside the Control Subsystem. We need to route inputs to it.

% Route inputs to the Control Subsystem Inputs
% The PMSM Controller block usually takes 'T_ref' or 'w_ref'.
% Since we set it to Torque Control, look for the Torque Reference input port.

% Find the input port of the Controller subsystem
% Note: You may need to manually connect these if port names differ by version
% We will add Outports to the top level to expose connections

% Connect Regen Switch
add_line(newModelName, 'Regen_Switch/1', 'Control/4', 'autorouting', 'on');

% We need to modify the Control Subsystem mask or inputs to accept Accel/Brake/Regen.
% This is complex via script. 
% ALTERNATIVE: Put the State Machine at the TOP LEVEL and feed TorqueRef INTO the controller.

% Moving State Machine to Top Level for easier integration
move_block(sfBlock.Path, [newModelName '/DriverLogic_StateMachine']);
set_param([newModelName '/DriverLogic_StateMachine'], 'Position', [500 50 700 200]);

% Connect inputs to State Machine
add_line(newModelName, 'Accel_In/1', 'DriverLogic_StateMachine/1');
add_line(newModelName, 'Brake_In/1', 'DriverLogic_StateMachine/2');
add_line(newModelName, 'Regen_Switch/1', 'DriverLogic_StateMachine/3');

% Now we need to feed the TorqueRef output into the Controller
% Identify the Controller block inside the subsystem
pmsmCtrlBlk = find_system(ctrlSubsys, 'BlockType', 'SubSystem', 'Name', '*/PMSM Controller*');
if ~isempty(pmsmCtrlBlk)
    ctrlInputPort = '1'; % Usually Port 1 is Reference
    % Add an input port to the Control Subsystem for External Torque
    add_block('simulink/Ports & Subsystems/In1', [ctrlSubsys '/T_ext_in'], ...
        'Position', [20 250 40 270]);
    
    % Connect Top Level State Machine -> Control Subsystem Input
    add_line(newModelName, 'DriverLogic_StateMachine/1', 'Control/5'); % Assuming 5th port is new
    
    % Connect Control Subsystem Input -> PMSM Controller Block
    % (This step often requires manual hookup inside the subsystem)
    disp('MANUAL STEP REQUIRED: Open the Control subsystem.');
    disp(['Connect the new input port T_ext_in to the Torque Reference input of ' pmsmCtrlBlk{1}]);
else
    disp('Could not find PMSM Controller block automatically.');
end

disp('Model Construction Complete.');