%% UAV Performance Analysis - Optimum Velocities and Power
%% AUTHOR = ABHIJEET MANGELA (AE25M034)

clear; clc;

%% Aircraft Parameters
Sref = 0.5307;          % Reference area (m^2)
AR = 7.53;              % Aspect ratio
b = 2;                  % Wingspan (m)
e = 0.8;                % Oswald efficiency factor
C_d_o = 0.02873;        % Zero-lift drag coefficient (given)
L_D_max = 12.55;        % Maximum lift-to-drag ratio

%% Flight Conditions
rho = 1.142;            % Air density at cruise altitude (200m)
g = 9.81;               % Gravitational acceleration

%% Weight Configurations
W1 = 6.5 * g;           % Weight before payload drop (N)
W2 = 5.3 * g;           % Weight after payload drop (N)

%% Initialize Results Storage
segments = {'Climb-1'; 'Cruise-1'; 'Loiter-1'; 'Climb-2'; 'Cruise-2'; 'Loiter-2'};
velocities = zeros(6,1);
thrusts = zeros(6,1);
powers = zeros(6,1);

%% Flight Phase Calculations - Before Payload Drop

% --- Climb-1 (gamma = 5 deg) ---
gamma_climb = 5;
W = W1;
L_climb = W * cosd(gamma_climb);
D_climb = 2 * L_climb / (sqrt(3) * L_D_max);
T_climb = D_climb + W * sind(gamma_climb);

% Find optimal velocity for climb
v0 = 20;  % Initial guess
C_L_climb_func = @(v) L_climb / (0.5 * rho * Sref * v^2);
C_D_climb_func = @(v) C_d_o + C_L_climb_func(v)^2 / (pi * e * AR);
D_climb_func = @(v) 0.5 * rho * Sref * v^2 * C_D_climb_func(v);
drag_error = @(v) D_climb_func(v) - D_climb;
V_climb1 = fsolve(drag_error, v0, optimset('Display','off'));
P_climb1 = T_climb * V_climb1;

velocities(1) = V_climb1;
thrusts(1) = T_climb;
powers(1) = P_climb1;

% --- Cruise-1 (gamma = 0 deg, level flight) ---
gamma_cruise = 0;
W = W1;
L_cruise = W * cosd(gamma_cruise);
D_cruise = L_cruise / L_D_max;
T_cruise = D_cruise + W * sind(gamma_cruise);

% Find optimal velocity for cruise
C_L_cruise_func = @(v) L_cruise / (0.5 * rho * Sref * v^2);
C_D_cruise_func = @(v) C_d_o + C_L_cruise_func(v)^2 / (pi * e * AR);
D_cruise_func = @(v) 0.5 * rho * Sref * v^2 * C_D_cruise_func(v);
drag_error = @(v) D_cruise_func(v) - D_cruise;
V_cruise1 = fsolve(drag_error, v0, optimset('Display','off'));
P_cruise1 = T_cruise * V_cruise1;

velocities(2) = V_cruise1;
thrusts(2) = T_cruise;
powers(2) = P_cruise1;

% --- Loiter-1 (gamma = 0 deg, minimum power) ---
gamma_loiter = 0;
W = W1;
L_loiter = W * cosd(gamma_loiter);
D_loiter = 2 * L_loiter / (sqrt(3) * L_D_max);
T_loiter = D_loiter + W * sind(gamma_loiter);

% Find optimal velocity for loiter
C_L_loiter_func = @(v) L_loiter / (0.5 * rho * Sref * v^2);
C_D_loiter_func = @(v) C_d_o + C_L_loiter_func(v)^2 / (pi * e * AR);
D_loiter_func = @(v) 0.5 * rho * Sref * v^2 * C_D_loiter_func(v);
drag_error = @(v) D_loiter_func(v) - D_loiter;
V_loiter1 = fsolve(drag_error, v0, optimset('Display','off'));
P_loiter1 = T_loiter * V_loiter1;

velocities(3) = V_loiter1;
thrusts(3) = T_loiter;
powers(3) = P_loiter1;

%% Flight Phase Calculations - After Payload Drop

% --- Climb-2 (gamma = 5 deg) ---
gamma_climb = 5;
W = W2;
L_climb = W * cosd(gamma_climb);
D_climb = 2 * L_climb / (sqrt(3) * L_D_max);
T_climb = D_climb + W * sind(gamma_climb);

% Find optimal velocity for climb
C_L_climb_func = @(v) L_climb / (0.5 * rho * Sref * v^2);
C_D_climb_func = @(v) C_d_o + C_L_climb_func(v)^2 / (pi * e * AR);
D_climb_func = @(v) 0.5 * rho * Sref * v^2 * C_D_climb_func(v);
drag_error = @(v) D_climb_func(v) - D_climb;
V_climb2 = fsolve(drag_error, v0, optimset('Display','off'));
P_climb2 = T_climb * V_climb2;

velocities(4) = V_climb2;
thrusts(4) = T_climb;
powers(4) = P_climb2;

% --- Cruise-2 (gamma = 0 deg, level flight) ---
gamma_cruise = 0;
W = W2;
L_cruise = W * cosd(gamma_cruise);
D_cruise = L_cruise / L_D_max;
T_cruise = D_cruise + W * sind(gamma_cruise);

% Find optimal velocity for cruise
C_L_cruise_func = @(v) L_cruise / (0.5 * rho * Sref * v^2);
C_D_cruise_func = @(v) C_d_o + C_L_cruise_func(v)^2 / (pi * e * AR);
D_cruise_func = @(v) 0.5 * rho * Sref * v^2 * C_D_cruise_func(v);
drag_error = @(v) D_cruise_func(v) - D_cruise;
V_cruise2 = fsolve(drag_error, v0, optimset('Display','off'));
P_cruise2 = T_cruise * V_cruise2;

velocities(5) = V_cruise2;
thrusts(5) = T_cruise;
powers(5) = P_cruise2;

% --- Loiter-2 (gamma = 0 deg, minimum power) ---
gamma_loiter = 0;
W = W2;
L_loiter = W * cosd(gamma_loiter);
D_loiter = 2 * L_loiter / (sqrt(3) * L_D_max);
T_loiter = D_loiter + W * sind(gamma_loiter);

% Find optimal velocity for loiter
C_L_loiter_func = @(v) L_loiter / (0.5 * rho * Sref * v^2);
C_D_loiter_func = @(v) C_d_o + C_L_loiter_func(v)^2 / (pi * e * AR);
D_loiter_func = @(v) 0.5 * rho * Sref * v^2 * C_D_loiter_func(v);
drag_error = @(v) D_loiter_func(v) - D_loiter;
V_loiter2 = fsolve(drag_error, v0, optimset('Display','off'));
P_loiter2 = T_loiter * V_loiter2;

velocities(6) = V_loiter2;
thrusts(6) = T_loiter;
powers(6) = P_loiter2;

%% Create Results Table
ResultsTable = table(segments, velocities, thrusts, powers, ...
    'VariableNames', {'Segment', 'Velocity_m_s', 'Thrust_N', 'Power_W'});

%% Display Results
fprintf('\n========================================\n');
fprintf('UAV PERFORMANCE ANALYSIS RESULTS\n');
fprintf('========================================\n\n');

fprintf('%-12s | %10s | %10s | %10s\n', 'Segment', 'Velocity', 'Thrust', 'Power');
fprintf('%-12s | %10s | %10s | %10s\n', '', '(m/s)', '(N)', '(W)');
fprintf('----------------------------------------------------------\n');

for i = 1:length(segments)
    fprintf('%-12s | %10.2f | %10.2f | %10.1f\n', ...
        segments{i}, velocities(i), thrusts(i), powers(i));
end

fprintf('----------------------------------------------------------\n\n');

% Display the table object
disp(ResultsTable);

%% Export Results to Excel
% Define output filename
output_filename = 'UAV_Performance_Results.xlsx';

% Write the main results table
writetable(ResultsTable, output_filename, 'Sheet', 'Performance Data');

% Create a summary sheet with aircraft parameters
summary_data = {
    'Parameter', 'Value', 'Unit';
    'Reference Area (Sref)', Sref, 'm^2';
    'Aspect Ratio (AR)', AR, '-';
    'Wingspan (b)', b, 'm';
    'Oswald Efficiency (e)', e, '-';
    'Zero-lift Drag Coeff (C_d_o)', C_d_o, '-';
    'Max L/D Ratio', L_D_max, '-';
    'Air Density (rho)', rho, 'kg/m^3';
    'Weight Before Payload Drop', W1, 'N';
    'Weight After Payload Drop', W2, 'N';
    'Payload Mass', (W1-W2)/g, 'kg'
};

% Write summary to second sheet
writecell(summary_data, output_filename, 'Sheet', 'Aircraft Parameters');

fprintf('\n*** Results exported to: %s ***\n\n', output_filename);
