clear;
clc;
close all;


%% Hardcoded Values
MTOW        = 6.9;         % kg
g           = 9.81;        % m/s^2
chord       = 0.265;       % m, wing MAC
b_wing      = 2.0;         % m, wing span
sigma_yield = 150e6;       % Pa, Al 1100
sigma_yield_gear = 276e6;  % Pa, Al 6061
E_gear      = 68.9e9;      % Pa, Al 6061 Young's modulus
n_gear      = 3.0;         % design landing load factor
t_leg       = 1e-3;        % m, Thickness
FOS         = 1.5;         % factor of safety
n_fatigue   = 1.5;         % fatigue factor
k_stress    = 1.5;         % stress concentration factor

CHANNEL_ASPECT_RATIO = 1; 
B_SEARCH_LO = 2.5*t_leg;   
B_SEARCH_HI = 0.10;        
BISECTION_TOL = 1e-9;      
BISECTION_MAX_ITER = 100;  
tipback_allow_deg   = 15;  % deg, standard GA minimum (Raymer)
turnover_allow_deg  = 63;  % deg, standard GA maximum (Raymer)
allow_defl_deg      = 5;   % deg, adopted allowable tip deflection angle


%% GEOMETRY

% Overall CG
X_cg            = 275.57e-3;  % m
Y_cg            = 13.03e-3;   % m
Z_cg            = 0.53e-3;    % m

% Nose landing gear
X_nose_attach   = 67.5e-3;    % m, fuselage attach point
Y_nose_attach   = -70.02e-3;  % m
Z_nose_attach   = 0;          % m

X_nose_wheel    = 56.85e-3;   % m, wheel centre
Y_nose_wheel    = -183.02e-3; % m
Z_nose_wheel    = 0;          % m

% Main landing gear
X_main_attach   = 369.31e-3;  % m, fuselage attach point
Y_main_attach   = -91.90e-3;  % m
Z_main_attach   = 0;          % m

X_main_wheel    = 369.31e-3;  % m, wheel centre
Y_main_wheel    = -183.02e-3; % m
Z_main_wheel    = 156.35e-3;  % m

% Wheels Specifications
d_wheel_nose    = 80e-3;      % m, nose wheel diameter
d_wheel_main    = 80e-3;      % m, main wheel diameter
m_wheel         = 170e-3;     % kg, mass per wheel

RAYMER_A_D      = 5.1;        % diameter constant
RAYMER_B_D      = 0.349;      % diameter exponent
RAYMER_A_B      = 2.3;        % width constant
RAYMER_B_B      = 0.312;      % width exponent

PLOT_N_PTS      = 100;       

theta_main_deg  = 35;         % deg, strut angle from vertical
MAX_BEAM_WIDTH  = 200e-3;     % m


%% ---------------- 3. DERIVED GEOMETRY -----------------------------------

pass_fail = {'FAIL', 'PASS'};

r_wheel_nose = d_wheel_nose/2;
r_wheel_main = d_wheel_main/2;

Y_ground = Y_main_wheel - r_wheel_main;

% Basic gear parameters
wheelbase = X_main_wheel - X_nose_wheel;                 % m
track     = 2*Z_main_wheel;                              % m, main gear track
h_cg      = Y_cg - Y_ground;                             % m, CG height above ground

leg_len_main = sqrt((X_main_wheel-X_main_attach)^2 + (Y_main_wheel-Y_main_attach)^2 ...
                     + (Z_main_wheel-Z_main_attach)^2);
leg_len_nose = sqrt((X_nose_wheel-X_nose_attach)^2 + (Y_nose_wheel-Y_nose_attach)^2 ...
                     + (Z_nose_wheel-Z_nose_attach)^2);

rake_main_frontview_deg = atand( (Z_main_wheel-Z_main_attach) / (Y_main_attach-Y_main_wheel) );

rake_nose_sideview_deg = atand( (X_nose_attach-X_nose_wheel) / (Y_nose_attach-Y_nose_wheel) );

height_main = Y_main_attach - Y_ground;   % m, vertical drop, attach -> ground
height_nose = Y_nose_attach - Y_ground;   % m, vertical drop, attach -> ground

fprintf('\n--- Basic gear parameters ---\n');
fprintf('Wheelbase            : %.2f mm\n', wheelbase*1e3);
fprintf('Main gear track      : %.2f mm\n', track*1e3);
fprintf('CG height above ground: %.2f mm\n', h_cg*1e3);
fprintf('Main leg length       : %.2f mm, front-view rake %.2f deg from vertical\n', ...
        leg_len_main*1e3, rake_main_frontview_deg);
fprintf('Nose leg length       : %.2f mm, side-view rake %.2f deg from vertical\n', ...
        leg_len_nose*1e3, rake_nose_sideview_deg);
fprintf('Main gear height (attach->ground): %.2f mm\n', height_main*1e3);
fprintf('Nose gear height (attach->ground): %.2f mm\n', height_nose*1e3);


%% ---------------- 4. STATIC LOAD DISTRIBUTION ---------------------------

W = MTOW*g;   % N, total static weight

W_main_total = W * (X_cg - X_nose_wheel) / wheelbase;   % both main legs combined
W_nose       = W * (X_main_wheel - X_cg) / wheelbase;   % nose gear

W_main_each  = W_main_total/2;                          % per main leg

fprintf('Total weight W        : %.2f N\n', W);
fprintf('Nose gear reaction     : %.2f N  (%.1f%% of W)\n', W_nose, 100*W_nose/W);
fprintf('Main gear reaction (tot): %.2f N  (%.1f%% of W)\n', W_main_total, 100*W_main_total/W);
fprintf('Main gear reaction (each leg): %.2f N\n', W_main_each);


%% ---------------- 5.1 TIRE SIZING -----------

W_nose_kg      = W_nose / g;          % kg, equivalent mass on nose wheel
W_main_each_kg = W_main_each / g;     % kg, equivalent mass on each main wheel

D_main_cm = RAYMER_A_D * W_main_each_kg^RAYMER_B_D;
b_main_cm = RAYMER_A_B * W_main_each_kg^RAYMER_B_B;

D_nose_cm = RAYMER_A_D * W_nose_kg^RAYMER_B_D;
b_nose_cm = RAYMER_A_B * W_nose_kg^RAYMER_B_B;

fprintf('Main wheel load (each): %.3f kg -> D_m = %.2f cm, b_m = %.2f cm\n', ...
        W_main_each_kg, D_main_cm, b_main_cm);
fprintf('Nose wheel load       : %.3f kg -> D_n = %.2f cm, b_n = %.2f cm\n', ...
        W_nose_kg, D_nose_cm, b_nose_cm);
fprintf('\nFor comparison, catalog wheels are %.1f mm dia (both) - the\n', d_wheel_main*1e3);
fprintf('calculated values above are what Raymer''s formula recommends as a\n');
fprintf('MINIMUM based on load alone; catalog wheels exceeding this is fine.\n');


%% ---------------- 5. DESIGN (IMPACT) LOADS ------------------------------

P_nose_design      = n_gear * W_nose;
P_main_design_each  = n_gear * W_main_each;

fprintf('\n--- Design (impact) loads, n_gear = %.2f ---\n', n_gear);
fprintf('Nose gear design load       : %.2f N\n', P_nose_design);
fprintf('Main gear design load (each): %.2f N\n', P_main_design_each);


%% ---------------- 6. STABILITY CHECKS (tip-back & turnover) -------------

tipback_deg = atand( (X_main_wheel - X_cg) / h_cg );

N_pt = [X_nose_wheel, 0];              % nose wheel, plan view
M_pt = [X_main_wheel, Z_main_wheel];   % one main wheel, plan view
C_pt = [X_cg, Z_cg];

line_vec = M_pt - N_pt;
NC_vec   = C_pt - N_pt;
d_perp   = abs(line_vec(1)*NC_vec(2) - line_vec(2)*NC_vec(1)) / norm(line_vec);

turnover_deg = atand(d_perp / h_cg);

fprintf('Tip-back angle  : %.2f deg  (allowable >= %.1f deg) -> %s\n', ...
        tipback_deg, tipback_allow_deg, pass_fail{(tipback_deg >= tipback_allow_deg)+1});
fprintf('Turnover angle  : %.2f deg  (allowable <= %.1f deg) -> %s\n', ...
        turnover_deg, turnover_allow_deg, pass_fail{(turnover_deg <= turnover_allow_deg)+1});

%% ---------------- 5.2 MAIN GEAR DESIGN - BMD -----

HLG_strut_main = Y_main_attach - Y_main_wheel;   % m
theta_main_rad = deg2rad(theta_main_deg);

L_strut_main     = HLG_strut_main / cos(theta_main_rad);      % m, inclined strut length
horiz_proj_main  = L_strut_main * sin(theta_main_rad);        % m, horizontal projection
beam_halfwidth_main = (track/2) - horiz_proj_main;              
sect2_len_main   = beam_halfwidth_main;                         

fprintf('Strut angle from vertical    : %.0f deg (design choice)\n', theta_main_deg);
fprintf('Section 1 (inclined strut)   : %.2f mm long\n', L_strut_main*1e3);
fprintf('Horizontal projection of strut: %.2f mm\n', horiz_proj_main*1e3);
fprintf('Section 2 (horizontal beam)  : %.2f mm HALF-width -> %.2f mm full beam width\n', ...
        beam_halfwidth_main*1e3, 2*beam_halfwidth_main*1e3);
fprintf('Fuselage width constraint     : %.2f mm max -> %s\n', ...
        MAX_BEAM_WIDTH*1e3, pass_fail{(2*beam_halfwidth_main <= MAX_BEAM_WIDTH)+1});

wt_main = b_main_cm / 100;   % m


F_normal = W_main_each;      % N, static reaction


F_worst = W_main_total;      % N, static reaction

M0_normal = F_normal * wt_main/2;
M0_worst  = F_worst  * wt_main/2;


fprintf('Tyre width : %.2f mm\n', wt_main*1e3);
fprintf('Normal case: F = %.2f N (per wheel), M0 = %.3f Nm\n', F_normal, M0_normal);
fprintf('Worst case : F = %.2f N (single wheel), M0 = %.3f Nm\n', F_worst, M0_worst);

x1 = linspace(0, L_strut_main, PLOT_N_PTS);             
x2 = linspace(0, sect2_len_main, PLOT_N_PTS);            

M1_normal = M0_normal + x1*F_normal*sin(theta_main_rad);
M2_normal = M0_normal + F_normal*(x2 + L_strut_main*sin(theta_main_rad));

M1_worst = M0_worst + x1*F_worst*sin(theta_main_rad);
M2_worst = M0_worst + F_worst*(x2 + L_strut_main*sin(theta_main_rad));

s1 = x1;                          
s2 = L_strut_main + x2;           

M_max_normal = max([M1_normal, M2_normal]);
M_max_worst  = max([M1_worst,  M2_worst]);

fprintf('M_max (normal case)  : %.3f Nm\n', M_max_normal);
fprintf('M_max (worst case)   : %.3f Nm  <-- use this for cross-section sizing\n', M_max_worst);
fprintf('(Note: M_max = M0 + F*(track/2) is independent of strut angle -\n');
fprintf(' the angle only reshapes the M(x) curve and affects axial/buckling.)\n');

try
    fig1 = figure('visible','off');
    plot([s1, s2]*1e3, [M1_normal, M2_normal], 'r-', 'LineWidth', 2);
    xlabel('Distance along gear contour (mm)');
    ylabel('Bending Moment (N\cdotm)');
    title('Main gear BMD - normal case (load split both wheels)');
    grid on;
    set(fig1, 'renderer', 'painters');
    exportgraphics(fig1, 'main_gear_BMD_normal.png', 'Resolution', 150);
    close(fig1);

    fig2 = figure('visible','off');
    plot([s1, s2]*1e3, [M1_worst, M2_worst], 'r-', 'LineWidth', 2);
    xlabel('Distance along gear contour (mm)');
    ylabel('Bending Moment (N\cdotm)');
    title('Main gear BMD - worst case (single wheel carries full load)');
    grid on;
    set(fig2, 'renderer', 'painters');
    exportgraphics(fig2, 'main_gear_BMD_worst.png', 'Resolution', 150);
    close(fig2);
catch ME
    fprintf('\n[Plot skipped - no graphics backend available: %s]\n', ME.message);
end


%% ---------------- 5.2c CROSS-SECTION SIZING --


sigma_allowable_gear = sigma_yield_gear / (n_fatigue * FOS * k_stress);
S_required = M_max_worst / sigma_allowable_gear;

fprintf('Allowable stress       : %.2f MPa\n', sigma_allowable_gear/1e6);
fprintf('Required section modulus S = M_max/sigma_allow: %.4e m^3\n', S_required);

b_lo = B_SEARCH_LO; b_hi = B_SEARCH_HI;
for iter = 1:BISECTION_MAX_ITER
    b_mid = (b_lo + b_hi)/2;
    h_mid = CHANNEL_ASPECT_RATIO * b_mid;
    [S_mid, ~, ~] = channel_section_props(b_mid, h_mid, t_leg);
    if S_mid < S_required
        b_lo = b_mid;
    else
        b_hi = b_mid;
    end
    if (b_hi - b_lo) < BISECTION_TOL
        break
    end
end
b_req = b_hi;                         
h_req = CHANNEL_ASPECT_RATIO * b_req;
[S_req_actual, y_max_req, I_req_actual] = channel_section_props(b_req, h_req, t_leg);

fprintf('Converged after %d bisection iterations\n', iter);
fprintf('Required b (flange width): %.2f mm\n', b_req*1e3);
fprintf('Required h (web height)  : %.2f mm  (h = %.1f * b)\n', h_req*1e3, CHANNEL_ASPECT_RATIO);
fprintf('Thickness t               : %.2f mm\n', t_leg*1e3);
fprintf('Resulting I                : %.4e m^4\n', I_req_actual);
fprintf('Resulting y_max             : %.4f mm\n', y_max_req*1e3);
fprintf('Resulting section modulus S: %.4e m^3  (required %.4e) -> %s\n', ...
        S_req_actual, S_required, pass_fail{(S_req_actual >= S_required)+1});


%% ---------------- 5.2d TIP DEFLECTION ------------------------------------

EI = E_gear * I_req_actual;

a1 = M0_worst;
b1 = F_worst*sin(theta_main_rad);
a2 = M0_worst + F_worst*L_strut_main*sin(theta_main_rad);
b2 = F_worst;

int_M_sec1  = a1*L_strut_main + b1*L_strut_main^2/2;
int_M_sec2  = a2*sect2_len_main + b2*sect2_len_main^2/2;
int_M2_sec1 = a1^2*L_strut_main + a1*b1*L_strut_main^2 + b1^2*L_strut_main^3/3;
int_M2_sec2 = a2^2*sect2_len_main + a2*b2*sect2_len_main^2 + b2^2*sect2_len_main^3/3;

int_M  = int_M_sec1 + int_M_sec2;
int_M2 = int_M2_sec1 + int_M2_sec2;

theta_tip_rad = int_M / EI;
theta_tip_deg = rad2deg(theta_tip_rad);
delta_tip     = int_M2 / (F_worst * EI);

fprintf('Tip deflection            : %.3f mm\n', delta_tip*1e3);
fprintf('Tip angular deflection    : %.3f deg  (allowable <= %.1f deg) -> %s\n', ...
        theta_tip_deg, allow_defl_deg, pass_fail{(theta_tip_deg <= allow_defl_deg)+1});


%% ---------------- LOCAL FUNCTIONS ----------------------------------------
function [S, y_max, I] = channel_section_props(b, h, t)

    A_web    = (h - 2*t) * t;
    A_flange = b * t;                   
    A_total  = A_web + 2*A_flange;

    x_web    = t/2;                       
    x_flange = b/2;                       

    x_bar = (A_web*x_web + 2*A_flange*x_flange) / A_total;   
    I_web_own    = (h - 2*t) * t^3 / 12;               
    I_flange_own = t * b^3 / 12;                       
    I = I_web_own + A_web*(x_bar - x_web)^2 ...
        + 2*(I_flange_own + A_flange*(x_bar - x_flange)^2);

    y_max = max(x_bar, b - x_bar);        
    S = I / y_max;                         
end