%% (L/D)max estimation of similar UAVs

%% AUTHOR = ABHIJEET MANGELA (AE25M034)

%ISA Parameters
h = hcruise;
rho = NaN(size(h));    
T0 = 288.15;             
P0 = 101325;              
rho0 = 1.225;             
L = 0.0065;               
R = 287.05;               
g = 9.80665;              



%------------------------------------------------------------------------
%for our UAV
Sref=0.5307
our_L_D_max = 12.55


% cruise
AR=7.53
MTOW=6.5*9.81; %
e=0.8 %for rectangular wing
rho=1.142; % @ cruising altitude = 200m
b=2;
s=0.5307
vcr=22;
Cl=(2*MTOW)/(rho*vcr^2*s);
k=1/(pi*AR*e);
v0=30; 
L_D_max = our_L_D_max
gamma = 8;
Lift_cruise = MTOW;
Drag_cruise = Lift_cruise/L_D_max;
thrust_required_cruise = Drag_cruise;
C_L_cruise = @(v) Lift_cruise/(0.5*rho*s*v^2);
C_d_i_cruise_func = @(v) C_L_cruise(v)^2/(pi*e*AR);
Drag_function_cruise = @(v)(0.5*rho*s*v^2*(2*C_d_i_cruise_func(v))) - Drag_cruise;
optimal_velocity_cruise = fsolve(Drag_function_cruise,v0);
Power_required_cruise = thrust_required_cruise*optimal_velocity_cruise;
T_W_cruise=thrust_required_cruise/MTOW;
P_W_cruse = Power_required_cruise/MTOW;
CDo = C_d_i_cruise_func(optimal_velocity_cruise);

% climb
gamma = 5;
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_climb = Drag + MTOW*sind(gamma);
C_L_climb = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_climb_func = @(v) C_L_climb(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_climb_func(v))) - Drag;
optimal_velocity_climb = fsolve(Drag_function,v0);
Power_required_climb = thrust_required_climb*optimal_velocity_climb
T_W_climb=thrust_required_climb/MTOW;
P_W_climb = Power_required_climb/MTOW;

% descent
gamma =-8;
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_descent = Drag + MTOW*sind(gamma);
C_L_descent = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_descent_func = @(v) C_L_descent(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_descent_func(v))) - Drag;
optimal_velocity_descent = fsolve(Drag_function,v0);
Power_required_descent = thrust_required_descent*optimal_velocity_descent;
T_W_descent=thrust_required_descent/MTOW;
P_W_climb = Power_required_descent/MTOW;

% loiter

gamma = 0
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_loiter = Drag + MTOW*sind(gamma);
C_L_loiter = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_loiter_func = @(v) C_L_loiter(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_loiter_func(v))) - Drag;
optimal_velocity_loiter= fsolve(Drag_function,v0);
Power_required_loiter = thrust_required_loiter*optimal_velocity_loiter;
T_W_loiter=thrust_required_loiter/MTOW;
P_W_loiter = Power_required_loiter/MTOW;

%MTOW AFTER PAYLOAD Drop

% cruise 2
AR=7.53
MTOW=5.3*9.81; %MTOW AFTER PAYLOAD Drop
e=0.8 %for rectangular wing
rho=1.142; % @ cruising altitude = 200m
b=2;
s=0.5307
vcr=22;
Cl=(2*MTOW)/(rho*vcr^2*s);
k=1/(pi*AR*e);
v0=30; 
L_D_max = our_L_D_max
gamma = 8;
Lift_cruise = MTOW;
Drag_cruise = Lift_cruise/L_D_max;
thrust_required_cruise = Drag_cruise;
C_L_cruise = @(v) Lift_cruise/(0.5*rho*s*v^2);
C_d_i_cruise_func = @(v) C_L_cruise(v)^2/(pi*e*AR);
Drag_function_cruise = @(v)(0.5*rho*s*v^2*(2*C_d_i_cruise_func(v))) - Drag_cruise;
optimal_velocity_cruise = fsolve(Drag_function_cruise,v0);
Power_required_cruise = thrust_required_cruise*optimal_velocity_cruise;
T_W_cruise=thrust_required_cruise/MTOW;
P_W_cruse = Power_required_cruise/MTOW;
CDo = C_d_i_cruise_func(optimal_velocity_cruise);

% climb 2
gamma = 5;
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_climb = Drag + MTOW*sind(gamma);
C_L_climb = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_climb_func = @(v) C_L_climb(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_climb_func(v))) - Drag;
optimal_velocity_climb = fsolve(Drag_function,v0);
Power_required_climb = thrust_required_climb*optimal_velocity_climb
T_W_climb=thrust_required_climb/MTOW;
P_W_climb = Power_required_climb/MTOW;

% descent 2
gamma =-8;
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_descent = Drag + MTOW*sind(gamma);
C_L_descent = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_descent_func = @(v) C_L_descent(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_descent_func(v))) - Drag;
optimal_velocity_descent = fsolve(Drag_function,v0);
Power_required_descent = thrust_required_descent*optimal_velocity_descent;
T_W_descent=thrust_required_descent/MTOW;
P_W_climb = Power_required_descent/MTOW;

% loiter 2

gamma = 0
Lift = MTOW*cosd(gamma);
Drag = 2*Lift/(sqrt(3)*L_D_max);
thrust_required_loiter = Drag + MTOW*sind(gamma);
C_L_loiter = @(v) Lift/(0.5*rho*s*v^2);
C_d_i_loiter_func = @(v) C_L_loiter(v)^2/(pi*e*AR);
Drag_function = @(v)(0.5*rho*s*v^2*(CDo + C_d_i_loiter_func(v))) - Drag;
optimal_velocity_loiter= fsolve(Drag_function,v0);
Power_required_loiter = thrust_required_loiter*optimal_velocity_loiter;
T_W_loiter=thrust_required_loiter/MTOW;
P_W_loiter = Power_required_loiter/MTOW;

%appox_c = b/AR