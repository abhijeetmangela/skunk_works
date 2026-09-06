% Author: Abhijeet Mangela

% cross section at spar location
h = 0.2 % m
b = 0.2 % m

t_skin = 0.5e-3
t_boom = 1e-3
V_y = -40
total_lift = 6.8*9.8

% We wiil assume distributed lift at the top
lift_distrib = (total_lift/2)/(b/2)

Area_boom = 20e-3*t_boom

I_skin = (1/12)*( b*h^3 - (b-2*t_skin)*(h-2*t_skin)^3)
I_boom = 4*Area_boom*(h/2)^2

I_total = I_skin + I_boom

% C_section bulkhead

length_middle = 15e-3
length_flange = 25e-3
As = (2*length_flange + length_middle)*t_boom

E = 70e9;       % Young's modulus [Pa]
G = 26e9;       % Shear modulus [Pa]



syms s n m f 

% n - Axial unknown
% f - Shear unknown
% m - moment unknown


% From 0 to 1
L1 = b/2

Axial_1 = n + (V_y*t_skin*h/(2*I_total))*(s^2/2) % Integral of Shear flow
Shear_1 = f 
Moment_1 = m - f*s

% loads at junction 1 

m_1 = m - f*b/2
n_1 = n + (V_y*t_skin*h/(2*I_total))*b^2/8
f_1 = f

% From 1 to 2

L2 = h

Axial_2 = f_1 + (V_y/I_total)*( (b*t_skin*h/4 + Area_boom*h/2)*s + h*t_skin*s^2/2 - t_skin*s^3/3 )
Shear_2 = -n_1
Moment_2 = m_1 + n_1*s

% loads at junction 2

L3 = b/2

m_2 = subs(Moment_2,s,h)
n_2 = subs(Axial_2,s,h)
f_2 = subs(Shear_2,s,h)

% From 2 to 3

Axial_3 = f_2 + V_y*t_skin*h/(4*I_total)*(b*s - s^2)
Shear_3 = lift_distrib*s - n_2
Moment_3 = lift_distrib*s^2/2 - m_2 - n_2*s


U1_bending = int( Moment_1^2/(2*E*I_total), s,0,L1)

U1_axial = int( Axial_1^2/(2*E*As), s,0,L1)

U1_Shear = int( Shear_1^2/(2*G*As), s,0,L1)

U1 = U1_bending + U1_axial + U1_Shear


U2_bending = int( Moment_2^2/(2*E*I_total), s,0,L2)

U2_axial = int( Axial_2^2/(2*E*As), s,0,L2)

U2_Shear = int( Shear_2^2/(2*G*As), s,0,L2)

U2 = U2_bending + U2_axial + U2_Shear

U3_bending = int( Moment_3^2/(2*E*I_total), s,0,L3)

U3_axial = int( Axial_3^2/(2*E*As), s,0,L3)

U3_Shear = int( Shear_3^2/(2*G*As), s,0,L3)

U3 = U3_bending + U3_axial + U3_Shear

U = simplify(U1 + U2 + U3)

% dU/dn = 0
% dU/df = 0
% dU/dm = 0

eq1 = diff(U,n) == 0
eq2 = diff(U,f) == 0
eq3 = diff(U,m) == 0

sol = solve([eq1,eq2,eq3],[n,f,m])


n_sol = double(sol.n)
f_sol = double(sol.f)
m_sol = double(sol.m)


N1 = subs(Axial_1,[n,f,m],[n_sol,f_sol,m_sol])
V1 = subs(Shear_1,[n,f,m],[n_sol,f_sol,m_sol])
M1 = subs(Moment_1,[n,f,m],[n_sol,f_sol,m_sol])

N2 = subs(Axial_2,[n,f,m],[n_sol,f_sol,m_sol])
V2 = subs(Shear_2,[n,f,m],[n_sol,f_sol,m_sol])
M2 = subs(Moment_2,[n,f,m],[n_sol,f_sol,m_sol])

N3 = subs(Axial_3,[n,f,m],[n_sol,f_sol,m_sol])
V3 = subs(Shear_3,[n,f,m],[n_sol,f_sol,m_sol])
M3 = subs(Moment_3,[n,f,m],[n_sol,f_sol,m_sol])

Npts = 1000

s1 = linspace(0,L1,Npts)
s2 = linspace(0,L2,Npts)
s3 = linspace(0,L3,Npts)

N1_num = double(subs(N1,s,s1));
V1_num = double(subs(V1,s,s1));
M1_num = double(subs(M1,s,s1));

N2_num = double(subs(N2,s,s2));
V2_num = double(subs(V2,s,s2));
M2_num = double(subs(M2,s,s2));

N3_num = double(subs(N3,s,s3));
V3_num = double(subs(V3,s,s3));
M3_num = double(subs(M3,s,s3));

x1 = s1;

x2 = L1 + s2;

x3 = L1 + L2 + s3;


x_all = [x1,x2,x3];

N_all = [N1_num,N2_num,N3_num];
V_all = [V1_num,V2_num,V3_num];
M_all = [M1_num,M2_num,M3_num];

N_max = max(abs(N_all))
V_max = max(abs(V_all))
M_max = max(abs(M_all))

fprintf('Maximum |N| = %.6f N\n',N_max);
fprintf('Maximum |V| = %.6f N\n',V_max);
fprintf('Maximum |M| = %.6f N-m\n',M_max);


% ============================================================
% PLOT AXIAL FORCE
% ============================================================

figure;

plot(x_all,N_all,'LineWidth',2);

xline(L1,'--','LineWidth',1.5);
xline(L1+L2,'--','LineWidth',1.5);

grid on;

xlabel('Global coordinate x [m]');
ylabel('Axial force N [N]');
title('Axial Force Distribution');

exportgraphics(gcf,"Axial_force.pdf",'ContentType','vector')



% ============================================================
% PLOT SHEAR FORCE
% ============================================================

figure;

plot(x_all,V_all,'LineWidth',2);

xline(L1,'--','LineWidth',1.5);
xline(L1+L2,'--','LineWidth',1.5);

grid on;

xlabel('Global coordinate x [m]');
ylabel('Shear force V [N]');
title('Shear Force Distribution');

exportgraphics(gcf,"Shear_force.pdf",'ContentType','vector')


% ============================================================
% PLOT BENDING MOMENT
% ============================================================

figure;

plot(x_all,M_all,'LineWidth',2);

xline(L1,'--','LineWidth',1.5);
xline(L1+L2,'--','LineWidth',1.5);

grid on;

xlabel('Global coordinate x [m]');
ylabel('Bending moment M [N-m]');
title('Bending Moment Distribution');

exportgraphics(gcf,"Bending_force.pdf",'ContentType','vector')



% STRESS CALCULATION

% Effective area used for axial stress
A_eff = As;

% Distance from neutral axis to extreme fibre
y_max = length_flange;

% Material yield strength
sigma_yield = 276e6;       % [Pa]

% NORMAL STRESS FROM AXIAL FORCE

sigma_N1 = N1_num/A_eff;
sigma_N2 = N2_num/A_eff;
sigma_N3 = N3_num/A_eff;

% BENDING STRESS
I_y_c_section = t_boom*length_flange^3/6

sigma_B1 = -M1_num*y_max/I_y_c_section;
sigma_B2 = -M2_num*y_max/I_y_c_section;
sigma_B3 = -M3_num*y_max/I_y_c_section;

% TOTAL NORMAL STRESS

sigma_1 = sigma_N1 + sigma_B1;
sigma_2 = sigma_N2 + sigma_B2;
sigma_3 = sigma_N3 + sigma_B3;

% SHEAR STRESS

% Approximate effective shear area
A_shear = As;

tau_1 = V1_num/A_shear;
tau_2 = V2_num/A_shear;
tau_3 = V3_num/A_shear;

% VON MISES STRESS

sigma_vm1 = sqrt(sigma_1.^2 + 3*tau_1.^2);
sigma_vm2 = sqrt(sigma_2.^2 + 3*tau_2.^2);
sigma_vm3 = sqrt(sigma_3.^2 + 3*tau_3.^2);

% COMBINE ALL SECTIONS

sigma_N_all = [sigma_N1 sigma_N2 sigma_N3];

sigma_B_all = [sigma_B1 sigma_B2 sigma_B3];

sigma_all = [sigma_1 sigma_2 sigma_3];

tau_all = [tau_1 tau_2 tau_3];

sigma_vm_all = [sigma_vm1 sigma_vm2 sigma_vm3];

% MAXIMUM STRESSES


sigma_N_max = max(abs(sigma_N_all));

sigma_B_max = max(abs(sigma_B_all));

sigma_max = max(abs(sigma_all));

tau_max = max(abs(tau_all));

sigma_vm_max = max(sigma_vm_all);


fprintf('STRESS RESULTS\n');


fprintf('Maximum axial stress       = %.6f MPa\n',sigma_N_max/1e6);

fprintf('Maximum bending stress     = %.6f MPa\n',sigma_B_max/1e6);

fprintf('Maximum normal stress      = %.6f MPa\n',sigma_max/1e6);

fprintf('Maximum shear stress       = %.6f MPa\n',tau_max/1e6);

fprintf('Maximum von Mises stress   = %.6f MPa\n',sigma_vm_max/1e6);

FoS = sigma_yield/sigma_vm_max;

fprintf('Yield strength             = %.6f MPa\n',sigma_yield/1e6);
fprintf('Factor of Safety            = %.4f\n',FoS);
