clc; clear; close all;

Wo = [15;25;7;9.5;5.1;7;17;26]; % MTOW data
We = [7;12.2;3.4;5.7;3.2;3.4;12.64;8.7]; % Empty weight data

% Borey 10
% Penguin C Mk2
% FLydragon
% Black swift s2
% Mq9
% Q200
% fire fly
% Borey 20

Y = We./Wo;
X = Wo; % To use power law i.e We/Wo = a*Wo^L

logX = log(X);
logY = log(Y);

A = [logX ones(length(X),1)];
b = logY;
Soln = A\b;

L = Soln(1);
loga = Soln(2);
a = exp(loga);

x = linspace(1.5,3.5,100);
y = linspace(0,30,200);

% parameters
rho = 1.225;
g = 9.81;
b = 2.25;
AR = 8;
S = b^2/AR;

V_cr1 = 16; V_cr2 = 16;
V_cl1 = 13.5; V_cl2 = 13.5;
V_de1 = 20; V_de2 = 20;
V_lo1 = 20; V_lo2 = 20; V_lo3 = 20;

t_cr1 = 2500; t_cr2 = 2500;
t_cl1 = 140;  t_cl2 = 115;
t_de1 = 100;  t_de2 = 130;
t_lo1 = 120;  t_lo2 = 120;  t_lo3 = 120;

gamma_cl = 15*pi/180;
gamma_de = 15*pi/180;

turn_rate = pi/60;
phi = atan(turn_rate*V_lo1/g);

CDo = 0.03;
e = 0.8;

Eta_p = 0.7;
Bat_eff = 0.95;
Bat_dis = 0.75;

mAh = 15000;
Vol = 22.2;
W_Ib = 1.3;

W_o = 7;
W_p = 2.13;

i = 1;
W_o_array = [];

while i < 100

    W_b = getbatteryweight( ...
        W_o,W_p,rho,g,S,V_cr1,V_cr2,V_cl1,V_cl2,V_de1,V_de2, ...
        V_lo1,V_lo2,V_lo3,t_cr1,t_cr2,t_cl1,t_cl2,t_de1,t_de2, ...
        t_lo1,t_lo2,t_lo3,gamma_cl,gamma_de,phi,CDo,AR,e, ...
        Eta_p,Bat_eff,Bat_dis,mAh,Vol,W_Ib);

    W_o = W_p/(1 - W_b/W_o - a*W_o^L);

    W_o_array(i) = W_o;
    iterations(i) = i;

    i = i + 1;
end


%% Figure 1 — log regression
figure(1)
scatter(logX,logY,36,"red","o","MarkerEdgeColor","flat",LineWidth=2)
hold on
plot(x,L*x+loga,LineStyle="-",Color="blue",LineWidth=3)

xlabel('log(Wo)','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
ylabel('log(We/Wo)','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
title('log(We/Wo) vs log(Wo)','FontSize',21,'FontWeight','bold','FontName',"FixedWidth")

legend('data','regression')
grid on
grid minor

ax = gca;
ax.LineWidth = 1.5;
ax.FontWeight = "bold";
box on


%% Figure 2 — iteration convergence
figure(2)
plot(iterations,W_o_array,LineStyle="-",Color="blue",LineWidth=3)

xlabel('no. of iterations','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
ylabel('Wo','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
title('MTOW vs no. of iterations','FontSize',21,'FontWeight','bold','FontName',"FixedWidth")

grid on
grid minor

ax = gca;
ax.LineWidth = 1.5;
ax.FontWeight = "bold";
box on


%% Figure 3 — power law fit
figure(3)
scatter(X,Y,36,"red","o","MarkerEdgeColor","flat",LineWidth=2)
hold on
plot(y,a*y.^L,LineStyle="-",Color="blue",LineWidth=3)

xlabel('Wo','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
ylabel('We/Wo','FontSize',21,'FontName',"FixedWidth","FontWeight","bold")
title('We/Wo vs Wo','FontSize',21,'FontWeight','bold','FontName',"FixedWidth")

legend('data','regression')
grid on
grid minor

ax = gca;
ax.LineWidth = 1.5;
ax.FontWeight = "bold";
box on
