% CONSTRAINT DIAGRAM TO FINDOUT WING LOADING

% AUTHOR: MONISHA VIJAYAN (AE25M037)
% Edited on 5/8/26 : Added proper service and absolute celings - Abhijeet.

rho_0=1.225;
rho_alt=0.9093; %Service celing = 3000 m
rho_abs=0.8194; %Absolute ceing = 4000 m
MTOW=68.516;
sigma_ceil=(rho_alt/rho_0);
sigma_ceil1=(rho_abs/rho_0);
g=9.81;
V_max=24;
CL_max=1.3;
S=0.413;
V_stall=sqrt((2*MTOW)/(rho_0*CL_max*S));
CD0=0.025;
AR=9.686;
e=0.8;
K=1/(pi*e*AR);
eta_p=0.77;
ROC=24.4720*sind(8);
L_D_max=15.54;
R_min=20;
S_TO=20;

ROC_ceil=0.508;
ROC_ceil1=0;
V_to= 1.2*V_stall;
CL_takeoff = (2*MTOW)/(rho_0*V_to^2*S);
mu=0.04;
WS=linspace(5,200,500);
% Stalling speed
WS_stall=0.5*rho_0*V_stall^2*CL_max;
% Maximum speed
WP_Vmax=eta_p./((0.5*rho_0*V_max^3*CD0./WS)+(2*K.*WS./(rho_0*V_max)));
% ROC
term_ROC=(1.155/(L_D_max))*sqrt((2.*WS)./(rho_alt*sqrt(3*CD0/K)));
WP_ROC=eta_p./(ROC+term_ROC);
% Takeoff ; CDG = CDO+KCL_to^2-muCL_to
CDG=CD0+K*(CL_takeoff)^2-mu*CL_takeoff;
exp_term=exp(0.6*rho_0*g*CDG*S_TO./WS);
WP_STO=((1-exp_term)./(mu-(mu+CDG/(CL_max/1.1^2)).*exp_term)).*(eta_p/V_to);
% Minimum Turn Radius
n=sqrt((V_max^2/(g*R_min))^2+1);
q_max=0.5*rho_0*V_max^2;
PW_Turn=(V_max/eta_p).*((q_max*CD0./WS)+(K*n^2.*WS./q_max));
WP_Turn=1./PW_Turn;
%Service ceiling
term_ceil=(1.155/(L_D_max*eta_p))*sqrt((2.*WS)./(rho_alt*sqrt(3*CD0/K)));
WP_Ceiling=sigma_ceil./((ROC_ceil/eta_p)+term_ceil);
%Absolute ceiling
term_ceil=(1.155/(L_D_max*eta_p))*sqrt((2.*WS)./(rho_abs*sqrt(3*CD0/K)));
WP_Ceiling1=sigma_ceil1./((ROC_ceil1/eta_p)+term_ceil);

%Maximum Stall Speed Constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_Vmax,'r','LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_Vmax(idx);
spacing = 5;  
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Maximum Speed Constraint','Color','k')
legend('Max Speed','Stall Limit')
ylim([0 0.5])
hold off

% Stalling Speed Constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
xline(WS_stall,'--r','LineWidth',2)
spacing = 5;
x = WS(WS<=WS_stall);
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 0.5], 'k','LineWidth',0.5)
end
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Stall Speed Constraint','Color','k')
legend('Stall Limit')
ylim([0 0.5])
xlim([0 max(WS)])
hold off

% ROC max constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_ROC,'b','LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_ROC(idx);
spacing = 5;
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Rate of Climb Constraint','Color','k')
legend('ROC','Stall Limit')
ylim([0 0.5])
hold off
% Takeoff constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_STO,'g','LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_STO(idx);
spacing = 5;
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Takeoff Constraint','Color','k')
legend('Takeoff','Stall Limit')
ylim([0 0.5])
hold off
% Turn radius constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_Turn,'m','LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_Turn(idx);
spacing = 5;
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Turn Radius Constraint','Color','k')
legend('Turn Radius','Stall Limit')
ylim([0 0.5])
hold off
% Service ceiling constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_Ceiling,'c','LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_Ceiling(idx);
spacing = 5;
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Service Ceiling Constraint','Color','k')
legend('Service Ceiling','Stall Limit')
ylim([0 2])
hold off
% Absolute ceiling constraint
figure('Color','w')
hold on
set(gca,'Color','w','XColor','k','YColor','k','GridColor','w'); grid on
plot(WS,WP_Ceiling1,'Color',[0.5 0 0.8],'LineWidth',2)
idx = WS <= WS_stall;
x = WS(idx);
y = WP_Ceiling1(idx);
spacing = 5;
for i = 1:spacing:length(x)
    plot([x(i) x(i)], [0 y(i)], 'k','LineWidth',0.5)
end
xline(WS_stall,'--k','LineWidth',2)
xlabel('Wing Loading W/S (N/m^2)','FontWeight','bold')
ylabel('Power Loading W/P (N/W)','FontWeight','bold')
title('Absolute Ceiling Constraint','Color','k')
legend('Absolute Ceiling','Stall Limit')
ylim([0 2])
hold off

%All Constraints
figure
fig=figure('Color','w');set(gcf,'InvertHardcopy','off');
ax=gca;ax.Color='w';ax.XColor='k';ax.YColor='k';ax.GridColor='k';ax.GridAlpha=0.2;ax.LineWidth=1.5;hold on;grid on;
plot(WS,WP_Vmax,'r','LineWidth',2,'DisplayName','MaxSpeed');
plot(WS,WP_ROC,'b','LineWidth',2,'DisplayName','ROC');
plot(WS,WP_STO,'g','LineWidth',2,'DisplayName','Take-off');
plot(WS,WP_Turn,'m','LineWidth',2,'DisplayName','TurnRadius');
plot(WS,WP_Ceiling,'c','LineWidth',2,'DisplayName','Service Ceiling');
plot(WS,WP_Ceiling1,'Color',[0.5 0 0.8],'LineWidth',2,'DisplayName','Absolute Ceiling');
xline(WS_stall,'--k','LineWidth',2,'DisplayName','StallLimit');
xlabel('Wing Loading (W/S) [N/m^2]','FontWeight','bold');
ylabel('Power Loading (W/P) [N/W]','FontWeight','bold');
title('Wing Loading Constraint Diagram','FontSize',14,'Fontweight','bold','Color','k');
ylim([0,0.5])
WS_design = 121.844;
WP_design = interp1(WS,WP_STO,WS_design);
plot(WS_design,WP_design,'ko','MarkerFaceColor','k','MarkerSize',8,'DisplayName','Design Point')
text(WS_design,WP_design+0.01,' Design Point','Color','k','FontWeight','bold')
legend('Location','northwest','TextColor','k','Color','w');
hold off;

% NEW SURFACE AREA
b=2;
W_S=121.844;
W=68.516;
S_new=W/(W_S)
AR=(b^2)/S_new