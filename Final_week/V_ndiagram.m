% FLIGHT ENVELOPE DIAGRAM
% AUTHOR : MONISHA VIJAYAN(AE25M037)
clear;clc;close all;
W=(6.592)*9.81; %(6.592-1.2) after payload drop
S=0.5307;
rho=1.220;
CL_max_pos=1.3;
CL_max_neg=-0.8;
n_max_pos=3.5;
n_max_neg=-1.5;
V_s=sqrt((2*W)/(rho*S*CL_max_pos));
V_max=24;
V_d=1.25*V_max;
V_a=V_s*sqrt(n_max_pos);
V_g=sqrt((2*abs(n_max_neg)*W)/(rho*S*abs(CL_max_neg)));
V=0:0.1:V_d;
n_curve_pos=(rho.*V.^2.*S.*CL_max_pos)/(2*W);
n_curve_neg=(rho.*V.^2.*S.*CL_max_neg)/(2*W);
n_curve_pos(n_curve_pos>n_max_pos)=n_max_pos;
n_curve_pos(V>V_a)=NaN;
n_curve_neg(V>V_g)=NaN;
figure('Color','w')
ax=gca;ax.Color='w';ax.XColor='k';ax.YColor='k';hold on;grid on;ax.GridColor='k';
plot(V,n_curve_pos,'b','LineWidth',2);
plot(V,n_curve_neg,'r','LineWidth',2);
line([V_a V_d],[n_max_pos n_max_pos],'Color','b','LineStyle','--','LineWidth',2);
line([V_g V_d],[n_max_neg n_max_neg],'Color','r','LineStyle','--','LineWidth',2);
line([V_max V_max],[n_max_neg n_max_pos],'Color','k','LineStyle','--');
line([V_d V_d],[n_max_neg n_max_pos],'Color','m','LineWidth',3);
plot(V_a,n_max_pos,'ko','MarkerFaceColor','g');
plot(V_g,n_max_neg,'ko','MarkerFaceColor','k');
plot(V_d,0,'ko','MarkerFaceColor','m');
xlabel('Velocity (V) [m/s]','Color','k','FontWeight','bold');
ylabel('Load Factor (n)','Color','k','FontWeight','bold');
title('UAV V-n Diagram','Color','k');
text(V_a,n_max_pos+0.2,'V_a','Color','k');
text(V_d,0.2,'V_d','Color','k');
text(V_g,n_max_neg+0.2,'V_g','Color','k')
axis([0 V_d+5 n_max_neg-0.5 n_max_pos+0.5]);
fprintf('Vs = %.2f m/s\n',V_s);
fprintf('Vmax = %.2f m/s\n',V_max);
fprintf('Vd = %.2f m/s\n',V_d);
fprintf('Va = %.2f m/s\n',V_a);
fprintf('Vg = %.2f m/s\n',V_g);
