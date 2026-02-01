function W_b = getbatteryweight(W_o,W_p,rho,g,S,V_cr1,V_cr2,V_cl1,V_cl2,V_de1,V_de2,V_lo1,V_lo2,V_lo3,t_cr1,t_cr2,t_cl1,t_cl2,t_de1,t_de2,t_lo1,t_lo2,t_lo3,gamma_cl,gamma_de,phi,CDo,AR,e,Eta_p,Bat_eff,Bat_dis,mAh,Vol,W_Ib)

CL_cr1 = W_o*g/(0.5*rho*V_cr1^2*S);
CL_cr2 = (W_o-W_p)*g/(0.5*rho*V_cr2^2*S);

CL_cl1 = W_o*g*cos(gamma_cl)/(0.5*rho*V_cl1^2*S);
CL_cl2 = (W_o - W_p)*g*cos(gamma_cl)/(0.5*rho*V_cl2^2*S);

CL_de1 = W_o*g*cos(gamma_de)/(0.5*rho*V_de1^2*S);
CL_de2 = (W_o - W_p)*g*cos(gamma_de)/(0.5*rho*V_de2^2*S);

CL_lo1 = W_o*g/(0.5*rho*V_lo1^2*S*cos(phi));
CL_lo2 = (W_o - W_p)*g/(0.5*rho*V_lo2^2*S*cos(phi));
CL_lo3 = (W_o - W_p)*g/(0.5*rho*V_lo3^2*S*cos(phi));

CD_cr1 = CDo + CL_cr1^2/(pi*AR*e);
CD_cr2 = CDo + CL_cr2^2/(pi*AR*e);

CD_cl1 = CDo + CL_cl1^2/(pi*AR*e);
CD_cl2 = CDo + CL_cl2^2/(pi*AR*e);

CD_de1 = CDo + CL_de1^2/(pi*AR*e);
CD_de2 = CDo + CL_de2^2/(pi*AR*e);

CD_lo1 = CDo + CL_lo1^2/(pi*AR*e);
CD_lo2 = CDo + CL_lo2^2/(pi*AR*e);
CD_lo3 = CDo + CL_lo3^2/(pi*AR*e);

D_cr1 = 0.5*rho*V_cr1^2*S*CD_cr1;
D_cr2 = 0.5*rho*V_cr2^2*S*CD_cr2;
D_cl1 = 0.5*rho*V_cl1^2*S*CD_cl1;
D_cl2 = 0.5*rho*V_cl2^2*S*CD_cl2;

D_de1 = 0.5*rho*V_de1^2*S*CD_de1;
D_de2 = 0.5*rho*V_de2^2*S*CD_de2;

D_lo1 = 0.5*rho*V_lo1^2*S*CD_lo1;
D_lo2 = 0.5*rho*V_lo2^2*S*CD_lo2;
D_lo3 = 0.5*rho*V_lo3^2*S*CD_lo3;

T_cr1 = D_cr1;
T_cr2 = D_cr2;

T_cl1 = D_cl1 + W_o*g*sin(gamma_cl);
T_cl2 = D_cl2 + (W_o - W_p)*g*sin(gamma_cl);

T_de1 = D_de1 - W_o*g*sin(gamma_de);
T_de2 = D_de2 - (W_o - W_p)*g*sin(gamma_de);

T_lo1 = D_lo1;
T_lo2 = D_lo2;
T_lo3 = D_lo3;

P_cr1 = T_cr1*V_cr1;
P_cr2 = T_cr2*V_cr2;

P_cl1 = T_cl1*V_cl1;
P_cl2 = T_cl2*V_cl2;

P_de1 = T_de1*V_de1;
P_de2 = T_de2*V_de2;

P_lo1 = T_lo1*V_lo1;
P_lo2 = T_lo2*V_lo2;
P_lo3 = T_lo3*V_lo3;

if P_de1 < 0
    P_de1 = 0;
end

if P_de2 < 0
    P_de2 = 0;
end

E_cr = P_cr2*t_cr2 + P_cr1*t_cr1;
E_cl = P_cl2*t_cl2 + P_cl1*t_cl1;
E_de = P_de2*t_de2 + P_de1*t_de1;
E_lo = P_lo1*t_lo1 + P_lo2*t_lo2 + P_lo3*t_lo3;

E_total = (E_cr + E_cl + E_de + E_lo)/Eta_p;
E_Bat = E_total/Bat_eff; % energy consumed from battery

Bat_cap = E_Bat/(3600*Bat_dis); % total battery capacity needed in Watt Hour

No_Bat = ceil(Bat_cap/(mAh*Vol*10^(-3)));

W_b = No_Bat*W_Ib;
