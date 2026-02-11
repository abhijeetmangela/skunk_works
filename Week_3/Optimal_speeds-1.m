 
%%      OPTIMAL SPEED ESTIMATION FOR DIFFERENT PHASES

%Cruise,Max Range,"Minimum Thrust (CDo=CDi​)"
%Climb,Max Angle,"Minimum Thrust (CDo=CDi​)"
%Loiter,Max Time,"Minimum Power (3CDo=CDi​)"
%Descent,Max Glide,"Minimum Thrust (CDo=CDi)"

% ==========================
% AUTHOR = MONISHA VIJAYAN (AE25M037)
%==========================



% BEFORE PAYLOAD DROP
% ==========================
% For our UAV
% ==========================
AR=9.686;
MTOW=68.516;
e=1.78*(1-0.045*AR^0.68)-0.64; %for rectangular wing
rho=1.225; % @ cruising altitude = 200m
b=2;
s=(b^2)/AR;
vcr=22;
Cl=(2*MTOW)/(rho*vcr^2*s);
k=1/(pi*AR*e); 
CDo=k*Cl^2; 
v0=30; 
%% =========================
%  Optimum velocity during CLIMB1  
% ==========================
gamma=8;
Clclimb=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Clclimb(v).^2);
v_optimal_climb=fsolve(fun,v0);
D=(0.5*rho*s*v_optimal_climb^2*2*k*Clclimb(v_optimal_climb)^2);
T_Wcl=(sind(gamma)+(D/MTOW));
T=(T_Wcl*MTOW);
Pclimb1=(T*v_optimal_climb);
fprintf('Optimum velocity during climb1 = %.3f\n',v_optimal_climb)
fprintf('T/W ratio for climb1 phase = %.3f\n',T_Wcl)
fprintf('Power required in Climb1 phase = %.3f\n',Pclimb1)
%% =========================
%  Optimum velocity during CRUISE1 
% ==========================
gamma=0;
Clcruise=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Clcruise(v).^2);
v_optimal_cruise=fsolve(fun,v0);
D=(0.5*rho*s*v_optimal_cruise^2*2*k*Clcruise(v_optimal_cruise)^2);
T_Wcr=(sind(gamma)+(D/MTOW));
T=(T_Wcr*MTOW);
Pcruise1=(T*v_optimal_cruise);
fprintf("Optimum velocity during cruise1 = %.3f\n",v_optimal_cruise)
fprintf("T/W ratio for cruise1 phase = %.3f\n",T_Wcr)
fprintf("Power required in cruise1 phase = %.3f\n",Pcruise1)
%% =========================
%  Optimum velocity during DESCENT1 
% ==========================
gamma=-8;
Cldescent=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Cldescent(v).^2);
v_optimal_descent=fsolve(fun,v0);
D=(0.5*rho*s*v_optimal_descent^2*2*k*Cldescent(v_optimal_descent)^2);
T_Wds=abs((sind(gamma)+(D/MTOW)));
T=(T_Wds*MTOW);
Pdescent1=(T*v_optimal_descent);
fprintf("Optimum velocity during descent1 = %.3f\n",v_optimal_descent)
fprintf("T/W ratio for descent1 phase = %.3f\n",T_Wds)
fprintf("Power required in descent1 phase = %.3f\n",Pdescent1)

%%==========================
% AFTER PAYLOAD DROP
% ==========================
% For our UAV
% ==========================
MTOW=57.94;
Cl_ref=(2*68.516)/(rho*vcr^2*s);
CDo=k*Cl_ref^2; 
L_Dm=1/(2*sqrt(k*CDo)); 
v0=30; 
%% =========================
%  Optimum velocity during LOITER
% ==========================
gamma=0;
Clloiter=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(3*0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Clloiter(v).^2);
v_optimal_loiter=fsolve(fun,v0);
D=0.5*rho*s*v_optimal_loiter^2*(CDo+k*Clloiter(v_optimal_loiter)^2);
T_Wl=(sind(gamma)+(D/MTOW));
T=(T_Wl*MTOW);
Ploiter=(T*v_optimal_loiter);
fprintf("Optimum velocity during loiter = %.3f\n",v_optimal_loiter)
fprintf("T/W ratio for loiter phase = %.3f\n",T_Wl)
fprintf("Power required in loiter phase = %.3f\n",Ploiter)
%% =========================
%  Optimum velocity during CLIMB2
% ==========================
gamma=8;
Clclimb=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Clclimb(v).^2);
v_optimal_climb=fsolve(fun,v0);
D=0.5*rho*s*v_optimal_climb^2*2*CDo;
T_Wcl=(sind(gamma)+(D/MTOW));
T=(T_Wcl*MTOW);
Pclimb2=(T*v_optimal_climb);
fprintf("Optimum velocity during climb2 = %.3f\n",v_optimal_climb)
fprintf("T/W ratio for climb2 phase = %.3f\n",T_Wcl)
fprintf("Power required in Climb2 phase = %.3f\n",Pclimb2)
%% =========================
%  Optimum velocity during CRUISE2
% ==========================
gamma=0;
Clcruise=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Clcruise(v).^2);
v_optimal_cruise=fsolve(fun,v0);
D=0.5*rho*s*v_optimal_cruise^2*2*CDo;
T_Wcr=(sind(gamma)+(D/MTOW));
T=(T_Wcr*MTOW);
Pcruise2=(T*v_optimal_cruise);
fprintf("Optimum velocity during cruise2 = %.3f\n",v_optimal_cruise)
fprintf("T/W ratio for cruise2 phase = %.3f\n",T_Wcr)
fprintf("Power required in cruise2 phase = %.3f\n",Pcruise2)
%% =========================
%  Optimum velocity during DESCENT2
% ==========================
gamma=-8;
Cldescent=@(v)(2*MTOW*cosd(gamma))/(rho*s*v.^2);
fun=@(v)(0.5*rho*s*v.^2*CDo)-(0.5*rho*s*v.^2.*k.*Cldescent(v).^2);
v_optimal_descent=fsolve(fun,v0);
D=0.5*rho*s*v_optimal_descent^2*2*CDo;
T_Wds=abs((sind(gamma)+(D/MTOW)));
T=(T_Wds*MTOW);
Pdescent2=(T*v_optimal_descent);
fprintf("Optimum velocity during descent2 = %.3f\n",v_optimal_descent)
fprintf("T/W ratio for descent2 phase = %.3f\n",T_Wds)
fprintf("Power required in descent2 phase = %.3f\n",Pdescent2)
