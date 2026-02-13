%% (L/D)max estimation of similar UAVs

%% AUTHOR = MONISHA VIJAYAN (AE25M037) 
% UAV Names
e=[];
Name = {
    'Yangda Mapird Pro'
    'Yangda Mapird Plus'
    'Yangda Nimbus Pro'
    'Yangda FW-250'
    'Avy Aera'
    'Foxtech Baby Shark 260'
    'P330 Pro'
    'Foxtech AYK-250'
    'Wingcopter 178'
    'FlyDragon FLY-2100'
    'Black Swift S2'
    'Q200'
    'SkyEye Sierra'
    'Foxtech Nimbus'
};
% MTOW (kg)
MTOW = [ ...
    7.5 11.5 8.6 12.0 12.0 13.0 14.0 15.0 ...
    18.0 7.0 9.5 7.0 12.5 6.0 ];
% Wing span(m)
b = [...
    2.1 2.43 1.95 2.5 2.4 2.5 2.53 2.5 ...
    1.78 2.1 3 1.95 2 1.8];
% Wing area (m^2)
Area = [ ...
    0.50 0.65 0.46 0.52 0.96 0.50 0.56 0.50 ...
    0.58 0.52 0.70 0.60 0.30 0.38 ];
% Aspect ratio
AR = [ ...
    8.82 9.08 8.27 12.02 6.00 12.50 11.43 12.50 ...
    5.46 8.48 12.86 6.34 13.33 8.53 ];
% Cruising Altitude (m)
hcruise = [ ...
    3500 3500 3500 3500 3100 3000 5000 2000 ...
    4500 500 4250 4000 3000 1500 ];
% Cruising speed (m/s)
vcruise = [ ...
    17.00 19.44 20.00 23.33 26.39 22 21.00 ...
    25.00 24.00 18.5 18.06 18 20 20];
%Swet/Sref 
sratio =[...
    3.7 3.7 3.7 3.4 2.20 2.3 2.20 2.20...
    3.6 3.7 2.3 3.6 3.6 3.7];
swet=[];
ARwet=[];
for i = 1:length(sratio)
    swt=sratio(i)*Area(i);
    swet(end+1)=swt;
    ARwt=(b(i)^2)/swet(i);
    ARwet(end+1)=ARwt;
end


%Planform type
Planform = {
    'Rectangular'
    'Rectangular'
    'Rectangular'
    'Tapered'
    'Tapered'
    'Tapered'
    'Tapered'
    'Tapered'
    'Rectangular'
    'Rectangular'
    'Tapered'
    'Rectangular'
    'Rectangular'
    'Rectangular'
};

% Taper Ratio
TR=[...
    1 1 1 0.141 0.287 0.119 0.149 0.3103 ...
    1 1 0.636 1 1 1];
%------------------------------------------------------------------------
%Oswald Efficiency Factor
for i = 1:length(MTOW)
    if strcmp(Planform{i},'Rectangular')
        eu(i) = 1.78*(1 - 0.045*AR(i)^0.68) - 0.64;

    elseif strcmp(Planform{i},'Tapered')
        d(i) = (1-TR(i))/(2*(1+TR(i)));
        eu(i) = 1/(1+d(i));
    end
    e(end+1)=eu(i);
end
e = max(e, 0.7);
e = min(e, 1.0);
%disp(e)
%------------------------------------------------------------------------
%ISA Parameters
h = hcruise
rho = NaN(size(h));    
T0 = 288.15;             
P0 = 101325;              
rho0 = 1.225;             
L = 0.0065;               
R = 287.05;               
g = 9.80665;              
%------------------------------------------------------------------------
%Density Calculation
for i = 1:length(h)
    if ~isnan(h(i))
        T = T0 - L*h(i);
        rho(i) = rho0*(T/T0)^(g/(R*L) - 1);
    end
end

%disp(rho)
%------------------------------------------------------------------------
%L/D max calculation
L_Dm = [];

for i=1:length(MTOW)
    Cl=(2*MTOW(i)*9.81)/(rho(i)*Area(i)*(vcruise(i)^2));
    L_D=(pi*e(i)*AR(i))/(2*Cl);
    L_Dm(end+1)=L_D;
end
%------------------------------------------------------------------------
%for our UAV
sratio1=3.6;
Sref=0.413;
swet=sratio1*Sref;
ARwsk =4/swet;
disp(sqrt(ARwsk))

%------------------------------------------------------------------------
%[AR_sorted, idx] = sort(sqrt(ARwet));
%L_dm_sorted = L_Dm(idx);
x= linspace(1,4,100);
fit =polyfit(ARwet,L_Dm,1)
y=polyval(fit,x)
figure
scatter(sqrt(ARwet),L_Dm,40,"green","o",'LineWidth',1.5)
hold on
plot(x,y)


