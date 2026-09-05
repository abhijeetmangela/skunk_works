clc;
clear;
close all;

%% ============================================================
% CLOSED RECTANGULAR SECTION WITH 4 CORNER LONGERONS
% BOOM CONTRIBUTION INCLUDED IN SHEAR FLOW
% ============================================================

%% ---------------- INPUTS ----------------

Vy = -50;              % Vertical shear force [N]

b = 0.180;             % Width [m]
h = 0.250;             % Height [m]

t = 1e-3;              % Skin thickness [m]

A_L = 10e-6;           % Area of EACH longeron [m^2]
                       % = 10 mm^2

N = 300;               % Points per skin panel


%% ============================================================
% LONGERON LOCATIONS
% ============================================================

% Start cut at middle of top panel and travel clockwise.
%
%             CUT
%              |
%       ----------------
%       |              |
%       |              |
%       |              |
%       ----------------
%
% Corners:
%
% 1 = top-left
% 2 = bottom-left
% 3 = bottom-right
% 4 = top-right

xB = [-b/2, -b/2,  b/2,  b/2];
yB = [ h/2, -h/2, -h/2,  h/2];


%% ============================================================
% SECOND MOMENT OF AREA
% ============================================================

% ---- Skin contribution ----

Ixx_skin = ...
      2*(b*t)*(h/2)^2 ...
    + 2*t*h^3/12;


% ---- Longerons contribution ----

Ixx_long = sum(A_L .* yB.^2);


% ---- Total ----

Ixx = Ixx_skin + Ixx_long;


%% ============================================================
% CREATE CLOSED SKIN CONTOUR
% CUT AT MIDDLE OF TOP PANEL
% ============================================================

% Top-left half
x1 = linspace(0, -b/2, N);
y1 = (h/2)*ones(1,N);

% Left wall
x2 = (-b/2)*ones(1,N);
y2 = linspace(h/2, -h/2, N);

% Bottom wall
x3 = linspace(-b/2, b/2, N);
y3 = (-h/2)*ones(1,N);

% Right wall
x4 = (b/2)*ones(1,N);
y4 = linspace(-h/2, h/2, N);

% Top-right half
x5 = linspace(b/2, 0, N);
y5 = (h/2)*ones(1,N);


% Combine
x = [x1, ...
     x2(2:end), ...
     x3(2:end), ...
     x4(2:end), ...
     x5(2:end)];

y = [y1, ...
     y2(2:end), ...
     y3(2:end), ...
     y4(2:end), ...
     y5(2:end)];

x = x(:);
y = y(:);


%% ============================================================
% SHEAR FLOW CALCULATION
%
% q_b = -(Vy/Ixx) * Q
%
% Q = integral(y*t*ds) + SUM(A_i*y_i)
%
% Longeron contributions are added when the contour passes
% each corner longeron.
% ============================================================

qb = zeros(size(x));

Q_skin = 0;
Q_boom = 0;

% Track which longerons have already been encountered
boom_used = false(1,4);


for i = 2:length(x)

    % --------------------------------------------------------
    % Skin contribution over current small segment
    % --------------------------------------------------------

    dx = x(i) - x(i-1);
    dy = y(i) - y(i-1);

    ds = sqrt(dx^2 + dy^2);

    y_avg = 0.5*(y(i-1) + y(i));

    dQ_skin = y_avg * t * ds;

    Q_skin = Q_skin + dQ_skin;


    % --------------------------------------------------------
    % Check whether a longeron is encountered
    % --------------------------------------------------------

    for k = 1:4

        if ~boom_used(k)

            distance_to_boom = ...
                sqrt((x(i)-xB(k))^2 + ...
                     (y(i)-yB(k))^2);

            if distance_to_boom < ds*1.5

                % Boom contribution:
                %
                % Q_boom += A_L * y_boom

                Q_boom = Q_boom + A_L*yB(k);

                boom_used(k) = true;

            end
        end
    end


    % --------------------------------------------------------
    % Total first moment
    % --------------------------------------------------------

    Q = Q_skin + Q_boom;


    % --------------------------------------------------------
    % Basic shear flow
    % --------------------------------------------------------

    qb(i) = -(Vy/Ixx)*Q;

end


%% ============================================================
% CLOSED CELL CONSTANT SHEAR FLOW q0
% ============================================================

% Compatibility:
%
%      integral(q/t ds) = 0
%
% q = qb + q0
%
% therefore
%
% q0 =
% - integral(qb/t ds)
% --------------------
%   integral(ds/t)

int_qb_t = 0;
int_1_t  = 0;

for i = 2:length(x)

    dx = x(i)-x(i-1);
    dy = y(i)-y(i-1);

    ds = sqrt(dx^2 + dy^2);

    qb_avg = 0.5*(qb(i-1)+qb(i));

    int_qb_t = int_qb_t + qb_avg/t*ds;

    int_1_t = int_1_t + ds/t;

end


q0 = -int_qb_t/int_1_t;


%% ============================================================
% TOTAL SHEAR FLOW
% ============================================================

q = qb + q0;


%% ============================================================
% RESULTS
% ============================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf(' CLOSED RECTANGULAR SECTION WITH LONGERONS\n');
fprintf('====================================================\n');

fprintf('Vy                    = %.3f N\n', Vy);

fprintf('\nGeometry:\n');
fprintf('Width                 = %.1f mm\n', b*1000);
fprintf('Height                = %.1f mm\n', h*1000);
fprintf('Skin thickness        = %.2f mm\n', t*1000);

fprintf('\nLongeron:\n');
fprintf('Area of each longeron = %.2f mm^2\n', A_L*1e6);
fprintf('Number of longerons   = 4\n');

fprintf('\nSecond moments:\n');
fprintf('Ixx skin              = %.6e m^4\n', Ixx_skin);
fprintf('Ixx longerons         = %.6e m^4\n', Ixx_long);
fprintf('Ixx total             = %.6e m^4\n', Ixx);

fprintf('\nClosed-cell correction:\n');
fprintf('q0                    = %.6f N/m\n', q0);

fprintf('\nShear flow:\n');
fprintf('qb,max                = %.3f N/m\n', max(qb));
fprintf('qb,min                = %.3f N/m\n', min(qb));

fprintf('q,max                 = %.3f N/m\n', max(q));
fprintf('q,min                 = %.3f N/m\n', min(q));

fprintf('abs(q)_max            = %.3f N/m\n', max(abs(q)));

fprintf('====================================================\n');


%% ============================================================
% DISTANCE ALONG PERIMETER
% ============================================================

s = zeros(size(x));

for i = 2:length(x)

    ds = sqrt( ...
        (x(i)-x(i-1))^2 + ...
        (y(i)-y(i-1))^2 );

    s(i) = s(i-1) + ds;

end


%% ============================================================
% SHEAR FLOW PLOT
% ============================================================

figure('Color','w');

plot(s,q,'b-','LineWidth',2);

grid on;

xlabel('Distance along perimeter, s [m]');
ylabel('Shear flow, q [N/m]');

title('Shear Flow with Four 10 mm^2 Longerons');

yline(0,'k--');


%% ============================================================
% SECTION + LONGERONS
% ============================================================

figure('Color','w');

plot(x,y,'k-','LineWidth',1.5);
hold on;

plot(xB,yB,'ro',...
     'MarkerSize',9,...
     'MarkerFaceColor','r');

axis equal;
grid on;

xlabel('x [m]');
ylabel('y [m]');

title('Rectangular Closed Section with Corner Longerons');

legend('Skin','Longerons');


%% ============================================================
% DISPLAY BOOM CONTRIBUTIONS
% ============================================================

fprintf('\nLongeron Q contributions:\n');

for k = 1:4

    Qk = A_L*yB(k);

    fprintf(['Longeron %d: x = %.1f mm, y = %.1f mm, ', ...
             'A*y = %.6e m^3\n'], ...
             k, ...
             xB(k)*1000, ...
             yB(k)*1000, ...
             Qk);

end