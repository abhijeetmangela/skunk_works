clc;
clear;
close all;

%% ============================================================
%      SHEAR FLOW OF A SINGLE-CELL CLOSED RECTANGULAR SECTION
% =============================================================

%% -------------------- INPUTS --------------------

Vy = -50;              % Vertical shear force [N]

b = 0.180;             % Width of section [m]
h = 0.250;             % Height of section [m]

t = 1e-3;              % Skin thickness [m]

N = 1000;              % Number of points per panel

%% ------------------------------------------------
% Coordinate system:
%
%        y
%        ^
%        |
%   (-b/2,h/2) -------- (b/2,h/2)
%        |                  |
%        |                  |
%        |                  |
%   (-b/2,-h/2) ------- (b/2,-h/2)
%
% x-axis -> right
% y-axis -> upward
%
% Shear force Vy acts vertically.
% -------------------------------------------------


%% -------------------- DISCRETIZE SECTION --------------------

% Top panel: Left -> Right
x_top = linspace(-b/2, b/2, N);
y_top = (h/2) * ones(1,N);

% Right panel: Top -> Bottom
x_right = (b/2) * ones(1,N);
y_right = linspace(h/2, -h/2, N);

% Bottom panel: Right -> Left
x_bottom = linspace(b/2, -b/2, N);
y_bottom = (-h/2) * ones(1,N);

% Left panel: Bottom -> Top
x_left = (-b/2) * ones(1,N);
y_left = linspace(-h/2, h/2, N);


%% -------------------- COMPLETE CLOSED CONTOUR --------------------

x = [x_top, ...
     x_right(2:end), ...
     x_bottom(2:end), ...
     x_left(2:end)];

y = [y_top, ...
     y_right(2:end), ...
     y_bottom(2:end), ...
     y_left(2:end)];

% Remove duplicate final point if required
x = x(:);
y = y(:);


%% -------------------- CALCULATE SECTION PROPERTIES --------------------

% Differential contour length
ds = sqrt(diff(x).^2 + diff(y).^2);

% Approximate enclosed-wall area
A = 0.5 * abs(sum(x(1:end-1).*y(2:end) - ...
                  x(2:end).*y(1:end-1)));

% Second moment of area about x-axis
% Thin-walled approximation:
Ixx = t * trapz(y.^2, ...
        [0; cumsum(ds)]);

% More convenient exact result for rectangular thin wall:
Ixx_exact = 2*(b*t)*(h/2)^2 + ...
            2*t*h^3/12;

Ixx = Ixx_exact;


%% -------------------- BASIC SHEAR FLOW --------------------
%
% For vertical shear:
%
%        qb = -(Vy / Ixx) * integral(y*t*ds)
%
% This is the thin-walled equivalent of
%
%        q = VQ/I
%
% with the sign determined by the chosen contour direction.
%
% We start the cut at the middle of the top panel.
% Due to symmetry, the compatibility constant q0 should
% evaluate to approximately zero for this particular section.
% -------------------------------------------------------------

qb = zeros(size(x));

for i = 2:length(x)

    ds_local = sqrt( ...
        (x(i)-x(i-1))^2 + ...
        (y(i)-y(i-1))^2 );

    % Trapezoidal integration of y*t*ds
    integral_ytds = 0;

    if i > 2
        y_prev = y(1:i-1);
        ds_prev = sqrt(diff(x(1:i-1)).^2 + ...
                       diff(y(1:i-1)).^2);

        integral_ytds = trapz( ...
            [0; cumsum(ds_prev)], ...
            y_prev * t );
    end

    % Add current segment contribution
    integral_ytds = integral_ytds + ...
        0.5*(y(i-1)+y(i))*t*ds_local;

    qb(i) = -(Vy/Ixx)*integral_ytds;
end


%% -------------------- CLOSED-CELL CONSTANT FLOW q0 --------------------
%
% Compatibility condition:
%
%        integral(q/t ds) = 0
%
% Therefore:
%
%        q0 = - integral(qb/t ds) / integral(ds/t)
%
% -------------------------------------------------------------

integral_qb = 0;
integral_inv_t = 0;

for i = 2:length(x)

    ds_local = sqrt( ...
        (x(i)-x(i-1))^2 + ...
        (y(i)-y(i-1))^2 );

    qb_avg = 0.5*(qb(i-1)+qb(i));

    integral_qb = integral_qb + ...
        qb_avg/t * ds_local;

    integral_inv_t = integral_inv_t + ...
        ds_local/t;
end

q0 = -integral_qb / integral_inv_t;


%% -------------------- TOTAL SHEAR FLOW --------------------

q = qb + q0;


%% -------------------- RESULTS --------------------

fprintf('\n============================================\n');
fprintf(' CLOSED RECTANGULAR SECTION SHEAR FLOW\n');
fprintf('============================================\n');

fprintf('Vy       = %.3f N\n', Vy);
fprintf('Width    = %.3f m\n', b);
fprintf('Height   = %.3f m\n', h);
fprintf('Thickness = %.3f mm\n', t*1000);
fprintf('Ixx      = %.6e m^4\n', Ixx);

fprintf('\nBasic shear flow q_b:\n');
fprintf('q_b,max  = %.3f N/m\n', max(qb));
fprintf('q_b,min  = %.3f N/m\n', min(qb));

fprintf('\nConstant closed-cell flow:\n');
fprintf('q0       = %.3f N/m\n', q0);

fprintf('\nTotal shear flow:\n');
fprintf('q_max    = %.3f N/m\n', max(q));
fprintf('q_min    = %.3f N/m\n', min(q));
fprintf('abs(q)_max = %.3f N/m\n', max(abs(q)));


%% -------------------- PLOT --------------------

figure('Color','w');

plot(x, y, 'k-', 'LineWidth', 1.5);
hold on;

% Scale factor for visualization
scale = 0.0001;

quiver(x(1:20:end), ...
       y(1:20:end), ...
       zeros(size(x(1:20:end))), ...
       q(1:20:end)*scale, ...
       0, ...
       'b');

axis equal;
grid on;

xlabel('x [m]');
ylabel('y [m]');
title('Shear Flow Distribution in Closed Rectangular Section');

legend('Section','Shear-flow direction');


%% -------------------- SHEAR FLOW VS PERIMETER --------------------

s = zeros(size(x));

for i = 2:length(x)
    s(i) = s(i-1) + ...
        sqrt((x(i)-x(i-1))^2 + ...
             (y(i)-y(i-1))^2);
end

figure('Color','w');

plot(s, q, 'LineWidth', 2);

grid on;

xlabel('Distance along perimeter, s [m]');
ylabel('Shear flow, q [N/m]');
title('Shear Flow Distribution Around Closed Section');

yline(0,'k--');