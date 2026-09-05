%% FLEXIBLE BULKHEAD DESIGN - NUMERICAL CASTIGLIANO METHOD
% Replacement for Appendix C.3 of the AS5100 report.
%
% The original report hard-coded the load polynomial, member lengths,
% material properties, cross-section dimensions, and several force terms.
% This version instead lets the user enter the physical inputs and performs
% the Castigliano calculation numerically.
%
% Model:
%   Three-member symmetric half-bulkhead: member 1 (vertical),
%   member 2 (horizontal), member 3 (vertical).
%   The three redundant/generalized forces at the cut are solved from
%   dU/d(F4,V4,M4) = 0.
%
% Load input:
%   - Vy = total vertical load applied to a selected member.
%   - The default distribution is uniform.
%   - Additional point loads can be entered in the POINT LOADS block.
%
% The code computes:
%   - Section properties A, I, As, ymax
%   - Internal axial force, shear force and bending moment distributions
%   - Castigliano redundants F4, V4, M4
%   - Maximum bending moment and stresses
%   - Required / selected section performance
%   - A dimension sweep for minimum channel or box dimensions
%
% Units: SI (m, N, Pa)

clc;
clear;
close all;

%% ========================================================================
% USER INPUTS
% ========================================================================

fprintf('\n===============================================================\n');
fprintf(' FLEXIBLE BULKHEAD DESIGN - NUMERICAL CASTIGLIANO METHOD\n');
fprintf('===============================================================\n\n');

% ---------------- Bulkhead geometry ----------------
B = input('Bulkhead width B [m]              = ');
H = input('Bulkhead height H [m]             = ');

% Number of points used on each member
nPerMember = input('Points per member [default 400] = ');
if isempty(nPerMember)
    nPerMember = 400;
end
nPerMember = max(20, round(nPerMember));

% ---------------- Applied loading ----------------
Vy = input('Total vertical load Vy [N]        = ');
Fx = input('Total horizontal load Fx [N]      = ');
Mz = input('Applied moment Mz [N m]            = ');

% Member receiving the distributed load
fprintf('\nMembers: 1 = left vertical, 2 = bottom horizontal, 3 = right vertical\n');
loadMember = input('Member for distributed load [1/2/3, default 2] = ');
if isempty(loadMember)
    loadMember = 2;
end
validateattributes(loadMember, {'numeric'}, {'scalar','integer','>=',1,'<=',3});

% Default distribution is uniform. The code normalizes it to the specified
% total resultant, so changing Vy does not require changing the load formula.
loadShapeType = lower(strtrim(input('Distributed load shape [uniform/triangular, default uniform]: ','s')));
if isempty(loadShapeType)
    loadShapeType = 'uniform';
end

% ---------------- Point loads ----------------
% Each row is: [member, local_s, Fx, Fy, Mz]
% local_s is measured from the start of that member.
usePointLoads = input('\nAdd point loads? [1=yes, 0=no, default 0] = ');
if isempty(usePointLoads)
    usePointLoads = 0;
end

pointLoads = zeros(0,5);
if usePointLoads == 1
    nPL = input('Number of point loads = ');
    pointLoads = zeros(nPL,5);
    fprintf('Enter each point load as: member, local_s[m], Fx[N], Fy[N], Mz[Nm]\n');
    for k = 1:nPL
        pointLoads(k,:) = input(sprintf('Point load %d = ',k));
        if numel(pointLoads(k,:)) ~= 5
            error('Each point load must contain 5 values.');
        end
    end
end

% ---------------- Material ----------------
E = input('Young''s modulus E [Pa]            = ');
nu = input('Poisson ratio nu                  = ');
if isempty(E)
    E = 70e9;
end
if isempty(nu)
    nu = 0.33;
end
G = E/(2*(1+nu));

% ---------------- Allowable stress ----------------
sigmaYield = input('Yield strength [Pa]               = ');
FOS = input('Factor of safety                  = ');
if isempty(sigmaYield)
    sigmaYield = 193e6;
end
if isempty(FOS)
    FOS = 8.5;
end
sigmaAllow = sigmaYield/FOS;

% ---------------- Cross-section type ----------------
sectionType = lower(strtrim(input('\nMember section [channel/box, default channel]: ','s')));
if isempty(sectionType)
    sectionType = 'channel';
end
if ~ismember(sectionType, {'channel','box'})
    error('Section type must be ''channel'' or ''box''.');
end

% Selected section dimensions for analysis
sec_b = input('Selected section width b [m]       = ');
sec_h = input('Selected section height h [m]     = ');
t = input('Selected section thickness t [m]   = ');

if isempty(sec_b) || isempty(sec_h) || isempty(t)
    error('Selected b, h and t are required.');
end

% Shear-area factor: As = factor*A. Leave editable rather than hiding it.
shearAreaFactor = input('Shear-area factor As/A [default 0.8] = ');
if isempty(shearAreaFactor)
    shearAreaFactor = 0.8;
end

% Optional automated sizing sweep
runSizing = input('\nRun automatic section-dimension sweep? [1=yes, 0=no, default 1] = ');
if isempty(runSizing)
    runSizing = 1;
end

if runSizing == 1
    bMin = input('Sizing sweep: minimum width b [m]    = ');
    bMax = input('Sizing sweep: maximum width b [m]    = ');
    hMin = input('Sizing sweep: minimum height h [m]   = ');
    hMax = input('Sizing sweep: maximum height h [m]   = ');
    nSize = input('Number of values per dimension [default 80] = ');
    if isempty(nSize)
        nSize = 80;
    end
else
    bMin = sec_b; bMax = sec_b;
    hMin = sec_h; hMax = sec_h;
    nSize = 1;
end

% Minimum practical thickness in the search (prevents invalid closed boxes)
if strcmp(sectionType,'box')
    fprintf('\nFor a box, valid dimensions require b > 2t and h > 2t.\n');
end

%% ========================================================================
% GEOMETRY AND LOAD DEFINITION
% ========================================================================

% Half-bulkhead path:
% member 1: down left web, member 2: across bottom, member 3: up right web.
memberLength = [H/2, B, H/2];
memberAngle  = [-pi/2, 0, pi/2];  % tangent angle in global x-y coordinates

% Start coordinates of each member
memberStart = [
     0,  H/2;
     0, -H/2;
     B, -H/2];

% Discretized contour
sLocal = cell(3,1);
xGlobal = cell(3,1);
yGlobal = cell(3,1);

for m = 1:3
    sLocal{m} = linspace(0, memberLength(m), nPerMember);
    c = cos(memberAngle(m));
    s = sin(memberAngle(m));
    xGlobal{m} = memberStart(m,1) + sLocal{m}*c;
    yGlobal{m} = memberStart(m,2) + sLocal{m}*s;
end

% Distributed load in global coordinates, force/length.
% We scale the requested resultant so that the numerical integral equals
% exactly the requested total load, independent of discretization.
wDist = cell(3,1);
for m = 1:3
    wDist{m} = zeros(2,numel(sLocal{m}));
end

if Fx ~= 0 || Vy ~= 0
    shape = buildLoadShape(loadShapeType, sLocal{loadMember}, memberLength(loadMember));
    shape = normalizeShape(shape, sLocal{loadMember});
    % Row 1 = global x-load/length, row 2 = global y-load/length.
    % The shape is normalized so its numerical integral gives the requested
    % total resultant exactly (to numerical integration accuracy).
    wDist{loadMember} = [Fx*shape; Vy*shape];
end

% Point loads are copied into a convenient structure.
pl = pointLoads;

%% ========================================================================
% CASTIGLIANO SOLUTION FOR SELECTED SECTION
% ========================================================================

props = sectionProperties(sectionType, sec_b, sec_h, t, shearAreaFactor);

fprintf('\nSelected section properties:\n');
fprintf('  A     = %.6e m^2\n', props.A);
fprintf('  I     = %.6e m^4\n', props.I);
fprintf('  As    = %.6e m^2\n', props.As);
fprintf('  ymax  = %.6e m\n', props.ymax);
fprintf('  G     = %.6e Pa\n', G);

sol = solveCastigliano(memberLength, memberAngle, memberStart, ...
    sLocal, xGlobal, yGlobal, wDist, pl, E, G, props, Fx, Vy, Mz);

F4 = sol.u(1);
V4 = sol.u(2);
M4 = sol.u(3);

fprintf('\nCastigliano solution:\n');
fprintf('  F4 = %+ .6f N\n', F4);
fprintf('  V4 = %+ .6f N\n', V4);
fprintf('  M4 = %+ .6f N m\n', M4);
fprintf('  Mmax = %.6f N m\n', sol.Mmax);

sigmaBendMax = sol.Mmax * props.ymax / props.I;
sigmaAxialMax = sol.Nmax / props.A;
tauMax = sol.Vmax / props.As;
sigmaCombinedConservative = sigmaBendMax + sigmaAxialMax;

fprintf('\nStress checks for selected section:\n');
fprintf('  Max |N|             = %.6f N\n', sol.Nmax);
fprintf('  Max |V|             = %.6f N\n', sol.Vmax);
fprintf('  Max bending stress  = %.6e Pa\n', sigmaBendMax);
fprintf('  Max axial stress    = %.6e Pa\n', sigmaAxialMax);
fprintf('  Max shear stress    = %.6e Pa\n', tauMax);
fprintf('  Allowable stress    = %.6e Pa\n', sigmaAllow);
fprintf('  Conservative N+M    = %.6e Pa\n', sigmaCombinedConservative);
if sigmaCombinedConservative <= sigmaAllow
    fprintf('  STATUS              = PASS\n');
else
    fprintf('  STATUS              = FAIL\n');
end

%% ========================================================================
% AUTOMATIC SECTION SIZING
% ========================================================================

if runSizing == 1
    bVals = linspace(bMin,bMax,nSize);
    hVals = linspace(hMin,hMax,nSize);

    bestFound = false;
    bestMassArea = inf;
    best = struct();

    for ib = 1:numel(bVals)
        for ih = 1:numel(hVals)
            bb = bVals(ib);
            hh = hVals(ih);

            % Geometry validity
            if t <= 0 || bb <= 0 || hh <= 0 || bb <= 2*t || hh <= 2*t
                continue;
            end

            p = sectionProperties(sectionType,bb,hh,t,shearAreaFactor);

            % Internal forces do not depend on cross-section stiffness only
            % under the present Castigliano model, so for each dimension the
            % correct compatibility solution must still be recomputed because
            % E*I, E*A and G*As change the redundants.
            testSol = solveCastigliano(memberLength, memberAngle, memberStart, ...
                sLocal, xGlobal, yGlobal, wDist, pl, E, G, p, Fx, Vy, Mz);

            sigM = testSol.Mmax*p.ymax/p.I;
            sigN = testSol.Nmax/p.A;
            sigCombined = sigM + sigN;

            if sigCombined <= sigmaAllow
                areaMetric = p.A;
                if areaMetric < bestMassArea
                    bestFound = true;
                    bestMassArea = areaMetric;
                    best.b = bb;
                    best.h = hh;
                    best.props = p;
                    best.sol = testSol;
                    best.sigM = sigM;
                    best.sigN = sigN;
                    best.sigCombined = sigCombined;
                end
            end
        end
    end

    if bestFound
        fprintf('\n===============================================================\n');
        fprintf(' AUTOMATIC MINIMUM-AREA SECTION FOUND\n');
        fprintf('===============================================================\n');
        fprintf('  Section type       = %s\n', sectionType);
        fprintf('  b                  = %.6f m (%.2f mm)\n', best.b, best.b*1e3);
        fprintf('  h                  = %.6f m (%.2f mm)\n', best.h, best.h*1e3);
        fprintf('  t                  = %.6f m (%.2f mm)\n', t, t*1e3);
        fprintf('  A                  = %.6e m^2\n', best.props.A);
        fprintf('  I                  = %.6e m^4\n', best.props.I);
        fprintf('  ymax               = %.6e m\n', best.props.ymax);
        fprintf('  Mmax               = %.6f N m\n', best.sol.Mmax);
        fprintf('  max bending stress = %.6e Pa\n', best.sigM);
        fprintf('  combined N+M       = %.6e Pa\n', best.sigCombined);
        fprintf('  allowable          = %.6e Pa\n', sigmaAllow);
    else
        fprintf('\nNo section in the requested sizing range satisfies the allowable stress.\n');
        fprintf('Increase bMax/hMax or reduce t/loading.\n');
    end
end

%% ========================================================================
% PLOTS
% ========================================================================

figure('Color','w','Name','Bulkhead internal force diagrams');

subplot(3,1,1); hold on;
for m = 1:3
    plot(sol.sGlobal{m}, sol.N{m}, 'LineWidth',1.5);
end
grid on;
xlabel('Contour coordinate s [m]');
ylabel('N [N]');
title('Axial force distribution');
legend('Member 1','Member 2','Member 3','Location','best');

subplot(3,1,2); hold on;
for m = 1:3
    plot(sol.sGlobal{m}, sol.V{m}, 'LineWidth',1.5);
end
grid on;
xlabel('Contour coordinate s [m]');
ylabel('V [N]');
title('Shear force distribution');

subplot(3,1,3); hold on;
for m = 1:3
    plot(sol.sGlobal{m}, sol.M{m}, 'LineWidth',1.5);
end
grid on;
xlabel('Contour coordinate s [m]');
ylabel('M [N m]');
title('Bending moment distribution');

figure('Color','w','Name','Bulkhead geometry');
hold on;
for m = 1:3
    plot(xGlobal{m}, yGlobal{m}, 'LineWidth',2);
end
scatter([0 B],[H/2 H/2],40,'filled');
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('Half-bulkhead structural idealization');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function shape = buildLoadShape(type, s, L)
    switch lower(type)
        case 'uniform'
            shape = ones(size(s));
        case 'triangular'
            % Zero at s=0 and maximum at s=L.
            shape = s/L;
        otherwise
            error('Unknown load shape. Use uniform or triangular.');
    end
end

function shape = normalizeShape(shape, s)
    area = trapz(s,shape);
    if abs(area) < eps
        error('Load shape has zero resultant.');
    end
    shape = shape/area;
end

function props = sectionProperties(type,b,h,t,shearAreaFactor)
    if t <= 0 || b <= 0 || h <= 0
        error('Section dimensions must be positive.');
    end

    switch lower(type)
        case 'channel'
            if h <= 2*t
                error('Channel height must be greater than 2t.');
            end
            % Non-overlapping rectangles: two flanges and one web.
            A1 = b*t;
            A2 = (h-2*t)*t;
            A = 2*A1 + A2;
            ybar = h/2;
            d = h/2 - t/2;
            I = 2*(b*t^3/12 + A1*d^2) + t*(h-2*t)^3/12;
            ymax = h/2;

        case 'box'
            if b <= 2*t || h <= 2*t
                error('Box dimensions must satisfy b>2t and h>2t.');
            end
            A = b*h - (b-2*t)*(h-2*t);
            I = (b*h^3 - (b-2*t)*(h-2*t)^3)/12;
            ymax = h/2;

        otherwise
            error('Unknown section type.');
    end

    As = shearAreaFactor*A;

    props.A = A;
    props.I = I;
    props.As = As;
    props.ymax = ymax;
end

function sol = solveCastigliano(Lm, ang, r0, sLoc, xg, yg, wDist, pointLoads, E, G, props, Fx, Vy, Mz)
    %#ok<INUSD>
    nMem = numel(Lm);

    % Generalized unknowns at the symmetry cut:
    % u = [F4; V4; M4].
    % Castigliano compatibility gives K*u = -r.
    [fixedN,fixedV,fixedM] = internalResponse([0;0;0]);

    influenceN = cell(3,1);
    influenceV = cell(3,1);
    influenceM = cell(3,1);

    for j = 1:3
        unit = zeros(3,1);
        unit(j) = 1;
        [influenceN{j}, influenceV{j}, influenceM{j}] = internalResponse(unit);
    end

    K = zeros(3,3);
    rhs = zeros(3,1);

    for m = 1:nMem
        s = sLoc{m};
        for j = 1:3
            nj = influenceN{j}{m};
            vj = influenceV{j}{m};
            mj = influenceM{j}{m};

            rhs(j) = rhs(j) + trapz(s, ...
                fixedN{m}.*nj/(E*props.A) + ...
                fixedV{m}.*vj/(G*props.As) + ...
                fixedM{m}.*mj/(E*props.I));

            for k = 1:3
                nk = influenceN{k}{m};
                vk = influenceV{k}{m};
                mk = influenceM{k}{m};
                K(j,k) = K(j,k) + trapz(s, ...
                    nj.*nk/(E*props.A) + ...
                    vj.*vk/(G*props.As) + ...
                    mj.*mk/(E*props.I));
            end
        end
    end

    if rcond(K) < 1e-12
        warning('Castigliano matrix is poorly conditioned; using pseudoinverse.');
        u = -pinv(K)*rhs;
    else
        u = -K\rhs;
    end

    [N,V,M] = internalResponse(u);

    sGlobal = cell(nMem,1);
    for m = 1:nMem
        sGlobal{m} = sum(Lm(1:m-1)) + sLoc{m};
    end

    allN = cell2mat(N(:));
    allV = cell2mat(V(:));
    allM = cell2mat(M(:));

    sol.u = u;
    sol.N = N;
    sol.V = V;
    sol.M = M;
    sol.sGlobal = sGlobal;
    sol.Nmax = max(abs(allN));
    sol.Vmax = max(abs(allV));
    sol.Mmax = max(abs(allM));

    function [Nout,Vout,Mout] = internalResponse(uvec)
        Nout = cell(nMem,1);
        Vout = cell(nMem,1);
        Mout = cell(nMem,1);

        % Applied Mz is treated as a moment at the starting cut.
        Rstart = uvec(1:2);
        Mstart = uvec(3) + Mz;

        for m = 1:nMem
            sl = sLoc{m};
            x = xg{m};
            y = yg{m};
            tHat = [cos(ang(m)); sin(ang(m))];
            nHat = [-sin(ang(m)); cos(ang(m))];

            Nvec = zeros(size(sl));
            Vvec = zeros(size(sl));
            Mvec = zeros(size(sl));

            for ii = 1:numel(sl)
                xCut = x(ii);
                yCut = y(ii);
                R = Rstart;
                MM = Mstart;

                % Distributed loads on all members traversed up to this cut.
                for p = 1:m
                    sp = sLoc{p};
                    xp = xg{p};
                    yp = yg{p};
                    wp = wDist{p};

                    if p < m
                        idx = 1:numel(sp);
                    else
                        idx = 1:ii;
                    end

                    if numel(idx) >= 2
                        spp = sp(idx);
                        xpp = xp(idx);
                        ypp = yp(idx);
                        fxp = wp(1,idx);
                        fyp = wp(2,idx);

                        R = R + [trapz(spp,fxp); trapz(spp,fyp)];
                        MM = MM + trapz(spp, ...
                            (xpp-xCut).*fyp - (ypp-yCut).*fxp);
                    end
                end

                % Point loads before/equal to the current cut.
                if ~isempty(pointLoads)
                    for r = 1:size(pointLoads,1)
                        mp = pointLoads(r,1);
                        sp = pointLoads(r,2);
                        if mp < m || (mp == m && sp <= sl(ii) + 1e-12)
                            FxP = pointLoads(r,3);
                            FyP = pointLoads(r,4);
                            MP = pointLoads(r,5);
                            xP = r0(mp,1) + sp*cos(ang(mp));
                            yP = r0(mp,2) + sp*sin(ang(mp));

                            R = R + [FxP;FyP];
                            MM = MM + (xP-xCut)*FyP - ...
                                  (yP-yCut)*FxP + MP;
                        end
                    end
                end

                % Equilibrium of the cut free body.
                Rcut = -R;
                Mcut = -MM;

                Nvec(ii) = dot(Rcut,tHat);
                Vvec(ii) = dot(Rcut,nHat);
                Mvec(ii) = Mcut;
            end

            Nout{m} = Nvec;
            Vout{m} = Vvec;
            Mout{m} = Mvec;
        end
    end
end
