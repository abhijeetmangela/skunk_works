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
%   - Additional point loads can be entered in the POINT LOADS block below.
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
%
% ALL INPUTS LIVE IN THE "USER INPUTS" BLOCK RIGHT BELOW. Nothing else in
% the file needs to be touched to change a run: just edit the values there
% and run the script (no interactive input() prompts).

clc;
clear;
close all;

%% ========================================================================
% USER INPUTS -- edit the values below, then run the script
% ========================================================================

% ---------------- Bulkhead geometry ----------------
B = 0.05;                   % Bulkhead width B [m]
H = 0.05;                   % Bulkhead height H [m]
nPerMember = 400;          % Discretization points per member (>= 20)

% ---------------- Applied loading ----------------
Vy = -50;                   % Total vertical load Vy [N]
Fx = 0;                    % Total horizontal load Fx [N]
Mz = 10;                    % Applied moment Mz [N m]

% Members: 1 = left vertical, 2 = bottom horizontal, 3 = right vertical
loadMember = 2;             % Member that carries the distributed load
loadShapeType = 'uniform';  % 'uniform' or 'triangular'
                            % (normalized so the integral still equals Vy/Fx
                            % exactly, so changing Vy needs no other edits)

% ---------------- Point loads ----------------
% Each row is: [member, local_s, Fx, Fy, Mz], local_s measured from the
% start of that member. Leave as zeros(0,5) for no point loads. Example:
%   pointLoads = [2, 0.5, 0, -20, 0;   % 20 N downward at mid-span of member 2
%                 1, 0.3, 5,   0, 0];  % 5 N axial load on member 1
pointLoads = zeros(0,5);

% ---------------- Material: Aluminium (6061-T6) ----------------
E = 69e9;                   % Young's modulus E [Pa]
nu = 0.33;                  % Poisson ratio
G = E/(2*(1+nu));

% ---------------- Allowable stress ----------------
sigmaYield = 170e6;         % Yield strength [Pa] (Aluminium 6061-T6)
FOS = 8.5;                  % Factor of safety
sigmaAllow = sigmaYield/FOS;

% ---------------- Cross-section type ----------------
sectionType = 'channel';    % 'channel' or 'box'
sec_b = 0.05;               % Selected section width b [m]
sec_h = 0.08;               % Selected section height h [m]
t = 0.001;                  % Selected section thickness t [m]
shearAreaFactor = 0.8;      % As = shearAreaFactor*A

% ---------------- Automated sizing sweep ----------------
% NOTE ON RUNTIME: the sweep re-solves the full numerical Castigliano
% problem (a fresh solve over all nPerMember points) for every (b,h) pair,
% i.e. roughly nSize^2 full solves. With the defaults below (nSize=20,
% nPerMember=400) that is 400 solves and can take a noticeable amount of
% time in MATLAB, and be very slow in Octave. Reduce nSize and/or
% nPerMember first if you just want a quick check.
runSizing = 1;              % 1 = also sweep dimensions for a minimum-area
                             % section that satisfies sigmaAllow, 0 = skip
bMin = 0.01;  bMax = 0.08;   % Sizing sweep: width b range [m]
hMin = 0.03;  hMax = 0.12;   % Sizing sweep: height h range [m]
nSize = 20;                  % Number of values per dimension

% ---------------- Input validation (leave as-is) ----------------
nPerMember = max(20, round(nPerMember));
validateattributes(loadMember, {'numeric'}, {'scalar','integer','>=',1,'<=',3});

if isempty(sec_b) || isempty(sec_h) || isempty(t)
    error('Selected b, h and t are required.');
end

if ~ismember(sectionType, {'channel','box'})
    error('Section type must be ''channel'' or ''box''.');
end

if strcmp(sectionType,'box')
    fprintf('\nFor a box, valid dimensions require b > 2t and h > 2t.\n');
end

%% ========================================================================
% GEOMETRY AND LOAD DEFINITION
% ========================================================================

memberLength = [H/2, B, H/2];
memberAngle  = [-pi/2, 0, pi/2];

memberStart = [
     0,  H/2;
     0, -H/2;
     B, -H/2];

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

wDist = cell(3,1);
for m = 1:3
    wDist{m} = zeros(2,numel(sLocal{m}));
end

if Fx ~= 0 || Vy ~= 0
    shape = buildLoadShape(loadShapeType, sLocal{loadMember}, memberLength(loadMember));
    shape = normalizeShape(shape, sLocal{loadMember});
    wDist{loadMember} = [Fx*shape; Vy*shape];
end

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

            if t <= 0 || bb <= 0 || hh <= 0 || bb <= 2*t || hh <= 2*t
                continue;
            end

            p = sectionProperties(sectionType,bb,hh,t,shearAreaFactor);

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
    nMem = numel(Lm);

    [fixedN,fixedV,fixedM] = internalResponse(Lm, ang, r0, sLoc, xg, yg, wDist, pointLoads, Mz, nMem, [0;0;0]);

    influenceN = cell(3,1);
    influenceV = cell(3,1);
    influenceM = cell(3,1);

    for j = 1:3
        unit = zeros(3,1);
        unit(j) = 1;
        [influenceN{j}, influenceV{j}, influenceM{j}] = internalResponse(Lm, ang, r0, sLoc, xg, yg, wDist, pointLoads, Mz, nMem, unit);
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

    [N,V,M] = internalResponse(Lm, ang, r0, sLoc, xg, yg, wDist, pointLoads, Mz, nMem, u);

    sGlobal = cell(nMem,1);
    for m = 1:nMem
        sGlobal{m} = sum(Lm(1:m-1)) + sLoc{m};
    end

    % Horizontal concatenation gives one flat vector across all members and
    % all points, so max(abs(...)) below is a true scalar peak. (The original
    % cell2mat(N(:)) vertically stacked the members into a 3-by-nPerMember
    % matrix, so max() returned a 1-by-nPerMember row of per-column maxima
    % instead of the single overall maximum.)
    allN = [N{:}];
    allV = [V{:}];
    allM = [M{:}];

    sol.u = u;
    sol.N = N;
    sol.V = V;
    sol.M = M;
    sol.sGlobal = sGlobal;
    sol.Nmax = max(abs(allN));
    sol.Vmax = max(abs(allV));
    sol.Mmax = max(abs(allM));
end

function [Nout,Vout,Mout] = internalResponse(Lm, ang, r0, sLoc, xg, yg, wDist, pointLoads, Mz, nMem, uvec)
        Nout = cell(nMem,1);
        Vout = cell(nMem,1);
        Mout = cell(nMem,1);

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
