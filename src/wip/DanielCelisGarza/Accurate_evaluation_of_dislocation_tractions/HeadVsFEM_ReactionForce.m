%%% Script to compare Head's analytical solutions to those obtained via FEM using analytic and numeric tractions.
close all
clear all
% TODO #45
amag = 3.18e-4;

CRYSTAL_STRUCTURE = 'bcc';
% FEM PARAMETERS
%Cantilever
simTime = 0;
use_gpu = 0;
para_scheme = 0;
n_threads = 0;

%MU = 160E9; %160GPa
MU = 1;
NU = 0.28;

% Lattice constants
a = 10;
bVec = [[1 0 0]; [0 1 0]];

planes = [1; 2; 3; 4; 5; 6];
dx = 1000;
dy = 1000;
dz = 1000;

figCounter = 0;
cntr = 0;

for j = 20
    mx = j;

    gridSize = mx;
    x = linspace(0, dx, gridSize);
    y = linspace(0, dy, gridSize);
    z = linspace(0.5 * dz, 0.5 * dz, gridSize);
    [X, Y] = meshgrid(x, y);
    Z = meshgrid(z);

    clear x y z;

    loading = 1;
    vertices = [0, 0, 0; ...
                dx, 0, 0; ...
                0, dy, 0; ...
                dx, dy, 0; ...
                0, 0, dz; ...
                dx, 0, dz; ...
                0, dy, dz; ...
                dx, dy, dz];

    plim = 12 / amag; %12microns
    [xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, Sbot, ...
            Sfront, Sback, Smixed, gammat, gammau, gammaMixed, fixedDofs, freeDofs, ...
            w, h, d, my, mz, mel] = STATIC_finiteElement3D(dx, dy, dz, mx, MU, NU, loading);

    gamma_dln = [gammat(:, 1)];

    % Set surface node labels for surface node extraction.
    n_nodes = 4;
    surf_node_util = zeros(n_nodes + 2, 6);
    xy = mx * my;
    xz = mx * mz;
    yz = my * mz;

    f_hat = zeros(3 * mno, 1);
    [x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx; my; mz], ...
        planes, 4);
    [f_dln_node, f_dln_se, ...
            f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
    tolerance = dx / 10^6;

    figCounter = figCounter + 1;
    figure(figCounter)
    clf; hold on; view(3)
    xlabel('x'); ylabel('y'); zlabel('z')

    %     plot3(xnodes(Stop(:, 1), 1), xnodes(Stop(:, 1), 2), xnodes(Stop(:, 1), 3), 'b*')
    %     plot3(xnodes(Sbot(:, 1), 1), xnodes(Sbot(:, 1), 2), xnodes(Sbot(:, 1), 3), 'b*')
    %     plot3(xnodes(Sright(:, 1), 1), xnodes(Sright(:, 1), 2), xnodes(Sright(:, 1), 3), 'r.')
    plot3(xnodes(Sleft(:, 1), 1), xnodes(Sleft(:, 1), 2), xnodes(Sleft(:, 1), 3), 'r.')
    %     plot3(xnodes(Sfront(:, 1), 1), xnodes(Sfront(:, 1), 2), xnodes(Sfront(:, 1), 3), 'k*')
    %     plot3(xnodes(Sback(:, 1), 1), xnodes(Sback(:, 1), 2), xnodes(Sback(:, 1), 3), 'k*')
    %     plot3(xnodes(Smixed(:, 1), 1), xnodes(Smixed(:, 1), 2), xnodes(Smixed(:, 1), 3), 'g*')
    plot3(xnodes(gammat(:, 1), 1), xnodes(gammat(:, 1), 2), xnodes(gammat(:, 1), 3), 'ko')
    axis('equal')
    hold off

    len = 100;
    xcoord = linspace(0, dx, j);
    ycoord = linspace(0, dy, j);
    % TODO #44
    xcoord = xcoord(2) / 2; % middle of first element.
    ycoord = (ycoord(floor(j / 2)) + ycoord(floor(j / 2) + 1)) / 2; % middle of the domain
    x = linspace(xcoord, xcoord, len);
    y = linspace(ycoord, ycoord, len);
    z = linspace(0, dz, len);
    x1 = x(1);
    y1 = y(1);
    t = [0 0 1];
    n = [1 0 0];
    rn = zeros(len, 3);
    rn(:, 1) = x;
    rn(:, 2) = y;
    rn(:, 3) = z;
    links = zeros(len - 1, 8);
    hold on
    plot3(rn(:, 1), rn(:, 2), rn(:, 3), 'r.')
    plot3(X, Y, Z, 'k.')
    hold off
end

doSave = true;

close all;

if doSave
    save('mesh')
end

%%
% TODO #34 Couplers to use features in to v2.0
addpath 'D:\DPhil\OneDrive - Nexus365\EasyDD\src'

for k = 1:2
    close all
    cntr = cntr + 1;
    b = bVec(k, :);

    for i = 1:len - 1
        links(i, :) = [i, i + 1, b, n];
    end

    [uhat, fend, Ubar, fan] = STATIC_analytic_FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
        0, 0, gammaMixed, fixedDofs, freeDofs, dx, simTime, ...
        gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node, ...
        f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);

    [uhat2, fend2, Ubar2, fnum] = STATIC_FEMcoupler(rn, links, 0, a, MU, NU, xnodes, mno, kg, L, U, ...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);

    sigmaA = hatStressSurf(uhat, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
    sxxA = squeeze(sigmaA(1, 1, :, :));
    syyA = squeeze(sigmaA(2, 2, :, :));
    sxyA = squeeze(sigmaA(1, 2, :, :));
    sigmaN = hatStressSurf(uhat2, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
    sxxN = squeeze(sigmaN(1, 1, :, :));
    syyN = squeeze(sigmaN(2, 2, :, :));
    sxyN = squeeze(sigmaN(1, 2, :, :));

    segments = constructsegmentlist(rn, links, doSave);
    p1 = [segments(:, 6) segments(:, 7) segments(:, 8)];
    p2 = [segments(:, 9) segments(:, 10) segments(:, 11)];

    sigmaFP = FieldPointStressSurf(X, Y, Z, p1, p2, b, a, MU, NU);
    sxxFP = squeeze(sigmaFP(1, 1, :, :));
    syyFP = squeeze(sigmaFP(2, 2, :, :));
    sxyFP = squeeze(sigmaFP(1, 2, :, :));

    x = linspace(0, dx, gridSize);
    y = linspace(0, dy, gridSize);
    b = 1; %sqrt(3) / 2;
    [X, Y] = meshgrid(x, y);

    if k == 1
        [txx, tyy, txy] = imageStressAnalyticEdgePerp(MU, b, NU, X, Y, x1, y1);
        [txxFP, tyyFP, txyFP] = FPStressAnalyticEdgePerp(MU, b, NU, X, Y, x1, y1);
        orientationB = 'Eperp';

        sxxAperp = sxxA;
        syyAperp = syyA;
        sxyAperp = sxyA;

        sxxNperp = sxxN;
        syyNperp = syyN;
        sxyNperp = sxyN;

        txxPerp = txx;
        tyyPerp = tyy;
        txyPerp = txy;
    else
        [txx, tyy, txy] = imageStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
        [txxFP, tyyFP, txyFP] = FPStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
        orientationB = 'Epar';

        sxxApar = sxxA;
        syyApar = syyA;
        sxyApar = sxyA;

        sxxNpar = sxxN;
        syyNpar = syyN;
        sxyNpar = sxyN;

        txxPar = txx;
        tyyPar = tyy;
        txyPar = txy;
    end
    
    % Edge
    % Image stresses
    symbol = '\hat{\sigma}';
    % Head
    [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txx, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyy, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    % FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA, x1, y1, orientationB, symbol, 'xx', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyA, x1, y1, orientationB, symbol, 'yy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyA, x1, y1, orientationB, symbol, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
    % FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN, x1, y1, orientationB, symbol, 'xx', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyN, x1, y1, orientationB, symbol, 'yy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyN, x1, y1, orientationB, symbol, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
%     % Rel Err FEM + Analytic tractions
%     plotCountourfSigmaHat(X, Y, sxxA ./ txx - 1, x1, y1, orientationB, symbol, 'xx', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyA ./ tyy - 1, x1, y1, orientationB, symbol, 'yy', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyA ./ txy - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     % Rel Err FEM + Numeric tractions
%     plotCountourfSigmaHat(X, Y, sxxN ./ txx - 1, x1, y1, orientationB, symbol, 'xx', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyN ./ tyy - 1, x1, y1, orientationB, symbol, 'yy', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyN ./ txy - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
%     % Abs Err FEM + Analytic tractions
%     plotCountourfSigmaHat(X, Y, sxxA - txx, x1, y1, orientationB, symbol, 'xx', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, syyA - tyy, x1, y1, orientationB, symbol, 'yy', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyA - txy, x1, y1, orientationB, symbol, 'xy', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
%     % Abs Err FEM + Numeric tractions
%     plotCountourfSigmaHat(X, Y, sxxN - txx, x1, y1, orientationB, symbol, 'xx', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, syyN - tyy, x1, y1, orientationB, symbol, 'yy', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
%     plotCountourfSigmaHat(X, Y, sxyN - txy, x1, y1, orientationB, symbol, 'xy', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)

%     % Real stresses
%     symbol = '\tilde{\sigma}';
%     % Head
%     orientationB = strcat(orientationB, '_real');
%     [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txxFP, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyyFP, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txyFP, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     % FEM + Analytic Fieldpoints
%     plotCountourfSigmaHat(X, Y, sxxFP, x1, y1, orientationB, symbol, 'xx', 'FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
%     plotCountourfSigmaHat(X, Y, syyFP, x1, y1, orientationB, symbol, 'yy', 'FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
%     plotCountourfSigmaHat(X, Y, sxyFP, x1, y1, orientationB, symbol, 'xy', 'FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
%     % Rel Err
%     plotCountourfSigmaHat(X, Y, sxxFP ./ txxFP - 1, x1, y1, orientationB, symbol, 'xx', '$\eta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
%     plotCountourfSigmaHat(X, Y, syyFP ./ syyFP - 1, x1, y1, orientationB, symbol, 'yy', '$\eta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
%     plotCountourfSigmaHat(X, Y, sxyFP ./ sxyFP - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
%     % Abs err
%     plotCountourfSigmaHat(X, Y, sxxFP - txxFP, x1, y1, orientationB, symbol, 'xx', '$\Delta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
%     plotCountourfSigmaHat(X, Y, syyFP - syyFP, x1, y1, orientationB, symbol, 'yy', '$\Delta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
%     plotCountourfSigmaHat(X, Y, sxyFP - sxyFP, x1, y1, orientationB, symbol, 'xy', '$\Delta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% 
%     % Total stresses
%     symbol = '\sigma';
%     % Head
%     orientationB = strcat(orientationB, '_total');
%     [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txx + txxFP, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyy + tyyFP, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy + txyFP, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
%     % FEM + Analytic tractions
%     plotCountourfSigmaHat(X, Y, sxxA + sxxFP, x1, y1, orientationB, symbol, 'xx', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
%     plotCountourfSigmaHat(X, Y, syyA + syyFP, x1, y1, orientationB, symbol, 'yy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
%     plotCountourfSigmaHat(X, Y, sxyA + sxyFP, x1, y1, orientationB, symbol, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
%     % FEM + Numeric tractions
%     plotCountourfSigmaHat(X, Y, sxxN + sxxFP, x1, y1, orientationB, symbol, 'xx', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
%     plotCountourfSigmaHat(X, Y, syyN + syyFP, x1, y1, orientationB, symbol, 'yy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
%     plotCountourfSigmaHat(X, Y, sxyN + sxyFP, x1, y1, orientationB, symbol, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
%     % Rel Err FEM + Analytic tractions
%     txxT = txx + txxFP;
%     tyyT = tyy + tyyFP;
%     txyT = txy + txyFP;
%     txxT(abs(txxT) < eps) = eps;
%     tyyT(abs(tyyT) < eps) = eps;
%     txyT(abs(txyT) < eps) = eps;
%     plotCountourfSigmaHat(X, Y, (sxxA + sxxFP) ./ (txxT) - 1, x1, y1, orientationB, symbol, 'xx', '$\eta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (syyA + syyFP) ./ (tyyT) - 1, x1, y1, orientationB, symbol, 'yy', '$\eta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (sxyA + sxyFP) ./ (txyT) - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     % Rel Err FEM + Numeric tractions
%     plotCountourfSigmaHat(X, Y, (sxxN + sxxFP) ./ (txxT) - 1, x1, y1, orientationB, symbol, 'xx', '$\eta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (syyN + syyFP) ./ (tyyT) - 1, x1, y1, orientationB, symbol, 'yy', '$\eta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (sxyN + sxyFP) ./ (txyT) - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     % Abs Err FEM + Analytic tractions
%     plotCountourfSigmaHat(X, Y, (sxxA + sxxFP) - (txxT), x1, y1, orientationB, symbol, 'xx', '$\Delta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (syyA + syyFP) - (tyyT), x1, y1, orientationB, symbol, 'yy', '$\Delta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (sxyA + sxyFP) - (txyT), x1, y1, orientationB, symbol, 'xy', '$\Delta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     % Abs Err FEM + Numeric tractions
%     plotCountourfSigmaHat(X, Y, (sxxN + sxxFP) - (txxT), x1, y1, orientationB, symbol, 'xx', '$\Delta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (syyN + syyFP) - (tyyT), x1, y1, orientationB, symbol, 'yy', '$\Delta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
%     plotCountourfSigmaHat(X, Y, (sxyN + sxyFP) - (txyT), x1, y1, orientationB, symbol, 'xy', '$\Delta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)

end

%%

% Screw
close all
b = [0 0 1];
orientationB = 'screw';

for i = 1:len - 1
    links(i, :) = [i, i + 1, b, n];
end

[uhat, fend, Ubar, fan] = STATIC_analytic_FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
    0, 0, gammaMixed, fixedDofs, freeDofs, dx, simTime, ...
    gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node, ...
    f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);
[uhat2, fend2, Ubar2, fnum] = STATIC_FEMcoupler(rn, links, 0, a, MU, NU, xnodes, mno, kg, L, U, ...
    gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);

segments = constructsegmentlist(rn, links, doSave);
p1 = [segments(:, 6) segments(:, 7) segments(:, 8)];
p2 = [segments(:, 9) segments(:, 10) segments(:, 11)];

sigmaFP = FieldPointStressSurf(X, Y, Z, p1, p2, b, a, MU, NU);
sxzFP = squeeze(sigmaFP(1, 3, :, :));
syzFP = squeeze(sigmaFP(2, 3, :, :));

sigmaA = hatStressSurf(uhat, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
sigmaN = hatStressSurf(uhat2, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
sxzA = squeeze(sigmaA(1, 3, :, :));
syzA = squeeze(sigmaA(2, 3, :, :));
sxzN = squeeze(sigmaN(1, 3, :, :));
syzN = squeeze(sigmaN(2, 3, :, :));

[txz, tyz] = imageStressAnalyticScrew(MU, 1, NU, X, Y, x1, y1);
[txzFP, tyzFP] = FPStressAnalyticScrew(MU, 1, NU, X, Y, x1, y1);
txzT = txz + txzFP;
txyT = txy + txyFP;
txzT(abs(txzT) < eps) = eps;
txyT(abs(txyT) < eps) = eps;

% Image stress
symbol = '\hat{\sigma}';
% Head
[~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txz, x1, y1, orientationB, symbol, 'xz', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
[~, meanvalyz, stddevyz] = plotCountourfSigmaHat(X, Y, tyz, x1, y1, orientationB, symbol, 'yz', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
% Analytic
plotCountourfSigmaHat(X, Y, sxzA, x1, y1, orientationB, symbol, 'xz', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
plotCountourfSigmaHat(X, Y, syzA, x1, y1, orientationB, symbol, 'yz', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyz, stddevyz)
% Numeric
plotCountourfSigmaHat(X, Y, sxzN, x1, y1, orientationB, symbol, 'xz', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
plotCountourfSigmaHat(X, Y, syzN, x1, y1, orientationB, symbol, 'yz', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyz, stddevyz)
% % Rel err FEM + analytic
% plotCountourfSigmaHat(X, Y, sxzA ./ txz - 1, x1, y1, orientationB, symbol, 'xz', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
% plotCountourfSigmaHat(X, Y, syzA ./ tyz - 1, x1, y1, orientationB, symbol, 'yz', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
% % Rel err FEM + numeric
% plotCountourfSigmaHat(X, Y, sxzN ./ txz - 1, x1, y1, orientationB, symbol, 'xz', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
% plotCountourfSigmaHat(X, Y, syzN ./ tyz - 1, x1, y1, orientationB, symbol, 'yz', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 30, doSave)
% % Abs err FEM + analytic
% plotCountourfSigmaHat(X, Y, sxzA - txz, x1, y1, orientationB, symbol, 'xz', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
% plotCountourfSigmaHat(X, Y, syzA - tyz, x1, y1, orientationB, symbol, 'yz', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
% % Abs err FEM + numeric
% plotCountourfSigmaHat(X, Y, sxzN - txz, x1, y1, orientationB, symbol, 'xz', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)
% plotCountourfSigmaHat(X, Y, syzN - tyz, x1, y1, orientationB, symbol, 'yz', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 30, doSave)

% % Real stresses
% symbol = '\tilde{\sigma}';
% % Head
% orientationB = strcat(orientationB, '_real');
% [~, meanvalxz, stddevxz] = plotCountourfSigmaHat(X, Y, txzFP, x1, y1, orientationB, symbol, 'xz', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
% [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txyFP, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
% % FEM + Analytic tractions
% plotCountourfSigmaHat(X, Y, sxzFP, x1, y1, orientationB, symbol, 'xz', 'FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxz, stddevxz)
% plotCountourfSigmaHat(X, Y, sxyFP, x1, y1, orientationB, symbol, 'xy', 'FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% % Rel Err
% plotCountourfSigmaHat(X, Y, sxzFP ./ txzFP - 1, x1, y1, orientationB, symbol, 'xz', '$\eta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxz, stddevxz)
% plotCountourfSigmaHat(X, Y, sxyFP ./ txyFP - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% % Abs Err
% plotCountourfSigmaHat(X, Y, sxzFP - txzFP, x1, y1, orientationB, symbol, 'xz', '$\Delta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxz, stddevxz)
% plotCountourfSigmaHat(X, Y, sxyFP - txyFP, x1, y1, orientationB, symbol, 'xy', '$\Delta$FP', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% 
% % Total stresses
% symbol = '\sigma';
% % Head
% orientationB = strcat(orientationB, '_total');
% [~, meanvalxz, stddevxz] = plotCountourfSigmaHat(X, Y, txz + txzFP, x1, y1, orientationB, symbol, 'xz', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
% [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy + txyFP, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
% % FEM + Analytic tractions
% plotCountourfSigmaHat(X, Y, sxzA + sxzFP, x1, y1, orientationB, symbol, 'xz', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxz, stddevxz)
% plotCountourfSigmaHat(X, Y, sxyA + sxyFP, x1, y1, orientationB, symbol, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% % FEM + Numeric tractions
% plotCountourfSigmaHat(X, Y, sxzN + sxzFP, x1, y1, orientationB, symbol, 'xz', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxz, stddevxz)
% plotCountourfSigmaHat(X, Y, sxyN + sxyFP, x1, y1, orientationB, symbol, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
% % Rel Err FEM + Analytic tractions
% plotCountourfSigmaHat(X, Y, (sxzA + sxzFP) ./ (txzT) - 1, x1, y1, orientationB, symbol, 'xz', '$\eta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, (sxyA + sxyFP) ./ (txyT) - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% % Rel Err FEM + Numeric tractions
% plotCountourfSigmaHat(X, Y, (sxzN + sxzFP) ./ (txzT) - 1, x1, y1, orientationB, symbol, 'xz', '$\eta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, (sxyN + sxyFP) ./ (txyT) - 1, x1, y1, orientationB, symbol, 'xy', '$\eta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% % Abs Err FEM + Analytic tractions
% plotCountourfSigmaHat(X, Y, (sxzA + sxzFP) - (txzT), x1, y1, orientationB, symbol, 'xz', '$\Delta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, (sxyA + sxyFP) - (txyT), x1, y1, orientationB, symbol, 'xy', '$\Delta$A', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% % Abs Err FEM + Numeric tractions
% plotCountourfSigmaHat(X, Y, (sxzN + sxzFP) - (txzT), x1, y1, orientationB, symbol, 'xz', '$\Delta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, (sxyN + sxyFP) - (txyT), x1, y1, orientationB, symbol, 'xy', '$\Delta$N', 'x,~b', 'y,~b', '$\mu$', 30, doSave)

% Line plots
% Edge perp
close all
i = 2;
symbol = '\hat{\sigma}';
orientationB = 'Eperp';
linePlot(sxxAperp(:, i), sxxNperp(:, i), txxPerp(:, i), orientationB, symbol, 'xx', 'Grid Point', '$\mu$', 30, doSave)
linePlot(syyAperp(:, i), syyNperp(:, i), tyyPerp(:, i), orientationB, symbol, 'yy', 'Grid Point', '$\mu$', 30, doSave)
linePlot(sxyAperp(:, i), sxyNperp(:, i), txyPerp(:, i), orientationB, symbol, 'xy', 'Grid Point', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, txxPerp, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% x = linspace(xcoord, xcoord, j);
% y = linspace(0, dy, j);
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off
% plotCountourfSigmaHat(X, Y, tyyPerp, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off
% plotCountourfSigmaHat(X, Y, txyPerp, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off

% Edge par
i = 2;
orientationB = 'Epar';
linePlot(sxxApar(:, i), sxxNpar(:, i), txxPar(:, i), orientationB, symbol, 'xx', 'Grid Point', '$\mu$', 30, doSave)
linePlot(syyApar(:, i), syyNpar(:, i), tyyPar(:, i), orientationB, symbol, 'yy', 'Grid Point', '$\mu$', 30, doSave)
linePlot(sxyApar(:, i), sxyNpar(:, i), txyPar(:, i), orientationB, symbol, 'xy', 'Grid Point', '$\mu$', 30, doSave)
% plotCountourfSigmaHat(X, Y, txxPar, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off
% plotCountourfSigmaHat(X, Y, tyyPar, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off
% plotCountourfSigmaHat(X, Y, txyPar, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off

% Screw
i = 2;
orientationB = 'screw';
linePlot(sxzA(:, i), sxzN(:, i), txz(:, i), orientationB, symbol, 'xz', 'Grid Point', '$\mu$', 30, doSave)
linePlot(syzA(:, i), syzN(:, i), tyz(:, i), orientationB, symbol, 'yz', 'Grid Point', '$\mu$', 30, doSave)
plotCountourfSigmaHat(X, Y, txz, x1, y1, orientationB, symbol, 'xz', '', 'x,~b', 'y,~b', '$\mu$', 30, false)
hold on
plot(linspace(x(2), x(2), j), y, 'LineWidth', 2, 'LineStyle', '--')
hold off
set(gcf(), 'Units', 'Inches');
pos = get(gcf(), 'Position');
set(gcf(), 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3)+pos(3)*0.1, pos(4)])
print(gcf(), sprintf('./paper/images/contourLine.pdf'), '-dpdf', '-r0')
% plotCountourfSigmaHat(X, Y, tyz, x1, y1, orientationB, symbol, 'yz', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave)
% hold on
% plot(x, y, 'LineWidth', 2, 'LineStyle', '--')
% hold off
close all
%%
plotCountourfSigmaHat(X, Y, txz, x1, y1, orientationB, symbol, 'xz', '', 'x,~b', 'y,~b', '$\mu$', 30, false)
linePlot(sxzA(:, i), sxzN(:, i), txz(:, i), orientationB, symbol, 'xz', 'Grid Point', '$\mu$', 30, doSave)
%%
addpath 'D:\DPhil\OneDrive - Nexus365\EasyDD\src'

for node = [2, 5, 11, 18]
    close all
    len = 100;
    xcoord = linspace(0, dx, j);
    ycoord = linspace(0, dy, j);
    % TODO #44
    xcoord = (xcoord(node-1) + xcoord(node)) / 2; % middle of first element.
    ycoord = (ycoord(floor(j / 2)) + ycoord(floor(j / 2) + 1)) / 2; % middle of the domain
    x = linspace(xcoord, xcoord, len);
    y = linspace(ycoord, ycoord, len);
    z = linspace(0, dz, len);
    x1 = x(1);
    y1 = y(1);
    t = [0 0 1];
    n = [1 0 0];
    rn = zeros(len, 3);
    rn(:, 1) = x;
    rn(:, 2) = y;
    rn(:, 3) = z;
    links = zeros(len - 1, 8);
    hold on
    plot3(rn(:, 1), rn(:, 2), rn(:, 3), 'r.')
    plot3(X, Y, Z, 'k.')
    hold off
    
    b = bVec(2, :);
    
    for i = 1:len - 1
        links(i, :) = [i, i + 1, b, n];
    end
    
    [uhat, fend, Ubar, fan] = STATIC_analytic_FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
        0, 0, gammaMixed, fixedDofs, freeDofs, dx, simTime, ...
        gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node, ...
        f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);
    
    [uhat2, fend2, Ubar2, fnum] = STATIC_FEMcoupler(rn, links, 0, a, MU, NU, xnodes, mno, kg, L, U, ...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);
    
    sigmaA = hatStressSurf(uhat, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
    sxxA = squeeze(sigmaA(1, 1, :, :));
    syyA = squeeze(sigmaA(2, 2, :, :));
    sxyA = squeeze(sigmaA(1, 2, :, :));
    sigmaN = hatStressSurf(uhat2, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
    sxxN = squeeze(sigmaN(1, 1, :, :));
    syyN = squeeze(sigmaN(2, 2, :, :));
    sxyN = squeeze(sigmaN(1, 2, :, :));
    
    segments = constructsegmentlist(rn, links, doSave);
    p1 = [segments(:, 6) segments(:, 7) segments(:, 8)];
    p2 = [segments(:, 9) segments(:, 10) segments(:, 11)];
    
    sigmaFP = FieldPointStressSurf(X, Y, Z, p1, p2, b, a, MU, NU);
    sxxFP = squeeze(sigmaFP(1, 1, :, :));
    syyFP = squeeze(sigmaFP(2, 2, :, :));
    sxyFP = squeeze(sigmaFP(1, 2, :, :));
    
    x = linspace(0, dx, gridSize);
    y = linspace(0, dy, gridSize);
    b = 1; %sqrt(3) / 2;
    [X, Y] = meshgrid(x, y);
    
    [txx, tyy, txy] = imageStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
    [txxFP, tyyFP, txyFP] = FPStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
    orientationB = 'EparQuarter';
    
    sxxApar = sxxA;
    syyApar = syyA;
    sxyApar = sxyA;
    
    sxxNpar = sxxN;
    syyNpar = syyN;
    sxyNpar = sxyN;
    
    txxPar = txx;
    tyyPar = tyy;
    txyPar = txy;
    
    % Edge
    % Image stresses
    symbol = '\hat{\sigma}';
    % Head
    [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txx, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyy, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    % FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA, x1, y1, orientationB, symbol, 'xx', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyA, x1, y1, orientationB, symbol, 'yy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyA, x1, y1, orientationB, symbol, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
    % FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN, x1, y1, orientationB, symbol, 'xx', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyN, x1, y1, orientationB, symbol, 'yy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyN, x1, y1, orientationB, symbol, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
    
    orientationB = sprintf('Epar%d', node);
    linePlot(sxxApar(:, node), sxxNpar(:, node), txxPar(:, node), orientationB, symbol, 'xx', 'Grid Point', '$\mu$', 30, doSave)
    linePlot(syyApar(:, node), syyNpar(:, node), tyyPar(:, node), orientationB, symbol, 'yy', 'Grid Point', '$\mu$', 30, doSave)
    linePlot(sxyApar(:, node), sxyNpar(:, node), txyPar(:, node), orientationB, symbol, 'xy', 'Grid Point', '$\mu$', 30, doSave)
    
    
    % Edge
    % Total stresses
    symbol = '\sigma';
    % Head
    [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txx+sxxFP, x1, y1, orientationB, symbol, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyy+syyFP, x1, y1, orientationB, symbol, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy+sxyFP, x1, y1, orientationB, symbol, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 30, doSave);
    % FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA+sxxFP, x1, y1, orientationB, symbol, 'xx', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyA+syyFP, x1, y1, orientationB, symbol, 'yy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyA+sxyFP, x1, y1, orientationB, symbol, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
    % FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN+sxxFP, x1, y1, orientationB, symbol, 'xx', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyN+syyFP, x1, y1, orientationB, symbol, 'yy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyN+sxyFP, x1, y1, orientationB, symbol, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 30, doSave, meanvalxy, stddevxy)
    
    orientationB = sprintf('Epar%d', node);
    linePlot(sxxApar(:, node)+sxxFP(:, node), sxxNpar(:, node)+sxxFP(:, node), txxPar(:, node)+sxxFP(:, node), orientationB, symbol, 'xx', 'Grid Point', '$\mu$', 30, doSave)
    linePlot(syyApar(:, node)+syyFP(:, node), syyNpar(:, node)+syyFP(:, node), tyyPar(:, node)+syyFP(:, node), orientationB, symbol, 'yy', 'Grid Point', '$\mu$', 30, doSave)
    linePlot(sxyApar(:, node)+sxyFP(:, node), sxyNpar(:, node)+sxyFP(:, node), txyPar(:, node)+sxyFP(:, node), orientationB, symbol, 'xy', 'Grid Point', '$\mu$', 30, doSave)
end
function fig = linePlot(analytic, numeric, head, orientationB, stress, component, xaxis, yaxis, fontSize, save)
    fig = figure();
    hold on
    plot(analytic, 'LineWidth', 2)
    plot(numeric, '--', 'LineWidth', 2)
    plot(head, ':', 'LineWidth', 2)
    hold off
    %     title(sprintf('$%s_{%s}^{\\textrm{%s}}$', stress, component, equation), 'Interpreter', 'latex', 'FontSize', fontSize)
    xlabel(sprintf('%s', xaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    ylabel(sprintf('$%s$', yaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    legend(sprintf('$%s_{%s}^{\\textrm{A}}$', stress, component), sprintf('$%s_{%s}^{\\textrm{N}}$', stress, component), sprintf('$%s_{%s}$', stress, component), 'Interpreter', 'latex', 'FontSize', fontSize)
    
    set(gca,'TickLabelInterpreter', 'latex');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel', a, 'fontSize', floor(fontSize*0.7), 'FontWeight', 'bold')
    set(gca,'XTickLabelMode','auto')
    xlim([1,size(analytic,1)])
    
    if save
        name = erase(sprintf('line_s%s%s', component, orientationB), ["\", "$"]);
        set(fig, 'Units', 'Inches');
        pos = get(fig, 'Position');
        set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3)+pos(3)*0.1, pos(4)])
        print(fig, sprintf('./paper/images/%s.pdf', name), '-dpdf', '-r0')
    end

end

function [fig, meanval, stddev] = plotCountourfSigmaHat(X, Y, Z, x0, y0, orientationB, stress, component, equation, xaxis, yaxis, units, fontSize, save, meanval, stddev, levels)
    fig = figure();

    if ~exist('levels', 'var')
        levels = 20;
    end

    contourf(X, Y, Z, levels);
    colormap(parula)
    cb = colorbar;
    set(cb,'FontSize', floor(fontSize*0.7), 'TickLabelInterpreter','latex')
  
    

    if ~exist('meanval', 'var') ||~exist('stddev', 'var')
        meanval = mean(Z, 'all');
        stddev = std(Z, 0, 'all');
    end

    displace = 5 * stddev;
    limits = [meanval - displace, meanval + displace];
    caxis(limits)
    title(sprintf('$%s_{%s}^{\\textrm{%s}}$', stress, component, equation), 'Interpreter', 'latex', 'FontSize', fontSize)
    xlabel(sprintf('$%s$', xaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    ylabel(sprintf('$%s$', yaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    ylabel(cb, sprintf('%s', units), 'Interpreter', 'latex', 'FontSize', fontSize)
    
    hold on
    plot(x0, y0, '.', 'color', 'white', 'MarkerSize', fontSize)
    hold off
    set(gca,'TickLabelInterpreter', 'latex');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel', a, 'FontWeight', 'bold', 'FontSize', floor(fontSize*0.7))
    set(gca,'XTickLabelMode','auto')

    if save
        name = erase(sprintf('s%s%s%s', component, equation, orientationB), ["\", "$"]);
        set(fig, 'Units', 'Inches');
        pos = get(fig, 'Position');
        set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3)+pos(3)*0.1, pos(4)])
        print(fig, sprintf('./paper/images/%s.pdf', name), '-dpdf', '-r0')
    end

end

function sigma = hatStressSurf(uhat, nc, x, D, mx, mz, w, h, d, X, Y, Z)

    gridSize = size(X);

    sigma = zeros(3, 3, gridSize(1), gridSize(2));
    x0 = zeros(3, 1);

    for col = 1:gridSize(2)

        for row = 1:gridSize(1)
            x0(1) = X(row, col);
            x0(2) = Y(row, col);
            x0(3) = Z(row, col);
            sigma(:, :, row, col) = hatStress(uhat, nc, x, D, mx, mz, w, h, d, x0);
        end

    end

end

function sigma = FieldPointStressSurf(X, Y, Z, x1, x2, b, a, mu, nu)
    gridSize = size(X);
    sigma = zeros(3, 3, gridSize(1), gridSize(2));
    x0 = zeros(1, 3);

    for col = 1:gridSize(2)

        for row = 1:gridSize(1)
            x0(1) = X(row, col);
            x0(2) = Y(row, col);
            x0(3) = Z(row, col);
            stress = FieldPointStress(x0, x1, x2, b, a, mu, nu);

            sigma(1, 1, row, col) = stress(:, 1);
            sigma(2, 2, row, col) = stress(:, 2);
            sigma(3, 3, row, col) = stress(:, 3);
            sigma(1, 2, row, col) = stress(:, 4);
            sigma(2, 1, row, col) = stress(:, 4);
            sigma(2, 3, row, col) = stress(:, 5);
            sigma(3, 2, row, col) = stress(:, 5);
            sigma(1, 3, row, col) = stress(:, 6);
            sigma(3, 1, row, col) = stress(:, 6);
        end

    end

end

function [txz, tyz] = imageStressAnalyticScrew(mu, b, nu, x, y, a, c)
    E = (2 * (1 + nu)) * mu;
    D = E .* b ./ (4 .* pi .* (1 - nu.^2));

    xpa = x + a;
    ymc = y - c;

    txz = -ymc ./ (xpa.^2 + ymc.^2);
    tyz = -xpa ./ (xpa.^2 + ymc.^2);
    txz = -D .* txz;
    tyz = D .* tyz;

end

function [txz, tyz] = FPStressAnalyticScrew(mu, b, nu, x, y, a, c)
    E = (2 * (1 + nu)) * mu;
    D = E .* b ./ (4 .* pi .* (1 - nu.^2));

    xma = x - a;
    xpa = x + a;
    ymc = y - c;

    txz = ymc ./ (xma.^2 + ymc.^2); %- ymc ./ (xpa.^2 + ymc.^2);
    tyz = xma ./ (xma.^2 + ymc.^2); %- xpa ./ (xpa.^2 + ymc.^2);
    txz = -D .* txz;
    tyz = D .* tyz;

end

function [txx, tyy, txy] = imageStressAnalyticEdgePerp(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (x, y) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (a, c).
    % b perpendicular to surface.
    %%%
    E = (2 * (1 + nu)) * mu;

    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx = ... - ymc .* (3 .* xma2 + ymc2) ./ den1 + ...
        ymc .* (3 .* xpa2 + ymc2) ./ den2 + ...
        4 .* a .* x .* ymc .* (3 .* xpa2 - ymc2) ./ den3;
    txx = D .* txx;

    tyy = ...%ymc .* (xma2 - ymc2)./den1 + ...
    -ymc .* (xpa2 - ymc2) ./ den2 + ...
        4 .* a .* ymc .* ((2 .* a - x) .* xpa2 + (3 .* x + 2 .* a) .* ymc2) ./ den3;
    tyy = D .* tyy;

    txy = ...%xma .* (xma2 - ymc2)./den1 + ...
    -xpa .* (xpa2 - ymc2) ./ den2 + ...
        2 .* a .* (-xma .* xpa .* xpa2 + 6 .* x .* xpa .* ymc2 - ymc2 .* ymc2) ./ den3;
    txy = D .* txy;
end

function [txx, tyy, txy] = totalStressAnalyticEdgePerp(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (x, y) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (a, c).
    % b perpendicular to surface.
    %%%
    E = (2 * (1 + nu)) * mu;

    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx =- ymc .* (3 .* xma2 + ymc2) ./ den1 + ...
        ymc .* (3 .* xpa2 + ymc2) ./ den2 + ...
        4 .* a .* x .* ymc .* (3 .* xpa2 - ymc2) ./ den3;
    txx = D .* txx;

    tyy = ymc .* (xma2 - ymc2) ./ den1 + ...
        -ymc .* (xpa2 - ymc2) ./ den2 + ...
        4 .* a .* ymc .* ((2 .* a - x) .* xpa2 + (3 .* x + 2 .* a) .* ymc2) ./ den3;
    tyy = D .* tyy;

    txy = xma .* (xma2 - ymc2) ./ den1 + ...
        -xpa .* (xpa2 - ymc2) ./ den2 + ...
        2 .* a .* (-xma .* xpa .* xpa2 + 6 .* x .* xpa .* ymc2 - ymc2 .* ymc2) ./ den3;
    txy = D .* txy;
end

function [txx, tyy, txy] = FPStressAnalyticEdgePerp(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (x, y) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (a, c).
    % b perpendicular to surface.
    %%%
    E = (2 * (1 + nu)) * mu;

    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx = -ymc .* (3 .* xma2 + ymc2) ./ den1;
    %                 + ...
    %                 ymc .* (3.*xpa2 + ymc2)./den2 + ...
    %       4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
    txx = D .* txx;

    tyy = ymc .* (xma2 - ymc2) ./ den1;
    %               + ...
    %              -ymc .* (xpa2 - ymc2)./den2 + ...
    %       4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
    tyy = D .* tyy;

    txy = xma .* (xma2 - ymc2) ./ den1;
    %                 + ...
    %                -xpa .* (xpa2 - ymc2)./den2 + ...
    %       2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
    txy = D .* txy;
end

function [txx, tyy, txy] = imageStressAnalyticEdgePar(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (a, c) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (x, y).
    % p parallel to surface.
    %%%
    E = (2 * (1 + nu)) * mu;
    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx = ...%xma .* (xma2 - ymc2) ./ den1 + ...
    -xpa .* (xpa2 - ymc2) ./ den2 + ...
        2 .* a .* ((3 * x + a) .* xpa2 .* xpa - 6 .* x .* xpa .* ymc2 - ymc2 .* ymc2) ./ den3;
    %%%
    txx = D .* txx;

    tyy = ...%xma .* (xma2 + 3.*ymc2) ./ den1 + ...
    -xpa .* (xpa2 + 3 .* ymc2) ./ den2 + ...
        -2 .* a .* (xma .* xpa .* xpa2 - 6 .* x .* xpa .* ymc2 + ymc2 .* ymc2) ./ den3;
    tyy = D .* tyy;

    txy = ...%ymc .* (xma2 - ymc2) ./ den1 + ...
    -ymc .* (xpa2 - ymc2) ./ den2 + ...
        4 .* a .* x .* ymc .* (3 .* xpa2 - ymc2) ./ den3;
    txy = D .* txy;
end

function [txx, tyy, txy] = totalStressAnalyticEdgePar(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (a, c) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (x, y).
    % p parallel to surface.
    %%%
    E = (2 * (1 + nu)) * mu;
    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx = xma .* (xma2 - ymc2) ./ den1 + ...
        -xpa .* (xpa2 - ymc2) ./ den2 + ...
        2 .* a .* ((3 * x + a) .* xpa2 .* xpa - 6 .* x .* xpa .* ymc2 - ymc2 .* ymc2) ./ den3;
    txx = D .* txx;

    tyy = xma .* (xma2 + 3 .* ymc2) ./ den1 + ...
        -xpa .* (xpa2 + 3 .* ymc2) ./ den2 + ...
        -2 .* a .* (xma .* xpa .* xpa2 - 6 .* x .* xpa .* ymc2 + ymc2 .* ymc2) ./ den3;
    tyy = D .* tyy;

    txy = ymc .* (xma2 - ymc2) ./ den1 + ...
        -ymc .* (xpa2 - ymc2) ./ den2 + ...
        4 .* a .* x .* ymc .* (3 .* xpa2 - ymc2) ./ den3;
    txy = D .* txy;
end

function [txx, tyy, txy] = FPStressAnalyticEdgePar(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (a, c) induced by edge dislocation parallel to the
    % surface at x = 0. Dislocation coordinates are (x, y).
    % p parallel to surface.
    %%%
    E = (2 * (1 + nu)) * mu;
    D = E .* b ./ (4 .* pi .* (1 - nu.^2));
    ymc = y - c;
    ymc2 = ymc.^2;
    xma = x - a;
    xma2 = xma.^2;
    xpa = x + a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx = xma .* (xma2 - ymc2) ./ den1;
    %             + ...
    %              -xpa .* (xpa2 - ymc2) ./ den2 + ...
    %       2.*a .* ymc .* ((3.*x+a) .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 - ymc2.*ymc2) ./ den3;
    txx = D .* txx;

    tyy = xma .* (xma2 + 3 .* ymc2) ./ den1;
    %                 + ...
    %                -xma .* (xpa2 + 3.*ymc2) ./ den2 + ...
    %       -2.*a .* (xma .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 + ymc2.*ymc2) ./ den3;
    tyy = D .* tyy;

    txy = ymc .* (xma2 - ymc2) ./ den1;
    %                 + ...
    %                 -ymc .* (xpa2 - ymc2) ./ den2 + ...
    %       4.*a.*x .* ymc .* (3.*xpa2 - ymc2) ./ den3;
    txy = D .* txy;
end
