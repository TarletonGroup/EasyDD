%%% Script to compare Head's analytical solutions to those obtained via FEM using analytic and numeric tractions.
close all
clear all

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
% TODO #36
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
    x = linspace(0.1 * dx, 0.1 * dx, len);
    y = linspace(0.5 * dy, 0.5 * dy, len);
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
% TODO #38
% TODO #39
% TODO #40
addpath 'D:\DPhil\OneDrive - Nexus365\EasyDD\src'

for k = 1:2
    %     close all
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
    else
        [txx, tyy, txy] = imageStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
        [txxFP, tyyFP, txyFP] = FPStressAnalyticEdgePar(MU, b, NU, X, Y, x1, y1);
        orientationB = 'Epar';
    end

    % Head
    [~, meanvalxx, stddevxx] = plotCountourfSigmaHat(X, Y, txx, x1, y1, orientationB, 'xx', '', 'x,~b', 'y,~b', '$\mu$', 15, doSave);
    [~, meanvalyy, stddevyy] = plotCountourfSigmaHat(X, Y, tyy, x1, y1, orientationB, 'yy', '', 'x,~b', 'y,~b', '$\mu$', 15, doSave);
    [~, meanvalxy, stddevxy] = plotCountourfSigmaHat(X, Y, txy, x1, y1, orientationB, 'xy', '', 'x,~b', 'y,~b', '$\mu$', 15, doSave);
    % FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA, x1, y1, orientationB, 'xx', 'A', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyA, x1, y1, orientationB, 'yy', 'A', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyA, x1, y1, orientationB, 'xy', 'A', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalxy, stddevxy)
    % FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN, x1, y1, orientationB, 'xx', 'N', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalxx, stddevxx)
    plotCountourfSigmaHat(X, Y, syyN, x1, y1, orientationB, 'yy', 'N', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalyy, stddevyy)
    plotCountourfSigmaHat(X, Y, sxyN, x1, y1, orientationB, 'xy', 'N', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanvalxy, stddevxy)
    % Rel Err FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA ./ txx - 1, x1, y1, orientationB, 'xx', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyA ./ tyy - 1, x1, y1, orientationB, 'yy', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyA ./ txy - 1, x1, y1, orientationB, 'xy', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    % Rel Err FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN ./ txx - 1, x1, y1, orientationB, 'xx', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyN ./ tyy - 1, x1, y1, orientationB, 'yy', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyN ./ txy - 1, x1, y1, orientationB, 'xy', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
    % Abs Err FEM + Analytic tractions
    plotCountourfSigmaHat(X, Y, sxxA - txx, x1, y1, orientationB, 'xx', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, syyA - tyy, x1, y1, orientationB, 'yy', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyA - txy, x1, y1, orientationB, 'xy', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
    % Abs Err FEM + Numeric tractions
    plotCountourfSigmaHat(X, Y, sxxN - txx, x1, y1, orientationB, 'xx', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, syyN - tyy, x1, y1, orientationB, 'yy', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
    plotCountourfSigmaHat(X, Y, sxyN - txy, x1, y1, orientationB, 'xy', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
end

% %% Screw
b = [0 0 1];
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
sigmaN = hatStressSurf(uhat2, nc, xnodes, D, mx, mz, w, h, d, X, Y, Z);
syzA = squeeze(sigmaA(2, 3, :, :));
syzN = squeeze(sigmaN(2, 3, :, :));
[txz, tyz] = imageStressAnalyticScrew(MU, 1, NU, X, Y, x1, y1);
[~, meanval, stddev] = plotCountourfSigmaHat(X, Y, tyz, x1, y1, 'screw', 'yz', '', 'x,~b', 'y,~b', '$\mu$', 15, doSave);
plotCountourfSigmaHat(X, Y, syzA, x1, y1, 'screw', 'yz', 'A', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanval, stddev)
plotCountourfSigmaHat(X, Y, syzN, x1, y1, 'screw', 'yz', 'N', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanval, stddev)
plotCountourfSigmaHat(X, Y, syzA ./ tyz - 1, x1, y1, 'screw', 'yz', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
plotCountourfSigmaHat(X, Y, syzN ./ tyz - 1, x1, y1, 'screw', 'yz', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
plotCountourfSigmaHat(X, Y, syzA - tyz, x1, y1, 'screw', 'yz', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
plotCountourfSigmaHat(X, Y, syzN - tyz, x1, y1, 'screw', 'yz', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)

sxzA = squeeze(sigmaA(1, 3, :, :));
sxzN = squeeze(sigmaN(1, 3, :, :));
[~, meanval, stddev] = plotCountourfSigmaHat(X, Y, txz, x1, y1, 'screw', 'xz', '', 'x,~b', 'y,~b', '$\mu$', 15, doSave);
plotCountourfSigmaHat(X, Y, sxzA, x1, y1, 'screw', 'xz', 'A', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanval, stddev)
plotCountourfSigmaHat(X, Y, sxzN, x1, y1, 'screw', 'xz', 'N', 'x,~b', 'y,~b', '$\mu$', 15, doSave, meanval, stddev)
plotCountourfSigmaHat(X, Y, sxzA ./ txz - 1, x1, y1, 'screw', 'xz', '$\eta$A', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
plotCountourfSigmaHat(X, Y, sxzN ./ txz - 1, x1, y1, 'screw', 'xz', '$\eta$N', 'x,~b', 'y,~b', 'Rel Err', 15, doSave)
plotCountourfSigmaHat(X, Y, sxzA - txz, x1, y1, 'screw', 'xz', '$\Delta$A', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)
plotCountourfSigmaHat(X, Y, sxzN - txz, x1, y1, 'screw', 'xz', '$\Delta$N', 'x,~b', 'y,~b', 'Abs Err', 15, doSave)

% figCounter = figCounter + 1;
% figure(figCounter)
% txx = txx; %./norm(txx);
% meantxx = mean(txx, 'all');
% stdtxx = std(txx, 0, 'all');
% displace = 5 * stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X, Z, txx);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' sxx'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% tyy = tyy; %./norm(tyy);
% meantzz = mean(tyy, 'all');
% stdtzz = std(tyy, 0, 'all');
% displace = 5 * stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X, Z, tyy);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' szz'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% txy = txy; %./norm(txy);
% meantxz = mean(txy, 'all');
% stdtxz = std(txy, 0, 'all');
% displace = 5 * stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X, Z, txy);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' sxz'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% contourf(X, Z, sxxFP);
% colormap(parula)
% colorbar
% title('FP FE sxx')
% xlabel('x')
% ylabel('z')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% contourf(X, Z, szzFP);
% colormap(parula)
% colorbar
% title('FP FE szz')
% xlabel('x')
% ylabel('z')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% contourf(X, Z, sxzFP);
% colormap(parula)
% colorbar
% title('FP FE sxz')
% xlabel('x')
% ylabel('z')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% txxFP = txxFP; %./norm(txx);
% meantxx = mean(txxFP, 'all');
% stdtxx = std(txxFP, 0, 'all');
% displace = 5 * stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X, Z, txxFP);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' FP sxx'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% tzzFP = tzzFP; %./norm(tyy);
% meantzz = mean(tzzFP, 'all');
% stdtzz = std(tzzFP, 0, 'all');
% displace = 5 * stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X, Z, tzzFP);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' FP szz'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off
%
% figCounter = figCounter + 1;
% figure(figCounter)
% txzFP = txzFP; %./norm(txy);
% meantxz = mean(txzFP, 'all');
% stdtxz = std(txzFP, 0, 'all');
% displace = 5 * stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X, Z, txzFP);
% colormap(parula)
% colorbar
% caxis(limits)
% title(strcat('b', name, ' FP sxz'))
% xlabel('b')
% ylabel('b')
% hold on
% plot(x1, y1, '.', 'color', 'black', 'MarkerSize', 10)
% hold off

function [fig, meanval, stddev] = plotCountourfSigmaHat(X, Y, Z, x0, y0, orientationB, component, equation, xaxis, yaxis, units, fontSize, save, meanval, stddev, levels)
    fig = figure();
    if ~exist('levels', 'var')
        levels = 20;
    end
    contourf(X, Y, Z, levels);
    colormap(parula)
    cb = colorbar;
    if ~exist('meanval', 'var') || ~exist('stddev', 'var')
        meanval = mean(Z, 'all');
        stddev = std(Z, 0, 'all');
    end
    displace = 5 * stddev;
    limits = [meanval - displace, meanval + displace];
    caxis(limits)
    title(sprintf('$\\hat{\\sigma}_{%s}^{\\textrm{%s}}$', component, equation), 'Interpreter', 'latex', 'FontSize', fontSize)
    xlabel(sprintf('$%s$', xaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    ylabel(sprintf('$%s$', yaxis), 'Interpreter', 'latex', 'FontSize', fontSize)
    ylabel(cb, sprintf('%s', units), 'Interpreter', 'latex', 'FontSize', fontSize)
    hold on
    plot(x0, y0, '.', 'color', 'black', 'MarkerSize', fontSize)
    hold off

    if save
        name = erase(sprintf('s%s%s%s', component, equation, orientationB), ["\", "$"]);
        set(fig, 'Units', 'Inches');
        pos = get(fig, 'Position');
        set(fig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
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

    xma = x - a;
    xpa = x + a;
    ymc = y - c;

    txz = -ymc ./ (xpa.^2 + ymc.^2);
    tyz = -xpa ./ (xpa.^2 + ymc.^2);
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
