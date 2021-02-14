close all
clear all

CRYSTAL_STRUCTURE = 'bcc';
amag = 3.18e-4;
simName = date;

simTime = 0;
use_gpu = 0;
para_scheme = 0;
n_threads = 0;

MU = 1;
NU = 0.28;
a = 5;

% No rotation matrix.
% bVec = [1 1 1];
% nVec = [1 -1 0];

% With rotation matrix.
n1 = [1 0 0]; e1 = [1; 1; 1]; e1 = e1 / norm(e1);
n2 = [0 1 0]; e2 = [1; -1; 0]; e2 = e2 / norm(e2);
n3 = cross(n1, n2); e3 = cross(e1, e2); e3 = e3 / norm(e3);
rotMatrix = [e1 e2 e3]';
bVec = [1 0 0] * sqrt(3) / 2;
nVec = [0 1 0];

bVec * rotMatrix
nVec * rotMatrix

planes = [1; 2; 3; 4; 5; 6];
dx = 2000;
dy = 2000;
dz = 2000;
mx = 20;
my = 20;
mz = 20;

% mobility = @mobbcc1;
mobility = @mobbcc_bb1b;
simName = strcat(simName, '_bb');
simType = @NumTracVsAnaTrac;

Fsim = [];
Usim = [];
t = [];

vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

plim = 12 / amag; %12microns

gridSize = mx;
x = linspace(0, dx, gridSize);
y = linspace(0, dy, gridSize);
z = linspace(0.5 * dz, 0.5 * dz, gridSize);
[X, Y] = meshgrid(x, y);
Z = meshgrid(z);

clear x y z;
len = 100;
xcoord = linspace(0, dx, gridSize);
ycoord = linspace(0, dy, gridSize);
%  xcoord = xcoord(2) / 4 % there were substantial differences
% xcoord = xcoord(2)/2; % middle of first element.
xcoord = dx / 32; % eigth of the domain.
% xcoord = 3*dx / 32; % eigth of the domain.
% xcoord = dy / 16; % middle of third element.
%% This makes the numeric tractions exit but not the analytic. Haven't found another case.
% ycoord = dy / 2 + 20; % middle of the domain
% xcoord = dx / 16; % eigth of the domain.
%%
ycoord = dy / 2; % middle of the domain
% ycoord = ycoord(2)/2; % middle of first element
% simName = strcat(simName, 'corner');
x = linspace(xcoord, xcoord, len);
y = linspace(ycoord, ycoord, len);
z = linspace(0, dz, len);
x1 = x(1);
y1 = y(1);
t = [0 0 1];
n = [1 0 0];
rn = zeros(len, 4);
rn(:, 1) = x;
rn(:, 2) = y;
rn(:, 3) = z;
links = zeros(len - 1, 8);
rn(1, 4) = 6; %61
rn(end, 4) = 6; %61
rn = [rn; x(1) y(1) -1000 * dz 67; x(1) y(1) 1000 * dz 67];

for i = 1:len - 1
    links(i, :) = [i, i + 1, bVec, nVec];
end

links = [links; len + 1, 1, bVec, nVec; len + 2, len, bVec, nVec];
% rntol = 45;

% %%
% xcoord = linspace(0, dx, gridSize);
% ycoord = linspace(0, dy, gridSize);
% %  xcoord = xcoord(2) / 4 % there were substantial differences
% % xcoord = xcoord(2)/2; % middle of first element.
% xcoord =  2*dx / 32; % eigth of the domain.
% % xcoord = dy / 16; % middle of third element.
% ycoord = 2*dy / 32; % middle of the domain
% % ycoord = ycoord(2)/2; % middle of first element
% simName = strcat(simName, 'corner');
% x = linspace(xcoord, xcoord, len);
% y = linspace(ycoord, ycoord, len);
% z = linspace(0, dz, len);
% x1 = x(1);
% y1 = y(1);
% t = [0 0 1];
% n = [1 0 0];
% rn2 = zeros(len, 4);
% rn2(:, 1) = x;
% rn2(:, 2) = y;
% rn2(:, 3) = z;
% links2 = zeros(len - 1, 8);
% rn2(1, 4) = 6;
% rn2(end, 4) = 6;
% rn2 = [rn2; x(1) y(1) -1000*dz 67; x(1) y(1) 1000*dz 67];
%
% for i = len+2+1:len+2+len - 1
%     links2(i - (len+2), :) = [i, i + 1, bVec, nVec];
% end
% links2 = [links2; len+2+len + 1, len+2+1, bVec, nVec; len+2+len + 2, len+2+len, bVec, nVec];
% links = [links; links2];
% rn = [rn; rn2];
%
% %%
% xcoord = linspace(0, dx, gridSize);
% ycoord = linspace(0, dy, gridSize);
% %  xcoord = xcoord(2) / 4 % there were substantial differences
% % xcoord = xcoord(2)/2; % middle of first element.
% xcoord =  1*dx / 32; % eigth of the domain.
% % xcoord = dy / 16; % middle of third element.
% ycoord = 1*dy / 32; % middle of the domain
% % ycoord = ycoord(2)/2; % middle of first element
% simName = strcat(simName, 'corner');
% x = linspace(xcoord, xcoord, len);
% y = linspace(ycoord, ycoord, len);
% z = linspace(0, dz, len);
% x1 = x(1);
% y1 = y(1);
% t = [0 0 1];
% n = [1 0 0];
% rn3 = zeros(len, 4);
% rn3(:, 1) = x;
% rn3(:, 2) = y;
% rn3(:, 3) = z;
% links3 = zeros(len - 1, 8);
% rn3(1, 4) = 6;
% rn3(end, 4) = 6;
% rn3 = [rn3; x(1) y(1) -1000*dz 67; x(1) y(1) 1000*dz 67];
%
% for i = len+2+len+2+1:len+2+len+2+len - 1
%     links3(i - (len+2+len+2), :) = [i, i + 1, bVec, nVec];
% end
% links3 = [links3; len+2+len+2+len + 1, len+2+len+2+1, bVec, nVec; len+2+len+2+len + 2, len+2+len+2+len, bVec, nVec];
% links = [links; links3];
% rn = [rn; rn3];

hold on
plot3(rn(:, 1), rn(:, 2), rn(:, 3), 'r.')
plot3(X, Y, Z, 'k.')
hold off

u_dot = 0;
f_dot = 0;

loading = @staticSim;

CUDA_flag = false;
para_scheme = 1;

% calculateTractions = @calculateNumericTractions;
% simName = strcat('numeric_', simName);
calculateTractions = @calculateAnalyticTractions;
simName = strcat('analytic_', simName);
plotFreq = 10;
saveFreq = 2; %4*plotFreq;

lmin = 10 * a;
lmax = 2.5 * lmin;

addpath '../../../'
Bcoeff = struct('screw', 10, 'edge', 1, 'climb', 1e10, 'line', 1e-4);
% dt0 = 1e9;
% totalSimTime = 1e12;
