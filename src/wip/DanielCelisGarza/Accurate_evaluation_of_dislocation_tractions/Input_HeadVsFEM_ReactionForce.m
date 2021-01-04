close all
clear all

CRYSTAL_STRUCTURE = 'bcc';
amag = 3.18e-4;

simTime = 0;
use_gpu = 0;
para_scheme = 0;
n_threads = 0;

MU = 1;
NU = 0.28;
a = 5;
bVec = [0 1 0];
nVec = [1 0 0];

planes = [1; 2; 3; 4; 5; 6];
dx = 2000;
dy = 2000;
dz = 2000;
mx = 40;
my = 40;
mz = 40;

mobility = @mobbcc1;
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
xcoord = xcoord(2) / 2; % middle of first element.
ycoord = dy / 2; % middle of the domain
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
rn(1, 4) = 7;
rn(end, 4) = 7;

for i = 1:len - 1
    links(i, :) = [i, i + 1, bVec, nVec];
end

hold on
plot3(rn(:, 1), rn(:, 2), rn(:, 3), 'r.')
plot3(X, Y, Z, 'k.')
hold off

u_dot = 0;
f_dot = 0;

loading = @staticSim;

CUDA_flag = false;
para_scheme = 2;

calculateTractions = @calculateNumericTractions;
plotFreq = 10;

addpath '../../../'