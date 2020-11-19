%%
% Daniel Celis Garza
% 19/11/2018
% Cubic domain. Single fixed end. No loading.
clear all
close all
amag = 1;
CRYSTAL_STRUCTURE = 'bcc';
simTime = 0;
totalSimTime = 1e14;
dt0 = 1e9;
use_gpu = 0;
para_scheme = 0;
n_threads = 0;

MU = 1;
NU = 0.28;

a = 10;
lmin = 2 * a;
lmax = 10 * a;
lmean = (lmin + lmax) / 2;

planes = [1; 2; 3; 4; 5; 6];
dx = 200;
dy = 200;
dz = 200;
mx = 10;
my = 10;
mz = 10;

b = [1 1 -1] / sqrt(3);
n = [1 0 1] / sqrt(2);

centre = [dx / 2 dy / 2 dz / 2];
edgeVec = cross(n, b) * lmean;
screwVec = b * lmean;

r1 = centre - edgeVec - screwVec;
r2 = r1 + edgeVec;
r3 = r2 + edgeVec;
r4 = r3 + screwVec;
r5 = r4 + screwVec;
r6 = r5 - edgeVec;
r7 = r6 - edgeVec;
r8 = r7 - screwVec;

rn = zeros(8, 4);

rn(1, 1:3) = r1;
rn(2, 1:3) = r2;
rn(3, 1:3) = r3;
rn(4, 1:3) = r4;
rn(5, 1:3) = r5;
rn(6, 1:3) = r6;
rn(7, 1:3) = r7;
rn(8, 1:3) = r8;

rn(:,1) = rn(:,1) + dx/5;
rn(:,2) = rn(:,2) + dy/15;
rn(:,3) = rn(:,3) + dz/5;

links = [
    1 2 b n;
    2 3 b n;
    3 4 b n;
    4 5 b n;
    5 6 b n;
    6 7 b n;
    7 8 b n;
    8 1 b n;
    ];

simType = 'atRest';
loading = 'noLoadNoDisp';
saveFreq = 5000;
plotFreq = 1000;
a_trac = 1;
simName = sprintf('%d_%s', a_trac, date);
