%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that includes checking some simple cases
% in a microcantilever of bcc tungsten.
%
% Sections are the following:
% 1) SOURCE GENERATION PARAMETERS
% 2) FEM PARAMETERS
% 3) MATERIAL CONSTANTS
% 4) DDLab PARAMETERS
%       - Dislocation nodes (rn) and segments (links) generator
%       - Edge and Screw Glide and Climb Mobility Parameters
%       - Meshing
%       - Simulation time
%       - Integrator
%       - Plotting
%
% Note SI units
% Distances: microns
% Time: s
% Force: N
% Pressure: Pa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% SOURCE GENERATION PARAMETERS
bmag = 1e-4;
mumag = 1; %160E3;% MPa only used for plotting

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 1;
DIST_SOURCE = 2.5 / bmag;

%% FEM PARAMETERS
%Cantilever

dx = 30 / bmag; %10micron
dy = 5 / bmag; %2micron
dz = 5 / bmag; %2micron

mx = 50; % number of elements along beam length
loading = 1;
vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

%% MATERIAL CONSTANTS

%MU = 160E9; %160GPa
MU = 1;
NU = 0.305;

%% DDLab PARAMETERS

%Dislocation nodes and segments generator
[rn, links] = checkGenerator(NUM_SOURCES, DIST_SOURCE, CRYSTAL_STRUCTURE, dx, dy, dz);

midPTS(1, :) = [0.1 * dx, 0.2 * dy, 0.5 * dz];
normal = [0.01 0.99 0];
% generate random <111> b vector
% b_vec = pmone(normal,NUM_SOURCES);
b_vec = [0 0.99 0]; % use [111] for debugging

%We have thus defined slip-plane and b-vector of the loops.
%We try a loop placed within the slip-plane (shear loop)...

% seg_vec = cross(normal',b_vec');
edge = 0.5 * [0.01 0 -0.99] * DIST_SOURCE;
screw = 0.5 * [0.99 0.01 0] * DIST_SOURCE;

rn = zeros(3, 4);
links = zeros(2, 8);

% pure edge segment
%     r1 = midPTS(p,:)
%     r2 = midPTS(p,:) + 0.5*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%     r3 = midPTS(p,:) + 1*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%
% pure screw segment
%     r1 = midPTS(p,:) - 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:));
%     r2 = midPTS(p,:);
%     r3 = midPTS(p,:) + 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:));

r1 = midPTS(1, :) - edge;
r2 = midPTS(1, :);
r3 = midPTS(1, :) + edge;
r4 = r3 + screw;
r5 = r4 + screw;
r6 = r5 - edge;
r7 = r6 - edge;
r8 = r7 - screw;

rn(1, :) = [r1 7];
rn(2, :) = [r2 0];
rn(3, :) = [r3 7];
rn(4, :) = [r4 0];
rn(5, :) = [r5 7];
rn(6, :) = [r6 0];
rn(7, :) = [r7 7];
rn(8, :) = [r8 0];

for p = 1:7
    links(p, 1:2) = [p, p + 1];
end

links(8, 1:2) = [8, 1];

links(1:8, 3:5) = repmat(b_vec, 8, 1);
links(1:8, 6:8) = repmat(normal, 8, 1);

%%
%Edge and screw glide and climb mobility parameters
mobility = 'mobbcc1';
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge = 1;
Bscrew = 10;
Beclimb = 1e10;
Bline = 1e-4 * min(Bscrew, Bedge);

%Meshing
maxconnections = 4;
lmax = 0.5 / bmag; %1 microns
lmin = lmax / 10; % 0.05 microns
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;
doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
%dt0=10E-9; %1 ns
dt0 = 1E9;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot = 1E9;
doplot = 1; % frame recording: 1 == on, 0 == off
totalSimTime = 1E12;
curstep = 0;
simTime = 0;

%Integrator
integrator = 'int_trapezoid';
%integrator='int_trapezoid_stoc'; %in development
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
rann = 0.5 * a;
rntol = 0.5 * rann; % need to do convergence studies on all these parameters
rmax = lmax;

%Plotting
plotFreq = 10;
plim = 12 / bmag; %12microns
viewangle = [-35, 15];
printfreq = 1;
printnode = 2;
