%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that includes NUM_SOURCES number of sources in the form of pinned
% b=<111> sources of *constant* size, determined by parameter DIST_SOURCE
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
% Distances: m
% Time: s
% Force: N
% Pressure: Pa
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOURCE GENERATION PARAMETERS

amag = 3.14E-10;

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 1;
DIST_SOURCE = 35 * 10E-9 / amag;

%% FEM PARAMETERS
%Cantilever

dx = 12E-6 / amag; %12micron
dy = 2E-6 / amag; %2micron
dz = 2E-6 / amag; %2micron

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
%[rn,links] = DislocationSourceGenerator(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz,5);
[rn, links] = checkGenerator101(NUM_SOURCES, DIST_SOURCE, CRYSTAL_STRUCTURE, dx, dy, dz);

%Edge and screw glide and climb mobility parameters
mobility = 'mobbcc0';
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge = 1;
Bscrew = 0.1;
Beclimb = 1e15;
Bline = 1e-4 * min(Bscrew, Bedge);

%Meshing
maxconnections = 4;
lmax = 20 * 10E-9 / amag; %20nm
lmin = 1 * 10E-9 / amag; %1nm
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;
doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
dt0 = 10E-9;
%dt0=10E4;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot = 1E5;
doplot = 1; % frame recording: 1 == on, 0 == off
totalSimTime = 10E12;
curstep = 0;
simTime = 0;

%Integrator
integrator = 'int_trapezoid';
%integrator='int_trapezoid_stoc'; %in development
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
rann = 0.5 * a;
rntol = 0.5 * rann;
rmax = lmin;

%Plotting
plotFreq = 10;
plim = 12E-6 / amag; %12microns
viewangle = [-35, 15];
printfreq = 1;
printnode = 4;
