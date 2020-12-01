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
close all
clear all
%% SOURCE GENERATION PARAMETERS
amag = 3.18e-4;
mumag = 160E3; % MPa only used for plotting

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 1;
DIST_SOURCE = 0.5 / amag;

%% FEM PARAMETERS
%Cantilever

dx = 30 / amag; %10micron
dy = 5 / amag; %2micron
dz = 5 / amag; %2micron

mx = 20; % number of elements along beam length
% loading = @displacementControl;
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

c = 1 / amag;
x = 0.5 / (amag * sqrt(6));
y = 0.5 / (amag * sqrt(3));

rn = [6 * c 2.5 * c + 2 * x - y 0.5 * c 0;
    6 * c + x 2.5 * c - y 0.5 * c + x 0;
    6 * c + 2 * x 2.5 * c - 2 * x - y 0.5 * c + 2 * x 0;
    6 * c + 2 * x + y 2.5 * c - 2 * x 0.5 * c + 2 * x + y 0;
    6 * c + 2 * x + 2 * y 2.5 * c - 2 * x + y 0.5 * c + 2 * x + 2 * y 0;
    6 * c + x + 2 * y 2.5 * c + y 0.5 * c + x + 2 * y 0;
    6 * c + 2 * y 2.5 * c + 2 * x + y 0.5 * c + 2 * y 0;
    6 * c + y 2.5 * c + 2 * x 0.5 * c + y 0;
    8 * c 2.5 * c + 2 * x - y 0.5 * c 0;
    8 * c - x 2.5 * c - y 0.5 * c + x 0;
    8 * c - 2 * x 2.5 * c - 2 * x - y 0.5 * c + 2 * x 0;
    8 * c - 2 * x - y 2.5 * c - 2 * x 0.5 * c + 2 * x + y 0;
    8 * c - 2 * x - 2 * y 2.5 * c - 2 * x + y 0.5 * c + 2 * x + 2 * y 0;
    8 * c - x - 2 * y 2.5 * c + y 0.5 * c + x + 2 * y 0;
    8 * c - 2 * y 2.5 * c + 2 * x + y 0.5 * c + 2 * y 0;
    8 * c - y 2.5 * c + 2 * x 0.5 * c + y 0];

b1 = [-1 -1 -1] / 2;
b2 = [-1 1 1] / 2;
n1 = [-1 0 1] / sqrt(2);
n2 = [1 0 1] / sqrt(2);
% a = 5*norm(b1); % Default value of a has num tractions freak out, a = 5*||b|| is even worse.
a = 10 * norm(b1);
calculateTractions = @calculateAnalyticTractions;
% calculateTractions = @calculateNumericTractions;

links = [1 2 b1 n1;
    2 3 b1 n1;
    3 4 b1 n1;
    4 5 b1 n1;
    5 6 b1 n1;
    6 7 b1 n1;
    7 8 b1 n1;
    8 1 b1 n1;
    9 10 b2 n2;
    10 11 b2 n2;
    11 12 b2 n2;
    12 13 b2 n2;
    13 14 b2 n2;
    14 15 b2 n2;
    15 16 b2 n2;
    16 9 b2 n2];

plotFreq = 5;
saveFreq = 1e9;
u_dot = dx / 160E6;

calculateLoading = @sixStageDisplacementByEndLoad;
% calculateLoading = @constantLoading;
calculateLoadingFunctionArgs = struct('u_dot_0', u_dot, ...
    'u_bar_crit', [75; 115; 125; 145; 165], 'scaleFactor', [1/2; 1/5; 1/10; 1/25; 1/50]);

%create mirror of prismatic loops (outside boundary)
% rn_mirror = [rn(:,1)+dx , rn(:,2) , rn(:,3)+dx , zeros(length(rn(:,1)),1)+67];
% links_mirror = links;
% links_mirror(:,1) = links(:,2) + max(max(links(:,1:2)));
% links_mirror(:,2) = links(:,1) + max(max(links(:,1:2)));
%
% rn = vertcat(rn,rn_mirror);
% links = vertcat(links,links_mirror);

% %%
% %Edge and screw glide and climb mobility parameters
% mobility = 'mobbcc1';
% global Bscrew Bedge Beclimb Bline
% %Bedge=1e-4; %Pa s
% %Bscrew=1e-5; %Pa s
% %Beclimb=1e5; %Pa s - really big
% %Bline=1e-4*min(Bscrew,Bedge);
% Bedge = 1;
% Bscrew = 10;
% Beclimb = 1e10;
% Bline = 1e-4 * min(Bscrew, Bedge);

% %Meshing
% maxconnections = 4;
% lmax = 0.25 / amag;
% lmin = 0.1 / amag;
% areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
% areamax = 20 * areamin;
% doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
% docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
% doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
% dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

% %Simulation time
% dt0 = 1E5;

% intSimTime = 0;
% simTime = 0;
% %dtplot=2E-9; %2ns
% dtplot = 3E4;
% doplot = 1; % frame recording: 1 == on, 0 == off
% totalSimTime = 6E6;
% curstep = 0;
% simTime = 0;

% %Integrator
% integrator = 'int_trapezoid';
% %integrator='int_trapezoid_stoc'; %in development
% a = lmin / sqrt(3) * 0.5;
% Ec = MU / (4 * pi) * log(a / 0.1);
% rann = 0.5 * a;
% rntol = 0.5 * rann; % need to do convergence studies on all these parameters
% rmax = lmax;

% %Plotting
% plotFreq = 20;
% plim = 12 / amag; %12microns
% viewangle = [-35, 15];
% printfreq = 1;
% printnode = 2;

% rann = 0.9 * lmin;
% rmax = 0.5 * rann;
% integrator = 'int_trapezoid';
% a = 11.3473;
% totalSimTime = 1e9;
% dt0 = 1e60;
% printfreq = 100000;
% plotFreq = 1;
% Bedge = 1;
% Bscrew = 10;
% mobility = 'mobbcc_bb1b';
% a_trac = 1;

my = round(mx * dy / dx); % # elements in y direction
my = max(my, 1);
mz = round(mx * dz / dx); % elements in z direction
mz = max(mz, 1);
