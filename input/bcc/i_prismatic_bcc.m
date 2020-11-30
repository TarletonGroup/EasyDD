%=========================================================================%
%-------------------------------------------------------------------------%
% Daniel Celis Garza
% 04/07/2018
%-------------------------------------------------------------------------%
%
% Input for prismatic bcc loops.
%
% Sections are the following:
% 1) Source generation
% 2) Material constants
% 3) FEM parameters
% 4) DDLab variables
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
%
%=========================================================================%
% Variables
%=========================================================================%
%
% Global variables are defined in the DDLab Manual.
%
%-------------------------------------------------------------------------%
% Source generation.
%-------------------------------------------------------------------------%
%
% a_m  := magnitude of dislocation core.
% mu_m := magnitude of mu in MPa. Only used in plotting.
% b_m  := magnitude of Burgers vector.
% crys_struct := crystal structure.
% n_source := number of dislocation sources.
% d_source := distance between sources.
% l_source := length of source.
%
%-------------------------------------------------------------------------%
% FEM parameters.
%-------------------------------------------------------------------------%
%
% dx, dy, dz := dimensions in the x, y, z dimensions.
% mx         := number of elements along the x dimension.
% vertices   := vertices of the cantilever.
% loading    := loading type.
%
%-------------------------------------------------------------------------%
% DDLab variables.
%-------------------------------------------------------------------------%
%
% mobility   := dislocation mobility.
%
% Plotting.
%
% plotFreq  := plot frequency.
% plim      := plot limit.
% printfreq := print frequency.
% printnode :=
% viewangle := view angle for plot in deg.
%
% Meshing.
%
% maxconnections := maximum connectivity.
% lmax  := maximum segment length.
% lmin  := minimum segment length.
% areamin  := minimum area enclosed by two connected segments.
% areamax  := minimum area enclosed by two connected segments.
%
% Runtime flags. 1 turns flag on, 0 turns flag off.
%
% doplot       := plot while the simulation runs.
% dosave       := save .mat file at every multiple of savefreq iterations.
% doremesh     := dislocation network remesh.
% docollision  := collision detection between dislocations.
% doseparation := splitting dislocation network.
% dovirtmesh   := remeshing of virtual nodes.
%
% Runtime parameters.
%
% savefreq := save frequency of the .mat file.
%=========================================================================%

clear all;
close all;

%% Source generation.
% Values are for tungsten.
a_m = 2.856e-4;
mu_m = 82E3;
b_m = sqrt(3) / 2;

crys_struct = 'bcc';

n_source = 10;
d_source = 0.2 / a_m;
l_source = 80 * b_m;

%% Material constants.
global MU NU

% Values are for tungsten.
MU = 1;
NU = 0.29;

%% FEM parameters.
% Rectangular prism cantilever.

% Dimensions.
dx = 12 / a_m;
dy = 3 / a_m;
dz = 3 / a_m;

vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

fem_n = [-1 0 0; ...% min(x), yz-plane, face 1
    1 0 0; ...% max(x), yz-plane, face 2
    0 -1 0; ...% min(y), xz-plane, face 4
    0 1 0; ...% max(y), xz-plane, face 3
    0 0 -1; ...% min(z), xy-plane, face 5
    0 0 1]; % max(z), xy-plane, face 6
% Loading type.
loading = 1;

%% DDLab variables.
global Bscrew Bedge Beclimb Bline
global cOVERa
global maxSeg
global burgsref planesref edgesref bpiref ppbref planestyperef burgnumbers planenumbers
global node_NRF flag_NRF flag_Junction node_Junction1 node_Junction2
global unloadflag unloadcount
global Ubarglobal

% Flags.
flag_NRF = 0;
flag_Junction = 0;
unloadflag = 0;
unloadcount = 0;
Ubarglobal = 0;

% Mobility law.
mobility = 'Hmobbcc1';

% Values are for tungsten.
Bedge = 5E-4;
Bscrew = 100E-4;
Beclimb = 1e5;
Bline = 1e-4 * min(Bscrew, Bedge);

% Plotting parameters.
plotFreq = 1000;
plim = 12 / a_m;
printfreq = 1000;
printnode = 2;
viewangle = [30, 60];

% Dislocation node and segment generator.
slips = [3];
[rn, links] = prismatic_bcc_generator(slips, dx, dy, dz, vertices, fem_n);

% Meshing parameters.
maxconnections = 4;
lmax = 0.1 / a_m;
lmin = 0.02 / a_m;
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;

% DDLab flags.
doplot = 1;
dosave = 0;
doremesh = 1;
docollision = 1;
doseparation = 1;
dovirtmesh = 1;

% Runtime parameters.
savefreq = 10;

% Temporal variables.
dt0 = 1E8;
dt = 1E8;
dtplot = 0;

intSimTime = 0;
sinTime = 0;

totalSimTime = 1E15;
curstep = 0;
simTime = 0;

% Integrator. All these require convergence studies.
integrator = 'int_trapezoid';
a = 40 * b_m;
Ec = MU / (4 * pi) * log(a / 0.1);
rann = 0.5 * a;
rntol = 0.5 * rann;
rmax = lmax;
maxSeg = lmax;

%plotnodes(rn,links,0,vertices)
