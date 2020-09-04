clear all

global Bscrew Bedge Beclimb Bline

bmag = 2.74E-4;

l = 0.5 / bmag;
b = [1 0 1];

side1 = [-1 0 1];
side2 = [0 -1 0];
side3 = -side1;
side4 = -side2;

n = cross(side2, b);

%start from center of grid

%P1 = [0.52 0.52  1.7]/bmag; %exiting top
%P1 = [0.52 0.1 0.52]/bmag; %lateral

%Loop n.1, inside [101]
P1 = [1.52 1.2 1.2] / bmag;
P2m = P1 + l * side1 * 0.5;
P2 = P1 + l * side1;
P3m = P2 + l * side2 * 0.5;
P3 = P2 + l * side2;
P4 = P3 + l * side3;
P4m = P3 + l * side3 * 0.5;
P1m = P4 + l * side4 * 0.5;

%Loop n.2, inside [101]
P5 = [2.52 1.2 1.2] / bmag;
P6m = P5 + l * side1 * 0.5;
P6 = P5 + l * side1;
P7m = P6 + l * side2 * 0.5;
P7 = P6 + l * side2;
P8 = P7 + l * side3;
P8m = P7 + l * side3 * 0.5;
P5m = P8 + l * side4 * 0.5;

%Loop n.3, inside [101]
P9 = [4.52 1.2 1.2] / bmag;
P10m = P9 + l * side1 * 0.5;
P10 = P9 + l * side1;
P11m = P10 + l * side2 * 0.5;
P11 = P10 + l * side2;
P12 = P11 + l * side3;
P12m = P11 + l * side3 * 0.5;
P9m = P9 + l * side4 * 0.5;

rn = [
    P1, 7;
    P2m, 0;
    P2, 7;
    P3m, 0;
    P3, 7;
    P4m, 0;
    P4, 7;
    P1m, 0;

    P5, 7;
    P6m, 0;
    P6, 7;
    P7m, 0;
    P7, 7;
    P8m, 0;
    P8, 7;
    P5m, 0;

    P9, 7;
    P10m, 0;
    P10, 7;
    P11m, 0;
    P11, 7;
    P12m, 0;
    P12, 7;
    P9m, 0;
    ];

links = [2 1 b n;
    3 2 b n;
    4 3 b n;
    5 4 b n;
    6 5 b n;
    7 6 b n;
    8 7 b n;
    1 8 b n;

    2 + 8 1 + 8 b n;
    3 + 8 2 + 8 b n;
    4 + 8 3 + 8 b n;
    5 + 8 4 + 8 b n;
    6 + 8 5 + 8 b n;
    7 + 8 6 + 8 b n;
    8 + 8 7 + 8 b n;
    1 + 8 8 + 8 b n;

    2 + 16 1 + 16 b n;
    3 + 16 2 + 16 b n;
    4 + 16 3 + 16 b n;
    5 + 16 4 + 16 b n;
    6 + 16 5 + 16 b n;
    7 + 16 6 + 16 b n;
    8 + 16 7 + 16 b n;
    1 + 16 8 + 16 b n;
    ];

MU = 1;
NU = 0.305;
maxconnections = 8;
lmax = 500;
lmin = 100;
coplanar_tol = 10;
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
dt0 = 1e5;
mobility = 'mobbcc0';

intSimTime = 0;
sinTime = 0;
dtplot = 5e5;
doplot = 1; % frame recording: 1 == on, 0 == off
totalSimTime = 0.3e12;
curstep = 0;
simTime = 0;

Bscrew = 1e0;
Bedge = 1e0;
Beclimb = 1e10;
Bline = 1.0e-4 * min(Bscrew, Bedge);

integrator = 'int_trapezoid';
rann = 0.5 * a;
rntol = 0.5 * rann;
doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on
plotFreq = 1;
plim = 7000;
viewangle = [-35, 15];
printfreq = 1;
printnode = 4;
rmax = 100; %maximum distance a node may travel in one cycle

%FEM CANTILEVER PARAMETERS

dx = 12; %microns
dy = 2; %microns
dz = 2; %microns
%tungsten "a" is 0.000274
dx = dx / bmag;
dy = dy / bmag;
dz = dz / bmag;

mx = 30; % number of elements along beam length
loading = 1;
vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

%plotnodes(rn,links,plim,vertices);
hold on;
plot3(rn(:, 1), rn(:, 2), rn(:, 3), 'ko'); %nodes
