clear all

global Bscrew Bedge Beclimb Bline

b = [0 0 1];
b1 = [0 0 1];

P1 = [(7.2993e+03 - 1000) (1.8248e+03 - 1000) 1.8248e+03 + 1000];
P2 = [(7.2993e+03 + 1000) (1.8248e+03 - 1000) 1.8248e+03 + 1000];
P3 = [(7.2993e+03 + 1000) (1.8248e+03 + 1000) 1.8248e+03 + 1000];
P4 = [(7.2993e+03 - 1000) (1.8248e+03 + 1000) 1.8248e+03 + 1000];

P5 = [(7.2993e+03 - 1000) (1.8248e+03 - 1000) 1.8248e+03 - 1000];
P6 = [(7.2993e+03 + 1000) (1.8248e+03 - 1000) 1.8248e+03 - 1000];
P7 = [(7.2993e+03 + 1000) (1.8248e+03 + 1000) 1.8248e+03 - 1000];
P8 = [(7.2993e+03 - 1000) (1.8248e+03 + 1000) 1.8248e+03 - 1000];

rn = [
    P1, 0;
    P2, 0;
    P3, 0;
    P4, 0;
    P5, 0;
    P6, 0;
    P7, 0;
    P8, 0;
    ];

n1 = cross(P2 - P1, b);
n2 = cross(P3 - P2, b);
n3 = cross(P4 - P3, b);
n4 = cross(P4 - P1, b);

links = [1 4 b n1;
    4 3 b n2;
    3 2 b n3;
    2 1 b n4;
    8 5 b1 n1;
    5 6 b1 n4;
    6 7 b1 n3;
    7 8 b1 n2;

    ];

MU = 1;
NU = 0.305;
maxconnections = 8;
lmax = 1000;
lmin = 50;
coplanar_tol = lmin * 0.5;
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
dt0 = 1e7;
mobility = 'mobbcc0';

intSimTime = 0;
sinTime = 0;
dtplot = 3e6; % 3 ms (code in units of ns)
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
appliedstress = 0 * [2 0 1; 0 2 -1; 1 -1 0];
viewangle = [-35, 15];
printfreq = 1;
printnode = 3;
rmax = 100; %maximum distance a node may travel in one cycle

%FEM CANTILEVER PARAMETERS

dx = 6; %microns
dy = 1; %microns
dz = 1; %microns
%tungsten "a" is 0.000274
dx = dx / 0.000274;
dy = dy / 0.000274;
dz = dz / 0.000274;

mx = 12; % number of elements along beam length
loading = 1;
vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];
