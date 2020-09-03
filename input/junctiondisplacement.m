clear all

global Bscrew Bedge Beclimb Bline

l = 250;
b = 0.5 * [1 1 1];

side1 = [1 -2 1];
side2 = [2 -1 -1];
side3 = [1 1 -2];
side4 = [-1 2 -1];
side5 = [-2 1 1];
side6 = [-1 -1 2];

n1 = cross(side1, b);
n2 = cross(side2, b);
n3 = cross(side3, b);
n4 = cross(side4, b);
n5 = cross(side5, b);
n6 = cross(side6, b);

%start from center of grid

% Loop 1
A1 = [2500 2500 900];
A2 = A1 + l * side1;
A3 = A2 + l * side2;
A4 = A3 + l * side3;
A5 = A4 + l * side4;
A6 = A5 + l * side5;
A7 = A6 + l * side6; %should equal P1, consistency check

% Loop 2
B1 = A1 - (side1 + side2) * l;
B2 = B1 + l * side1;
B3 = B2 + l * side2;
B4 = B3 + l * side3;
B5 = B4 + l * side4;
B6 = B5 + l * side5;
B7 = B6 + l * side6; %should equal P1, consistency check

rn = [
    A1, 0;
    A2, 0;
    A3, 0;
    A4, 0;
    A5, 0;
    A6, 0;

    B1, 0;
    B2, 0;
    B3, 0;
    B4, 0;
    B5, 0;
    B6, 0;
    ];

links = [1 2 b n1;
    2 3 b n2;
    3 4 b n3;
    4 5 b n4;
    5 6 b n5;
    6 1 b n6;

    7 8 b n1;
    8 9 b n2;
    9 10 b n3;
    10 11 b n4;
    11 12 b n5;
    12 7 b n6;
    ];

MU = 1;
NU = 0.305;
maxconnections = 8;
lmax = 1000;
lmin = 100;
areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
areamax = 20 * areamin;
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
dt0 = 1;
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

global USING_GPU;
USING_GPU = 0; %0 if CPU only.

integrator = 'int_trapezoid';
rann = 0.5 * a;
rntol = 0.5 * rann;
doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on
plotFreq = 1;
plim = 7000;
appliedstress = 10^ - 2 * [2 0 1; 0 2 -1; 1 -1 0];
viewangle = [90 0];
printfreq = 1;
printnode = 3;
rmax = 100; %maximum distance a node may travel in one cycle

%FEM CANTILEVER PARAMETERS
a_mag = 3.18e-4;
dx = 8; %microns
dy = 2; %microns
dz = 2; %microns
%tungsten "a" is 0.000274
dx = dx / a_mag;
dy = dy / a_mag;
dz = dz / a_mag;
DIST_SOURCE = 0.5 / a_mag;

mx = 20;
loading = 1;
vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

if (~exist('dt', 'var'))
    dt = dt0;
end

dt = min(dt, dt0);

disp('Initiliazing connectivity list. Please wait.');
[connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
[B, xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, Sbot, ...
        Sfront, Sback, gammat, gammau, gammaMixed, fixedDofs, freeDofs, ...
        w, h, d, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, MU, NU, loading);

disp('Creating surface mesh. Please wait.');
[TriangleCentroids, TriangleNormals, tri, Xb] = ...
    MeshSurfaceTriangulation(xnodes, Stop, Sbot, Sfront, Sback, Sleft, Sright);

disp('Calculating displacements from segments. Please wait.');
[uhat, ~, ~] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
    gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);
%fprintf('fend = %d, Ubar = %d, simTime = %d \n',fend,Ubar,simTime);

disp('Initiliazing motion...');
[rn, vn, dt, fn, fseg] = feval(integrator, rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
    rmax, rntol, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);

simTime = simTime + dt;
rnnew = [rn(:, 1:3) vn rn(:, 4)];
linksnew = links;
connectivitynew = connectivity;
linksinconnectnew = linksinconnect;
fsegnew = fseg;

if (doseparation)
    %spliting of nodes with 4 or more connections
    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
        separation(rnnew, linksnew, connectivitynew, linksinconnectnew, ...
        fsegnew, mobility, MU, NU, a, Ec, 2 * rann, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
end

%save restart.mat
if (docollision)
    %collision detection and handling
    [colliding_segments] = CollisionCheckerMex(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
        rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann);

    if colliding_segments == 1%scan and update dislocation structure.
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
            collision(rnnew, linksnew, connectivitynew, linksinconnectnew, ...
            fsegnew, rann, MU, NU, a, Ec, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
    end

end

if (doremesh)%do virtual re-meshing first
    %remeshing virtual dislocation structures
    if (dovirtmesh)
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = virtualmeshcoarsen_mex(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, DIST_SOURCE * 0.49, dx, MU, NU, a, Ec);
    end

    %remeshing internal dislocation structures
    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, mobility, doremesh, dovirtmesh, vertices, ...
        uhat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, TriangleNormals);
end

rnnew = rnnew(:, [1, 2, 3, 7]);

%% Compare displacement after remesh

disp('Comparing pre-meshed and post-meshed displacements...');
disp('Pre-meshed:');
[uhat, fend, Ubar] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
    gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);
fprintf('fend = %d, Ubar = %d, simTime = %d \n', fend, Ubar, simTime);

disp('Post-meshed');
[uhat_rem, fend_rem, Ubar_rem] = FEMcoupler(rnnew, linksnew, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
    gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);
fprintf('fend = %d, Ubar = %d, simTime = %d \n', fend_rem, Ubar_rem, simTime);

fprintf('\n \n Error in fend = %d, Error in Ubar = %d\n', norm(fend_rem - fend), norm(Ubar_rem - Ubar));
