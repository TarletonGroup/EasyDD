clear all;
close all;
CRYSTAL_STRUCTURE = 'fcc';

slipPlanes = [
        1.0 1.0 1.0;
        1.0 1.0 1.0;
        1.0 1.0 1.0;
        -1.0 1.0 1.0;
        -1.0 1.0 1.0;
        -1.0 1.0 1.0;
        1.0 -1.0 1.0;
        1.0 -1.0 1.0;
        1.0 -1.0 1.0;
        1.0 1.0 -1.0;
        1.0 1.0 -1.0;
        1.0 1.0 -1.0
        ];

bVec = [
    1.0 -1.0 0.0;
    1.0 0.0 -1.0;
    0.0 1.0 -1.0;
    1.0 1.0 0.0;
    1.0 0.0 1.0;
    0.0 1.0 -1.0;
    1.0 1.0 0.0;
    0.0 1.0 1.0;
    1.0 0.0 -1.0;
    1.0 0.0 1.0;
    0.0 1.0 1.0;
    1.0 -1.0 0.0
    ];

% Values from https://www.azom.com/properties.aspx?ArticleID=2193
% Nickel lattice parameter: 3.499 Angstroms
amag = 3.499 * 1e-4; % microns * burgers vector

% Nickel shear modulus: 72 - 86 GPa
mumag = 79e3; % MPa
MU = 1.0;

% Nickel poisson's ratio 0.305 - 0.315.
NU = 0.31;

% x = <100>, y = <010>, z = <001>
% FE domain dimensions: x/amag := x microns.
% The size effect is given by the size disparity between dislocations and
% volume. We can either change the size of the volume or the size of the
% dislocations. 2^3, 4^3, 8^3, 16^3, 32^3 micron domain.
dx = 32 / amag; % This gives microns.
dy = dx;
dz = dx;
mx = 20;
my = 20;
mz = 20;

% Dislocation segment lengths and integration tolerances.
segLen = 0.5 / amag; % Source segment length.
lmin = segLen / 10; % Minimum allowed segment length.
lmax = segLen / 5; % Maximum allowed segment length.
a = lmin / 20; % Dislocation core radius (for non-singular stresses).
rann = lmin;%lmin/2; % Collision distance between two dislocations.
rntol = lmin/2; % Error tolerance between one step to the next.
rmax = lmin/2; % Maximum change in position between one step to the next.

vertices = [0, 0, 0; ...
            dx, 0, 0; ...
            0, dy, 0; ...
            dx, dy, 0; ...
            0, 0, dz; ...
            dx, 0, dz; ...
            0, dy, dz; ...
            dx, dy, dz];

% Loading rate.
% time_real = time_EasyDD / mumag / 1e6;
% If we were to use the same loading rate we'd be loading at
% 5e-3/mumag/1e6 \approx 6.3291e-14. This is untractable. Instead we scale
% with a heuristic and using beam theory.
% timeUnit = 5e-3 * mumag * 1e6;
timeUnit = 5e-3*mumag * 1e6;
u_dot = dx / timeUnit;
% 1 := tensile, -1 := compressive
sign_u_dot = 1;

% Sets boundary conditions and simulation type.
loading = @displacementControlMicropillarTensile;
simType = @micropillarTensile;

% Precomputed FCC loops with 4 sides.
run fccLoops
% Scaling the Burgers vectors they're normalised.
% Prismatic loops.
prismbVec(:, :) = prismbVec(:, :) / max(abs(prismbVec(1, :)));
prismbVec(:, :) = prismbVec(:, :) * norm(prismbVec(1, :));
% Shear loops.
shearbVec(:, :) = shearbVec(:, :) / max(abs(shearbVec(1, :)));
shearbVec(:, :) = shearbVec(:, :) * norm(shearbVec(1, :));

% Range where the loops will be centred. In this case in the middle of the
% domain. If we set the loop right in the middle with the same number of
% FE nodes, the displacements diverge so we displace the centre scaled
% according to the volume in number of elements times the minimum value
% double precision can distinguish between two numbers. This as is close to
% the minimum displacement away from the centre that we can get without
% running into numerical under or overflow.
xmin = (0.5 + (mx * my * mz)^3 * eps) * dx;
xmax = (0.5 + (mx * my * mz)^3 * eps) * dx;
ymin = 0.5 * dy;
ymax = 0.5 * dy;
zmin = 0.5 * dz;
zmax = 0.5 * dz;

distRange = [xmin ymin zmin; xmax ymax zmax];
displacement = distRange(1, :) + (distRange(2, :) - distRange(1, :)) .* rand(12, 3);
links = [];
rn = [];

% We simulate a FR source by pinning all nodes of the loop but one. We also
% only want to see the effect a single dislocation in the middle of the
% domain has on the behaviour.
for i = 1:1
    idx = (i - 1) * 8;
    links = [links; (prismLinks((1:8) + idx, :) + idx) prismbVec((1:8) + idx, :) prismSlipPlane((1:8) + idx, :)];
    displacedCoord = prismCoord((1:8) + idx, :) * segLen + displacement(i, :);
    rn = [rn; displacedCoord [7; 7; 7; 7; 0; 7; 7; 7]];
end

plotnodes(rn, links, dx, vertices);
dt0 = timeUnit;
dtMin = 10 * eps; % Minimum allowed timestep.
totalSimTime = timeUnit * 1e4; % Total simulated time.
mobility = @mobfcc0; % Mobility law.
saveFreq = 200; % Saving frequency (steps between saves).
plotFreq = 10; % Plotting frequency (steps between plotting).

% Set scaling factors for plotting the displacements and force.
plotArgs = struct("factDisp", amag, "factForce", amag^2*mumag);
% What plots to show.
plotFlags = struct('nodes', true, 'secondary', true);

% Mobility in different directions. [Pa][s]/[b]
Bcoeff = struct('screw', 1, 'edge', 1, 'climb', 1e6, 'line', 1e-4);

% Use analytic tractions.
calculateTractions = @calculateAnalyticTractions;

simName = date;
simName = strcat(simName, sprintf('_Tensile_Ni_100_size_effect_%d', dx / segLen));
