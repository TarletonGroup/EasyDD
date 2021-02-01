%% Dual time-step input parameters

%% Preliminaries (DHT)

pause on
format compact

%% Initial globals

global Bscrew Bedge Beclimb Bline

%% Initialise new globals (DHT)

global k eV
global MU_SI a_SI Bglide_SI
global atomVol vacFormEnergy vacMigrEnergy diffCoeffPreExpo
global radiusCore radiusInfinity
global temperature cInfinityRatio
global osmoticForce
global inclusions impenetrable

%% Generate dislocation loop (DHT)

[rn, links] = generateCircular100ShearLoop(50, 5);

%% Inclusion input parameters (DHT)

%  inclusions(:, 1) - radius
%  inclusions(:, 2:4) - position
%  inclusions(:, 5) - volumetric misfit

%inclusions = []; % Empty array for no inclusions
inclusions = [20 0 0 0 -0.1]; 
%inclusions = [20  50  50 0 -0.1;
%              20  40 -50 0 -0.1;
%              20 -50  40 0 -0.1;
%              20 -40 -40 0 -0.1];

impenetrable = 1; % == 1 for impenetrable inclusions,
                  % == anything else for penetrable ones

%% Define new globals (DHT)

% SI Constants
k = 1.3806488e-23;    % Boltzmann's constant k
eV = 1.602176565e-19; % 1 eV

% SI properties of bcc iron for use converting between SI and DDLab units
MU_SI = 82e9;      % Shear modulus, 82 GPa
a_SI = 2.8665e-10; % Lattice parameter, 2.886 Angstroms
Bglide_SI = 1e-4;  % Glide drag factor, 1e-4 Pa s (guestimate)
                   % so [T] = Bglide_SI/MU_SI = 1.22e-15 s i.e. 1.22 fs

% Other SI properties of bbc iron
atomVol = 11.8e-30;        % Atomic volume, 11.8 cubic Angstroms
vacFormEnergy = 2 * eV;    % Vacancy formation energy, 2 eV
vacMigrEnergy = 0.6 * eV;  % Vacancy migration energy, 0.6 eV
diffCoeffPreExpo = 2e-4;   % Vacancy diffusion coefficient pre-exponent

% Simulation parameters
radiusCore = 5;         % Dislocation core radius in Burgers vectors
radiusInfinity = 10000; % Outer cut-off radius in Burgers vectors

%% Define parameters

MU = 1;
NU = 0.305;
maxconnections = 8;
lmax = 100; % Default value = 1000
lmin = 5; % Default value = 200
areamin = lmin * lmin * sin(60 / 180 * pi) * 0.5; 
areamax = 20 * areamin;
a = lmin / sqrt(3) * 0.5;
Ec = MU / (4 * pi) * log(a / 0.1);
totalsteps = 1001; % Put 1 more than you actually want!
dt0 = 1e15; % Maximum time step = 1.22 s

mobility = 'mobbcc0';

%Drag (Mobility) parameters
Bscrew = 1;
Bedge = 1;

%% Temperature dependent climb (DHT)

temperature = 973; % Temperature in Kelvins, 700 degress Celcius
cInfinityRatio = 1;
Beclimb = calcClimbDrag;
osmoticForce = calcOsmoticForce;
Bline = 1e-4*min(Bscrew,Bedge); % Default
%Bline = 1e-4*Beclimb; % For climb only motion

%% Time-step switching (DHT)

doGlide = 1; % Start with glide
minStrainRatio = 1e-3;
maxclimbsteps = 1;

%% Movie frame recording (DHT)

doMovie = 0; % On == 1, Off == anything else
frameCounter = 0; % Initialise
frameFreq = 1;

%% Define more parameters

integrator = 'int_eulerforward';
rann = 0.5 * a;       
rntol = 0.5 * rann;
doremesh = 1;
docollision = 1;
doseparation = 1;
plotfreq = 1;      
plim = 60; % 60 for most loop and inclusion sims
           % 200 when necessary (certain glide + climb sims)

%% Applied stress (DHT)

%appliedstress = 1e-3 .* [2 0 1; 0 2 -1; 1 -1 0]; % Default value
appliedstress = zeros(3,3);                      % Zero applied stress
%appliedstress = 1e-3 .* [0 5 0; 5 0 0; 0 0 0];   % 12 Shear load
%appliedstress = 1e-3 .* [0 0 5; 0 0 0; 5 0 0];   % 13 Shear load
%appliedstress = 1e-3 .* [0 0 0; 0 0 5; 0 5 0];   % 23 Shear load
%appliedstress = 1e-3 .* [10 0 0; 0 0 0; 0 0 0];  % Uniaxial 11
                                                  % tension/compression
                                                  % In units of MU = 82 GPa
                                               
%% View angle (DHT)

%viewangle = [45 -45]; % Default value
%viewangle = [90 0];   % Down x axis for b=[100] prismatic loop
viewangle = [0 0];    % Down y axis for rubber-banding b=[100] prismatic loop
%viewangle = [0 90];   % Down z axis
%viewangle = [45 45];  % 45 degrees to both x and y axes,
                       % and elevated 45 degrees out of xy plane
                       % for rubber-banding b=[100] prismatic loop
%viewangle = [45 0];   % 45 degrees to both x and y axes 
                       % for rubber-banding b=[100] prismatic loop
                      
%% Define final parameters

printfreq = 1;      
printnode = 3;
rmax = 1; % Default value = 100

% simulation box size (for paradis)
L = 40000;
minSideX = -L/2;
minSideY = -L/2;
minSideZ = -L/2;
maxSideX =  L/2;
maxSideY =  L/2;
maxSideZ =  L/2;

% boundary conditions (for paradis)
xBoundType = 1; % free
yBoundType = 1; % free
zBoundType = 1; % free

% output and communication settings
paradis_dir = '../..';
paradis_input_dir = strcat(paradis_dir, '/Runs');
paradis_output_dir = strcat('Outputs/frank_read_results'); 
