%% Input file created by Daniel Hortelano Roig (14/02/2021) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sections are the following:
% 1) SIMULATION META-INFORMATION
% 2) MATERIAL PARAMETERS
% 3) FEM PARAMETERS
% 4) EasyDD PARAMETERS
%       - Parallelisation
%       - Meshing
%       - Time integrator
%       - Simulation time
%       - Plotting
%       - Mobility law
%       - Dislocation nodes (rn) and segments (links)
% 5) INFO & CHECKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear % Removes all items from workspace
dbstop if error % Applies debugging breakpoint instead of error
%dbstop if warning % Applies debugging breakpoint instead of warning


%% SIMULATION META-INFORMATION %%


%%% Simulation name

INPUTNAME = mfilename; % Input file name excluding '.m' extension
simName = append(INPUTNAME,'_',date,'_',char(timeofday(datetime)));

%%% Input paths & folders

% Functions designed for use in the input file:
if ~exist('./inputscripts/', 'dir')
       mkdir('inputscripts/');
end
addpath(genpath('inputscripts/'));

% Initial nodal network figures:
if ~exist('./figures/', 'dir')
       mkdir('figures/');
end


%% MATERIAL PARAMETERS %%

materialname = "alpha-Zr"; % alpha-Zirconium (temperature below ~900K)
CRYSTAL_STRUCTURE = 'hcp';

%%% Simulation units (SI)

lengthSI = 3.233e-10; % Length: [m]
pressureSI = 37.1E9; % Pressure: [Pa]
dragSI = 1.0e0; % Drag coefficient: [Pa s]
timeSI = dragSI/pressureSI; % Time: [s]
velocitySI = lengthSI/timeSI; % Velocity: [m/s]
forceSI = pressureSI * lengthSI^2; % Force: [Pa m^2]
nodalforceSI = pressureSI * lengthSI; % Nodal force: [Pa m]
temperatureSI = 1.0e0; % Temperature: [K]

%%% Crystallographic parameters

temperature = 600 / temperatureSI; % Units: [K]

% Lattice parameter magnitude:
length_inputunits = 1e-6; % Units: [m]
amag = lengthSI / length_inputunits; % Units: [microns]

% HCP a,c lattice parameters:
HCPa = 1.0; % Units: lengthSI
HCPc = 1.593; % Units: lengthSI

% Shear modulus:
pressure_inputunits = 1e6; % Units: [Pa]
mumag = pressureSI / pressure_inputunits; % Units: [MPa]
MU = 1; % Units: pressureSI

% Poisson's ratio:
NU = 0.32; % Units: unitless


%% FEM PARAMETERS %%


%%% Canilever beam discretisation

% Cantilever axes:
% ^z
% |     (y going into page)
% |
% |------>x

% Cantilever beam dimensions:
dx = 30/amag;
dy = 5/amag;
dz = 5/amag;

mx = 40; % Number of FEs along x-direction
my = max(round(mx * dy / dx),1);
mz = max(round(mx * dz / dx),1);

%%% Geometry and boundary conditions

simType = @cantileverBending;
loading = @displacementControl;
calculateLoading = @constantLoading;
calculateTractions = @calculateAnalyticTractions;

strain_dot = 1e4*timeSI; % Strain rate
u_dot = strain_dot*dx; % Displacement rate
f_dot = 0; % Force rate


%% DDLab PARAMETERS %%


%% Parallelisation & processes flags %%

% Parallelisation:
n_threads = 256; % Number of threads per GPU block
para_scheme = 1; % Parallelisation scheme: =0, none
%                                          =1, over dislocations
%                                          =2, over surface elements

% Flags:
CUDA_flag = false; % Compile CUDA code
CUDA_segseg_flag = false; % Parallelise segseg forces
doremesh = 1;
docollision = 1;
doseparation = 1;
dovirtmesh = 1;

%% Integrator %%

integrator = @int_trapezoid; % Time integrator function

a = 90; % Isotropic spreading radius
Ec = (MU/(4*pi))*log(a/0.1); % Segment core energy per unit length
rann = 0.5*a; % Segment annihilation dist
rntol = 0.5*rann; % Tolerance for max dist travelled by segment
rmax = 2.5*a*(sqrt(3)*2); % Max dist a node may move during time incr
mindist = 2*rann;

%% Meshing %%

maxconnections = 8;
lmax = rmax;
lmin = lmax/2.5;
areamin = 2*rntol*lmax;
areamax = 2*areamin + lmax^2*sqrt(3)/8;

%% Simulation time %%

% Permitted absolute time increment:
dt0 = 1E-7/timeSI; % Maximum
dt = 0.5E-9/timeSI; % Minimum

% Absolute time duration of simulation:
totalSimTime = 1E-5/timeSI;

% Initial time of simulation:
simTime = 0/timeSI;

% Initial time-step:
curstep = 0;

%% Plotting %%

% Plotting parameters:
plotFreq = 0; % =0,1: plot at each time-step, =Inf: never plot
viewangle = [0,0]; % view(azimuth,elevation)
saveFreq = 100; % Frequency to save workspace variables per time-step
printfreq = 100; % =0,1: print info at each time-step, =Inf: never print

%% Mobility %%

%%% Referential information

mobility = @mobhcp0;

% HCP lattice vectors:
refHCPa1 = HCPa * [-1/2 sqrt(3)/2 0];
refHCPa2 = HCPa * [-1/2 -sqrt(3)/2 0];
refHCPa3 = HCPa * [1 0 0];
refHCPa4 = HCPc * [0 0 1];

% Number of Burgers vector indices of each type:
numburgsA = 3; % Type <a>
numburgsCA = 6; % Type <c+a>
numburgsC = 1; % Type <c>
numburgs = numburgsA + numburgsCA + numburgsC;

numplanetypes = 5; % Number of types of glissile planes
% NOTE: slip plane types are indexed as:
% Ba = 1, Pr = 2, PyI = 3, PyII = 4, sP = 5, Sessile = 6


%%% Drag coefficients

% Plane Type -- Burgers Vector Index -- Orientation

glidecoeffsSTR = struct; % Stores mobility coefficients

% Glide drag coefficients for <a> segments (Burgers vector indices 1-3):
glidecoeffsSTR.Ba1Screw = 180e-6; glidecoeffsSTR.Ba1Edge = 90e-6;
glidecoeffsSTR.Pr1Screw = 85e-6; glidecoeffsSTR.Pr1Edge = 39e-6;
glidecoeffsSTR.PyI1Screw = 85e-6; glidecoeffsSTR.PyI1Edge = 39e-6;
glidecoeffsSTR.PyII1Screw = 85e-6; glidecoeffsSTR.PyII1Edge = 39e-6;

% Glide drag coefficients for <c+a> segments (4-9):
glidecoeffsSTR.Pr4Screw = 500e-6; glidecoeffsSTR.Pr4Edge = 100e-6;
glidecoeffsSTR.PyI4Screw = 500e-6; glidecoeffsSTR.PyI4Edge = 100e-6;
glidecoeffsSTR.PyII4Screw = 500e-6; glidecoeffsSTR.PyII4Edge = 100e-6;
glidecoeffsSTR.sP4Screw = 500e-6; glidecoeffsSTR.sP4Edge = 100e-6;

% Glide drag coefficients for <c> segments (10):
glidecoeffsSTR.Pr10Screw = 85e-6; glidecoeffsSTR.Pr10Edge = 39e-6;

% Glide drag for junction segments:
glidecoeffsSTR.dragsessile = 1e1;

% Organize glide drag coefficients in matrix:
glidecoefficients = OrganizeGlideDragCoefficients( ... % size (numburgs,numplanetypes,2)
numburgsA,numburgsCA,numburgsC,numburgs,numplanetypes, ... % Organization-related inputs
glidecoeffsSTR);                                        % Glide coefficient structure

% Climb drag:
dragclimb = glidecoeffsSTR.dragsessile;

% Line drag:
dragline = min(1e-1 * min(glidecoefficients,[],'all'), 1e-1 / areamin);


%%% Glissile slip systems

slipsystems = SlipSystemsDef(numburgs, numplanetypes);

%%% Transformation of crystal grain slip systems

% Convert slip systems to cartesian:
slipsystemscart = SlipSystemsToCartesian(slipsystems, refHCPa1,refHCPa2,refHCPa3,refHCPa4);

% Prepare to rotate slip systems:
ssangx = -90/180 * pi;
ssangy = -45/180 * pi;
%ssangz = 0/180 * pi;
ssRx = [1 0 0; 0 cos(ssangx) -sin(ssangx); 0 sin(ssangx) cos(ssangx)];
ssRy = [cos(ssangy) 0 sin(ssangy); 0 1 0; -sin(ssangy) 0 cos(ssangy)];
%ssRz = [cos(ssangz) -sin(ssangz) 0; sin(ssangz) cos(ssangz) 0; 0 0 1];

ssrotmatrix = ssRy * ssRx; % Rotation matrix

% Rotate reference Burgers, normal, and HCP lattice vectors:
[slipsystemscartrot,HCPa1,HCPa2,HCPa3,HCPa4] = ...
    RotateCartesianSlipSystems(slipsystemscart, refHCPa1,refHCPa1,refHCPa1,refHCPa1, ssrotmatrix);

% Store rotated reference cartesian Burgers vectors:
burgscart = cat(1,slipsystemscartrot{:,1,1});


%%% Critical resolved shear stress (CRSS)

enableCRSS = 0; % If = 0, all CRSS values are set to zero

crsscoeffsSTR = struct; % Stores CRSS coefficients

% Macro CRSS values (see [3])
% <a>:
crsscoeffsSTR.a_Ba = 204e6 / pressureSI; % 204 MPa
crsscoeffsSTR.a_Pr = 153e6 / pressureSI; % 153 MPa
crsscoeffsSTR.a_PyI = crsscoeffsSTR.a_Pr;
crsscoeffsSTR.a_PyII = crsscoeffsSTR.a_Pr;
% <c+a>:
crsscoeffsSTR.ca_PyI = 532e6 / pressureSI; % 532 MPa
crsscoeffsSTR.ca_PyII = crsscoeffsSTR.ca_PyI;
crsscoeffsSTR.ca_sP = crsscoeffsSTR.ca_PyI;
crsscoeffsSTR.ca_Pr = crsscoeffsSTR.ca_PyI;
% <c>:
crsscoeffsSTR.c_Pr = crsscoeffsSTR.a_Pr;
% Sessile:
crsssessile = 10e9 / pressureSI; % 10 GPa
crsscoeffsSTR.sessile = crsssessile;

% Organize CRSS coefficients in matrix:
crsscoefficients = OrganizeCoefficients( ... % size (numburgs,numplanetypes + 1)
numburgsA,numburgsCA,numburgsC,numburgs,numplanetypes, ... % Organization-related inputs
crsscoeffsSTR);                                               % Coefficients structure

if enableCRSS ~= 1
    enableCRSS = 0;
    crsscoefficients(:) = 0;
end


%%% Store mobility information inside structure

Bcoeff = struct;

Bcoeff.slipsystemsref = slipsystemscartrot; % Cell: size (numburgs,numplanetypes,2)
Bcoeff.burgscart = burgscart; % Matrix: size (numburgs,3)
Bcoeff.HCPa1 = HCPa1;
Bcoeff.HCPa2 = HCPa2;
Bcoeff.HCPa3 = HCPa3;
Bcoeff.HCPa4 = HCPa4;
Bcoeff.dragline = dragline;
Bcoeff.dragclimb = dragclimb;
Bcoeff.glidecoefficients = glidecoefficients; % Matrix: size (numburgs,numplanetypes,2)
Bcoeff.dragsessile = glidecoeffsSTR.dragsessile;
Bcoeff.crsscoefficients = crsscoefficients; % Matrix: size (numburgs,numplanetypes+1)
Bcoeff.crsssessile = crsssessile;

%% Initial nodes and links %%

%%% Nodes %%%

%%%
% Frank-Read source:
rnfrref =   [0  -0.5   0 7;
             0   0     0 0;
             0   0.5   0 7];
% Neutron irradiation-induced prismatic loop (square):
rnvacref =  [0   0.5  0.5 0; % 8 nodes
             0   0.5  0   0;
             0   0.5 -0.5 0;
             0   0   -0.5 0;
             0  -0.5 -0.5 0;
             0  -0.5  0   0;
             0  -0.5  0.5 0;
             0   0    0.5 0];
%%%

% Frank-Read source: 2Ba
b_frref = slipsystemscartrot{2,1,1};
b_frmag = norm(b_frref);
b_fr = b_frmag * (b_frref ./ norm(b_frref));
b_frunit = b_fr / norm(b_fr);
n_fr = slipsystemscartrot{2,2,2};

% Vacancy loop: 2Pr
b_vref = slipsystemscartrot{2,2,1};
b_vmag = -norm(b_vref);
b_v = b_vmag * (b_vref ./ norm(b_vref));
b_vunit = b_v / norm(b_v);
n_1v = slipsystemscartrot{2,1,2}; % Ba plane
n_2v = slipsystemscartrot{2,2,2}; % Pr plane

%%% Scale and rotate nodes

% Prepare rotation:
%rnangx = 0/180 * pi;
rnangy = -45/180 * pi;
%rnangz = 0/180 * pi;
%rnRx = [1 0 0; 0 cos(rnangx) -sin(rnangx); 0 sin(rnangx) cos(rnangx)];
rnRy = [cos(rnangy) 0 sin(rnangy); 0 1 0; -sin(rnangy) 0 cos(rnangy)];
%rnRz = [cos(rnangz) -sin(rnangz) 0; sin(rnangz) cos(rnangz) 0; 0 0 1];

% FR source
rnfrscaled = rnfrref;
% Scale Frank-Read source:
frlen = 350e-9 / lengthSI; % Frank-Read source length
rnfrscaled(:,1:3) = rnfrscaled(:,1:3) .* frlen;

% Vacancy loop
rnvacscaled = rnvacref;
% Scale vacancy loop:
vaclen = 40 * 10e-9 / lengthSI; % Square side length
rnvacscaled(:,1:3) = rnvacscaled(:,1:3) .* vaclen;
% Rotate vacancy loop:
rnrotmatrix = rnRy;
rnvaccent = sum(rnvacscaled(:,1:3),1) / size(rnvacscaled,1);
rnvacscaled(:,1:3) = (rnrotmatrix * (rnvacscaled(:,1:3) - rnvaccent)')' + rnvaccent;

%%% Translate nodes

% Prepare translation:
centpos_lo = [dx/4 dy/2 dz/10];
centpos_hi = [dx dy dz] - centpos_lo;

% Higher FR source:
rnfr_hi = rnfrscaled;
rnfr_hi(:,1:3) = rnfr_hi(:,1:3) + centpos_hi;

% Lower FR source:
rnfr_lo = rnfrscaled;
rnfr_lo(:,1:3) = rnfr_lo(:,1:3) + centpos_lo;

% High vacancy loop:
rnvac_hi = rnvacscaled;
rnvac_hi(:,1:3) = rnvac_hi(:,1:3) + centpos_hi;

% Lower vacancy loop:
rnvac_lo = rnvacscaled;
rnvac_lo(:,1:3) = rnvac_lo(:,1:3) + centpos_lo;

% Higher node networks

rnv1_hi = rnvac_hi;
rnv1_hi(:,1:3) = rnvac_hi(:,1:3) + (dz*sqrt(2)*1/8) * (b_vunit);

rnv2_hi = rnvac_hi;
scalecentrnv3_hi = (dz/2 - centpos_hi(3))/b_vunit(3);
rnv2_hi(:,1:3) = rnvac_hi(:,1:3) + scalecentrnv3_hi * (b_vunit);
%rnv3(:,1:3) = rnvac(:,1:3) + 25 * (110e-9 / lengthSI) * (-b_vunit);

rnv3_hi = rnvac_hi;
rnv3_hi(:,1:3) = rnvac_hi(:,1:3) + (dz*sqrt(2)*6/8) * (b_vunit);

% Lower node networks

% rnv1 = rnvac;
% rnv1(:,1:3) = rnvac(:,1:3) + 40 * (10e-9 / lengthSI) * (-b_vunit);

rnv1_lo = rnvac_lo;
rnv1_lo(:,1:3) = rnvac_lo(:,1:3) + (dz*sqrt(2)*1/8) * (-b_vunit);

rnv2_lo = rnvac_lo;
scalecentrnv3_lo = (dz/2 - centpos_lo(3))/b_vunit(3);
rnv2_lo(:,1:3) = rnvac_lo(:,1:3) + scalecentrnv3_lo * (b_vunit);
%rnv3(:,1:3) = rnvac(:,1:3) + 25 * (110e-9 / lengthSI) * (-b_vunit);

rnv3_lo = rnvac_lo;
rnv3_lo(:,1:3) = rnvac_lo(:,1:3) + (dz*sqrt(2)*6/8) * (-b_vunit);

% rnv5 = rnvac;
% rnv5(:,1:3) = rnvac(:,1:3) + 25 * (210e-9 / lengthSI) * (-b_vunit);

%%% Collect nodes

% Higher:
rn = [      rnfr_hi;    rnv1_hi; rnv2_hi; rnv3_hi];
% Lower:
rn = [rn;   rnfr_lo;    rnv1_lo; rnv2_lo; rnv3_lo];


%%% Links %%%

% Reference
linksfr =    [3  2   b_fr n_fr; % FR source
              2  1   b_fr n_fr];
linksvac =   [1  2   b_v n_1v; % Vac loop
              2  3   b_v n_1v;
              3  4   b_v n_2v;
              4  5   b_v n_2v;
              5  6   b_v n_1v;
              6  7   b_v n_1v;
              7  8   b_v n_2v;
              8  1   b_v n_2v];
%

%%% Higher

linksfr_hi = linksfr;

linksv1_hi = linksvac;
v1hi_idx = max(linksfr_hi(:,1:2),[],'all');
linksv1_hi(:,1:2) = linksvac(:,1:2) + v1hi_idx;

linksv2_hi = linksvac;
v2hi_idx = max(linksv1_hi(:,1:2),[],'all');
linksv2_hi(:,1:2) = linksvac(:,1:2) + v2hi_idx;

linksv3_hi = linksvac;
v3hi_idx = max(linksv2_hi(:,1:2),[],'all');
linksv3_hi(:,1:2) = linksvac(:,1:2) + v3hi_idx;

% linksv4_hi = linksvac;
% linksv4_hi(:,1:2) = linksvac(:,1:2) + 27;

% linksv5_hi = linksvac;
% linksv5_hi(:,1:2) = linksvac(:,1:2) + 35;

%%% Lower

linksfr_lo = linksfr;
frlo_idx = max(linksv3_hi(:,1:2),[],'all');
linksfr_lo(:,1:2) = linksfr(:,1:2) + frlo_idx;

linksv1_lo = linksvac;
v1lo_idx = max(linksfr_lo(:,1:2),[],'all');
linksv1_lo(:,1:2) = linksvac(:,1:2) + v1lo_idx;

linksv2_lo = linksvac;
v2lo_idx = max(linksv1_lo(:,1:2),[],'all');
linksv2_lo(:,1:2) = linksvac(:,1:2) + v2lo_idx;

linksv3_lo = linksvac;
v3lo_idx = max(linksv2_lo(:,1:2),[],'all');
linksv3_lo(:,1:2) = linksvac(:,1:2) + v3lo_idx;

% linksv4_lo = linksvac;
% linksv4_lo(:,1:2) = linksvac(:,1:2) + 27;

% linksv5_lo = linksvac;
% linksv5_lo(:,1:2) = linksvac(:,1:2) + 35;

%%% Collect links

% Higher:
links = [           linksfr_hi;     linksv1_hi; linksv2_hi; linksv3_hi];
% Lower:
links = [links;     linksfr_lo;     linksv1_lo; linksv2_lo; linksv3_lo];

%% INFO & CHECKS %%


run PrintInformation.m % Print simulation input information

run PlotInitialNodalNetwork.m % Save figure of nodal network
% To open saved figure, run: openfig(['figures/' INPUTNAME '.fig'])

run InputConsistencyCheck.m % Checks and warnings

run OrganizeInStructs.m % Store inessential variables in structures


%% REFERENCES %%


% [1]: Tarleton 2015: http://dx.doi.org/10.1016/j.actamat.2015.01.030
% [2]: Drouet 2014: https://doi.org/10.1016/j.jnucmat.2013.11.049
% [3]: Gong 2015: https://doi.org/10.1016/j.actamat.2015.06.020
% [4]: Introduction to Dislocations. Hull and Bacon, (2011)
% [5]: Khater 2010: https://doi.org/10.1016/j.actamat.2010.01.028
