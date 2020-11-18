%% Input file created by Daniel Hortelano Roig (17/07/2020) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sections are the following:
% 0) SIMULATION META-INFORMATION
% 1) MATERIAL CONSTANTS
% 2) SOURCE GENERATION PARAMETERS
% 3) FEM PARAMETERS
% 4) DDLab PARAMETERS
%       - Meshing
%       - Integrator
%       - Simulation time
%       - Plotting
%       - Mobility module
%       - Dislocation nodes (rn) and segments (links)
% 5) PRINT INFORMATION
% 6) REFERENCES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
dbstop if error % Applies debugging breakpoint instead of error





%% SIMULATION META-INFORMATION %%


%%% Input name

INPUTFILENAME = mfilename;
INPUTNAME = INPUTFILENAME(7:end);

%%% MATLAB version year

% Release version of MATLAB client currently running:
MATLABrelease = ver('MATLAB').Release;

%%% Add EasyDD subfolders to MATLAB path

% Main 'EasyDD' folder:
mainpath = '../../';
addpath(mainpath);

% Subprocesses folders:
addpath([mainpath 'Mobility']); % Mobility laws
addpath([mainpath 'Integrator']); % Time integrators

% Functions designed for use in the input file:
addpath([mainpath 'Inputs/HCP/InputFunctions']);

% Input figures:
%addpath([mainpath 'Inputs/HCP/Figures']);

% Output figures:
frameext = '.jpg'; % File extension of saved frames
movieext = '.mp4'; % File extension of saved movies
outputspath = append('Outputs/',INPUTNAME,'/'); % Path to saved figures
if ~exist([mainpath outputspath], 'dir')
       mkdir([mainpath outputspath]);
end
outputsframespath = append(outputspath,'Frames/');
if ~exist([mainpath outputsframespath], 'dir')
       mkdir([mainpath outputsframespath]);
end

%%% Simulation units in SI

lengthSI = 3.233e-10; % Length: [m]
pressureSI = 37.1E9; % Pressure: [Pa]
dragSI = 1.0e0; % Drag coefficient: [Pa s]
timeSI = dragSI/pressureSI; % Time: [s]
velocitySI = lengthSI/timeSI; % Velocity: [m/s]
forceSI = pressureSI * lengthSI^2; % Force: [Pa m^2]
nodalforceSI = pressureSI * lengthSI; % Nodal force: [Pa m]
temperatureSI = 1.0e0; % Temperature: [K]





%% MATERIAL CONSTANTS %%


% Material name:
materialname = "alpha-Zr"; % alpha-Zirconium (temperature below ~900K)

% Temperature:
temperature = 600; % Units: [K]

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





%% SOURCE GENERATION PARAMETERS %%


DIST_SOURCE = 0.5/amag;





%% FEM PARAMETERS %%


% Cantilever Axes:
% ^z
% |     (y going into page)
% |
% |------>x

% Cantilever beam dimensions:
dx = 30/amag; % NOTE: 1/amag = 1 micron
dy = 5/amag;
dz = 5/amag;

% Number of finite elements along beam x-direction:
mx = 40;

% Vertices of cantilever beam:
vertices = [ 0,  0,  0;...
            dx,  0,  0;...
             0, dy,  0;...
            dx, dy,  0;...
             0,  0, dz;...
            dx,  0, dz;...
             0, dy, dz;...
            dx, dy, dz];

% Linear loading:
loading = 1; % 1 = on, 0 = off
straindot = 10e5*timeSI; % Strain rate
Udot = straindot*dx; % Loading rate





%% DDLab PARAMETERS %%


%%%%% Integrator %%%%%

% Time integrator function:
integrator = 'int_trapezoid';

% Time integration variables:
a = 90; % Isotropic spreading radius
Ec = (MU/(4*pi))*log(a/0.1); % Segment core energy per unit length
rann = 0.5*a; % Segment annihilation dist
rntol = 0.5*rann; % Tolerance for max dist travelled by segment
rmax = 2.5*a*(sqrt(3)*2); % Max dist a node may move during time incr





%%%%% Meshing %%%%%

maxconnections = 8;
lmax = rmax;
lmin = lmax/2.5;
areamin = 2*rntol*lmax;
areamax = 2*areamin + lmax^2*sqrt(3)/8;
doremesh = 1;
docollision = 1;
doseparation = 1;
dovirtmesh = 1;





%%%%% Simulation time %%%%%

% Permitted absolute time increment:
dt0 = 1e-5/timeSI; % Maximum
dt = 0.5e-9/timeSI; % Minimum

% Absolute time duration of simulation:
totalSimTime = 1/timeSI;

% Initial time of simulation:
initSimTime = 0/timeSI;

% Initial time-step:
curstep = 0;





%%%%% Plotting %%%%%

% Plotting parameters:
plotfreq = 0; % =1: plot at each time-step, =Inf: never plot
pauseplot = 1e-3; % Time paused after plot is generated (Units: [s])
plim = 12/amag; % Does nothing...
viewangle = [0,0]; % view(azimuth,elevation)
printfreq = Inf; % =0: print info at each time-step, =Inf: never print
printnode = 1;

% Frame plotting:
doframeplot = 1; % Frame recording: =1: on, =0: off
alwaysframeplot = 1; % Frame record each time-step
intSimTime = 0/timeSI; % Initial counting time before first frame
numframes = 1000; % Total number of recorded frames
dtplot = floor(totalSimTime/numframes); % Time increment between frames
viewangleframes = [-45,-5];
%viewaxisframes = [dx*5/32 dx*7/32 dy*12/32 dy*20/32 dz*12/32 dz*20/32];
%viewaxisframes = [dx*2/16 dx*3/16 dy*3/8 dy*5/8 dz*1/4 dz*9/16];
%viewaxisframes = [dx*2/16 dx*3/16 dy*3/8 dy*5/8 dz*1/4 dz*9/16];
%viewaxisframes = [dx*2/16 dx*3/16 dy*1/4 dy*3/4 dz*1/4 dz*2/3]; % axis(viewaxisframes)





%%%%% Mobility %%%%%

%%% Referential information

mobility = 'mobhcp0_DHR';

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
glidecoeffsstruct = struct;

% Glide drag coefficients for <a> segments (Burgers vector indices 1-3):
glidecoeffsstruct.Ba1Screw = 180e-6; glidecoeffsstruct.Ba1Edge = 90e-6;
glidecoeffsstruct.Pr1Screw = 85e-6; glidecoeffsstruct.Pr1Edge = 39e-6;
glidecoeffsstruct.PyI1Screw = 85e-6; glidecoeffsstruct.PyI1Edge = 39e-6;
glidecoeffsstruct.PyII1Screw = 85e-6; glidecoeffsstruct.PyII1Edge = 39e-6;

% Glide drag coefficients for <c+a> segments (4-9):
glidecoeffsstruct.Pr4Screw = 500e-6; glidecoeffsstruct.Pr4Edge = 100e-6;
glidecoeffsstruct.PyI4Screw = 500e-6; glidecoeffsstruct.PyI4Edge = 100e-6;
glidecoeffsstruct.PyII4Screw = 500e-6; glidecoeffsstruct.PyII4Edge = 100e-6;
glidecoeffsstruct.sP4Screw = 500e-6; glidecoeffsstruct.sP4Edge = 100e-6;

% Glide drag coefficients for <c> segments (10):
glidecoeffsstruct.Pr10Screw = 85e-6; glidecoeffsstruct.Pr10Edge = 39e-6;

% Glide drag for junction segments:
glidecoeffsstruct.dragsessile = 1e1;

% Organize glide drag coefficients in matrix:
glidecoefficients = OrganizeGlideDragCoefficients( ... % size (numburgs,numplanetypes,2)
numburgsA,numburgsCA,numburgsC,numburgs,numplanetypes, ... % Organization-related inputs
glidecoeffsstruct);                                        % Glide coefficient structure

% Climb drag:
dragclimb = glidecoeffsstruct.dragsessile;

% Line drag:
dragline = min(1e-1 * min(glidecoefficients,[],'all'), 1e-1 / areamin);


%%% Glissile slip systems

slipsystemsref = cell(numburgs,numplanetypes,2);

% Type <a> Burgers vectors
% [b1](p)
slipsystemsref(1,1:5,1) = {1/3 * [2,-1,-1,0]};         % b1
slipsystemsref(1,:,2) = {[0,0,0,1]; ...  % 1Ba
					  [0,1,-1,0]; ... % 1Pr
					  [0,1,-1,1]; ... % 1PyI
					  [0,-1,1,1]; ... % 1PyII
					  []; ...         % 1sP (empty)
					  };
% [b2](p)
slipsystemsref(2,1:5,1) = {1/3 * [-1,-1,2,0]};         % b2
slipsystemsref(2,:,2) = {[0,0,0,1]; ...  % 2Ba
					  [1,-1,0,0]; ... % 2Pr
					  [1,-1,0,1]; ... % 2PyI
					  [-1,1,0,1]; ... % 2PyII
					  []; ...         % 2sP (empty)
					  };
% [b3](p)
slipsystemsref(3,1:5,1) = {1/3 * [-1,2,-1,0]};         % b3
slipsystemsref(3,:,2) = {[0,0,0,1]; ...  % 3Ba
					  [-1,0,1,0]; ... % 3Pr
					  [-1,0,1,1]; ... % 3PyI
					  [1,0,-1,1]; ... % 3PyII
					  []; ... 	      % 3sP (empty)
					  };

% Type <c+a> Burgers vectors
% [b4](p)
slipsystemsref(4,1:5,1) = {1/3 * [2,-1,-1,3]};         % b4
slipsystemsref(4,:,2) = {[]; ...         % 4Ba (empty)
					  [0,1,-1,0]; ... % 4Pr
					  [-1,0,1,1]; ... % 4PyI
					  [-1,1,0,1]; ... % 4PyII
					  [-2,1,1,2]; ... % 4sP
					  };
% [b5](p)
slipsystemsref(5,1:5,1) = {1/3 * [2,-1,-1,-3]};        % b5
slipsystemsref(5,:,2) = {[]; ... 		   % 5Ba (empty)
					  [0,1,-1,0]; ...  % 5Pr
					  [1,0,-1,1]; ...  % 5PyI
					  [1,-1,0,1]; ...  % 5PyII
					  [2,-1,-1,2]; ... % 5sP
					  };
% [b6](p)
slipsystemsref(6,1:5,1) = {1/3 * [-1,-1,2,3]};         % b6
slipsystemsref(6,:,2) = {[]; ...         % 6Ba (empty)
					  [1,-1,0,0]; ... % 6Pr
					  [1,0,-1,1]; ... % 6PyI
					  [0,1,-1,1]; ... % 6PyII
					  [1,1,-2,2]; ... % 6sP
					  };
% [b7](p)
slipsystemsref(7,1:5,1) = {1/3 * [-1,-1,2,-3]};        % b7
slipsystemsref(7,:,2) = {[]; ...         % 7Ba (empty)
					  [1,-1,0,0]; ... % 7Pr
					  [1,0,-1,1]; ... % 7PyI
					  [0,1,-1,1]; ... % 7PyII
					  [1,1,-2,2]; ... % 7sP
					  };
% [b8](p)
slipsystemsref(8,1:5,1) = {1/3 * [-1,2,-1,3]};         % b8
slipsystemsref(8,:,2) = {[]; ...         % 8Ba (empty)
					  [-1,0,1,0]; ... % 8Pr
					  [1,-1,0,1]; ... % 8PyI
					  [0,-1,1,1]; ... % 8PyII
					  [1,-2,1,2]; ... % 8sP
					  };
% [b9](p)
slipsystemsref(9,1:5,1) = {1/3 * [-1,2,-1,-3]};        % b9
slipsystemsref(9,:,2) = {[]; ... 		   % 9Ba (empty)
					  [-1,0,1,0]; ...  % 9Pr
					  [0,1,-1,1]; ...  % 9PyI
					  [-1,1,0,1]; ...  % 9PyII
					  [-1,2,-1,2]; ... % 9sP
					  };

% Type <c> Burgers vectors
% [b10](p)
slipsystemsref(10,1:5,1) = {1/3 * [0,0,0,3]};          % b10
slipsystemsref(10,:,2) = {[]; ...         % 10Ba (empty)
					   [0,1,-1,0; ...  % 10PrI
                       -1,0,1,0; ...   % 10PrII
                        1,-1,0,0]; ... % 10PrIII
					   []; ...         % 10PyI (empty)
					   []; ...         % 10PyII (empty)
					   []; ...         % 10sP (empty)
					   };
                   
% Maximum number of possible data directions:
maxdatadirs = CalculateMaxDataDirs(slipsystemsref);
maxedgedatadirs = maxdatadirs - 1;


%%% Transformation of crystal grain slip systems

% Convert slip systems to cartesian:
slipsystemscart = SlipSystemsToCartesian(slipsystemsref, refHCPa1,refHCPa2,refHCPa3,refHCPa4);

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

% Macro CRSS values (see [3])
crsscoeffsstruct = struct;
% <a>:
crsscoeffsstruct.a_Ba = 204e6 / pressureSI; % 204 MPa
crsscoeffsstruct.a_Pr = 153e6 / pressureSI; % 153 MPa
crsscoeffsstruct.a_PyI = crsscoeffsstruct.a_Pr;
crsscoeffsstruct.a_PyII = crsscoeffsstruct.a_Pr;
% <c+a>:
crsscoeffsstruct.ca_PyI = 532e6 / pressureSI; % 532 MPa
crsscoeffsstruct.ca_PyII = crsscoeffsstruct.ca_PyI;
crsscoeffsstruct.ca_sP = crsscoeffsstruct.ca_PyI;
crsscoeffsstruct.ca_Pr = crsscoeffsstruct.ca_PyI;
% <c>:
crsscoeffsstruct.c_Pr = crsscoeffsstruct.a_Pr;
% Sessile:
crsscoeffsstruct.sessile = 10e9 / pressureSI; % 10 GPa

% Organize CRSS coefficients in matrix:
crsscoefficients = OrganizeCoefficients( ... % size (numburgs,numplanetypes + 1)
numburgsA,numburgsCA,numburgsC,numburgs,numplanetypes, ... % Organization-related inputs
crsscoeffsstruct);                                         % Coefficients structure

if enableCRSS ~= 1
    enableCRSS = 0;
    crsscoefficients(:) = 0;
end


%%% Store mobility information inside structure

mobilitystruct = struct;

mobilitystruct.slipsystemsref = slipsystemscartrot; % Cell: size (numburgs,numplanetypes,2)
mobilitystruct.burgscart = burgscart; % Matrix: size (numburgs,3)
mobilitystruct.HCPa1 = HCPa1;
mobilitystruct.HCPa2 = HCPa2;
mobilitystruct.HCPa3 = HCPa3;
mobilitystruct.HCPa4 = HCPa4;
mobilitystruct.dragline = dragline;
mobilitystruct.dragclimb = dragclimb;
mobilitystruct.glidecoefficients = glidecoefficients; % Matrix: size (numburgs,numplanetypes,2)
mobilitystruct.dragsessile = glidecoeffsstruct.dragsessile;
mobilitystruct.crsscoefficients = crsscoefficients; % Matrix: size (numburgs,numplanetypes+1)
mobilitystruct.maxedgedatadirs = maxedgedatadirs;
mobilitystruct.MATLABrelease = MATLABrelease;





%%%%% Initial nodes and links %%%%%

%%% Nodes %%%

%%%
% Frank-Read source:
rnfrref =   [0  -0.5   0 7;
             0   0     0 0;
             0   0.5   0 7];
% Neutron irradiation-induced vacancy loop (square):
rnvacref =  [0   0.5  0.5 0; % 8 nodes
             0   0.5  0   0;
             0   0.5 -0.5 0;
             0   0   -0.5 0;
             0  -0.5 -0.5 0;
             0  -0.5  0   0;
             0  -0.5  0.5 0;
             0   0    0.5 0];
%%%

% Frank-Read source: 2Pr
b_frref = slipsystemscartrot{2,2,1};
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
centpos_lo = [dx/16 dy/2 dz/4];
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

%%% Higher

rnv1_hi = rnvac_hi;
rnv1_hi(:,1:3) = rnvac_hi(:,1:3) + (dz*sqrt(2)*1/8) * (b_vunit);

rnv2_hi = rnvac_hi;
scalecentrnv2_hi = (dz/2 - centpos_hi(3))/b_vunit(3);
rnv2_hi(:,1:3) = rnvac_hi(:,1:3) + scalecentrnv2_hi * (b_vunit);
%rnv3(:,1:3) = rnvac(:,1:3) + 25 * (110e-9 / lengthSI) * (-b_vunit);

rnv3_hi = rnvac_hi;
rnv3_hi(:,1:3) = rnvac_hi(:,1:3) + (dz*sqrt(2)*6/8) * (b_vunit);

%%% Lower

% rnv1 = rnvac;
% rnv1(:,1:3) = rnvac(:,1:3) + 40 * (10e-9 / lengthSI) * (-b_vunit);

rnv1_lo = rnvac_lo;
rnv1_lo(:,1:3) = rnvac_lo(:,1:3) + (dz*sqrt(2)*1/8) * (-b_vunit);

rnv2_lo = rnvac_lo;
scalecentrnv2_lo = (dz*1/2 - centpos_lo(3))/b_vunit(3);
rnv2_lo(:,1:3) = rnvac_lo(:,1:3) + scalecentrnv2_lo * (b_vunit);
% centrnv2_lo = sum(rnv2_lo(:,1:3),1) / size(rnv2_lo,1);
% rnv2delta = [(vaclen*4) (vaclen*4) (vaclen*4)];
% viewaxisframes = [(centrnv2_lo(1) - rnv2delta(1)/2) (centrnv2_lo(1) + rnv2delta(1)/2) ...
%                   (centrnv2_lo(2) - rnv2delta(2)/2) (centrnv2_lo(2) + rnv2delta(2)/2) ...
%                   (centrnv2_lo(3) - rnv2delta(3)/2) (centrnv2_lo(3) + rnv2delta(3)/2)];
%rnv3(:,1:3) = rnvac(:,1:3) + 25 * (110e-9 / lengthSI) * (-b_vunit);
%viewaxisframes = [dx*1/10 dx*5/32 dy*3/8 dy*5/8 dz*1/4 dz*9/16];
%viewaxisframes = [dx*0.98/16 dx*6/32 dy*3/8 dy*5/8 dz*0.98/16 dz*9/16];
viewaxisframes = [dx*0.98/16 dx*3.7/32 dy*3/8 dy*5/8 dz*0.98/4 dz*9/16];

rnv3_lo = rnvac_lo;
rnv3_lo(:,1:3) = rnvac_lo(:,1:3) + (dz*sqrt(2)*6/8) * (-b_vunit);

% rnv5 = rnvac;
% rnv5(:,1:3) = rnvac(:,1:3) + 25 * (210e-9 / lengthSI) * (-b_vunit);

%%% Collect nodes and plot

% Higher:
%rn = [      rnfr_hi;    rnv1_hi; rnv2_hi; rnv3_hi];
% Lower:
%rn = [rn;   rnfr_lo;    rnv1_lo; rnv2_lo; rnv3_lo];
rn = [rnfr_lo; rnv2_lo;];

% Plot nodes:
rnfigname = append('rn_input_',INPUTNAME);
rnfig = figure('Name',rnfigname,'visible','off');
grid on, hold on
for n = 1:size(rn,1)
    plot3(rn(n,1),rn(n,2),rn(n,3),'.','MarkerSize',5,'Color','b')
end
ax = gca;
ax.XLim = [0 dx];
ax.YLim = [0 dy];
ax.ZLim = [0 dz];
set(rnfig, 'visible', 'on');
savefilefig = append('Figures/rn_input_',INPUTNAME,'.fig');
savefig(rnfig,savefilefig,'compact');
close(rnfig);
% To open saved figure: openfig('Figures/rn_input_xxx.fig')


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

% linksfr_hi = linksfr;
% 
% linksv1_hi = linksvac;
% v1hi_idx = max(linksfr_hi(:,1:2),[],'all');
% linksv1_hi(:,1:2) = linksvac(:,1:2) + v1hi_idx;
% 
% linksv2_hi = linksvac;
% v2hi_idx = max(linksv1_hi(:,1:2),[],'all');
% linksv2_hi(:,1:2) = linksvac(:,1:2) + v2hi_idx;
% 
% linksv3_hi = linksvac;
% v3hi_idx = max(linksv2_hi(:,1:2),[],'all');
% linksv3_hi(:,1:2) = linksvac(:,1:2) + v3hi_idx;

% linksv4_hi = linksvac;
% linksv4_hi(:,1:2) = linksvac(:,1:2) + 27;

% linksv5_hi = linksvac;
% linksv5_hi(:,1:2) = linksvac(:,1:2) + 35;

%%% Lower

linksfr_lo = linksfr;
frlo_idx = 0;
linksfr_lo(:,1:2) = linksfr(:,1:2) + frlo_idx;

linksv1_lo = linksvac;
v1lo_idx = max(linksfr_lo(:,1:2),[],'all');
linksv1_lo(:,1:2) = linksvac(:,1:2) + v1lo_idx;

% linksv2_lo = linksvac;
% v2lo_idx = max(linksv1_lo(:,1:2),[],'all');
% linksv2_lo(:,1:2) = linksvac(:,1:2) + v2lo_idx;
% 
% linksv3_lo = linksvac;
% v3lo_idx = max(linksv2_lo(:,1:2),[],'all');
% linksv3_lo(:,1:2) = linksvac(:,1:2) + v3lo_idx;

% linksv4_lo = linksvac;
% linksv4_lo(:,1:2) = linksvac(:,1:2) + 27;

% linksv5_lo = linksvac;
% linksv5_lo(:,1:2) = linksvac(:,1:2) + 35;

%%% Collect links

% Higher:
%links = [           linksfr_hi;     linksv1_hi; linksv2_hi; linksv3_hi];
% Lower:
%links = [links;     linksfr_lo;     linksv1_lo; linksv2_lo; linksv3_lo];
links = [linksfr_lo; linksv1_lo];





%% PRINT INFORMATION %%


%%% Material
fprintf('================ %s ================\n\n', mfilename);
fprintf('     Material: %s\n\n',        materialname);
fprintf('  Temperature: %0.2e K\n',     temperature);
fprintf('Shear modulus: %0.2e GPa\n',   pressureSI/1e9);
fprintf('Poisson ratio: %0.2e\n',       NU);
fprintf('c/a ratio: %0.2e\n',           HCPc);

%%% Selected subprocess functions
fprintf('\n--------------   Subprocesses   --------------\n\n');
fprintf('Mobility function: %s\n',      mobility);
fprintf('Integrator function: %s\n',    integrator);

%%% Simulation units in SI
fprintf('\n-------------- Simulation Units --------------\n\n');
fprintf('Distance: %0.2e m\n',          lengthSI);
fprintf('Time: %0.2e s\n',              timeSI);
fprintf('Velocity: %0.2e m/s\n',        velocitySI);
fprintf('Pressure: %0.2e Pa\n',         pressureSI);
fprintf('Force: %0.2e Pa m^2\n',        forceSI);
fprintf('Nodal force: %0.2e Pa m\n',    nodalforceSI);
fprintf('Drag: %0.2e Pa s\n',           dragSI);
fprintf('Temperature: %0.2e K\n',       temperatureSI);

%%% Warnings
% Nodal network:
if any(rn(:,1) > dx) || any(rn(:,2) > dy) || any(rn(:,3) > dz)
    warning('Material volume does not fully contain the nodal network.')
end
% Drag system:
if min(glidecoefficients,[],'all') / max(glidecoefficients,[],'all') < eps(1e1)
    warning('Very large disparity of drag coefficients.')
end
% Topological constraints:
if ~(2*lmin < lmax)
    warning('Recommended constraints on lmax and lmin not imposed.')
end
if ~(0 <= areamin && 4*areamin < areamax && areamax <= lmax^2*sqrt(3)/4)
    warning('Recommended topological constraints on areamax and areamin not imposed.')
end





%% REFERENCES %%


% [1]: Tarleton 2015: http://dx.doi.org/10.1016/j.actamat.2015.01.030
% [2]: Drouet 2014: https://doi.org/10.1016/j.jnucmat.2013.11.049
% [3]: Gong 2015: https://doi.org/10.1016/j.actamat.2015.06.020
% [4]: Introduction to Dislocations. Hull and Bacon, (2011)
% [5]: Khater 2010: https://doi.org/10.1016/j.actamat.2010.01.028
