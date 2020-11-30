%% Input file created by Daniel Hortelano Roig (17/07/2020) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sections are the following:
% 0) SIMULATION PATHS AND UNITS
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

% MATLAB version year:
MATLABversion = ver('MATLAB');
MATLABversionyear = str2double(MATLABversion.Release(3:6));





%% SIMULATION PATHS AND UNITS %%


%%% Add EasyDD subfolders to MATLAB path

% Main 'EasyDD' folder:
mainPath = '../../';
addpath(mainPath);

% Subprocesses folders:
addpath([mainPath 'Mobility']); % Mobility laws
addpath([mainPath 'Integrator']); % Time integrators

% Folder with functions designed for use in the input file:
addpath("InputFunctions");


%%% Simulation units in SI

lengthSI = 2.29e-10; % Length: [m]
pressureSI = 144E9; % Pressure: [Pa]
dragSI = 1.0e0; % Drag coefficient: [Pa s]
timeSI = dragSI/pressureSI; % Time: [s]
velocitySI = lengthSI/timeSI; % Velocity: [m/s]
forceSI = pressureSI * lengthSI^2; % Force: [Pa m^2]
nodalforceSI = pressureSI * lengthSI; % Nodal force: [Pa m]





%% MATERIAL CONSTANTS %%


% Material name:
materialname = "Be"; % Beryllium

% Lattice parameter magnitude:
length_inputunits = 1e-6; % Units: [m]
amag = lengthSI / length_inputunits; % Units: [microns]

% HCP a,c lattice parameters:
HCPa = 1.0; % Units: lengthSI
HCPc = 1.568; % Units: lengthSI

% Shear modulus:
pressure_inputunits = 1e6; % Units: [Pa]
mumag = pressureSI / pressure_inputunits; % Units: [MPa]
MU = 1; % Units: pressureSI

% Poisson's ratio:
NU = 0.064; % Units: unitless





%% SOURCE GENERATION PARAMETERS %%


DIST_SOURCE = 0.5/amag;





%% FEM PARAMETERS %%


% Cantilever Axes:
% ^z
% |     (y going into page)
% |
% |------>x

% Cantilever beam dimensions:
dx = 20/amag; % NOTE: 1/amag = 1 micron
dy = 7/amag;
dz = 7/amag;

% Number of finite elements along beam x-direction:
mx = 40;

% Vertices of beam:
vertices = [ 0,  0,  0;...
            dx,  0,  0;...
             0, dy,  0;...
            dx, dy,  0;...
             0,  0, dz;...
            dx,  0, dz;...
             0, dy, dz;...
            dx, dy, dz];

% Linear loading rate:
loading = 1; % = 1: on, = 0: off
Udot = (dx*lengthSI)*10e3*(timeSI/lengthSI);





%% DDLab PARAMETERS %%


%%%%% Meshing %%%%%

maxconnections = 8;
lmax = 0.25/amag;
lmin = 0.1/amag;
areamin = lmin*lmin*sin(60/180*pi)*0.5;
areamax = 20*areamin;
doremesh = 1; % flat set to 0 or 1 that turns the remesh functions off or on
docollision = 1; % flat set to 0 or 1 that turns collision detection off or on
doseparation = 1; % flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh = 1; % flat set to 0 or 1 that turns remeshing of virtual nodes off or on





%%%%% Integrator %%%%%


%%% Time integrator function

integrator = 'int_trapezoid';


%%% Time integration variables

a = lmin/sqrt(3)*0.5; % Isotropic spreading radius
Ec = MU/(4*pi)*log(a/0.1); % Segment core energy per unit length
rann = 0.5*a; % Segment annihilation distance
rntol = 0.5*rann; % Time integrator segment distance error tolerance
rmax = lmax; % Max distance a node may move during a time increment





%%%%% Simulation time %%%%%


%%% Nodal time

% Permitted absolute time increment:
dt0 = 1E-7/timeSI; % Maximum
dt = 1E-10/timeSI; % Minimum

totalSimTime = 6.25E-5/timeSI; % Absolute time duration of simulation

curstep = 0; % Initial time-step


%%% Plotting time

intSimTime = 0/timeSI; % Initial absolute time for plotting
dtplot = floor(totalSimTime/5); % Absolute time increment for plotting
doplot = 1; % Frame recording: 1 == on, 0 == off





%%%%% Plotting %%%%%

plotfreq = 1;
plim = 12/amag;
viewangle = [-35,15];
printfreq = 500;
printnode = 2;





%%%%% Mobility %%%%%


%%% Required inputs

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
% Units: Pa s
% Plane Type -- Burgers Vector Index -- Orientation

% Glide drag coefficients for <a> segments (Burgers vector indices 1-3):
Ba1Screw = 1.65e-3; Ba1Edge = 1.18e-3;
Pr1Screw = 7.15e-4; Pr1Edge = 5.7e-4;
PyI1Screw = 1.05e-3; PyI1Edge = 7.28e-4;
PyII1Screw = 1.05e-3; PyII1Edge = 7.28e-4;

% Glide drag coefficients for <c+a> segments (4-9):
Pr4Screw = 2.37e-2; Pr4Edge = 2.37e-2;
PyI4Screw = 2.37e-2; PyI4Edge = 2.37e-2;
PyII4Screw = 2.37e-2; PyII4Edge = 2.37e-2;
% WARNING! This sP data is not validated yet:
sP4Screw = 2.37e-2; sP4Edge = 2.37e-2;

% Glide drag coefficients for <c> segments (10):
Pr10Screw = 2.37e-2; Pr10Edge = 2.37e-2;

% Glide drag for junction segments:
dragsessile = 2.37e1;

% Climb drag:
dragclimb = dragsessile;

% Line drag:
dragline = 5.0e-5;

% Organize glide drag coefficients:
glidecoefficients = OrganizeGlideDragCoefficients( ...
numburgsA,numburgsCA,numburgsC,numburgs,numplanetypes, ... % Organization-related inputs
Ba1Screw,Ba1Edge,Pr1Screw,Pr1Edge,PyI1Screw,PyI1Edge,PyII1Screw,PyII1Edge, ... % <a> glide
Pr4Screw,Pr4Edge,PyI4Screw,PyI4Edge,PyII4Screw,PyII4Edge,sP4Screw,sP4Edge, ... % <c+a> glide
Pr10Screw,Pr10Edge, ...                                                        % <c> glide
dragsessile);                                                                  % Junction glide


%%% Glissile slip systems

slipsystemsref = cell(numburgs,numplanetypes,2);

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


%%% Slip systems transformation

% Convert slip systems to cartesian:
slipsystemscart = ...
    SlipSystemsToCartesian(slipsystemsref, refHCPa1,refHCPa2,refHCPa3,refHCPa4);

% Prepare to rotate slip systems:
%angx = 0/180 * pi;
angy = -45/180 * pi;
angz = 270/180 * pi;
%Rx = [1 0 0; 0 cos(angx) -sin(angx); 0 sin(angx) cos(angx)];
Ry = [cos(angy) 0 sin(angy); 0 1 0; -sin(angy) 0 cos(angy)];
Rz = [cos(angz) -sin(angz) 0; sin(angz) cos(angz) 0; 0 0 1];

rotmatrix = Ry * Rz; % Rotation matrix

% Rotate reference Burgers, normal, and HCP lattice vectors:
[slipsystemscartrot,HCPa1,HCPa2,HCPa3,HCPa4] = ...
    RotateCartesianSlipSystems(slipsystemscart, refHCPa1,refHCPa1,refHCPa1,refHCPa1, rotmatrix);

% Rotated reference Cartesian Burgers vectors:
burgscart = cat(1,slipsystemscartrot{:,1,1});


%%% Critical resolved shear stress (CRSS)

crss = 0;


%%% Store mobility information inside of a structure

mobilitystruct = struct; % Empty structure

mobilitystruct.slipsystems = slipsystemscartrot;
mobilitystruct.burgscart = burgscart;
mobilitystruct.HCPa1 = HCPa1;
mobilitystruct.HCPa2 = HCPa2;
mobilitystruct.HCPa3 = HCPa3;
mobilitystruct.HCPa4 = HCPa4;
mobilitystruct.dragline = dragline;
mobilitystruct.dragclimb = dragclimb;
mobilitystruct.glidecoefficients = glidecoefficients;
mobilitystruct.dragsessile = dragsessile;
mobilitystruct.crss = crss;
mobilitystruct.maxedgedatadirs = maxedgedatadirs;
mobilitystruct.MATLABversionyear = MATLABversionyear;





%%%%% Initial nodes and links %%%%%


%%% Nodes %%%

cn = 1/amag;
xn = 0.5/(amag*sqrt(6));
yn = 0.5/(amag*sqrt(3));

rn = [          6*cn   2.5*cn+2*xn-yn           0.5*cn   7; % Loop 1
               6*cn+xn       2.5*cn-yn         0.5*cn+xn   7;
             6*cn+2*xn   2.5*cn-2*xn-yn       0.5*cn+2*xn   0;
           6*cn+2*xn+yn     2.5*cn-2*xn     0.5*cn+2*xn+yn   0;
         6*cn+2*xn+2*yn   2.5*cn-2*xn+yn   0.5*cn+2*xn+2*yn   0;
           6*cn+xn+2*yn       2.5*cn+yn     0.5*cn+xn+2*yn   7;
             6*cn+2*yn   2.5*cn+2*xn+yn       0.5*cn+2*yn   7;
               6*cn+yn     2.5*cn+2*xn         0.5*cn+yn   0;
                 8*cn   2.5*cn+2*xn-yn           0.5*cn   7; % Loop 2
               8*cn-xn       2.5*cn-yn         0.5*cn+xn   0;
             8*cn-2*xn   2.5*cn-2*xn-yn       0.5*cn+2*xn   0;
           8*cn-2*xn-yn     2.5*cn-2*xn     0.5*cn+2*xn+yn   0;
         8*cn-2*xn-2*yn   2.5*cn-2*xn+yn   0.5*cn+2*xn+2*yn   0;
           8*cn-xn-2*yn       2.5*cn+yn     0.5*cn+xn+2*yn   0;
             8*cn-2*yn   2.5*cn+2*xn+yn       0.5*cn+2*yn   7;
               8*cn-yn     2.5*cn+2*xn         0.5*cn+yn   0];

rn(1:8,1:3) = rn(1:8,1:3) * 2;


%%% Links %%%

% (1): 2Ba
b_1def = slipsystemscartrot{2,1,1};
b_1mag = a/10;
b_1 = b_1mag * (b_1def ./ norm(b_1def));
n_1 = slipsystemscartrot{2,1,2};


% (2): 10PrIII
b_2def = slipsystemscartrot{10,2,1};
b_2mag = a/5;
b_2 = b_2mag * (b_2def ./ norm(b_2def));
n_2 = slipsystemscartrot{10,2,2}(3,:);

links = [   1  2   b_1 n_1;
            2  3   b_1 n_1;
            3  4   b_1 n_1;
            4  5   b_1 n_1;
            5  6   b_1 n_1;
            6  7   b_1 n_1;
            7  8   b_1 n_1;
            8  1   b_1 n_1;
            9 10   b_2 n_2;
           10 11   b_2 n_2;
           11 12   b_2 n_2;
           12 13   b_2 n_2;
           13 14   b_2 n_2;
           14 15   b_2 n_2;
           15 16   b_2 n_2;
           16  9   b_2 n_2];





%% PRINT INFORMATION %%


% Input file:
fprintf('============== %s ==============\n\n', mfilename);
fprintf('Material: %s\n', materialname);

% Selected subprocess functions:
fprintf('\n-------------- Subprocesses --------------\n\n');
fprintf('Mobility function: %s\n', mobility);
fprintf('Integrator function: %s\n', integrator);

% Simulation units in SI:
fprintf('\n-------------- Simulation Units --------------\n\n');
fprintf('   Distance: %0.2s m\n',       lengthSI);
fprintf('       Time: %0.2s s\n',       timeSI);
fprintf('   Velocity: %0.2s m/s\n',     velocitySI);
fprintf('   Pressure: %0.2s Pa\n',      pressureSI);
fprintf('      Force: %0.2s Pa m^2\n',  forceSI);
fprintf('Nodal force: %0.2s Pa m\n',    nodalforceSI);
fprintf('       Drag: %f Pa s\n',       dragSI);





%% REFERENCES %%


% [1]: Aubry 2016: http://dx.doi.org/10.1016/j.jmps.2016.04.019
