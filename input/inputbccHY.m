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
%HY20180414: modified significantly by HY to consider bcc iron for HHDDD
%simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% SOURCE GENERATION PARAMETERS
amag = 2.856e-4; 
LENGTH_b = sqrt(3)/2;
mumag = 82E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 10;
DIST_SOURCE = 0.2/amag;
LENGTH_SOURCE = 500*LENGTH_b;

%% FEM PARAMETERS
%Cantilever

dx=12/amag; %10micron=10/amag
dy=3/amag; %
dz=3/amag; %

% dx=30/amag; %10micron=10/amag
% dy=5/amag; %
% dz=5/amag; %

mx=50; % number of elements along beam length
loading=1; 
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz]; 

%% DDLab PARAMETERS

%Edge and screw glide and climb mobility parameters
mobility='mobbcc1'; 
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge=5E-4;
Bscrew=100e-4; 
Beclimb=1e5;
Bline=1e-4*min(Bscrew,Bedge);

global cOVERa
global MU NU maxSeg
global burgsref planesref edgesref bpiref ppbref planestyperef burgnumbers planenumbers

global node_NRF flag_NRF flag_Junction node_Junction1 node_Junction2

global unloadflag unloadcount

global Ubarglobal

flag_NRF = 0;
flag_Junction = 0;
unloadflag = 0;
unloadcount = 0;
Ubarglobal = 0;

%% MATERIAL CONSTANTS

MU = 1;
NU = 0.29;

%Plotting
plotfreq=500; 
savefreq=10; 
plim=12/amag; %12microns
% viewangle=[0,0]; 
viewangle=[30,60]; 
printfreq=1000; 
printnode=2; 

%Dislocation nodes and segments generator
[rn,links] = checkGeneratorbccHY(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz,plim,LENGTH_SOURCE);


%create mirror of prismatic loops (outside boundary)
% rn_mirror = [rn(:,1)+dx , rn(:,2) , rn(:,3)+dx , zeros(length(rn(:,1)),1)+67];
% links_mirror = links;
% links_mirror(:,1) = links(:,2) + max(max(links(:,1:2)));
% links_mirror(:,2) = links(:,1) + max(max(links(:,1:2)));
% 
% rn = vertcat(rn,rn_mirror);
% links = vertcat(links,links_mirror);
%HY20180112: ^^^^initial dislocation structures^^^^

%Meshing
maxconnections=4; 
lmax =0.1/amag;
lmin = 0.02/amag;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
dosave = 0;
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
dt0=1E8;
dt=1E8;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot=0;
doplot=1; % frame recording: 1 == on, 0 == off
totalSimTime = 1E15;
curstep = 0;
simTime = 0;

%Integrator
integrator='int_trapezoid'; 
%integrator='int_trapezoid_stoc'; %in development
% a=lmin/sqrt(3)*0.5; 
a = 40*LENGTH_b;
Ec = MU/(4*pi)*log(a/0.1); 
rann = 0.5*a; 
rntol = 0.5*rann; % need to do convergence studies on all these parameters
rmax = lmax;
maxSeg = lmax;

% %HY20180316: hydrgoen related parameters
% delta_VH =1;
% k_BH = 1;
% TH = 400;
% c_max = 2E-7;
% % c_max = 0E-7;

%HY20180320: hydrogen related parameters
%HY20180409: corrected hydrogen related parameters
% delta_VH =0.0907;
delta_VH = 1.4E-12/amag/amag/amag;
k_BH = 1.38E-23*1E6/amag*1E6/(amag*amag*mumag);
TH = 300;
c_max = 8.46E28*1E-18/(1/amag)^3;%HY20180409: N_L=8.46E28/m^3 in bcc steel
% kappa_remote = 1E-6;%HY20180408: corresponding to 1appm in bcc steel
kappa_remote = 1E-7;%HY20180409: corresponding to 0.1appm in bcc steel
% kappa_remote = 0;
potential_remote = log(1/kappa_remote-1)*k_BH*TH;

% tpause = 1E-6*1E-3/(Bedge/mumag);

rn00 = rn;
links00 = links;


tpause = 1.8E10;
use_gpu = 1;
n_threads = 256;
para_scheme = 1;

save('dcg_20_07_18')
