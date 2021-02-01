%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that includes NUM_SOURCES number of sources in the form of pinned
% b=<111> sources of *constant* size, determined by parameter DIST_SOURCE
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% SOURCE GENERATION PARAMETERS
a_mag=3.14E-4; 
mumag = 160E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 30;
DIST_SOURCE = 0.5/a_mag;

%% FEM PARAMETERS
%Cantilever

dx=30/a_mag; %10micron
dy=5/a_mag; %2micron
dz=5/a_mag; %2micron

mx=30; % number of elements along beam length
loading=1; 
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];

%% MATERIAL CONSTANTS

%MU = 160E9; %160GPa
MU = 1;
NU = 0.305; 

%% DDLab PARAMETERS

%Dislocation nodes and segments generator
[rn,links] = FRSGenerator(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,0.5*dx,dy,dz);

%Edge and screw glide and climb mobility parameters
mobility='mobbcc0'; 
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge=1;
Bscrew=1; 
Beclimb=1e10;
Bline=1e-4*min(Bscrew,Bedge);

%Meshing
maxconnections=4; 
lmax =0.5/a_mag; %1 microns 
lmin = lmax/5; % 0.1 microns
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
%dt0=10E-9; %1 ns
dt0=1E7;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot=1E9;
doplot=1; % frame recording: 1 == on, 0 == off
totalSimTime = 1E12;
curstep = 0;
simTime = 0;

%Integrator
integrator='int_trapezoid'; 
%integrator='int_trapezoid_stoc'; %in development
a=lmin/sqrt(3)*0.5; 
Ec = MU/(4*pi)*log(a/0.1); 
rann = 0.5*a; 
rntol = 0.5*rann; % need to do convergence studies on all these parameters
rmax = lmax;

%Plotting
plotfreq=100; 
plim=12/a_mag; %12microns
viewangle=[-35,15]; 
printfreq=1; 
printnode=2; 

