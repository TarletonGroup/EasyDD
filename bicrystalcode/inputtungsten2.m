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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% SOURCE GENERATION PARAMETERS
amag=3.18e-4;
mumag = 135E3; % MPa only used for plotting

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 4;
DIST_SOURCE = 0.1/amag;

%% FEM PARAMETERS
%Cantilever

dx=19.5/amag; %10micron
dy=5.4/amag; %2micron
dz=5.3/amag; %2micron

mx=40; % number of elements along beam length
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
NU = 0.28;

%% DDLab PARAMETERS

% c=1/amag;
% x=0.25/(amag*sqrt(6));
% y=0.25/(amag*sqrt(3));
%
% rn=[         3.7*c   2.7*c+2*x-y           0.25*c   0;
%            3.7*c+x       2.7*c-y         0.25*c+x   0;
%          3.7*c+2*x   2.7*c-2*x-y       0.25*c+2*x   0;
%        3.7*c+2*x+y     2.7*c-2*x     0.25*c+2*x+y   0;
%      3.7*c+2*x+2*y   2.7*c-2*x+y   0.25*c+2*x+2*y   0;
%        3.7*c+x+2*y       2.7*c+y     0.25*c+x+2*y   0;
%          3.7*c+2*y   2.7*c+2*x+y       0.25*c+2*y   0;
%            3.7*c+y     2.7*c+2*x         0.25*c+y   0];
%
% b1=[-1 -1 -1]/sqrt(3);
% n1=[-1 0 1]/sqrt(2);
% n2=[1 0 -1]/sqrt(2);
%
% links=[1 2   b1 n1;
%        2 3   b1 n1;
%        3 4   b1 n1;
%        4 5   b1 n1;
%        5 6   b1 n1;
%        6 7   b1 n1;
%        7 8   b1 n1;
%        8 1   b1 n1];

    [rn,links] = bccsurfsourcegen(NUM_SOURCES,DIST_SOURCE,dx,dy,dz);
%     [rn]=[dx-0.5*dx,dy-0.5*dy,dz-0.5*dz,7;
%           dx-0.5*dx,dy-0.5*dy,dz-0.55*dz,7;
%           dx-0.49*dx,dy-0.5*dy,dz-0.55*dz,7];
%     [links]=[1,2,1,0,0,1,0,0;
%              2,3,1,0,0,1,0,0;
%              3,1,1,0,0,1,0,0];

%create mirror of prismatic loops (outside boundary)
% rn_mirror = [rn(:,1)+dx , rn(:,2) , rn(:,3)+dx , zeros(length(rn(:,1)),1)+67];
% links_mirror = links;
% links_mirror(:,1) = links(:,2) + max(max(links(:,1:2)));
% links_mirror(:,2) = links(:,1) + max(max(links(:,1:2)));
%
% rn = vertcat(rn,rn_mirror);
% links = vertcat(links,links_mirror);

%%
%Edge and screw glide and climb mobility parameters
mobility='mobbcc_bb';
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge=20;
Bscrew=100;
Beclimb=1e12;
Bline=1e-4*min(Bscrew,Bedge);

%Meshing
maxconnections=4;
lmax =0.25/amag;
lmin = 0.05/amag;
areamin=lmin*lmin*sin(60/180*pi)*0.5;
areamax=20*areamin;
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
dt0=1E6;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot=5E4;
doplot=1; % frame recording: 1 == on, 0 == off
totalSimTime = (2/amag)/(100*1E3*dx*(1E-4/160E9));
curstep = 0;
simTime = 0;

%Integrator
integrator='int_trapezoid';
%integrator='int_trapezoid_stoc'; %in development
a=lmin/sqrt(3)*0.5;
Ec = MU/(4*pi)*log(a/0.1);
rann = 0.5*a; 
rntol = 50*rann; % need to do convergence studies on all these parameters
rmax = lmax;

%Plotting
plotfreq=10E7;
plim=12/amag; %12microns
viewangle=[-35,15];
printfreq=500;
printnode=2;

%GPU Setup

n_threads=512;
