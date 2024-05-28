%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that includes checking some simple cases
% in a compression of bcc iron.
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
% Distances: microns  %%%amag
% Time: s  
% Force: N     
% Pressure: Pa     %%% mu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear;
clear all
global tem amag RR QQ do_inclusion mumag SimBox strain0 rntol rmax do_twosource
 
%% SOURCE GENERATION PARAMETERS
amag= 2.87e-4; 
% bur=sqrt(3)/2*amag;
amag=sqrt(3)/2*amag; % bur as the length unit
mumag = 83E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES = 2;
%DIST_SOURCE = 0.18/amag;
DIST_SOURCE = 0.5/amag;


%% FEM PARAMETERS
%pillar

dx=4.0/amag; %1.2micron
dy=2.0/amag; %0.6micron
dz=4.0/amag; %0.6micron
SimBox=[dx dy dz];

mz=40; % number of elements along beam length
loading=1; % 1 for compression , 0 for bending
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
MU = 1; % 0.83e11Pa
NU = 0.29; 
tem = 900; %K

%%
%Edge and screw glide and climb mobility parameters
mobility='mobbcc_bb1b'; 
global Bscrew Bedge Beclimb Bline
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge=1;   %1e-2 Pa s
Bscrew=1; 
Beclimb=1e10;
Bline=1e-4*min(Bscrew,Bedge);

%Meshing
maxconnections=4; 
lmax =0.1/amag;
lmin = 0.04/amag;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
dt0=1E2;
dt0=1e6*2048*100;

intSimTime = 0;
simTime = 0;
%dtplot=2E-9; %2ns
% dtplot=5E4;
% doplot=1; % frame recording: 1 == on, 0 == off
%totalSimTime = (2/amag)/(100*1E3*dx*(1E-4/160E9));
totalSimTime = 10E-6;

curstep = 0;

%Integrator
integrator='int_trapezoid_bb_old'; 
%integrator='int_trapezoid_stoc'; %in development
% a=lmin/sqrt(3)*0.5;
a=5;
Ec = MU/(4*pi)*log(a/0.1); %%Daniel has found papers that explicitly give this quantity for Zr, might exist for other metals
rann = 4*a; 
rntol = 2*rann; % need to do convergence studies on all these parameters
rmax = 2*lmin;
% rmax=4*lmin;

%Plotting
plotfreq=20; 
savefreq=5;
% plim=12/amag; %12microns
viewangle=[-35,15]; 
% printfreq=500; 
% printnode=1000; 

%GPU Setup

n_threads=512;

% initiation: check the mechanisms
do_cross_slip=0;  % consider crossslip
L_cross_crit=30;  % critical length of the screw dislcoation to allow cross slip

do_inclusion =0;  % consider inclusion
do_twosource = 1;

volume_frac = 0.02; % volume frace of the inclusion 

do_climb =0;  % consder climb motion 

do_visualisestress = 0; % visualise the stress field
do_nucleation =0; % consider nucleation during the calculation

do_regeneration = 1 ;  % generate a new inital config

do_virtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

passivate_sur=0; %1 passivated surface; 0 free surface 
strain0= 0.02; % dialatational misfit strain;
 
%% DDlab parameters
% [rn,links] = bccsourcegen(NUM_SOURCES,DIST_SOURCE,dx,dy,dz);     
if do_regeneration 
    [rn,links] = bccsourcegen_bicrystal(NUM_SOURCES,DIST_SOURCE,dx,dy,dz);
    % ----------- generate a new structure
    %         [rn,links]=inicon_random_FRs(lmax,lmin,SimBox);
    
% rn=[-100 -100 dz/2 67;
%     dx+100 -100 dz/2 67;
%     dx+100 dy+100 dz/2 67;
%     -100 dy+100  dz/2 67;];
% links=[1 2 0 0 1 0 1 0;
%     2 3 0 0 1 0 1 0;
%     3 4 0 0 1 0 1 0;
%     4 1  0 0 1 0 1 0;];
        
else
    %-----------read from a given initial structure
    fid =   fopen('.\input\initial_structure','r');
    ll1=fread(fid,1,'integer*4');
    rn   =   fread(fid,[ll1 4],'real*8');
    ll2=fread(fid,1,'integer*4');
    links   =   fread(fid,[ll2 8],'real*8');
    fclose(fid);
    
end
    
%parameters for the change between climb and glide
flagc=0;
num_glide=0;  % number of steps which are steady during glide motion
num_climb=0;  % number of steps in which the climb distance is larger than the critical value
numofstep=20; % number of step for the swich between glide and climb
crit_glide=20; % switch to climb when strain rate < crit_glide
crit_climb=20;   % switch to glide when max(vn_c*dt) > crit_climb
strat_change=0; % flag for the begining of the change


num_fig=0; % number of saved figures

%% draw the inital figure
plim =SimBox;
PlotBox_b(SimBox);
plotnodes_climb1(rn,links,plim,0,0)


% save initial2

