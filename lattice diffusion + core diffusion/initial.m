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
clc; clear;

global tem bur 


%% MATERIAL CONSTANTS
MU = 0.83e11; % Pa
NU = 0.29; 
bur=0.2473; %nm
tem = 1300; %K

CRYSTAL_STRUCTURE = 'bcc';
%% FEM PARAMETERS
dx=500/bur; %
dy=dx;
dz=dx/2;
SimBox=[dx dy dz];

mx=10; % number of elements along beam length
loading=0; % 1 for compression , 0 for bending
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];

        %%
%Edge and screw glide and climb mobility parameters
global flag_c flag_b doclimb doglide
mobility='mobbcc_full'; 
doclimb=1; %flag for climb motion
flag_b=1; % flag for bulk diffusion controlled climb
flag_c=1; % flag for core diffusion controlled climb
doglide=0; %flag for glide motion


% ------------parameters for climb, diffusivity
if doclimb
    
    R = 8.314; % J/mol/K
    Boltz    = 1.3807e-23; % J/K
    Atomvol  = (2.856)^3*1e-30;  %volum of the atom, m^3/s
    rc       = 5*bur*1e-9;  % dislocation core radius, m
    
    if flag_b  % diffusivity for bulk diffusion  (https://www.tf.uni-kiel.de/matwis/amat/iss/kap_5/illustr/s5_2_3d.html)
        pre_exponential = 1e-4; % m^2/s
        active_energy = 174*1e3; % kJ/mol  1.8034eV
        D_bulk = pre_exponential*exp(-active_energy/0.6/R/tem);  %  bulk diffusivity, m^2/s
        Dv=D_bulk*Atomvol/Boltz/tem;  % effective diffusivity for bulk diffusion (D_b V/kTb) m^3/J.s
    end
    
    if flag_c  
        
        pre_exponential_c = 1e-23; % ac*Dc0  m^4/s
        active_energy = 174*1e3; % kJ/mol  1.8034eV
        D_core = pre_exponential_c*exp(-active_energy/R/tem);  % core diffusivity
        Dc=D_core*Atomvol/Boltz/tem;
    end
end

%---------------parameters for glide
global Bscrew Bedge Beclimb Bline
Bedge=1;
Bscrew=2; 
Beclimb=1e10;
Bline=1e-4*min(Bscrew,Bedge);

%Meshing
maxconnections=4; 
lmax = 50;
lmin = 20;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=100*areamin; 

doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on

Boundary=[0 0 0];

%Simulation time
dt0=1e-4;
totalsteps =1e6;
curstep = 0;

%Integrator
integrator='int_trapezoid_bb_old'; 
a=5; 
Ec = MU/(4*pi)*log(a/0.1); 
rann = 2*a; 
rntol = 2; % need to do convergence studies on all these parameters
rmax = lmax;

%Plotting
plotfreq=100; 
savefreq=100;

do_regeneration=1; 

%GPU Setup
n_threads=512;
 
%% DDlab parameters   
if do_regeneration == 1
    
    % ----------- generate a new structure
    [rn,links]=inicon_loop();
        
else
    %-----------read from a given initial structure
    fid =   fopen('.\input\initial_structure','r');
    l1=fread(fid,1,'integer*4');
    rn   =   fread(fid,[l1 4],'real*4');
    l2=fread(fid,1,'integer*4');
    links   =   fread(fid,[l2 8],'real*4');
    fclose(fid);

end
    
%parameters for the switch between climb and glide
flagc=0;
num_glide=0;  % number of steps which are steady during glide motion
num_climb=0;  % number of steps in which the climb distance is larger than the critical value
numofstep=20; % number of step for the swich between glide and climb
crit_glide=2e-2; % switch to climb when strain rate < crit_glide
crit_climb=1;   % switch to glide when max(vn_c*dt) > crit_climb
strat_change=0; % flag for the begining of the change

num_fig=0; % number of saved figures


appliedstress = -2e8 * [ 0  0  0
                         0  0  0
                         0  0  0];   % sigma                     
deltStress =0;




%% draw the inital figure
hold on
plim =SimBox(1);
viewangle = [135 45];
plotnodes(rn,links,plim); 
view(viewangle); 
hold on
PlotBox_b(SimBox);
   
hold on
plim =5000;
viewangle = [135 45];
plotnodes(rn,links,2*plim); 
view(viewangle);

hold on