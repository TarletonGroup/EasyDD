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

%HY20190705: simplest cases of "nanoindentation"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% SOURCE GENERATION PARAMETERS
amag = 2.856e-4; 
LENGTH_b = sqrt(3)/2;
mumag = 82E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
NUM_SOURCES =1;
DIST_SOURCE = 0.2/amag;
LENGTH_SOURCE = 20*LENGTH_b;

%% FEM PARAMETERS
%Cantilever

dx=1/amag*1; %10micron=10/amag
dy=1/amag*1; %
dz=1/amag*0.2*1; %

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
% mobility='mobbcc_bb'; 
% mobility='Hmobbcc_peierls';
mobility='Hmobbcc7rotation';
global Bscrew Bedge Beclimb Bline Bclimb
%Bedge=1e-4; %Pa s
%Bscrew=1e-5; %Pa s
%Beclimb=1e5; %Pa s - really big
%Bline=1e-4*min(Bscrew,Bedge);
Bedge=1E-4*1E0;
Bscrew=1e-4*1E0; 
Beclimb=1e4*min(Bscrew,Bedge);
Bclimb = Beclimb;
Bline=1e-4*min(Bscrew,Bedge);

global unloadflag unloadcount

global Ubarglobal

unloadflag = 0;
unloadcount = 0;
Ubarglobal = 0;

%% MATERIAL CONSTANTS

MU = 1;
NU = 0.3;

eps = 1E-6;

global burgsref planesref burgnumbers planenumbers rotationBCC;

% planesref = [0 1 1; 1 -1 0; 1 1 0]/sqrt(2);
% burgsref = [-1 -1 1; 1 1 -1; 1 -1 1]/2;

% planesref = [1 0 -1; 1 -1 0; 1 1 0]/sqrt(2);
planesref = [1 -1 0; 1 -1 0; 1 1 0]/sqrt(2);
burgsref = [1 1 1; 1 1 -1; 1 -1 1]/2;

burgnumbers = [1 2 3];
planenumbers = [1 2 3];

%HY20181010: calculate the rotation matrix.
% crossV = cross(planesref(1,:),burgsref(1,:));
% crossV = [-1 -1 -1]/sqrt(3);
% rotV = vrrotvec(crossV,[0 0 -1]);
% W = [0 -rotV(3) rotV(2);rotV(3) 0 -rotV(1);-rotV(2) rotV(1) 0];
% Q1 = [eye(3)+sin(rotV(4))*W+(2*(sin(rotV(4)/2))^2)*W*W]'

%HY20191001: calculate the rotation matrix.

n1=[1 0 0];n2=[0 1 0];n3=[0 0 1];
e1=[0 0 1];e2=[1 -1 0];e3=[1 1 0];%(110)
% e1=[1 1 -2];e2=[1 -1 0];e3=[1 1 1];%(111)
e1 = e1/norm(e1);e2 = e2/norm(e2);e3 = e3/norm(e3);%e1, e2, e3 must be normalized firstly
Q1=[e1;e2;e3];  % your general formulas R

% %% http://www.meshola.com/Articles/converting-between-coordinate-systems
%        Vector3d X1 = XAxisWorld;            %  This is vector (1,0,0)
%        Vector3d X2 = YAxisWorld;            %  This is vector (0,1,0)
%        Vector3d X3 = ZAxisWorld;            %  This is vector (0,0,1)
%        
%        %  These vectors are the local X,Y,Z of the rotated object
%        Vector3d X1Prime = XAxisLocal;
%        Vector3d X2Prime = YAxisLocal;
%        Vector3d X3Prime = ZAxisLocal;
%        %  This matrix will transform points from the rotated axis to the world
%        LocalToWorldTransform = new Matrix3x3()
%        {
%             M11 = (float)Vector3d.DotProduct(X1, X1Prime),
%             M12 = (float)Vector3d.DotProduct(X1, X2Prime),
%             M13 = (float)Vector3d.DotProduct(X1, X3Prime),
%             M21 = (float)Vector3d.DotProduct(X2, X1Prime),
%             M22 = (float)Vector3d.DotProduct(X2, X2Prime),
%             M23 = (float)Vector3d.DotProduct(X2, X3Prime),
%             M31 = (float)Vector3d.DotProduct(X3, X1Prime),
%             M32 = (float)Vector3d.DotProduct(X3, X2Prime),
%             M33 = (float)Vector3d.DotProduct(X3, X3Prime),
%        };
%        %  This matrix will transform points from the world back to the rotated axis
%        WorldToLocalTransform = new Matrix3x3()
%        {
%             M11 = (float)Vector3d.DotProduct(X1Prime, X1),
%             M12 = (float)Vector3d.DotProduct(X1Prime, X2),
%             M13 = (float)Vector3d.DotProduct(X1Prime, X3),
%             M21 = (float)Vector3d.DotProduct(X2Prime, X1),
%             M22 = (float)Vector3d.DotProduct(X2Prime, X2),
%             M23 = (float)Vector3d.DotProduct(X2Prime, X3),
%             M31 = (float)Vector3d.DotProduct(X3Prime, X1),
%             M32 = (float)Vector3d.DotProduct(X3Prime, X2),
%             M33 = (float)Vector3d.DotProduct(X3Prime, X3),
%        };
   
   
rotationBCC = Q1';

burgsref = burgsref*rotationBCC;
planesref = planesref*rotationBCC;

%HY20180205
planesref(find(abs(planesref)<eps)) = 0;
burgsref(find(abs(burgsref)<eps)) = 0;

%Plotting
plotfreq=100; 
savefreq=50; 
plim=12/amag; %12microns
viewangle2=[-81,81]; 
viewangle=[-25,40]; 
printfreq=plotfreq; 
printnode=2; 

%Dislocation nodes and segments generator
[rn,links,normals,b_vecs,midpoints] = checkGeneratorbccHYrotate(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz,plim,LENGTH_SOURCE,burgsref,planesref,burgnumbers,planenumbers);

%Meshing
maxconnections=4; 
% lmax =0.1/amag*0.1.5;
% lmin = 0.02/amag*0.1.5;
% lmax =LENGTH_SOURCE/2;
% lmax =dy/50;
lmax = LENGTH_SOURCE/2;
lmin = lmax/5;
global areamin;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
dosave = 0;
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on

%Simulation time
dt0=1E-1;
dtmax=1E3;

intSimTime = 0;
sinTime = 0;
%dtplot=2E-9; %2ns
dtplot=0;
doplot=1; % frame recording: 1 == on, 0 == off
totalSimTime = 1E15;
curstep = 0;
simTime = 0;

%Integrator
% integrator='int_trapezoid_bbab'; 
% integrator='int_trapezoid_bb'; 
integrator='int_trapezoid'; 
%integrator='int_trapezoid_stoc'; %in development
% a=lmin/sqrt(3)*0.5; 
% a = 10*LENGTH_b;
a=lmin/sqrt(3)*0.5; 
Ec = MU/(4*pi)*log(a/0.1); 
rann = 0.5*a; 
% % rann = 0.001/amag;%HY20190509: accordin to Ed, we don't care about the interaction outside a radius of 100nm
% rann = lmin/2;%HY20190509: accordin to Ed, we don't care about the interaction outside a radius of 100nm
rntol = 0.5*rann*100; % need to do convergence studies on all these parameters
rmax = rann/2;
maxSeg = lmax;

% %HY20180316: hydrgoen related parameters
% delta_VH =1;
% k_BH = 1;
% TH = 400;
% c_max = 2E-7;
% % c_max = 0E-7;

%HY20190404:
global delta_VH k_BH TH c_max potential_remote

%HY20180320: hydrogen related parameters
%HY20180409: corrected hydrogen related parameters
% delta_VH =0.0907;
delta_VH = 1.4E-12/amag/amag/amag;
k_BH = 1.38E-23*1E6/amag*0.1E6/(amag*amag*mumag);
TH = 300;
c_max = 8.46E28*1E-18/(1/amag)^3;%HY20180409: N_L=8.46E28/m^3 in bcc steel
kappa_remote = 0E-6;%HY20180408: corresponding to 1appm in bcc steel
% kappa_remote = 1E-7;%HY20180409: corresponding to 0.1appm in bcc steel
% kappa_remote = 0;
potential_remote = log(1/kappa_remote-1)*k_BH*TH;

% tpause = 1E-6*1E-3/(Bedge/mumag);

rn00 = rn;
links00 = links;


tpause = 1.8E6;

global load_rate Uchange

load_rate = dz*1e-6;
Uchange = 0.0;


global unique_errmag_point unique_distmag_point

unique_errmag_point = 0;
unique_distmag_point = 0;

activation = 0;
PKactiv = 0.02;

global do_cross_slip do_reset last_step L_cross_crit
do_cross_slip=1
do_reset = 1;
last_step = 0;
checkcrossfreq = 1;
addsourcefreq = 5;
L_cross_crit = 4*lmax;

combinations110 = 1/sqrt(2) * [1 1 0 ; -1 -1 0 ; 1 -1 0 ; -1 1 0 ; ...
                               1 0 1 ; -1 0 -1 ; 1 0 -1 ; -1 0 1 ; ...
                               0 1 1 ; 0 -1 -1 ; 0 1 -1 ; 0 -1 1 ];
keepgoing = 0;

global vdotl

output = 'C:\Users\Haiyang\Desktop\DDD\HY\indentationHHDDLabFEMab20191013\temp\temp';
save(fullfile(output));

global vdotl

run dd3drotatebcc