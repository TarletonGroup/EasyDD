clear all

global Bscrew Bedge Beclimb Bline

bmag = 2.74E-4; %microns
mumag = 160e3;  % MPa
l = 0.5/bmag;
b = 1/2*[1 1 1];


side1 = [1  1  1]/sqrt(3);


n = [-1 0 1] %cross(side2,b);

%start from center of grid

%P1 = [0.52 0.52  1.7]/bmag; %exiting top
%P1 = [0.52 0.1 0.52]/bmag; %lateral
dx=5; %microns
dy=2; %microns 
dz=2; %microns
ns= 10 %number of sources
 
P1 = [dx 0.5*dy, 0.8*dz]/bmag;
P2m = P1 + l*side1*0.5;
P2 = P1 + l*side1;


P3 = [0.05*dx 0.5*dy, 0.2*dz]/bmag;
P3m = P3 + l*side1*0.5;
P4 = P3 + l*side1;


rn = [ 
       P1, 7;
       P2m, 0;
       P2, 7;
       P3, 7;
       P3m, 0;
       P4, 7;    
     ];
 
links = [ 1 2 b n;
          2 3 b n;          
          4 5 b n;
          5 6 b n;
         ];
           
MU = 1; 
NU = 0.305; 
maxconnections=4; 
lmax = 1000; 
lmin = 100; 
coplanar_tol = 10;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
a=lmin/sqrt(3)*0.5; 
Ec = MU/(4*pi)*log(a/0.1); 
dt0=1e7;  
mobility='mobbcc0'; 

intSimTime = 0;
sinTime = 0;
dtplot=5e6;
doplot=1;% frame recording: 1 == on, 0 == off
totalSimTime = 0.3e12;
curstep = 0;
simTime = 0;

Bscrew=1e0;
Bedge=1e0;
Beclimb=1e10;
Bline=1.0e-4*min(Bscrew,Bedge);

integrator='int_trapezoid'; 
rann = 0.5*a; 
rntol = 0.5*rann; 
doremesh=1; %flat set to 0 or 1 that turns the remesh functions off or on
docollision=1; %flat set to 0 or 1 that turns collision detection off or on
doseparation=1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
dovirtmesh=1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on
plotfreq=100; 
plim=7000;
viewangle=[-35,15]; 
printfreq=1; 
printnode=2; 
rmax=100; %maximum distance a node may travel in one cycle

%FEM CANTILEVER PARAMETERS

% dx=12; %microns
% dy=2; %microns 
% dz=2; %microns
%tungsten "a" is 0.000274
dx = dx/bmag;
dy = dy/bmag;
dz = dz/bmag;

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

%plotnodes(rn,links,plim,vertices);
hold on;
plot3(rn(:,1), rn(:,2), rn(:,3),'ko'); %nodes
