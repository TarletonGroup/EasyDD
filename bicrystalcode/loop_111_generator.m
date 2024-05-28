clear all

global Bscrew Bedge Beclimb Bline

l = 300;
b = [1 1 1];

side1 = [1 -2 1];
side2 = [2 -1 -1];
side3 = [1 1 -2];
side4 = [-1 2 -1];
side5 = [-2 1 1];
side6 = [-1 -1 2];

n1 = cross(side1,b);
n2 = cross(side2,b);
n3 = cross(side3,b);
n4 = cross(side4,b);
n5 = cross(side5,b);
n6 = cross(side6,b);

%start from center of grid

P1 = [10000 2000 500];
P2 = P1 + l*side1;
P3 = P2 + l*side2;
P4 = P3 + l*side3;
P5 = P4 + l*side4;
P6 = P5 + l*side5;

P7 = P6 + l*side6; %should equal P1, consistency check

rn = [ 
       P1, 0;
       P2, 0;
       P3, 0;
       P4, 0;
       P5, 0;
       P6, 0;
     ];
 
links = [ 1 2 b n1;
          2 3 b n2;
          3 4 b n3;
          4 5 b n4;
          5 6 b n5;
          6 1 b n6;
         ];
           
MU = 1; 
NU = 0.305; 
maxconnections=8; 
lmax = 1000; 
lmin = 100; 
coplanar_tol = 10;
areamin=lmin*lmin*sin(60/180*pi)*0.5; 
areamax=20*areamin; 
a=lmin/sqrt(3)*0.5; 
Ec = MU/(4*pi)*log(a/0.1); 
dt0=1e5;  
mobility='mobbcc0'; 

intSimTime = 0;
sinTime = 0;
dtplot=3e6;% 3 ms (code in units of ns)
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
plotfreq=1; 
plim=7000;
appliedstress =10^-2*[2 0 1; 0 2 -1; 1 -1 0];
viewangle=[90 0 ]; 
printfreq=1; 
printnode=3; 
rmax=100; %maximum distance a node may travel in one cycle

%FEM CANTILEVER PARAMETERS

dx=4; %microns
dy=1; %microns 
dz=1; %microns
%tungsten "a" is 0.000274
dx = dx/0.000274;
dy = dy/0.000274;
dz = dz/0.000274;

mx=6;
loading=1; 
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];
