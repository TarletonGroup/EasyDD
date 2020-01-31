%%
% Static comparison of tractions+FEM vs analytic image forces.
%
%
%% SOURCE GENERATION PARAMETERS
amag=3.18e-4; 
mumag = 145E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
%% FEM PARAMETERS
%Cantilever
simTime = 0;
use_gpu=0;
para_scheme = 0;
n_threads = 0;
dx=15/amag;
dy=15/amag;
dz=15/amag;

mx=5; % number of elements along beam length
loading=1; 
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];


%MU = 160E9; %160GPa
MU = 1;
NU = 0.28; 
lmin = 0.1/amag;
a=lmin/sqrt(3)*0.5; 
rn = [0, 0.5, 0.5;
      1, 0.5, 0.5]*dx;
links = [1, 2, 0, 1, 1,  1, -1,  1];
l = (rn(2,:)-rn(1,:))/norm(rn(2,:)-rn(1,:));



plim=12/amag; %12microns
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,Smixed,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = STATIC_finiteElement3D(dx,dy,dz,mx,MU,NU,loading);

 gammad_dln=[Stop;Sbot;Sright;Sleft;Sfront;Sback;Smixed]; % t=0         
    
Dofs = [3*gammad_dln(:,1)-2; 3*gammad_dln(:,1)-1; 3*gammad_dln(:,1)];

%% Addition by Daniel Celis Garza 12/19/2017.
% Precomputing variables needed to couple tractions induced by dislocations
% to FEM. Generate the list of planes to be extracted. The FEM coupler is 
% hardcoded so that the min(x) yz-plane does not have traction boundary
% conditions.
f_hat = zeros(3*mno, 1);

planes = (1:1:6)';
[x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
                                               planes, 4);
gamma_dln = Dofs;%[gammat(:,1); gammaMixed(:,1)];
[f_dln_node, f_dln_se,...
 f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
tolerance = dx/10^7;


[uhat,fend,Ubar] = STATIC_analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
                       Dofs, Dofs, Dofs,Dofs,Dofs,dx,simTime ,...
                       gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
                       f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);
segments=constructsegmentlist(rn,links);
image_force = pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);

function f=pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments)
%nodal force on dislocation segments due to applied stress sigext
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,nx,ny,nz)
[nseg,~]=size(segments);
f=zeros(nseg,3);

b=segments(:,3:5);
r0=segments(:,6:8);
r1=segments(:,9:11);
r01=r1-r0;
rmid = 0.5*(r0+r1);
for i =1:nseg % evaluate sigmahat at midpoint of each segment (do something better later)
    xi = rmid(i,1:3);

    if any(isnan(xi))
        disp('Segment midpoint is undefined! See segforcevec.m');
        pause;
    end

    sigext = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,xi); %must not pass virtsegs!

    sigb=sigext*b(i,1:3)';
    l = r01(i,:);
    f(i,:)=cross(sigb',l);
end
end