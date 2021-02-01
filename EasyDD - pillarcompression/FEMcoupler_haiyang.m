function [uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
gammau,gammat, gammaMixed,fixedDofs,freeDofs,unfixedDofs,dx,dy,dz,t,mx,my,mz,utilda_0)

%Coupling of FEM and DDD
% u = uhat + utilda
% f = fhat + ftilda

segments = constructsegmentlist(rn,links);

% Udot = 1E3*dx*(1E-4/160E9); %test by BB
Udot = 100*1E3*dx*(1E-4/160E9); %for tungsten...
% Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
% Udot = dx*1e-5; %HY

Ubar = Udot*t;
%Ubar = 0.1*1E4; for debuggin
u=zeros(3*(mno),1);
gamma=[gammau;gammaMixed];

u(3*gammaMixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

uhat=zeros(3*mno,1);
utilda=zeros(3*mno,1);

gn = gamma(:,1); % global node number
% x0 = xnodes(gn,1:3); % field point
% point_array_length = size(x0,1);
% segments_array_length = size(segments,1);

%Matlab wrapper
% tic;
% displacements = displacement_fivel(x0,segments,NU); %utilda array ordered by x0
% toc;

%Full C version (including combinatorials, checks, corrections etc.)
% tic;
% [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
%                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
%                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
%                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
%                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
%                        NU,point_array_length,segments_array_length);
% displacements = horzcat(Ux,Uy,Uz);
% toc;
% disp(displacementsMEX-displacements);
%  [Ux, Uy, Uz]  = displacement_fivel(x0,segments,NU);

% % Displacement calculation Bruce Bromage
% [Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);
% 
% utilda(3*gn -2) = Ux;
% utilda(3*gn -1) = Uy;
% utilda(3*gn   ) = Uz;
% 
% utilda = utilda - utilda_0;
% 
% if any(isnan(utilda))
%     disp('some entries of utilda are NaN')
%     pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
% end
% 
uhat(fixedDofs) = u(fixedDofs);% - utilda(fixedDofs);


fhat=zeros(3*(mno),1);

gamma=[gammat;gammaMixed];

%ftilda = zeros(3*mno,1);
ftilda = traction(gamma,segments,xnodes, mno, a, MU, NU);

%%
%ftilda=zeros(3*mno,1); %ET uncomment later!

%% we think the plane normals are wrong, lets experiment.
%fhat(freeDofs) = -ftilda(freeDofs);% no applied forces
fhat(freeDofs) = ftilda(freeDofs);% no applied forces
%% 

%f=zeros(2*(mno),1);
f=fhat-kg(:,fixedDofs)*uhat(fixedDofs);

%HY20171206:********************************************************
%HY20171206: modified by HY to make the code cleaner by removing the
%equations related to the fixedDofs; since FreeDofs has been used to
%represent the free boundary nodes, a new term, unfixedDofs, is used to
%represent all the nodes other than the fixed ones. i.e.
%unfixedDofs=allDofs - fixedDofs
% 
uhat(unfixedDofs) = U\(L\f(unfixedDofs));

%HY20171206:********************************************************

%HY20171206: commented by HY
% bcwt=mean(diag(kg));%=trace(K)/length(K)
% bcwt = full(bcwt);
%
% f(fixedDofs) = bcwt*uhat(fixedDofs);
% uhat = U\(L\f); %using LU decomposition
% uhat2=K\f;

rhat=kg*uhat; % reaction force

fend = rhat(3*gammaMixed(:,1))+ftilda(3*gammaMixed(:,1));
fend = sum(fend);

end
