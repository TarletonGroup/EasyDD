% TODO: add Udot as an input. Optional to turn displacements on or off. Add option for non-zero initial load.
function [uhat,Ftop,Uend,aveuhat,aveutilda] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
   Stop,gammau,gammat, gammaMixed,fixedDofs,freeDofs,dx,dy,dz,curstep,mx,my,mz,unfixedDofs,Kred,Lred,Ured,utilda_0)

%Coupling of FEM and DDD
% u = uhat + utilda
% f = fhat + ftilda

global mumag
segments = constructsegmentlist(rn,links);

% Udot = 1E3*dx*(1E-4/160E9); %test by BB
%Udot = (1/2048)*1000000*1E3*dx*(1E-4/160E9); %for tungsten...
% Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle 
 %Udot = dx*1e-10; %HY
% Fdot = -600/mumag*1e-3;
% t= 2*10^1*1e3* t*(1/mumag); % B=10^3 Pa.s
% 
% 
% Fbar = Fdot*t;
 Fbar = -200/mumag*dx*dy/size(Stop,1)/100*curstep;

fbar=zeros(3*(mno),1);
fbar(Stop(:,1)*3) = Fbar;
%Ftop = sum(fbar(Stop(:,1)*3).*Stop(:,2)/max(Stop(:,2)));
Ftop = sum(fbar(Stop(:,1)*3));


%Ubar = Udot*t;
%Ubar = 0.1*1E4; for debuggin
u=zeros(3*(mno),1);
gamma=[gammau;gammat];

%u(3*gammaMixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes  

uhat=zeros(3*mno,1);
utilda=zeros(3*mno,1); 

gn = gamma(:,1); % global node number
% gn(:,1)= 1:mno;

[Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);

utilda(3*gn -2) = Ux;
utilda(3*gn -1) = Uy;
utilda(3*gn   ) = Uz;

utilda = utilda - utilda_0;

if any(isnan(utilda))
    disp('some entries of utilda are NaN')
    pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
end

uhat(fixedDofs) = u(fixedDofs) - utilda(fixedDofs);

fhat=zeros(3*(mno),1);

gamma=[gammat;gammaMixed];

%HY20200602:
%find segments that are virtual with a surface node.
% [virtindex,~,realindex,realID,~] = virtualsegfinder(rn,links,segments);
% virtsegs = segments(virtindex,:);

S=size(segments,1);
%get only "real" segments, not virtual ones.
index=true(S,1);
for i=1:S
    if rn(segments(i,1),4) == 67 || rn(segments(i,2),4) == 67
        index(i) = false;
    end
end
segments = segments(index,:);

%ftilda = zeros(3*mno,1);
ftilda = traction(gamma,segments,xnodes, mno, a, MU, NU);   
%ftilda=zeros(3*mno,1); %ET uncomment later!

fhat(freeDofs) = fbar(freeDofs)-ftilda(freeDofs);% no applied forces

%f=zeros(2*(mno),1);
f=fhat-kg(:,fixedDofs)*uhat(fixedDofs);

%%%%%%%%%%%%%%%%%%%%%%%%%

%HY20171206: modified by HY to make the code cleaner by removing the
%equations related to the fixedDofs; since FreeDofs has been used to
%represent the free boundary nodes, a new term, unfixedDofs, is used to
%represent all the nodes other than the fixed ones. i.e.
%unfixedDofs=allDofs - fixedDofs
% tic;
% disp('setdiff')
% uhat = uhat;
%uhatHY = uhat;
fred = f(unfixedDofs);
u_new = Ured\(Lred\fred);
uhat(unfixedDofs) = u_new;

rhat=kg*uhat; % reaction force

 uend = uhat(3*Stop(:,1))+utilda(3*Stop(:,1));
% uend = uhat(3*Stop(:,1));
Uend = mean(uend);
%fend = sum(fend);
aveuhat = mean(uhat(3*Stop(:,1)));
aveutilda = mean(utilda(3*Stop(:,1)));
end
