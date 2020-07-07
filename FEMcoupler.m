% TODO: add Udot as an input. Optional to turn displacements on or off. Add option for non-zero initial load.
function [uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
    gammau,gammat, gammaMixed,fixedDofs,freeDofs,dx,dy,dz,t,mx,my,mz,utilda_0)

%Coupling of FEM and DDD
% u = uhat + utilda
% f = fhat + ftilda

segments = constructsegmentlist(rn,links,1);

% Udot = 1E3*dx*(1E-4/160E9); %test by BB
Udot = (1/2048)*100*1E3*dx*(1E-4/160E9)*2048*100; %for tungsten...
% Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
% Udot = dx*1e-5; %HY

Ubar = Udot*t;
u=zeros(3*(mno),1);
gamma=[gammau;gammaMixed];

u(3*gammaMixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

uhat=zeros(3*mno,1);
utilda=zeros(3*mno,1);

gn = gamma(:,1); % global node number

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

ftilda = traction(gamma,segments,xnodes, mno, a, MU, NU);

fhat(freeDofs) = -ftilda(freeDofs);% no applied forces

f=fhat-kg(:,fixedDofs)*uhat(fixedDofs);

bcwt = mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

f(fixedDofs) = bcwt*uhat(fixedDofs);
uhat = U\(L\f); %using LU decomposition

rhat=kg*uhat; % reaction force

fend = rhat(3*gammaMixed(:,1))+ftilda(3*gammaMixed(:,1));
fend = sum(fend);

end
