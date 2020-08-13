% TODO: add Udot as an input. Optional to turn displacements on or off. Add option for non-zero initial load.
function [u_hat,fend,Ubar] = FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
    gammau,gammat, gammaMixed,fixedDofs,freeDofs,dx,dy,dz,t,mx,my,mz,u_tilda_0)

%Coupling of FEM and DDD
% u = u_hat + u_tilda
% f = fhat + ftilda

[segments,~] = constructsegmentlist(rn,links,1);

% Udot = 1E3*dx*(1E-4/160E9); %test by BB
Udot = (1/2048)*100*1E3*dx*(1E-4/160E9)*2048*100; %for tungsten...
% Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
% Udot = dx*1e-5; %HY

Ubar = Udot*t;
u=zeros(3*(mno),1);
gamma=[gammau;gammaMixed];

u(3*gammaMixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

u_hat=zeros(3*mno,1);
u_tilda=zeros(3*mno,1);

gn = gamma(:,1); % global node number

[Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);

u_tilda(3*gn -2) = Ux;
u_tilda(3*gn -1) = Uy;
u_tilda(3*gn   ) = Uz;

u_tilda = u_tilda - u_tilda_0;

if any(isnan(u_tilda))
    disp('some entries of u_tilda are NaN')
    pause; %some entries of u_tilda are NaN -> leads to issues in hatStress routine
end

u_hat(fixedDofs) = u(fixedDofs) - u_tilda(fixedDofs);

fhat=zeros(3*(mno),1);

gamma=[gammat;gammaMixed];

ftilda = traction(gamma,segments,xnodes, mno, a, MU, NU);

fhat(freeDofs) = -ftilda(freeDofs);% no applied forces

f=fhat-kg(:,fixedDofs)*u_hat(fixedDofs);

bcwt = mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

f(fixedDofs) = bcwt*u_hat(fixedDofs);
u_hat = U\(L\f); %using LU decomposition

rhat=kg*u_hat; % reaction force

fend = rhat(3*gammaMixed(:,1))+ftilda(3*gammaMixed(:,1));
fend = sum(fend);

end
