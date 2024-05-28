% TODO: add Udot as an input. Optional to turn displacements on or off. Add option for non-zero initial load.

function [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
    gamma_disp, gammat, gamma_mixed, fixedDofs,freeDofs,dx,dy,dz,t,mx,my,mz,utilda_0,...
    gamma_dln, x3x6, n_nodes, n_nodes_t, n_se, idxi, ...
    f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, para_tol)

% % % 
% Coupling of FEM and DDD
% u = uhat + utilda
% f = f_hat + ftilda
% segments = constructsegmentlist(rn,links);
% Udot = (1/2048)*100*1E3*dx*(1E-4/160E9)*2048*100; %for tungsten...
% Udot = 100*1E3*dx*(1E-4/160E9); %for tungsten...
Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
Ubar = Udot*t;

% Ubar0 = 4.25e2;
% Udot = (1/2048)*100*1E3*dx*(1E-4/160E9); %for tungsten...
% Ubar = Udot*t + Ubar0;

%Ubar = 0.1*1E4; for debucontourfggin
u=zeros(3*(mno),1);

u(3*gamma_mixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

uhat=zeros(3*mno,1);
utilda=zeros(3*mno,1);

gn = gamma_disp(:,1); % global node number

[Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz); %must be used with virtual segments projected normal to the surface

utilda(3*gn-2) = Ux;
utilda(3*gn-1) = Uy;
utilda(3*gn  ) = Uz;

utilda = utilda - utilda_0;

if any(isnan(utilda))
    disp('some entries of utilda are NaN')
    pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
end

uhat(fixedDofs) = u(fixedDofs) - utilda(fixedDofs);

[x1x2, b, n_dln] = extract_dislocation_nodes(rn, links);
f_dln(:,1) = 0;
f_hat(:,1) = 0;
[f_dln,~] = analytic_traction(x3x6 , x1x2, b, n_nodes, n_nodes_t,...
                        n_se, n_dln, 3*gamma_dln, idxi,...
                        f_dln_node, f_dln_se, f_dln,...
                        MU, NU, a, use_gpu, n_threads, para_scheme, para_tol);

f_hat(freeDofs) = -f_dln(freeDofs)';% no applied forces

f    = f_hat-kg(:,fixedDofs)*uhat(fixedDofs);

bcwt = mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

f(fixedDofs) = bcwt*uhat(fixedDofs);
uhat = U\(L\f); %using LU decomposition

rhat=kg*uhat;

fend = rhat(3*gamma_mixed(:,1))+f_dln(3*gamma_mixed(:,1));
fend = sum(fend);
end
