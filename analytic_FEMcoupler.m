function [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
    gamma_disp, gamma_mixed, fixedDofs,freeDofs,dx,dy,dz,mx,my,mz,utilda_0,t,...
    gamma_dln, x3x6, n_nodes, n_nodes_t, n_se, idxi, ...
    f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme)
    
%Coupling of FEM and DDD
% u = uhat + utilda
% f = f_hat + ftilda

%segments = constructsegmentlist(rn,links);

Udot = 100*1E3*dx*(1E-4/160E9); %for tungsten...
% Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle 

Ubar = Udot*t;
%Ubar = 0.1*1E4; for debucontourfggin
u=zeros(3*(mno),1);

u(3*gamma_mixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes  

uhat=zeros(3*mno,1);
utilda=zeros(3*mno,1); 

gn = gamma_disp(:,1); % global node number
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
[Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz); %must be used with virtual segments projected normal to the surface

utilda(3*gn -2) = Ux;
utilda(3*gn -1) = Uy;
utilda(3*gn   ) = Uz;

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
                        MU, NU, a, use_gpu, n_threads, para_scheme);
                    
%f_dln = - f_dln;
% any(isnan(f_dln))
% any(isinf(f_dln))
%fhat=zeros(3*mno,1);
%% fhat2 replaced by f_hat in the rest of the code.
% fhatn=zeros(3*mno,1);

% ftilda = traction([gammat; gamma_mixed],segments,xnodes, mno, a, MU, NU);

%%
% f_hat(freeDofs) = -f_dln(freeDofs);% no applied forces
f_hat(freeDofs) = f_dln(freeDofs)';% no applied forces
% fhatn(freeDofs) = -ftilda(freeDofs);% no applied forces

% plot(f_hat - fhatn)

% plot(rel_err, '.')
%f=zeros(2*(mno),1);
% f_hat = f_hat-kg(:,fixedDofs)*uhat(fixedDofs);
f2    = f_hat-kg(:,fixedDofs)*uhat(fixedDofs);
%f     = fhat -kg(:,fixedDofs)*uhat(fixedDofs);

bcwt =mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

% f_hat(fixedDofs) = bcwt*uhat(fixedDofs);
% uhat = U\(L\f_hat); %using LU decomposition
% uhat2=K\f;
f2(fixedDofs) = bcwt*uhat(fixedDofs);
uhat = U\(L\f2); %using LU decomposition
%f(fixedDofs) = bcwt*uhat(fixedDofs);
%uhat = U\(L\f); %using LU decomposition

rhat2=kg*uhat;
%rhat=kg*uhat; % reaction force

%% The next two lines used to be fend2 instead of fend. Changed it so Bruce
%% can run simulations with my traction code. This whole function needs to 
%% be edited further to remove unused variables and optimise it.
fend = rhat2(3*gamma_mixed(:,1))-f_dln(3*gamma_mixed(:,1));
fend = sum(fend);
%%
%fend = rhat(3*gamma_mixed(:,1))+ftilda(3*gamma_mixed(:,1));
%fend = sum(fend);

%% for plotting diferences between analtical and numerical traction
% idx = f_dln ~= 0;
% tf_dln = f_dln(idx);
% tftilda = ftilda(idx);
% idx2 = tftilda ~= 0;
% tf_dln2 = tf_dln(idx2);
% tftilda2 = tftilda(idx2);
% abs_err = tftilda2 + tf_dln2;
% rel_err = (tftilda2 + tf_dln2)./tf_dln2;
% min(rel_err)
% max(rel_err)
% mean(rel_err)
 
% fig1 = figure;
% subplot(1,1,1)
% plot(rel_err, '.')
%  xlabel('Surface Node, arbitrary ordering' , 'Interpreter', 'latex');
%  ylabel('Force Relative Error, $\frac{F_{\mathrm{n}}-F_{\mathrm{a}}}{F_{\mathrm{a}}}$', 'Interpreter', 'latex');
end