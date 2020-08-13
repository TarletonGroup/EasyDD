function [f_hat, u_hat, r_hat] = FEM_DDD_Superposition(rn, links, a, MU, NU,...
    xnodes, mno, kg, L, U, gamma_disp, gammat, gamma_mixed, fixedDofs,...
    freeDofs, dx, dy, dz, t, mx, my, mz, u_dot, f_dot, u_tilda_0, u,...
    u_hat, u_tilda, simType, a_trac, gamma_dln, x3x6, n_nodes,...
    n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, f_hat, use_gpu,...
    n_threads, para_scheme, para_tol)


u(:,1) = 0;
u_hat(:,1) = 0;
u_tilda(:,1) = 0;
f(:,1) = 0;
f_hat(:,1) = 0;
f_tilda(:,1) = 0;

% Displacement control.
if simType == 1
    u_bar = u_dot*t;
    u(3*gamma_mixed(:,1)) = -u_bar;
    % Force control.
elseif simType == 2
    f_bar = f_dot*t;
    f(3*gamma_mixed(:,1)) = -f_bar;
else
    u_bar = u_dot*t;
    u(3*gamma_mixed(:,1)) = -u_bar;
    f_bar = f_dot*t;
    f(3*gamma_mixed(:,1)) = -f_bar;
end

% Calculate adjusted U_tilda.
u_tilda = calculateUtilda(rn, links, gamma_disp, NU, xnodes, dx,...
    dy, dz, mx, my, mz, u_tilda) - u_tilda_0;

u_hat(fixedDofs) = u(fixedDofs) - u_tilda(fixedDofs);

if a_trac == true
    [x1x2, b, n_dln] = extract_dislocation_nodes(rn, links);
    [f_tilda,~] = analytic_traction(x3x6 , x1x2, b, n_nodes, n_nodes_t,...
        n_se, n_dln, 3*gamma_dln(:,1), idxi,...
        f_tilda_node, f_tilda_se, f_tilda,...
        MU, NU, a, use_gpu, n_threads, para_scheme, para_tol);
else
    f_tilda = traction(gamma_dln,segments,xnodes, a, MU, NU, f_tilda);  
end

f_hat(freeDofs) = f(freeDofs) - f_tilda(freeDofs);% no applied forces
f    = f_hat-kg(:,fixedDofs)*u_hat(fixedDofs);

bcwt = mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

f(fixedDofs) = bcwt*u_hat(fixedDofs);
u_hat = U\(L\f); %using LU decomposition

r_hat=kg*u_hat;

end