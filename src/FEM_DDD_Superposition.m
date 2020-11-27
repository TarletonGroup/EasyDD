function [f_bar, f_hat, f_tilda, u_bar, u_hat, u_tilda, r_hat] = FEM_DDD_Superposition(...
        rn, links, a, MU, NU, xnodes, kg, K, L, U, P_l, P_u, gamma_disp, gamma_mixed, fixedDofs, ...
        freeDofs, dx, dy, dz, simTime, mx, my, mz, sign_u_dot, u_dot, sign_f_dot, ...
        f_dot, u_tilda_0, u, u_hat, u_tilda, loading, calculateTractions, gamma_dln, ...
        x3x6, n_nodes, n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, ...
        f_hat, use_gpu, n_threads, para_scheme, para_tol)

    u(:, 1) = 0;
    u_hat(:, 1) = 0;
    u_tilda(:, 1) = 0;
    f(:, 1) = 0;
    f_hat(:, 1) = 0;
    f_tilda(:, 1) = 0;

    [u, u_bar, f, f_bar, calculateR_hat] = loading(u, u_dot, sign_u_dot, f, f_dot, sign_f_dot, simTime, gamma_mixed);

    % Calculate adjusted U_tilda.
    u_tilda = calculateUtilda(rn, links, gamma_disp, NU, xnodes, dx, ...
        dy, dz, mx, my, mz, u_tilda) - u_tilda_0;

    u_hat(fixedDofs) = u(fixedDofs) - u_tilda(fixedDofs);

    f_tilda = calculateTractions(rn, links, gamma_dln, xnodes, a, MU, NU, ...
        f_tilda, x3x6, n_nodes, n_nodes_t, n_se, idxi, f_tilda_node, ...
        f_tilda_se, use_gpu, n_threads, para_scheme, para_tol);

    f_hat(freeDofs) = f(freeDofs) - f_tilda(freeDofs); % no applied forces
    f = f_hat - kg(:, fixedDofs) * u_hat(fixedDofs);

    bcwt = mean(diag(kg)); %=trace(K)/length(K)
    bcwt = full(bcwt);

    f(fixedDofs) = bcwt * u_hat(fixedDofs);

    if isempty(U)
        u_hat = K \ f;
    else
        u_hat = P_l \ (U \ (L \ (P_u \ f))); % using LU decomposition for sparse matrices
    end

    % If not using force control, calculate reaction force.
    if calculateR_hat
        r_hat = kg * u_hat;
        % If using force control, don't calculate reaction force.
    else
        r_hat = 0;
    end

end
