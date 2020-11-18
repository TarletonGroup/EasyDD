function [f, f_hat, f_tilda, u, u_hat, u_tilda, r_hat, BCstore] = FEM_DDD_Superposition(...
        rn, links, a, MU, NU, xnodes, kg, L, U, bcwt, gamma_disp, gamma, S, ...
        fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, diffBC, ...
        u_tilda_0, u, u_hat, u_tilda, loading, BCstore, a_trac, gamma_dln, x3x6, n_nodes, ...
        n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, f_hat, use_gpu, ...
        n_threads, para_scheme, para_tol)
    %=====================================================================%
    % Oxford Materials (as of 11/11/2020)

    % Couples FEM and DDD using the superposition principle method (SPM).
    %=====================================================================%
    
    %% Reset fields
    
    % Reset displacements and forces.
    u_hat(:) = 0; u_tilda(:) = 0;
    f_hat(:) = 0; f_tilda(:) = 0;
    
    %% Compute dislocation fields
    
    % Calculate adjusted U_tilda.
    u_tilda = calculateUtilda(rn, links, gamma_disp, NU, xnodes, dx, ...
        dy, dz, mx, my, mz, u_tilda) - u_tilda_0;
    
    % Calculate dislocation forces on FE nodes.
    if a_trac
        [x1x2, b, n_dln] = extract_dislocation_nodes(rn, links);
        [f_tilda, ~] = analytic_traction(x3x6, x1x2, b, n_nodes, n_nodes_t, ...
            n_se, n_dln, 3 * gamma_dln(:, 1), idxi, ...
            f_tilda_node, f_tilda_se, f_tilda, ...
            MU, NU, a, use_gpu, n_threads, para_scheme, para_tol);
    else
        [segments, ~] = constructsegmentlist(rn, links, true);
        f_tilda = traction(gamma_dln, segments, xnodes, a, MU, NU, f_tilda);
    end
    
    %% Specify domain boundary conditions
    
    [u, f] = feval(loading, ...
        u, f, diffBC, simTime, gamma, S, fixedDofs, freeDofs);
    
    %% Compute corrective fields
    
    % Calculate corrective displacements and forces.
    u_hat(fixedDofs) = u(fixedDofs) - u_tilda(fixedDofs);
    f_hat(freeDofs) = f(freeDofs) - f_tilda(freeDofs);
    
    % Modify corrective forces to apply BCs on fixed dofs
    f_temp = f_hat - kg(:, fixedDofs) * u_hat(fixedDofs);
    f_temp(fixedDofs) = bcwt * u_hat(fixedDofs);
    
    % Solve for corrective displacements.
    u_hat = U \ (L \ f_temp); % Using LU decomposition
    
    %% Compute reactive fields
    
    % If not using force control, calculate reaction force.
    if calculateR_hat
        r_hat = kg * u_hat;
    else % If using force control, don't calculate reaction force.
        r_hat = 0;
    end
    
    % Print and store relevant BCs for plotting and information.
    [BCstore] = feval(processForceDisp, ...
        BCstore, f, f_hat, f_tilda, u, u_hat, u_tilda, r_hat, ...
        gamma, fixedDofs, freeDofs, curstep, simTime);
    
end
