function [u, u_hat, u_tilda, f, f_hat, f_tilda, r_hat] = FEM_DDD_Superposition(...
        rn, links, ... % Network variables
        simTime, ... % Time
        u, u_hat, u_tilda, u_tilda_0, ... % Displacements
        f, f_hat, f_tilda, ... % Forces
        matpara, mods, flags, decompK, FEM, gamma, dofs, tract, diffBC) % Structures
    %=====================================================================%
    % Oxford Materials (11/11/2020)

    % Couples FEM and DDD using the superposition principle method (SPM).
    %=====================================================================%
    
    %% Extraction
    
    % simpara:
    a = matpara.a;
    MU = matpara.MU;
    NU = matpara.NU;
    
    % flags:
    a_trac = flags.a_trac;
    CUDA_flag = flags.CUDA_flag;
    calculateR_hat = flags.calculateR_hat;
    
    % Kdecomp (global stiffness matrix decompositions):
    Lchol = decompK.Lchol;
    Uchol = decompK.Uchol;
    bcwt = decompK.bcwt;
    
    % DoFs:
    fixedDofs = dofs.fixedDofs;
    freeDofs = dofs.freeDofs;
    
    % FEM:
    kg = FEM.kg;
    xnodes = FEM.xnodes;
    
    %% Compute dislocation fields
    
    % Reset displacements and forces.
    u_hat(:) = 0; u_tilda(:) = 0;
    f_hat(:) = 0; f_tilda(:) = 0;
    
    % Calculate adjusted U_tilda.
    u_tilda = calculateUtilda(rn, links, gamma, FEM, NU, u_tilda) - u_tilda_0;
    
    % Calculate dislocation forces on FE nodes.
    if a_trac
        [dln_node_coord, b, n_dln] = extract_dislocation_nodes(rn, links);
        [f_tilda, ~] = analytic_traction(tract, gamma, ...
            f_tilda, dln_node_coord, b, n_dln, MU, NU, a, CUDA_flag);
    else
        [segments, ~] = constructsegmentlist(rn, links, true);
        f_tilda = traction(gamma.dln, segments, xnodes, a, MU, NU, f_tilda);
    end
    
    %% Specify domain boundary conditions
    
    [u, f] = feval(mods.boundaryConditions, ...
        u, f, diffBC, simTime, gamma, dofs);
    
    %% Compute corrective fields
    
    % Calculate corrective displacements and forces.
    u_hat(fixedDofs) = u(fixedDofs) - u_tilda(fixedDofs);
    f_hat(freeDofs) = f(freeDofs) - f_tilda(freeDofs);
    
    % Modify corrective forces to apply BCs on fixed Dofs.
    f_temp = f_hat - kg(:, fixedDofs) * u_hat(fixedDofs);
    f_temp(fixedDofs) = bcwt * u_hat(fixedDofs);
    
    % Solve for corrective displacements.
    u_hat = Uchol \ (Lchol \ f_temp); % Using LU decomposition
    
    %% Compute reactive fields
    
    % If not using force control, calculate reaction force.
    if calculateR_hat
        r_hat = kg * u_hat;
    else % Else, don't calculate reaction force.
        r_hat = 0;
    end
        
end
