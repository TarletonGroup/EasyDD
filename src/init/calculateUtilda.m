function u_tilda = calculateUtilda(rn, links, gamma, FEM, NU, u_tilda)
    %=========================================================================%
    % Wrapper for Utilda calculation by Bruce Bromage.
    %=========================================================================%
    
    %% Extraction
    
    % gamma:
    gamma_disp = gamma.disp;
    
    %% Calculation
    
    [Ux, Uy, Uz] = Utilda_bb_vec(rn, links, gamma, FEM, NU);

    u_tilda(3 * gamma_disp(:,1) - 2) = Ux;
    u_tilda(3 * gamma_disp(:,1) - 1) = Uy;
    u_tilda(3 * gamma_disp(:,1)) = Uz;
end
