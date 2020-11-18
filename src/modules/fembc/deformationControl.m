function [u, f] = deformationControl(...
    u, f, diffBC, simTime, gamma, S, fixedDofs, freeDofs)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Specifies (possibly time-dependent) applied boundary conditions
    % given the fixed and free degrees of freedom and their
    % corresponding surface sets.
    %===============================================================%
    
    %% Initialise BCs
    
    u(:) = 0; f(:) = 0;
    
    %% Force conditions
    
    f(freeDofs) = 0; % All freeDofs have zero applied force BCs
    
    %% Displacement conditions
    
    % FixedDofs corresponding to zero displacement:
    fixedDofs_0 = [3*gamma.u(:,1)-2; 3*gamma.u(:,1)-1; 3*gamma.u(:,1)];
    u(fixedDofs_0) = 0;
    
    % FixedDofs corresponding to affinely increasing displacement:
    fixedDofs_U = 3*gamma.Mixed(:, 1);
    u(fixedDofs_U) = diffBC.u_dot * simTime;
    
end
