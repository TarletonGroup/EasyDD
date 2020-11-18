function [gamma, fixedDofs, freeDofs, t, Usim, Fsim] = ...
    cantileverBC(S, maxstepsBC)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Specifies the fixed and free degrees of freedom of the surface
    % nodes of a deforming cantilever (cuboid).
    %===============================================================%
    
    %% Specify BCs
    
    %%% Store surfaces
    
    gamma = struct; % Stores surfaces with displacement/traction BCs
    
    % Traction BCs:
    gamma.t = [S.top; S.bot; S.right; S.front; S.back; ... % Faces
        S.botright; S.frontright; S.backright; ...
        S.topfront; S.topback; S.botfront; S.botback; ... % Edges
        reshape(S.corners(2,1,1,:),1,5); reshape(S.corners(2,2,1,:),1,5)]; % Corners
    % Displacement BCs:
    gamma.u = [S.left; ... % Faces
        S.topleft; S.botleft; S.frontleft; S.backleft; ... % Edges
        reshape(S.corners(1,1,1,:),1,5); reshape(S.corners(1,2,1,:),1,5); ...
        reshape(S.corners(1,1,2,:),1,5); reshape(S.corners(1,2,2,:),1,5)]; % Corners
    % Mixed traction and displacement BCs:
    gamma.Mixed = [S.topright; ... % Edges
        reshape(S.corners(2,1,2,:),1,5); reshape(S.corners(2,2,2,:),1,5)]; % Corners
    
    %%% Specify fixed Dofs

    dofUx = 3*gamma.u(:, 1)-2;
    dofUy = 3*gamma.u(:, 1)-1;
    dofUz = [3*gamma.u(:, 1); 3*gamma.Mixed(:, 1)];
    
    fixedDofs = [dofUx; dofUy; dofUz];
    
    %%% Specify free Dofs
    
    dofTx = [3*gamma.t(:, 1)-2; 3*gamma.Mixed(:, 1)-2];
    dofTy = [3*gamma.t(:, 1)-1; 3*gamma.Mixed(:, 1)-1];
    dofTz = 3*gamma.t(:, 1);
    
    freeDofs = [dofTx; dofTy; dofTz]; 
    
    %% Initialise BC storage
    
    t = zeros(maxstepsBC, 1);
    Usim = zeros(maxstepsBC, 1);
    Fsim = zeros(maxstepsBC, 1);
    
end
