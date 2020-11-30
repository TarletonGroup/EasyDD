function [gamma, dofs] = cantilever_bending(S)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk
    
    % Specifies the fixed and free degrees of freedom of the surface
    % nodes of a deforming cantilever (cuboid). In addition,
    % complimentary storing and plotting modules are defined here.
    %===============================================================%
    
    %% Specify BCs
    
    %%% Specify surfaces
    
    gamma = struct; % Stores boundary node information
    
    % Traction-only BCs:
    gamma.t = [S.top; S.bot; S.right; S.front; S.back; ... % Faces
        S.botright; S.frontright; S.backright; ...
        S.topfront; S.topback; S.botfront; S.botback; ... % Edges
        S.corners(2,:); S.corners(4,:)]; % Corners
    % Displacement-only BCs:
    gamma.u = [S.left; ... % Faces
        S.topleft; S.botleft; S.frontleft; S.backleft; ... % Edges
        S.corners(1,:); S.corners(3,:); ...
        S.corners(5,:); S.corners(7,:)]; % Corners
    % Mixed traction and displacement BCs:
    gamma.Mixed = [S.topright; ... % Edges
        S.corners(6,:); S.corners(8,:)]; % Corners
    
    % Nodes with displacement BCs:
    gamma.disp = [gamma.u; gamma.Mixed];
    
    % Nodes with traction BCs:
    gamma.dln = [gamma.t; gamma.Mixed];
    
    %%% Specify DoFs
    
    dofs = struct; % Stores DoFs
    
    % Fixed DoFs:
    dofUx = 3*gamma.u(:, 1)-2;
    dofUy = 3*gamma.u(:, 1)-1;
    dofUz = [3*gamma.u(:, 1); 3*gamma.Mixed(:, 1)];
    dofs.fixedDofs = [dofUx; dofUy; dofUz];
    
    % Free DoFs:
    dofTx = [3*gamma.t(:, 1)-2; 3*gamma.Mixed(:, 1)-2];
    dofTy = [3*gamma.t(:, 1)-1; 3*gamma.Mixed(:, 1)-1];
    dofTz = 3*gamma.t(:, 1);
    dofs.freeDofs = [dofTx; dofTy; dofTz];
    
end
