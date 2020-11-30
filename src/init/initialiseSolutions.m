function [u, u_hat, u_tilda, u_tilda_0, f, f_hat, f_tilda, fseg] = ...
    initialiseSolutions(FEM)
    %===============================================================%
    % Daniel Hortelano Roig (22/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Initialise main simulation data (not FEM data).
    %===============================================================%
    
    %% Extraction
    
    ndofs = FEM.ndofs;
    
    %% Initialise data
    
    % Nodal displacements:
    u = zeros(ndofs, 1);
    u_hat = zeros(ndofs, 1);
    u_tilda = zeros(ndofs, 1);
    u_tilda_0 = zeros(ndofs, 1);
    
    % Nodal forces:
    f = zeros(ndofs, 1);
    f_hat = zeros(ndofs, 1);
    f_tilda = zeros(ndofs, 1);
    
    % Segment forces:
    fseg = [];
    
end
