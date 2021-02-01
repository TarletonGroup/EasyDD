function [coords, burgers, n_dln] = extract_dislocation_nodes(...
                                        dln_nodes, dln_node_cnct)
    %%===================================================================%%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 11/12/2017--15/12/2017.
    % Refactored into independent subroutines in 11/12/2017--15/12/2017.
    %---------------------------------------------------------------------%
    %
    % This subroutine extracts the node coordinates of individual 
    % dislocation line segments for the analytical traction calculation
    % found in:
    % S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014)
    %
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % dln_nodes = dimension(:, 4). Assumed shape 2D array with 4 columns.
    %   Describes the nodes of the dislocation ensemble and mobility.
    %
    % dln_node_cnct = dimension(:, 8). Assumed shape 2D array with 8
    %   columns. Describes the node connectivity of the dislocation ensemble,
    %   associated Burgers vectors and slip planes.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % dln := dislocation line segment list. Used to reshape it into the
    %   form that we need it.
    %
    %=====================================================================%
    % Outputs
    %=====================================================================%
    %
    % coords := dimension(n_dln*3, 2). Contains the Cartesian coordinates
    %   of the dislocation line segments' nodes.
    %
    % burgers := dimension(n_dln). Contains the Burgers vectors
    %   corresponding to the dislocation line segments.
    %
    % n_dln := number of dislocation line segments.
    %
    %%===================================================================%%
    %% Generate dislocation line segments.
    %dln   = constructsegmentlist(dln_nodes, dln_node_cnct)';
    dln = construct_segment_list_analytic_tractions(dln_nodes, dln_node_cnct)';
    n_dln = size(dln, 2);

    %% Extract dislocation line nodes and burgers vectors.
    % Preallocate the coordinate array for both nodes.
    coords = zeros(3 * n_dln, 2);
    % Make sure the arrays are arranged in a way that can be used by the
    % traction calculation.
    coords(:, 1) = reshape(dln(6:8 , :), 3 * n_dln, 1);
    coords(:, 2) = reshape(dln(9:11, :), 3 * n_dln, 1);
    burgers      = reshape(dln(3:5,  :), 3 * n_dln, 1);
    
    %% Cleanup.
    clear dln;
end %function