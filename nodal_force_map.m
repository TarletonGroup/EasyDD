function [idxf, idxi, f_dln, n_nodes_t] = nodal_force_map(...
                                            se_node_label, gamma, n_neighbours)
    %%===================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 11/12/2017--15/12/2017.
    % Refactored into independent subroutines in 11/12/2017--15/12/2017.
    %---------------------------------------------------------------------%
    %
    % This subroutine generates the indices which map the nodal forces used
    % in the analytical traction calculation to an array with perfect
    % congruency to the array containing the finite element forces
    % acting on the nodes subjected to traction boundary conditions. The
    % map is non-trivial and specific to a simulation for three reasons:
    %   1) the nodes are arranged depending on how gamma is arranged,
    %   2) surface element nodes can be shared by more than one surface
    %       element (up to n_neighbours),
    %   3) some nodes in a surface element might not have traction boundary
    %       conditions but still have to be accounted for in the force
    %       calculation because the analytical solution requires knowledge
    %       of all nodes in each surface element.
    %
    %=====================================================================%
    % Inputs
    %=====================================================================%
    % 
    % se_node_label := dimension(n_elem, n_nodes). 2D array that holds the
    %   global node labels of the nodes that make up the surface elements
    %   that are to be used in the analytical traction calculation.
    %
    % gamma := dimension(:). Assumed shape array that holds the global node
    %   labels for the nodes with traction boundary conditions.
    %
    % n_neighbours := the maximum number of surface elements that can share
    %   a node (in a rectangular grid 4 elements can share a node).
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % tmp := temporary array that stores the indices where a particular
    %   node label is found in the nodal force array of the analytical
    %   traction calculation.
    %
    % j := makes sure that the output index advances 4 places onto the
    %   block of 4 entries corresponding to the next node of gamma.
    %
    %=====================================================================%
    % Outputs
    %=====================================================================%
    %
    % idxf := dimension(n_nodes_t). Precomputed array of indices
    %   corresponding to the z-coordinates of the finite element node force 
    %   array (f_dln is ordered the same).
    %
    % idxi := dimension(n_nodes*n_neighbours). Array with the indices where
    %   a given node in gamma is found in the nodal force array of the
    %   analytical traction solution.
    %
    % f_dln := dimension(3*n_nodes_t). Preallocated array with perfect
    %   congruency to the finite element node force array.
    %
    % n_nodes_t := total number of nodes with traction boundary conditions.
    %   This is not necessarily the same as the number of elements times
    %   the number of nodes per element because not all nodes in an element
    %   may have traction boundary conditions, but still need to be 
    %   included in the calculation due to the problem's formulation.
    %
    %%===================================================================%%
    
    %% Map analytical nodal forces into a useful form for the force superposition scheme.
    % Find the number of nodes for which tractions need to be calculated.
    n_nodes_t = size(gamma, 1);
    % Allocate the force vector for forces induced by dislocations.
    f_dln = zeros(3*n_nodes_t, 1);
    % The maximum number of elements that share a node.
    idxi = zeros(n_nodes_t * n_neighbours, 1);
    % The FEM forces are given in the same node order as gamma, so this
    % finds the index of the nodal forces corresponding to gamma(i).
    % Multiplied by 3 because the nodes have 3 coordinates. This lands 
    % idxi on the index which corresponds to the z-coordinate of node 
    % gamma(i).
    idxf = (1:1:n_nodes_t)'*3;
    % Loop through the number of nodes.
    j = 0;
    for i = 1: n_nodes_t
        % Find the indices of the node labels corresponding to the node
        % gamma(i). Multiplied by three because each node has 3
        % coordinates. This lands idxi on the index which corresponds to 
        % the z-coordinate of node gamma(i).
        tmp = 3*find(se_node_label == gamma(i));
        idxi(1+j: j + size(tmp, 1), 1) = tmp;
        j = j + 4;
    end %for
    
    %% Cleanup.
    clear tmp; clear j;
end %function