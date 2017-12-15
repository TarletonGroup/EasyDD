function [idxf, idxi, f_dln, n_nodes_t] = nodal_force_map(se_node_label, gamma, n_neighbours)
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