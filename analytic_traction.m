function [nodal_force, total_force] = analytic_traction(              ...
                                        se_node_coord , se_node_label,...
                                        dln_node_coord, burgers      ,...
                                        n_nodes, n_se, n_dln, gamma  ,...
                                        mu, nu, a, use_gpu, n_threads,...
                                        para_scheme)
    %% Force calculation.
    % Allocating nodal and total force arrays
    nodal_force = zeros(3*n_se, n_nodes);
    total_force = zeros(3*n_se, 1);

    %% CUDA C calculation
    if use_gpu == 1
        % Provide a default number of threads in case none is given.
        if ~exists(n_threads)
            n_threads = 512;
        end %if
        % Provide a default parallelisaion scheme in case none is given.
        if ~exists(para_scheme)
            % Parallelise over dislocations.
            para_scheme = 1;
        end %if
        [nodal_force(:, 1), nodal_force(:, 2),...
         nodal_force(:, 3), nodal_force(:, 4),...
         total_force] = nodal_surface_force_linear_rectangle_mex_cuda(        ...
                                dln_node_coord(:, 1), dln_node_coord(:, 2)   ,...
                                se_node_coord (:, 1), se_node_coord (:, 2)   ,...
                                se_node_coord (:, 3), se_node_coord (:, 4)   ,...
                                burgers(:), mu, nu, a, n_se, n_dln, n_threads,...
                                para_scheme);
    % Serial force calculation in C.
    elseif use_gpu == 0
        [nodal_force(:, 1), nodal_force(:, 2),...
         nodal_force(:, 3), nodal_force(:, 4),...
         total_force] = nodal_surface_force_linear_rectangle_mex(          ...
                                dln_node_coord(:, 1), dln_node_coord(:, 2),...
                                se_node_coord (:, 1), se_node_coord (:, 2),...
                                se_node_coord (:, 3), se_node_coord (:, 4),...
                                burgers(:), mu, nu, a, n_se, n_dln);
    % Matlab version.
    else
        [nodal_force(:, 1), nodal_force(:, 2),...
         nodal_force(:, 3), nodal_force(:, 4),...
         total_force] = NodalSurfForceLinearRectangle2(                    ...
                                dln_node_coord(:, 1), dln_node_coord(:, 2),...
                                se_node_coord (:, 1), se_node_coord (:, 2),...
                                se_node_coord (:, 3), se_node_coord (:, 4),...
                                burgers(:), mu, nu, a, n_se, n_dln);
    end %if
    
    %% Map analytical nodal forces into a useful form for the force superposition scheme.
    % Allocate f_dln.
    % Find the number of nodes for which tractions need to be calculated.
    n_nodes_t = size(gamma, 1);
    % Allocate the force vector for forces induced by dislocations.
    f_dln = zeros(3*n_nodes_t, 1);
    % Loop through the number of nodes.
    for i = 1: n_nodes_t
        % Find the indices of the node labels corresponding to the node
        % gamma(i). Multiplied by three because each node has 3
        % coordinates. This lands idxi on the index which corresponds to 
        % the z-coordinate of node gamma(i).
        idxi = 3*find(se_node_label == gamma(i));
        % The FEM forces are given in the same node order as gamma, so this
        % finds the index of the nodal forces corresponding to gamma(i).
        % Multiplied by 3 because the nodes have 3 coordinates. This lands 
        % idxi on the index which corresponds to the z-coordinate of node 
        % gamma(i).
        idxf = i*3;
        % Loop through coordinates.
        for j = 2:-1:0
            % The index is displaced by -2, -1, 0, corresponding to the
            % indices of the x, y, z coordinate respectively.
            % We add the force contributions from all surface elements a node
            % is part of. This gives us the total x,y,z forces on each node.
            f_dln(idxf-j) = sum(nodal_force(idxi-j));
        end %for
    end %for
    
    %% Cleanup.
    clear n_nodes_t; clear idxi; clear idxf;
end % function