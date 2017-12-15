function [labels, coords, n_se] = extract_surface_nodes(...
                                      fem_nodes     , fem_node_cnct,...
                                      fem_dim       , fem_planes   ,...
                                      n_nodes       , surf_node_util)
    %% Calculate face dimensions.
    xy = fem_dim(1)*fem_dim(2);
    xz = fem_dim(1)*fem_dim(3);
    yz = fem_dim(2)*fem_dim(3);

    %% Node selection utility matrix.
    % If the matrix is undefined then provide a default for a cubic prism.
    util_flag = ~true(1);
    if ~exist('surf_node_util','var')
        util_flag = true(1);
        % Set surface node labels for surface node extraction.
        surf_node_util = zeros(n_nodes+2, 6);
        % For rectangular surface elements. Add an if statement and redifine
        % matrix to generalise it for triangular surface elements.
        surf_node_util(1:6, 1) = [5, 1, 8, 4, yz, 1]; % min(x), yz-plane, face 1
        surf_node_util(1:6, 2) = [2, 6, 3, 7, yz, 1]; % max(x), yz-plane, face 2
        surf_node_util(1:6, 3) = [6, 5, 7, 8, xz, 2]; % min(y), xz-plane, face 4
        surf_node_util(1:6, 4) = [1, 2, 4, 3, xz, 2]; % max(y), xz-plane, face 3
        surf_node_util(1:6, 5) = [5, 6, 1, 2, xy, 3]; % min(z), xy-plane, face 5
        surf_node_util(1:6, 6) = [4, 3, 8, 7, xy, 3]; % max(z), xy-plane, face 6
    end %if
    
    %% Calculate number of surface elements.
    n_se = 0;
    n_nodes_p1 = n_nodes + 1;
    for i = 1: size(fem_planes, 1)
        n_se = n_se + surf_node_util(n_nodes_p1, fem_planes(i));
    end %for    
    
    %% Extract relevant node planes in the correct order.
    [labels, coords] = extract_node_planes(fem_nodes, fem_node_cnct, surf_node_util, fem_planes, n_se, n_nodes);
    
    %% Cleanup.
    if util_flag
        clear surf_node_util;
    end %if
    clear util_flag; clear xy; clear xz; clear yz; clear n_nodes_p1;
end % function