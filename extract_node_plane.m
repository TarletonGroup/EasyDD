function [node_plane, idxi]= extract_node_plane(fem_nodes, fem_node_cnct, node_labels, dim, filter, idxi, node_plane)
    %%
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % fem_nodes = dimension(:,3). Assumed shape 2D array with 3 columns.
    % Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    % columns. Describes the node connectivity of the ft_inite elements in
    % the model.
    %%
    n_nodes     = size(node_labels, 2);
    if(n_nodes ~= size(node_plane , 2))
        error('extract_node_plane:: n_nodes = %d, size(node_plane, 2) = %d; they must be the same i.e. the number of nodes per planar element', n_nodes, size(node_plane,2))
    end %if
    idxf = idxi + 3*dim - 1;
    tmp_nodes  = fem_nodes(fem_node_cnct(:, node_labels), 1:3);    
    mask       = fem_nodes(fem_node_cnct(:, node_labels), filter(1)) == filter(2);
    tmp_nodes  = reshape(tmp_nodes(mask, :)', 3 * n_nodes * dim, 1);
    
    step = 0;
    for i = 1: n_nodes
        node_plane(idxi: idxf, i) = tmp_nodes(idxi + step: idxf + step);
        step = step + 3*dim;
    end %for
    idxi = idxf + 1;
    
    clear tmp_nodes; clear mask; clear n_nodes; clear step; clear idxf;
end %function




function node_plane = extract_node_plane(fem_nodes, fem_node_cnct, node_labels, fem_faces, n_plane_elem, n_nodes)
    %%
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % fem_nodes = dimension(:,3). Assumed shape 2D array with 3 columns.
    % Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    % columns. Describes the node connectivity of the ft_inite elements in
    % the model.
    %%
    node_plane = zeros(3*n_plane_elem, n_nodes);
    idxi = 1;
    for i = 1: size(fem_faces, 1)
        face_idx = fem_faces(i);
        idxf = idxi + 3*dim - 1;
        tmp_nodes  = fem_nodes(fem_node_cnct(:, node_labels(1:n_nodes, face_idx), 1:3));
        % If face is odd it is a minimum value face.
        if(mod(face_idx, 2) ~= 0)
            filter = min(tmp_nodes(:, node_labels(6, face_idx)));
        % If face is even it is a maximum value face.
        else
            filter = max(tmp_nodes(:, node_labels(6, face_idx)));
        end %if
        mask       = fem_nodes(fem_node_cnct(:, node_labels(1:n_nodes, face_idx)), node_labels(6, face_idx)) == filter;
        tmp_nodes  = reshape(tmp_nodes(mask, :)', 3 * n_nodes * dim, 1);

        step = 0;
        for i = 1: n_nodes
            node_plane(idxi: idxf, i) = tmp_nodes(idxi + step: idxf + step);
            step = step + 3*dim;
        end %for
        idxi = idxf + 1;
    end %for
    
    clear tmp_nodes; clear mask; clear n_nodes; clear step; clear idxf;
end %function