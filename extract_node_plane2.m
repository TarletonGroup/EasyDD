function node_plane = extract_node_plane2(fem_nodes, fem_node_cnct, node_labels, fem_faces, n_elem, n_nodes_elem)
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
    node_plane = zeros(3*n_elem, n_nodes_elem);
    idxi = 1;
    for i = 1: size(fem_faces, 1)
        face_idx = fem_faces(i);
        dim_3 = 3 * node_labels(5, face_idx);
        coord = node_labels(6, face_idx);
        idxf  = idxi + dim_3 - 1;
        tmp_nodes = fem_nodes(fem_node_cnct(:, node_labels(1:n_nodes_elem, face_idx)), 1:3);  
        % If face is odd it is a minimum value face.
        if(mod(face_idx, 2) ~= 0)
            filter = min(tmp_nodes(:, coord));
        % If face is even it is a maximum value face.
        else
            filter = max(tmp_nodes(:, coord));
        end %if
        mask       = tmp_nodes(:, coord) == filter;
        tmp_nodes  = reshape(tmp_nodes(mask,:)', n_nodes_elem * dim_3, 1);

        step = 0;
        for j = 1: n_nodes_elem
            node_plane(idxi: idxf, j) = tmp_nodes(1: dim_3);
            step = step + dim_3;
        end %for
        idxi = idxf + 1;
    end %for
    
    clear tmp_nodes; clear mask; clear n_nodes; clear step; clear idxf;
end %function