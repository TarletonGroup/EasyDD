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
    dim_3 = dim*3;
    idxf = idxi + dim_3 - 1;
    tmp_nodes  = fem_nodes(fem_node_cnct(:, node_labels), 1:3);
    mask       = tmp_nodes(:, filter(1)) == filter(2);
    tmp_nodes  = reshape(tmp_nodes(mask, :)', n_nodes * dim_3, 1);
    step = 0;
    disp(idxi)
    for i = 1: n_nodes
        node_plane(idxi: idxf, i) = tmp_nodes(1: dim_3);
        step = step + dim_3;
    end %for
    idxi = idxf + 1;
    
    clear tmp_nodes; clear mask; clear n_nodes; clear step; clear idxf;
end %function