function node_plane = extract_node_plane(fem_nodes  , fem_node_cnct,...
                                         node_labels, dim, filter  ,...
                                         node_plane)
    %%===================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 13/11/17. 
    %
    % This code has been optimised as much as I can and has been thoroughly
    % tested as of 15/11/17.
    %---------------------------------------------------------------------%
    % 
    % Extracts a single node plane according to given filtering criteria
    % and puts it into a 2D array arranged in xyz-coordinates like so.
    %
    %     Node 1    Node 2    Node 3 ... Node N
    %
    %      x_11      x_12      x_13       x_1N   |
    %      y_11      y_12      y_13       y_1N   | Element 1
    %      z_11      z_12      z_13  ...  z_1N   |___________
    %      x_21      x_22      x_23       x_2N   |
    %      y_21      y_22      y_23       y_2N   | Element 2
    %      z_21      z_22      z_23       z_2N   |___________
    %             .                                    .
    %                       .                          .
    %                                 .                .
    %                                             ___________
    %      x_E1      x_E2      x_E3       x_EN   |
    %      y_E1      y_E2      y_E3  ...  y_EN   | Element E
    %      z_E1      z_E2      z_E3       z_EN   |
    %
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % fem_nodes := dimension(:, 3). Assumed shape 2D array with 3 columns.
    % Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct := dimension(:, 8). Assumed shape 2D array with 8
    % columns. Describes the node connectivity of the finite elements in
    % the model.
    %
    % node_labels := dimension(:). The labels of the nodes that are to be
    % extracted from fem_nodes.
    %
    % dim := three times the area of the plane to be extracted (3
    % dimensions).
    %
    % filter := the selection criteria used to filter which nodes we want
    % to obtain.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % n_nodes := number of nodes per planar element.
    %
    % tmp_nodes := dimension(n_nodes, 3). Contains the xyz-coordinates of 
    % all nodes with the labels specified by node_label.
    %
    % mask := same dimension as tmp_nodes. Logical mask whose entries are 1
    % only when a given condition is met. It is used to index tmp_nodes to 
    % find only the nodes that correspond to a given plane.
    %
    %=====================================================================%
    % Input/output variables
    %=====================================================================%
    %
    % node_plane := dimension(dim/3, n_nodes). It is the slice of the array
    % that holds the nodes of the plane. Only pass the array slice
    % that corresponds to the current plane of nodes to be extracted,
    % a la Fortran intent(in out).
    %
    %%===================================================================%%
    
    %% Error checking
    n_nodes     = size(node_labels, 2);
    if(n_nodes ~= size(node_plane , 2))
        error('extract_node_plane:: n_nodes = %d, size(node_plane, 2) = %d; they must be the same i.e. the number of nodes per planar element', n_nodes, size(node_plane,2))
    end %if
    
    %% Finding plane nodes
    % Find the all the nodes labelled as stated in node_labels 
    % (corresponding to a plane family)
    tmp_nodes = fem_nodes(fem_node_cnct(:, node_labels), 1:3);
    % Find an array of entries that meet the filter criteria in tmp_nodes
    mask = tmp_nodes(:, filter(1)) == filter(2);
    % Reshape the array to a column array containing only the nodes which
    % meet the filter criteria.
    tmp_nodes = reshape(tmp_nodes(mask, :)', n_nodes * dim, 1);
    
    %% Passing data to the slice of the array that contains the nodes for every plane
    step = 0;
    for i = 1: n_nodes
        node_plane(1: dim, i) = tmp_nodes(1 + step: dim + step);
        step = step + dim_3;
    end %for
    
    %% Cleanup
    clear dim_3; clear idxf; clear tmp_nodes; clear mask; clear step;
end %function