function [node_plane_lbl, node_plane] = extract_node_planes(            ...
                                          fem_nodes     , fem_node_cnct,...
                                          surf_node_util, fem_planes   ,...
                                          n_elem        , n_nodes)
    %%====================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 13/11/17.
    %
    % This code has been optimised as much as I can and has been thoroughly
    % tested as of 15/11/17.
    %---------------------------------------------------------------------%
    %
    % _lbl := variables suffixed with it are the global label equivalent
    %   of those without.
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
    %   Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct := dimension(:, 8). Assumed shape 2D array with 8
    %   columns. Describes the node connectivity of the finite elements in
    %   the model.
    %
    % surf_node_util := dimension(n_nodes+2, 6). Each column is laid out as
    %   follows:
    %       [node_label_initial; ... ; node_label_final;
    %       plane_area; coordinate_number].
    %   Each row corresponds to one plane.
    %
    % fem_planes := dimension(:). Contains the planes whose nodes are to be
    %   extracted.
    %
    % n_elem := total number of elements to be extracted.
    %
    % n_nodes := number of nodes per planar element.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % plane_idx := the index of the fem plane to be extracted during the
    %   current iteration.
    %
    % dim := plane_area, there are this many surface elements per plane.
    %
    % dim_3 := 3*plane_area, there are three coordinates per plane
    %   node, plane_area is the total number of nodes in a plane and is found
    %   in surf_node_utils.
    %
    % tmp_nodes := dimension(n_nodes, 3). Contains the xyz-coordinates of
    %   all nodes with the labels specified by node_label.
    %
    % mask := same dimension as tmp_nodes. Logical mask whose entries are 1
    %   only when a given condition is met. It is used to index tmp_nodes to
    %   find only the nodes that correspond to a given plane.
    %
    % coord := orthogonal coordinate to the plane (for node filtering
    %   purposes).
    %
    % idxi, idxf := initial and final index that slice the coordinate
    %   output array into sections for each extracted plane.
    %
    % idxl, idxm := initial and final index that slice the label
    %   output array into sections for each extracted plane.
    %
    % filter := value for filtering data.
    %
    % step := step that ensures proper array indexing for each section of
    %   the linear array that contains the nodes of a given plane.
    %
    % n_nodes_p1 := n_nodes + 1
    %
    % n_nodes_p2 := n_nodes + 2
    %
    %=====================================================================%
    % Output variables
    %=====================================================================%
    %
    % node_plane := dimension(3*n_elem, n_nodes). It is the array that
    %   holds the planes' nodes' Cartesian coordinates.
    %
    %%===================================================================%%

    %% Allocating node plane coordinates and labels
    node_plane     = zeros(3 * n_elem, n_nodes);
    node_plane_lbl = zeros(    n_elem, n_nodes);

    %% Looping through the planes to be extracted
    idxl = 1;
    idxi = 1;
    n_nodes_p_1 = n_nodes + 1;
    n_nodes_p_2 = n_nodes + 2;
    for i = 1: size(fem_planes, 1)
        %% Assigning control variables for the iteration
        plane_idx = fem_planes(i);
        dim   = surf_node_util(n_nodes_p_1, plane_idx);
        dim_3 = 3 * dim;
        coord = surf_node_util(n_nodes_p_2, plane_idx);
        % Set final index of the output array slice containing the nodes
        % for this iteration's plane
        idxm = idxl + dim   - 1;
        idxf = idxi + dim_3 - 1;
        %% Finding global node labels and node coordinates of the plane
        % Extracting all nodes whose labels correspond to the plane family
        % to be extracted this iteration
        tmp_nodes_lbl = fem_node_cnct(:, surf_node_util(1:n_nodes, plane_idx));
        tmp_nodes     = fem_nodes(tmp_nodes_lbl, 1:3);
        % If the plane index is odd it is a face where its orthogonal
        % coordinate is at a minimum.
        if(mod(plane_idx, 2) ~= 0)
            filter = min(tmp_nodes(:, coord));
        % Otherwise find the face whose orthogonal coordinate is at a maximum
        else
            filter = max(tmp_nodes(:, coord));
        end %if
        % Find an array of entries that meet the filter criteria in tmp_nodes
        mask = tmp_nodes(:, coord) == filter;
        % Reshape the array to a column array containing only the nodes which
        % meet the filter criteria.
        tmp_nodes = reshape(tmp_nodes(mask,:)', n_nodes * dim_3, 1);
        % Reshape the array with global node labels so it can be filtered
        % using the same mask.
        tmp_nodes_lbl_size = size(tmp_nodes_lbl, 1);
        tmp_nodes_lbl      = reshape(tmp_nodes_lbl ,...
                                     tmp_nodes_lbl_size * n_nodes, 1);
        tmp_nodes_lbl = tmp_nodes_lbl(mask);

        %% Passing data to the slice of the array that contains the nodes labels and coords for every plane
        % Pass the node labels of this iteration to the global array.
        % Set the slice of the tmp_nodes array to be indexed
        step     = 0;
        step_lbl = 0;
        for j = 1: n_nodes
            % Grab the array slice that corresponds to node j in the
            % tmp_nodes array and send it off to the corresponding column
            % of the slice of the output array containing all extracted nodes
            node_plane(idxi: idxf, j) = tmp_nodes(...
                                              1 + step: dim_3 + step...
                                            );
            node_plane_lbl(idxl: idxm, j) = tmp_nodes_lbl(...
                                              1 + step_lbl: dim + step_lbl...
                                            );
            % Advance the step to move onto the array slice that contains
            % node j+1 in tmp_nodes
            step     = step     + dim_3;
            step_lbl = step_lbl + dim  ;
        end %for
        % Set the initial index to shift the output array slice onto the
        % next plane to be extracted.
        idxi = idxf + 1;
        idxl = idxm + 1;
        clear tmp_nodes; clear tmp_nodes_lbl; clear mask;
    end %for

    %% Cleanup
    clear idxl; clear idxi; clear plane_idx; clear dim; clear dim_3;
    clear coord; clear idxm; clear idxf; clear filter; clear tmp_nodes_lbl_size;
    clear step; clear step_lbl; clear n_nodes_p1; clear n_nodes_p2;
end %function
