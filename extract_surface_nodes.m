function [labels, coords, n_se] = extract_surface_nodes(            ...
                                      fem_nodes     , fem_node_cnct,...
                                      fem_dim       , fem_planes   ,...
                                      n_nodes       , surf_node_util)
    %%===================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 11/12/2017--15/12/2017.
    % Refactored into independent subroutines in 11/12/2017--15/12/2017.
    %---------------------------------------------------------------------%
    % Couples analytical forces to the surfaces of a rectangular
    % cantilever whose nodes are labelled as follows:
    %
    % 4. ----- .3
    %  |\       |\
    %  | \      | \
    % 1. -\---- .2 \
    %   \  \     \  \
    %    \ 8. ----\- .7
    %     \ |      \ |
    %      \|       \|
    %      5. ----- .6
    % ^z
    % |  y
    % | /
    % |/
    % |------>x
    %
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % fem_nodes = dimension(:,3). Assumed shape 2D array with 3 columns.
    %   Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    %   columns. Describes the node connectivity of the finite elements in
    %   the model.
    %
    % dln_nodes = dimension(:,4). Assumed shape 2D array with 4 columns.
    %   Describes the nodes of the dislocation ensemble and mobility.
    %
    % dln_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    %   columns. Describes the node connectivity of the dislocation ensemble,
    %   associated Burgers' vectors and slip planes.
    %
    % fem_dim := [dim_x; dim_y; dim_z], finite element model dimensions.
    %   dim_x := dimension in x-direction, dim_y := dimension in y-direction,
    %   dim_z := dimension in z-direction.
    %
    % fem_planes := dimension(:). Assumed shape 1D array (maximum length 6)
    %   with the faces for which tractions are to be calculated. In order to
    %   make it consistent with the FEMcoupler and the analytical forces
    %   exerted by dislocations on free surfaces. The faces are by default
    %   defined by the areas made up of the following nodes in the diagram:
    %   min(x), yz-plane, nodes [5, 1, 8, 4], face 1
    %   max(x), yz-plane, nodes [2, 6, 3, 7], face 2
    %   min(y), xz-plane, nodes [6, 5, 7, 8], face 3
    %   max(y), xz-plane, nodes [1, 2, 4, 3], face 4
    %   min(z), xy-plane, nodes [5, 6, 1, 2], face 5
    %   max(z), xy-plane, nodes [4, 3, 8, 7], face 6.
    %   However, should the need arise, the node labels and face numbering
    %   can be changed to accomodate externally generated meshes.
    %   Alternatively this subroutine need not be used if only surface
    %   meshes are imported. However, the node ordering is crucial for the
    %   correct calculation of the problem as it must be consistent with
    %   S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).
    %
    % n_nodes := number of nodes per element.
    %
    %---------------------------------------------------------------------%
    %
    % Optional parameters:
    %
    % surf_node_util := dimension(n_nodes+2, 6). Each column is laid out as
    %   follows:
    %       [node_label_initial; ... ; node_label_final; plane_area; coordinate_number].
    %   Each row corresponds to a plane.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % util_flag := logical flag used to store whether surf_node_util was
    %   provided. If it wasn't, surf_node_util is deleted at the end of the
    %   subroutine to free memory.
    %
    % n_nodes_p1 := n_nodes + 1. Used when surf_node_util.
    %
    %=====================================================================%
    % Outputs
    %=====================================================================%
    %
    % labels := dimension(n_se, n_nodes). 2D array containing the global
    %   node labels of the relevant surface nodes.
    %
    % coords := dimension(3*n_se, n_nodes). 2D array containing the
    % cartesian coordinates of the relevant surface element nodes.
    %
    % n_se := number of surface elements in the calculation.
    %
    %%===================================================================%%

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
        surf_node_util(1:6, 1) = [1, 5, 4, 8, yz, 1]; % min(x), yz-plane, face 1 ~Sleft
        surf_node_util(1:6, 2) = [6, 2, 7, 3, yz, 1]; % max(x), yz-plane, face 2 ~Sright
        surf_node_util(1:6, 3) = [5, 6, 8, 7, xz, 2]; % min(y), xz-plane, face 3 ~Sfront
        surf_node_util(1:6, 4) = [2, 1, 3, 4, xz, 2]; % max(y), xz-plane, face 4 ~Sback
        surf_node_util(1:6, 5) = [6, 5, 2, 1, xy, 3]; % min(z), xy-plane, face 5 ~Sbot
        surf_node_util(1:6, 6) = [3, 4, 7, 8, xy, 3]; % max(z), xy-plane, face 6 ~Stop
    end %if

    %% Calculate number of surface elements.
    n_se = 0;
    n_nodes_p1 = n_nodes + 1;
    for i = 1: size(fem_planes, 1)
        n_se = n_se + surf_node_util(n_nodes_p1, fem_planes(i));
    end %for

    %% Extract relevant node planes in the correct order.
    [labels, coords] = extract_node_planes(                           ...
                            fem_nodes , fem_node_cnct, surf_node_util,...
                            fem_planes, n_se         , n_nodes);

    %% Cleanup.
    if util_flag
        clear surf_node_util;
    end %if
    clear util_flag; clear xy; clear xz; clear yz; clear n_nodes_p1;
end % function
