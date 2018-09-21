function [f_dln] = analytic_traction(...
                fem_nodes, fem_node_cnct, dln_nodes  , dln_node_cnct,...
                fem_dim  , fem_planes   , n_nodes    , gamma        ,...
                mu, nu, a, use_gpu      , para_scheme, surf_node_util)
    %%===================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 12/11/17--16/11/17.
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
    % Describes the coordinates of the nodes of the ft_inite element model.
    %
    % fem_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    % columns. Describes the node connectivity of the finite elements in
    % the model.
    %
    % dln_nodes = dimension(:,4). Assumed shape 2D array with 4 columns.
    % Describes the nodes of the dislocation ensemble and mobility.
    %
    % dln_node_cnct = dimension(:,8). Assumed shape 2D array with 8
    % columns. Describes the node connectivity of the dislocation ensemble,
    % associated Burgers' vectors and slip planes.
    %
    % fem_dim := [dim_x; dim_y; dim_z], finite element model dimensions.
    % dim_x := dimension in x-direction, dim_y := dimension in y-direction,
    % dim_z := dimension in z-direction.
    %
    % fem_planes := dimension(:). Assumed shape 1D array (maximum length 6)
    % with the faces for which tractions are to be calculated. In order to
    % make it consistent with the FEMcoupler and the analytical forces
    % exerted by dislocations on free surfaces. The faces re det_fined by the
    % areas made up of the following nodes in the diagram:
    % min(x), yz-plane, nodes [5, 1, 8, 4], face 1
    % max(x), yz-plane, nodes [2, 6, 3, 7], face 2
    % min(y), xz-plane, nodes [6, 5, 7, 8], face 3
    % max(y), xz-plane, nodes [1, 2, 4, 3], face 4
    % min(z), xy-plane, nodes [5, 6, 1, 2], face 5
    % max(z), xy-plane, nodes [4, 3, 8, 7], face 6.
    % The node ordering for each face ensures self-consistency in the force
    % calculation according to the labelling scheme in S. Queyreau,
    % J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).
    %
    % n_nodes := number of nodes per element.
    %
    % mu := shear modulus of material.
    %
    % nu := poisson's ratio.
    %
    % a := dislocation core size parameter (non-singular dislocations).
    %
    %---------------------------------------------------------------------%
    %
    % Flags:
    %
    % use_gpu := Flag to use the Graphics Processing Unit (GPU) in the
    % calculation. If (use_gpu == 1) {use CUDA C to run on an NVidia
    % GPU} else {use serial code in C}.
    % code.
    %
    %---------------------------------------------------------------------%
    %
    % Optional parameters:
    %
    % para_scheme := parallelisation scheme used when use_gpu == 1. If
    % (para_scheme == 1) {parallelise over dislocation network} elseif
    % (para_scheme == 2) {parallelise over surface elements}. Defaults to
    % 1.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % n_se := number of surface elements.
    %
    % n_dln := number of dislocation line segments.
    %
    % dln := dislocation line segments with Burgers' vector.
    %
    % surf_node_util := dimension(n_nodes+2, 6). Each column is laid out as
    % follows:
    %   [node_label_initial; ... ; node_label_final;
    %    plane_area; coordinate_number].
    % Each row corresponds to one plane.
    %
    % x1x2 := dimension(3*n_dln, 2). The first column has the
    % xyz-coordinates of the x1 nodes (dislocation line nodes) in the force
    % calculation, the second contains those of the x2 nodes.
    %
    % x3x6 := dimension(3*n_se, 4). Same thing as x1x2, only for the surface
    % elements.
    %
    %=====================================================================%
    % Outputs
    %=====================================================================%
    %
    % nodal_force := dimension(4*n_se, 3). 3D force on each
    % surface element node.
    %
    % total_force := dimension(n_se, 3). Total 3D force on each surface
    % element (sum over all nodes of an element).
    %
    %%===================================================================%%

    %% Generate dislocation line segments.
    dln   = constructsegmentlist(dln_nodes, dln_node_cnct)';
    n_dln = size(dln, 2);

    % Extract dislocation line nodes and burgers vectors.
    x1x2 = zeros(3 * n_dln, 2);
    x1x2(:, 1) = reshape(dln(6:8 , :), 3 * n_dln, 1);
    x1x2(:, 2) = reshape(dln(9:11, :), 3 * n_dln, 1);
    b  = reshape(dln(3:5,  :), 3 * n_dln, 1);
    clear dln;

    %% Generate surface element nodes.
    % Calculate face dimensions.
    xy = fem_dim(1)*fem_dim(2);
    xz = fem_dim(1)*fem_dim(3);
    yz = fem_dim(2)*fem_dim(3);

    % If the matrix is undefined then provide a default.
    if ~exist('surf_node_util','var')
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
    % Calculate number of surface elements.
    n_se = 0;
    for i = 1: size(fem_planes, 1)
        n_se = n_se + surf_node_util(5, fem_planes(i));
    end %for    

    [x3x6_lbl, x3x6] = extract_node_planes(fem_nodes, fem_node_cnct, surf_node_util, fem_planes, n_se, n_nodes);
    clear surf_node_util;

    %% Force calculation.
    % Allocating nodal and total force arrays
    nodal_force = zeros(3*n_se, n_nodes);
    total_force = zeros(3*n_se, 1);

    % CUDA C calculation
    if (use_gpu == 1)
        n_threads = 512;
        if (~exists(para_scheme))
            para_scheme = 1;
        end %if
        [nodal_force(:, 1), nodal_force(:, 2),...
         nodal_force(:, 3), nodal_force(:, 4),...
         total_force] = nodal_surface_force_linear_rectangle_mex(           ...
                                x1x2(:, 1), x1x2(:, 2),                        ...
                                x3x6(:, 1), x3x6(:, 2), x3x6(:, 3), x3x6(:, 4),...
                                b(:), mu, nu, a, n_se, n_dln,...
                                n_threads, para_scheme);
    % Serial force calculation in C
    else
        [nodal_force(:, 1), nodal_force(:, 2),...
         nodal_force(:, 3), nodal_force(:, 4),...
         total_force] = nodal_surface_force_linear_rectangle_arr(           ...
                                x1x2(:, 1), x1x2(:, 2),                        ...
                                x3x6(:, 1), x3x6(:, 2), x3x6(:, 3), x3x6(:, 4),...
                                b(:), mu, nu, a, n_se, n_dln);
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
        idxi = 3*find(x3x6_lbl == gamma(i));
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
end %function


%%=======================================================================%%
% This is faster if only extracting 1 to 3 planes. Replace the
% "%% Generate surface element nodes" block with this one if that is the case.
% This way should always be slower but MATLAB is dumb like that.
%%-----------------------------------------------------------------------%%
% % Set surface node labels for surface node extraction.
%     surf_node_util         = zeros(n_nodes+3, 6);
%     surf_node_util(1:6, 1) = [5, 1, 8, 4, yz, 1]; % min(x), yz-plane, face 1
%     surf_node_util(1:6, 2) = [2, 6, 3, 7, yz, 1]; % max(x), yz-plane, face 2
%     surf_node_util(1:6, 3) = [6, 5, 7, 8, xz, 2]; % min(y), xz-plane, face 4
%     surf_node_util(1:6, 4) = [1, 2, 4, 3, xz, 2]; % max(y), xz-plane, face 3
%     surf_node_util(1:6, 5) = [5, 6, 1, 2, xy, 3]; % min(z), xy-plane, face 5
%     surf_node_util(1:6, 6) = [4, 3, 8, 7, xy, 3]; % max(z), xy-plane, face 6
%
%     clear xy; clear xz; clear yz;
%     % Allocate x3 to x6.
%     x3x6 = zeros(3*n_se, 4);
%     % Extremal values for x, y, z.
%     for i = 1: 2: size(surf_node_util, 2)
%         surf_node_util(7, i  ) = min(fem_nodes(fem_node_cnct(:, surf_node_util(1:n_nodes, i  )), surf_node_util(6, i  )));
%         surf_node_util(7, i+1) = max(fem_nodes(fem_node_cnct(:, surf_node_util(1:n_nodes, i+1)), surf_node_util(6, i+1)));
%     end %for
%     % Indices for vectorised code.
%     idxi = 1;
%      for i = 1: size(fem_planes, 1)
%         plane_idx = fem_planes(i);
%         dim_3 = surf_node_util(n_nodes+1, plane_idx)*3;
%         idxf = idxi + dim_3 - 1;
%         x3x6(idxi: idxf, :) = extract_node_plane(fem_nodes, fem_node_cnct,...
%                             surf_node_util(1:n_nodes          , plane_idx)',...
%                             dim_3 ,...
%                             surf_node_util(n_nodes+2:n_nodes+3, plane_idx),...
%                             x3x6(idxi: idxf, :));
%         idxi = idxf + 1;
%      end %for
%
%     test = extract_node_planes(fem_nodes, fem_node_cnct, surf_node_util_test, fem_planes, n_se, n_nodes);
%     clear surf_node_util; clear plane_idx;
