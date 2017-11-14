function x3x6 = analytic_traction(                   ...
                                        fem_nodes, fem_node_cnct, ...
                                        dln_nodes, dln_node_cnct, ...
                                        fem_dim  , fem_faces    , ...
                                        mu,    nu,     a        , ...
                                        use_gpu  , para_scheme    ...
                                       )
    %[nodal_force, total_force]
    %%===================================================================%%
    % Couples analytical forces to the surfaces of a rectangular
    % cantilever.
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
    % columns. Describes the node connectivity of the ft_inite elements in
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
    % fem_faces := dimension(:). Assumed shape 1D array (maximum length 6) 
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
    % surf_nodes := dimension(7, 6). Each column is laid out as follows:
    % [node_label_1; node_label_2; node_label_3; node_label_4; 
    %  plane_dimension; coordinate_number; plane_selection_criteria]
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
    
    % Generate dislocation line segments.
    dln   = constructsegmentlist(dln_nodes, dln_node_cnct)';
    n_dln = size(dln, 2);
    
    % Generate dislocation line nodes and burgers vectors.
    x1 = reshape(dln(6:8 , :), 3 * n_dln, 1);
    x2 = reshape(dln(9:11, :), 3 * n_dln, 1);
    b  = reshape(dln(3:5,  :), 3 * n_dln, 1);
    clear dln;
    
    % Generate surface element nodes.
    xy = fem_dim(1)*fem_dim(2);
    xz = fem_dim(1)*fem_dim(3);
    yz = fem_dim(2)*fem_dim(3);
    
    % Calculate the number of surface elements.
    n_se = 0;
    for i = 1: size(fem_faces, 1)
        fem_face = fem_faces(i);
        % If we have an xz face we add mx*mz elements.
        if(fem_face == 1 || fem_face == 2)
            n_se = n_se + xz;
        % If we have a yz face we add my*mz elements. 
        elseif(fem_face == 3 || fem_face == 4)
            n_se = n_se + yz;
   %     min(fem_node_cnct(:, surf_nodes(1:4, 1)),...
%                                             surf_nodes(6  , 1))
     % Otherwise we have an xy face and we add mx*my elements.
        else
            n_se = n_se + xy;
        end %if
    end %for
    % Allocate x3 to x6.
    x3x6 = zeros(3*n_se, 4);
    
    % Set surface node labels for surface node extraction.
    surf_nodes         = zeros(7, 6);
    surf_nodes(1:6, 1) = [5, 1, 8, 4, yz, 1]; % min(x), yz-plane, face 1
    surf_nodes(1:6, 2) = [2, 6, 3, 7, yz, 1]; % max(x), yz-plane, face 2
    surf_nodes(1:6, 3) = [6, 5, 7, 8, xz, 2]; % min(y), xz-plane, face 4
    surf_nodes(1:6, 4) = [1, 2, 4, 3, xz, 2]; % max(y), xz-plane, face 3
    surf_nodes(1:6, 5) = [5, 6, 1, 2, xy, 3]; % min(z), xy-plane, face 5
    surf_nodes(1:6, 6) = [4, 3, 8, 7, xy, 3]; % max(z), xy-plane, face 6
    
    clear xy; clear xz; clear yz;
    % Extremal values for x, y, z.
    for i = 1: 2: size(surf_nodes, 2)
        surf_nodes(7, i  ) = min(fem_nodes(fem_node_cnct(:, surf_nodes(1:4, i  )), surf_nodes(6, i  )));
        surf_nodes(7, i+1) = max(fem_nodes(fem_node_cnct(:, surf_nodes(1:4, i+1)), surf_nodes(6, i+1)));
    end %for
      
    % Indices for vectorised code.
    idxi = 1;
     for i = 1: size(fem_faces, 1)
        face_idx = fem_faces(i);
        [x3x6, idxi] = extract_node_plane(fem_nodes, fem_node_cnct,...
                            surf_nodes(1:4, face_idx)',...
                            surf_nodes(5  , face_idx) ,...
                            surf_nodes(6:7, face_idx) , idxi, x3x6);
     end %for
    clear surf_nodes; clear face_idx;
    

end %function
