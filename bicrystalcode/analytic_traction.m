function [f_dln, f_dln_se, f_dln_node] = analytic_traction(                                    ...
                     se_node_coord, dln_node_coord, burgers    , n_nodes,...
                     n_nodes_t    , n_se          , n_dln      , idxf   ,...
                     idxi         , f_dln_node    , f_dln_se   , f_dln  ,...
                     mu, nu, a    , use_gpu       , n_threads  , para_scheme, ...
                     para_tol)
    %%===================================================================%%
    %---------------------------------------------------------------------%
    % Written by famed MATLAB hater and fan of compiled languages,
    % Daniel Celis Garza, in 12/11/17--16/11/17.
    % Refactored into independent subroutines in 11/12/2017--15/12/2017.
    %---------------------------------------------------------------------%
    % Calculates analytical forces to the surfaces of a rectangular
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
    % Everything has been arranged to conform with the analytical solutions
    % for forces exerted on linear rectangular surface elements published by:
    % S. Queyreau, J. Marian, B.D. Wirth, A. Arsenlis, MSMSE, 22(3):035004, (2014).
    %
    %=====================================================================%
    % Inputs
    %=====================================================================%
    %
    % se_node_coord := dimension(n_se*n_nodes, 3). 2D array
    %   with the cartesian coordinates of the finite element nodes subject
    %   to traction boundary conditions. The analytical solution requires
    %   they be arranged element by element in a specific order.
    %
    % dln_node_coord := dimension(n_dln*2, 3). 2D array with
    %   the cartesian coordinates of the dislocation line nodes. The
    %   analytical solution requires they be arranged line segment by
    %   line segment.
    %
    % burgers := dimension(n_dln, 3). 2D array with the
    %   individual dislocation line segments' Burgers vector. The
    %   analytical solution needs the Burgers vector for each line segment.
    %
    % n_nodes := number of nodes per surface element.
    %
    % n_nodes_t := total number of nodes with traction boundary conditions.
    %   It is needed to map the analytical solutions' forces to the force
    %   array used to account for the forces exerted by the dislocation
    %   ensemble on the finite element nodes. This does not necessarily
    %   contain a factor of n_nodes less entries than se_node_coord because
    %   not all nodes in an element may have traction boundary conditions,
    %   but still need to be included in the calculation due to the
    %   problem's formulation.
    %
    % n_se := number of surface elements
    %
    % n_dln := number of dislocation line segments
    %
    % idxf := dimension(:). 1D array containing the indices of the force
    %   array used to account for the forces exerted by the dislocation
    %   ensemble on the finite element nodes. It is needed to map the
    %   analytical solutions' forces to the right index. It is usually
    %   equal to gamma*3 because only the nodes with traction boundary
    %   conditions are needed.
    %
    % idxi := dimension(n_nodes_t*n_nodes). 1D assumed shape array
    %   containing the global label indices of the nodes whose tractions
    %   were calculated by the analytical solution.
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
    %   calculation. If use_gpu == 1, the code is run by a CUDA-enabled
    %   NVidia GPU. The code will crash if there is none, or if there is no
    %   compiled MEX CUDA file. If use_gpu == 0, a serial C code will be
    %   executed. The code will crash if there is no compiled MEX C file.
    %   If use_gpu is any other number, a MATLAB version will be used. This
    %   is much slower than either of the others.
    %
    %---------------------------------------------------------------------%
    %
    % Optional parameters:
    %
    % n_threads := when use_gpu ==1, it is the number of threads per block
    %   to be used in the parallelisation. This may require some
    %   experimentation to optimise. Defaults to 512. Does nothing if a
    %   GPU isn't used.
    %
    % para_scheme := when use_gpu is enabled, it dicates the
    %   parallelisation scheme used. If para_scheme == 1 parallelise over
    %   dislocation line segments. If para_scheme == 0 parallelise over
    %   surface elements. Other parallelisation schemes can be added in the
    %   CUDA C code.
    %
    %=====================================================================%
    % Dummy variables
    %=====================================================================%
    %
    % f_dln_node := dimension(3*n_se, 4). 3D force on each
    %   surface element node.
    %
    % f_dln_se := dimension(3*n_se, 1). Total 3D force on each surface
    %   element (sum over all nodes of an element). This isn't used to
    %   couple the forces to the finite element forces but can be used to
    %   provide the evolution of the forces exerted by the dislocations on
    %   the surface elements.
    %
    % tmp := dimension(n_nodes). 1D array containing the at-most 4
    %   instances where a surface element node appears (a single node can
    %   be shared by up to 4 surface elements).
    %
    % tmp2 := dimension(:). Assumed shape 1D array containing the at-most 4
    %   instances where a surface element node appears (a single node can
    %   be shared by up to 4 surface elements) but without any zero-valued
    %   entries because 0 indices crash MATLAB.
    %
    % k := used to traverse idxi to the next block of 4 entries which
    %   contain data relevant to the node in the
    %
    %=====================================================================%
    % Inputs/Outputs
    %=====================================================================%
    %
    % f_dln := dimension(n_nodes_t*3). 1D array containing the forces
    %   exerted by the dislocations on the surface nodes. It is used to
    %   correct the finite element forces.
    %
    %%===================================================================%%

    %% Analytical force calculation.
    % Parallel CUDA C calculation.
    % serial last dislocation and surface element
    if use_gpu == 1
        % Provide a default number of threads in case none is given.
        if ~exist('n_threads', 'var')
            n_threads = ceil(mod(n_dln,256)/32)*32;
        end %if
        % Provide a default parallelisaion scheme in case none is given.
        if ~exist('para_scheme', 'var')
            % Parallelise over dislocations.
            para_scheme = 1;
        end %if

        [f_dln_node(:, 1), f_dln_node(:, 2), f_dln_node(:, 3), f_dln_node(:, 4),...
         f_dln_se] = nodal_surface_force_linear_rectangle_mex_cuda(        ...
                                dln_node_coord(:, 1), dln_node_coord(:, 2)   ,...
                                se_node_coord (:, 1), se_node_coord (:, 2)   ,...
                                se_node_coord (:, 3), se_node_coord (:, 4)   ,...
                                burgers(:), mu, nu, a, n_se, n_dln, n_threads,...
                                para_scheme, para_tol);
    % Serial force calculation in C.
    elseif use_gpu == 0
        [f_dln_node(:, 1), f_dln_node(:, 2), f_dln_node(:, 3), f_dln_node(:, 4), ...
         f_dln_se] = nodal_surface_force_linear_rectangle_mex(          ...
                                dln_node_coord(:, 1), dln_node_coord(:, 2),...
                                se_node_coord (:, 1), se_node_coord (:, 2),...
                                se_node_coord (:, 3), se_node_coord (:, 4),...
                                burgers(:), mu, nu, a, n_se, n_dln, para_tol);
    end %if

%     f_dln_node(:, 1) = f_dln_se*0.25;
%     f_dln_node(:, 2) = f_dln_se*0.25;
%     f_dln_node(:, 3) = f_dln_se*0.25;
%     f_dln_node(:, 4) = f_dln_se*0.25;

    %% Map analytical nodal forces into a useful form for the force superposition scheme.
    % Loop through the number of nodes.
    k = 0;
    tmp = zeros(n_nodes, 1);
    for i = 1: n_nodes_t
        % Populate tmp array with the indices corresponding to nodes of
        % the surface relevant surface element.
        tmp  = idxi(1 + k: k + n_nodes);
        % Obtain only the indices which are non-zero. Zero indices mean
        % those nodes are not under traction boundary conditions.
        tmp2 = tmp(idxi(1+k: k + n_nodes, 1) ~= 0);
        for j = 2:-1:0
            % The index is displaced by -2, -1, 0, corresponding to the
            % indices of the x, y, z coordinate respectively.
            % We add the force contributions from all surface elements a node
            % is part of. This gives us the total x,y,z forces on each node.
            f_dln(idxf(i) - j) = sum(f_dln_node(tmp2 - j));
        end %for
        % Step to the block in idxi where the next surface node's indices
        % are found.
        k = k + 4;
    end %for
    clear tmp; clear tmp2; clear k;

end % function
