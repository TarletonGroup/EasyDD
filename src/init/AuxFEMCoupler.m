function [f_tilda, tract] = ...
    AuxFEMCoupler(FEM, gamma, tract, flags)
    %=========================================================================%
    % Sets up auxiliary data structures for analytic traction calculations.
    %
    % Daniel Celis Garza, Aug 2020
    % daniel.celisgarza@materials.ox.ac.uk
    %-------------------------------------------------------------------------%
    % Inputs
    % mno := total number of FE nodes.
    % dx, dy, dz := dimensions in x, y, z coordinates.
    % mx, my, mz := number of nodes in x, y, z dimension.
    % xnodes := coordinates and labels of FE nodes
    % nc := FE node connectivity matrix
    % gamma := traction, displacement, mixed boundary
    %   conditions.
    % a_trac := flag for analytic tractions.
    % CUDA_flag := flag in case CUDA codes required. If true compile, else do
    %   not compile.
    %-------------------------------------------------------------------------%
    % Local variables
    % x3x6_lbl := node labels of surface elements.
    %-------------------------------------------------------------------------%
    % Outputs
    % f_hat := tractions on FE nodes
    % para_tol := tolerance for calling segments parallel to the surface
    % x3x6 := coordinates of surface element nodes
    % n_se := number of surface elements
    % gamma_dln := nodes where analytic tractions need to be calculated
    % f_tilda_node := dislocation forces on nodes of surface elements
    % f_tilda_se := dislocation forces on single surface elements
    % f_tilda := dislocation forces on FE nodes
    % idxi := index for adding forces from nodes shared by different surface
    %   elements
    % n_nodes := number of nodes per surface element
    % n_nodes_t := number of nodes with traction boundary conditions
    % n_threads := number of threads per GPU block
    % para_scheme := parallelisation scheme, 1 parallelises over dislocations
    %   2 parallelises over surface elements
    % gamma_disp := nodes with displacement boundary conditions
    %=========================================================================%
    
    %% Extraction
    
    % flags:
    a_trac = flags.a_trac;
    CUDA_flag = flags.CUDA_flag;
    
    % tract:
    n_threads = tract.n_threads;
    para_scheme = tract.para_scheme;
    
    % FEM:
    ndofs = FEM.ndofs;
    mno = FEM.mno;
    dx = FEM.dx; dy = FEM.dy; dz = FEM.dz;
    mx = FEM.mx; my = FEM.my; mz = FEM.mz;
    n_nodes = FEM.n_nodes;
    xnodes = FEM.xnodes;
    nc = FEM.nc;
    
    %% Set up auxiliary data
    
    if (~exist('CUDA_flag', 'var'))
        CUDA_flag = false;
    end

    % Parallel CUDA C flags.
    if CUDA_flag == true
        % Provide a default number of threads in case none is given.
        if ~exist('n_threads', 'var')
            n_threads = 256;
        end %if

        % Provide a default parallelisaion scheme in case none is given.
        if ~exist('para_scheme', 'var')
            % Parallelise over dislocations.
            para_scheme = 1;
        end %if

    else
        n_threads = 0;
        para_scheme = 0;
    end %if

    if (~exist('a_trac', 'var'))
        a_trac = false;
    end
    
    if a_trac == false % Analytic tractions disabled
        
        para_tol = 0;
        x3x6 = 0;
        n_se = 0;
        f_tilda_node = 0;
        f_tilda_se = 0;
        idxi = 0;
        n_nodes_t = 0;
        f_tilda = zeros(ndofs, 1);
        
    else % Analytic tractions enabled
        
        dimension = sqrt(dx * dx + dy * dy + dz * dz);
        para_tol = dimension / 1e7;

        planes = (1:1:6)';
        [x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx; my; mz], ...
            planes, n_nodes);

        [f_tilda_node, f_tilda_se, ...
                f_tilda, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma.dln, n_nodes, n_se, mno);
    end
        
    %% Storage
    
    tract = struct; % Stores auxillary variables
    
    tract.para_tol = para_tol;
    tract.x3x6 = x3x6;
    tract.n_se = n_se;
    tract.idxi = idxi;
    tract.n_nodes = n_nodes;
    tract.n_nodes_t = n_nodes_t;
    tract.n_threads = n_threads;
    tract.para_scheme = para_scheme;
    tract.f_tilda_node = f_tilda_node;
    tract.f_tilda_se = f_tilda_se;
    
end
