function [f_hat, para_tol, x3x6, n_se, gamma_dln, f_dln_node, f_dln_se,...
    f_dln, idxi, n_nodes_t, n_threads, para_scheme, gamma_disp] = ...
    analyticTractionAux(mno, dx, dy, dz, mx, my, mz, xnodes, nc, gammat,...
    gammau, gammaMixed, a_trac, CUDA_flag, n_threads, para_scheme)
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
% gammat, gammau, gammaMixed := traction, displacement, mixed boundary
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
% gamma_dln := nodes where analytic tractions need to be calcualted
% f_dln_node := dislocation forces on nodes of surface elements
% f_dln_se := dislocation forces on single surface elements
% f_dln := dislocation forces on FE nodes
% idxi := index for adding forces from nodes shared by different surface 
%   elements 
% n_nodes_t := number of nodes with traction boundary conditions
% n_threads := number of threads per GPU block
% para_scheme := parallelisation scheme, 1 parallelises over dislocations
%   2 parallelises over surface elements
% gamma_disp := nodes with displacement boundary conditions
%=========================================================================%

f_hat = zeros(3*mno, 1);

if(~exist('a_trac','var'))
    a_trac = false;
end
if a_trac == false
    return
end

if (~exist('CUDA_flag', 'var'))
    CUDA_flag = false;
end


dimension=sqrt(dx*dx+dy*dy+dz*dz);
para_tol = dimension/1e7;

planes = (1:1:6)';
[x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
    planes, 4);
gamma_dln = [gammat(:,1); gammaMixed(:,1)];
[f_dln_node, f_dln_se,...
    f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);

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

gamma_disp = [gammau; gammaMixed];
end