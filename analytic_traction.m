function [nodal_force, total_force] = traction(                   ...
                                        fem_nodes, fem_node_cnct, ...
                                        dln_nodes, dln_node_cnct, ...
                                        fem_dim  , fem_faces    , ...
                                        mu,    nu,     a        , ...
                                        use_gpu  , para_scheme    ...
                                       )
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
    % fem_dim := [mx my mz], ft_inite element model dimensions. 
    % mx := dimension in x-direction, my := dimension in y-direction, 
    % mz := dimension in z-direction.
    %
    % fem_faces := dimension(:). Assumed shape 1D array (maximum length 6) 
    % with the faces for which tractions are to be calculated. In order to 
    % make it consistent with the FEMcoupler and the analytical forces 
    % exerted by dislocations on free surfaces. The faces re det_fined by the
    % areas made up of the following nodes in the diagram:
    % min(y), xz-plane, nodes [6, 5, 7, 8], face 1
	% max(y), xz-plane, nodes [1, 2, 4, 3], face 2
	% min(x), yz-plane, nodes [5, 1, 8, 4], face 3
	% max(x), yz-plane, nodes [2, 6, 3, 7], face 4
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
    n_dln = size(dln,1);
    % Generate dislocation line nodes and burgers vectors.
    x1 = reshape(dln(6:8 , :), 3 * n_dln, 1);
    x2 = reshape(dln(9:11, :), 3 * n_dln, 1);
    b  = reshape(dln(3:5 , :), 3 * n_dln, 1);
    clear dln;
   
    % Generate surface element nodes.
    mx_mz = fem_dim(1)*fem_dim(3);
    my_mz = fem_dim(2)*fem_dim(3);
    mx_my = fem_dim(1)*fem_dim(2);
    % Calculate the number of surface elements.
    n_se = 0;
    for i = 1: size(fem_faces)
        fem_face = fem_faces(i);
        % If we have an xz face we add mx*mz elements.
        if(fem_face == 1 || fem_face == 2)
            n_se = n_se + mx_mz;
        % If we have a yz face we add my*mz elements. 
        elseif(fem_face == 3 || fem_face == 4)
            n_se = n_se + my_mz;
        % Otherwise we have an xy face and we add mx*my elements.
        else
            n_se = n_se + mx_my;
        end %if
    end %for
    
    % Allocate x3 to x6.
    x3   = zeros( 3 * n_se, 1);
    x4   = zeros( 3 * n_se, 1);
    x5   = zeros( 3 * n_se, 1);
    x6   = zeros( 3 * n_se, 1);
    x3x6 = zeros(12 * n_se, 1);
    % Indices for vectorised code.
    t_init = 0;
    t_fin  = 0;
    n_init = 0;
    n_fin  = 0;
    % Extremal values for x, y, z.
    max_x = max(fem_node_cnct(:, [1, 2, 4, 3]), 1);
    min_x = min(fem_node_cnct(:, [5, 1, 8, 4]), 1);
    
    max_y = max(fem_node_cnct(:, [2, 6, 3, 7]), 2);
    min_y = min(fem_node_cnct(:, [6, 5, 7, 8]), 2);
    
    max_z = max(fem_node_cnct(:, [4, 3, 8, 7]), 3);
    min_z = min(fem_node_cnct(:, [5, 6, 1, 2]), 3);
    
    for i = 1: size(fem_faces)
        fem_face = fem_faces(i);
        % If we have an xz face at minimum y.
        if(fem_face == 1)
            % Increase t_initial index to start after the end of the previous block.
            t_init = t_fin;
            t_fin  = t_fin + 12 * mx_mz;
            n_init = n_fin;
            n_fin  = n_fin + 3 * mx_mz;
            % This can be made into a function, x3, x4, x5, x6 can be made
            % into a 2D array of dimension(3*n_se, 4) and we can loop
            % through rows while vectorising columns. Set the dimensions of
            % x3x6, x3, x4, x5, x6 outside the function and assign stuff in
            % the function.
            nodes  = fem_nodes(fem_node_cnct(:, [6, 5, 7, 8]), 1:3);
            mask   = fem_nodes(fem_node_cnct(:, [6, 5, 7, 8]), 2) == min_y;
            x3x6(1 + t_init: t_fin) = reshape(nodes(mask, :)', 12 * mx_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*mx_mz: n_fin + 3*mx_mz);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*mx_mz: n_fin + 6*mx_mz);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*mx_mz: n_fin + 9*mx_mz);
        elseif(fem_face == 2)
            t_init = t_fin;
            t_fin  = t_fin + 12 * mx_mz;
            n_init = n_fin;
            n_fin  = n_fin + 3 * mx_mz;
            nodes  = fem_nodes(fem_node_cnct(:, [2, 6, 3, 7]), 1:3);
            mask   = fem_nodes(fem_node_cnct(:, [2, 6, 3, 7]), 2) == max_y;
            x3x6(1+t_init:t_fin) = reshape(nodes(mask, :)', 12 * mx_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*mx_mz: n_fin + 3*mx_mz);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*mx_mz: n_fin + 6*mx_mz);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*mx_mz: n_fin + 9*mx_mz);
        elseif(fem_face == 3)
            t_init = t_fin;
            t_fin  = t_fin + 12 * my_mz;
            n_init = n_fin;
            n_fin  = n_fin + 3 * my_mz;
            nodes  = fem_nodes(fem_node_cnct(:, [5, 1, 8, 4]), 1:3);
            mask   = fem_nodes(fem_node_cnct(:, [5, 1, 8, 4]), 1) == min_x;
            x3x6(1+t_init:t_fin) = reshape(nodes(mask, :)', 12 * my_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*my_mz: n_fin + 3*my_mz);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*my_mz: n_fin + 6*my_mz);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*my_mz: n_fin + 9*my_mz);
        elseif(fem_face == 4)
            t_init = t_fin;
            t_fin  = t_fin + 12 * my_mz;
            n_init = n_fin;
            n_fin  = n_fin + 3 * my_mz;
            nodes  = fem_nodes(fem_node_cnct(:, [1, 2, 4, 3]), 1:3);
            mask   = fem_nodes(fem_node_cnct(:, [1, 2, 4, 3]), 1) == max_x;
            x3x6(1+t_init:t_fin) = reshape(nodes(mask, :)', 12 * my_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*my_mz: n_fin + 3*my_mz);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*my_mz: n_fin + 6*my_mz);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*my_mz: n_fin + 9*my_mz);
        elseif(fem_face == 5)
            t_init = t_fin;
            t_fin  = t_fin + 12 * mx_my;
            n_init = n_fin;
            n_fin  = n_fin + 3 * mx_my;
            nodes  = fem_nodes(fem_node_cnct(:, [5, 6, 1, 2]), 1:3);
            mask   = fem_nodes(fem_node_cnct(:, [5, 6, 1, 2]), 3) == min_z;
            x3x6(1+t_init:t_fin) = reshape(nodes(mask, :)', 12 * my_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*mx_my: n_fin + 3*mx_my);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*mx_my: n_fin + 6*mx_my);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*mx_my: n_fin + 9*mx_my);
        else
            t_init = t_fin;
            t_fin  = t_fin + 12 * mx_my;
            n_init = n_fin;
            n_fin  = n_fin + 3 * mx_my;
            nodes = fem_nodes(fem_node_cnct(:, [4, 3, 8, 7]), 1:3);
            mask  = fem_nodes(fem_node_cnct(:, [4, 3, 8, 7]), 3) == max_z;
            x3x6(1+t_init:t_fin) = reshape(nodes(mask, :)', 12 * my_mz, 1);
            x3  (1 + n_init: n_fin) = x3x6(1 + n_init          : n_fin          );
            x4  (1 + n_init: n_fin) = x3x6(1 + n_init + 3*mx_my: n_fin + 3*mx_my);
            x5  (1 + n_init: n_fin) = x3x6(1 + n_init + 6*mx_my: n_fin + 6*mx_my);
            x6  (1 + n_init: n_fin) = x3x6(1 + n_init + 9*mx_my: n_fin + 9*mx_my);
        end %if
    end %for
end %function

% gamma = [gammat;gammaMixed]; ndis = 1; xdis=dx/2,ydis=dy/2,zdis=dz/2;

% function ftilda = traction(gamma,segments,xnodes, mno, a, MU, NU)                                 
%             
% % t_find traction on boundary due to dislocations
% % 
% %
% % rectangulare domain. Note for shared nodes. Set priority is 
% % 1) left (x=0)/right  (x=dx)
% % 2) top (z=dz)/bottom  (z=0)
% % 3) front (y=0)/back (y=dy)
% %
% %-----------------------------------------------
% %           
% %                   (mx,my)
% %                 
% %-----------------------------------------------
% 
% % modify  boundary conditions for concave domain
% 
% ftilda = zeros(mno*3,1);
% 
% nodes = gamma(:,1);
% area = gamma(:,2);
% normal = gamma(:,3:5);
% 
% lseg = size(segments,1);
% lgrid = size(nodes,1);
% 
% p1x = segments(:,6);
% p1y = segments(:,7);
% p1z = segments(:,8);
% 
% p2x = segments(:,9);
% p2y = segments(:,10);
% p2z = segments(:,11);
% 
% bx = segments(:,3);
% by = segments(:,4);
% bz = segments(:,5);
% 
% x = xnodes(nodes,1);
% y = xnodes(nodes,2);
% z = xnodes(nodes,3);
% 
% [sxx, syy, szz, sxy, syz, sxz] = StressDueToSegs(lgrid, lseg,...
%                                               x, y, z,...
%                                               p1x,p1y,p1z,...
%                                               p2x,p2y,p2z,...
%                                               bx,by,bz,...
%                                               a,MU,NU); 
%                                           
% Tx = sxx.*normal(:,1) + sxy.*normal(:,2) + sxz.*normal(:,3);
% Ty = sxy.*normal(:,1) + syy.*normal(:,2) + syz.*normal(:,3);
% Tz = sxz.*normal(:,1) + syz.*normal(:,2) + szz.*normal(:,3);
% 
% ATx = area.*Tx;
% ATy = area.*Ty;
% ATz = area.*Tz;
% 
% %populate ftilda
% ftilda(3*nodes-2)=ATx;
% ftilda(3*nodes-1)=ATy;
% ftilda(3*nodes)=ATz;
% % for j=1:lgrid
% %     gn=nodes(j);
% %     ftilda(3*gn-2) = ATx(j);
% %     ftilda(3*gn-1) = ATy(j);
% %     ftilda(3*gn) = ATz(j);
% % end
%     
% end
% 
% 
%  
