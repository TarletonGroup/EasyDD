% Remesh rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. This includes the original remesh rule.
% 2. Topological changes on free surface do not include force/velocity
% calculation. (Thus, just rearrangements of rn and links matrix)
%
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% January 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rn, links, ...
        connectivity, linksinconnect, fseg, lmin, lmax, areamin, areamax, MU, NU, a, Ec, mobility, rotMatrix, ...
        doremesh, dovirtmesh, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d, P, fn, CUDA_flag, Bcoeff)

    if dovirtmesh
        % Beginning of surface remeshing for surface node.
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_surf(rn, links, ...
            connectivity, linksinconnect, fseg, vertices, P, fn);
    else
        rnnew = rn;
        linksnew = links;
        connectivitynew = connectivity;
        linksinconnectnew = linksinconnect;
        fsegnew = fseg;
    end

    if doremesh
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh(rnnew, linksnew, ...
            connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, ...
            Ec, mobility, vertices, rotMatrix, uhat, nc, xnodes, D, mx, mz, w, h, d, CUDA_flag, Bcoeff);
    end

end
