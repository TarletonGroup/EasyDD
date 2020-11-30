function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remeshPreCollision(...
    rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
    u_hat, ...
    matpara, mods, flags, surfmesh, FEM, Bcoeff)
    %% Extraction
    
    % flags:
    doremesh = flags.doremesh;
    dovirtmesh = flags.dovirtmesh;
    
    %% Function

    if doremesh% do virtual re-meshing first
        %remeshing virtual dislocation structures
        if dovirtmesh
            [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = virtualmeshcoarsen(...
                rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
                matpara);
        end

        %remeshing internal dislocation structures
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(...
            rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
            u_hat, ...
            matpara, mods, flags, surfmesh, FEM, Bcoeff);
    end

end
