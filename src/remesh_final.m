function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_final(...
    rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
    u_hat, ...
    matpara, mods, flags, FEM, Bcoeff)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remesh rule
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
    
    %% Extraction
    
    % flags:
    doremesh = flags.doremesh;
    
    %% Function
    
    if doremesh
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh(...
            rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
            u_hat, ...
            matpara, mods, flags, FEM, Bcoeff);
    end
    
end
