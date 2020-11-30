function plotSimulation(rn, links, plotFreq, viewangle, curstep, ...
    FEM, saveBC)
    %===============================================================%
    % Oxford Materials (11/11/2020)
    
    % Constructs plots of the simulation.
    %===============================================================%
    
    %%
    
    if (mod(curstep, plotFreq) == 0)
        
        % Geometry and/or dislocation network plot:
        plotnodes(rn, links, FEM, viewangle);
        
        % Boundary conditions plot:
        plotBCs(saveBC, curstep);
    end
end
