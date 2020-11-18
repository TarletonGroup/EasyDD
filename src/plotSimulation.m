function plotSimulation(rn, links, plim, vertices, plotFreq, viewangle, plotForceDisp, BCstore, simscale, curstep)
    
    if (mod(curstep, plotFreq) == 0)
        figure(1)
        plotnodes(rn, links, plim, vertices);
        view(viewangle);
        drawnow
        
        feval(plotForceDisp, BCstore, simscale, curstep);
        pause(0.01);
    end

end
