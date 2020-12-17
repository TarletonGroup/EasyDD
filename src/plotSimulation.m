function plotSimulation(Usim, Fsim, rn, links, plim, vertices, plotFreq, viewangle, plotForceDisp, amag, mumag, curstep, plotFlags)

    figNum = 0;
    
    if (mod(curstep, plotFreq) == 0)
        if plotFlags.nodes
            figNum = figNum + 1;
            figure(figNum)
            plotnodes(rn, links, plim, vertices);
            view(viewangle);
            drawnow
        end
        
        if plotFlags.secondary
            figNum = figNum + 1;
            figure(figNum)
            plotForceDisp(Usim, Fsim, amag, mumag, curstep);
        end
    end

end
