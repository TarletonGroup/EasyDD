function plotSimulation(Usim, Fsim, rn, links, plim, vertices, plotFreq, viewangle, plotForceDisp, amag, mumag, curstep)

    if (mod(curstep, plotFreq) == 0)
        figure(1)
        plotnodes(rn, links, plim, vertices);
        view(viewangle);
        drawnow

        feval(plotForceDisp, Usim, Fsim, amag, mumag, curstep);
        pause(0.01);
    end

end
