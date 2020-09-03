function cantileverBendingPlot(Usim, Fsim, amag, mumag, curstep)
    figure(2)
    plot(Usim(1:curstep) * amag, -Fsim(1:curstep) * amag^2 * mumag);
    drawnow
end
