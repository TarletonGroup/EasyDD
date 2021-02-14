function micropillarTensilePlot(Usim, Fsim, amag, mumag, curstep, args)
    factDisp = args.factDisp;
    factForce = args.factForce;
    figure(2)
    plot(Usim(1:curstep-1)*factDisp,Fsim(1:curstep-1)*factForce)
    drawnow
end
