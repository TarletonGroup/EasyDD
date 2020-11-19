function [Fsim, Usim, t] = atRestForceDisp(Fsim, f_bar, f_hat, f_tilda, Usim, u_bar, u_hat, u_tilda, r_hat, ...
        gammaMixed, fixedDofs, freeDofs, curstep, simTime)

    t(curstep + 1) = simTime;
    fprintf('simTime = %d\n', simTime);

end
