function [Fsim, Usim, t] = cantileverBendingForceDisp(Fsim, f_bar, f_hat, f_tilda, Usim, u_bar, u_hat, u_tilda, r_hat, ...
        gammaMixed, fixedDofs, freeDofs, curstep, simTime)

    f_out = sum(r_hat(3 * gammaMixed(:, 1)) + f_tilda(3 * gammaMixed(:, 1)));
    u_out = u_bar;

    Fsim(curstep + 1) = f_out;
    Usim(curstep + 1) = u_out;
    t(curstep + 1) = simTime;

    fprintf('f_out = %d, u_out = %d, simTime = %d\n', f_out, u_out, simTime);

end