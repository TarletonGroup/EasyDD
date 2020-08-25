function [f_out, u_out] = processForceDisp(f_bar, f_hat, f_tilda, u_bar, u_hat, u_tilda, r_hat, gamma_disp, ...
        gamma_mixed, fixedDofs, freeDofs, simType)

    if simType == 1
        f_out = sum(r_hat(3 * gamma_mixed(:, 1)) + f_tilda(3 * gamma_mixed(:, 1)));
        u_out = u_bar;
    end

end
