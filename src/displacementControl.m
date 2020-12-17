function [u, u_bar, f, f_bar, calculateR_hat] = displacementControl(u, u_dot, sign_u_dot, u_bar_0, f, f_dot, sign_f_dot, f_bar_0, simTime, holdingTime, gamma_mixed)
    f_bar = 0;
    if simTime < holdingTime
        u_bar = u_bar_0;
    else
        u_bar = u_dot * (simTime - holdingTime) + u_bar_0;
    end
    u(3 * gamma_mixed(:, 1)) = sign_u_dot * u_bar;
    calculateR_hat = true;
end