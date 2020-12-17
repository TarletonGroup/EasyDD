function [u, u_bar, f, f_bar, calculateR_hat] = displacementControl(u, u_dot, sign_u_dot, u_bar_0, f, f_dot, sign_f_dot, f_bar_0, simTime, gamma_mixed)
    f_bar = 0;
    u_bar = u_dot * simTime + u_bar_0;
    u(3 * gamma_mixed(:, 1)) = sign_u_dot * u_bar;
    calculateR_hat = true;
end
