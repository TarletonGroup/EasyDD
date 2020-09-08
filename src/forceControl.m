function [u, u_bar, f, f_bar, calculateR_hat] = forceControl(u, u_dot, sign_u_dot, f, f_dot, sign_f_dot, simTime, gamma_mixed)
    u_bar = 0;
    f_bar = f_dot * simTime;
    f(3 * gamma_mixed(:, 1)) = sign_f_dot * f_bar;
    calculateR_hat = false;
end
