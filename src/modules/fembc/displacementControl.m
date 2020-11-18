function [u, u_bar, f, f_bar, calculateR_hat] = displacementControl(...
    u, u_dot, sign_u_dot, f, f_dot, sign_f_dot, simTime, ...
    gamma)
    
    % Deformation conditions:
    f_bar = 0;
    u_bar = u_dot * simTime;
    u(3 * gamma_mixed(:, 1)) = sign_u_dot * u_bar;
    
    calculateR_hat = true;
end
