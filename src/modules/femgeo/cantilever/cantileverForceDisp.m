function [Fsim, Usim, t] = cantileverForceDisp(BCstore, ...
        f, f_hat, f_tilda, u, u_hat, u_tilda, r_hat, ...
        gamma, fixedDofs, freeDofs, curstep, simTime)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk
    
    % Prints and stores applied forces and displacements into the Fsim and
    % Usim objects. These are later plotted.
    %===============================================================%
    
    %%
    
    % Sum concentrated forces:
    f_out = sum(r_hat(3*gamma.Mixed(:, 1)) + f_tilda(3*gamma.Mixed(:, 1)));
    
    % FixedDofs corresponding to affinely increasing displacement:
    fixedDofs_U = 3*gamma.Mixed(:, 1);
    u_out = abs(u(fixedDofs_U(1)));

    BCstore.Fsim(curstep + 1) = f_out;
    BCstore.Usim(curstep + 1) = u_out;
    BCstore.t(curstep + 1) = simTime;

    fprintf('f_out = %d, u_out = %d, simTime = %d\n', f_out, u_out, simTime);

end