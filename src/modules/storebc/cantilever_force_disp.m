function [saveBC] = cantilever_force_disp(...
        u, u_hat, u_tilda, f, f_hat, f_tilda, r_hat, ...
        curstep, simTime, ...
        saveBC, gamma, dofs)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk
    
    % Prints and stores boundary conditions into the Ysim and Xsim
    % objects. This prepares a plot to be made for Xsim vs Ysim.
    %===============================================================%
        
    %% Calculate saved variables
    
    % Sum concentrated forces on cantilever edge:
    f_out = sum(r_hat(3*gamma.Mixed(:, 1)) + f_tilda(3*gamma.Mixed(:, 1)));
    
    % FixedDofs corresponding to affinely increasing displacement:
    fixedDofs_U = 3*gamma.Mixed(:, 1);
    u_out = abs(u(fixedDofs_U(1)));
    
    %% Store variables
    
    % Initialisation:
    if curstep == 0
        
        saveBC.t = zeros(1e6, 1);
        saveBC.Xsim = zeros(1e6, 1);
        saveBC.Ysim = zeros(1e6, 1);
        
        saveBC.Xlabel = 'Displacement';
        saveBC.Ylabel = 'Force';
        saveBC.figname = 'Bending cantilever edge BCs';
    end
    
    saveBC.t(curstep + 1) = simTime;
    saveBC.Xsim(curstep + 1) = u_out;
    saveBC.Ysim(curstep + 1) = f_out;
    
    %% Print information
    
    fprintf('f_out = %d, u_out = %d, simTime = %d\n', f_out, u_out, simTime);

end