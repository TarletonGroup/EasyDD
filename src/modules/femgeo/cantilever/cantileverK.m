function [K, L, U, bcwt] = cantileverK(kg, fixedDofs, freeDofs)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Reformats stiffness matrix based on the prescribed
    % boundary conditions on the degrees of freedom.
    %===============================================================%
    
    %%
    
    fprintf('Cantilever boundary conditions: reformatting K.\n');
    
    % {freeDofs,fixedDofs} should contain every degree of freedom on
    % boundary + any internal Dofs with a force or displacement specified.
    if length([fixedDofs; freeDofs]) > length(unique([fixedDofs; freeDofs]))
        fprintf('error\n')
        pause
    end
    
    K = kg; % Reallocate stiffness matrix
    
    % Weighing factor to retain conditioning of K:
    bcwt = mean(diag(kg)); % = trace(K)/length(K)
    bcwt = full(bcwt); % Convert sparse matrix --> full matrix
    
    for m = 1:length(fixedDofs)
        i = fixedDofs(m);
        K(:, i) = 0;
        K(i, :) = 0;
        K(i, i) = bcwt;
    end
    % K should still be symmetric!
    
    fprintf('Cholesky factorisation of K...\n');
    tic;
    U = chol(K);
    L = U';
    toc;
    fprintf('Finished Cholesky factorisation of K.\n');
    
    processForceDisp = 'cantileverBendingForceDisp';
    plotForceDisp = 'cantileverBendingPlot';

end
