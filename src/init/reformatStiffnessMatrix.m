function [Kdecomp] = reformatStiffnessMatrix(FEM, dofs)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Reformats stiffness matrix based on the prescribed
    % boundary conditions on the degrees of freedom.
    %===============================================================%
    
    %% Extraction
    
    % FEM:
    kg = FEM.kg;
    
    % DoFs:
    fixedDofs = dofs.fixedDofs;
    freeDofs = dofs.freeDofs;
    
    %% Reformatting
    
    fprintf('Cantilever boundary conditions: reformatting K.\n');
    
    % {freeDofs,fixedDofs} should contain every degree of freedom on
    % boundary + any internal Dofs with a force or displacement specified.
    if length([fixedDofs; freeDofs]) > length(unique([fixedDofs; freeDofs]))
        error('length([fixedDofs; freeDofs]) > length(unique([fixedDofs; freeDofs]))\n')
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
    Uchol = chol(K);
    Lchol = Uchol';
    toc;
    fprintf('Finished Cholesky factorisation of K.\n');
    
    %% Data storage
    
    Kdecomp = struct; % Stores stiffness matrix decomposition data
    
    Kdecomp.bcwt = bcwt;
    Kdecomp.K = K;
    Kdecomp.Uchol = Uchol;
    Kdecomp.Lchol = Lchol;
    
end
