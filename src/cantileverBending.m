function [K, L, U, P_l, P_u, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, gammau,...
    gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = cantileverBending(...
    kg, w, h, d, mx, my, mz)

    fprintf('Cantilever bending boundary conditions: reformatting K\n');

    K = kg;
    bcwt = mean(diag(kg)); %=trace(K)/length(K)
    bcwt = full(bcwt);

    % define empty sets (not used but allows same code as for bending)
    % priority for inclusion of shared nodes is in order below
    %S = [node #, nodal Area, outward normal]
    Smixed = zeros(my + 1, 5); % top right edge (u3=U) first priority

    Sleft = zeros((my + 1) * (mz + 1), 5); % second priority
    Sright = zeros((my + 1) * mz, 3);

    Stop = zeros((mx - 1) * (my + 1), 5); % third priority
    Sbot = zeros((mx - 1) * (my + 1), 5);

    Sfront = zeros((mz - 1) * (mx - 1), 5); % fourth priority
    Sback = zeros((mz - 1) * (mx - 1), 5);

    m = 1;

    for j = 1:my + 1

        for k = 1:mz + 1
            % left surface of beam (u=0)
            gn = 1 + (k - 1) * (mx + 1) + (mx + 1) * (mz + 1) * (j - 1);
            Sleft(m, 1) = gn;

            if k == 1 || k == mz + 1
                L3 = 0.5 * d;
            else
                L3 = d;
            end

            if j == 1 || j == my + 1
                L2 = 0.5 * h;
            else
                L2 = h;
            end

            Sleft(m, 2) = L2 * L3;
            Sleft(m, 3:5) = [-1, 0, 0];
            m = m + 1;
        end

    end

    m = 1;

    for j = 1:my + 1

        for k = 1:mz
            % right surface of beam (t=0)
            gn = k * (mx + 1) + (mx + 1) * (mz + 1) * (j - 1);
            Sright(m, 1) = gn;

            if k == 1
                L3 = 0.5 * d;
            else
                L3 = d;
            end

            if j == 1 || j == my + 1
                L2 = 0.5 * h;
            else
                L2 = h;
            end

            Sright(m, 2) = L2 * L3;
            Sright(m, 3:5) = [1, 0, 0];
            m = m + 1;
        end

    end

    m = 1;

    for j = 1:my + 1
        % top right loaded edge (t1=0, u2=-U) mixed surface
        gn = (mx + 1) * (mz + 1) * j;
        Smixed(m, 1) = gn;

        if j == 1 || j == my + 1
            L2 = 0.5 * h;
        else
            L2 = h;
        end

        L3 = 0.5 * d;
        Smixed(m, 2) = L2 * L3;
        Smixed(m, 3:5) = [1, 0, 0];
        m = m + 1;
    end

    m = 1;

    for j = 1:my + 1

        for i = 2:mx
            % bottom surface of beam (t=0)
            gn = i + (mx + 1) * (mz + 1) * (j - 1);
            Sbot(m, 1) = gn;
            % top surface of beam (t=0)
            gn = i + (mx + 1) * mz + (mx + 1) * (mz + 1) * (j - 1);
            Stop(m, 1) = gn;

            L1 = w;

            if j == 1 || j == my + 1
                L2 = 0.5 * h;
            else
                L2 = h;
            end

            Sbot(m, 2) = L1 * L2;
            Sbot(m, 3:5) = [0, 0, -1];
            Stop(m, 2) = L1 * L2;
            Stop(m, 3:5) = [0, 0, 1];
            m = m + 1;
        end

    end

    m = 1;

    for k = 2:mz

        for i = 2:mx
            %front surface of beam (t=0)
            gn = i + (k - 1) * (mx + 1);
            Sfront(m) = gn; %(y=0)

            gn = i + (k - 1) * (mx + 1) + (mx + 1) * (mz + 1) * my;
            Sback(m) = gn; %(y=dy)

            L1 = w;
            L3 = d;
            Sfront(m, 2) = L1 * L3;
            Sfront(m, 3:5) = [0, -1, 0];
            Sback(m, 2) = L1 * L3;
            Sback(m, 3:5) = [0, 1, 0];

            m = m + 1;
        end

    end

    gammat = [Stop; Sbot; Sright; Sfront; Sback]; % t=0
    gammau = Sleft; %  set containg fixed nodes u=0
    gammaMixed = Smixed; % t1=t2=0, u3= U

    fixedDofs = [3 * gammau(:, 1) - 2; 3 * gammau(:, 1) - 1; 3 * gammau(:, 1); 3 * gammaMixed(:, 1)];
    freeDofs = [3 * gammat(:, 1) - 2; 3 * gammat(:, 1) - 1; 3 * gammat(:, 1); 3 * gammaMixed(:, 1) - 2; ...
                3 * gammaMixed(:, 1) - 1];

    K(:, fixedDofs) = 0;
    K(fixedDofs, :) = 0;
    idx = logical(speye(size(K)));
    diagonal = K(idx);
    diagonal(fixedDofs) = bcwt;
    K(idx) = diagonal;

    % {freeDofs,fixedDofs} should contain every degree of freedom on boundary +
    % any internal Dofs with a force or displacement specified.
    if length([fixedDofs; freeDofs]) > length(unique([fixedDofs; freeDofs]))
        fprintf('error\n')
        pause
    end

    try
        fprintf('Cholesky Factorization of K...\n'); %should be symmetric!
        % Special algorithm for sparse matrices
        % [R, flag, P] = chol(S)
        % R'*R = P'*S*P -> P*R'*R*P' = S
        tic;
            [U, ~, P_u] = chol(K);
            L = U';
            P_l = P_u';
        toc
    catch
        sprintf('Ran out of memory in cholesky factorisation, use explicit K.\n')
        U = [];
        L = [];
        P_u = [];
        P_l = [];
    end
    
    processForceDisp = @cantileverBendingForceDisp;
    plotForceDisp = @cantileverBendingPlot;

    fprintf('finished FEM\n')

end
