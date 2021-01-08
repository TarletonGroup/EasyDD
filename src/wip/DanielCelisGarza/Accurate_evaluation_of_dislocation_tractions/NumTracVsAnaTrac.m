function [K, L, U, P_l, P_u, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, gammau,...
    gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = NumTracVsAnaTrac(...
    kg, w, h, d, mno, mx, my, mz)

    fprintf('Static simulation for numeric vs analytic traction comparison\n');
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

    gammat = Sleft; % t=0
    gammau = [Stop; Sbot; Sright; Sfront; Sback; Smixed]; %  set containg fixed nodes u=0
    gammaMixed = []; % t1=t2=0, u3= U

    fixedDofs = [3 * gammau(:, 1) - 2; 3 * gammau(:, 1) - 1; 3 * gammau(:, 1)];
    freeDofs = setdiff([1:3 * mno], fixedDofs);

    K = K(freeDofs, freeDofs);

%     figure(2); clf; hold on; view(3)
%     xlabel('x'); ylabel('y'); zlabel('z')
%     plot3(xnodes(Stop(:, 1), 1), xnodes(Stop(:, 1), 2), xnodes(Stop(:, 1), 3), 'r*')
%     plot3(xnodes(Sbot(:, 1), 1), xnodes(Sbot(:, 1), 2), xnodes(Sbot(:, 1), 3), 'r*')
%     plot3(xnodes(Sright(:, 1), 1), xnodes(Sright(:, 1), 2), xnodes(Sright(:, 1), 3), 'b.')
%     plot3(xnodes(Sleft(:, 1), 1), xnodes(Sleft(:, 1), 2), xnodes(Sleft(:, 1), 3), 'b.')
%     plot3(xnodes(Sfront(:, 1), 1), xnodes(Sfront(:, 1), 2), xnodes(Sfront(:, 1), 3), 'k*')
%     plot3(xnodes(Sback(:, 1), 1), xnodes(Sback(:, 1), 2), xnodes(Sback(:, 1), 3), 'k*')
%     plot3(xnodes(Smixed(:, 1), 1), xnodes(Smixed(:, 1), 2), xnodes(Smixed(:, 1), 3), 'g*')
%     axis('equal')
%     hold off

    disp('Cholesky Factorization of K...'); %should be symmetric!
    % Special algorithm for sparse matrices
    % [R, flag, P] = chol(S)
    % R'*R = P'*S*P -> P*R'*R*P' = S
    try
        [U, ~, P_u] = chol(K);
        L = U';
        P_l = P_u';
    catch
        sprintf('Ran out of memory in cholesky factorisation, use explicit K.')
        U = [];
        L = [];
        P_u = [];
        P_l = [];
    end
    processForceDisp = @processNothing;
    plotForceDisp = @plotNoForceDisp;
end