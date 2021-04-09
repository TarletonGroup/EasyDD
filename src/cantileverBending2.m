function [K, L, U, P_l, P_u, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, gammau,...
    gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = cantileverBending2(kg, ~, ~, ~, mno, ~, ~, ~, S)

    Sleft = [S.left; S.topleft; S.botleft; S.frontleft; S.backleft; S.corners([1,3,5,7],:)];
    for i = 1:size(Sleft,1)
        Sleft(i, 3:5) = [-1 0 0];
    end

    Smixed = [S.topright; S.corners([6, 8], :)];
    for i = 1:size(Smixed,1)
        Smixed(i, 3:5) = [1 0 0];
    end

    for i = [2,4,6,8]
        S.corners(i,3:5) = [1 0 0];
    end

    Stop = [S.top; S.topfront; S.topback];
    for i = 1:size(Stop,1)
        Stop(i, 3:5) = [0 0 1];
    end

    Sbot = [S.bot; S.botfront; S.botback];
    for i = 1:size(Sbot,1)
        Sbot(i, 3:5) = [0 0 -1];
    end

    Sfront = [S.front];
    Sback = [S.back];
    Sright = [S.right; S.botright; S.frontright; S.backright; S.corners([2, 4], :)];
    for i = 1:size(Sright)
        Sright(i, 3:5) = [1 0 0];
    end

    srfSet = cat(1, Sleft, Smixed, Sright, Stop, Sbot, Sfront, Sback);

    gammau = Sleft;
    gammaMixed = Smixed;
    [~, gammatIdx] = setdiff(srfSet(:,1),[gammau(:,1); gammaMixed(:,1)]);
    gammat = srfSet(gammatIdx, :);

    % Sleft = [S.left(:,1);
    % S.botleft(:,1);
    % S.topleft(:,1);
    % S.frontleft(:,1);
    % S.backleft(:,1);
    % S.corners([1,3,5,7], 1)];
    % 
    % Sback = [S.backleft(:,1);
    %     S.corners([3,7], 1)];
    % 
    % Sfront = [S.botleft(:,1);
    %     S.corners([1,3], 1)];
    % 
    % Sright = [S.corners([2,4,6,8], 1)];

    fixedDofs = [
        3 * S.left(:, 1) - 2;
        3 * S.botleft(:, 1) - 2;
        3 * S.topleft(:, 1) - 2;
        3 * S.frontleft(:, 1) - 2;
        3 * S.backleft(:, 1) - 2;
        3 * S.corners([1, 3, 5, 7], 1) - 2;

        3 * S.left(:, 1) - 1;
        3 * S.botleft(:, 1) - 1;
        3 * S.topleft(:, 1) - 1;
        3 * S.frontleft(:, 1) - 1;
        3 * S.backleft(:, 1) - 1;
        3 * S.corners([1, 3, 5, 7], 1) - 1;

        3 * S.left(:, 1);
        3 * S.botleft(:, 1);
        3 * S.topleft(:, 1);
        3 * S.frontleft(:, 1);
        3 * S.backleft(:, 1);
        3 * S.corners([1, 3, 5, 7], 1);

        3 * S.topright(:, 1);
        3 * S.corners([6, 8], 1)
    ];


    % fixedDofs = [3*gammau(:, 1) - 2; 3*gammau(:, 1) - 1; 3*gammau(:, 1); 3*gammaMixed(:,1) - 2];
    freeDofs = setdiff([1:3*mno], fixedDofs);

    K = kg(freeDofs, freeDofs);

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
end