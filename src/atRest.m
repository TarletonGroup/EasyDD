function [K, L, U, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, gammau,...
    gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = atRest(...
    kg, w, h, d, mx, my, mz)

    fprintf('At rest on a surface boundary conditions: reformatting K\n');

    K = kg;
    bcwt = mean(diag(kg)); %=trace(K)/length(K)
    bcwt = full(bcwt);

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
    
    gammau = Sleft; % u=0
    gammat = [Stop; Sbot; Sright; Sfront; Sback; Smixed]; %  set containg fixed nodes t=0
    gammaMixed = []; % t1=t2=0, u3= U

    fixedDofs = [3 * gammau(:, 1) - 2; 3 * gammau(:, 1) - 1; 3 * gammau(:, 1)];
    freeDofs = [3 * gammat(:, 1) - 2; 3 * gammat(:, 1) - 1; 3 * gammat(:, 1)];
    
for m = 1:length(fixedDofs)
    i = fixedDofs(m);
    K(:, i) = 0;
    K(i, :) = 0;
    K(i, i) = bcwt;
end

if length([fixedDofs; freeDofs]) > length(unique([fixedDofs; freeDofs]))
    disp('error')
    pause
end


% disp('LU decomposition of K...')
% %L=0;U=0;
% tic;
% [L,U] = lu(K);
% toc;

disp('Cholesky Factorization of K...'); %should be symmetric!
tic;
U = chol(K);
L = U';
toc;

processForceDisp = 'atRestForceDisp';
plotForceDisp = 'atRestPlot';

fprintf('finished FEM\n')
end