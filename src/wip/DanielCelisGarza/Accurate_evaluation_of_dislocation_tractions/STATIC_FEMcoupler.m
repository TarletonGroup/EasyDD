function [uhat, fend, Ubar, ftilda] = STATIC_FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, K, L, U, P_l, P_u,...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, t)

    %Coupling of FEM and DDD
    % u = uhat + utilda
    % f = fhat + ftilda

    segments = constructsegmentlist(rn, links, true);

    Udot = 0; %0.01*1E3*dx*(1E-4/160E9);%for tungsten...
    % Udot = 100*1E3*dx*(1E-4/160E9); %for tungsten...
    %Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle

    Ubar = Udot * t;
    %Ubar = 0.1*1E4; for debuggin
    u = zeros(3 * (mno), 1);
    % gamma=[gammau;gammaMixed];

    % u(3*gammaMixed(:,1)) = -Ubar;  %applied displacements in z at right edge nodes

    uhat = zeros(3 * mno, 1);
    % utilda=zeros(3*mno,1);

    % gn = gamma(:,1); % global node number
    % x0 = xnodes(gn,1:3); % field point
    % point_array_length = size(x0,1);
    % segments_array_length = size(segments,1);

    %Matlab wrapper
    % tic;
    % displacements = displacement_fivel(x0,segments,NU); %utilda array ordered by x0
    % toc;

    %Full C version (including combinatorials, checks, corrections etc.)
    % tic;
    % [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
    %                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
    %                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
    %                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
    %                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
    %                        NU,point_array_length,segments_array_length);
    % % displacements = horzcat(Ux,Uy,Uz);
    % % toc;
    % % disp(displacementsMEX-displacements);
    % %  [Ux, Uy, Uz]  = displacement_fivel(x0,segments,NU);
    %
    % utilda(3*gn -2) = Ux;
    % utilda(3*gn -1) = Uy;
    % utilda(3*gn   ) = Uz;
    %
    % if any(isnan(utilda))
    %     disp('some entries of utilda are NaN')
    %     pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
    % end
    % Commenting out utilda for the simulation.
    uhat(fixedDofs) = u(fixedDofs); % - utilda(fixedDofs);

    fhat = zeros(3 * (mno), 1);
    ftilda = zeros(3 * (mno), 1);

    % gamma=[gammat;gammaMixed];
    gamma = gammat;

    %ftilda = zeros(3*mno,1);
    ftilda = numeric_traction(gamma, segments, xnodes, a, MU, NU, ftilda);

    %%
    %ftilda=zeros(3*mno,1); %ET uncomment later!

    fhat(freeDofs) = -ftilda(freeDofs); % no applied forces

    %f=zeros(2*(mno),1);
    f = fhat - kg(:, fixedDofs) * uhat(fixedDofs);

    bcwt = mean(diag(kg)); %=trace(K)/length(K)
    bcwt = full(bcwt);

    % f(fixedDofs) = bcwt*uhat(fixedDofs);
    if isempty(U)
        uhat = K \ f;
    else
        uhat = P_l \ (U \ (L \ (P_u \ f))); % using LU decomposition for sparse matrices
    end

%     uhat(freeDofs) = K(freeDofs,freeDofs)\f(freeDofs);
    

    % uhat2=K\f;

    rhat = kg * uhat; % reaction force

    % fend = rhat(3*gammaMixed(:,1))+ftilda(3*gammaMixed(:,1));
    fend = rhat;
    % fend = sum(fend);

end
