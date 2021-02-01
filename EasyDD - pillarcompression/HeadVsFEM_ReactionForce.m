%%% Script to compare Head's analytical solutions to those obtained via FEM using analytic and numeric tractions.
close all
amag=3.18e-4; 

CRYSTAL_STRUCTURE = 'bcc';
% FEM PARAMETERS
%Cantilever
simTime = 0;
use_gpu=0;
para_scheme = 0;
n_threads = 0;

%MU = 160E9; %160GPa
MU = 1;
NU = 0.28;

% Lattice constants
a = 10;
bVec = [-[1 0 0]; -[0 0 1]];


planes = [1;2;3;4;5;6];
dx = 1000;
dy = 1000;
dz = 1000;

figCounter = 0;
cntr = 0;
    
for j = 20
    mx=j;

    gridSize = mx;
    x = linspace(0, dx, gridSize);
    y = linspace(0.5*dy, 0.5*dy, gridSize);
    z = linspace(0, dz, gridSize);
    [X,Z] = meshgrid(x,z);
    Y = meshgrid(y);

    clear x y z;

    loading=1; 
    vertices = [0,0,0;...
                dx,0,0;...
                0,dy,0;...
                dx,dy,0;...
                0,0,dz;...
                dx,0,dz;...
                0,dy,dz;...
                dx,dy,dz];

    plim=12/amag; %12microns
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,Smixed,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,my,mz,mel] = STATIC_finiteElement3D(dx,dy,dz,mx,MU,NU,loading);

    gamma_dln = [gammat(:,1)];

     % Set surface node labels for surface node extraction.
    n_nodes = 4;
    surf_node_util = zeros(n_nodes+2, 6);
    xy = mx*my;
    xz = mx*mz;
    yz = my*mz;

    f_hat = zeros(3*mno, 1);
    [x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
                                                   planes, 4);
    [f_dln_node, f_dln_se,...
     f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
    tolerance = dx/10^6;

    figCounter = figCounter + 1;
    figure(figCounter)
    clf;hold on;view(3)
    xlabel('x');ylabel('y');zlabel('z')

    plot3(xnodes(Stop(:,1),1),xnodes(Stop(:,1),2),xnodes(Stop(:,1),3),'r*')
    plot3(xnodes(Sbot(:,1),1),xnodes(Sbot(:,1),2),xnodes(Sbot(:,1),3),'r*')
    plot3(xnodes(Sright(:,1),1),xnodes(Sright(:,1),2),xnodes(Sright(:,1),3),'b.')
    plot3(xnodes(Sleft(:,1),1),xnodes(Sleft(:,1),2),xnodes(Sleft(:,1),3),'b.')
    plot3(xnodes(Sfront(:,1),1),xnodes(Sfront(:,1),2),xnodes(Sfront(:,1),3),'k*')
    plot3(xnodes(Sback(:,1),1),xnodes(Sback(:,1),2),xnodes(Sback(:,1),3),'k*')
    plot3(xnodes(Smixed(:,1),1),xnodes(Smixed(:,1),2),xnodes(Smixed(:,1),3),'g*')
    axis('equal')
    hold off

    len = 100;
    x = linspace(0.1*dx,0.1*dx,len);
    y = linspace(0,dy,len);
    z = linspace(0.5*dz,0.5*dz,len);
    x1 = x(1);
    z1 = z(1);
    t = [0 1 0];
    n = [0 0 1];
    rn = zeros(len,3);
    rn(:,1) = x;
    rn(:,2) = y;
    rn(:,3) = z;
    links = zeros(len-1,8);
    hold on
        plot3(rn(:,1),rn(:,2),rn(:,3),'r.')
    hold off
end
%%
    for k = 2
        cntr = cntr + 1;
        b = bVec(k, :);
        for i = 1:len-1
            links(i,:) = [i, i+1, b,  n];
        end

        [uhat,fend,Ubar,fan] = STATIC_analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
                   0, 0, gammaMixed,fixedDofs,freeDofs,dx,simTime ,...
                   gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
                   f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);

        [uhat2,fend2,Ubar2,fnum] = STATIC_FEMcoupler(rn,links,0,a,MU,NU,xnodes,mno,kg,L,U,...
                             gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,simTime);

        sigmaA = hatStressSurf(uhat,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);
        sxxA = squeeze(sigmaA(1,1,:,:));
        szzA = squeeze(sigmaA(3,3,:,:));
        sxzA = squeeze(sigmaA(1,3,:,:));
        sigmaN = hatStressSurf(uhat2,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);
        sxxN = squeeze(sigmaN(1,1,:,:));
        szzN = squeeze(sigmaN(3,3,:,:));
        sxzN = squeeze(sigmaN(1,3,:,:));
        
        segments = constructsegmentlist(rn,links);
        p1 = [segments(:,6) segments(:,7) segments(:,8)];
        p2 = [segments(:,9) segments(:,10) segments(:,11)];

        sigmaFP = FieldPointStressSurf(X,Y,Z,p1,p2,b,a,MU,NU);
        sxxFP = squeeze(sigmaFP(1,1,:,:));
        szzFP = squeeze(sigmaFP(3,3,:,:));
        sxzFP = squeeze(sigmaFP(1,3,:,:));

        x = linspace(0,dx,gridSize);
        z = linspace(0,dz,gridSize);
        b = sqrt(3)/2;
        [X,Z] = meshgrid(x,z);
        if k == 1
            [txx,tzz,txz] = imageStressAnalyticEdgePerp(MU, b, NU, X, Z, x1, z1);
            [txxFP, tzzFP, txzFP] = FPStressAnalyticEdgePerp(MU, b, NU, X, Z, x1, z1);
            name = ' perp';
        else
            [txx,tzz,txz] = imageStressAnalyticEdgePar(MU, b, NU , X, Z, x1, z1);
            [txxFP, tzzFP, txzFP] = FPStressAnalyticEdgePar(MU, b, NU, X, Z, x1, z1);
            name = ' par';
        end
        close all

%         save(sprintf('./mat_files/headVsFEM_%d', cntr))

%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,sxxA);
%             colormap(parula)
%             colorbar
%             title('A sxx')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%     
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,szzA);
%             colormap(parula)
%             colorbar
%             title('A szz')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%     
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,sxzA);
%             colormap(parula)
%             colorbar
%             title('A sxz')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%             
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,sxxN);
%             colormap(parula)
%             colorbar
%             title('N sxx')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%     
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,szzN);
%             colormap(parula)
%             colorbar
%             title('N szz')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%             
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             contourf(X,Z,sxzN);
%             colormap(parula)
%             colorbar
%             title('N sxz')
%             xlabel('x')
%             ylabel('z')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%             
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             txx = txx;%./norm(txx);
%             meantxx = mean(txx,'all');
%             stdtxx = std(txx,0,'all');
%             displace = 5*stdtxx;
%             limits = [meantxx - displace, meantxx + displace];
%             contourf(X,Z,txx);
%             colormap(parula)
%             colorbar
%             caxis(limits)
%             title(strcat('b', name, ' sxx'))
%             xlabel('b')
%             ylabel('b')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
%             
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             tzz = tzz;%./norm(tzz);
%             meantzz = mean(tzz,'all');
%             stdtzz = std(tzz,0,'all');
%             displace = 5*stdtzz;
%             limits = [meantzz - displace, meantzz + displace];
%             contourf(X,Z,tzz);
%             colormap(parula)
%             colorbar
%             caxis(limits)
%             title(strcat('b', name, ' szz'))
%             xlabel('b')
%             ylabel('b')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
% 
%             figCounter = figCounter + 1;
%             figure(figCounter)
%             txz = txz;%./norm(txz);
%             meantxz = mean(txz,'all');
%             stdtxz = std(txz,0,'all');
%             displace = 5*stdtxz;
%             limits = [meantxz - displace, meantxz + displace];
%             contourf(X,Z,txz);
%             colormap(parula)
%             colorbar
%             caxis(limits)
%             title(strcat('b', name, ' sxz'))
%             xlabel('b')
%             ylabel('b')
%             hold on
%             plot(x1,z1,'.','color','black','MarkerSize',10)
%             hold off
end
% end

%%

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxxA);
colormap(parula)
colorbar
title('A sxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,szzA);
colormap(parula)
colorbar
title('A szz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxzA);
colormap(parula)
colorbar
title('A sxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxxN);
colormap(parula)
colorbar
title('N sxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,szzN);
colormap(parula)
colorbar
title('N szz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxzN);
colormap(parula)
colorbar
title('N sxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
txx = txx;%./norm(txx);
meantxx = mean(txx,'all');
stdtxx = std(txx,0,'all');
displace = 5*stdtxx;
limits = [meantxx - displace, meantxx + displace];
contourf(X,Z,txx);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' sxx'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
tzz = tzz;%./norm(tzz);
meantzz = mean(tzz,'all');
stdtzz = std(tzz,0,'all');
displace = 5*stdtzz;
limits = [meantzz - displace, meantzz + displace];
contourf(X,Z,tzz);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' szz'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
txz = txz;%./norm(txz);
meantxz = mean(txz,'all');
stdtxz = std(txz,0,'all');
displace = 5*stdtxz;
limits = [meantxz - displace, meantxz + displace];
contourf(X,Z,txz);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' sxz'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off
%%

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxxFP);
colormap(parula)
colorbar
title('FP FE sxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,szzFP);
colormap(parula)
colorbar
title('FP FE szz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
contourf(X,Z,sxzFP);
colormap(parula)
colorbar
title('FP FE sxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off


figCounter = figCounter + 1;
figure(figCounter)
txxFP = txxFP;%./norm(txx);
meantxx = mean(txxFP,'all');
stdtxx = std(txxFP,0,'all');
displace = 5*stdtxx;
limits = [meantxx - displace, meantxx + displace];
contourf(X,Z,txxFP);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' FP sxx'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
tzzFP = tzzFP;%./norm(tzz);
meantzz = mean(tzzFP,'all');
stdtzz = std(tzzFP,0,'all');
displace = 5*stdtzz;
limits = [meantzz - displace, meantzz + displace];
contourf(X,Z,tzzFP);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' FP szz'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figCounter = figCounter + 1;
figure(figCounter)
txzFP = txzFP;%./norm(txz);
meantxz = mean(txzFP,'all');
stdtxz = std(txzFP,0,'all');
displace = 5*stdtxz;
limits = [meantxz - displace, meantxz + displace];
contourf(X,Z,txzFP);
colormap(parula)
colorbar
caxis(limits)
title(strcat('b', name, ' FP sxz'))
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off


function sigma = hatStressSurf(uhat,nc,x,D,mx,mz,w,h,d,X,Y,Z)

    gridSize = size(X);

    sigma = zeros(3,3,gridSize(1),gridSize(2));
    x0 = zeros(3,1);
    for col = 1:gridSize(2)
        for row = 1:gridSize(1)
            x0(1) = X(row,col);
            x0(2) = Y(row,col);
            x0(3) = Z(row,col);
            sigma(:,:,row,col) = hatStress(uhat,nc,x,D,mx,mz,w,h,d,x0);
        end
    end
end

function sigma = FieldPointStressSurf(X,Y,Z,x1,x2,b,a,mu,nu)
    gridSize = size(X);
    sigma = zeros(3,3,gridSize(1),gridSize(2));
    x0 = zeros(1,3);
    for col = 1:gridSize(2)
        for row = 1:gridSize(1)
            x0(1) = X(row,col);
            x0(2) = Y(row,col);
            x0(3) = Z(row,col);
            stress = FieldPointStress(x0,x1,x2,b,a,mu,nu);
            
            sigma(1,1,row,col) = stress(:, 1);
            sigma(2,2,row,col) = stress(:, 2);
            sigma(3,3,row,col) = stress(:, 3);
            sigma(1,2,row,col) = stress(:, 4);
            sigma(2,1,row,col) = stress(:, 4);
            sigma(2,3,row,col) = stress(:, 5);
            sigma(3,2,row,col) = stress(:, 5);
            sigma(1,3,row,col) = stress(:, 6);
            sigma(3,1,row,col) = stress(:, 6);
        end
    end
    
end

function [txx, tyy, txy] = imageStressAnalyticEdgePerp(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (x, y) induced by edge dislocation parallel to the 
    % surface at x = 0. Dislocation coordinates are (a, c).
    % b perpendicular to surface.
    %%%
    E = (2*(1+nu)) * mu;
    
    D = E.*b./(4.*pi.*(1-nu.^2));
    ymc = y-c;
    ymc2 = ymc.^2;
    xma = x-a;
    xma2 = xma.^2;
    xpa = x+a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx =       ...-ymc .* (3.*xma2 + ymc2)./den1 + ...
                ymc .* (3.*xpa2 + ymc2)./den2 + ...
      4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
    txx = D.*txx;

    tyy =     ... %ymc .* (xma2 - ymc2)./den1 + ...
             -ymc .* (xpa2 - ymc2)./den2 + ...
      4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
    tyy = D.*tyy;

    txy =       ... %xma .* (xma2 - ymc2)./den1 + ...
               -xpa .* (xpa2 - ymc2)./den2 + ...
      2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
    txy = D.*txy;
end

function [txx, tyy, txy] = FPStressAnalyticEdgePerp(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (x, y) induced by edge dislocation parallel to the 
    % surface at x = 0. Dislocation coordinates are (a, c).
    % b perpendicular to surface.
    %%%
    E = (2*(1+nu)) * mu;
    
    D = E.*b./(4.*pi.*(1-nu.^2));
    ymc = y-c;
    ymc2 = ymc.^2;
    xma = x-a;
    xma2 = xma.^2;
    xpa = x+a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx =       -ymc .* (3.*xma2 + ymc2)./den1;
%                 + ...
%                 ymc .* (3.*xpa2 + ymc2)./den2 + ...
%       4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
    txx = D.*txx;

    tyy =     ymc .* (xma2 - ymc2)./den1;
%               + ...
%              -ymc .* (xpa2 - ymc2)./den2 + ...
%       4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
    tyy = D.*tyy;

    txy =       xma .* (xma2 - ymc2)./den1;
%                 + ...
%                -xpa .* (xpa2 - ymc2)./den2 + ...
%       2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
    txy = D.*txy;
end

function [txx, tyy, txy] = imageStressAnalyticEdgePar(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (a, c) induced by edge dislocation parallel to the 
    % surface at x = 0. Dislocation coordinates are (x, y).
    % p parallel to surface.
    %%%
    E = (2*(1+nu)) * mu;
    D = E.*b./(4.*pi.*(1-nu.^2));
    ymc = y-c;
    ymc2 = ymc.^2;
    xma = x-a;
    xma2 = xma.^2;
    xpa = x+a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx =     ...%xma .* (xma2 - ymc2) ./ den1 + ...
             -xpa .* (xpa2 - ymc2) ./ den2 + ...
      2.*a .* (3 .* xpa2 .* xpa2 - 6.*x.*xpa.*ymc2 - ymc2.*ymc2) ./ den3;
    txx = D.*txx;

    tyy =       ...%xma .* (xma2 + 3.*ymc2) ./ den1 + ...
               -xpa .* (xpa2 + 3.*ymc2) ./ den2 + ...
      -2.*a .* (xma .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 + ymc2.*ymc2) ./ den3;
    tyy = D.*tyy;

    txy =        ...%ymc .* (xma2 - ymc2) ./ den1 + ...
                -ymc .* (xpa2 - ymc2) ./ den2 + ...
      4.*a.*x .* ymc .* (3.*xpa2 - ymc2) ./ den3;
    txy = D.*txy;
end

function [txx, tyy, txy] = FPStressAnalyticEdgePar(mu, b, nu, x, y, a, c)
    %%%
    % Stress on point (a, c) induced by edge dislocation parallel to the 
    % surface at x = 0. Dislocation coordinates are (x, y).
    % p parallel to surface.
    %%%
    E = (2*(1+nu)) * mu;
    D = E.*b./(4.*pi.*(1-nu.^2));
    ymc = y-c;
    ymc2 = ymc.^2;
    xma = x-a;
    xma2 = xma.^2;
    xpa = x+a;
    xpa2 = xpa.^2;
    den1 = (xma2 + ymc2).^2;
    den2 = (xpa2 + ymc2).^2;
    den3 = den2 .* (xpa2 + ymc2);

    txx =     xma .* (xma2 - ymc2) ./ den1; 
%             + ...
%              -xpa .* (xpa2 - ymc2) ./ den2 + ...
%       2.*a .* ymc .* ((3.*x+a) .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 - ymc2.*ymc2) ./ den3;
    txx = D.*txx;

    tyy =       xma .* (xma2 + 3.*ymc2) ./ den1;
%                 + ...
%                -xma .* (xpa2 + 3.*ymc2) ./ den2 + ...
%       -2.*a .* (xma .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 + ymc2.*ymc2) ./ den3;
    tyy = D.*tyy;

    txy =        ymc .* (xma2 - ymc2) ./ den1;
%                 + ...
%                 -ymc .* (xpa2 - ymc2) ./ den2 + ...
%       4.*a.*x .* ymc .* (3.*xpa2 - ymc2) ./ den3;
    txy = D.*txy;
end

