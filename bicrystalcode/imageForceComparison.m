%%
% Static comparison of tractions+FEM vs analytic image forces.
%
%
%% SOURCE GENERATION PARAMETERS
close all
amag=3.18e-4; 
mumag = 145E3; % MPa only used for plotting  

CRYSTAL_STRUCTURE = 'bcc';
% FEM PARAMETERS
%Cantilever
simTime = 0;
use_gpu=0;
para_scheme = 0;
n_threads = 0;
dx=25;
dy=50;
dz=50;

mx=10; % number of elements along beam length
loading=1; 
vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];


%MU = 160E9; %160GPa
MU = 1;
NU = 0.28; 
lmin = 0.1/amag;
a=lmin/sqrt(3)*0.5;



plim=12/amag; %12microns
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,Smixed,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = STATIC_finiteElement3D(dx,dy,dz,mx,MU,NU,loading);
% save('./mat_files/StaticMesh')
%%
% yz plane with x=0
% planes = 1;
% gammat = Sleft;
% gammau = Sleft;
% gammaMixed = Sleft;

planes = [1;2;3;4;5;6];
% % gammat = [Sleft];%[Stop;Sbot;Sright;Sfront;Sback];%[Sleft;Sright;Stop;Sbot;Sfront;Sback;Smixed];
% % gammau = [Stop;Sright;Sbot;Smixed;Sfront;Sback];
% % gammaMixed = [];%Smixed

gamma_dln = [gammat(:,1)];%, gammaMixed(:,1)];
% % fixedDofs =[3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1)];
% % freeDofs = [3*gammat(:,1)-2; 3*gammat(:,1)-1; 3*gammat(:,1)];


len = 100; % good behaviour 1000
x = linspace(0.1*dx,0.1*dx,len);
y = linspace(0,dy,len);
z = linspace(0.5*dz,0.5*dz,len);
x1 = x(1);
z1 = z(1);
t = [0 1 0];
b = -[1 0 0];
% b = -[0 0 1];
n = [0 0 1];
rn = zeros(len,3);
rn(:,1) = x;
rn(:,2) = y;
rn(:,3) = z;
links = zeros(len-1,8);
for i = 1:len-1
    links(i,:) = [i, i+1, b,  n];
end

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


figure(25);clf;hold on;view(3)
xlabel('x');ylabel('y');zlabel('z')
%
% plot3(x3x6(1:3:end,1), x3x6(2:3:end,1), x3x6(3:3:end,1),'ko')
% plot3(x3x6(1:3:end,2), x3x6(2:3:end,2), x3x6(3:3:end,2),'ko')
% plot3(x3x6(1:3:end,3), x3x6(2:3:end,3), x3x6(3:3:end,3),'ko')
% plot3(x3x6(1:3:end,4), x3x6(2:3:end,4), x3x6(3:3:end,4),'ko')
plot3(xnodes(Stop(:,1),1),xnodes(Stop(:,1),2),xnodes(Stop(:,1),3),'r*')
plot3(xnodes(Sbot(:,1),1),xnodes(Sbot(:,1),2),xnodes(Sbot(:,1),3),'r*')
plot3(xnodes(Sright(:,1),1),xnodes(Sright(:,1),2),xnodes(Sright(:,1),3),'b.')
plot3(xnodes(Sleft(:,1),1),xnodes(Sleft(:,1),2),xnodes(Sleft(:,1),3),'b.')
plot3(xnodes(Sfront(:,1),1),xnodes(Sfront(:,1),2),xnodes(Sfront(:,1),3),'k*')
plot3(xnodes(Sback(:,1),1),xnodes(Sback(:,1),2),xnodes(Sback(:,1),3),'k*')
plot3(xnodes(Smixed(:,1),1),xnodes(Smixed(:,1),2),xnodes(Smixed(:,1),3),'g*')
plot3(rn(:,1),rn(:,2),rn(:,3),'r.')
axis('equal')
hold off

% surf_node_util(1:6, 1) = [1, 5, 4, 8, yz, 1]; % min(x), yz-plane, face 1 ~Sleft
% surf_node_util(1:6, 2) = [6, 2, 7, 3, yz, 1]; % max(x), yz-plane, face 2 ~Sright
% surf_node_util(1:6, 3) = [5, 6, 8, 7, xz, 2]; % min(y), xz-plane, face 3 ~Sfront
% surf_node_util(1:6, 4) = [2, 1, 3, 4, xz, 2]; % max(y), xz-plane, face 4 ~Sback
% surf_node_util(1:6, 5) = [6, 5, 2, 1, xy, 3]; % min(z), xy-plane, face 5 ~Sbot
% surf_node_util(1:6, 6) = [3, 4, 7, 8, xy, 3]; % max(z), xy-plane, face 6 ~Stop
% 4. ----- .3
%  |\       |\
%  | \      | \
% 1. -\---- .2 \
%   \  \     \  \
%    \ 8. ----\- .7
%     \ |      \ |
%      \|       \|
%      5. ----- .6
% ^z
% |  y
% | /
% |/
% |------>x
%        


[uhat,fend,Ubar,fan] = STATIC_analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
                       0, 0, gammaMixed,fixedDofs,freeDofs,dx,simTime ,...
                       gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
                       f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);
                   
[uhat2,fend2,Ubar2,fnum] = STATIC_FEMcoupler(rn,links,0,a,MU,NU,xnodes,mno,kg,L,U,...
                     gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,simTime);
% [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
%                        gamma_disp, gammat, gammaMixed,fixedDofs,freeDofs,dx,simTime ,...
%                        gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
%                        f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance, Ubar, dt);
% segments=constructsegmentlist(rn,links);
% image_force = pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);
% image_force2 = pkforcevec(uhat2,nc,xnodes,D,mx,mz,w,h,d,segments);

%%
% gridSize = mx;
% X = linspace(0,dx,gridSize)';
% Y = linspace(0.5*dy,0.5*dy,gridSize)';
% Z = linspace(0,dz,gridSize)';
% 
% X_size = length(X);
% Y_size = length(Y);
% Z_size = length(Z);
% 
% Sxxu = zeros(X_size,Y_size);
% Syyu = zeros(X_size,Y_size);
% Szzu = zeros(X_size,Y_size);
% Sxyu = zeros(X_size,Y_size);
% Sxzu = zeros(X_size,Y_size);
% Syzu = zeros(X_size,Y_size);
% 
% 
% for i= 1:X_size
%         for j = 1:Z_size
% %             x0 = [X(i) 0.5*dy Z(j)]; % field point
%             x0 = [X(i) 0.5*dy Z(j)]; % field point
%             sigmahat = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
%             Sxxu(i,j) = sigmahat(1,1);
%             Syyu(i,j) = sigmahat(2,2);
%             Szzu(i,j) = sigmahat(3,3);
% 
%             Sxyu(i,j) = sigmahat(1,2); %isotropic
%             Sxzu(i,j) = sigmahat(1,3); %isotropic
%             Syzu(i,j) = sigmahat(2,3); %isotropic
%         end
% end
% figure(1)
% surf(X,Z,Sxxu','EdgeColor','none');
%
gridSize = mx;
x = linspace(0, dx, gridSize);
y = linspace(0.5*dy, 0.5*dy, gridSize);
z = linspace(0, dz, gridSize);
[X,Z] = meshgrid(x,z);
Y = meshgrid(y);

clear x y z;

sigma = hatStressSurf(uhat,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);

sxxA = squeeze(sigma(1,1,:,:));
figure(1)
contourf(X,Z,sxxA);
colormap(parula)
colorbar
title('A sxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

szzA = squeeze(sigma(3,3,:,:));
figure(2)
contourf(X,Z,szzA);
colormap(parula)
colorbar
title('A szz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

sxzA = squeeze(sigma(1,3,:,:));
figure(3)
contourf(X,Z,sxzA);
colormap(parula)
colorbar
title('A sxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off


sigma = hatStressSurf(uhat2,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);

sxxN = squeeze(sigma(1,1,:,:));
figure(4)
contourf(X,Z,sxxN);
colormap(parula)
colorbar
title('N sxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

szzN = squeeze(sigma(3,3,:,:));
figure(5)
contourf(X,Z,szzN);
colormap(parula)
colorbar
title('N szz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

sxzN = squeeze(sigma(1,3,:,:));
figure(6)
contourf(X,Z,sxzN);
colormap(parula)
colorbar
title('N sxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off


x = linspace(0,dx,gridSize);
z = linspace(0,dz,gridSize);
b = sqrt(3)/2;
[X,Z] = meshgrid(x,z);
[txx,tzz,txz] = imageStressAnalyticEdgePerp(1, b, 0.28, X, Z, x1, z1);
% [txx,tzz,txz] = imageStressAnalyticEdgePar(1, b, 0.28 , X, Z, x1, z1);
figure(7)
txx = txx;%./norm(txx);
meantxx = mean(txx,'all');
stdtxx = std(txx,0,'all');
displace = 5*stdtxx;
limits = [meantxx - displace, meantxx + displace];
contourf(X,Z,txx);
colormap(parula)
colorbar
caxis(limits)
title('b perp sxx')
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figure(8)
tzz = tzz;%./norm(tzz);
meantzz = mean(tzz,'all');
stdtzz = std(tzz,0,'all');
displace = 5*stdtzz;
limits = [meantzz - displace, meantzz + displace];
contourf(X,Z,tzz);
colormap(parula)
colorbar
caxis(limits)
title('b perp szz')
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

figure(9)
txz = txz;%./norm(txz);
meantxz = mean(txz,'all');
stdtxz = std(txz,0,'all');
displace = 5*stdtxz;
limits = [meantxz - displace, meantxz + displace];
contourf(X,Z,txz);
colormap(parula)
colorbar
caxis(limits)
title('b perp sxz')
xlabel('b')
ylabel('b')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

%%
segments=constructsegmentlist(rn,links);
FPStress = gridFieldPointStress(segments, a, MU, NU, X,Y,Z);
sxx = squeeze(FPStress(1,:,:));
figure(1)
contourf(X,Z,sxx);
colormap(parula)
colorbar
title('A fpxx')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

szz = squeeze(FPStress(3,:,:));
figure(2)
contourf(X,Z,szz);
colormap(parula)
colorbar
title('A fpzz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off

sxz = squeeze(FPStress(6,:,:));
figure(3)
contourf(X,Z,sxz);
colormap(parula)
colorbar
title('A fpxz')
xlabel('x')
ylabel('z')
hold on
plot(x1,z1,'.','color','black','MarkerSize',10)
hold off


x = linspace(0,dx,gridSize);
z = linspace(0,dz,gridSize);
b = sqrt(3)/2;
[X,Z] = meshgrid(x,z);
[txx,tzz,txz] = DlnStressAnalyticEdgePerp(1, b, 0.28, X, Z, x1, z1);
figure(7)
txx = txx;%./norm(txx);
meantxx = mean(txx,'all');
stdtxx = std(txx,0,'all');
displace = 5*stdtxx;
limits = [meantxx - displace, meantxx + displace];
contourf(X,Z,txx,30);
colormap(parula)
colorbar
caxis(limits)
title('b perp sxx')
xlabel('b')
ylabel('b')

figure(8)
tzz = tzz;%./norm(tzz);
meantzz = mean(tzz,'all');
stdtzz = std(tzz,0,'all');
displace = 5*stdtzz;
limits = [meantzz - displace, meantzz + displace];
contourf(X,Z,tzz,30);
colormap(parula)
colorbar
caxis(limits)
title('b perp szz')
xlabel('b')
ylabel('b')

figure(9)
txz = txz;%./norm(txz);
meantxz = mean(txz,'all');
stdtxz = std(txz,0,'all');
displace = 5*stdtxz;
limits = [meantxz - displace, meantxz + displace];
contourf(X,Z,txz,30);
colormap(parula)
colorbar
caxis(limits)
title('b perp sxz')
xlabel('b')
ylabel('b')

function [txx, tyy, txy] = imageStressAnalyticEdgePerp(E, b, nu, x, y, a, c)
%%%
% Stress on point (x, y) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (a, c).
% b perpendicular to surface.
%%%
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


function [txx, tyy, txy] = imageStressAnalyticEdgePar(E, b, nu, x, y, a, c)
%%%
% Stress on point (a, c) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (x, y).
% p parallel to surface.
%%%
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
  2.*a .* ymc .* ((3.*x+a) .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 - ymc2.*ymc2) ./ den3;
txx = D.*txx;
   
tyy =       ...%xma .* (xma2 + 3.*ymc2) ./ den1 + ...
           -xma .* (xpa2 + 3.*ymc2) ./ den2 + ...
  -2.*a .* (xma .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 + ymc2.*ymc2) ./ den3;
tyy = D.*tyy;

txy =        ...%ymc .* (xma2 - ymc2) ./ den1 + ...
            -ymc .* (xpa2 - ymc2) ./ den2 + ...
  4.*a.*x .* ymc .* (3.*xpa2 - ymc2) ./ den3;
txy = D.*txy;
end


function [txx, tyy, txy] = DlnStressAnalyticEdgePerp(E, b, nu, x, y, a, c)
%%%
% Stress on point (x, y) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (a, c).
% b perpendicular to surface.
%%%
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

txx =       -ymc .* (3.*xma2 + ymc2)./den1;% + ...
%             ymc .* (3.*xpa2 + ymc2)./den2 + ...
%   4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
txx = D.*txx;

tyy =     ymc .* (xma2 - ymc2)./den1;% + ...
%          -ymc .* (xpa2 - ymc2)./den2 + ...
%   4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
tyy = D.*tyy;

txy =       xma .* (xma2 - ymc2)./den1;% + ...
%            -xpa .* (xpa2 - ymc2)./den2 + ...
%   2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
txy = D.*txy;
end


% [sxx2, syy, szz2, sxy, sxz2, syz] = hatStress(uhat2,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);
% figure(7)
% contourf(X,Z,sxx2);
% colormap(parula)
% colorbar
% title('N sxx2')
% xlabel('x')
% ylabel('z')

% figure(8)
% contourf(X,Z,syy);
% colormap(parula)
% colorbar
% title('N syy')
% xlabel('x')
% ylabel('z')

% figure(9)
% contourf(X,Z,szz2);
% colormap(parula)
% colorbar
% title('N szz2')
% xlabel('x')
% ylabel('z')
% 
% figure(10)
% contourf(X,Z,sxy);
% colormap(parula)
% colorbar
% title('N sxy')
% xlabel('x')
% ylabel('z')

% figure(11)
% contourf(X,Z,sxz2);
% colormap(parula)
% colorbar
% title('N sxz2')
% xlabel('x')
% ylabel('z')
% save('./mat_files/stresses_bperp_surf', 'uhat', 'uhat2', 'fend', 'fend2', 'sxx', 'szz', 'sxz', 'sxx2', 'szz2', 'sxz2')

% figure(12)
% contourf(X,Z,syz);
% colormap(parula)
% colorbar
% title('N syz')
% xlabel('x')
% ylabel('z')

% gridSize = 40;
% 
% x = linspace(0,5e2,gridSize);
% z = linspace(0,5e2,gridSize);
% b = sqrt(3)/2;
% [X,Z] = meshgrid(x,z);
% [txx,tzz,txz] = imageStressAnalyticEdgePerp(145, b, 0.28, X, Z, 1e2, 2.5e2);
% figure(13)
% txx = txx;%./norm(txx);
% meantxx = mean(txx,'all');
% stdtxx = std(txx,0,'all');
% displace = 5*stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X,Z,txx);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b perp sxx')
% xlabel('b')
% ylabel('b')
% 
% figure(14)
% tzz = tzz;%./norm(tzz);
% meantzz = mean(tzz,'all');
% stdtzz = std(tzz,0,'all');
% displace = 5*stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X,Z,tzz);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b perp szz')
% xlabel('b')
% ylabel('b')
% 
% figure(15)
% txz = txz;%./norm(txz);
% meantxz = mean(txz,'all');
% stdtxz = std(txz,0,'all');
% displace = 5*stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X,Z,txz);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b perp sxz')
% xlabel('b')
% ylabel('b')
% 
% [txx,tzz,txz] = imageStressAnalyticEdgePar(145, b, 0.28, X, Z, 1e2, 2.5e2);
% figure(16)
% txx = txx;%./norm(txx);
% meantxx = mean(txx,'all');
% stdtxx = std(txx,0,'all');
% displace = 5*stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X,Z,txx);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b par sxx')
% xlabel('b')
% ylabel('b')
% 
% figure(17)
% tzz = tzz;%./norm(tzz);
% meantzz = mean(tzz,'all');
% stdtzz = std(tzz,0,'all');
% displace = 5*stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X,Z,tzz);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b par szz')
% xlabel('b')
% ylabel('b')
% 
% figure(18)
% txz = txz;%./norm(txz);
% meantxz = mean(txz,'all');
% stdtxz = std(txz,0,'all');
% displace = 5*stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X,Z,txz);
% colormap(parula)
% colorbar
% caxis(limits)
% title('b par sxz')
% xlabel('b')
% ylabel('b')

% sigma2 = hatStressSurf(uhat,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);
% surf(X*amag,Y*amag,mumag*PKresolved1','EdgeColor','none'); 

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

function sigma = gridFieldPointStress(segments, a, mu, nu, X,Y,Z)

gridSize = size(X);

sigma = zeros(6, gridSize(1),gridSize(2));

x1 = [segments(:,6), segments(:,7), segments(:,8)];
x2 = [segments(:,9), segments(:,10), segments(:,11)];

b = [segments(:,3), segments(:,4), segments(:,5)];

x0 = zeros(1,3);
for col = 1:gridSize(2)
    for row = 1:gridSize(1)
        x0(1,1) = X(row,col);
        x0(1,2) = Y(row,col);
        x0(1,3) = Z(row,col);
        sigma(:,row,col)=FieldPointStress(x0,x1,x2,b,a,mu,nu);
    end
end


end

% function [txx, tyy, txy] = imageStressAnalyticEdgePerp(E, b, nu, x, y, a, c)
% %%%
% % Stress on point (x, y) induced by edge dislocation parallel to the 
% % surface at x = 0. Dislocation coordinates are (a, c).
% % b perpendicular to surface.
% %%%
% D = E.*b./(4.*pi.*(1-nu.^2));
% ymc = y-c;
% ymc2 = ymc.^2;
% xma = x-a;
% xma2 = xma.^2;
% xpa = x+a;
% xpa2 = xpa.^2;
% den1 = (xma2 + ymc2).^2;
% den2 = (xpa2 + ymc2).^2;
% den3 = den2 .* (xpa2 + ymc2);
% 
% txx =       -ymc .* (3.*xma2 + ymc2)./den1 + ...
%              ymc .* (3.*xpa2 + ymc2)./den2 + ...
%   4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
% txx = D.*txx;
% 
% tyy =     ymc .* (xma2 - ymc2)./den1 + ...
%          -ymc .* (xpa2 - ymc2)./den2 + ...
%   4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
% tyy = D.*tyy;
% 
% txy =       xma .* (xma2 - ymc2)./den1 + ...
%            -xpa .* (xpa2 - ymc2)./den2 + ...
%   2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
% txy = D.*txy;
% end

function [txx, tyy, txy] = StressAnalyticEdgePar(E, b, nu, x, y, a, c)
%%%
% Stress on point (a, c) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (x, y).
% p parallel to surface.
%%%
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

txx =     xma .* (xma2 - ymc2) ./ den1 + ...
         -xpa .* (xpa2 - ymc2) ./ den2 + ...
  2.*a .* ymc .* ((3.*x+a) .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 - ymc2.*ymc2) ./ den3;
txx = D.*txx;
   
tyy =       xma .* (xma2 + 3.*ymc2) ./ den1 + ...
           -xma .* (xpa2 + 3.*ymc2) ./ den2 + ...
  -2.*a .* (xma .* xpa .* xpa2 - 6.*x.*xpa.*ymc2 + ymc2.*ymc2) ./ den3;
tyy = D.*tyy;

txy =        ymc .* (xma2 - ymc2) ./ den1 + ...
            -ymc .* (xpa2 - ymc2) ./ den2 + ...
  4.*a.*x .* ymc .* (3.*xpa2 - ymc2) ./ den3;
txy = D.*txy;
end

function f=pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments)
%nodal force on dislocation segments due to applied stress sigext
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,nx,ny,nz)
[nseg,~]=size(segments);
f=zeros(nseg,3);

b=segments(:,3:5);
r0=segments(:,6:8);
r1=segments(:,9:11);
r01=r1-r0;
rmid = 0.5*(r0+r1);
for i =1:nseg % evaluate sigmahat at midpoint of each segment (do something better later)
    xi = rmid(i,1:3);

    if any(isnan(xi))
        disp('Segment midpoint is undefined! See segforcevec.m');
        pause;
    end

    sigext = hatStressStatic(uhat,nc,xnodes,D,mx,mz,w,h,d,xi); %must not pass virtsegs!

    sigb=sigext*b(i,1:3)';
    l = r01(i,:);
    f(i,:)=cross(sigb',l);
end
end

