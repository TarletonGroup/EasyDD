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
dx=1/amag;
dy=5/amag;
dz=5/amag;

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

links = [1, 2, 0, 1, 1,  1, -1,  1];



plim=12/amag; %12microns
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,Smixed,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = STATIC_finiteElement3D(dx,dy,dz,mx,MU,NU,loading);
%%
planes = 0;
length = 100; % good behaviour 1000
range = linspace(0,dy,length);
lvec = range(2)-range(1);
scale = 1; % good behaviour 10
% Parameters with nice fend plots length = 1000, scale = 100
scenario = 2;
if scenario == 1
%     range = flip(range);
    gammat = Sleft;
    gammau = Sleft;
    gammaMixed = Sleft;
    midPoints = [0 - scale; 0.5*dy; 0.5*dz];%[0.5*dx; 0.5*dy; 0.5*dz];%
    planes = 1;    
elseif scenario == 2
    gammat = Sleft;
    gammau = Sleft;
    gammaMixed = Sleft;
    midPoints = [0 + lvec*scale; 0.5*dy; 0.5*dz];
    planes = 1;
elseif scenario == 3
    gammat = [Sright;Smixed];
    gammau = [Sright;Smixed];
    gammaMixed = [Sright;Smixed];
    midPoints = [dx - lvec*scale; 0.5*dy; 0.5*dz];%[0.5*dx; 0.5*dy; 0.5*dz];%
    planes = 2;
elseif scenario == 4
    gammat = [Sright;Smixed];
    gammau = [Sright;Smixed];
    gammaMixed = [Sright;Smixed];
    midPoints = [dx + lvec*scale; 0.5*dy; 0.5*dz];
    planes = 2;
end
gamma_dln = gammat(:,1);
fixedDofs =[3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1)];
freeDofs = [3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1)];
t = [0 1 0];
b = [1 0 0];
n = [0 0 1];
rn = zeros(length,3);
rn(:,1) = midPoints(1);
rn(:,2) = range;
rn(:,3) = midPoints(3);
links = zeros(length-1,8);
for i = 1:length-1
    links(i,:) = [i, i+1, b,  n];
end

 % Set surface node labels for surface node extraction.
        n_nodes = 4;
        surf_node_util = zeros(n_nodes+2, 6);
        xy = mx*my;
        xz = mx*mz;
        yz = my*mz;
        % For rectangular surface elements. Add an if statement and redifine
        % matrix to generalise it for triangular surface elements.
        surf_node_util(1:6, 1) = [5, 1, 8, 4, yz, 1]; % min(x), yz-plane, face 1 Sleft
%         surf_node_util(1:6, 1) = [1, 5, 4, 8, yz, 1]; % min(x), yz-plane, face 1 Sleft
        surf_node_util(1:6, 2) = [2, 6, 3, 7, yz, 1]; % max(x), yz-plane, face 2 Sright+SMixed
%         surf_node_util(1:6, 2) = [6, 2, 7, 3, yz, 1]; % max(x), yz-plane, face 2 Sright+SMixed
        surf_node_util(1:6, 3) = [6, 5, 7, 8, xz, 2]; % min(y), xz-plane, face 4 Sfront
        surf_node_util(1:6, 4) = [1, 2, 4, 3, xz, 2]; % max(y), xz-plane, face 3 Sback
        surf_node_util(1:6, 5) = [5, 6, 1, 2, xy, 3]; % min(z), xy-plane, face 5 Sbot
        surf_node_util(1:6, 6) = [4, 3, 8, 7, xy, 3]; % max(z), xy-plane, face 6 Stop

f_hat = zeros(3*mno, 1);
[x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
                                               planes, 4, surf_node_util);
[f_dln_node, f_dln_se,...
 f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
tolerance = dx/10^6;

figure(1)
hold on
plot3(x3x6(1:3:end,1), x3x6(2:3:end,1), x3x6(3:3:end,1),'ko')
plot3(x3x6(1:3:end,2), x3x6(2:3:end,2), x3x6(3:3:end,2),'ko')
plot3(x3x6(1:3:end,3), x3x6(2:3:end,3), x3x6(3:3:end,3),'ko')
plot3(x3x6(1:3:end,4), x3x6(2:3:end,4), x3x6(3:3:end,4),'ko')
plot3(xnodes(Sleft(:,1),1),xnodes(Sleft(:,1),2),xnodes(Sleft(:,1),3),'b.')
plot3(xnodes(Sright(:,1),1),xnodes(Sright(:,1),2),xnodes(Sright(:,1),3),'b.')
plot3(xnodes(Smixed(:,1),1),xnodes(Smixed(:,1),2),xnodes(Smixed(:,1),3),'b.')
plot3(rn(:,1),rn(:,2),rn(:,3),'r.')
hold off

[uhat,fend,Ubar] = STATIC_analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
                       0, 0, gammaMixed,fixedDofs,freeDofs,dx,simTime ,...
                       gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
                       f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance);

% [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
%                        gamma_disp, gammat, gammaMixed,fixedDofs,freeDofs,dx,simTime ,...
%                        gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f_dln_node,...
%                        f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, tolerance, Ubar, dt);

segments=constructsegmentlist(rn,links);
image_force = pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);
%%
gridSize = 3;
x = linspace(0, 2*lvec*scale, gridSize);
y = linspace(0.5*dy, 0.5*dy, gridSize);
z = linspace(0.5*dz-lvec*scale, 0.5*dz+lvec*scale, gridSize);
[X,Z] = meshgrid(x,z);
Y = meshgrid(y);

% r0=segments(50,6:8);
% r1=segments(55,9:11);
% r01=r1-r0;
% rmid = 0.5*(r0+r1);
% xi = rmid(1:3);
% X = xi(1);
% Y = xi(2);
% Z = xi(3);

[sxx, syy, szz, sxy, sxz, syz] = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,X,Y,Z);

% xi(1),xi(2),xi(3)); test that it works


% sigma2 = hatStressStatic(uhat,nc,xnodes,D,mx,mz,w,h,d,[X,Y,Z]);
% % xi); test that it works
% difference = max(abs(sigma2 - [sxx; syy; szz; sxy; sxz; syz]))

% pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);
function [sxx, syy, szz, sxy, sxz, syz] = hatStress(uhat,nc,x,D,mx,mz,w,h,d,X,Y,Z)

gridSize = size(X,1);
i = ceil(X/w);
j = ceil(Y/h);
k = ceil(Z/d);

i0 = find(i==0);
j0 = find(j==0);
k0 = find(k==0);
if any(i0)
    i(i0) = 1;
end
if any(j0)
    i(j0) = 1;
end
if any(k0)
    i(k0) = 1;
end
p = i + (k-1)*mx + (j-1)*mx*mz;


xc = [0.5*(x(nc(p,1),1) + x(nc(p,7),1)),... % midpoint x
      0.5*(x(nc(p,1),2) + x(nc(p,7),2)),... % midpoint y
      0.5*(x(nc(p,1),3) + x(nc(p,7),3))];   % midpoint z

Xc = reshape(xc(:,1), gridSize, gridSize);
Yc = reshape(xc(:,2), gridSize, gridSize);
Zc = reshape(xc(:,3), gridSize, gridSize);

a = 0.5*w;
b = 0.5*h;
c = 0.5*d;

s1 = (X-Xc)/a;
s2 = (Y-Yc)/b;
s3 = (Z-Zc)/c;

ds1dx=1/a;
ds2dy=1/b;
ds3dz=1/c;

pm1 =[-1  1  1 -1 -1  1  1 -1];
pm2 =[ 1  1  1  1 -1 -1 -1 -1];
pm3 =[-1 -1  1  1 -1 -1  1  1];

dNds1 = zeros(gridSize, gridSize, 8);
dNds2 = zeros(gridSize, gridSize, 8);
dNds3 = zeros(gridSize, gridSize, 8);
B = zeros(gridSize, gridSize, 6, 24);
U = zeros(gridSize, gridSize, 24);
sigma = zeros(6, gridSize, gridSize);
for i = 1:gridSize
    for j = 1:gridSize
        for a = 1:8
            dNds1(i,j,a) = 1/8*pm1(a)*(1+pm2(a)*s2(i,j))*(1+pm3(a)*s3(i,j));
            dNds2(i,j,a) = 1/8*(1+pm1(a)*s1(i,j))*pm2(a)*(1+pm3(a)*s3(i,j));
            dNds3(i,j,a) = 1/8*(1+pm1(a)*s1(i,j))*(1+pm2(a)*s2(i,j))*pm3(a);

            B(i,j,1,3*(a-1)+1) = dNds1(i,j,a)*ds1dx;
            B(i,j,2,3*(a-1)+2) = dNds2(i,j,a)*ds2dy;
            B(i,j,3,3*(a-1)+3) = dNds3(i,j,a)*ds3dz;

            B(i,j,4,3*(a-1)+1) = B(i,j,2,3*(a-1)+2);
            B(i,j,4,3*(a-1)+2) = B(i,j,1,3*(a-1)+1);

            B(i,j,5,3*(a-1)+1) = B(i,j,3,3*a);
            B(i,j,5,3*a)   = B(i,j,1,3*(a-1)+1);

            B(i,j,6,3*(a-1)+2) = B(i,j,3,3*a);
            B(i,j,6,3*a)   = B(i,j,2,3*(a-1)+2);
            
            U(i,j, 3*a-2) = uhat(3*nc(p(i,j),a)-2);
            U(i,j, 3*a-1) = uhat(3*nc(p(i,j),a)-1);
            U(i,j, 3*a  ) = uhat(3*nc(p(i,j),a)  );
        end
    end
end

for i = 1:gridSize % Traverse first dimension
    for j = 1:gridSize % Traverse second dimension
        Utmp = squeeze(U(i,j,:));
        Btmp = squeeze(B(i,j,:,:));
        sigma(:,i,j) = D*(Btmp*Utmp);
    end
end

sxx = squeeze(sigma(1,:,:));
syy = squeeze(sigma(2,:,:));
szz = squeeze(sigma(3,:,:));
sxy = squeeze(sigma(4,:,:));
sxz = squeeze(sigma(5,:,:));
syz = squeeze(sigma(6,:,:));

% sigmaA(:,:)=D*(B(:,:)*U(:));
% 
% sigma=zeros(3);
% sigma(1,1)=sigmaA(1); %11
% sigma(2,2)=sigmaA(2); %22
% sigma(3,3)=sigmaA(3); % 33
% sigma(1,2)=sigmaA(4);% 12
% sigma(1,3)=sigmaA(5);% 13
% sigma(2,3)=sigmaA(6);% 23
% sigma(2,1)=sigma(1,2);
% sigma(3,1)=sigma(1,3);
% sigma(3,2)=sigma(2,3);

end

function sigma=hatStressStatic(uhat,nc,x,D,mx,mz,w,h,d,x0)
                              
% returns stress tensor at a point x=x0 by constructing B matrix B(x)
% notes with mx=90 dx=6,dy=dz=1, agrees with S(1:5) agree with Abaqus to 1% 
% S(6) is approx zero so the error is huge!

i=ceil(x0(1)/w);
i = max(i,1);
j=ceil(x0(2)/h);
j =max(j,1);
k=ceil(x0(3)/d);
k=max(k,1);
p = i + (k-1)*mx + (j-1)*mx*mz;
% if any(x0<0) || p > mx*mz*mz || isnan(p) 
%     %disp('Node outside domain! See hatStress.m')
%     %disp(strcat('x0 = [',num2str(x0),']')); 
%     %Usually stemming from NaN in utilda!
%     %pause;
%     sigma = zeros(3); % do something better like remeshing later...
%     %Sort of have to do this when using int_trapezium.m because real
%     %segment will move out of domain during trial before it can be
%     %remeshed?
%     return;
% end

% 4. ----- .3
%  �\       �\
%  � \      � \
% 1. -\---- .2 \ 
%   \  \     \  \
%    \ 8. ----\- .7
%     \ �      \ �
%      \�       \�
%      5. ----- .6


%  s2   s3    
%    \  ^ 
%     \ �      
%      \�       
%       . -----> s1
% redefine local corordinate system (s1,s2,s3) to have same orientation as
% the global system (x,y,z) this should save calcualting Jij and inv(J)
xc = zeros(3,1);
for i=1:3
    xc(i)= 0.5*(x(nc(p,1),i)+x(nc(p,7),i)); %xc is center of element p
end

a = w/2; % element dx
b = h/2; % element dy
c = d/2; % element dz

s1 =  (x0(1) - xc(1))/a;
s2 = (x0(2) - xc(2))/b;
s3 = (x0(3) - xc(3))/c;

ds1dx=1/a;
ds2dy=1/b;
ds3dz=1/c;

pm1 =[-1  1  1 -1 -1  1  1 -1];
pm2 =[ 1  1  1  1 -1 -1 -1 -1];
pm3 =[-1 -1  1  1 -1 -1  1  1];
%eg shape function a: Na = 1/8*(1+pm1(a)*s1)(1+pm2(a)*s2)(1+pm3(a)*s3)
% dNa/ds1 = 1/8* pm1(a)*(1+pm2(a)*s2)(1+pm3(a)*s3)
% dNa/dx = (dNa/ds1)(ds1/dx) where ds1/dx = 1/a
% dN1/dx = 1/a(dNa/ds1) = -b*c*(1+s2)(1-s3)/(abc)

B=zeros(6,24);
for a= 1:8
    dNds1(a) = 1/8*pm1(a)*(1+pm2(a)*s2)*(1+pm3(a)*s3);
    dNds2(a) = 1/8*(1+pm1(a)*s1)*pm2(a)*(1+pm3(a)*s3);
    dNds3(a) = 1/8*(1+pm1(a)*s1)*(1+pm2(a)*s2)*pm3(a);
    
    B(1,3*(a-1)+1) = dNds1(a)*ds1dx;
    B(2,3*(a-1)+2) = dNds2(a)*ds2dy;
    B(3,3*(a-1)+3) = dNds3(a)*ds3dz;
    
    B(4,3*(a-1)+1) = B(2,3*(a-1)+2);
    B(4,3*(a-1)+2) = B(1,3*(a-1)+1);
    
    B(5,3*(a-1)+1) = B(3,3*a);
    B(5,3*a)   = B(1,3*(a-1)+1);
    
    B(6,3*(a-1)+2) = B(3,3*a);
    B(6,3*a)   = B(2,3*(a-1)+2);
end
%
%----------------------------------------------------------------------
U=zeros(24,1); sigmaA = zeros(6,1);
for a=1:8
    U(3*a-2)=uhat(3*nc(p,a)-2);
    U(3*a-1)=uhat(3*nc(p,a)-1);
    U(3*a)=uhat(3*nc(p,a));
end

sigmaA=D*(B*U); % 

sigma = sigmaA;

% sigma=zeros(3);
% sigma(1,1)=sigmaA(1); %11
% sigma(2,2)=sigmaA(2); %22
% sigma(3,3)=sigmaA(3); % 33
% sigma(1,2)=sigmaA(4);% 12
% sigma(1,3)=sigmaA(5);% 13
% sigma(2,3)=sigmaA(6);% 23
% sigma(2,1)=sigma(1,2);
% sigma(3,1)=sigma(1,3);
% sigma(3,2)=sigma(2,3);

end


% %----------------------------------------------------------------------
% U=zeros(24,1); sigmaA = zeros(6,1);
% for a=1:8
%     U(3*a-2)=uhat(3*nc(p,a)-2);
%     U(3*a-1)=uhat(3*nc(p,a)-1);
%     U(3*a)=uhat(3*nc(p,a));
% end
% 
% sigmaA=D*(B*U); % 
% 
% sigma=zeros(3);
% sigma(1,1)=sigmaA(1); %11
% sigma(2,2)=sigmaA(2); %22
% sigma(3,3)=sigmaA(3); % 33
% sigma(1,2)=sigmaA(4);% 12
% sigma(1,3)=sigmaA(5);% 13
% sigma(2,3)=sigmaA(6);% 23
% sigma(2,1)=sigma(1,2);
% sigma(3,1)=sigma(1,3);
% sigma(3,2)=sigma(2,3);
% 
% end


% x = linspace(0,5e2,1000);
% z = linspace(0,5e2,1000);
% [X,Z] = meshgrid(x,z);
% 
% [txx,tzz,txz] = imageStressAnalyticEdge1(200e3, 1, 0.28, X, Z, 1e2, 2.5e2);
% figure(1)
% txx = txx;%./norm(txx);
% meantxx = mean(txx,'all');
% stdtxx = std(txx,0,'all');
% displace = 5*stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X,Z,txx,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% xlabel('b')
% ylabel('b')
% 
% figure(2)
% tzz = tzz;%./norm(tzz);
% meantzz = mean(tzz,'all');
% stdtzz = std(tzz,0,'all');
% displace = 5*stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X,Z,tzz,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% xlabel('b')
% ylabel('b')
% 
% figure(3)
% txz = txz;%./norm(txz);
% meantxz = mean(txz,'all');
% stdtxz = std(txz,0,'all');
% displace = 5*stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X,Z,txz,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% xlabel('b')
% ylabel('b')
% 
% 
% %
% [txx,tzz,txz] = imageStressAnalyticEdge2(1e8, 0.001, 0.28, X, Z, 1e-3, 1e-3);
% figure(4)
% txx = txx./norm(txx);
% meantxx = mean(txx,'all');
% stdtxx = std(txx,0,'all');
% displace = 5*stdtxx;
% limits = [meantxx - displace, meantxx + displace];
% contourf(X,Z,txx,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% 
% figure(5)
% tzz = tzz./norm(tzz);
% meantzz = mean(tzz,'all');
% stdtzz = std(tzz,0,'all');
% displace = 5*stdtzz;
% limits = [meantzz - displace, meantzz + displace];
% contourf(X,Z,tzz,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% 
% figure(6)
% txz = txz./norm(txz);
% meantxz = mean(txz,'all');
% stdtxz = std(txz,0,'all');
% displace = 5*stdtxz;
% limits = [meantxz - displace, meantxz + displace];
% contourf(X,Z,txz,400);
% colormap(parula(400))
% colorbar
% caxis(limits)
% 
% %
% % xnodes2 = xnodes(gammau(:,1),1:3);
% % uhat2 = uhat(fixedDofs);
% pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);


function [txx, tyy, txy] = imageStressAnalyticEdge1(E, b, nu, x, y, a, c)
%%%
% Stress on point (x, y) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (a, c).
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

txx =       -ymc .* (3.*xma2 + ymc2)./den1 + ...
             ymc .* (3.*xpa2 + ymc2)./den2 + ...
  4.*a.*x .* ymc .* (3.*xpa2 - ymc2)./den3;
txx = D.*txx;

tyy =     ymc .* (xma2 - ymc2)./den1 + ...
         -ymc .* (xpa2 - ymc2)./den2 + ...
  4.*a .* ymc .* ((2.*a - x) .* xpa2 + (3.*x + 2.*a) .* ymc2)./den3;
tyy = D.*tyy;

txy =       xma .* (xma2 - ymc2)./den1 + ...
           -xpa .* (xpa2 - ymc2)./den2 + ...
  2.*a .* (-xma .* xpa .* xpa2 + 6.*x.*xpa.*ymc2 - ymc2.*ymc2)./den3;
txy = D.*txy;
end
function [txx, tyy, txy] = imageStressAnalyticEdge2(E, b, nu, x, y, a, c)
%%%
% Stress on point (a, c) induced by edge dislocation parallel to the 
% surface at x = 0. Dislocation coordinates are (x, y).
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

txx =     xma .* (xma2 - ymc2) ./ den1 + ... % + ymc2
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

