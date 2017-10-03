%load a test condition
run inputCheck.m

segments = constructsegmentlist(rn,links);

%construct finite element arrays
%construct stiffeness matrix K and pre-compute L,U decompositions.
disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);    

X = linspace(0,dx,2*mx)';
Z = linspace(0,dy,2*my)';
Y = 0.5*dy; 

X_size = length(X);
Y_size = length(Y);
Z_size = length(Z);

Sxxuhat = zeros(X_size,Z_size);
Syyuhat = zeros(X_size,Z_size);
Szzuhat = zeros(X_size,Z_size);
Sxyuhat = zeros(X_size,Z_size);
Sxzuhat = zeros(X_size,Z_size);
Syzuhat = zeros(X_size,Z_size);

Sxxutilde = zeros(X_size,Z_size);
Syyutilde = zeros(X_size,Z_size);
Szzutilde = zeros(X_size,Z_size);
Sxyutilde = zeros(X_size,Z_size);
Sxzutilde = zeros(X_size,Z_size);
Syzutilde = zeros(X_size,Z_size);

Sxx = zeros(X_size,Z_size);
Syy = zeros(X_size,Z_size);
Szz = zeros(X_size,Z_size);
Sxy = zeros(X_size,Z_size);
Sxz = zeros(X_size,Z_size);
Syz = zeros(X_size,Z_size);

%%
%---------------------------------------------------------------
% evaluate tractions on gammat and solve for uhat  to recover
% dislocation stress field 

gn = 1:mno; %gammau(:,1);  % global node number
x0 = xnodes(gn,1:3); % field point
point_array_length = size(x0,1);
segments_array_length = size(segments,1);

[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);

utilda(3*gn - 2) = Ux;
utilda(3*gn - 1) = Uy;
utilda(3*gn   ) = Uz;

%% ---------------------------------------------------------------
for i= 1:X_size;
    
    for j = 1:Z_size
        
        x0 = [X(i) Y Z(j)]; % field point
        
        %sigmahatFEM = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
        %Sxxuhat(i,j) = sigmahatFEM(1,1);
        %Syyuhat(i,j) = sigmahatFEM(2,2);
        %Szzuhat(i,j) = sigmahatFEM(3,3);
        %Sxyuhat(i,j) = sigmahatFEM(1,2); %isotropic
        %Sxzuhat(i,j) = sigmahatFEM(1,3); %isotropic
        %Syzuhat(i,j) = sigmahatFEM(2,3); %isotropic
        
        sigmatildeFEM = hatStress(utilda,nc,xnodes,D,mx,mz,w,h,d,x0);
        Sxxutilde(i,j) = sigmatildeFEM(1,1);
        Syyutilde(i,j) = sigmatildeFEM(2,2);
        Szzutilde(i,j) = sigmatildeFEM(3,3);
        Sxyutilde(i,j) = sigmatildeFEM(1,2); %isotropic
        Sxzutilde(i,j) = sigmatildeFEM(1,3); %isotropic
        Syzutilde(i,j) = sigmatildeFEM(2,3); %isotropic
        
        x1=segments(:,6:8);
        x2=segments(:,9:11);
        b=segments(:,3:5);
        sigmatilde=FieldPointStress(x0,x1,x2,b,a,MU,NU);
        Sxx(i,j) = sigmatilde(1);
        Syy(i,j) = sigmatilde(2);
        Szz(i,j) = sigmatilde(3);
        Sxy(i,j) = sigmatilde(4);
        Syz(i,j) = sigmatilde(5);
        Sxz(i,j) = sigmatilde(6);           
    end
end

%%
subplot(5,2,1)
contourf(X,Z,Sxx'); 
axis([0 dx/2 0 dz])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xx}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);
subplot(5,2,2)
contourf(X,Z,Sxxutilde'); 
axis([0 dx/2 0 dz])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xx}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,3)
contourf(X,Z,Syy'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{yy}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,4)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Syyutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{yy}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,5)
contourf(X,Z,Szz'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{zz}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,6)
%c=1E-4*linspace(-1,1,100);
contourf(X,Z,Szzutilde'); 
%surf(X,Z,Szzutilde','EdgeColor','none'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{zz}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,7)
contourf(X,Z,Sxy'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xy}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,8)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Sxyutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xy}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,9)
contourf(X,Z,Sxz'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xz}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,10)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Sxzutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xz}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);