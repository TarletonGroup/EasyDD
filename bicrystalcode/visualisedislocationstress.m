% % code below is just to visulise the defomed mesh

% clear all
% close all
% % disp('Loading restart file');
% load matlabstep24319.mat

plotmesh=0;
plotstress=1;

%refine mesh
mx = 20*mx;

w = dx/mx; % elements width

my = round(mx*dy/dx); % # elements in y direction
my=max(my,1);
h=dy/my; % element height

mz = round(mx*dz/dx); % elements in z direction
mz=max(mz,1);
d = dz/mz; % element depth

mel = mx*my*mz;
mno=(mx+1)*(my+1)*(mz+1);

xnodes= zeros(mno,3); %global nodal coordinates gn

%update results
% disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
% % [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
% %      Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
% %      w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);    
% [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
%     Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
%     w,h,d,my,mz,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);
% disp('Done!');
% 
% [TriangleCentroids,TriangleNormals,tri,Xb] = ...
%      MeshSurfaceTriangulation(xnodes,Stop,Sbot,Sfront,Sback,Sleft,Sright);
% 
% % [uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
% %                     gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,simTime);
% [uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
%     gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,simTime,unfixedDofs,Kred,Lred,Ured);

segments = constructsegmentlist(rn,links);
utilda = zeros(3*mno,1);
gn = 1:mno; % global node number
x0 = xnodes(gn,1:3); % field point
point_array_length = size(x0,1);
segments_array_length = size(segments,1);
%Full C version (including combinatorials, checks, corrections etc.)

% if plotmesh
% disp('Calculating displacements');
% [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
%                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
%                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
%                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
%                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
%                        NU,point_array_length,segments_array_length);                       
% %[Uxf, Uyf, Uzf] =displacement_fivel(x0,segments,NU); %gives same answer
% 
% utilda(3*gn -2) = Ux;
% utilda(3*gn -1) = Uy;
% utilda(3*gn   ) = Uz;
% 
% 
%     disp('Plotting mesh');
%     xp = zeros(mno,3);
% %     for j =1:mno;
% %          xp(j,1:3) = xnodes(j,1:3)...
% %              + 5e4*utilda(3*j-2:3*j)'+ 0e4*uhat(3*j-2:3*j)';
% %     end
% %HY20190125
%     for j =1:mno;
%          xp(j,1:3) = xnodes(j,1:3)...
%              + 5E2 * utilda(3*j-2:3*j)'+ 1 * uhat(3*j-2:3*j)';
%     end
%     xp = xp*amag;
%     scale = 1;
%     xp = scale*xp;
%     figure;clf;hold on;view(0,0)
%     xlabel('x');ylabel('y');zlabel('z')
%     style='-k';
%     for p =1:mel
%     %      plot3(xp(:,1),xp(:,2),xp(:,3),'.') % plot nodes
%         % plot elements
%         plot3(xp(nc(p,[1:4,1]),1),xp(nc(p,[1:4,1]),2),xp(nc(p,[1:4,1]),3),style) 
%         plot3(xp(nc(p,[5:8,5]),1),xp(nc(p,[5:8,5]),2),xp(nc(p,[5:8,5]),3),style) 
%         plot3(xp(nc(p,[1,5]),1),xp(nc(p,[1,5]),2),xp(nc(p,[1,5]),3),style) % 
%         plot3(xp(nc(p,[2,6]),1),xp(nc(p,[2,6]),2),xp(nc(p,[2,6]),3),style) % 
%         plot3(xp(nc(p,[3,7]),1),xp(nc(p,[3,7]),2),xp(nc(p,[3,7]),3),style) % 
%         plot3(xp(nc(p,[4,8]),1),xp(nc(p,[4,8]),2),xp(nc(p,[4,8]),3),style) % 
%     end
%     axis equal
%     zlabel('z-direction (\mum)');
%     xlabel('x-direction (\mum)');
% %     zlim([-6 6])
% %     xlim([0 31])
%     title('$\tilde{u}\ scaled\ by\ 20\ +\ \hat{u}$','FontSize',14,'Interpreter','Latex');
% end

%-------------------generate stress--------------------------
if plotstress
    X = linspace(0,dx,mx)';
    Z = linspace(0,dy,my)';
    Y = 0.5*dy; %middle slide

    X_size = length(X);
    Y_size = length(Y);
    Z_size = length(Z);

    Sxxu = zeros(X_size,Y_size);
    Syyu = zeros(X_size,Y_size);
    Szzu = zeros(X_size,Y_size);
    Sxyu = zeros(X_size,Y_size);
    Sxzu = zeros(X_size,Y_size);
    Syzu = zeros(X_size,Y_size);
    Sxx = zeros(X_size,Y_size);
    Syy = zeros(X_size,Y_size);
    Szz = zeros(X_size,Y_size);
    Sxy = zeros(X_size,Y_size);
    Sxz = zeros(X_size,Y_size);
    Syz = zeros(X_size,Y_size);
    p1x = segments(:,6);
    p1y = segments(:,7);
    p1z = segments(:,8);
    p2x = segments(:,9);
    p2y = segments(:,10);
    p2z = segments(:,11);
    bx = segments(:,3);
    by = segments(:,4);
    bz = segments(:,5);
    x1 =  [ p1x,p1y,p1z];
    x2 = [p2x,p2y,p2z];
    b=[bx, by, bz];
    
    disp('Calculating stresses');
    for i= 1:X_size
        for j = 1:Z_size
            xx0 = [X(i) Y Z(j)]; % field point
            sigmahat = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,xx0);
            Sxxu(i,j) = sigmahat(1,1);
            Syyu(i,j) = sigmahat(2,2);
            Szzu(i,j) = sigmahat(3,3);

            Sxyu(i,j) = sigmahat(1,2); %isotropic
            Sxzu(i,j) = sigmahat(1,3); %isotropic
            Syzu(i,j) = sigmahat(2,3); %isotropic

            sigmatilde=FieldPointStress(xx0,x1,x2,b,a,MU,NU);
            Sxx(i,j) = sigmatilde(1);
            Syy(i,j) = sigmatilde(2);
            Szz(i,j) = sigmatilde(3);
            Sxy(i,j) = sigmatilde(4);
            Syz(i,j) = sigmatilde(5);
            Sxz(i,j) = sigmatilde(6);
            
%             sigmatildeH=HFieldPointStress(xx0,x1,x2,b,a,MU,NU);
%             SxxH(i,j) = sigmatildeH(1);
%             SyyH(i,j) = sigmatildeH(2);
%             SzzH(i,j) = sigmatildeH(3);
%             SxyH(i,j) = sigmatildeH(4);
%             SyzH(i,j) = sigmatildeH(5);
%             SxzH(i,j) = sigmatildeH(6);
            
            progress= [i,j]

        end
    end
    
%     figure; clf
%     subplot(3,1,1)
%     surf(X*amag,Z*amag,mumag*(Sxxu+Sxx)','EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)','FontSize',24)
%     ylabel('z-direction (\mum)','FontSize',24)
%     title('$$\hat{\sigma}_{xx}$$+$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     hh=colorbar;
% %     caxis([-100,100])
%     xlabel(hh,'MPa','FontSize',24);
%     ax = gca;
%     ax.FontSize = 20;
%     
%     subplot(3,1,2)
%     surf(X*amag,Z*amag,mumag*Sxx','EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)','FontSize',24)
%     ylabel('z-direction (\mum)','FontSize',24)
%     title('$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     hh=colorbar;
% %     caxis([-200,200])
%     xlabel(hh,'MPa','FontSize',24);
%     ax = gca;
%     ax.FontSize = 20;
% 
%     subplot(3,1,3)
%     surf(X*amag,Z*amag,mumag*Sxxu','EdgeColor','none'); 
%     %HY20180125
% %     isovalues = (-500:250:1500);
% %     contour(X*amag,Z*amag,mumag*(Sxxu+Sxx)',isovalues); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)','FontSize',24)
%     ylabel('z-direction (\mum)','FontSize',24)
%     title('$$\hat{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     hh=colorbar;
% %     caxis([-50,50])
%     xlabel(hh,'MPa','FontSize',24); 
%     ax = gca;
%     ax.FontSize = 20;
    %saveas(gcf,'sxx','epsc')
    
    %HY20180125: 3D stress contour
    disp('3D stress contour');
    Nnodes = size(x0,1);
    for i= 1:Nnodes
            x0 = xnodes(i,:); % field point
            sigmahat = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
            Sxxu(i) = sigmahat(1,1);
            Syyu(i) = sigmahat(2,2);
            Szzu(i) = sigmahat(3,3);

            Sxyu(i) = sigmahat(1,2); %isotropic
            Sxzu(i) = sigmahat(1,3); %isotropic
            Syzu(i) = sigmahat(2,3); %isotropic

            sigmatilde=FieldPointStress(x0,x1,x2,b,a,MU,NU);
            Sxx(i) = sigmatilde(1);
            Syy(i) = sigmatilde(2);
            Szz(i) = sigmatilde(3);
            Sxy(i) = sigmatilde(4);
            Syz(i) = sigmatilde(5);
            Sxz(i) = sigmatilde(6);

            i
    end
%     
%     figure; clf
%     surf(xnodes(:,1)*amag,xnodes(:,2)*amag,xnodes(:,3)*amag,mumag*(Sxxu+Sxx),'EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)')
%     ylabel('z-direction (\mum)')
%     title('$$\hat{\sigma}_{xx}$$+$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     h=colorbar;
%     xlabel(h,'MPa');
end
