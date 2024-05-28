%%



% % code below is just to visulise the defomed mesh

% clear all
% close all
% disp('Loading restart file');
% load 53201.mat %restart.50sources.75744.mat

plotmesh=0;
plotstress=1;

%refine mesh
mx=mx*1;

%update results
disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
     Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
     w,h,d,mx,my,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D(dx,dy,dz,mz,MU,NU);    
disp('Done!');

[TriangleCentroids,TriangleNormals,tri,Xb] = ...
     MeshSurfaceTriangulation(xnodes,Stop,Sbot,Sfront,Sback,Sleft,Sright);

[uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
                    Stop,gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,dy,dz,simTime0,mx,my,mz,unfixedDofs,Kred,Lred,Ured,utilda_0);

segments = constructsegmentlist(rn,links);
utilda = zeros(3*mno,1);
gn = 1:mno; % global node number
x0 = xnodes(gn,1:3); % field point
point_array_length = size(x0,1);
segments_array_length = size(segments,1);
%Full C version (including combinatorials, checks, corrections etc.)

disp('Calculating displacements');
[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);                       
%[Uxf, Uyf, Uzf] =displacement_fivel(x0,segments,NU); %gives same answer

utilda(3*gn -2) = Ux;
utilda(3*gn -1) = Uy;
utilda(3*gn   ) = Uz;

if plotmesh
    disp('Plotting mesh');
    xp = zeros(mno,3);
    for j =1:mno;
         xp(j,1:3) = xnodes(j,1:3)...
             + 1e2*utilda(3*j-2:3*j)'; %+ 0e4*uhat(3*j-2:3*j)';
    end
%     amag=1;
%     xp = amag*xp;
    figure;clf;hold on;view(0,0)
    xlabel('x');ylabel('y');zlabel('z')
    style='-k';
    for p =1:mel
    %      plot3(xp(:,1),xp(:,2),xp(:,3),'.') % plot nodes
        % plot elements
        plot3(xp(nc(p,[1:4,1]),1),xp(nc(p,[1:4,1]),2),xp(nc(p,[1:4,1]),3),style) 
        plot3(xp(nc(p,[5:8,5]),1),xp(nc(p,[5:8,5]),2),xp(nc(p,[5:8,5]),3),style) 
        plot3(xp(nc(p,[1,5]),1),xp(nc(p,[1,5]),2),xp(nc(p,[1,5]),3),style) % 
        plot3(xp(nc(p,[2,6]),1),xp(nc(p,[2,6]),2),xp(nc(p,[2,6]),3),style) % 
        plot3(xp(nc(p,[3,7]),1),xp(nc(p,[3,7]),2),xp(nc(p,[3,7]),3),style) % 
        plot3(xp(nc(p,[4,8]),1),xp(nc(p,[4,8]),2),xp(nc(p,[4,8]),3),style) % 
    end
    axis equal
    zlabel('z-direction (\mum)');
    xlabel('x-direction (\mum)');
%     zlim([-6 6])
%     xlim([0 9e4])
    title('$\tilde{u}$ scaled','FontSize',14,'Interpreter','Latex');
end

%-------------------generate stress--------------------------
if plotstress
    X = linspace(0,dx,2*mx)';
    Z = linspace(0,dx,2*mx)';
    Y = 0.5*dy; %middle slide

    X_size = length(X);
    Y_size = length(Y);
    Z_size = length(Z);
    
    RSSu = zeros(X_size,Y_size);
    RSS = zeros(X_size,Y_size);
    
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
    bur=[1 1 1]/sqrt(3);
    normal=[0 -1 1]/sqrt(2);
    cosphi=dot(normal,[0 0 1])/norm(normal);
    coslamda=dot(bur,[0 0 1])/norm(bur);

    
    disp('Calculating stresses');
    for i= 1:X_size;
        for j = 1:Z_size
            x0 = [X(i) Y Z(j)]; % field point
            sigmahat = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
            
            
%             rssu=sigmahat*bur';
%             RSSu(i,j)=norm(rssu-rssu'*normal'*normal');
%             RSSu(i,j)=norm(rss);
            
            
            Sxxu(i,j) = sigmahat(1,1);
            Syyu(i,j) = sigmahat(2,2);
            Szzu(i,j) = sigmahat(3,3);

            Sxyu(i,j) = sigmahat(1,2); %isotropic
            Sxzu(i,j) = sigmahat(1,3); %isotropic
            Syzu(i,j) = sigmahat(2,3); %isotropic
            
            RSSu(i,j)=Szzu(i,j)*cosphi*coslamda;

            sigmatilde=FieldPointStress(x0,x1,x2,b,a,MU,NU);
            
%             sigmatilde1=[sigmatilde(1) sigmatilde(4) sigmatilde(5); sigmatilde(4) sigmatilde(2) sigmatilde(6); sigmatilde(5) sigmatilde(6) sigmatilde(3) ];
%             rss=sigmatilde1*bur';

            
            Sxx(i,j) = sigmatilde(1);
            Syy(i,j) = sigmatilde(2);
            Szz(i,j) = sigmatilde(3);
            Sxy(i,j) = sigmatilde(4);
            Syz(i,j) = sigmatilde(5);
            Sxz(i,j) = sigmatilde(6);
            
            RSS(i,j)=Szz(i,j)*cosphi*coslamda;

        end
    end
    
    figure; clf
    surf(X*amag,Z*amag,mumag*(RSS)','EdgeColor','none');
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    box on
    set(gca,'LineWidth',4,'FontSize',16)
    xlabel('x-direction (\mum)','FontSize',18)
    ylabel('z-direction (\mum)','FontSize',18)
    title('Resolved shear stress','Interpreter','Latex','FontSize',22);
    grid off
    h=colorbar;
    xlabel(h,'MPa','FontSize',18);
    clim([-400, 400]);
    h.FontSize = 16;
    h.LineWidth = 2;
    strtime = ['Time = ', num2str(simTime),'s' ];
    text(-1.5,4,strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',18,'Color','b');

    figure; clf
    surf(X*amag,Z*amag,mumag*(Sxzu+Sxz)','EdgeColor','none'); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    box on
    set(gca,'LineWidth',4,'FontSize',16)
    xlabel('x-direction (\mum)','FontSize',18)
    ylabel('z-direction (\mum)','FontSize',18)
    title('$$\hat{\sigma}_{xz}$$+$$\tilde{\sigma}_{xz}$$','Interpreter','Latex','FontSize',22);
    grid off
    h=colorbar;
    xlabel(h,'MPa','FontSize',18);
    clim([-800, 1200]);
    h.FontSize = 16;
    h.LineWidth = 2;
    strtime = ['Time = ', num2str(simTime),'s' ];
    text(-1.5,4,strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',18,'Color','b');
        
%     figure; clf
%     subplot(2,1,1)
%     surf(X*amag,Z*amag,mumag*Sxz','EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)')
%     ylabel('z-direction (\mum)')
%     title('$$\tilde{\sigma}_{xz}$$','Interpreter','Latex');
%     grid off
%     h=colorbar;
%     xlabel(h,'MPa');
% 
%     subplot(2,1,2)
%     surf(X*amag,Z*amag,mumag*Sxxu','EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)')
%     ylabel('z-direction (\mum)')
%     title('$$\hat{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%         h=colorbar;
%     xlabel(h,'MPa');
%     %saveas(gcf,'sxx','epsc')
end
