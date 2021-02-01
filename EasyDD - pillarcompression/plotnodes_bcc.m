function [output]=plotnodes_bcc(rn,links,plim,vertices)
%plot nodes
%only those nodes within [-plim,plim] in x,y,z directions are plotted
figure(1);
clf
amag = 3.18e-4; %lattice vector BCC W, in microns
% amag=1;
plot3(0,0,0); hold on;
LINKMAX=length(links(:,1));
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    %to skip external nodes...
    if rn(n0,end)==67||rn(n1,end)==67
       continue;
    end
    lvec = amag*(rn(n1,1:3)-rn(n0,1:3));
    plane_n = links(i,6:8);
    bvec = links(i,3:5);
    normb=bvec/norm(bvec);
    plane_n = cross(lvec/norm(lvec),bvec/norm(bvec));
    check1pos=norm(normb-[1 1 1]*0.5);
    check1neg=norm(normb+[1 1 1]*0.5);
    check2pos=norm(normb-[-1 1 1]*0.5);
    check2neg=norm(normb+[-1 1 1]*0.5);
    check3pos=norm(normb-[1 -1 1]*0.5);
    check3neg=norm(normb+[1 -1 1]*0.5);
    check4pos=norm(normb-[1 1 -1]*0.5);
    check4neg=norm(normb+[1 1 -1]*0.5);
%     if plane_n == [-1 0 1]/sqrt(2)
    r0 = rn(n0,1:3)*amag;
%      if rn(n0,4)==0
        %filter out "infinity" lines
        if check1pos<eps || check1neg<eps
       p1=plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'r','LineWidth',2,'DisplayName','[1 1 1]');
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'r','LineWidth',1);
        elseif check2pos<eps || check2neg<eps
            p2=plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'g','LineWidth',2,'DisplayName','[-1 1 1]');
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'g','LineWidth',1);
        elseif check3pos<eps || check3neg<eps
            p3=plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'b','LineWidth',2,'DisplayName','[1 -1 1]');
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'b','LineWidth',1);
        elseif check4pos<eps || check4neg<eps
            p4=plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'y','LineWidth',2,'DisplayName','[1 1 -1]');
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'y','LineWidth',1);
        else
            p5=plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'m','LineWidth',2,'DisplayName','Other');
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'m','LineWidth',1);
        end
        %plot3(rn(n0,1)*bmag,rn(n0,2)*bmag,rn(n0,3)*bmag,'k.');
        %plot3(rn(n1,1)*bmag,rn(n1,2)*bmag,rn(n1,3)*bmag,'k.');
%      end

end

%plot film film location
%dt = delaunayTriangulation(vertices*bmag);
% dt = delaunayTri(vertices);
%[tri, Xb] = freeBoundary(dt);
%trisurf(tri, Xb(:,1), Xb(:,2), Xb(:,3), 'FaceColor', 'white','FaceAlpha', 0.1);
%plot bounding box
face1=[1 2 4 3 1];
face2=[5 6 8 7 5];
vertices_scaled=vertices*amag;
surf1=vertices_scaled(face1,:);
surf2=vertices_scaled(face2,:);

plot3(surf1(:,1),surf1(:,2),surf1(:,3),'k','LineWidth',2);
hold on;
plot3(surf2(:,1),surf2(:,2),surf2(:,3),'k','LineWidth',2);

side = vertices_scaled([1 5],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([2 6],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([3 7],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([4 8],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);

% plot virtual segments
% if isempty(virtual_seg)
%     %do nothing
% else
% for j=1:size(virtual_seg,1)% (node_ID_int, node_ID_inf, bx,by,bz,x_int,y_int,z_int,x_inf,y_inf,z_inf,nx,ny,nz)
% plot3([virtual_seg(j,6) virtual_seg(j,9)] , [virtual_seg(j,7) virtual_seg(j,10)], [virtual_seg(j,8) virtual_seg(j,11)],'k--');
% end

plotHandle=gcf;
hold off
axis equal
% grid
xlabel('x-direction (\mu m)','FontSize',10);
ylabel('y-direction (\mu m)','FontSize',10);
zlabel('z-direction (\mu m)','FontSize',10);

 xlim([0 vertices(2,1)]*amag);
 ylim([0 vertices(3,2)]*amag);
 zlim([0 vertices(5,3)]*amag);
 legend([p1 p2 p3 p4 p5],{'[1 1 1]','[-1 1 1]','[1 -1 1]','[1 1 -1]','Other'})

%axis equal;
output=plotHandle;
