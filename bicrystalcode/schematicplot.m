function [output]=schematicplot(rn,links,vertices,U_bar,Fend,amag,dx,totalSimTime)
figure(1);
clf
%amag=1;
Udot =100*1E3*dx*(1E-4/160E9)*100;
xmax=Udot*totalSimTime;
ymax=1000*xmax;
subplot(2,3,[1 2]);
flatplot(rn,links,vertices,amag,1,3);
subplot(2,3,3);
flatplot(rn,links,vertices,amag,2,3);
subplot(2,3,[4 5]);
flatplot(rn,links,vertices,amag,1,2);
ax4=subplot(2,3,6);
plot(U_bar,-Fend);
xlim(ax4,[0 xmax]);
ylim(ax4,[0 ymax]);

plotHandle=gcf;
output=plotHandle;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [graph]=flatplot(rn,links,vertices,amag,xind,yind)
hold on;
LINKMAX=length(links(:,1));
for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    %to skip external nodes...
    if rn(n0,end)==67||rn(n1,end)==67
       continue;
    end
    lvec = amag*(rn(n1,[xind,yind])-rn(n0,[xind,yind]));
    r0 = rn(n0,[xind,yind])*amag;
    if isequal(links(i,3:5),[-1 -1 -1]/sqrt(3))
    plot(rn([n0,n1],xind)*amag,rn([n0,n1],yind)*amag,'r','LineWidth',2);  
        quiver(r0(1),r0(2),lvec(1),lvec(2),'r','LineWidth',1); 
    elseif isequal(links(i,3:5),[-1 1 1]/sqrt(3))
    plot(rn([n0,n1],xind)*amag,rn([n0,n1],yind)*amag,'y','LineWidth',2);  
        quiver(r0(1),r0(2),lvec(1),lvec(2),'y','LineWidth',1); 
    else
        plot(rn([n0,n1],xind)*amag,rn([n0,n1],yind)*amag,'b','LineWidth',2);  
        quiver(r0(1),r0(2),lvec(1),lvec(2),'b','LineWidth',1);
    end
end

face1=[1 2 4 3 1];
face2=[5 6 8 7 5];
vertices_scaled=vertices*amag;
surf1=vertices_scaled(face1,:);
surf2=vertices_scaled(face2,:);
plot(surf1(:,xind),surf1(:,yind),'k','LineWidth',2);
hold on;
plot(surf2(:,xind),surf2(:,yind),'k','LineWidth',2);

side = vertices_scaled([1 5],[xind yind]);
plot(side(:,1),side(:,2),'k','LineWidth',2);
side = vertices_scaled([2 6],[xind yind]);
plot(side(:,1),side(:,2),'k','LineWidth',2);
side = vertices_scaled([3 7],[xind yind]);
plot(side(:,1),side(:,2),'k','LineWidth',2);
side = vertices_scaled([4 8],[xind yind]);
plot(side(:,1),side(:,2),'k','LineWidth',2);

plothandle=gcf;
hold off
axis equal
% grid
if xind==1
    xlabel('x-direction (\mu m)','FontSize',10);
elseif xind==2
    xlabel('y-direction (\mu m)','FontSize',10);
elseif xind==3
    xlabel('z-direction (\mu m)','FontSize',10);
end
    
if yind==1
    ylabel('x-direction (\mu m)','FontSize',10);
elseif yind==2
    ylabel('y-direction (\mu m)','FontSize',10);
elseif yind==3
    ylabel('z-direction (\mu m)','FontSize',10);
end

if xind==1
    xind2=2;
elseif xind==2
    xind2=3;
elseif xind==3
    xind2=5;
end
if yind==1
    yind2=2;
elseif yind==2
    yind2=3;
elseif yind==3
    yind2=5;
end

 xlim([0 vertices(xind2,xind)]*amag);
 ylim([0 vertices(yind2,yind)]*amag);

%axis equal;
graph=plothandle;
end
