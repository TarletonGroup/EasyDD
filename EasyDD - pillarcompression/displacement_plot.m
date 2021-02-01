function [output]=displacement_plot(rn,links,NU,xnodes,dx,dy,dz,mx,my,mz,Sleft,Sright,Stop,Sbot,Sfront,Sback)

Sall=[Sleft;Sright;Stop;Sbot;Sfront;Sback];
gnl=Sall(:,1);
[Ux, Uy, Uz] = Utilda_bb3(rn,links,gnl,NU,xnodes,dx,dy,dz,mx,my,mz);
surfnodes=xnodes(gnl,:);
quiver3(surfnodes(:,1),surfnodes(:,2),surfnodes(:,3),Ux,Uy,Uz)
plotHandle=gcf;
axis equal

output=plotHandle;