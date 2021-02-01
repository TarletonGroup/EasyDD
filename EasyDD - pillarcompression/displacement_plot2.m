function [output]=displacement_plot2(mx2,my2,mz2,rn,links,NU,dx,dy,dz,mx,my,mz)

[xsn,snc]=gensurfnodemesh(mx2,my2,mz2,dx,dy,dz);
gnl=1:size(xsn,1);
gnl=gnl';
mel=size(snc,1);

[Ux, Uy, Uz] = Utilda_bb3(rn,links,gnl,NU,xsn,dx,dy,dz,mx,my,mz);
utilda(:,1)=Ux;
utilda(:,2)=Uy;
utilda(:,3)=Uz;
scalefactor=10;
xp=surfnodes+utilda*scalefactor;

swap=[snc(:,4),snc(:,3)];
snc(:,3:4)=swap;

figure;clf;hold on;
style='-k';
for p =1:mel
    plot3(xp(snc(p,[1:4,1]),1),xp(snc(p,[1:4,1]),2),xp(snc(p,[1:4,1]),3),style)
end
plotHandle=gcf;
axis equal

output=plotHandle;