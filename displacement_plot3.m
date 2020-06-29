%load('C:\Users\Edmund.OUMS-TARLETONHP\Dropbox\Bruce\DDLabFEM\Images\Displacement Code Simulations\Paper Figure\demorun2ET.mat')

mx2 = mx*4;
my2 = my*4;
mz2 = mz*4;

[xsn,snc]=gensurfnodemesh(mx2,my2,mz2,dx,dy,dz);
gnl=1:size(xsn,1);
gnl=gnl';
mel=size(snc,1);

tic
%[Ux, Uy, Uz] = Utilda_bb3(rn,links,gnl,NU,xsn,dx,dy,dz,mx2,my2,mz2);
[Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gnl,NU,xsn,dx,dy,dz,mx2,my2,mz2);
toc
utilda(:,1)=Ux;
utilda(:,2)=Uy;
utilda(:,3)=Uz;

scalefactor=100;

xp=xsn+utilda*scalefactor;

swap=[snc(:,4),snc(:,3)];
snc(:,3:4)=swap;

%%
h=figure(2);
% clf
% set(axes,'position',[0 0 1 1])
utildaNorm = sqrt(utilda(:,1).^2 + utilda(:,2).^2 + utilda(:,3).^2);
col = utildaNorm;
f=patch('Faces',snc,'Vertices',xp,'FaceVertexCData',col,'FaceColor','flat','FaceAlpha',1, 'LineWidth',0.01);
view(3)
axis image
axis off
colorbar
title('$\tilde{u}/b$', 'interpreter','latex')

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print('utildaLowRes','-dpdf','-r600')
print('utilda','-dpng','-r600')
%%
% figure(2);
% clf;
% hold on;
% style='-k';
% for p =1:mel
%     plot3(xp(snc(p,[1:4,1]),1),xp(snc(p,[1:4,1]),2),xp(snc(p,[1:4,1]),3),style)
% end
plotHandle=gcf;
% axis equal
%%
output=plotHandle;