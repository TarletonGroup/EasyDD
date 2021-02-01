clear all
close all
screen_size = get(0, 'ScreenSize');

figure(21)
clf
figure(22)
clf

load('withH-8sources-withOB-R5E-7-U0.25.mat')
% close(1)


visualisedislocationstress
Xplot = X;
Xplotsize = size(Xplot,1);
Zplot = Z;
Zplotsize = size(Zplot,1);
Sxxplot = Sxx;
SxxHplot = SxxH;
mx = mx/20;

% HY20181201: interpolation
vars = {'X','Y','Z'};
clear(vars{:})
visualisehatstress

Xsize = size(X,1);
Zsize = size(Z,1);
if (size(Zplot,1)/size(Z,1)~=size(Xplot,1)/size(X,1))
    display('mesh number not exactly scaled')
%     pause
end

%HY20181201: make square matrice
Xresize = X;
Z2 = (vertices(8,3)+(vertices(8,1)-vertices(8,3))/(Xsize-Zsize)):((vertices(8,1)-vertices(8,3))/(Xsize-Zsize)):vertices(8,1);
Zresize = [Z;Z2'];
Sxxuresize = zeros(Xsize,Xsize);
Sxxuresize(1:Xsize,1:Zsize) = Sxxu;
Xplotresize = Xplot;
Zplot2 = (vertices(8,3)+(vertices(8,1)-vertices(8,3))/(Xplotsize-Zplotsize)):((vertices(8,1)-vertices(8,3))/(Xplotsize-Zplotsize)):vertices(8,1);
Zplotresize = [Zplot;Zplot2'];

[Xmesh,Zmesh] = meshgrid(Xresize,Zresize);
[Xplotmesh,Zplotmesh] = meshgrid(Xplotresize,Zplotresize);

Sxxumesh = interp2(Xmesh,Zmesh,Sxxuresize,Xplotmesh,Zplotmesh,'linear');
Sxxuplot = Sxxumesh(1:Xplotsize,1:Zplotsize);

% close(2,3,22);

screen_size = get(0, 'ScreenSize');
figure(21)
clf
figure(22)
clf
% visualise
f21 = figure(21)
hold on
    subplot(2,1,1)
    surf(Xplot*amag,Zplot*amag,mumag*(Sxxplot+SxxHplot)','EdgeColor','none'); 
%     surf(Xplot*amag,Zplot*amag,mumag*(SxxHplot)','EdgeColor','none');
%     surf(X*amag,Z*amag,mumag*(Sxx)','EdgeColor','none'); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)','FontSize',24)
    ylabel('z-direction (\mum)','FontSize',24)
    title('$$\tilde{\sigma}_{xx}+\bar{\sigma}_{xx}$$','Interpreter','Latex');
%     title('$$\bar{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    hh=colorbar;
    caxis([-3000,3000])
%     xlim([0 5])
    xlabel(hh,'MPa','FontSize',24);
    ax = gca;
    ax.FontSize = 20;

    subplot(2,1,2)
    surf(Xplot*amag,Zplot*amag,mumag*(Sxxuplot+Sxxplot+SxxHplot)','EdgeColor','none'); 
%     surf(Xplot*amag,Zplot*amag,mumag*(SxxHplot)','EdgeColor','none'); 
    %HY20180125
%     isovalues = (-500:250:1500);
%     contour(X*amag,Z*amag,mumag*(Sxxu+Sxx)',isovalues); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)','FontSize',24)
    ylabel('z-direction (\mum)','FontSize',24)
    title('$$\hat{\sigma}_{xx}+\tilde{\sigma}_{xx}+\bar{\sigma}_{xx}$$','Interpreter','Latex');
%     title('$$\bar{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    hh=colorbar;
    caxis([-3000,3000])
%     xlim([0 5])
    xlabel(hh,'MPa','FontSize',24); 
    ax = gca;
    ax.FontSize = 20;
    
%     subplot(3,1,3)
%     surf(Xplot*amag,Zplot*amag,mumag*(Sxxuplot+Sxxplot+SxxHplot)','EdgeColor','none'); 
%     %HY20180125
% %     isovalues = (-500:250:1500);
% %     contour(X*amag,Z*amag,mumag*(Sxxu+Sxx)',isovalues); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)','FontSize',24)
%     ylabel('z-direction (\mum)','FontSize',24)
%     title('$$\hat{\sigma}_{xx}+\tilde{\sigma}_{xx}+\bar{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     hh=colorbar;
%     caxis([-2000,2000])
%     xlim([0 5])
% %     zlim([2 2.5])
%     xlabel(hh,'MPa','FontSize',24); 
%     ax = gca;
%     ax.FontSize = 20;

set(f21,'Units','Inches');
pos = get(f21,'Position');
set(f21,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print('stress-withH-withOB','-dpdf','-r1200')
print('stress-withH-withOB-Rev-refined','-dtiff','-r1200')


f22 = figure(22)
hold on
    subplot(1,2,1)
    surf(Xplot*amag,Zplot*amag,mumag*(Sxxplot+SxxHplot)','EdgeColor','none'); 
%     surf(Xplot*amag,Zplot*amag,mumag*(SxxHplot)','EdgeColor','none');
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)','FontSize',24)
    ylabel('z-direction (\mum)','FontSize',24)
    title('$$\tilde{\sigma}_{xx}+\bar{\sigma}_{xx}$$','Interpreter','Latex');
%     title('$$\bar{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    hh=colorbar;
    caxis([-3000,3000])
    xlim([1 3])
%     xlabel(hh,'MPa','FontSize',24);
    ax = gca;
    ax.FontSize = 20;

    subplot(1,2,2)
    surf(Xplot*amag,Zplot*amag,mumag*(Sxxuplot+Sxxplot+SxxHplot)','EdgeColor','none'); 
%     surf(Xplot*amag,Zplot*amag,mumag*(SxxHplot)','EdgeColor','none'); 
    %HY20180125
%     isovalues = (-500:250:1500);
%     contour(X*amag,Z*amag,mumag*(Sxxu+Sxx)',isovalues); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)','FontSize',24)
    ylabel('z-direction (\mum)','FontSize',24)
    title('$$\hat{\sigma}_{xx}+\tilde{\sigma}_{xx}+\bar{\sigma}_{xx}$$','Interpreter','Latex');
%     title('$$\bar{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    hh=colorbar;
    caxis([-3000,3000])
    xlim([1 3])
    xlabel(hh,'MPa','FontSize',24); 
    ax = gca;
    ax.FontSize = 20;
% 
set(f22,'Units','Inches');
pos = get(f22,'Position');
set(f22,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print('stress-withH-withOB-zoom','-dpdf','-r1200')
print('stress-withH-withOB-zoom-Rev-refined','-dtiff','-r1200')













% clear all
% screen_size = get(0, 'ScreenSize');
% 
% load('withH-8sources-withOB-R5E-7-U0.2.mat')
% % close(1)
% 
% visualise
% 
% f22 = figure(22)
% hold on
%     subplot(2,1,1)
%     surf(X*amag,Z*amag,mumag*Sxx','EdgeColor','none'); 
%     view(2)
%     axis equal;
%     axis([0 dx*amag 0 dz*amag])
%     xlabel('x-direction (\mum)','FontSize',24)
%     ylabel('z-direction (\mum)','FontSize',24)
%     title('$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
%     grid off
%     hh=colorbar;
%     caxis([-200,200])
%     xlabel(hh,'MPa','FontSize',24);
%     ax = gca;
%     ax.FontSize = 20;
% 
%     subplot(2,1,2)
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
%     caxis([-2000,2000])
%     xlabel(hh,'MPa','FontSize',24); 
%     ax = gca;
%     ax.FontSize = 20;
%    
% set(f22,'Units','Inches');
% pos = get(f22,'Position');
% set(f22,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print('stress-withH-withOB','-dpdf','-r1200')