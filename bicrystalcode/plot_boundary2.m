
% r1=50; r2=50;
% centerx=[75 75]/0.2473;
% centery=[-110 110]/0.2473;
% centerz=[0 0]/0.2473;
% r=[r1/0.2473 r2/0.2473];
% color=['k' 'k'];
% alpha=[0.2 0.2];
% [x, y, z]=sphere(20);
% for k=1:2
%     surf(r(k)*x+centerx(k),r(k)*y+centery(k),r(k)*z+centerz(k),'FaceColor', color(k), 'LineStyle', 'none', 'FaceAlpha',alpha(k))
%     hold on
% end
% hold on
function plot_boundary()

global SimBox
dim=1;%in x direction
loc_boundary=SimBox(1,dim)/2;% location of the boundary is x=dx/2

[y,z]=meshgrid(0:1000:SimBox(2),0:2000:SimBox(3));
x=ones(13).*loc_boundary;
hold on
surf(x,y,z,'FaceColor', 'y', 'EdgeColor','y');
hold off
end