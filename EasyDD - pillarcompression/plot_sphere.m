
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
function plot_sphere(RR,QQ,bur)

hold on

centerx=QQ(:,1);
centery=QQ(:,2);
centerz=QQ(:,3);
r=[RR ];
color=['k'];
alpha=[0.3];
[x, y, z]=sphere(20);
for k=1:size(RR,2)

    surf(r(k)*x+centerx(k),r(k)*y+centery(k),r(k)*z+centerz(k),'FaceColor', color, 'LineStyle', 'none', 'FaceAlpha',alpha)
    hold on
end
hold off
end