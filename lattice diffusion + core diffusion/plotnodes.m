function plotnodes(rn,links,plim)
%plot nodes
%only those nodes within [-plim,plim] in x,y,z directions are plotted
%hold on;

for i=1:size(links,1)
    n0 = links(i,1);
    n1 = links(i,2);
    
    if n0~=0 && n1~=0 && max(max(abs(rn([n0,n1],:))))<=plim
        %filter out "infinity" lines
        plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'r','linewidth',2.0);
        hold on
    end 
end
% 
% hold off
axis equal
grid

%set(h,'EraseMode','xor','MarkerSize',18)
%xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);
%view(viewangle);
%saveas(h,strcat('G:\DDlab\DDLab\image\',num2str(curstep)),'bmp');
%drawnow
%pause(0.01);