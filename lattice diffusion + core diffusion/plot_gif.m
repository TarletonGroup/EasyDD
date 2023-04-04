load('.\output\input_area.dat')
x=input_area(:,1);
y=input_area(:,2);
yy=max(y);
% j=0;
% for i=1:10:length(yy)
%     j=j+1;
%     x(j)=xx(i);
%     y(j)=yy(i);
%     
% end
 plot(x,y);
    
%  plot(x,y,'LineWidth',3);
% xlim([0 inf]);
% ylim([0 yy*1.1]);

xlabel('Time(s)');
ylabel('Enclosed area (b^2)');

%???
i=0;
%?????????????'xdata',??????'ydata',??????30
h = line('xdata',[],'ydata',[],'color','r','marker','.','markersize',30);
for ii=1:1000:length(y)
    %??????
    set(h,'xdata',x(ii),'ydata',y(ii));
    %?????
    figure();
    %??
    drawnow
    %??????????
    pause(0.02)
    hold on
    i=i+1;
    %????????:/image????i.bmp
     print(1,'-dbmp',sprintf('%d',i))
    %     print(i,'-dbmp')
    %??figure()
    close;
end

for j=1:i
    %??????
    A=imread(sprintf('%d.bmp',j));
    [I,map]=rgb2ind(A,256);
    %??gif????
    if(j==1)
        imwrite(I,map,'movefig_ellipse1.gif','DelayTime',0.1,'LoopCount',Inf)
    else
        imwrite(I,map,'movefig_ellipse1.gif','WriteMode','append','DelayTime',0.1)    
    end
end

