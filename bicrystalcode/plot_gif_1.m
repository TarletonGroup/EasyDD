
% %mmmm1: the number of figures
num_fig=292;
for j=1:1:num_fig
    A=imread(sprintf('%d.bmp',j-1));
    [I,map]=rgb2ind(A,256);
    %??gif????
    if(j==1)
        imwrite(I,map,'01movefig_400MPa.gif','DelayTime',0.2,'LoopCount',1)
    else
        imwrite(I,map,'01movefig_400MPa.gif','WriteMode','append','DelayTime',0.2)    
    end
    

end
    

% plot .avi
writerObj = VideoWriter('01_bicrystal_crossslip.avi'); % output the name +.avi
writerObj.FrameRate = 50; % ???25fps
writerObj.Quality = 90;   % ???????90%
open(writerObj);
for i = 11:1:num_fig % ?100???
    img = imread(sprintf('%d.bmp',i-1)); %????????img????1.bmp 2.bmp ...
    writeVideo(writerObj, img);
end
close(writerObj);