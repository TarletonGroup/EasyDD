
%mmmm1: the number of figures

for j=1:1512
    A=imread(sprintf('%d.bmp',j-1));
    [I,map]=rgb2ind(A,256);
    %??gif????
    if(j==1)
        imwrite(I,map,'01movefig_glide_climb.gif','DelayTime',0.2,'LoopCount',1)
    else
        imwrite(I,map,'01movefig_glide_climb.gif','WriteMode','append','DelayTime',0.2)    
    end
    

end
    

% plot .avi
writerObj = VideoWriter('3loops_bulk+core_100+60_1300K.avi'); % output the name +.avi
writerObj.FrameRate = 25; % ???25fps
writerObj.Quality = 90;   % ???????90%
open(writerObj);
for i = 1:2:1512 % ?100???
    img = imread(sprintf('%d.bmp',i-1)); %????????img????1.bmp 2.bmp ...
    writeVideo(writerObj, img);
end
close(writerObj);