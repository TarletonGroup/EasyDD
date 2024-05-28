clear 
load restart.50sources.start.mat

for i=228:-1:211
    
    finish=0;
    num=num2str(i,'%03d');
    fullFileName=strcat(num,'.mat');
    
    while finish~=1
        if exist(fullFileName,'file')
        load(fullFileName);
        plotHandle=plotnodes(rn,links,plim,vertices);
        saveas(plotHandle,num,'epsc');
        finish=1;
        else
            pause(60);
        end
    end
end
