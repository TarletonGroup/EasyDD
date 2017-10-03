 load plim.mat

for i=2041
num=num2str(i,'%03d');
%load(strcat(num,'.mat'));

%load(strcat(num,'-100sources.mat'));
load(strcat(num,'-100sources.mat'));

plotHandle = plotnodes(rn,links,plim,vertices);
saveas(plotHandle,strcat(num,'-50sources'),'epsc');
end
