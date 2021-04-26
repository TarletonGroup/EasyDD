%===============================================================%
% Daniel Hortelano Roig (11/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk 

% Saves a figure of the initial input nodal network.
%===============================================================%
%%
rnfig = figure('Name',INPUTNAME,'visible','off');
grid on, hold on
for n = 1:size(rn,1)
    plot3(rn(n,1),rn(n,2),rn(n,3),'.','MarkerSize',5,'Color','blue')
end
ax = gca;
ax.XLim = [0 dx];
ax.YLim = [0 dy];
ax.ZLim = [0 dz];
set(rnfig, 'visible', 'on');
savefilefig = append('figures/',INPUTNAME,'.fig');
savefig(rnfig,savefilefig,'compact');
close(rnfig);
% To open saved figure: openfig(['figures/' INPUTNAME '.fig'])