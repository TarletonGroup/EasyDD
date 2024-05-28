
hold on
figure
if do_visualisestress==1
    visualisedislocationstress(SimBox,curstep,rnnew,links,a,MU,NU);
end

PlotBox_b(SimBox);
plotnodes_climb(connectivity,rnnew,links,SimBox,simTime,dt)

%         real_simTime=simTime*1e3/mumag/1e6; %s
%         time=real_simTime;
% strtime = ['Time = ', num2str(simTime),'s' ];
% text(-1*plim,0.5*plim, 2.0*plim,strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'Color','b');
if do_inclusion
    plot_boundary()
end
if do_twosource
    plot_boundary()
end
view([0 -1 0])


hold on
n=1;
burgg = links(n,3:5);
normplane = links(n,6:8);

m=links(n,2);
quiver3(rn(m,1),rn(m,2),rn(m,3),burgg(1),burgg(2),burgg(3),2e3,'-b','LineWidth',2); % burgers
quiver3(rn(m,1),rn(m,2),rn(m,3),normplane(1),normplane(2),normplane(3),2e3,'-r','LineWidth',2);