hold on;view(3)
xlabel('x');ylabel('y');zlabel('z')
for i=1:size(links,1)
    if isempty( gammaMiu_el(i,:))==0
        p1=gammaMiu_el(i,gammaMiu_el(i,:)~=0);
        %         plot3(xnodes(p,1),xnodes(p,2),xnodes(p,3),'r*')
        for n=1:size(p1,2)
            p=p1(n);
            plot3(xnodes(nc(p,[1:4,1]),1),xnodes(nc(p,[1:4,1]),2),xnodes(nc(p,[1:4,1]),3),'r-')
            plot3(xnodes(nc(p,[5:8,5]),1),xnodes(nc(p,[5:8,5]),2),xnodes(nc(p,[5:8,5]),3),'r-')
            plot3(xnodes(nc(p,[1,5]),1),xnodes(nc(p,[1,5]),2),xnodes(nc(p,[1,5]),3),'r-') %
            plot3(xnodes(nc(p,[2,6]),1),xnodes(nc(p,[2,6]),2),xnodes(nc(p,[2,6]),3),'r-') %
            plot3(xnodes(nc(p,[3,7]),1),xnodes(nc(p,[3,7]),2),xnodes(nc(p,[3,7]),3),'r-') %
            plot3(xnodes(nc(p,[4,8]),1),xnodes(nc(p,[4,8]),2),xnodes(nc(p,[4,8]),3),'r-') %
        end
    end
    
end
figure(2);clf;hold on;view(3)
xlabel('x');ylabel('y');zlabel('z')
plot3(xnodes(gammaU(:,1),1),xnodes(gammaU(:,1),2),xnodes(gammaU(:,1),3),'r*')
% plot3(xnodes(Sbot(:,1),1),xnodes(Sbot(:,1),2),xnodes(Sbot(:,1),3),'r*')
% plot3(xnodes(Sright(:,1),1),xnodes(Sright(:,1),2),xnodes(Sright(:,1),3),'b.')
% plot3(xnodes(Sleft(:,1),1),xnodes(Sleft(:,1),2),xnodes(Sleft(:,1),3),'b.')
% plot3(xnodes(Sfront(:,1),1),xnodes(Sfront(:,1),2),xnodes(Sfront(:,1),3),'k*')
% plot3(xnodes(Sback(:,1),1),xnodes(Sback(:,1),2),xnodes(Sback(:,1),3),'k*')
% plot3(xnodes(Smixed(:,1),1),xnodes(Smixed(:,1),2),xnodes(Smixed(:,1),3),'g*')