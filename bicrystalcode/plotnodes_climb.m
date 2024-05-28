function plotnodes_climb(connectivity, rn,links,plim,tottime,dt)
%plot nodes
%only those nodes within [-plim,plim] in x,y,z directions are plotted

eps=4E-2;
% eps=7.6154e-07;
plot3(0,0,0); hold on;


LINKMAX=length(links(:,1));

colo=zeros(length(rn(:,1)),1);

for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    if(n0~=0)&& (n1~=0) 
        lvec = (rn(n1,1:3)-rn(n0,1:3));
        bvec = links(i,3:5);
%         color = 'r';
        dline = norm(lvec);
        linedir=lvec./dline;
        costh2=(linedir*bvec')^2/(bvec*bvec');
        sinth2 = 1-costh2;
        if sinth2 < eps%HY20180503: +-5 degrees
            color = 'r';%HY20180504: plot pure screw components with color red
            plane = links(i,6:8);
%             pause
        else
            plane=cross(lvec,bvec);
            plane=plane/norm(plane);
        end
        %filter out "infinity" lines
%         plot3(rn(n0,1),rn(n0,2),rn(n0,3),'bo','LineWidth',0.5);

%%%%--------------plot the climb nodes
        if abs(rn(n0,3))> 5 
            
             colo(n1)=1;
%             plot3(rn(n0,1),rn(n0,2),rn(n0,3),'bo','LineWidth',1);
%         elseif rn(n0,3)< -5 
%             plot3(rn(n0,1),rn(n0,2),rn(n0,3),'r.','LineWidth',1);
        end
        
        if abs(rn(n1,3))> 5
            %             plot3(rn(n1,1),rn(n1,2),rn(n1,3),'bo','LineWidth',1);
            colo(n1)=1;
            
            %         elseif rn(n1,3)< -5
            %             plot3(rn(n1,1),rn(n1,2),rn(n1,3),'r.','LineWidth',1);
        end
        
        
%         plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'k-','LineWidth',2);
        %HY20180728: only used for plotting LC fcc locks
        bref1 = [  -1   1    1   ]/2;
        bref2 = [  1   -1    1   ]/2;
        bref3 = [  1    1   -1   ]/2;
        
%         plane = links(i,6:8);

        %
 %%%%%%   plot seg with different normal planes       
        if norm(plane-1/sqrt(2)*[1 -1 0])<eps || norm(plane+1/sqrt(2)*[1 -1 0])<eps %HY20191006:(-1-11)
            color = [0 1 0]; %green
        elseif norm(plane-1/sqrt(2)*[1  0 -1])<eps || norm(plane+1/sqrt(2)*[1  0 -1])<eps 
            color = [1 0.2 0.9];%HY20191006:(111) %pink
        elseif norm(plane-1/sqrt(2)*[0  -1 1])<eps || norm(plane+1/sqrt(2)*[0  -1 1])<eps 
            color = [0.5 0 0.2];
        elseif norm(plane-1/sqrt(2)*[0  1 1])<eps || norm(plane+1/sqrt(2)*[0  1 1])<eps 
            color = [0 0.6 0.7]; %green
        elseif norm(plane-1/sqrt(2)*[1 0 1])<eps || norm(plane+1/sqrt(2)*[1 0 1])<eps
            color = [0.4 0 1];
        elseif norm(plane-1/sqrt(2)*[1 1 0])<eps || norm(plane+1/sqrt(2)*[1 1 0])<eps
            color = [0 0 1]; %blue
        else
            color=[0 0 0];
        end

%%%%%  plot seg with different motion mode_glide or climb
%         if colo(n0)==1 || colo(n1)==1  %HY20191006:(-1-11)
%             color = 'r';
% %         elseif norm(plane-1/sqrt(2)*[1 1 0])<eps || norm(plane+1/sqrt(2)*[1 1 0])<eps 
% %             color = 'b';%HY20191006:(111)
% %         elseif norm(plane-1/sqrt(2)*[ 0 1 -1])<eps || norm(plane+1/sqrt(2)*[0 1 -1])<eps 
% %             color = 'r';
%         else
%             color = 'k';            
%         end
%%%%%  plot junction seg
        if connectivity(n0)>=3 || connectivity(n1)>=3  %HY20191006:(-1-11)
            color = 'r';
            plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'color',color,'LineWidth',2);
        else
%             color = 'k';    
             plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'color',color,'LineWidth',1);
        end


%         if sinth2 < sind(5)^2%HY20180503: +-5 degrees
%             color = 'r';%HY20180504: plot pure screw components with color red
%             %             pause
%         end
%         plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'color',color,'LineWidth',1);
%         plot3(rn(n0,1),rn(n0,2),rn(n0,3),'g.','markersize',12);
%         if n0==4
%             plot3(rn(n0,1),rn(n0,2),rn(n0,3),'ro','markersize',12);
%         end
%         if n0==2
%             plot3(rn(n0,1),rn(n0,2),rn(n0,3),'mo','markersize',12);
%         end
    end 
end

% hold on 
% axis equal
% % grid
% xlabel('x-direction (\mu m)','FontSize',10);
% ylabel('y-direction (\mu m)','FontSize',10);
% zlabel('z-direction (\mu m)','FontSize',10);
% 
%  xlim([0 vertices(2,1)]*bur);
%  ylim([0 vertices(3,2)]*bur);
%  zlim([0 vertices(5,3)]*bur);
% plot3([-plim, plim]',[0, 0]',[0, 0]','r--','LineWidth',3);

%  time = 1E3*tottime*(1/mumag);
% % time=tottime;
% 
% strtime = ['Time = ', num2str(time),'s' ];
% % strtime1= [' dt  = ', num2str(dt),'s'];
% % strtime2={strtime;strtime1};
% text(-0.3*plim(1),-0.5*plim(2),1.4*plim(3),strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'Color','b');
% 


% text(-1*plim,0.7*plim,strtime1,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'Color','b');


% r = 150; % pixels per inch
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1 1]);
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1600 900]/r);

set(gcf,'Position',[100 100 1300 700])
hold off
axis equal

%grid
% set(gca,'fontsize',20)
set(gca,'LineWidth',2,'FontSize',16)
xlabel('x / b','FontSize',18);
ylabel('y / b','FontSize',18);
zlabel('z / b','FontSize',18);
xlim([0 plim(1)]); ylim([0 plim(2)]); zlim([0 plim(3)]);