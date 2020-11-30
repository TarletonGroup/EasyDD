function plotnodes(rn, links, FEM, viewangle)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Plots the dislocation network nodes.
    %===============================================================%

%% Extraction

vertices = FEM.vertices;

%% Plot nodes

%%% Setup figure

figname = "Dislocaton network";
figure('Name',figname);
grid on, hold on

%%% Plot links

plot3(0,0,0);

for i = 1:length(links(:,1))
    
    n0 = links(i,1);
    n1 = links(i,2);
    
    r0 = rn(n0,1:3);
    r1 = rn(n1,1:3);
    
    midseg = (r1 + r0)/2;
    
    %to skip external nodes...
    if rn(n0,end) == 67 || rn(n1,end) == 67
       continue;
    end
    
    lvec = rn(n1,1:3)-rn(n0,1:3);
    lunit = lvec / norm(lvec);
    bvec = links(i,3:5);
    bunit = bvec / norm(bvec);
    bplot = 2e2 * bvec;
    
    screworient = abs(dot(bunit,lunit));
    color = [1 screworient 0]; % =red: edge, =yellow: screw
    
	plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3), ...
        'Color',color,'LineWidth',2);
    quiver3(midseg(1),midseg(2),midseg(3),bplot(1),bplot(2),bplot(3), ...
        'Color','blue','LineWidth',1, ...
        'ShowArrowHead','on','MaxHeadSize',1e2);
end

%%% Plot nodes

for j = 1:size(rn,1)
    
    % Assign different colour to pinned nodes:
    if rn(j,4) == 7
        nodemarkertype = '*';
        nodemarkersize = 3;
        nodelinewidth = 2;
        nodecolor = [0.6350 0.0780 0.1840];
    elseif rn(j,4) == 67
        continue; % Skip virtual nodes
    else
        nodemarkertype = 'o';
        nodemarkersize = 1;
        nodelinewidth = 1;
        nodecolor = 'black';
    end
    
	plot3(rn(j,1),rn(j,2),rn(j,3), nodemarkertype, 'Color', nodecolor, ...
        'MarkerSize', nodemarkersize, 'LineWidth', nodelinewidth);
end

%%% Plot bounding box

face1 = [1 2 4 3 1];
face2 = [5 6 8 7 5];
vertices_scaled = vertices;
surf1=vertices_scaled(face1,:);
surf2=vertices_scaled(face2,:);

plot3(surf1(:,1),surf1(:,2),surf1(:,3),'k','LineWidth',2);
plot3(surf2(:,1),surf2(:,2),surf2(:,3),'k','LineWidth',2);

side = vertices_scaled([1 5],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([2 6],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([3 7],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);
side = vertices_scaled([4 8],:);
plot3(side(:,1),side(:,2),side(:,3),'k','LineWidth',2);

%%% Axes and labels

axis equal

xlabel('x','FontSize',12);
ylabel('y','FontSize',12);
zlabel('z','FontSize',12);
xlim([0 vertices(8,1)]);
ylim([0 vertices(8,2)]);
zlim([0 vertices(8,3)]);

title(figname,'FontSize',18); % Title

view(viewangle); % View angle
end