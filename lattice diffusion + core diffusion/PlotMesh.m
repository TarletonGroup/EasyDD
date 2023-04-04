function PlotMesh(coordinates,nodes,show)
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [X Y Z] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%           show - to dispaly nodal and element numbers 
%                  0 (default) - do not display 
%                  1           - display
%
% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   
%
% version 1: 28 August 2011
% Version 2: 16 September 2016
%--------------------------------------------------------------------------
if nargin == 2
    show = 0 ;
end
dimension = size(coordinates,2) ;  % Dimension of the mesh
nel = length(nodes) ;                  % number of elements
nnode = length(coordinates) ;          % total number of nodes in system
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
Z = zeros(nnel,nel) ;
if dimension == 3   % For 3D plots
    if nnel==4 % surface in 3D
        for iel=1:nel   
            nd = nodes(iel,:) ;
            X(:,iel)=coordinates(nd,1);    % extract x value of the node
            Y(:,iel)=coordinates(nd,2);    % extract y value of the node
            Z(:,iel)=coordinates(nd,3) ;   % extract z value of the node 
        end    
        % Plotting the FEM mesh
        figure
        fill3(X,Y,Z,'w')
        rotate3d ;
        title('Finite Element Mesh') ;
        axis off ;
        % display Node numbers and Element numbers
        if show ~= 0
            k = 1:nnode ;
            nd = k' ;
            for i = 1:nel
                text(X(:,i),Y(:,i),Z(:,i),int2str(nd(i)),....
                    'fontsize',8,'color','k');
                text(sum(X(:,i))/4,sum(Y(:,i))/4,sum(Z(:,i))/4,int2str(i),.....
                    'fontsize',10,'color','r') ;
            end
        end
        
    elseif nnel==8  % solid in 3D
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,nel) ;
        for e=1:nel
            nd=nodes(e,:);
            X(:,e) = coordinates(nd,1) ;
            Y(:,e) = coordinates(nd,2) ;
            Z(:,e) = coordinates(nd,3) ;
            XYZ{e} = [X(:,e)  Y(:,e) Z(:,e)] ;
        end
        % Plot FEM mesh 
        figure
        set(gcf,'color','w')
        axis off 
        cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,.......
            repmat({'Faces'},1,nel),repmat({fm},1,nel),......
            repmat({'FaceColor'},1,nel),repmat({'w'},1,nel));
        view(3)
        set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;
        % display Node numbers and Element numbers
        if show ~= 0
            k = 1:nnode ;
            nd = k' ;
            for i = 1:nel
                text(X(:,i),Y(:,i),Z(:,i),int2str(nd(i)),....
                    'fontsize',8,'color','k');
                text(sum(X(:,i))/8,sum(Y(:,i))/8,sum(Z(:,i))/8,int2str(i),.....
                    'fontsize',10,'color','r') ;
            end
        end
    end
    
elseif dimension == 2           % For 2D plots
    for iel=1:nel   
        nd = nodes(iel,:) ;
        X(:,iel)=coordinates(nd,1);    % extract x value of the node
        Y(:,iel)=coordinates(nd,2);    % extract y value of the node
    end
    
    % Plotting the FEM mesh, diaplay Node numbers and Element numbers
    figure
    fill(X,Y,'w')
    title('Finite Element Mesh') ;
    axis off ;
    if show ~= 0
        k = 1:nnode ;
        nd = k' ;
        for i = 1:nel
            text(X(:,i),Y(:,i),int2str(nd(i)),'fontsize',8,'color','k');
            text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
        end
    end
end