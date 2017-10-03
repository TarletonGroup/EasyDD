function [rn,links] = FRSGenerator(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz)

%dx,dy,dz: size of cantilever
%bufferfactor: factor normalised by dist_source indicating minimum distance
%from surfaces in which sources can be generated.

if strcmp(CRYSTAL_STRUCTURE,'bcc')
    disp('Crystal structure recognized');
else
    disp('Crystal structure not recognized. Aborting');
    return;
end
bufferfactor = 0.6% > 1/2 
%NB Sources are idealised as squares...
Xmin = 0+bufferfactor*DIST_SOURCE;
Xmax = dx-bufferfactor*DIST_SOURCE;
Ymin = 0+bufferfactor*DIST_SOURCE;
Ymax = dy-bufferfactor*DIST_SOURCE;
Zmin = 0+bufferfactor*DIST_SOURCE;
Zmax = dz-bufferfactor*DIST_SOURCE;

%Generate midpoints of sources
midX = Xmin + (Xmax - Xmin).*rand(NUM_SOURCES,1);
midY = Ymin + (Ymax - Ymin).*rand(NUM_SOURCES,1);
midZ = Zmin + (Zmax - Zmin).*rand(NUM_SOURCES,1);
midPTS = horzcat(midX,midY,midZ);
% midPTS(1,:) = [0.2*dx,0.5*dy,0.9*dz];
% for p = 2:NUM_SOURCES
%     midPTS(p,:) = midPTS(1,:) + (p-1)*[0.02*dx,0 0 ];
% end
%Generate random {110} type habit planes for each source
% normal = slipplane(NUM_SOURCES);
normal = repmat([-1 0 1],NUM_SOURCES,1);
% generate random <111> b vector
% b_vec = pmone(normal,NUM_SOURCES);
b_vec = repmat([1 0 1],NUM_SOURCES,1); % use [111] for debugging

%We have thus defined slip-plane and b-vector of the loops.
%We try a loop placed within the slip-plane (shear loop)...

% seg_vec = cross(normal,b_vec,2);
seg_vec = [0 0 0];

rn = zeros(3*NUM_SOURCES,4);
links = zeros(2*NUM_SOURCES,8);

for p=1:NUM_SOURCES %for each source...

    % pure edge segment
%     r1 = midPTS(p,:)  
%     r2 = midPTS(p,:) + 0.5*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%     r3 = midPTS(p,:) + 1*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%    
    % pure screw segment
%     r1 = midPTS(p,:) - 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:)); 
%     r2 = midPTS(p,:);
%     r3 = midPTS(p,:) + 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:));
    
    r1 = midPTS(p,:)-0.5*DIST_SOURCE*seg_vec;
    r2 = midPTS(p,:);
    r3 = midPTS(p,:) + 0.5*DIST_SOURCE*seg_vec;


    rn(3*p-2,:) = [r1 7]; 
    rn(3*p-1,:) = [r2 0]; 
    rn(3*p,:) = [r3 7]; 

  
    % this is wrong, fix it this afternoon
    links(2*p-1,1:2) = [3*p-2, 3*p-1];
    links(2*p,1:2) = [3*p-1, 3*p];
   
    
    links(2*p-1:2*p,3:5) = repmat(b_vec(p,:),2,1);
    links(2*p-1:2*p,6:8) = repmat(normal(p,:),2,1);
end

vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];
% plotnodes(rn,links,0,vertices);
hold on;
plot3(rn(:,1), rn(:,2), rn(:,3),'ko'); %nodes
plot3(midPTS(:,1),midPTS(:,2),midPTS(:,3),'b*'); %midpoints
quiver3(midPTS(:,1),midPTS(:,2),midPTS(:,3), b_vec(:,1),b_vec(:,2),b_vec(:,3))
end

function normal = slipplane(NUM_SOURCES)
    
    normal = zeros(NUM_SOURCES,3);
    for i=1:NUM_SOURCES
        
        random = rand;
        
        if random<0.3333
            normal(i,:) = [1 1 0];
            if rand < 0.5
                normal(i,:) = [1 -1 0];
            end
        elseif random > 0.3333 && random < 0.6666
            normal(i,:) = [1 0 1];
            if rand < 0.5
                normal(i,:) = [1 0 -1];
            end
        else
            normal(i,:) = [0 1 1];
            if rand < 0.5
                normal(i,:) = [0 1 -1];
            end
        end
    end
        
end

function b_vec = pmone(normal,NUM_SOURCES)

    b_vec = zeros(NUM_SOURCES,3);
    b_combinations = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; 1 -1 -1; -1 1 -1; -1 -1 -1];
    
    for i = 1:NUM_SOURCES
        for j=1:8;
            if dot( normal(i,:) , b_combinations(j,:) ) == 0
                b_vec(i,:) = b_combinations(j,:);
                continue;
            end
        end
    end  
end