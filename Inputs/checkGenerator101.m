function [rn,links] = checkGenerator101(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz)

%dx,dy,dz: size of cantilever
%bufferfactor: factor normalised by dist_source indicating minimum distance
%from surfaces in which sources can be generated.

% if strcmp(CRYSTAL_STRUCTURE,'bcc')
%     disp('Crystal structure recognized');
% else
%     disp('Crystal structure not recognized. Aborting');
%     return;
% end
bufferfactor = 0.65;% > 1/2 
%NB Sources are idealised as squares...
Xmin = 0+bufferfactor*DIST_SOURCE;
Xmax = dx-bufferfactor*DIST_SOURCE;
Ymin = 0+bufferfactor*DIST_SOURCE;
Ymax = dy-bufferfactor*DIST_SOURCE;
Zmin = 0+bufferfactor*DIST_SOURCE;
Zmax = dz-bufferfactor*DIST_SOURCE;

%Generate midpoints of sources
 midX = Xmin + (Xmax - Xmin).*rand(NUM_SOURCES,1); %*0.75;
 midY = Ymin + (Ymax - Ymin).*rand(NUM_SOURCES,1); %*0.1; 
 midZ = Zmin + (Zmax - Zmin).*rand(NUM_SOURCES,1); %*0.45;
 midPTS = horzcat(midX,midY,midZ);
%midPTS(1,:) = [0.1*dx,0.5*dy,0.8*dz];

%Generate random {101} type habit planes for each source
% normal = slipplane(NUM_SOURCES);
normal = [0 1 0];
% generate random <101> b vector
% b_vec = pmone(normal,NUM_SOURCES);
b_vec = [-0.2 0.5 -0.35]; % use [101] for debugging

%We have thus defined slip-plane and b-vector of the loops.
%We try a loop placed within the slip-plane (shear loop)...

% seg_vec = cross(normal',b_vec');
edge = 0.5*[-0.35 .2 0.4]*DIST_SOURCE;%[1 0 1]*DIST_SOURCE;
screw = 0.25*[0.2 -0.5 0.3]*DIST_SOURCE;%[-1 0 1]*DIST_SOURCE;


rn = zeros(8*NUM_SOURCES,4);
links = zeros(8*NUM_SOURCES,8);

for p=1:NUM_SOURCES
    % pure edge segment
%     r1 = midPTS(p,:)  
%     r2 = midPTS(p,:) + 0.5*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%     r3 = midPTS(p,:) + 1*DIST_SOURCE*seg_vec(p,:)/norm(seg_vec(p,:));
%    
    % pure screw segment
%     r1 = midPTS(p,:) - 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:)); 
%     r2 = midPTS(p,:);
%     r3 = midPTS(p,:) + 0.5*DIST_SOURCE*b_vec(p,:)/norm(b_vec(p,:));
    
    r1 = midPTS(p,:)-edge;
    r2 = midPTS(p,:);
    r3 = midPTS(p,:) + edge;
    r4 = r3 + screw;
    r5 = r4 + screw;
    r6 = r5 - edge;
    r7 = r6 - edge;
    r8 = r7 - screw;

    rn((p-1)*8+1,:) = [r1 7]; 
    rn((p-1)*8+2,:) = [r2 0]; 
    rn((p-1)*8+3,:) = [r3 7]; 
    rn((p-1)*8+4,:) = [r4 0];
    rn((p-1)*8+5,:) = [r5 7];
    rn((p-1)*8+6,:) = [r6 0];
    rn((p-1)*8+7,:) = [r7 7];
    rn((p-1)*8+8,:) = [r8 0];
  
    for m = 1:7
        links((p-1)*8+m,1:2) = [(p-1)*8+m, (p-1)*8+(m+1)];
    end
    links((p-1)*8+8,1:2) = [(p-1)*8+8,(p-1)*8+1];
    
    links(((p-1)*8+1):((p-1)*8+8),3:5) = repmat(b_vec,8,1);
    links(((p-1)*8+1):((p-1)*8+8),6:8) = repmat(normal,8,1);
    
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
plot3(rn(:,1), rn(:,2), rn(:,3),'b-'); %nodes
plot3(rn(:,1), rn(:,2), rn(:,3),'b.'); %nodes
plot3(midPTS(:,1),midPTS(:,2),midPTS(:,3),'b*'); %midpoints
bvec_plot = repmat(b_vec,NUM_SOURCES,1);
quiver3(midPTS(:,1),midPTS(:,2),midPTS(:,3), bvec_plot(:,1),bvec_plot(:,2),bvec_plot(:,3))
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