function [rn,links] = DislocationSourceGenerator(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz,bufferfactor)

%dx,dy,dz: size of cantilever
%bufferfactor: factor normalised by dist_source indicating minimum distance
%from surfaces in which sources can be generated.

if strcmp(CRYSTAL_STRUCTURE,'bcc')
    disp('Crystal structure recognized');
else
    disp('Crystal structure not recognized. Aborting');
    return;
end

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

%Generate random {110} type habit planes for each source
%slip_array = slipplane(NUM_SOURCES);
slip_array = repmat([-1 0 1],NUM_SOURCES,1);

%b_array = pmone(slip_array,NUM_SOURCES);
b_array = repmat([1 0 1],NUM_SOURCES,1);

%We have thus defined slip-plane and b-vector of the loops.
%We try a loop placed within the slip-plane (shear loop)...

seg_vec = cross(slip_array,b_array,2);
%seg_vec=[0 0 0];

rn = zeros(8*NUM_SOURCES,4);
links = zeros(8*NUM_SOURCES,8);

for p=1:NUM_SOURCES %for each source...

%     %inplane
%     r1 = midPTS(p,:) + 0.5*DIST_SOURCE*(seg_vec(p,:) + b_array(p,:));
%     r2 = midPTS(p,:) + 0.5*DIST_SOURCE*(-seg_vec(p,:) + b_array(p,:));
%     r3 = midPTS(p,:) + 0.5*DIST_SOURCE*(-seg_vec(p,:) - b_array(p,:));
%     r4 = midPTS(p,:) + 0.5*DIST_SOURCE*(seg_vec(p,:) - b_array(p,:));
    
    r1 = midPTS(p,:) + 0.5*DIST_SOURCE*(seg_vec(p,:) + slip_array(p,:))/norm(seg_vec(p,:) + slip_array(p,:));
    r2 = midPTS(p,:) + 0.5*DIST_SOURCE*(-seg_vec(p,:) + slip_array(p,:))/norm(-seg_vec(p,:) + slip_array(p,:));
    r3 = midPTS(p,:) + 0.5*DIST_SOURCE*(-seg_vec(p,:) - slip_array(p,:))/norm(-seg_vec(p,:) - slip_array(p,:));
    r4 = midPTS(p,:) + 0.5*DIST_SOURCE*(seg_vec(p,:) - slip_array(p,:))/norm(seg_vec(p,:) - slip_array(p,:));

    rn((p-1)*8+1,:) = [r1 7]; 
    rn((p-1)*8+3,:) = [r2 7]; 
    rn((p-1)*8+5,:) = [r3 7]; 
    rn((p-1)*8+7,:) = [r4 7]; 
    
    rn((p-1)*8+2,:) = [0.5*(r1+r2) 0]; 
    rn((p-1)*8+4,:) = [0.5*(r2+r3) 0]; 
    rn((p-1)*8+6,:) = [0.5*(r3+r4) 0]; 
    rn((p-1)*8+8,:) = [0.5*(r4+r1) 0]; 
    
    links((p-1)*8+1,1:2) = [(p-1)*8+1, (p-1)*8+2];
    links((p-1)*8+2,1:2) = [(p-1)*8+2, (p-1)*8+3];
    links((p-1)*8+3,1:2) = [(p-1)*8+3, (p-1)*8+4];
    links((p-1)*8+4,1:2) = [(p-1)*8+4, (p-1)*8+5];
    links((p-1)*8+5,1:2) = [(p-1)*8+5, (p-1)*8+6];
    links((p-1)*8+6,1:2) = [(p-1)*8+6, (p-1)*8+7];
    links((p-1)*8+7,1:2) = [(p-1)*8+7, (p-1)*8+8];
    links((p-1)*8+8,1:2) = [(p-1)*8+8, (p-1)*8+1];
    
    links(((p-1)*8+1:(p-1)*8+8),3:5) = repmat(b_array(p,:),8,1);
    links(((p-1)*8+1:(p-1)*8+8),6:8) = repmat(slip_array(p,:),8,1);
end

vertices = [0,0,0;...
            dx,0,0;...
            0,dy,0;...
            dx,dy,0;...
            0,0,dz;...
            dx,0,dz;...
            0,dy,dz;...
            dx,dy,dz];
plotnodes(rn,links,0,vertices);
hold on;
plot3(rn(:,1), rn(:,2), rn(:,3),'ko'); %nodes
plot3(midPTS(:,1),midPTS(:,2),midPTS(:,3),'b*'); %midpoints
quiver3(midPTS(:,1),midPTS(:,2),midPTS(:,3), b_array(:,1),b_array(:,2),b_array(:,3),0.1*dy)
end

function slip_array = slipplane(NUM_SOURCES)
    
    slip_array = zeros(NUM_SOURCES,3);
    for i=1:NUM_SOURCES
        
        random = rand;
        
        if random<0.3333
            slip_array(i,:) = [1 1 0];
            if rand < 0.5
                slip_array(i,:) = [1 -1 0];
            end
        elseif random > 0.3333 && random < 0.6666
            slip_array(i,:) = [1 0 1];
            if rand < 0.5
                slip_array(i,:) = [1 0 -1];
            end
        else
            slip_array(i,:) = [0 1 1];
            if rand < 0.5
                slip_array(i,:) = [0 1 -1];
            end
        end
    end
        
end

function b_array = pmone(slip_array,NUM_SOURCES)

    b_array = zeros(NUM_SOURCES,3);
    b_combinations = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; 1 -1 -1; -1 1 -1; -1 -1 -1];
    
    for i = 1:NUM_SOURCES
        for j=1:8;
            if dot( slip_array(i,:) , b_combinations(j,:) ) == 0
                b_array(i,:) = b_combinations(j,:);
                continue;
            end
        end
    end  
end