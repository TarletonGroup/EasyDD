function [rn,links] = bccsourcegen_bicrystal(NUM_SOURCES,DIST_SOURCE,dx,dy,dz)

slipsys=[
%      1, 1, 0, -1,  1,  1;
%      1, 1, 0,  1, -1,  1;
%       -1, 1, 0,  1,  1,  1;
%      -1, 1, 0,  1,  1, -1;
%    1, 0, 1,  1,  1, -1;
       1, 0, 1, -1,  1,  1;
%      -1, 0, 1,  1,  1,  1;
%       -1, 0, 1,  1, -1,  1;
      0, 1, 1,  1, -1,  1;
%         0, 1, 1,  1,  1, -1;
%       0,-1, 1,  1,  1,  1;
%     0,-1, 1, -1,  1,  1
    ]; %plane normal & burgers vector
%% plane normal = glide plane

bufferfactor = 2;
%bufferfactor = 4;
%NB Sources are idealised as squares...
Xmin = 0+bufferfactor*DIST_SOURCE;
%Xmax = dx-Xmin;%dx*0.75-bufferfactor*DIST_SOURCE;
Xmax = 1.025*dx+bufferfactor*DIST_SOURCE;
Ymin = 0+bufferfactor*DIST_SOURCE;
Ymax = dy-bufferfactor*DIST_SOURCE;
Zmin = 0+bufferfactor*DIST_SOURCE;
Zmax = 0.78*dz-bufferfactor*DIST_SOURCE;
edgevecs(:,1:3)=cross(slipsys(:,1:3),slipsys(:,4:6));


rn = zeros(size(slipsys,1)*8*(NUM_SOURCES-1),4);
links = zeros(size(slipsys,1)*8*(NUM_SOURCES-1),8);

for i=1:size(slipsys,1)
    normal=slipsys(i,1:3);
    normal=normal/norm(normal);
%     screw=slipsys(1,4:6);
%     screw=screw/norm(screw);
    edge=edgevecs(i,1:3);
    edge=edge/norm(edge);
    b_vec=slipsys(i,4:6);
    b_vec=b_vec/norm(b_vec);
    
    mobvec=DIST_SOURCE*edge;
    fixvec=DIST_SOURCE*normal;
    %Generate midpoints of sources
    midX = Xmin + (Xmax - Xmin);
    midY = Ymin + (Ymax - Ymin);
    midZ = Zmin + (Zmax - Zmin);
    midPTS = horzcat(midX,midY,midZ);
    
    for p=1:NUM_SOURCES
        
          midPTS(1,1) = midPTS(1,1)-2*bufferfactor*DIST_SOURCE;
        r1=midPTS-mobvec-fixvec;
        r2=r1+mobvec;
        r3=r2+mobvec;
        r4=r3+fixvec;
        r5=r4+fixvec;
        r6=r5-mobvec;
        r7=r6-mobvec;
        r8=r7-fixvec;
        
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+1,:)=[r1,7];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+2,:)=[r2,0];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+3,:)=[r3,7];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+4,:)=[r4,7];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+5,:)=[r5,7];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+6,:)=[r6,0];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+7,:)=[r7,7];
        rn((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+8,:)=[r8,7];
        
        for m = 1:7
            links((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+m,1:2) = [(i-1)*(NUM_SOURCES-1)*8+(p-1)*8+m, (i-1)*(NUM_SOURCES-1)*8+(p-1)*8+(m+1)];
        end
        links((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+8,1:2) = [(i-1)*(NUM_SOURCES-1)*8+(p-1)*8+8,(i-1)*(NUM_SOURCES-1)*8+(p-1)*8+1];
    
        links(((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+1):((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+8),3:5) = repmat(b_vec,8,1);
        links(((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+1):((i-1)*(NUM_SOURCES-1)*8+(p-1)*8+8),6:8) = repmat(normal,8,1);
        
%     midPTS(1,1) = midPTS(1,1)-2*bufferfactor*DIST_SOURCE;
        
    end


end


flag_screw=0;
for i=1:size(links,1)
    
    linevec=rn(links(i,2),1:3)-rn(links(i,1),1:3);
    linevec=linevec/norm(linevec);
    burvec=links(i,3:5);
    normplane=cross(burvec,linevec);
    if norm(normplane)>1e-6
        normplane=normplane/norm(normplane);
    else
        fprintf('pure screw seg in initial config, see bccsourcegen.m')
        flag_screw=1;
    end
    if ~flag_screw
        links(i,6:8)=normplane;
    end
end
   

rn2(1:8,:) = rn(1:8,:);
rn2(9:16,:) = rn(17:24,:);
links2 = zeros(16,8);
links2(1:8,:) = links(1:8,:);
links2(9:16,1:2) = links(9:16,1:2);
links2(9:16,3:8) = links(17:24,3:8);
clear links
clear rn
rn = rn2;
links = links2;

fid =   fopen('.\input\initial_structure','w');
ll1=size(rn,1);
fwrite(fid,ll1,'integer*4');
fwrite(fid,rn,'real*8');
ll2=size(links,1);
fwrite(fid,ll2,'integer*4');
fwrite(fid,links,'real*8');
fclose(fid);
end