function [rn,links] = bccsurfsourcegen(NUM_SOURCES,DIST_SOURCE,dx,dy,dz)

slipsys=[1, 1, 0, -1,  1,  1;
         1, 1, 0,  1, -1,  1;
        -1, 1, 0,  1,  1,  1;
        -1, 1, 0,  1,  1, -1;
         1, 0, 1,  1,  1, -1;
         1, 0, 1, -1,  1,  1;
        -1, 0, 1,  1,  1,  1;
        -1, 0, 1,  1, -1,  1;
         0, 1, 1,  1, -1,  1;
         0, 1, 1,  1,  1, -1;
         0,-1, 1,  1,  1,  1;
         0,-1, 1, -1,  1,  1];
%          1, 1, 2,  1,  1, -1;
%          1, 1,-2,  1,  1,  1;
%          1,-1, 2, -1,  1,  1;
%         -1, 1, 2,  1, -1,  1;
%          1, 2, 1,  1, -1,  1;
%          1, 2,-1, -1,  1,  1;
%          1,-2, 1,  1,  1,  1;
%         -1, 2, 1,  1,  1, -1;
%          2, 1, 1, -1,  1,  1;
%          2, 1,-1,  1, -1,  1;
%          2,-1, 1,  1,  1, -1;
%         -2, 1, 1,  1,  1,  1;
%          1, 2, 3,  1,  1, -1;
%          1, 2,-3,  1,  1,  1;
%          1,-2, 3, -1,  1,  1;
%         -1, 2, 3,  1, -1,  1;
%          2, 1, 3,  1,  1, -1;
%          2, 1,-3,  1,  1,  1;
%          2,-1, 3, -1,  1,  1;
%         -2, 1, 3,  1, -1,  1;
%          1, 3, 2,  1, -1,  1;
%          1, 3,-2, -1,  1,  1;
%          1,-3, 2,  1,  1,  1;
%         -1, 3, 2,  1,  1, -1;
%          2, 3, 1,  1, -1,  1;
%          2, 3,-1, -1,  1,  1;
%          2,-3, 1,  1,  1,  1;
%         -2, 3, 1,  1,  1, -1;
%          3, 1, 2, -1,  1,  1;
%          3, 1,-2,  1, -1,  1;
%          3,-1, 2,  1,  1, -1;
%         -3, 1, 2,  1,  1,  1;
%          3, 2, 1, -1,  1,  1;
%          3, 2,-1,  1, -1,  1;
%          3,-2, 1,  1,  1, -1;
%         -3, 2, 1,  1,  1,  1;
%          1, 1, 0,  1, -1, -1;
%          1, 1, 0, -1,  1, -1;
%         -1, 1, 0, -1, -1, -1;
%         -1, 1, 0, -1, -1,  1;
%          1, 0, 1, -1, -1,  1;
%          1, 0, 1,  1, -1, -1;
%         -1, 0, 1, -1, -1, -1;
%         -1, 0, 1, -1,  1, -1;
%          0, 1, 1, -1,  1, -1;
%          0, 1, 1, -1, -1,  1;
%          0,-1, 1, -1, -1, -1;
%          0,-1, 1,  1, -1, -1;
%          1, 1, 2, -1, -1,  1;
%          1, 1,-2, -1, -1, -1;
%          1,-1, 2,  1, -1, -1;
%         -1, 1, 2, -1,  1, -1;
%          1, 2, 1, -1,  1, -1;
%          1, 2,-1,  1, -1, -1;
%          1,-2, 1, -1, -1, -1;
%         -1, 2, 1, -1, -1,  1;
%          2, 1, 1,  1, -1, -1;
%          2, 1,-1, -1,  1, -1;
%          2,-1, 1, -1, -1,  1;
%         -2, 1, 1, -1, -1, -1;
%          1, 2, 3, -1, -1,  1;
%          1, 2,-3, -1, -1, -1;
%          1,-2, 3,  1, -1, -1;
%         -1, 2, 3, -1,  1, -1;
%          2, 1, 3, -1, -1,  1;
%          2, 1,-3, -1, -1, -1;
%          2,-1, 3,  1, -1, -1;
%         -2, 1, 3, -1,  1, -1;
%          1, 3, 2, -1,  1, -1;
%          1, 3,-2,  1, -1, -1;
%          1,-3, 2, -1, -1, -1;
%         -1, 3, 2, -1, -1,  1;
%          2, 3, 1, -1,  1, -1;
%          2, 3,-1,  1, -1, -1;
%          2,-3, 1, -1, -1, -1;
%         -2, 3, 1, -1, -1,  1;
%          3, 1, 2,  1, -1, -1;
%          3, 1,-2, -1,  1, -1;
%          3,-1, 2, -1, -1,  1;
%         -3, 1, 2, -1, -1, -1;
%          3, 2, 1,  1, -1, -1;
%          3, 2,-1, -1,  1, -1;
%          3,-2, 1, -1, -1,  1;
%         -3, 2, 1, -1, -1, -1];

bufferfactor = 1.5;
%NB Sources are idealised as squares...
Xmin = 0+bufferfactor*DIST_SOURCE+(0.8/3)*dx;
Xmax = Xmin+(dx-Xmin)/5;%dx*0.75-bufferfactor*DIST_SOURCE;
Ymin = 0+bufferfactor*DIST_SOURCE;
Ymax = dy-bufferfactor*DIST_SOURCE;
Zmin = 0+bufferfactor*DIST_SOURCE;
Zmin2 = Zmin+dz/5;
Zmax = dz-bufferfactor*DIST_SOURCE;
Zmax2 = Zmax-dz/5;
edgevecs(:,1:3)=cross(slipsys(:,1:3),slipsys(:,4:6));

if mod(NUM_SOURCES,2)==1
    NUM_SOURCES=NUM_SOURCES+1;
end
NUM_SOURCES=NUM_SOURCES/2;

rn = zeros(size(slipsys,1)*8*NUM_SOURCES,4);
links = zeros(size(slipsys,1)*8*NUM_SOURCES,8);

for surf=1:2
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
    midX = Xmin + (Xmax - Xmin).*rand(NUM_SOURCES,1);
    midY = Ymin + (Ymax - Ymin).*rand(NUM_SOURCES,1);

    if surf==1
    midZ = Zmin + (Zmin2 - Zmin).*rand(NUM_SOURCES,1);
    else
    midZ = Zmax2 + (Zmax - Zmax2).*rand(NUM_SOURCES,1);
    end

    midPTS = horzcat(midX,midY,midZ);

    for p=1:NUM_SOURCES

        r1=midPTS(p,:)-mobvec-fixvec;
        r2=r1+mobvec;
        r3=r2+mobvec;
        r4=r3+fixvec;
        r5=r4+fixvec;
        r6=r5-mobvec;
        r7=r6-mobvec;
        r8=r7-fixvec;

        rn((i-1)*NUM_SOURCES*8+(p-1)*8+1,:)=[r1,7];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+2,:)=[r2,0];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+3,:)=[r3,7];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+4,:)=[r4,7];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+5,:)=[r5,7];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+6,:)=[r6,0];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+7,:)=[r7,7];
        rn((i-1)*NUM_SOURCES*8+(p-1)*8+8,:)=[r8,7];

        for m = 1:7
            links((i-1)*NUM_SOURCES*8+(p-1)*8+m,1:2) = [(i-1)*NUM_SOURCES*8+(p-1)*8+m, (i-1)*NUM_SOURCES*8+(p-1)*8+(m+1)];
        end
        links((i-1)*NUM_SOURCES*8+(p-1)*8+8,1:2) = [(i-1)*NUM_SOURCES*8+(p-1)*8+8,(i-1)*NUM_SOURCES*8+(p-1)*8+1];

        links(((i-1)*NUM_SOURCES*8+(p-1)*8+1):((i-1)*NUM_SOURCES*8+(p-1)*8+8),3:5) = repmat(b_vec,8,1);
        links(((i-1)*NUM_SOURCES*8+(p-1)*8+1):((i-1)*NUM_SOURCES*8+(p-1)*8+8),6:8) = repmat(normal,8,1);

    end


end

if surf==1
    rn1=rn;
    links1=links;
else
    rn2=rn;
    links2=links;
    links2(:,1)=links2(:,1)+size(slipsys,1)*8*NUM_SOURCES;
    links2(:,2)=links2(:,2)+size(slipsys,1)*8*NUM_SOURCES;
end

end

rn=[rn1;rn2];
links=[links1;links2];

end
