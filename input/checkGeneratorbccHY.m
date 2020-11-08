function [rn,links,normal,b_vec] = checkGeneratorbccHY(NUM_SOURCES,DIST_SOURCE,CRYSTAL_STRUCTURE,dx,dy,dz,plim,LENGTH_SOURCE)
%HY20171128: modified to check plane norm and b_vec.
%HY20180114: modified by HY to consider all 12 bcc slip systems, for HHDDD
%simulation on bcc iron
%dx,dy,dz: size of cantilever
%bufferfactor: factor normalised by dist_source indicating minimum distance
%from surfaces in which sources can be generated.

eps = 1e-12;
amag = 2.856e-4; 

if strcmp(CRYSTAL_STRUCTURE,'bcc')
    disp('Crystal structure recognized');
else
    disp('Crystal structure not recognized. Aborting');
    return;
end

ALLslips = [1 0 1 -1 1 1;
            1 0 1 1 1 -1;
            1 0 -1 1 1 1;
            1 0 -1 1 -1 1;
            0 1 1 1 -1 1;
            0 1 1 1 1 -1;
            0 1 -1 1 1 1;
            0 1 -1 -1 1 1;
            1 1 0 -1 1 1;
            1 1 0 1 -1 1;
            1 -1 0 1 1 1;
            1 -1 0 1 1 -1];

bufferfactor = 1% > 1/2 

%NB Sources are idealised as squares...
Xmin = dx*0.1+bufferfactor*DIST_SOURCE;
Xmax = dx*0.6-bufferfactor*DIST_SOURCE;
Ymin = dy*0.5+bufferfactor*DIST_SOURCE;
Ymax = dy*0.5-bufferfactor*DIST_SOURCE;
Zmin = dz*0.1+bufferfactor*DIST_SOURCE;
Zmax = dz-bufferfactor*DIST_SOURCE;

%Generate midpoints of sources
 midX = Xmin + (Xmax - Xmin).*rand(NUM_SOURCES,1);
 midY = Ymin + (Ymax - Ymin).*rand(NUM_SOURCES,1);
 midZ = Zmin + (Zmax - Zmin).*rand(NUM_SOURCES,1);
 midPTS = horzcat(midX,midY,midZ);
%midPTS(1,:) = [0.1*dx,0.5*dy,0.8*dz];
%              Screw
%         --------------
%         |            |
%         |      b     |
%   Edge  |     -->    | Edge
%         |            |
%         |------------|
%              Screw
%HY20171128: This is the in-plane case: edge line direction is always perpendicular to b_vec; screw,
%parallel. so, [1 0 1] is the burgers vector, [0 1 0] is selected as
%perpendicular to burgers vector but remaining inside the slip plane. So
%how to determine this a=[0 1 0]? - by these two conditions: a*burgers=0
%and a*normal=0
rn = zeros(8*NUM_SOURCES,4);
links = zeros(8*NUM_SOURCES,8);

for p=1:NUM_SOURCES
    
    index = randi([1,3]);
    
    normal = ALLslips(index,1:3);
    b_vec = ALLslips(index,4:6);
    
    b_vec = b_vec/norm(b_vec);
    b_vec = sqrt(3)/2*b_vec/norm(b_vec);
% %HY20180415: Rotate the loop around global x axis by alpha, y axis by beta and z axis by
% %theta 
% 
% alpha = 45/180*pi;
% beta = (90-asin(sqrt(2)/sqrt(3))*180/3.14)/180*pi;
% theta = (atan(0.7071/1.2247)*180/3.14)/180*pi;
% 
% Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
% Rz = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% 
% %HY20180124: NOTE! The sequence of rotation is important!!!!
% Q = Rx'*Ry'*Rz';
% 
% normal = normal*Q;
% b_vec = b_vec*Q;
% 
% 
% %HY20180415: rotate to the maximum shear direction in the xz plane
% 
% 
% alpha = 0/180*pi;
% beta = -45/180*pi;
% theta = 0/180*pi;
% 
% Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
% Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
% Rz = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% 
% Q = Rx'*Ry'*Rz';
% 
% normal = normal*Q;
% b_vec = b_vec*Q;
    
    bvec_plot(p,:) = b_vec;
    normal_plot(p,:) = normal;
    
    %HY20180414: prismatic loops
    seg_vec = cross(normal,b_vec);
    screw = 1*normal/norm(normal)*LENGTH_SOURCE/2;
    edge =  1*(seg_vec/norm(seg_vec))*LENGTH_SOURCE/2;

    temp=screw;
    screw=edge;
    edge=temp;
    
    r2 = midPTS(p,:)-screw;
    r1 = r2-edge;
    r2 = r2;
    r3 = r2 + edge;
    r4 = r3 + screw;
    r5 = r4 + screw;
    r6 = r5 - edge;
    r7 = r6 - edge;
    r8 = r7 - screw;

    %HY20180124: "4-8"or"2-6"
    rn((p-1)*8+1,:) = [r1 7]; 
    rn((p-1)*8+2,:) = [r2 7]; 
    rn((p-1)*8+3,:) = [r3 7]; 
    rn((p-1)*8+4,:) = [r4 0];
    rn((p-1)*8+5,:) = [r5 7];
    rn((p-1)*8+6,:) = [r6 7];
    rn((p-1)*8+7,:) = [r7 7];
    rn((p-1)*8+8,:) = [r8 0];

%     rn((p-1)*8+1,:) = [r1 0]; 
%     rn((p-1)*8+2,:) = [r2 0]; 
%     rn((p-1)*8+3,:) = [r3 0]; 
%     rn((p-1)*8+4,:) = [r4 0];
%     rn((p-1)*8+5,:) = [r5 0];
%     rn((p-1)*8+6,:) = [r6 0];
%     rn((p-1)*8+7,:) = [r7 0];
%     rn((p-1)*8+8,:) = [r8 0];
  
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
viewangle=[0,0];
clf;
figure(1); 
plotnodes(rn*amag,links,plim*amag,vertices*amag); 
view(viewangle);
hold on;
% plot3(rn(:,1), rn(:,2), rn(:,3),'b-'); %nodes
% plot3(rn(:,1), rn(:,2), rn(:,3),'b.'); %nodes
plot3(midPTS(:,1)*amag,midPTS(:,2)*amag,midPTS(:,3)*amag,'b*'); %midpoints
% bvec_plot = repmat(b_vec,NUM_SOURCES,1);
% normal_plot = repmat(normal,NUM_SOURCES,1);
bscaler = 1e3;
nscaler = 1e3;
quiver3(midPTS(:,1)*amag,midPTS(:,2)*amag,midPTS(:,3)*amag, bscaler*bvec_plot(:,1)*amag,bscaler*bvec_plot(:,2)*amag,bscaler*bvec_plot(:,3)*amag, 'g','LineWidth',1.5)
quiver3(midPTS(:,1)*amag,midPTS(:,2)*amag,midPTS(:,3)*amag, nscaler*normal_plot(:,1)*amag,nscaler*normal_plot(:,2)*amag,nscaler*normal_plot(:,3)*amag,'r','LineWidth',1.5)
% quiver3(midPTS(p,1),midPTS(p,2),midPTS(p,3), 1e4*b_vec(1),1e4*b_vec(2),1e4*b_vec(3), 'g')
% quiver3(midPTS(p,1),midPTS(p,2),midPTS(p,3), 1e4*normal(1),1e4*normal(2),1e4*normal(3),'r')
% view([0 0])

%HY20180122: expansion direction of the loop, judging from PK force
%HY20180122: loop number ii
ii = 1;
tensor = [1 0 0;0 0 0;0 0 0];
burgers = bvec_plot(ii,:);
normalll = normal_plot(ii,:);
rnexpansion = rn((ii-1)*8+1:ii*8,1:3);

figure(2); 
% plotnodes(rn,links,plim,vertices); 
view(viewangle);
hold on;
plot3(rnexpansion(1:8,1), rnexpansion(1:8,2), rnexpansion(1:8,3),'b-'); %nodes
plot3(rnexpansion(1:8,1), rnexpansion(1:8,2), rnexpansion(1:8,3),'b.'); %nodes
plot3(midPTS(ii,1),midPTS(ii,2),midPTS(ii,3),'b*'); %midpoints
for i = 1:3
    seg(i,:) = rnexpansion(2*i+1,:)-rnexpansion(2*i-1,:);
    burgers;
    PK1 = tensor*burgers';
    PK2 = seg(i,:)';
    PK(i,:) = cross(PK1,PK2);
    %PK(i,:) = PK(i,:)/norm(PK(i,:));
%     PK(i,:)
    quiver3(rnexpansion(2*i,1),rnexpansion(2*i,2),rnexpansion(2*i,3), PK(i,1),PK(i,2),PK(i,3),'b')
end

bscaler2 =1E3;
nscaler2 =1E3;

burgers = burgers*bscaler2;
normalll = normalll*nscaler2;
quiver3(midPTS(ii,1),midPTS(ii,2),midPTS(ii,3), burgers(1),burgers(2),burgers(3), 'g')
quiver3(midPTS(ii,1),midPTS(ii,2),midPTS(ii,3), normalll(1),normalll(2),normalll(3), 'r')
plot3(rn(4,1), rn(4,2), rn(4,3),'r*'); %node no.4
% axis equal

end
