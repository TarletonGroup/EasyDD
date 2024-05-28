function [Ux,Uy,Uz]=Utilda_bb3_vec(rn,links,gnl,NU,xnodes,dx,dy,dz,mx,my,mz)

C=[dx/mx,dy/my,dz/mz];
nodenum=size(gnl,1);                                        %Number of FE nodes
segnum=size(links,1);                                       %Number of dislocation segments
Utilda=zeros(nodenum,3);                                    %Preallocate the displacement
Utilda=Utilda';

nodepoints=zeros(nodenum,3);
for j=1:nodenum
        nodepoints(j,:)=xnodes((gnl(j)),1:3);
end

vertices=[0,0,0;                                            %Vertices of cuboid
          dx,0,0;
          dx,dy,0;
          dx,dy,dz;
          dx,0,dz;
          0,0,dz;
          0,dy,dz;
          0,dy,0];
faces=[1,2,3,8;                                             %Faces of cuboid as defined by vertices
       1,2,5,6;
       1,6,7,8;
       2,3,4,5;
       3,4,7,8;
       4,5,6,7];
normals=[1,0,0;                                             %Normalised surface normal vectors for all faces
         0,1,0;
         0,0,1];
surfnodecons=zeros(size(rn,1),1);                           %Preallocate for connections to surface nodes

for i=1:segnum                                              %This loop finds virtual nodes connected to surface nodes and records the connections in surfnodecons
    
    if rn(links(i,1),end) == 6 && rn(links(i,2),end) == 67
        surfnodecons(links(i,2))=links(i,1);
    elseif rn(links(i,2),end) == 6 && rn(links(i,1),end) == 67
        surfnodecons(links(i,1))=links(i,2);
    end
    
end

for i=1:segnum
    
    if rn(links(i,1),4) == 6 && rn(links(i,2),4) == 67      %Ignore all segments connecting surface nodes to virtual nodes
        continue
    elseif rn(links(i,2),4) == 6 && rn(links(i,1),4) == 67
        continue
    end
    
   A=rn(links(i,1),1:3);                                   %Coordinates of first node connected to current segment
    B=rn(links(i,2),1:3);                                   %Coordinates of second node connected to current segment
    b=links(i,3:5);                                         %Burgers vector of AB
    
    if rn(links(i,1),4) == 67 ||  rn(links(i,2),4) == 67                    %Finds virtual segments
        out=1;                                                              %Flags virtual segments as external
        surfptsA=[];
        surfptsB=[];
        
        for v=1:size(normals,1)                                             %Finds the points on the surface which intersect virtual nodes along the surface normals
            lineA=[A,normals(v,:)];
            lineB=[B,normals(v,:)];
            surfptsA=[surfptsA;intersectLineMesh3d(lineA,vertices,faces)];
            surfptsB=[surfptsB;intersectLineMesh3d(lineB,vertices,faces)];
        end
        
        if isempty(surfptsA)                                                %If no surface points were found for A, perturb A and try again (check for corners)
            perturb=800.*normals;                                           %The factor for perturb should be >= rmax
            
            for p=1:size(perturb,1)
                 Atemp1=A+perturb(p,:);
                 Atemp2=A-perturb(p,:);
                
                for v=1:size(normals,1)
                     lineA1=[Atemp1,normals(v,:)];
                     lineA2=[Atemp2,normals(v,:)];
                     surfptsA=[surfptsA;intersectLineMesh3d(lineA1,vertices,faces)];
                     surfptsA=[surfptsA;intersectLineMesh3d(lineA2,vertices,faces)];
                end
                
            end
            
        end
        
        if isempty(surfptsB)                                                %If no surface points were found for B, perturb B and try again
            perturb=800.*normals;                                           %The factor for perturb should be >= rmax
            
            for p=1:size(perturb,1)
                 Btemp1=B+perturb(p,:);
                 Btemp2=B-perturb(p,:);
                
                for v=1:size(normals,1)
                     lineB1=[Btemp1,normals(v,:)];
                     lineB2=[Btemp2,normals(v,:)];
                     surfptsB=[surfptsB;intersectLineMesh3d(lineB1,vertices,faces)];
                     surfptsB=[surfptsB;intersectLineMesh3d(lineB2,vertices,faces)]; %#ok<*AGROW>
                end
                
            end
            
        end
        
        [~,Aindex]=min((surfptsA(:,1)-A(1)).^2+(surfptsA(:,2)-A(2)).^2+(surfptsA(:,3)-A(3)).^2);    %Find surface point closest to A
        [~,Bindex]=min((surfptsB(:,1)-B(1)).^2+(surfptsB(:,2)-B(2)).^2+(surfptsB(:,3)-B(3)).^2);    %Find surface point closest to B
        Aprime=surfptsA(Aindex,:);                                                                  %Point where A was projected from
        Bprime=surfptsB(Bindex,:);                                                                  %Point where B was projected from
        
        if surfnodecons(links(i,1))~=0                                                              %If A is connected to a surface node, change surface point to coordinates of that node
            Aprime=rn(surfnodecons(links(i,1)),1:3);
        end
        
        if surfnodecons(links(i,2))~=0                                                              %If B is connected to a surface node, change surface point to coordinates of that node
            Bprime=rn(surfnodecons(links(i,2)),1:3);
        end
            
        Ctemp=A+0.5.*(Bprime-A);                                                                    %Temporary closure point for external loop
        
    else
        out=0;                                                    % Flags internal segments as internal
    end
    
    
          
          if out==0                                                   %For internal segments, calculate displacements using Barnett triangles with closure point
              Utilda=Utilda+disp_dislo_tri_ABC_gen(A',B',C',nodepoints',b',NU);
          else                                                        %For external segments, close internal loop with normal closure point and calculate displacement for external loop with temporary closure point
              Utilda=Utilda+disp_dislo_tri_ABC_plas(Aprime',Bprime',C',nodepoints',b')+disp_dislo_tri_ABC_plas(Aprime',A',Ctemp',nodepoints',b')+disp_dislo_tri_ABC_plas(A',B',Ctemp',nodepoints',b')+disp_dislo_tri_ABC_plas(B',Bprime',Ctemp',nodepoints',b')+disp_dislo_tri_ABC_plas(Bprime',Aprime',Ctemp',nodepoints',b');
          end
          
          if any(isnan(Utilda))
              fprintf('Error in line 130 of Utilda_bb3_vec: Some output values are NaN')
              pause
          end
    
end

    Utilda=Utilda';
    Ux=Utilda(:,1);                                                     %Organises outputs
    Uy=Utilda(:,2);
    Uz=Utilda(:,3);
end
