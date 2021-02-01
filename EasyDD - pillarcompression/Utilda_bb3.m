function [Ux,Uy,Uz]=Utilda_bb3(rn,links,gnl,NU,xnodes,dx,dy,dz,mx,my,mz)

% C=[dx/mx,dy/my,dz/mz];
C=[dx/2,dy/2,dz/2];
nodenum=size(gnl,1);                                        %Number of FE nodes
segnum=size(links,1);                                       %Number of dislocation segments
Utilda=zeros(nodenum,3);                                    %Preallocate the displacement
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
    
    for j=1:nodenum
        nodepoint=xnodes((gnl(j)),1:3);                           %Coordinates of the FE node
          p=nodepoint';                                             %Column vector for use in displacement_et
          
          if out==0                                                   %For internal segments, calculate displacements using Barnett triangles with closure point
              Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',NU)+displacement_et_plas(p,A',B',C',b');
          else                                                        %For external segments, close internal loop with normal closure point and calculate displacement for external loop with temporary closure point
              Utilda(j,:)=Utilda(j,:)+displacement_et_plas(p,Aprime',Bprime',C',b')+displacement_et_plas(p,Aprime',A',Ctemp',b')+displacement_et_plas(p,A',B',Ctemp',b')+displacement_et_plas(p,B',Bprime',Ctemp',b')+displacement_et_plas(p,Bprime',Aprime',Ctemp',b');
          end
          
          if any(isnan(Utilda(j,:)))
              pause
          end
    end
    
end

    Ux=Utilda(:,1);                                                     %Organises outputs
    Uy=Utilda(:,2);
    Uz=Utilda(:,3);
end

function [u] = displacement_et_el(p,A,B,b,nu)

%Calculates only elastic components of displacement for segment AB using
%code from displacement_et

con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));


% R vectors (from P to nodes)
RA = A - p;
RB = B - p;

lamA = safenorm(RA);
lamB = safenorm(RB);

vecAB = B - A;

  
tAB = safenorm(vecAB);

    

% f = fab
f = fab(b, tAB, lamA, lamB, RA, RB);

% g = gab
g = gab(b, lamA, lamB);

% update displacement inc. solid angle
u = - con1.*f + con2.*g;

u = u';

end

function [u] = displacement_et_plas(p,A,B,C,b)
% Calculates the plastic (solid angle) components of displacement for
% triangle ABC using code from displacement_et


% R vectors (from P to nodes)
RA = A - p;
RB = B - p;
RC = C - p;

lamA = safenorm(RA);
lamB = safenorm(RB);
lamC = safenorm(RC);

vecAB = B - A;
vecBC = C - B;
  
tAB = safenorm(vecAB);
tBC = safenorm(vecBC);

% calculate slip plane normal
% vec_a = segs(1,:)/norm(segs(1,:));
%     vec_b = -segs(end,:)/norm(segs(end,:));

    n = cross(tAB,tBC);
    
omega  = solang(lamA, lamB, lamC,n);
% update displacement inc. solid angle
u =  - b*omega/(4*pi);

u = u';

end

function [ out ] = fab( b, tAB, lamA, lamB, RA, RB )
%calculates f vector

numerator = norm(RB)*(1 + dot(lamB,tAB));
denominator = norm(RA)*(1 + dot(lamA,tAB));

   % this can happen if field point is on the dislocation segment
    if abs(denominator) < eps
        logarithm  = 0;
    elseif abs(numerator) < eps
        logarithm = 0;
    else
        logarithm = real(log(numerator/denominator)); % for large step calculations, log will produce an imaginary 
    end    

out = cross(b,tAB)*logarithm;

end

function [ out ] = gab( b, lamA, lamB )
%calculates g vector

numerator = dot(b,cross(lamA,lamB)) * (lamA + lamB);
denominator = 1 + dot(lamA,lamB);

if abs(denominator) < eps

    out = [0; 0; 0];

else
    out = numerator/denominator;
end

end

function [ omega ] = solang( lamA, lamB, lamC, n)
%calculates solid angle due to triangular loop ABC
% 1) this may not be robust, check different formulae for the 
% spherical excess (wiki spherical geometry)
% 2) it could be a lot faster to calculate a,b,c in the main loop 
% for each pair of lambdas, then use a,b,c as input arguments 

a = acos(dot(lamB,lamC));
b = acos(dot(lamC,lamA));
c = acos(dot(lamA,lamB));

s = (a + b + c)/2;

svec = [s, s-a, s-b, s-c]/2;

temp = tan(svec(1))*tan(svec(2))*tan(svec(3))*tan(svec(4));

if temp < eps
    temp = 0;
end
omega = 4*real(atan(sqrt(temp)));

sgn = sign(dot(lamA,n));

omega = -sgn*omega;


%-------- check sign explicitly for c code ------
% if dot(lamA,n) > eps
%     sgn2 = 1;
% elseif dot(lamA,n) < -eps    
%     sgn2 = -1;
% else
%     sgn2 = 0;
% end
% if sgn2 ~= sgn
%     disp('sign error')
%    % pause
% end
%--------------------------

end

function unitR = safenorm(R)

normR = sqrt(R(1)*R(1) + R(2)*R(2)+R(3)*R(3));

if normR > eps
    unitR = R/normR;
else
    unitR = zeros(3,1);
end

end