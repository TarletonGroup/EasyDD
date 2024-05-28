%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Designed to calculate the displacement vectors on a cuboid convex hull due to the
% presence of a dislocation structure where any virtual nodes have been projected
% away from the surface along the line of the Burgers vector except when
% the Burgers vector is perpendicular to the local surface normal
% 
% Bruce Bromage
% Michromechanical Testing Group
% Department of Materials, University of Oxford
% bruce.bromage@materials.ox.ac.uk
% October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ux, Uy, Uz] = Utilda_bb2(rn,links,gnl,nu,xnodes,dx,dy,dz,mx,my,mz)

%This function finds the displacement caused by all dislocation segements
%in links on all of the FE nodes in gnl.


nodenum=size(gnl,1);                                        %Number of FE nodes
segnum=size(links,1);                                       %Number of dislocation segments
Utilda=zeros(nodenum,3);                                    %Preallocate the displacement
x=dx/(mx);                                                  %x coordinate inside simulated cuboid
y=dy/(my);                                                  %y coordinate inside simulated cuboid
z=dz/(mz);                                                  %z coordinate inside simulated cuboid
C=[x,y,z];                                                  %Closure point inside simulated cuboid
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

for i=1:segnum
    
    if rn(links(i,1),4) == 6 && rn(links(i,2),4) == 67      %Ignore all segments connecting surface nodes to virtual nodes
        continue
    elseif rn(links(i,2),4) == 6 && rn(links(i,1),4) == 67
        continue
    end
    
    A=rn(links(i,1),1:3);                                   %Coordinates of first node connected to current segment
    B=rn(links(i,2),1:3);                                   %Coordinates of second node connected to current segment
    b=links(i,3:5);                                         %Burgers vector of AB
    bA=[links(links(:,1)==links(i,1),3:5);links(links(:,2)==links(i,1),3:5)];  %Average Burgers vector at node A
    bA=sum(bA,1);
    bA=bA/norm(bA);
    bB=[links(links(:,1)==links(i,2),3:5);links(links(:,2)==links(i,2),3:5)];  %Average Burgers vector at node B
    bB=sum(bB,1);
    bB=bB/norm(bB);
    proj_vecA=1e7*bA;                                                          %Backup projection vector for node A, should be matched with remesh_surf
    lineA=[A,bA];                                                              %Line of Burgers vector from node A for use inintersectLineMesh3d
    proj_vecB=1e7*bB;                                                          %Backup projection vector for node B, should be matched with remesh_surf
    lineB=[B,bB];                                                              %Line of Burgers vector from node B for use inintersectLineMesh3d
    AB=B-A;                                                                    %Vector for segment AB
    planenormal=cross(AB,b);                                                   %Slip plane normal for segment AB
    screwcheck=norm(planenormal);                                              %Check if AB is screw type
    
    if A == B                                                                  %Sanity check in case of faulty remeshing
        continue
    end
    
     CA=A-C;                                                                   %vector from node A to closure point
     CB=B-C;                                                                   %vector from node B to closure point
     
      if rn(links(i,1),4) == 67 ||  rn(links(i,2),4) == 67                     %Finds virtual segments
          out=1;                                                               %Flags virtual segments as external
          surfptsA = intersectLineMesh3d(lineA, vertices, faces);              %Finds points on the surface where Burgers vector from node A intersects with the surface
          surfptsB = intersectLineMesh3d(lineB, vertices, faces);              %Finds points on the surface where Burgers vector from node B intersects with the surface
                    
          if isempty(surfptsA)                                      %Projects node A along the Burgers vector in the case where a face has been fully exited
              Aprime=A-proj_vecA;
              CAprime=Aprime-C;
              normCA=norm(CA);
              normCAprime=norm(CAprime);
              
              if normCAprime >= normCA                              %Ensures projected node is closer to the volume than original virtual node
                  Aprime=A+proj_vecA;
              end

          else
              Aprime=surfptsA(1,:);                                 %Designates equivalent surface node for node A
              Aprime2=surfptsA(2,:);
              AAprime=Aprime-A;
              AAprime2=Aprime2-A;
              normAAprime=norm(AAprime);
              normAAprime2=norm(AAprime2);
          
          
              if normAAprime > normAAprime2                         %Ensures surface node on the correct face has been chosen
                  Aprime=Aprime2;
              end 
              
          end
          
          if isempty(surfptsB)
              Bprime=B-proj_vecB;                                   %Projects node B along the Burgers vector in the case where a face has been fully exited
              CBprime=Bprime-C;
              normCB=norm(CB);
              normCBprime=norm(CBprime);
              
              if normCBprime >= normCB                              %Ensures projected node is closer to the volume than original virtual node
                  Bprime=B+proj_vecB;
              end
              
          else
              Bprime=surfptsB(1,:);                                 %Designates equivalent surface node for node B
              Bprime2=surfptsB(2,:);
              BBprime=Bprime-B;
              BBprime2=Bprime2-B;
              normBBprime=norm(BBprime);
              normBBprime2=norm(BBprime2);
          
              if normBBprime > normBBprime2                         %Ensures surface node on the correct face has been chosen
                  Bprime=Bprime2; 
              end
              
          end
          
          planecheck=dot(planenormal,CA);                           %Checks if closure point is in the plane of the loop
          
          if screwcheck < eps                                       %Treat unprojected screw type virtual segments as internal
              out=0;
          elseif planecheck == 0                                    %Do not generate convex hull if closure point is in the plane of the loop
              checker=zeros(size(xnodes,1),1);
          else
              pyramid=[A;B;C;Aprime;Bprime];                        %Designates original virtual nodes, closure point and projected nodes as vertices
              tess=convhulln(pyramid,{'QJ','Pp'});                  %Finds the convex hull of the designated vertices
              tol=1e1;
              checker=inhull(xnodes(:,1:3),pyramid,tess,tol);       %Flags FE nodes inside the convex hull
          end
      else
          out=0;                                                    %Flags non-virtual segments
      end
            
      for j=1:nodenum
          
          nodepoint=xnodes((gnl(j)),1:3);                           %Coordinates of the FE node
          p=nodepoint';                                             %Column vector for use in displacement_et
          
          if out == 0                                               %Calculates displacement vector caused by internal segments using Barnett triangles
              Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',nu)+displacement_et_plas(p,A',B',C',b');
              
          elseif out == 1                                           %Identifies external virtual segments
              
              if checker(gnl(j)) == 0                               %Calculates displacement vector caused  by virtual segments and the segments connecting them to the surface using Barnett triangles
                  Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',nu)+displacement_et_el(p,B',Bprime',b',nu)+displacement_et_el(p,Aprime',A',b',nu)+displacement_et_plas(p,A',B',C',b')+displacement_et_plas(p,B',Bprime',C',b')+displacement_et_plas(p,Aprime',A',C',b');
                  
              elseif checker(gnl(j)) == 1                           %Identifies FE nodes inside the pyramidal convex hull and calcultes the correct displacement accordingly
                  Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',nu)+displacement_et_el(p,B',Bprime',b',nu)+displacement_et_el(p,Aprime',A',b',nu)+solang_correction(p,A',B',Bprime',Aprime',C',b');
              end
              
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

function [u] = solang_correction(p,A,B,Bprime,Aprime,C,b)

%Calculates the solid angle contribution dor the special case where the
%field point is inside the pyramid created when the loop connects to the
%closure point. Uses code from displacement_et

segs=[Aprime',A',; A',B'; B',Bprime'; Bprime',Aprime'];
num=size(segs,1);
omega_big=0;

for i=1:num
  
RA = segs(i,1:3)' - p;
RB = segs(i,4:6)' - p;
RC = C - p;

lamA = safenorm(RA);
lamB = safenorm(RB);
lamC = safenorm(RC);

vecAB = B - A;
vecBC = C - B;

tAB = safenorm(vecAB);
tBC = safenorm(vecBC);

n = cross(tAB,tBC);

addon = solang(lamA, lamB, lamC,n);

omega_big  = omega_big + addon;
    
end

RA = segs(num,4:6)' - p;
RB = segs(num,1:3)' - p;
RC = C - p;

lamA = safenorm(RA);
lamB = safenorm(RB);
lamC = safenorm(RC);

vecAB = B - A;
vecBC = C - B;

tAB = safenorm(vecAB);
tBC = safenorm(vecBC);

n = cross(tAB,tBC);

omega_small  = solang(lamA, lamB, lamC,n);

omega_corr = -sign(omega_big)*2*pi + omega_big + omega_small;

u =  - b*omega_corr/(4*pi);

u = u';

end

function unitR = safenorm(R)

normR = sqrt(R(1)*R(1) + R(2)*R(2)+R(3)*R(3));

if normR > eps
    unitR = R/normR;
else
    unitR = zeros(3,1);
end

end