%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Designed to calculate the displacement vectors on a cuboid convex hull due to the
% presence of a dislocation structure where any virtual nodes have been projected
% away from the surface along the line of the Burgers vector 
% 
% Bruce Bromage
% Michromechanical Testing Group
% Department of Materials, University of Oxford
% bruce.bromage@materials.ox.ac.uk
% October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    bA=[links(links(:,1)==links(i,1),3:5);links(links(:,2)==links(i,1),3:5)];  %Average burgers vector at node A
    bA=sum(bA,1);
    bA=bA/norm(bA);
    bB=[links(links(:,1)==links(i,2),3:5);links(links(:,2)==links(i,2),3:5)];  %Average burgers vector at node B
    bB=sum(bB,1);
    bB=bB/norm(bB);
%     proj_vecA=1e5*bA;
    lineA=[A,bA];
%     proj_vecB=1e5*bB;
    lineB=[B,bB];
    AB=B-A;
    screwcheck=dot(AB,b);
    planenormal=cross(AB,b);
    
    if A == B
        continue
    end
      
      if rn(links(i,1),4) == 67 ||  rn(links(i,2),4) == 67    %Finds virtual segments
          out=1;
          
          surfptsA = intersectLineMesh3d(lineA, vertices, faces);
%           Aprime=A-proj_vecA;                                  %Projects first node on to the surface along the Burgers vector
          Aprime=surfptsA(1,:);
          Aprime2=surfptsA(2,:);
          CA=A-C;
%           CAprime=Aprime-C;
          AAprime=Aprime-A;
          AAprime2=Aprime2-A;
%           normCA=norm(CA);
%           normCAprime=norm(CAprime);
          normAAprime=norm(AAprime);
          normAAprime2=norm(AAprime2);
          
%           if normCAprime >= normCA                              %Ensures projected node is closer to the volume than original virtual node
%              Aprime=A+proj_vecA;
%           end
          
          if normAAprime > normAAprime2
              Aprime=Aprime2;
          end
          
          surfptsB = intersectLineMesh3d(lineB, vertices, faces);
%           Bprime=B-proj_vecB;                                  %Projects second node on to the surface along the Burgers vector
          Bprime=surfptsB(1,:);
          Bprime2=surfptsB(2,:);
%           CB=B-C;
%           CBprime=Bprime-C;
          BBprime=Bprime-B;
          BBprime2=Bprime2-B;
%           normCB=norm(CB);
          normBBprime=norm(BBprime);
          normBBprime2=norm(BBprime2);
          
%           if normCBprime >= normCB                              %Ensures projected node is closer to the volume than original virtual node
%              Bprime=B+proj_vecB;
%           end
          
          if normBBprime > normBBprime2
              Bprime=Bprime2; 
          end
          
          planecheck=dot(planenormal,CA);
          
          if screwcheck == 0
              out=0;
          elseif planecheck == 0
              checker=zeros(size(xnodes,1),1);
          else
              pyramid=[A;B;C;Aprime;Bprime];                    %Designates original virtual nodes, closure point and projected nodes as vertices
              tess=convhulln(pyramid,{'QJ','Pp'});                          %Finds the convex hull of the designated vertices
              tol=1e1;
              checker=inhull(xnodes(:,1:3),pyramid,tess,tol);   %Flags FE nodes inside the convex hull
          end
      else
          out=0;                                           %Flags non-virtual segments
      end
            
      for j=1:nodenum
          
          nodepoint=xnodes((gnl(j)),1:3);                   %Coordinates of the FE node
          p=nodepoint';                                     %Column vector for use in displacement_et
          
          if out == 0                                      %Calculates displacement vector caused by internal segments using Barnett triangles
              Utilda(j,:)=Utilda(j,:)+displacement_et(p,A',B',C',b',nu);
          elseif out == 1                                  %Identifies virtual segments
              
              if checker(gnl(j)) == 0                            %Calculates displacement vector caused  by virtual segments and the segments connecting them to the surface using Barnett triangles
                  Utilda(j,:)=Utilda(j,:)+displacement_et(p,A',B',C',b',nu)+displacement_et(p,B',Bprime',C',b',nu)+displacement_et(p,Aprime',A',C',b',nu);
              elseif checker(gnl(j)) == 1                        %Identifies FE nodes inside the pyramidal convex hull and calcultes the correct displacement accordingly
%                   Ctemp=A+0.5*AB+10*bA;
%                   Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',nu)+displacement_et_el(p,B',Bprime',b',nu)+displacement_et_el(p,Aprime',A',b',nu)+displacement_et_plas(p,A',B',Ctemp',b')+displacement_et_plas(p,Aprime',A',Ctemp',b')+displacement_et_plas(p,B',Bprime',Ctemp',b')+displacement_et_plas(p,Bprime',Aprime',Ctemp',b')+displacement_et_plas(p,Aprime',Bprime',C',b');
%                   Utilda(j,:)=Utilda(j,:)+displacement_et(p,A',B',Ctemp',b',nu)+displacement_et(p,B',Bprime',Ctemp',b',nu)+displacement_et(p,Aprime',A',Ctemp',b',nu)+displacement_et(p,Bprime',Aprime',Ctemp',b',nu)+displacement_et(p,Aprime',Bprime',C',b',nu);
                  Utilda(j,:)=Utilda(j,:)+displacement_et_el(p,A',B',b',nu)+displacement_et_el(p,B',Bprime',b',nu)+displacement_et_el(p,Aprime',A',b',nu)+solang_correction(p,A',B',Bprime',Aprime',C',b');
              end
              
          end
          
      end
      
end

Ux=Utilda(:,1);                                             %Organises outputs
Uy=Utilda(:,2);
Uz=Utilda(:,3);

end

function [u] = displacement_et(p,A,B,C,b,nu)
% Code to calculate the displacement vector at a point p 
% due to a loop specified by a set of nodes. Loop must be 
% closed. 
% Steve Fitzgerald 16 Oct 2013
% Ed Tarleton 05 Dec 2014 -- used for debugging the solid angle/slip plane
% normal 
% This is equivalent a matlab version of displacmentmex.c but without bugs
% hopefully!

% p is field point
% ABC define triangle 
% n is slip plane normal defined by closing the original loop with RH sense
con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));


% R vectors (from P to nodes)
RA = A - p;
RB = B - p;
RC = C - p;

lamA = safenorm(RA);
lamB = safenorm(RB);
lamC = safenorm(RC);

vecAB = B - A;
vecBC = C - B;
vecCA = A - C;

  
tAB = safenorm(vecAB);
tBC = safenorm(vecBC);
tCA = safenorm(vecCA);

% calculate slip plane normal
% vec_a = segs(1,:)/norm(segs(1,:));
%     vec_b = -segs(end,:)/norm(segs(end,:));

    n = cross(tAB,tBC);
    

% f = fab+fbc+fca
f = fab(b, tAB, lamA, lamB, RA, RB) ...
    + fab(b, tBC, lamB, lamC, RB, RC) ...
    + fab(b, tCA, lamC, lamA, RC, RA);

% g = gab+ gbc + gca
g = gab(b, lamA, lamB) + gab(b, lamB, lamC) + gab(b, lamC, lamA);

omega  = solang(lamA, lamB, lamC,n);
% update displacement inc. solid angle
u =  - b*omega/(4*pi) - con1.*f + con2.*g;

u = u';

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