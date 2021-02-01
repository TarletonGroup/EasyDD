function [utilda] = displacement(p,segments,loop_list,nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to calculate the displacement vector at a point p 
% ue to a dislocation structure (closed loops / networks)
%
% Steve Fitzgerald, 16 Oct 2013
% Francesco Ferroni, 28 Oct 2013
%
% University of Oxford, Department of Materials, Parks Road, United Kingdom
% Defects Group & Materials for Fission and Fusion Power
%
% p: field point [x y z]
% u: displacement [ux uy uz]
% tol: tolerance to check whether node in loop is coplanar (normal distance
% to plane, in units of b). 5-10b would be acceptable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

utilda = [0 0 0]; %initialize
con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));

n_loops = max(loop_list(:,2)); 
size_list = size(loop_list,1);
size_seglist = size(segments,1);

for z=1:n_loops
    
    %isolate single "loop" from "loop_list"
    index=false(size_list,1);
    for d=1:size_list
        if z==loop_list(d,2);
            index(d)=true;
        end
    end
    loop=loop_list(index,1);
    n=size(loop,1);
       
    %generate node pairs in loop
    id_pairs = zeros(n,2);
    for i=1:n-1
        id_pairs(i,1) = loop(i);
        id_pairs(i,2) = loop(i+1);
    end
    id_pairs(n,1) = loop(n);
    id_pairs(n,2) = loop(1);
    
    index=false(size_seglist,1);
    for d=1:size_seglist
        for i=1:n 
            if id_pairs(i,1)==segments(d,1) && id_pairs(i,2)==segments(d,2)
                index(d)=true;
            elseif id_pairs(i,1)==segments(d,2) && id_pairs(i,2)==segments(d,1)
                index(d)=true;
            end
        end
    end
    %make 'nodes' and 'segs' list
    nodes=segments(index,6:8);
    b=segments(index,3:5);
    segs=segments(index,9:11)-segments(index,6:8);
    
    %check burgers consistency. they should be all the same since they are
    %in a single loop, I think.
    for d=1:n
        if not ( all(b(:,1) == b(1,1)) && all(b(:,2) == b(1,2)) && all(b(:,3) == b(1,3)) )
            disp('Not-identical Burgers vector for loop number ',num2str(z))
            pause;
        end
    end
    b=b(1,:);
    
    %calculate loop plane normal [a b c d], (using 2 segments of loop)
    vec_a = segs(1,:)/norm(segs(1,:));
    vec_b = -segs(end,:)/norm(segs(end,:));
    plane_n = cross(vec_a,vec_b);
    d = -dot(nodes(1,:),plane_n);
    plane_n(4) = d;
    
%     %check if loop nodes are coplanar. not necessary: no cross-slip.
%     cros
%     vrts=rn(loop(:),1:3);
%     coplanar = vrts(:,1)*plane_n(1) + vrts(:,2)*plane_n(2) + vrts(:,3)*plane_n(3) + plane_n(4);
%     test=any(coplanar>tol | coplanar<-tol);
%     if test
%         disp('Loop is not coplanar')
%     end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPF

    %centroid
    C = sum(nodes)./n;

    % t vectors (unit tangents to segments)
    t = zeros(size(nodes));
    for i=1:n
        t(i,:)=segs(i,:)./norm(segs(i,:));
    end

    % R vectors (from P to nodes)
    R = zeros(size(nodes));
    for i=1:n
        R(i,:) = nodes(i,:) - p;
    end
    RC = C - p;

    % absolute values
    modR = zeros(n,1);
    for i=1:n
        modR(i) = norm(R(i,:));
    end
    modRC = norm(RC);

    % lambda vectors (unit vectors of the R's)
    lam = zeros(size(nodes));
    for i=1:n
        lam(i,:)=R(i,:)./modR(i); %if point p lies on the lam vector, this creates NaN term! - FF
    end
    lamC = RC./modRC;

    % extra segments FROM nodes TO centroid
    extras = zeros(size(nodes));
    for i=1:n
        extras(i,:) = C - nodes(i,:);
    end

    % tangents for these -- FROM nodes TO centroid
    tex = zeros(size(nodes));
    for i=1:n
        tex(i,:) = extras(i,:)./norm(extras(i,:));
    end

    % displacement vector
    u = [0 0 0];

    % start from first node, make a triangle ABC

    for i=1:n-1;
        ab = t(i,:);
        bc = tex(i+1,:);
        ca = -tex(i,:);%opposite sense
        lamA = lam(i,:);
        lamB = lam(i+1,:);
        % already know lamC
        RA = modR(i);
        RB = modR(i+1);
        % already know modRC

        % get the three f's
        f = fab(b, ab, lamA, lamB, RA, RB) ...
            + fab(b, bc, lamB, lamC, RB, modRC) ...
            + fab(b, ca, lamC, lamA, modRC, RA);

        % and the three g's
        g = gab(b, lamA, lamB) + gab(b, lamB, lamC) + gab(b, lamC, lamA);

        % update displacement inc. solid angle
        u = u - b.*(solang(lamA, lamB, lamC , p , plane_n))./(4*pi) - con1.*f + con2.*g;
     end

        % last one
        ab = t(n,:);
        bc = tex(1,:);
        ca = -tex(n,:);%opposite sense
        lamA = lam(n,:);
        lamB = lam(1,:);
        % already know lamC
        RA = modR(n);
        RB = modR(1);
        % already know modRC

        % get the three f's
        f = fab(b, ab, lamA, lamB, RA, RB) ...
            + fab(b, bc, lamB, lamC, RB, modRC) ...
            + fab(b, ca, lamC, lamA, modRC, RA);

        % and the three g's
        g = gab(b, lamA, lamB) + gab(b, lamB, lamC) + gab(b, lamC, lamA);

        % update displacement inc. solid angle
        u = u - b.*(solang(lamA, lamB, lamC , p , plane_n))./(4*pi) - con1.*f + con2.*g;
        
        % add onto final displacement
        utilda = utilda + u;
        
         %disp('Oosterom Method (new)');
         %disp(solang2(lamA, lamB, lamC , p , plane_n));
       
         %disp('Huilier Method (original)');
         %disp(solang(lamA, lamB, lamC , p , plane_n));
end

utilda = utilda'; %make it a column vector

end

function [ out ] = fab( b, t, lamA, lamB, RA, RB )
%calculates f vector

numerator = RB*(1 + dot(lamB,t));
denominator = RA*(1 + dot(lamA,t));

logarithm = log(numerator/denominator);

out = cross(b,t).*logarithm;

end

function [ out ] = gab( b, lamA, lamB )
%calculates g vector

numerator = dot(b,cross(lamA,lamB)).*(lamA + lamB);
denominator = 1 + dot(lamA,lamB);

out = numerator./denominator;

end

function [ omega ] = solang( lamA, lamB, lamC , p , plane_n )
%calculates solid angle due to triangular loop ABC
% 1) this may not be robust, check different formulae for the 
% spherical excess (wiki spherical geometry)
% 2) it could be a lot faster to calculate a,b,c in the main loop 
% for each pair of lambdas, then use a,b,c as input arguments 
a = acos(dot(lamB,lamC));
b = acos(dot(lamC,lamA));
c = acos(dot(lamA,lamB));

s = (a + b + c)/2;

svec = [s s-a s-b s-c]./2;

temp = tan(svec(1))*tan(svec(2))*tan(svec(3))*tan(svec(4));

omega = 4*atan(sqrt(abs(temp)));

if dot(lamA,plane_n(1:3)) > 0
    sign = 1;
else
    sign = -1;
end

omega = -sign*omega;

%Eq(18) Barnett 1985
%test whether point is above or below plane or on it.
% x = dot(p,plane_n(1:3)) + plane_n(4);
% if x < 0
%     omega = -omega;
% end

end

function [ omega ] = solang2( lamA, lamB, lamC , p , plane_n )
%alternative solid angle formula to see if edge dislocation miscalculation
%is stemming from some error in solang function.

%given 3 input vectors: lamA, lamB, lamC

determ = det([lamA;lamB;lamC]);
 
al = norm(lamA);
bl = norm(lamB);
cl = norm(lamC);

div = al*bl*cl + dot(lamA,lamB)*cl + dot(lamA,lamC)*bl + dot(lamB,lamC)*al;
at = atan2(determ, div);

if at<0
    at = at + pi;
end

omega = 2*at;

if dot(lamA,plane_n(1:3)) > 0
    sign = 1;
else
    sign = -1;
end

omega = -sign*omega;

%Eq(18) Barnett 1985
%test whether point is above or below plane
% x = dot(p,plane_n(1:3)) + plane_n(4);
% if x < 0
%     omega = -omega;
% end

end

