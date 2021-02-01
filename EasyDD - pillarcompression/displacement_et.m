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
omega = 4*atan(sqrt(temp));

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









