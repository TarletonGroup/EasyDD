function [u] = detplas(p,A,B,C,b)
% Calculates the plastic (solid angle) components of displacement for
% segment AB using code from displacement_et


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