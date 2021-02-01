function [u] = displacement_bb(p,A,B,cm,b,nu)
% Code to calculate displacement vector contribution at point p due to a
% segmaent between A and B. Solid angle contribution is measured by forming
% a triangle with point cm. Based on Hirth & Lothe, 1982, and
% displacement_et.m, 2014.
%N.B. This is not a physical diplacement vector. Actual displacement can
%only be found by summing the contributions of all segments in a
%dislocation network compromised of one or mor closed loops.
%N.B.2 There are inconsitencies in notation between H&L and the work by
%Barnett on which displacement_et.m is primarily based. Variables here have
%been named as is convenient.
%Bruce Bromage 26 June 2017

%constants in elastic equations
con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));

% R vectors (from p to nodes)
RA = A - p;
RB = B - p;
Rcm = cm - p;  %R' in H&L

%Normalised R vectors
lamA = safenorm(RA);
lamB = safenorm(RB);
LA = A-cm;  %lamA in H&L 3rd edition p 132
LB = B-cm;  %lamB in H&L

%Line vector and line direction
vecAB = B - A;
tAB = safenorm(vecAB);

% calculate ABcm normal
n = cross(LA, tAB);
n = safenorm(n);

%calculate height of ABcm
eAB = cross(n, tAB);
eAB = safenorm(eAB);
hAB = dot(-LA, eAB);


%calculate perpendicular distance from ABcm to p
Z = -dot(Rcm, n);

%calculate equation terms
f = fab(b, tAB, lamA, lamB, RA, RB);
g = gab(b, lamA, lamB);
omega = solang2(LA, LB, hAB, Z);

%combine terms for total displacement contribution
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

function [ omega ] = solang2(LA, LB, hAB, Z)

lA=norm(LA);  %magnitude of LA
lB=norm(LB);  %magnitude of LB
Acm=LA/lA;  %direction of LA
Bcm=LB/lB;  %direction of LB
phiAB=acos(dot(Acm,Bcm));  %angle in ABcm at cm
sgn=sign(Z);  %sign of Z

term1=(phiAB-pi)*sgn;

term2=atan2((hAB)*sqrt(lA^2+Z^2),(Z)*sqrt(lA^2-hAB^2));

term3=atan2((hAB)*sqrt(lB^2+Z^2),(Z)*sqrt(lB^2-hAB^2)); %check atan2 and atan

if abs(hAB)<eps && abs(Z)<eps
    term2=0;
    term3=0;
end

omega=term1+term2+term3;
omega=atan(tan(omega));

end

function unitR = safenorm(R)

normR = sqrt(R(1)*R(1) + R(2)*R(2)+R(3)*R(3));

if normR > eps
    unitR = R/normR;
else
    unitR = zeros(3,1);
end

end
