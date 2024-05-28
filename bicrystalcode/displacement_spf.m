function [u] = displacement_spf(p,nodes,b,nu)
% Code to calculate the displacement vector at a point p 
% due to a loop specified by a set of nodes. Loop must be 
% closed. 
% Steve Fitzgerald 16 Oct 2013

% nodes
% nodes = [10.4082   10.4082    9.1835
%      9.7011   11.1154    9.1835
%      9.2929   10.7071   10.0000
%      8.8846   10.2989   10.8165
%      9.5918    9.5918   10.8165
%      10.2989    8.8846   10.8165
%      10.7071    9.2929   10.0000
%      11.1154    9.7011    9.1835];
 
%field point
%p=[0 0 0];

% burgers' vector
%b = [0 1 0];

%nu = 0.305;
con1 = (1-2*nu)/(8*pi*(1-nu));
con2 = 1/(8*pi*(1-nu));

% segments
n = size(nodes,1);
segs = zeros(size(nodes));
for i=1:n-1
    segs(i,:) = nodes(i+1,:) - nodes(i,:);
end
segs(n,:) = nodes(1,:) - nodes(n,:);

%centroid
C = sum(nodes)./n;
%C = [10 10.01 0];

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
    lam(i,:)=R(i,:)./modR(i);
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
    u = u - b.*(solang(lamA, lamB, lamC))./(4*pi) - con1.*f + con2.*g;

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
    u = u - b.*(solang(lamA, lamB, lamC))./(4*pi) - con1.*f + con2.*g;
    %disp(u)


% finish

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

function [ omega ] = solang( lamA, lamB, lamC )
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

omega = 4*atan(sqrt(temp));



end













