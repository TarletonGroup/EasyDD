function [dist2,ddist2dt,L1,L2,Rori]=mindistcalc(x0vx0,x1vx1,y0vy0,y1vy1)
% this function finds the minimum distance bewtween two line segments
% seg1=x0->x1 seg2=y0->y1
% dist2 = square of the minimum distance between the two points
% L1 = normalize position on seg1 that is closest to seg2
% L2 = normalized position on seg2 that is closest to seg1
% ddist2dt = time rate of change of the distance between L1 and L2
x0=x0vx0(1:3);
x1=x1vx1(1:3);
y0=y0vy0(1:3);
y1=y1vy1(1:3);
if length(x0vx0)==6
    vx0=x0vx0(4:6);
    vx1=x1vx1(4:6);
    vy0=y0vy0(4:6);
    vy1=y1vy1(4:6);
else
    vx1=zeros(1,3);
    vx0=zeros(1,3);
    vy1=zeros(1,3);
    vy0=zeros(1,3);
end

seg1=x1-x0;
seg2=y1-y0;
vseg1=vx1-vx0;
vseg2=vy1-vy0;

A=seg1*seg1';
B=2*seg1*(x0'-y0');
C=2*seg1*seg2';
D=2*seg2*(y0'-x0');
E=seg2*seg2';
F=x0*x0'+y0*y0';
G=C*C-4*A*E;


eps=1e-12;
if A<eps % seg1 is just a point
    L1=0;
    if E<eps
        L2=0;
    else
        L2=-0.5*D/E;
    end
elseif E<eps % seg2 is just a point
    L2=0;
    if A<eps
        L1=0;
    else
        L1=-0.5*B/A;
    end
elseif abs(G)<eps % lines are parallel
    dist2=[(y0-x0)*(y0-x0)' (y1-x0)*(y1-x0)' (y0-x1)*(y0-x1)' (y1-x1)*(y1-x1)'];
    [mindist2,pos]=min(dist2);
    L1=floor(pos/2);
    L2=mod(pos-1,2);
else
    L2=(2*A*D+B*C)/G;
    L1=0.5*(C*L2-B)/A;
end

% now check to make sure that L2 and L1 are betwen 0 and 1
L1=min(max([L1,0]),1);
L2=min(max([L2,0]),1);

% now calculate the distance^2 and the time rate of change of the distance between the points at L1 and L2
dist2=(x0+seg1.*L1-y0-seg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)';
ddist2dt=2*((vx0+vseg1.*L1-vy0-vseg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)');