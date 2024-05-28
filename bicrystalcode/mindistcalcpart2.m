%Marielle Thibault 
%Second part of mindistcalc function 

function [dist2,ddist2dt,L1,L2]=mindistcalcpart2(A,B,C,D,E,F,G)

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
% dist2=(x0+seg1.*L1-y0-seg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)';
% ddist2dt=2*((vx0+vseg1.*L1-vy0-vseg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)');

end