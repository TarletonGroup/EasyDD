%Marielle Thibault
%logical tests from mindist function

function [L1,L2,dist2,ddist2dt] = mindistcalcGPU2(seg11, seg12, seg13, seg21, seg22, seg23,A,B,C,D,E,F,G,x01,x02,x03,x11,x12,x13,y01,y02,y03,y11,y12,y13,...
                                            vx01,vx02,vx03,vx11,vx12,vx13,vy01,vy02,vy03,vy11,vy12,vy13)

eps=1e-12;
L1=0;
L2=0;
dist2=0;
ddist2dt=0;

%vseg1=vx1-vx0;
vseg11=vx11-vx01;
vseg12=vx12-vx02;
vseg13=vx13-vx03;

%vseg2=vy1-vy0;
vseg21=vy11-vy01;
vseg22=vy12-vy02;
vseg23=vy13-vy03;

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
    %dist2=[(y0-x0)*(y0-x0)' (y1-x0)*(y1-x0)' (y0-x1)*(y0-x1)' (y1-x1)*(y1-x1)'];
    part11=y01-x01; 
    part12=y02-x02;
    part13=y03-x03;
    part21=y11-x01;
    part22=y12-x02;
    part23=y13-x03;
    part31=y01-x11;
    part32=y02-x12;
    part33=y03-x13;
    part41=y11-x11;
    part42=y12-x12;
    part43=y13-x13;
    
    dist21=(part11*part11+part12*part12+part13*part13);
    dist22=(part21*part21+part22*part22+part23*part23);
    dist23=(part31*part31+part32*part32+part33*part33);
    dist24=(part41*part41+part42*part42+part43*part43);
    
    mindist2=dist21; pos=1;
    if dist22<mindist2
        mindist2=dist22; pos=2;
    end
    if dist23<mindist2
        mindist2=dist23; pos=3;
    end
    if dist24<mindist2
        mindist2=dist24; pos=4;
    end
%   mindist2=min(dist21,dist22,dist23,dist24);
    L1=floor(pos/2);
    L2=mod(pos-1,2);

else
    L2=(2*A*D+B*C)/G;
    L1=0.5*(C*L2-B)/A;
end

%Check L1,L2
if L1>0
    L1=L1;
else
    L1=0;
end
if L1<1
    L1=L1;
else
    L1=1;
end

if L2>0
    L2=L2;
else
    L2=0;
end
if L2<1
    L2=L2;
else
    L2=1;
end

% dist2=(x0+seg1.*L1-y0-seg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)';
Long1=x01+seg11*L1-y01-seg21*L2;
Long2=x02+seg12*L1-y02-seg22*L2;
Long3=x03+seg13*L1-y03-seg23*L2;
dist2=Long1*Long1+Long2*Long2+Long3*Long3;

%ddist2dt=2*((vx0+vseg1.*L1-vy0-vseg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)');
ddist11=vx01+vseg11*L1-vy01-vseg21*L2;
ddist12=vx02+vseg12*L1-vy02-vseg22*L2;
ddist13=vx03+vseg13*L1-vy03-vseg23*L2;

ddist21=x01+seg11*L1-y01-seg21*L2;
ddist22=x02+seg12*L1-y02-seg22*L2;
ddist23=x03+seg13*L1-y03-seg23*L2;

ddist2dt=2*(ddist11*ddist21+ddist12*ddist22+ddist13*ddist23);

end

