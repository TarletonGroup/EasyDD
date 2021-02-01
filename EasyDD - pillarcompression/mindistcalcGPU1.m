%Marielle Thibault 
%mindist_exemple is a function to improve the fist part of mindistcalc = algebraic computation

function [seg11, seg12, seg13, seg21, seg22, seg23, a, b,c,d,e,f,g]=mindistcalcGPU1(x01,x02,x03,x11,x12,x13,y01,y02,y03,y11,y12,y13)

%Not possible to deal with array
% x0=[x01 x02 x03];
% x1=[x11 x12 x13];
% y0=[y01 y02 y03];
% y1=[y11 y12 y13];

%3 coordinates of sg1
seg11=x11-x01;
seg12=x12-x02;
seg13=x13-x03;

%3 coordinates of seg2
seg21=y11-y01;
seg22=y12-y02;
seg23=y13-y03;

%expression for a
a=seg11*seg11+seg12*seg12 +seg13*seg13;
% expression for b
%xo-yo = xy
xy1=x01-y01;
xy2=x02-y02;
xy3=x03-y03;
b=2*(seg11*xy1+seg12*xy2+seg13*xy3);
%expression for c
c=2*(seg11*seg21+seg12*seg22+seg13*seg23);
%expression for d
%y0-x0
yx1=y01-x01;
yx2=y02-x02;
yx3=y03-x03;
d=2*(seg21*yx1+seg22*yx2+seg23*yx3);
%expression for e
e=seg21*seg21+seg22*seg22+seg23*seg23;
%expression for f 
f=(x01*x01+x02*x02+x03*x03)+(y01*y01+y02*y02+y03*y03);
% % %expression for g
g=(c*c)-4*a*e;


end