function [fx3,fx4,fx5,fx6, ftot]=NodalSurfForceLinearRectangle2(x1,x2,x3,x4,x5,x6,b,mu,nu,a,t)
% This routine calculates the Nodal Forces induced by a non-singular
% straigth segment of dislocation on a linear rectangular surface element.
% Inputs:x1(nseg,1:3) a list of starting positions of the dislocations segments
%        x2(nseg,1:3) a list of ending positions of the dislocation segments
%        x3-x6 Cartesian coordinates of the rectangular element
%        b(nseg,1:3) a list of burgers vectors associated with the segments going from x1 to x2
%        a=core width
%        mu= shear modulus
%        nu=poisson's ratio
%
% Notations and details can be found in S. Queyreau, J. Marian, B.D. Wirth,
% A. Arsenlis, MSMSE, 22(3):035004, (2014).
%
% Sylvain Queyreau 05/2014

% p and q are unit vectors
p= x4-x3 ;
L=norm(p);
p= p./L; %unit vector along one side of the rectawhyngular surface element
q= x5-x3;
Lprime=norm(q);
q=q./Lprime; %unit vector along the other side of the surface element

t=x2-x1;
t=t./norm(t);
%test0=x6-x3-(x5-x3+x4-x3) % this test passed to be complemented if not passed
n=cross(p,q);
pxt=cross(p,t);
qxt=cross(q,t);

% Warning: this formulation assumes x3-x6 to be diagonal

% r and s are two scalars that span the surface element
% y spans the dislocation segment
Rlim=[x3-x1;x6-x2]; %distance vector between the segment and the surface
y1=dot(Rlim(1,1:3),n)/dot(t,n);
y2=dot(Rlim(2,1:3),n)/dot(t,n);
r1=dot(Rlim(1,1:3),qxt)/dot(p,qxt);
r2=dot(Rlim(2,1:3),qxt)/dot(p,qxt);
s1=dot(Rlim(1,1:3),pxt)/dot(q,pxt);
s2=dot(Rlim(2,1:3),pxt)/dot(q,pxt);

rv= [r2;r2;r2;r2;r1;r1;r1;r1];
sv= [s2;s2;s1;s1;s2;s2;s1;s1];
yv= [y2;y1;y2;y1;y2;y1;y2;y1];

signv= [1 -1 -1 1 -1 1 1 -1];

pdott= dot(p,t);
qdott= dot(q,t);

R=zeros(8,3);
i=1:8;
j=1:3;
R(i,j)= rv(i)*p(j)+sv(i)*q(j)+yv(i)*t(j);

Ra=zeros(1,8);
i=1:8;
Ra(i)=((rv(i).*rv(i) + sv(i).*sv(i) + yv(i).*yv(i) + 2*pdott*rv(i).*yv(i) + 2*qdott*sv(i).*yv(i) + a^2).^0.5);
Rdotp= R(i,1)*p(1)+R(i,2)*p(2)+R(i,3)*p(3);
Rdotq= R(i,1)*q(1)+R(i,2)*q(2)+R(i,3)*q(3);
Rdott= R(i,1)*t(1)+R(i,2)*t(2)+R(i,3)*t(3);

% LINEAR SEED INTEGRALS
A0m1= log(Ra' + Rdotp) ; % m stands for '-'
B0m1= log(Ra' + Rdotq) ;
C0m1= log(Ra' + Rdott) ;

% LINEAR INTEGRALS
%A_{il} = \int r^i Ra^l dr % m stands for '-'
%B_{jl} = \int s^j Ra^l dr
%C_{kl} = \int y^k Ra^l dr
A1m1= Ra' - pdott.*yv.*A0m1;
B1m1= Ra' - qdott*yv.*B0m1;
C1m1= Ra' - (pdott*rv+qdott*sv).*C0m1;

%test1=[A0m1'*signv'; B0m1'*signv'; C0m1'*signv'; A1m1'*signv';B1m1'*signv'; C1m1'*signv'] % passed

A01= 0.5*(Rdotp.*Ra' + (Ra'.^2-(Rdotp.^2)).*A0m1);
B01= 0.5*(Rdotq.*Ra' + (Ra'.^2-(Rdotq.^2)).*B0m1);
C01= 0.5*(Rdott.*Ra' + (Ra'.^2-(Rdott.^2)).*C0m1);

A11= 1/3*Ra'.^3 - pdott*yv.*A01 ;
B11= 1/3*Ra'.^3 - qdott*yv.*B01 ;
%C11= 1/3*Ra'.^3 - (pdott*rv+qdott*sv).*C01 ;

%test2=[A01'*signv'; B01'*signv'; C01'*signv'; A11'*signv'; B11'*signv'; C11'*signv'] % passed

% DOUBLE SEED INTEGRALS:
udotu=(1-pdott^2)  ;
vdotv=(1-qdott^2)  ;
D00m3=[0 0 0 0 0 0 0 0]; E00m3=[0 0 0 0 0 0 0 0]; F00m3=[0 0 0 0 0 0 0 0];
for i=1:8
    if (a^2+(1-qdott^2-pdott^2)*yv(i)^2 >0)
        temp1= (a^2+(1-qdott^2-pdott^2)*yv(i)^2)^0.5 ;
        D00m3(i) = (2/temp1 * atan((Ra(i) - Rdotp(i) + Rdotq(i))/temp1)); % m stands for '-'
    elseif (a^2+(1-qdott^2-pdott^2)*yv(i)^2 < 0)
        temp1= (abs(a^2+(1-qdott^2-pdott^2)*yv(i)^2))^0.5 ;
        D00m3(i)= (-2/temp1 * atanh((Ra(i) - Rdotp(i) + Rdotq(i))/temp1));
    end
end
for i=1:8
    if ((1-pdott^2)*a^2+(1-qdott^2-pdott^2)*sv(i)^2 >0)
        temp1=((1-pdott^2)*a^2+(1-qdott^2-pdott^2)*sv(i)^2)^0.5 ;
        E00m3(i) = (2/temp1 * atan(((1-pdott)*(Ra(i)-rv(i)+yv(i))+qdott*sv(i))/temp1));
    elseif ((1-pdott^2)*a^2+(1-qdott^2-pdott^2)*sv(i)^2 < 0)
        temp1= (abs((1-pdott^2)*a^2+(1-qdott^2-pdott^2)*sv(i)^2))^0.5 ;
        E00m3(i) = (-2/temp1 * atanh(((1-pdott)*(Ra(i)-rv(i)+yv(i))+qdott*sv(i))/temp1));
    end
end
for i=1:8
    if ((1-qdott^2)*a^2+(1-qdott^2-pdott^2)*rv(i)^2 >0)
        temp1=((1-qdott^2)*a^2+(1-qdott^2-pdott^2)*rv(i)^2)^0.5 ;
        F00m3(i) = (2/temp1 * atan(((1-qdott)*(Ra(i)-sv(i)+yv(i))+pdott*rv(i))/temp1));
    elseif ((1-qdott^2)*a^2+(1-qdott^2-pdott^2)*rv(i)^2 < 0)
        temp1= (abs((1-qdott^2)*a^2+(1-qdott^2-pdott^2)*rv(i)^2))^0.5 ;
        F00m3(i) = (-2/temp1 * atanh(((1-qdott)*(Ra(i)-sv(i)+yv(i))+pdott*rv(i))/temp1)) ;
    end
end

%test3=[D00m3*signv'; E00m3*signv'; F00m3*signv']

% DOUBLE INTEGRALS:
% D_{ijl} = \iint r^i s^j Ra^l dr ds % m stands for '-'
% E_{ikl} = \iint r^i y^k Ra^l dr dy
% F_{jkl} = \iint s^j y^k Ra^l ds dy
D01m3= -A0m1 - qdott*yv.*D00m3' ;
D10m3= -B0m1 - pdott*yv.*D00m3' ;
E01m3= -1/udotu *(A0m1 - pdott*C0m1 + qdott*sv.*E00m3') ;
F01m3= -1/vdotv *(B0m1 - qdott*C0m1 + pdott*rv.*F00m3') ;
E10m3= -C0m1-pdott*E01m3;
F10m3= -C0m1-qdott*F01m3;

%test4=[D01m3'*signv'; D10m3'*signv'; E01m3'*signv'; F01m3'*signv'; E10m3'*signv'; F10m3'*signv'] %limit

D00m1= rv.*B0m1 + sv.*A0m1 - pdott*yv.*D10m3 - qdott*yv.*D01m3 - (a^2+yv.^2).*D00m3' ;
E00m1= rv.*C0m1 + yv.*A0m1 - qdott*sv.*E01m3 - (a^2+sv.^2).*E00m3' ;
F00m1= sv.*C0m1 + yv.*B0m1 - pdott*rv.*F01m3 - (a^2+rv.^2).*F00m3' ;
D11m3= -A1m1 - qdott*yv.*D10m3;
E11m3= 1/udotu *(-A1m1 + pdott*rv.*C0m1 - pdott*E00m1-qdott*sv.*E10m3);
F11m3= 1/vdotv *(-B1m1 + qdott*sv.*C0m1 - qdott*F00m1-pdott*rv.*F10m3);

%test5=[D00m1'*signv'; E00m1'*signv'; F00m1'*signv'; D11m3'*signv'; E11m3'*signv'; F11m3'*signv']

D20m3= D00m1 - pdott*yv.*D10m3 - rv.*B0m1;
D02m3= D00m1 - qdott*yv.*D01m3 - sv.*A0m1;
E20m3= E00m1 - pdott*E11m3 - rv.*C0m1;
F20m3= F00m1 - qdott*F11m3 - sv.*C0m1;
%E02m3= E00m1 - pdott*E11m3 - qdott*sv.*E01m3 - yv.*A1m1;
%F02m3= F00m1 - qdott*F11m3 - pdott*rv.*F01m3 - yv.*B1m1;

%test6=[D02m3'*signv'; D20m3'*signv'; E20m3'*signv'; F20m3'*signv';E02m3'*signv'; F02m3'*signv']

D01m1= A01 - qdott*yv.*D00m1;
D10m1= B01 - pdott*yv.*D00m1;
E01m1= 1/udotu * (A01 - pdott*C01 - qdott*sv.*E00m1) ;
F01m1= 1/vdotv * (B01 - qdott*C01 - pdott*rv.*F00m1) ;
E10m1= C01 - pdott*E01m1 ;
F10m1= C01 - qdott*F01m1 ;

%test7=[D10m1'*signv'; D01m1'*signv'; E01m1'*signv'; F01m1'*signv';E10m1'*signv'; F10m1'*signv'] %passed

D001= 1/3*( (yv.^2+a^2).*D00m1 + rv.*B01 + sv.*A01 + yv.*pdott.*D10m1 +yv.*qdott.*D01m1 ) ;
E001= 1/3*( (sv.^2+a^2).*E00m1 + rv.*C01 + yv.*A01 + qdott.*sv.*E01m1);
F001= 1/3*( (rv.^2+a^2).*F00m1 + sv.*C01 + yv.*B01 + pdott.*rv.*F01m1);

%test8=[D001'*signv'; E001'*signv'; F001'*signv'] %passed

D11m1=A11-qdott*yv.*D10m1;
E11m1=(1/udotu)*(A11-pdott*rv.*C01 - qdott*sv.*E10m1+pdott*E001);
F11m1=(1/vdotv)*(B11-qdott*sv.*C01 - pdott*rv.*F10m1+qdott*F001);
D02m1=sv.*A01-D001-qdott*yv.*D01m1;
D20m1=rv.*B01-D001-pdott*yv.*D10m1;
E20m1=rv.*C01-E001 -pdott*E11m1;
F20m1=sv.*C01-F001 -qdott*F11m1;

%test9=[D11m1'*signv'; E11m1'*signv'; F11m1'*signv'; D20m1'*signv'; D02m1'*signv'; E20m1'*signv'; F20m1'*signv'] %passed

D12m3= D10m1 - qdott*yv.*D11m3 - sv.*A1m1 ;
D21m3= D01m1 - pdott*yv.*D11m3 - rv.*B1m1 ;
E21m3= 1/udotu * (E01m1 - pdott*E10m1 - rv.*C1m1 + pdott*yv.*A1m1 + qdott*pdott*sv.*E11m3) ;
F21m3= 1/vdotv * (F01m1 - qdott*F10m1 - sv.*C1m1 + qdott*yv.*B1m1 + qdott*pdott*rv.*F11m3) ;
D22m3=-D001 - qdott*yv.*D01m1 - pdott*yv.*D10m1 + pdott*qdott*yv.*yv.*D11m3-rv.*sv.*Ra' + rv.*B01+sv.*A01 + qdott*rv.*yv.*B1m1 + pdott*sv.*yv.*A1m1;

%test10=[D12m3'*signv'; D21m3'*signv'; E21m3'*signv'; F21m3'*signv'; D22m3'*signv']

% TRIPLE INTEGRALS
% H_{ijkl} = \iiint r^i s^j y^k Ra^l drdsdy
% m stands for '-'
H001m3= -1 * (D00m1 - qdott*E00m1 - pdott*F00m1) /(1-pdott^2-qdott^2);
H100m3= -pdott*H001m3 - F00m1;
H010m3= -qdott*H001m3 - E00m1;

H001m1= 1/(1-pdott^2-qdott^2) * (D001 - pdott*F001 - qdott*E001);
H100m1= F001 - pdott*H001m1;
H010m1= E001 - qdott*H001m1;

H111m3=(1/(1-pdott^2-qdott^2))*(pdott*rv.*F10m1 + qdott*sv.*E10m1 - pdott*H010m1 - qdott*H100m1-D11m1);
H210m3= H010m1 - rv.*F10m1 - pdott*H111m3;
H120m3= H100m1 - sv.*E10m1 - qdott*H111m3;
H021m3=1/(1-pdott^2-qdott^2)*(-2*qdott*H010m1-D02m1+qdott*sv.*sv.*E00m1+pdott*F20m1);
H201m3=1/(1-pdott^2-qdott^2)*(-2*pdott*H100m1-D20m1+pdott*rv.*rv.*F00m1+qdott*E20m1);

%H000m3 is the only integral for which there is no analytical
%antiderivative. This integral could be very easily estimated by a Gauss
%quadrature. This is not actually required at all as all the H000m3 terms cancel
%out exactly in the final nodal force calculation. H000m3 is set to 0 for
%convenience.
H000m3= [0 0 0 0 0 0 0 0]';
%This means that all antiderivatives in relation with H000m3 are
%modified and their value would disagree with a Gauss quadrature if H000m3
%is not reset to its regular value. However the final sum and nodal forces
%remain unchanged.

H000m1= 1/2*(rv.*F00m1 + sv.*E00m1 + yv.*D00m1 - a^2*H000m3);

AA = 1-pdott*pdott-qdott*qdott ;

H101m3= (-D10m1+pdott*rv.*F00m1+qdott*E10m1-pdott*H000m1)/AA ;
H011m3= (-D01m1+qdott*sv.*E00m1+pdott*F10m1-qdott*H000m1)/AA ;
H110m3= -E10m1 - qdott*H101m3;
H200m3= H000m1 - rv.*F00m1 - pdott*H101m3;
H020m3= H000m1 - sv.*E00m1 - qdott*H011m3;
%H002m3= H000m1 - yv.*D00m1 - pdott*H101m3 - qdott*H011m3;

H001m5=  (pdott*F00m3 + qdott*E00m3 - D00m3) /3/AA ;
H100m5= -1/3*F00m3 - pdott*H001m5;
H010m5= -1/3*E00m3 - qdott*H001m5;

H101m5= (-pdott*H000m3-D10m3+pdott*rv.*F00m3'+qdott*E10m3) /3/AA;
H011m5= (-qdott*H000m3-D01m3+qdott*sv.*E00m3'+pdott*F10m3) /3/AA;
H110m5= -1/3*F10m3 - pdott*H011m5;
H200m5= 1/3*H000m3 - rv/3.*F00m3' - pdott*H101m5;
H020m5= 1/3*H000m3 - sv/3.*E00m3' - qdott*H011m5;

H111m5=  (-D11m3 + pdott*rv.*F10m3 + qdott*sv.*E10m3 - pdott*H010m3 - qdott*H100m3) /3/AA ;
H210m5= -rv/3 .*F10m3 + 1/3*H010m3 - pdott*H111m5;
H120m5= -sv/3 .*E10m3 + 1/3*H100m3 - qdott*H111m5;
H201m5= (-2*pdott*H100m3-D20m3+pdott*rv.*rv.*F00m3'+qdott*E20m3)/3/AA;
H021m5= (-2*qdott*H010m3-D02m3+qdott*sv.*sv.*E00m3'+pdott*F20m3)/3/AA ;
H012m5= (H010m3 -qdott*H001m3-yv.*D01m3+qdott*sv.*E01m3+pdott*F11m3)/3/AA;
H102m5= (H100m3 -pdott*H001m3-yv.*D10m3+pdott*rv.*F01m3+qdott*E11m3)/3/AA;

H112m5= (-yv.*D11m3 + pdott*rv.*F11m3 + qdott*sv.*E11m3 + H110m3 - pdott*H011m3 - qdott*H101m3) /3/AA;
H211m5= -rv/3.*F11m3 + 1/3*H011m3 - pdott*H112m5;
H121m5= -sv/3.*E11m3 + 1/3*H101m3 - qdott*H112m5;
H202m5= (H200m3-2*pdott*H101m3-yv.*D20m3+pdott*rv.*rv.*F01m3+qdott*E21m3) /3/AA;
H022m5= (H020m3-2*qdott*H011m3-yv.*D02m3+qdott*sv.*sv.*E01m3+pdott*F21m3) /3/AA;
H301m5= -1/3*rv.^2.*F01m3 + 2/3*H101m3 - pdott*H202m5;
H031m5= -1/3*sv.^2.*E01m3 + 2/3*H011m3 - qdott*H022m5;

H122m5= (H120m3-pdott*H021m3-2*qdott*H111m3-yv.*D12m3+pdott*rv.*F21m3+qdott*sv.*sv.*E11m3) /3/AA;
H212m5= (H210m3-qdott*H201m3-2*pdott*H111m3-yv.*D21m3+qdott*sv.*E21m3+pdott*rv.*rv.*F11m3) /3/AA;
H221m5= (pdott*rv.*rv.*F20m3+qdott*sv.*sv.*E20m3-2*pdott*H120m3-2*qdott*H210m3-D22m3) /3/AA;
H311m5= (2*(1-qdott*qdott)*H111m3+qdott*pdott*H201m3-pdott*H210m3...
    -(1-qdott*qdott)*rv.*rv.*F11m3-pdott*qdott*sv.*E21m3+pdott*yv.*D21m3) /3/AA;
H131m5= (2*(1-pdott*pdott)*H111m3+qdott*pdott*H021m3-qdott*H120m3...
    -(1-pdott*pdott)*sv.*sv.*E11m3-pdott*qdott*rv.*F21m3+qdott*yv.*D12m3)/3/AA;

%SCALARIZATION
%Every tripple integrals are now evaluated.
%if terms are not functions of r,s and y they disappear in the summation
%only functions of the three variable matter.
scH001m3= dot(H001m3,signv);
scH100m3= dot(H100m3,signv);
scH010m3= dot(H010m3,signv);
%scH001m1= dot(H001m1,signv); % intermediate integral, scalar evaluation not required
%scH100m1= dot(H100m1,signv); % intermediate integral, scalar evaluation not required
%scH010m1= dot(H010m1,signv); % intermediate integral, scalar evaluation not required
scH111m3= dot(H111m3,signv);
scH210m3= dot(H210m3,signv);
scH120m3= dot(H120m3,signv);
%scH021m3= dot(H021m3,signv); % intermediate integral, scalar evaluation not required
%scH201m3= dot(H201m3,signv); % intermediate integral, scalar evaluation not required
%scH000m3= dot(H000m3,signv); % intermediate integral, scalar evaluation not required
%scH000m1= dot(H000m1,signv); % intermediate integral, scalar evaluation not required
scH101m3= dot(H101m3,signv);
scH110m3= dot(H110m3,signv);
scH011m3= dot(H011m3,signv);
scH200m3= dot(H200m3,signv);
scH020m3= dot(H020m3,signv);
%scH002m3= dot(H002m3,signv); % intermediate integral, scalar evaluation not required
scH001m5= dot(H001m5,signv);
scH100m5= dot(H100m5,signv);
scH010m5= dot(H010m5,signv);
scH101m5= dot(H101m5,signv);
scH011m5= dot(H011m5,signv);
scH110m5= dot(H110m5,signv);
scH200m5= dot(H200m5,signv);
scH020m5= dot(H020m5,signv);
scH111m5= dot(H111m5,signv);
scH210m5= dot(H210m5,signv);
scH120m5= dot(H120m5,signv);
scH021m5= dot(H021m5,signv);
scH201m5= dot(H201m5,signv);
scH102m5= dot(H102m5,signv);
scH012m5= dot(H012m5,signv);
scH112m5= dot(H112m5,signv);
scH211m5= dot(H211m5,signv);
scH121m5= dot(H121m5,signv);
scH202m5= dot(H202m5,signv);
scH022m5= dot(H022m5,signv);
scH301m5= dot(H301m5,signv);
scH031m5= dot(H031m5,signv);
scH122m5= dot(H122m5,signv);
scH221m5= dot(H221m5,signv);
scH131m5= dot(H131m5,signv);
scH311m5= dot(H311m5,signv);
scH212m5= dot(H212m5,signv);

% FINAL SUM
% the force is now calculated for each of the 4 patch points
tcrossb= cross(t,b);
pcrossb= cross(p,b);
qcrossb= cross(q,b);
bcrosst= cross(b,t);
tdotn= dot(t,n);
ttbn= t*dot(tcrossb,n);
tbtn= t*dot(bcrosst,n);
tpbn= t*dot(pcrossb,n);
tqbn= t*dot(qcrossb,n);
factor=mu/4/pi/(1-nu)/L/Lprime;
a2=a*a;

% Force in x3
ri= r2;
si= s2;
% Vectors (identical for all xi nodes of the surface element)
I111m3=(tcrossb*tdotn+t*(dot(tcrossb,n)))*(1-nu) + (bcrosst*tdotn+tbtn);
I120m3=(qcrossb*tdotn+tqbn)*(1-nu)-(dot(qcrossb,t)*n)+(q*(dot(bcrosst,n)));
I210m3=(pcrossb*tdotn+tpbn)*(1-nu)-(dot(pcrossb,t)*n)+(p*(dot(bcrosst,n)));
I111m5=(tcrossb*tdotn+ttbn)*1.5*(1-nu)*a2;
I120m5=(qcrossb*tdotn+tqbn)*1.5*(1-nu)*a2-(dot(qcrossb,t)*n) * 3*a2;
I210m5=(pcrossb*tdotn+tpbn)*1.5*(1-nu)*a2-(dot(pcrossb,t)*n) * 3*a2;
I131m5=-(dot(qcrossb,t)*(tdotn))*q * 3;
I311m5=-(dot(pcrossb,t)*(tdotn))*p * 3;
I122m5=-(dot(qcrossb,t)*tdotn)*t * 3;
I212m5=-(dot(pcrossb,t)*tdotn)*t * 3;
I221m5=-((dot(pcrossb,t)*tdotn)*q +(dot(qcrossb,t)*tdotn)*p)*3;

% Final integrals required for x3
F111m3=scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3;
F120m3=scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3;
F210m3=scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3;
F111m5=scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5;
F120m5=scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5;
F210m5=scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5;
F131m5=scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5;
F311m5=scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5;
F122m5=scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5;
F212m5=scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5;
F221m5=scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5;

fLLprime= I111m3 * F111m3 + I210m3 * F210m3 + I120m3 * F120m3 ...
    + I111m5 * F111m5 + I210m5 * F210m5 + I120m5 * F120m5 ...
    + I212m5 * F212m5 + I122m5 * F122m5 + I221m5 * F221m5 ...
    + I311m5 * F311m5 + I131m5 * F131m5;

fx3= fLLprime*factor;

% Force in x4
ri= r1;
si= s2;
% Integrals required for x4
F111m3=-(scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3);
F120m3=-(scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3);
F210m3=-(scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3);
F111m5=-(scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5);
F120m5=-(scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5);
F210m5=-(scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5);
F131m5=-(scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5);
F311m5=-(scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5);
F122m5=-(scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5);
F212m5=-(scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5);
F221m5=-(scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5);

fLLprime= I111m3 * F111m3 + I210m3 * F210m3 + I120m3 * F120m3 ...
    + I111m5 * F111m5 + I210m5 * F210m5 + I120m5 * F120m5 ...
    + I212m5 * F212m5 + I122m5 * F122m5 + I221m5 * F221m5 ...
    + I311m5 * F311m5 + I131m5 * F131m5;

fx4= fLLprime*factor;

%Force in x5
ri= r2;
si= s1;
% Integrals required for x5
F111m3=-(scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3);
F120m3=-(scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3);
F210m3=-(scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3);
F111m5=-(scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5);
F120m5=-(scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5);
F210m5=-(scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5);
F131m5=-(scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5);
F311m5=-(scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5);
F122m5=-(scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5);
F212m5=-(scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5);
F221m5=-(scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5);

fLLprime= I111m3 * F111m3 + I210m3 * F210m3 + I120m3 * F120m3 ...
    + I111m5 * F111m5 + I210m5 * F210m5 + I120m5 * F120m5 ...
    + I212m5 * F212m5 + I122m5 * F122m5 + I221m5 * F221m5 ...
    + I311m5 * F311m5 + I131m5 * F131m5;

fx5= fLLprime*factor;

%Force in x6
ri= r1;
si= s1;
% Integrals required for x6
F111m3=scH111m3-si*scH101m3-ri*scH011m3+ri*si*scH001m3;
F120m3=scH120m3-si*scH110m3-ri*scH020m3+ri*si*scH010m3;
F210m3=scH210m3-si*scH200m3-ri*scH110m3+ri*si*scH100m3;
F111m5=scH111m5-si*scH101m5-ri*scH011m5+ri*si*scH001m5;
F120m5=scH120m5-si*scH110m5-ri*scH020m5+ri*si*scH010m5;
F210m5=scH210m5-si*scH200m5-ri*scH110m5+ri*si*scH100m5;
F131m5=scH131m5-si*scH121m5-ri*scH031m5+ri*si*scH021m5;
F311m5=scH311m5-si*scH301m5-ri*scH211m5+ri*si*scH201m5;
F122m5=scH122m5-si*scH112m5-ri*scH022m5+ri*si*scH012m5;
F212m5=scH212m5-si*scH202m5-ri*scH112m5+ri*si*scH102m5;
F221m5=scH221m5-si*scH211m5-ri*scH121m5+ri*si*scH111m5;

fLLprime= I111m3 * F111m3 + I210m3 * F210m3 + I120m3 * F120m3 ...
    + I111m5 * F111m5 + I210m5 * F210m5 + I120m5 * F120m5 ...
    + I212m5 * F212m5 + I122m5 * F122m5 + I221m5 * F221m5 ...
    + I311m5 * F311m5 + I131m5 * F131m5;

fx6= fLLprime*factor;

%Total force on the surface
ftot=fx3+fx4+fx5+fx6;

end
