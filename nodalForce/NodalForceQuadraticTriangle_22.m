function [fx3,fx4,fx5,fx6,fx7,fx8,ftot]=NodalForceQuadraticTriangle_22(x1,x2,x3,x4,x5,b,mu,nu,a)
% COPY OF NodalForceQuadraticTriangle_vf CREATED TO MAKE SOME TESTS
%[fx3,fx4,fx5,ftot,testA,testB,testCr,testCs,TestD,TestE,TestF,testH3,testH
%5]=NodalForceQuadraticTriangle_vf(x1,x2,x3,x4,x5,b,mu,nu,a,t)
% x1,x2 Cartesian coordinates of the dislocation segment
% x3,x4,x5 Cartesian coordinates of the triangular surface element

p=double(x4-x3);
p=p/norm(p);
q=double(x5-x3);
q=q/norm(q);

t=double(x2-x1);
t=t/norm(t);

tp1=dot(p,q);
tp1=acos(tp1);
alpha=sin(tp1);


n=cross(p,q);
n=n/norm(n);
pxt=cross(p,t);
qxt=cross(q,t);
tdotn= dot(t,n);

c=dot(p,t);
d=dot(q,t);
e=dot(p,q);
c2=c^2;
d2=d^2;
e2=e^2;
a2=a^2;

y1=dot(x3-x1,n)/tdotn;
y2=dot(x3-x2,n)/tdotn;
r1=dot(x3-x1,qxt)/dot(p,qxt);
r2=dot(x4-x1,qxt)/dot(p,qxt);
s1=dot(x3-x1,pxt)/dot(q,pxt);
s2=dot(x5-x1,pxt)/dot(q,pxt);

%y22= y2*y2; y12= y1*y1;
s12= s1*s1; r12= r1*r1;

%La=y2-y1;
Lq=s2-s1;
Lp=r2-r1;
lplq=Lp/Lq;
lqlp=Lq/Lp;
s0=s2+lqlp*r1;
r0=r2+lplq*s1;
f=-lplq;
g=r0;
h=-lqlp;
m=s0;
h2=h*h;f2=f*f;
g2=g*g;m2=m*m;
eh=e*h;ef=e*f;
dh=h*d;cf=f*c;
cd=c*d;
cg=c*g;
fg=f*g;
eg=e*g;
em=e*m;
hm=h*m;
dm=d*m;
de=d*e;ce=c*e;
f3=f2*f;h3=h2*h;
m3=m2*m;g3=g2*g;
h2m=h2*m;hm2=h*m2;
f2g=f2*g;fg2=f*g2;

%rf=-lplq*sv+r0
%sf=-lqlp*rv+s0

rv= [r2;r2;r2;r2;r1;r1;r1;r1];
sv= [s2;s2;s1;s1;s2;s2;s1;s1];
yv= [y2;y1;y2;y1;y2;y1;y2;y1];
rv2=rv.*rv;
sv2=sv.*sv;
yv2=yv.*yv;
% rfs sont les valeurs que prend rf=fct(s) quand l'integration sur r est
% faite avant l'integration sur s
rfs= [-lplq*(sv(1))+r0;-lplq*(sv(2))+r0;-lplq*(sv(3))+r0;-lplq*(sv(4))+r0;r1;r1;r1;r1];
% sfr sont les valeurs que prend sf=fct(r) quand l'integration sur s est
% faite avant l'integration sur r
sfr= [-lqlp*rv(1)+s0;-lqlp*rv(2)+s0;s1;s1;-lqlp*rv(5)+s0;-lqlp*rv(6)+s0;s1;s1];
signv= [1 -1 -1 1 -1 1 1 -1];
% signr2= [1; -1; -1; 1; 0; 0; 0; 0];
% signr1= [0; 0; 0; 0; -1; 1; 1; -1];
% signs2= [1; -1; 0; 0; -1; 1; 0; 0];
% signs1= [0; 0; -1; 1; 0; 0; 1; -1];
% signy2= [1; 0; -1; 0; -1; 0; 1; 0;];
% signy1= [0; -1; 0; 1; 0; 1; 0; -1;];
ID=[1; 1; 1; 1; 1; 1; 1; 1];

%rR=yv*t+rv*p+sfr*q;
%sR=yv*t+rfs*p+sv*q;
% rRa=sqrt(yv2+rv2+sfr.^2+2*c*rv.*yv+2*d*sfr.*yv+2*e*rv.*sfr+a2);
% sRa=sqrt(yv2+rfs.^2+sv2+2*c*rfs.*yv+2*d*sv.*yv+2*e*rfs.*sv+a2);
a2v = [a2; a2; a2; a2; a2; a2; a2; a2];
g2v = [g2; g2; g2; g2; g2; g2; g2; g2];
m2v = [m2; m2; m2; m2; m2; m2; m2; m2];
clear vtemp; vtemp = [yv2 rv2 sfr.^2 2*c*rv.*yv 2*d*sfr.*yv 2*e*rv.*sfr a2v];
rRa=sqrt(smartsum(vtemp));
clear vtemp; vtemp = [yv2 rfs.^2 sv2 2*c*rfs.*yv 2*d*sv.*yv 2*e*rfs.*sv a2v];
sRa=sqrt(smartsum(vtemp));

% Rar1=sqrt(yv2+r1^2+sv2+2*c*r1*yv+2*d*sv.*yv+2*e*r1*sv+a2);
% Rar2=sqrt(yv2+(1+f2+2*e*f)*sv2+2*sv.*yv*(cf+d)+2*sv*(fg+eg)+2*cg*yv+g2+a2);


Rar1=sqrt( smartsum([ yv2 r12*ID sv2 2*c*r1*yv 2*d*sv.*yv 2*e*r1*sv a2v ]) );
Rar2=sqrt( smartsum([ yv2 (1+f2+2*e*f)*sv2 2*sv.*yv*(cf+d) 2*sv*(fg+eg) 2*cg*yv g2v a2v]) );

Ras1=sqrt( smartsum([yv2 rv2 s12*ID 2*c*rv.*yv 2*d*s1*yv 2*e*rv*s1 a2v ]) );
Ras2=sqrt( smartsum([ yv2 (1+h2+2*e*h)*rv2 2*rv.*yv*(dh+c) 2*rv*(hm+em) 2*dm*yv m2v a2v ]) );

Ra=sqrt(smartsum([ yv2 rv2 sv2 2*c*rv.*yv 2*d*sv.*yv 2*e*rv.*sv a2v ]) );

Rdotp= rv+c*yv+e*sv;
Rdotq= sv+d*yv+e*rv;
%Rdott= yv+c*rv+d*sv;

%sRdotp= rfs+c*yv+e*sv;
%rRdotq= sfr+d*yv+e*rv;
rRdott= yv+c*rv+d*sfr;
sRdott= yv+c*rfs+d*sv;

RaRdotp= Ra + rv+c*yv+e*sv;
RaRdotq= Ra + sv+d*yv+e*rv;
%RaRdott=Ra + Rdott;

%sRaRdotp=sRa + sRdotp;
%rRaRdotq=rRa + rRdotq;
rRaRdott= rRa +rRdott;
sRaRdott= sRa +sRdott;

%INTEGRALES INITIALES
% all s2A integrals are s,y functions
% all r2B integrals are r,y functions
% all C integrals must be doubled: r,y functions
% and s,y functions
i=0;j=0;k=0;l=1  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%indices=[i j k l; ii jj kk ll]
As1(:,1,3)= log(RaRdotp);
Br1(:,1,3)= log(RaRdotq);
Cr1(:,1,3)= log(sRaRdott) ;
Cs1(:,1,3)= log(rRaRdott) ;

% AA=sqrt(1+h2+2*e*h);
% BB=(yv*(h*d+c)+m*(h+e))/AA;
% CC=sqrt(-BB.^2+(yv2+2*m*d*yv+m^2+a2));
% As2(:,1,3)=log(2*sqrt((AA*rv+BB).^2+CC.^2)+2*(AA*rv+BB))/AA;

AA=(1+h2+2*e*h);
BB=(yv*(h*d+c)+m*(h+e));
CC=(yv2+2*m*d*yv+m^2+a2);
As2(:,1,3)=log(2*sqrt(AA)*sqrt(AA*rv2+2*BB.*rv+CC)+2*(AA*rv+BB))/sqrt(AA);

% AA=sqrt(1+f2+2*e*f);
% BB=(yv*(f*c+d)+g*(f+e))/AA;
% CC=sqrt(-BB.^2+(yv2+2*g*c*yv+g^2+a2));
% Br2(:,1,3)=log(2*sqrt((AA*sv+BB).^2+CC.^2)+2*(AA*sv+BB))/AA;

AA=(1+f2+2*e*f);
BB=(yv*(cf+d)+g*(f+e));
CC=(yv2+2*cg*yv+g2+a2);
Br2(:,1,3)=log(2*sqrt(AA)*sqrt(AA*sv2+2*BB.*sv+CC)+2*(AA*sv+BB))/sqrt(AA);

BB=cg+(cf+d)*sv;
CC=sqrt(-BB.^2+sv2*(1+f2+2*ef)+2*sv*(fg+eg)+g2+a2);
Cr2(:,1,3)=log(smartsum([2*sqrt((yv+BB).^2+CC.^2)  2*(yv+BB) ]));

BB=dm+(dh+c)*rv;
CC=sqrt(-BB.^2+rv2*(1+h2+2*eh)+2*rv*(hm+em)+m2+a2);
Cs2(:,1,3)=log(smartsum([2*sqrt((yv+BB).^2+CC.^2)  2*(yv+BB) ]));

% %INTEGRALES SIMPLES
% % indice i,j,k ne peuvent pas demarrer a 0 dans matlab (erreur indice tableau)
% % donc tous les indices ii,jj,kk DANS LES TABLES SEULMT doivent etre +1
% % ll doit etre i+2 car i demarre a [-1:5]
% % Ail=r^i/Ra^l
% % positive l means a negative power of Ra
% i=0;j=0;k=0;l=1;
% As1(:,ii+1,ll)=1/(l-2)*(-(rv.^i)./Ras1.^(l-2)+i*As1(:,ii-1,ll-2)-(c*yv+e*s1)*(l-2).*As1(:,ii,ll));
% As2(:,ii+1,ll)=1/(l-2)/(1+h2+2*eh)*(-(rv.^i)./Ras2.^(l-2)+i*As2(:,ii-1,ll-2)-(yv*(dh+c)+m*(h+e))*(l-2).*As2(:,ii,ll));
% Br1(:,jj+1,ll)=1/(l-2)*(-(sv.^j)./Rar1.^(l-2)+j*Br1(:,jj-1,ll-2)-(d*yv+e*r1)*(l-2).*Br1(:,jj,ll));
% Br2(:,jj+1,ll)=1/(l-2)/(1+f2+2*ef)*(-(sv.^j)./Rar2.^(l-2)+j*Br2(:,jj-1,ll-2)-(yv*(cf+d)+g*(f+e))*(l-2).*Br2(:,jj,ll));
% rC(:,kk+1,ll+2)=1/l*(k*rC(:,kk-1,ll)-yv.^k./rRa.^l)-(c*rv+d*sfr).*rC(:,kk,ll+2)
% sC(:,kk+1,ll+2)=1/l*(k*sC(:,kk-1,ll)-yv.^k./sRa.^l)-(c*rfs+d*sv).*sC(:,kk,ll+2)
%
% As1(:,ii,ll+2)= 1/l/((1-c2)*yv2+(1-e2)*s1*s1+2*(d-ec)*yv*s1+a2)*((l-i-1)*As1(:,ii,ll)-i*(c*yv+e*s1)*As1(:,ii-1,ll)+ rv.^i.*(rv+c*yv+e*s1)./Ras1.^l)
% Br1(:,jj,ll+2)= 1/l/((1-d2)*yv2+(1-e2)*r1*r1+2*(c-ed)*yv*r1+a2)*((l-j-1)*Br1(:,jj,ll)-j*(d*yv+e*r1)*Br1(:,jj-1,ll)+ sv.^j.*(sv+d*yv+e*r1)./Rar1.^l)
% As1(:,ii,ll)=1/(l-i-1)*(l*((1-c2)*yv2+(1-e2)*s1*s1+2*(d-ec)*yv*s1+a2).*As1(:,ii,ll+2)+i*(c*yv+e*s1)*As1(:,ii-1,ll)-rv.^i.*(rv+c*yv+e*s1)./Ras1.^l)

% vtmp1(:)=((dh+c)*yv+hm+em)/(1+h2+2*eh)
% vtmp2(:)=((dh+c)*yv+hm+em).^2/(1+h2+2*eh)
% As2(:,ii,ll)=1/(l-i-1)*(i*vtmp1.*As2(:,ii-1,ll)-(rv+vtmp1).*rv.^i./Ras2.^l +l*(yv2+m2+2*dm*yv+a2-vtmp2).*As2(:,ii,ll+2))

% A(:,ii,ll+2)=(1/l*1./(yv2*(1-c2)+a2+sv2+2*d*sv.*yv-2*c*e*yv.*sv-e2*sv2)).*((l-i-1)*A(:,ii,ll)-i*(c*yv+e*sv).*A(:,ii-1,ll)+(rfs.^i./sRa.^l).*(rfs+c*yv+e*sv))
% B(:,jj,ll+2)=(1/l*1./(yv2*(1-d2)+a2+yv2+2*c*rv.*yv-2*d*e*yv.*rv-e2*rv2)).*((l-j-1)*B(:,jj,ll)-j*(d*yv+e*rv).*B(:,jj-1,ll)+(sfr.^j./rRa.^l).*(sfr+d*yv+e*rv))
% rC(:,kk,ll+2)=(1/l*1./(rv2*(1-c2)+a2+sfr.^2-2*d*c*sfr.*rv-d2*sfr.^2)).*((l-k-1)*rC(:,kk,ll)-k*(d*sfr+c*rv).*rC(:,kk-1,ll)+(yv.^k./rRa^l).*(sfr+d*yv+c*rv))
% sC(:,kk,ll+2)=(1/l*1./(rfs^2*(1-c2)+a2+sv2-2*c*d*sv.*rfs-d2*sv2)).*((l-k-1)*sC(:,kk,ll)-k*(d*sv+c*rfs).*sC(:,kk-1,ll)+(yv.^k./sRa^l).*(sv+d*yv+c*rfs))

% ABC[1 m1]
i=0;j=0;k=0;l=-1  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%indices=[i j k l; ii jj kk ll]
As1(:,2,3)=rRa-(e*s1+c*yv).*As1(:,1,3);
As2(:,2,3)=(rRa-(hm+em+yv*(dh+c)).*As2(:,1,3))*1/(1+h2+2*eh);

Br1(:,2,3)=sRa-(e*r1+d*yv).*Br1(:,1,3);
Br2(:,2,3)=(sRa-(fg+eg+yv*(cf+d)).*Br2(:,1,3))*1/(1+f2+2*ef);
 
Cr1(:,2,3)=sRa-(c*r1+d*sv).*Cr1(:,1,3);
Cr2(:,2,3)=sRa-(c*rfs+d*sv).*Cr2(:,1,3);

Cs1(:,2,3)=rRa-(d*s1+c*rv).*Cs1(:,1,3);
Cs2(:,2,3)=rRa-(d*sfr+c*rv).*Cs2(:,1,3);

% BB=2*(c*r1+d*sv);
% CC=2*e*r1*sv+r1^2+sv.^2+a2;
% Cr1(:,2,3)= sRa-(c*r1+d*sv).*log(2*(sRa+yv+c*r1+d*sv));
% 
% BB=2*(c*rfs+d*sv);
% CC=2*e*rfs.*sv+rfs.^2+sv.^2+a2;
% Cr2(:,2,3)= sRa-(c*rfs+d*sv).*log(2*(sRa+yv+c*rfs+d*sv));

%pause


%ABC[0 1] ABC[1 1]
i=0;j=0;k=0;l=-3  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%indices=[i j k l; ii jj kk ll]'
As1(:,ii,ll+2)=0.5*(Rdotp.*Ra + (a2+(1-e2)*s12+2*(d-ce)*s1*yv+(1-c2)*yv2).*As1(:,1,3));
As2(:,ii,ll+2)= 0.5*(rv.*rRa+(yv*(dh+c)+m*(h+e)).*As2(:,2,3)+(yv2+2*dm*yv+m2+a2).*As2(:,1,3));
As1(:,ii+1,ll+2)=1/l*(-1.*Ras1.^(-l))-(e*sv+c*yv).*As1(:,ii,ll+2);
As2(:,ii+1,ll+2)=1/(1+h2+2*eh)*(1/3*rRa.^3-(yv*(dh+c)+m*(h+e)).*As2(:,ii,ll+2));

Br1(:,jj,ll+2)=0.5*(Rdotq.*Ra + (a2+(1-e2)*r12+2*(c-de)*r1*yv+(1-d2)*yv2).*Br1(:,1,3));
Br2(:,jj,ll+2)= 0.5*(sv.*sRa+(yv*(cf+d)+g*(f+e)).*Br2(:,2,3)+(yv2+2*cg*yv+g2+a2).*Br2(:,1,3));
Br1(:,jj+1,ll+2)=1/l*(-1.*Rar1.^(-l))-(e*rv+d*yv).*Br1(:,jj,ll+2);
Br2(:,jj+1,ll+2)=1/(1+f2+2*ef)*(1/3*sRa.^3-(yv*(cf+d)+g*(f+e)).*Br2(:,jj,ll+2));

Cr1(:,kk,ll+2)=0.5*(rRdott.*rRa + (rRa.^2-(rRdott.^2)).*Cr1(:,1,3));
Cr2(:,kk,ll+2)=0.5*(+Rar2.*(yv+sv*(cf+d)+cg)+((1+f2+2*ef-(cf+d)^2)*sv2+2*sv*(fg+eg-cg*(cf+d))+(1-c2)*g2+a2).*Cr2(:,1,3));
Cr1(:,kk+1,ll+2)=1/l*(-1./sRa.^l)-(c*rfs+d*sv).*Cr1(:,kk,ll+2);
Cr2(:,kk+1,ll+2)=1/3*Rar2.^3-(sv*(cf+d)+cg).*Cr2(:,kk,ll+2);

Cs1(:,kk,ll+2)=0.5*(sRdott.*sRa + (sRa.^2-(sRdott.^2)).*Cs1(:,1,3));
Cs2(:,kk,ll+2)=0.5*(+Ras2.*(yv+rv*(dh+c)+dm)+((1+h2+2*eh-(dh+c)^2)*rv2+2*rv*(hm+em-dm*(dh+c))+(1-d2)*m2+a2).*Cs2(:,1,3));
Cs1(:,kk+1,ll+2)=1/l*(-1./rRa.^l)-(c*rv+d*sfr).*Cs1(:,kk,ll+2);
Cs2(:,kk+1,ll+2)=1/3*Ras2.^3-(rv*(dh+c)+dm).*Cs2(:,kk,ll+2);

% ABC[0 3]
i=0; j=0; k=0; l=-3  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%A03s1(1:8) =
%1/(l-i-1)*(l*((1-c2)*yv2+(1-e2)*s1*s1+2*(d-ce)*yv*s1+a2).*As1(:,ii,ll+2)-(rv+c*yv+e*s1)./Ras1.^l);
A03s1(1:8)=1/4*(Rdotp.*Ras1.^3 + 3*(a2+(1-e2)*s12+2*(d-ce)*s1*yv+(1-c2)*yv2).*As1(:,1,1));
vtmp1(1:8)=(hm+em+(dh+c)*yv)/(1+h2+2*eh);vtmp1=vtmp1';
vtmp2(1:8)=(hm+em+(dh+c)*yv).^2/(1+h2+2*eh); vtmp2= vtmp2';
A03s2(1:8)=1/(l-i-1)*(smartsum([ -(rv+vtmp1).*Ras2.^-l   +l*a2.*As2(1:8,ii,ll+2) +l*m2.*As2(1:8,ii,ll+2) +l*2*dm*yv.*As2(1:8,ii,ll+2) +l*(-vtmp2).*As2(1:8,ii,ll+2) +l*(yv2).*As2(1:8,ii,ll+2) ]));

B03r1(1:8)=1/4*(Rdotq.*Rar1.^3 + 3*(a2+(1-e2)*r12+2*(c-de)*r1*yv+(1-d2)*yv2).*Br1(:,1,1));
vtmp1(1:8)=(fg+eg+(cf+d)*yv)/(1+f2+2*ef);
vtmp2(1:8)=((cf+d)*yv+fg+eg).^2/(1+f2+2*ef);
B03r2(1:8)=1/(l-j-1)*(smartsum([ -(sv+vtmp1).*Rar2.^-l   +l*(+a2).*Br2(1:8,1,1) +l*(g2).*Br2(1:8,1,1) +l*(2*cg*yv).*Br2(1:8,1,1) +l*(-vtmp2).*Br2(1:8,1,1) +l*(yv2).*Br2(1:8,1,1)]));

C03s1(1:8)=1/4*(sRdott.*Ras1.^3 + 3*(a2*ID+(1-c2)*rv2+(1-d2)*s12*ID+2*rv*s1*(e-cd)).*Cs1(:,1,1));
C03s2(1:8)=1/4*(+Ras2.^3.*(yv+rv*(dh+c)+dm)+3*((1+h2+2*eh-(dh+c)^2)*rv2+2*rv*(hm+em-dm*(dh+c))+(1-d2)*m2+a2).*Cs2(:,1,1));
C03r1(1:8)=1/4*(rRdott.*Rar1.^3 + 3*(a2*ID+(1-d2)*sv2+(1-c2)*r12*ID+2*sv*r1*(e-cd)).*Cr1(:,1,1));
C03r2(1:8)=1/4*(+Rar2.^3.*(yv+sv*(cf+d)+cg)+3*((1+f2+2*ef-(cf+d)^2)*sv2+2*sv*(fg+eg-cg*(cf+d))+(1-c2)*g2+a2).*Cr2(:,1,1));

% bb= (2*c*rv +2*d*s1);
% cc= (rv.^2 +s1^2 +2*e*rv*s1 +a2);
% C03s1(1:8) = 1/64*(bb+2*yv).*(-3*bb.^2 +8*bb.*yv +20*cc +8*yv2).*sqrt(yv2+bb.*yv+cc)+3/128*(bb.^2 -4*cc).^2.*log(2*(sqrt(yv2+bb.*yv+cc)+yv)+bb);
% bb= (2*c*rv +2*d*sfr);
% cc= (rv.^2 +sfr.^2 +2*e*rv.*sfr +a2);
% C03s2(1:8) = 1/64*(bb+2*yv).*(-3*bb.^2 +8*bb.*yv +20*cc +8*yv2).*sqrt(yv2+bb.*yv+cc)+3/128*(bb.^2 -4*cc).^2.*log(2*(sqrt(yv2+bb.*yv+cc)+yv)+bb);


% ABC[1 3]
i=0; j=0; k=0; l=-3  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;

A13s1(1:8)=1/(l-2)*(-1./Ras1.^(l-2))-(c*yv+e*s1).*A03s1(:);
A13s2(1:8)=1/(l-2)/(1+h2+2*eh)*(-1./Ras2.^(l-2)-(yv*(dh+c)+m*(h+e))*(l-2).*A03s2(:));
B13r1(1:8)=1/(l-2)*(-1./Rar1.^(l-2))-(d*yv+e*r1).*B03r1(:);
B13r2(1:8)=1/(l-2)/(1+f2+2*ef)*(-1./Rar2.^(l-2)-(yv*(cf+d)+g*(f+e))*(l-2).*B03r2(:));

%ABC[2 1]
i=1;j=1;k=0;l=-1  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%indices=[i j k l; ii jj kk ll]
As1(:,ii+1,ll)=1/(l-2)*(-rv./Ras1.^(l-2)+A03s1(:)-(c*yv+e*s1)*(l-2).*As1(:,ii,ll));
As2(:,ii+1,ll)=1/(l-2)/(1+h2+2*eh)*(-rv./Ras2.^(l-2)+i*A03s2(:)-(yv*(dh+c)+m*(h+e))*(l-2).*As2(:,ii,ll));
Br1(:,jj+1,ll)=1/(l-2)*(-sv./Rar1.^(l-2)+B03r1(:)-(d*yv+e*r1)*(l-2).*Br1(:,jj,ll));
Br2(:,jj+1,ll)=1/(l-2)/(1+f2+2*ef)*(-sv./Rar2.^(l-2)+j*B03r2(:)-(yv*(cf+d)+g*(f+e))*(l-2).*Br2(:,jj,ll));

%ABC[2 m1]
i=1;j=1;k=0;l=1  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
As1(:,ii+1,ll)=1/(l-2)*(-rv./Ras1.^(l-2)+i*As1(:,ii-1,ll-2)-(c*yv+e*s1)*(l-2).*As1(:,ii,ll));
As2(:,ii+1,ll)=1/(l-2)/(1+h2+2*eh)*(-rv./Ras2.^(l-2)+i*As2(:,ii-1,ll-2)-(yv*(dh+c)+m*(h+e))*(l-2).*As2(:,ii,ll));
Br1(:,jj+1,ll)=1/(l-2)*(-sv./Rar1.^(l-2)+j*Br1(:,jj-1,ll-2)-(d*yv+e*r1)*(l-2).*Br1(:,jj,ll));
Br2(:,jj+1,ll)=1/(l-2)/(1+f2+2*ef)*(-sv./Rar2.^(l-2)+j*Br2(:,jj-1,ll-2)-(yv*(cf+d)+g*(f+e))*(l-2).*Br2(:,jj,ll));

%ABC[3 m1]
i=2;j=2;k=0;l=1  ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
As1(:,ii+1,ll)=1/(l-2)*(-(rv.^i)./Ras1.^(l-2)+i*As1(:,ii-1,ll-2)-(c*yv+e*s1)*(l-2).*As1(:,ii,ll));
As2(:,ii+1,ll)=1/(l-2)/(1+h2+2*eh)*(-(rv.^i)./Ras2.^(l-2)+i*As2(:,ii-1,ll-2)-(yv*(dh+c)+m*(h+e))*(l-2).*As2(:,ii,ll));
Br1(:,jj+1,ll)=1/(l-2)*(-(sv.^j)./Rar1.^(l-2)+j*Br1(:,jj-1,ll-2)-(d*yv+e*r1)*(l-2).*Br1(:,jj,ll));
Br2(:,jj+1,ll)=1/(l-2)/(1+f2+2*ef)*(-(sv.^j)./Rar2.^(l-2)+j*Br2(:,jj-1,ll-2)-(yv*(cf+d)+g*(f+e))*(l-2).*Br2(:,jj,ll));

As1(1:2,:,:)=0;As1(5:6,:,:)=0;
As2(3:4,:,:)=0;As2(7:8,:,:)=0;
A03s1(1:2)=0;A03s1(5:6)=0;
A03s2(3:4)=0;A03s2(7:8)=0;
A13s1(1:2)=0;A13s1(5:6)=0;
A13s2(3:4)=0;A13s2(7:8)=0;
B03r1(1:4)=0;B03r2(5:8)=0;
B13r1(1:4)=0;B13r2(5:8)=0;
Br1(1:4,:,:)=0;Br2(5:8,:,:)=0;
Cr1(1:4,:,:)=0;Cr2(5:8,:,:)=0;
Cs1(1:2,:,:)=0;Cs1(5:6,:,:)=0;
Cs2(3:4,:,:)=0;Cs2(7:8,:,:)=0;
C03s1(1:2)=0;C03s1(5:6)=0;
C03s2(3:4)=0;C03s2(7:8)=0;
C03r1(1:4)=0;C03r2(5:8)=0;

% % DEBUGGUNG A SINGLE INTEGRALS
% F22=@(r) 1.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) 1.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) 1.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) 1.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIA0m1=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) r.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) r.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) r.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIA1m1=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) 1.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(1/2);
% F21=@(r) 1.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(1/2);
% F12=@(r) 1.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(1/2);
% F11=@(r) 1.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(1/2);
% NIA01=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(1/2);
% F21=@(r) r.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(1/2);
% F12=@(r) r.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(1/2);
% F11=@(r) r.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(1/2);
% NIA11=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.*r.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) r.*r.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) r.*r.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) r.*r.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIA2m1=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) (r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(3/2);
% F21=@(r) (r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(3/2);
% F12=@(r) (r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(3/2);
% F11=@(r) (r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(3/2);
% NIA03=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.*r.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(1/2);
% F21=@(r) r.*r.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(1/2);
% F12=@(r) r.*r.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(1/2);
% F11=@(r) r.*r.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(1/2);
% NIA21=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.^3.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) r.^3.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) r.^3.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) r.^3.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIA3m1=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% F22=@(r) r.*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(3/2);
% F21=@(r) r.*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(3/2);
% F12=@(r) r.*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(3/2);
% F11=@(r) r.*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(3/2);
% NIA13=quad(F22,r1,r2,1e-16)-quad(F21,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16);
% 
% vA0m1(1:2)=As2(1:2,1,3);vA0m1(3:4)=As1(3:4,1,3);vA0m1(5:6)=As2(5:6,1,3);vA0m1(7:8)=As1(7:8,1,3);
% vA1m1(1:2)=As2(1:2,2,3);vA1m1(3:4)=As1(3:4,2,3);vA1m1(5:6)=As2(5:6,2,3);vA1m1(7:8)=As1(7:8,2,3);
% vA01(1:8)=As2(1:8,1,1)+As1(1:8,1,1);
% vA11(1:8)=As2(1:8,2,1)+As1(1:8,2,1);
% vA2m1(1:8)=As2(1:8,3,3)+As1(1:8,3,3);
% vA3m1(1:8)=As2(1:8,4,3)+As1(1:8,4,3);
% vA03(1:8)=A03s2(1:8)+A03s1(1:8);
% vA21(1:8)=As2(1:8,3,1)+As1(1:8,3,1);
% vA13(1:8)=A13s2(1:8)+A13s1(1:8);
% testA(:,1)=[SmartDot(vA0m1,signv) SmartDot(vA1m1,signv) SmartDot(vA01,signv) SmartDot(vA11,signv) SmartDot(vA2m1,signv) SmartDot(vA3m1,signv) SmartDot(vA03,signv) SmartDot(vA21,signv) SmartDot(vA13,signv) ];
% 
% testA(:,2)=[NIA0m1 NIA1m1  NIA01 NIA11 NIA2m1 NIA3m1 NIA03 NIA21 NIA13];
% testA(:,3)=(testA(:,1)-testA(:,2));
% nt=size(testA(:,1)); TEST_A(3)= 0;
% for i=1:nt
%     TEST_A(:)= testA(i,:)
% end

%%pause

% % DEBUGGUNG B SINGLE INTEGRALS
% F22=@(s) 1.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-1/2);
% F21=@(s) 1.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-1/2);
% F12=@(s) 1.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-1/2);
% F11=@(s) 1.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-1/2);
% NIB0m1=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-1/2);
% F21=@(s) s.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-1/2);
% F12=@(s) s.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-1/2);
% F11=@(s) s.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-1/2);
% NIB1m1=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.*s.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-1/2);
% F21=@(s) s.*s.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-1/2);
% F12=@(s) s.*s.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-1/2);
% F11=@(s) s.*s.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-1/2);
% NIB2m1=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.^3.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-1/2);
% F21=@(s) s.^3.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-1/2);
% F12=@(s) s.^3.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-1/2);
% F11=@(s) s.^3.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-1/2);
% NIB3m1=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) 1.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(1/2);
% F21=@(s) 1.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(1/2);
% F12=@(s) 1.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(1/2);
% F11=@(s) 1.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(1/2);
% NIB01=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) 1.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(3/2);
% F21=@(s) 1.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(3/2);
% F12=@(s) 1.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(3/2);
% F11=@(s) 1.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(3/2);
% NIB03=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(1/2);
% F21=@(s) s.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(1/2);
% F12=@(s) s.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(1/2);
% F11=@(s) s.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(1/2);
% NIB11=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.*s.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(1/2);
% F21=@(s) s.*s.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(1/2);
% F12=@(s) s.*s.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(1/2);
% F11=@(s) s.*s.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(1/2);
% NIB21=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% F22=@(s) s.*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(3/2);
% F21=@(s) s.*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(3/2);
% F12=@(s) s.*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(3/2);
% F11=@(s) s.*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(3/2);
% NIB13=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% vB0m1(1:8)=Br2(1:8,1,3)+Br1(1:8,1,3);
% vB1m1(1:8)=Br2(1:8,2,3)+Br1(1:8,2,3);
% vB01(1:8)=Br2(1:8,1,1)+Br1(1:8,1,1);
% vB03(1:8)=B03r2(1:8)+B03r1(1:8);
% vB13(1:8)=B13r2(1:8)+B13r1(1:8);
% vB11(1:8)=Br2(1:8,2,1)+Br1(1:8,2,1);
% vB21(1:8)=Br2(1:8,3,1)+Br1(1:8,3,1);
% vB2m1(1:8)=Br2(1:8,3,3)+Br1(1:8,3,3);
% vB3m1(1:8)=Br2(1:8,4,3)+Br1(1:8,4,3);
% testB(:,1)=[SmartDot(vB0m1,signv) SmartDot(vB1m1,signv) SmartDot(vB01,signv) SmartDot(vB11,signv) SmartDot(vB2m1,signv) SmartDot(vB3m1,signv) SmartDot(vB03,signv) SmartDot(vB21,signv) SmartDot(vB13,signv) ];
% testB(:,2)=[NIB0m1 NIB1m1 NIB01  NIB11 NIB2m1 NIB3m1 NIB03 NIB21  NIB13];
% testB(:,3)=(testB(:,1)-testB(:,2));
% nt=size(testB(:,1)); TEST_B(3)= 0;
% for i=1:nt
%     TEST_B(:)= testB(i,:)
% end
%pause

% DEBUGGUNG C SINGLE INTEGRALS
% vrC0m1(1:8)=Cr2(1:8,1,3)+Cr1(1:8,1,3);vsC0m1(1:8)=Cs2(1:8,1,3)+Cs1(1:8,1,3);
% vrC1m1(1:8)=Cr2(1:8,2,3)+Cr1(1:8,2,3);vsC1m1(1:8)=Cs2(1:8,2,3)+Cs1(1:8,2,3);
% vrC01(1:8)=Cr2(1:8,1,1)+Cr1(1:8,1,1);vsC01(1:8)=Cs2(1:8,1,1)+Cs1(1:8,1,1);
% vrC11(1:8)=Cr2(1:8,2,1)+Cr1(1:8,2,1);vsC11(1:8)=Cs2(1:8,2,1)+Cs1(1:8,2,1);
% vsC03(1:8)=C03s1(1:8)+C03s2(1:8); vrC03(1:8)=C03r1(1:8)+C03r2(1:8);

% F22=@(y) 1./(y.^2+(-lplq*s2+r0)^2+s2^2+2*c*(-lplq*s2+r0)*y+2*d*s2*y+2*e*(-lplq*s2+r0)*s2+a2).^(1/2);
% F21=@(y) 1./(y.^2+(-lplq*s1+r0)^2+s1^2+2*c*(-lplq*s1+r0)*y+2*d*s1*y+2*e*(-lplq*s1+r0)*s1+a2).^(1/2);
% F12=@(y) 1./(y.^2+r1^2+s2^2+2*c*(r1)*y+2*d*s2*y+2*e*(r1)*s2+a2).^(1/2);
% F11=@(y) 1./(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIrC0m1=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) y./(y.^2+(-lplq*s2+r0)^2+s2^2+2*c*(-lplq*s2+r0)*y+2*d*s2*y+2*e*(-lplq*s2+r0)*s2+a2).^(1/2);
% F21=@(y) y./(y.^2+(-lplq*s1+r0)^2+s1^2+2*c*(-lplq*s1+r0)*y+2*d*s1*y+2*e*(-lplq*s1+r0)*s1+a2).^(1/2);
% F12=@(y) y./(y.^2+r1^2+s2^2+2*c*(r1)*y+2*d*s2*y+2*e*(r1)*s2+a2).^(1/2);
% F11=@(y) y./(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIrC1m1=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) (y.^2+(-lplq*s2+r0)^2+s2^2+2*c*(-lplq*s2+r0)*y+2*d*s2*y+2*e*(-lplq*s2+r0)*s2+a2).^(1/2);
% F21=@(y) (y.^2+(-lplq*s1+r0)^2+s1^2+2*c*(-lplq*s1+r0)*y+2*d*s1*y+2*e*(-lplq*s1+r0)*s1+a2).^(1/2);
% F12=@(y) (y.^2+r1^2+s2^2+2*c*(r1)*y+2*d*s2*y+2*e*(r1)*s2+a2).^(1/2);
% F11=@(y) (y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIrC01=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) y.*(y.^2+(-lplq*s2+r0)^2+s2^2+2*c*(-lplq*s2+r0)*y+2*d*s2*y+2*e*(-lplq*s2+r0)*s2+a2).^(1/2);
% F21=@(y) y.*(y.^2+(-lplq*s1+r0)^2+s1^2+2*c*(-lplq*s1+r0)*y+2*d*s1*y+2*e*(-lplq*s1+r0)*s1+a2).^(1/2);
% F12=@(y) y.*(y.^2+r1^2+s2^2+2*c*(r1)*y+2*d*s2*y+2*e*(r1)*s2+a2).^(1/2);
% F11=@(y) y.*(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIrC11=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) (y.^2+(-lplq*s2+r0)^2+s2^2+2*c*(-lplq*s2+r0)*y+2*d*s2*y+2*e*(-lplq*s2+r0)*s2+a2).^(3/2);
% F21=@(y) (y.^2+(-lplq*s1+r0)^2+s1^2+2*c*(-lplq*s1+r0)*y+2*d*s1*y+2*e*(-lplq*s1+r0)*s1+a2).^(3/2);
% F12=@(y) (y.^2+r1^2+s2^2+2*c*(r1)*y+2*d*s2*y+2*e*(r1)*s2+a2).^(3/2);
% F11=@(y) (y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(3/2);
% NIrC03_1=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% 
% clear bb cc test
% bb= (2*c*rfs +2*d*sv);
% cc= (rfs.^2 +sv.^2 +2*e*rfs.*sv +a2);
% test = 1/64*(bb+2*yv).*(-3*bb.^2 +8*bb.*yv +20*cc +8*yv2).*sqrt(yv2+bb.*yv+cc)+3/128*(bb.^2 -4*cc).^2.*log(2*(sqrt(yv2+bb.*yv+cc)+yv)+bb);
% NIrC03_2 = SmartDot(test(:),signv(:))
% 
% testCr(:,1)=[vrC0m1*signv' vrC1m1*signv'  vrC01*signv' vrC11*signv' vrC03*signv' vrC03*signv'];
% testCr(:,2)=[NIrC0m1 NIrC1m1  NIrC01 NIrC11 NIrC03_1 NIrC03_2];
% testCr(:,3)=(testCr(:,1)-testCr(:,2));
% nt=size(testCr(:,1)); TEST_Cr(3)= 0;
% for i=1:nt
%     TEST_Cr(:)= testCr(i,:)
% end

% F22=@(y) 1./(y.^2+r2^2+(h*r2+m)^2+2*c*r2*y+2*d*(h*r2+m)*y+2*e*r2*(h*r2+m)+a2).^(1/2);
% F21=@(y) 1./(y.^2+r2^2+s1^2+2*c*r2*y+2*d*s1*y+2*e*r2*s1+a2).^(1/2);
% F12=@(y) 1./(y.^2+r1^2+(h*r1+m)^2+2*c*(r1)*y+2*d*(h*r1+m)*y+2*e*(r1)*(h*r1+m)+a2).^(1/2);
% F11=@(y) 1./(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIsC0m1=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) y./(y.^2+r2^2+(h*r2+m)^2+2*c*r2*y+2*d*(h*r2+m)*y+2*e*r2*(h*r2+m)+a2).^(1/2);
% F21=@(y) y./(y.^2+r2^2+s1^2+2*c*r2*y+2*d*s1*y+2*e*r2*s1+a2).^(1/2);
% F12=@(y) y./(y.^2+r1^2+(h*r1+m)^2+2*c*(r1)*y+2*d*(h*r1+m)*y+2*e*(r1)*(h*r1+m)+a2).^(1/2);
% F11=@(y) y./(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIsC1m1=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) 1.*(y.^2+r2^2+(h*r2+m)^2+2*c*r2*y+2*d*(h*r2+m)*y+2*e*r2*(h*r2+m)+a2).^(1/2);
% F21=@(y) 1.*(y.^2+r2^2+s1^2+2*c*r2*y+2*d*s1*y+2*e*r2*s1+a2).^(1/2);
% F12=@(y) 1.*(y.^2+r1^2+(h*r1+m)^2+2*c*(r1)*y+2*d*(h*r1+m)*y+2*e*(r1)*(h*r1+m)+a2).^(1/2);
% F11=@(y) 1.*(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIsC01=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) y.*(y.^2+r2^2+(h*r2+m)^2+2*c*r2*y+2*d*(h*r2+m)*y+2*e*r2*(h*r2+m)+a2).^(1/2);
% F21=@(y) y.*(y.^2+r2^2+s1^2+2*c*r2*y+2*d*s1*y+2*e*r2*s1+a2).^(1/2);
% F12=@(y) y.*(y.^2+r1^2+(h*r1+m)^2+2*c*(r1)*y+2*d*(h*r1+m)*y+2*e*(r1)*(h*r1+m)+a2).^(1/2);
% F11=@(y) y.*(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(1/2);
% NIsC11=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% F22=@(y) 1.*(y.^2+r2^2+(h*r2+m)^2+2*c*r2*y+2*d*(h*r2+m)*y+2*e*r2*(h*r2+m)+a2).^(3/2);
% F21=@(y) 1.*(y.^2+r2^2+s1^2+2*c*r2*y+2*d*s1*y+2*e*r2*s1+a2).^(3/2);
% F12=@(y) 1.*(y.^2+r1^2+(h*r1+m)^2+2*c*(r1)*y+2*d*(h*r1+m)*y+2*e*(r1)*(h*r1+m)+a2).^(3/2);
% F11=@(y) 1.*(y.^2+r1^2+s1^2+2*c*(r1)*y+2*d*s1*y+2*e*(r1)*s1+a2).^(3/2);
% NIsC03=quad(F22,y1,y2,1e-16)-quad(F12,y1,y2,1e-16)-quad(F21,y1,y2,1e-16)+quad(F11,y1,y2,1e-16);
% 
% testCs(:,1)=[vsC0m1*signv' vsC1m1*signv'  vsC01*signv' vsC11*signv' vsC03*signv'];
% testCs(:,2)=[NIsC0m1 NIsC1m1  NIsC01 NIsC11 NIsC03];
% testCs(:,3)=(testCs(:,1)-testCs(:,2));
% nt=size(testCs(:,1)); TEST_Cs(3)= 0;
% for i=1:nt
%     TEST_Cs(:)= testCs(i,:)
% end

%pause

% temp(1:8)=Asi(1:8,ii,jj,ll-2);
% testAsi=temp(1:8)*signv'
% F22=@(r) r.^(ii-1).*(h*r+m).^(jj-1).*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) r.^(ii-1).*(h*r+m).^(jj-1).*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) r.^(ii-1).*s1^(jj-1).*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) r.^(ii-1).*s1^(jj-1).*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIAsi=quad(F22,r1,r2)-quad(F21,r1,r2)-quad(F12,r1,r2)+quad(F11,r1,r2)
% temp(1:8)=Bri(1:8,ii,jj,ll-2);
% testBri=temp(1:8)*signv'
% F22=@(s) (f*s+g).^(ii-1).*s.^(jj-1).*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-(ll-2)/2);
% F21=@(s) (f*s+g).^(ii-1).*s.^(jj-1).*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-(ll-2)/2);
% F12=@(s) r1.^(ii-1).*s.^(jj-1).*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-(ll-2)/2);
% F11=@(s) r1.^(ii-1).*s.^(jj-1).*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-(ll-2)/2);
% NIBri=quad(F22,s1,s2,1e-16)-quad(F21,s1,s2,1e-16)-quad(F12,s1,s2,1e-16)+quad(F11,s1,s2,1e-16);
% %pause

A=As1+As2;
B=Br1+Br2;
rC=Cr2+Cr1;
sC=Cs2+Cs1;

%INTEGRALES DOUBLES
%Integrales E sont fonctions seulement de s,y
%Integrales F sont fonctions seulement de r,y
%rien de change pour le moment pour les integrales D
%SEED DOUBLE E[0,0,3]
%!!!!!!!!!!
% tous mes prefacteurs en r ou s devant une integral simple doivent rentrer
% sous l'integrals pour r2F et s2E, D et H
%!!!!!!!!!!
Fr1(1:8,1:7,1:4,1:5)=0;Fr2(1:8,1:7,1:4,1:5)=0;
Es1(1:8,1:7,1:4,1:5)=0;Es2(1:8,1:7,1:4,1:5)=0;
i=0 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%indices=[i j k l; ii jj kk ll]
root=1./sqrt((1-c2)*a2+((1-c2)*(1-e2)-(d-c*e)^2)*sv2);
Es1(:,ii,kk,ll)=2.*root.*atan(((1-c)*(sRa-rv-c*yv-e*sv)+(1-c2)*yv+(d-c*e)*sv).*root);
root=1./sqrt((1-d2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*rv2);
Fr1(:,jj,kk,ll)=2.*root.*atan(((1-d)*(rRa-sv-d*yv-e*rv)+(1-d2)*yv+(c-d*e)*rv).*root);

%root=1./sqrt((1-e2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*yv2);
%Dr1(:,jj,kk,ll)=2*root.*atan(((1-e)*(Rar1-sv-d*yv-e*rv)+(c-d*e)*yv+(1-e2)*rv).*root);
root=1./sqrt((1-e2)*a2+(1-c2-d2-e2+2*c*de)*yv2);
Dr1(:,jj,kk,ll)=2*root.*atan(((1-e)*(Rar1+rv-sv)+(c-d)*yv).*root); %more compact form

% tp1=sqrt(1+f2+2*ef);
% tp2=(yv*(cf+d)+fg+eg)/tp1;
% tp3=sqrt(tp1-f-e);
% tp4=(tp1*(c*yv+g)-tp2*(f+e))/tp3;
% tp5=sqrt(-tp4.^2-(tp1+f+e)*tp2.^2+(tp1+f+e)*(a2+g2+yv2+2*cg*yv));
% Dr2(:,jj,kk,ll)=2/tp3./tp5.*atan((tp3*(Rar2-tp1*sv-tp2)+tp4)./tp5);
AA=sqrt(1+f2+2*ef);
DD=sqrt(AA-f-e);
root=1./sqrt(yv2*(1-c2-d2-e2+2*c*de)+a2*(1-e2));
Dr2(:,jj,kk,ll)=2*root.*atan((DD^2*Rar2+(AA*c-cf-d)*yv+AA*g-fg-eg-AA*DD^2*sv).*root);
%more compact form

root=1./sqrt((1-e2)*a2+(1-c2-d2-e2+2*c*de)*yv2);
Ds1(:,jj,kk,ll)=2*root.*atan(((1-e)*(Ras1+sv-rv)+(d-c)*yv).*root); %more compact form
AA=sqrt(1+h2+2*eh);
DD=sqrt(AA-h-e);
root=1./sqrt(yv2*(1-c2-d2-e2+2*c*de)+a2*(1-e2));
Ds2(:,jj,kk,ll)=2*root.*atan((DD^2*Ras2+(AA*d-dh-c)*yv+AA*m-hm-em-AA*DD^2*rv).*root);
%more compact form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dr1(1:4,:,:)=0;Dr2(5:8,:,:)=0;
D00m3_1(1:8)=Dr1(1:8,1,1,5)+Dr2(1:8,1,1,5);
Ds1(1:2,:,:)=0;Ds1(5:6,:,:)=0;
Ds2(3:4,:,:)=0;Ds2(7:8,:,:)=0;
D00m3_2(1:8)=Ds1(1:8,1,1,5)+Ds2(1:8,1,1,5);
D(1:8,1,1,5)=0.5*D00m3_1(1:8) +0.5*D00m3_2(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp1=sqrt(1+h2+2*eh);
tp2=(yv*(dh+c)+hm+em)/tp1;
tp3=sqrt(tp1-dh-c);
tp4=(tp1*(yv+dm)-tp2*(dh+c))/tp3;
tp5=sqrt(-tp4.^2-(tp1+dh+c)*tp2.^2+(tp1+dh+c)*(a2+m2+yv2+2*dm*yv));
Es2(:,ii,kk,ll)=2/tp3./tp5.*atan((tp3*(Ras2-tp1*rv-tp2)+tp4)./tp5);

% AA=sqrt(1+h2+2*eh);
% DD=sqrt(AA-dh-c);
% root=1./sqrt(m2*(1-c2-d2-e2+2*c*de)+a2*(AA^2-(dh+c)^2));
% Es2(:,ii,kk,ll)=2*root.*atan((DD^2*Ras2+(AA-dh-c)*yv+AA*dm-hm-em-AA*DD^2*rv).*root);

tp1=sqrt(1+f2+2*ef);
tp2=(yv*(cf+d)+fg+eg)/tp1;
tp3=sqrt(tp1-cf-d);
tp4=(tp1*(yv+cg)-tp2*(cf+d))/tp3;
tp5=sqrt(-tp4.^2-(tp1+cf+d)*tp2.^2+(tp1+cf+d)*(a2+g2+yv2+2*cg*yv));
Fr2(:,jj,kk,ll)=2/tp3./tp5.*atan((tp3*(Rar2-tp1*sv-tp2)+tp4)./tp5);

% AA=sqrt(1+f2+2*ef);
% DD=sqrt(AA-cf-d);
% root=1./sqrt(g2*(1-c2-d2-e2+2*c*de)+a2*(AA^2-(cf+d)^2));
% Fr2(:,jj,kk,ll)=2*root.*atan((DD^2*Rar2+(AA-cf-d)*yv+AA*cg-fg-eg-AA*DD^2*
% sv).*root);

%DEF[1,0,3] DEF[0,1,3]
i=0 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll+2)=1/(1-c2)*(1/l*(c*As1(:,ii,ll)-Cs1(:,kk,ll))+(c*d-e)*sv.*Es1(:,ii,kk,ll+2));
Es1(:,ii,kk+1,ll+2)=1/(1-c2)*(1/l*(c*Cs1(:,kk,ll)-As1(:,ii,ll))+(c*e-d)*sv.*Es1(:,ii,kk,ll+2));
Es2(:,ii+1,kk,ll+2)=1/l/(1+h2+2*eh-(dh+c)^2)*((dh+c)*As2(:,ii,ll)-Cs2(:,kk,ll)-l*(hm+em-dm*(dh+c))*Es2(:,jj,kk,ll+2) );
Es2(:,ii,kk+1,ll+2)=1/l/(1+h2+2*eh-(dh+c)^2)*((dh+c)*Cs2(:,ii,ll)-(1+h2+2*eh)*As2(:,ii,ll)+l*((dh+c)*(hm+em)-(1+h2+2*eh)*dm)*Es2(:,jj,kk,ll+2));
Fr1(:,jj+1,kk,ll+2)=1/(1-d2)*(1/l*(d*Br1(:,jj,ll)-Cr1(:,kk,ll))+(c*d-e)*rv.*Fr1(:,jj,kk,ll+2));
Fr2(:,jj+1,kk,ll+2)=1/l/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Br2(:,jj,ll)-Cr2(:,kk,ll)-l*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll+2) );
Fr1(:,jj,kk+1,ll+2)=1/(1-d2)*(1/l*(d*Cr1(:,kk,ll)-Br1(:,jj,ll))+(d*e-c)*rv.*Fr1(:,jj,kk,ll+2));
Fr2(:,jj,kk+1,ll+2)=1/l/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Cr2(:,jj,ll)-(1+f2+2*ef)*Br2(:,jj,ll)+l*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll+2));
D(:,ii+1,jj,ll+2)=1/(1-e2)*(1/l*(+e*A(:,ii,ll)-B(:,jj,ll))-(c-de)*yv.*D(:,ii,jj,ll+2));
D(:,ii,jj+1,ll+2)=1/(1-e2)*(1/l*(+e*B(:,jj,ll)-A(:,ii,ll))-(d-ce)*yv.*D(:,ii,jj,ll+2));

%E[0,0,1] E[1,0,1] E[0,1,1]
i=0 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii,kk,ll) = 1/(1-c2)/(l-2-i-k)*((l*(a2+sv2)*(1-c2)+2*c*d*e*l*sv2-l*sv2*(e2+d2)).*Es1(:,ii,kk,ll+2)+(-rfs.^i.*(rfs*(1-c2)-d*sv*c+e*sv).*Cs1(:,kk,ll)-yv.^k.*(yv*(1-c2)+d*sv-e*c*sv).*As1(:,ii,ll)));
tp1=1+h2+2*eh-(dh+c)^2;
Es2(:,ii,kk,ll) = 1/(l-i-k-2)*(((hm+em)*(dh+c)-dm*(1+h2+2*eh))/tp1*As2(:,ii,ll)-yv.^(k+1).*As2(:,ii,ll) + (dm*(dh+c)-(hm+em))/tp1*Cs2(:,kk,ll)-rv.^(i+1).*Cs2(:,kk,ll) +l*(a2+m2)*Es2(:,ii,kk,ll+2) - l*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(1+h2+2*eh)-(dh+c)*(hm+em)))/tp1*Es2(:,ii,kk,ll+2));
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(c*As1(:,ii,ll-2)-Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es1(:,ii,kk+1,ll)=1/(1-c2)*(1/(l-2)*(c*Cs1(:,kk,ll-2)-As1(:,ii,ll-2))+(c*e-d)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*As2(:,ii,ll-2)-Cs2(:,kk,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );
Es2(:,ii,kk+1,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*Cs2(:,kk,ll-2)-(1+h2+2*eh)*As2(:,ii,ll-2)+(l-2)*((dh+c)*(hm+em)-(1+h2+2*eh)*dm)*Es2(:,ii,kk,ll));
Fr1(:,jj,kk,ll) = 1/(1-d2)/(l-2-j-k)*((l*(a2+rv2)*(1-d2)+2*c*d*e*l*rv2-l*rv2*(e2+c2)).*Fr1(:,jj,kk,ll+2)+(-sfr.^j.*(sfr*(1-d2)-d*rv*c+e*rv).*Cr1(:,kk,ll)-yv.^k.*(yv*(1-d2)+c*rv-e*d*rv).*Br1(:,jj,ll)));
tp1=1+f2+2*ef-(cf+d)^2;
Fr2(:,jj,kk,ll) = 1/(l-j-k-2)*(((fg+eg)*(cf+d)-cg*(1+f2+2*ef))/tp1*Br2(:,jj,ll)-yv.^(k+1).*Br2(:,jj,ll) + (cg*(cf+d)-(fg+eg))/tp1*Cr2(:,kk,ll)-sv.^(j+1).*Cr2(:,kk,ll) +l*(a2+g2)*Fr2(:,jj,kk,ll+2) - l*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(1+f2+2*ef)-(cf+d)*(fg+eg)))/tp1*Fr2(:,jj,kk,ll+2));
Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(d*Br1(:,jj,ll-2)-Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Br2(:,jj,ll-2)-Cr2(:,kk,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );
Fr1(:,jj,kk+1,ll)=1/(1-d2)*(1/(l-2)*(d*Cr1(:,kk,ll-2)-Br1(:,jj,ll-2))+(d*e-c)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj,kk+1,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Cr2(:,kk,ll-2)-(1+f2+2*ef)*Br2(:,jj,ll-2)+(l-2)*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll));

%D001
D(:,ii,jj,ll) =1/((1-e2)*(l-2-i-j))*((l*(a2+yv2)*(1-e2)+2*c*d*e*l*yv2-l*yv2*(c2+d2)).*D(:,ii,jj,ll+2)-(yv*(d-e*c)+m.^(j+1)*(1-e2)).*As2(:,ii,ll)-(h*(1-e2)).*As2(:,ii+1,ll)-(yv*(d-e*c)+sv.^(j+1)*(1-e2)).*As1(:,ii,ll)-(yv*(c-e*d)+g.^(i+1)*(1-e2)).*Br2(:,jj,ll)-(f*(1-e2)).*Br2(:,jj+1,ll)-(yv*(c-e*d)+rv.^(i+1)*(1-e2)).*Br1(:,jj,ll));

clear test0 test 
test0(1:8,1)=(l*(a2)*(1-e2)).*D(1:8,ii,jj,ll+2);
test0(1:8,2)=(+2*c*d*e*l*yv2(1:8)).*D(1:8,ii,jj,ll+2);
test0(1:8,3)=(-l*yv2(1:8)*(c2+d2)).*D(1:8,ii,jj,ll+2);
test0(1:8,4)=(l*(yv2(1:8))*(1-e2)).*D(1:8,ii,jj,ll+2);
test(1:8,1)= SmartSum(test0);
test(:,2)= -(yv*(d-ce)+m*(1-e2)).*As2(:,ii,ll);
test(:,3)= -(h*(1-e2)).*As2(:,ii+1,ll);
test(:,4)= -(yv*(d-ce)+sv*(1-e2)).*As1(:,ii,ll);
test(:,5)= -(yv*(c-de)+g*(1-e2)).*Br2(:,jj,ll);
test(:,6)= -(f*(1-e2)).*Br2(:,jj+1,ll);
test(:,7)= -(yv*(c-de)+r1*(1-e2)).*Br1(:,jj,ll);
vtest(1:8) = 1/((1-e2)*(l-2))*SmartSum(test);
D(:,ii,jj,ll) = vtest(1:8);

D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(+e*A(:,ii,ll-2)-B(:,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(+e*B(:,jj,ll-2)-A(:,ii,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));

%DEF[0,0,-1]
i=0 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii,kk,ll) = 1/(1-c2)/(l-2-i-k)*((l*(a2+sv2)*(1-c2)+2*c*d*e*l*sv2-l*sv2*(e2+d2)).*Es1(:,ii,kk,ll+2)+(-rfs.^i.*(rfs*(1-c2)-d*sv*c+e*sv).*Cs1(:,kk,ll)-yv.^k.*(yv*(1-c2)+d*sv-e*c*sv).*As1(:,ii,ll)));
tp1=1+h2+2*eh-(dh+c)^2;
Es2(:,ii,kk,ll) = 1/(l-i-k-2)*(((hm+em)*(dh+c)-dm*(1+h2+2*eh))/tp1*As2(:,ii,ll)-yv.^(k+1).*As2(:,ii,ll) + (dm*(dh+c)-(hm+em))/tp1*Cs2(:,kk,ll)-rv.^(i+1).*Cs2(:,kk,ll) +l*(a2+m2)*Es2(:,ii,kk,ll+2) - l*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(1+h2+2*eh)-(dh+c)*(hm+em)))/tp1*Es2(:,ii,kk,ll+2));
Fr1(:,jj,kk,ll) = 1/(1-d2)/(l-2-j-k)*((l*(a2+rv2)*(1-d2)+2*c*d*e*l*rv2-l*rv2*(e2+c2)).*Fr1(:,jj,kk,ll+2)+(-sfr.^j.*(sfr*(1-d2)-d*rv*c+e*rv).*Cr1(:,kk,ll)-yv.^k.*(yv*(1-d2)+c*rv-e*d*rv).*Br1(:,jj,ll)));
tp1=1+f2+2*ef-(cf+d)^2;
Fr2(:,jj,kk,ll) = 1/(l-j-k-2)*(((fg+eg)*(cf+d)-cg*(1+f2+2*ef))/tp1*Br2(:,jj,ll)-yv.^(k+1).*Br2(:,jj,ll) + (cg*(cf+d)-(fg+eg))/tp1*Cr2(:,kk,ll)-sv.^(j+1).*Cr2(:,kk,ll) +l*(a2+g2)*Fr2(:,jj,kk,ll+2) - l*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(1+f2+2*ef)-(cf+d)*(fg+eg)))/tp1*Fr2(:,jj,kk,ll+2));

D(:,ii,jj,ll) =1/((1-e2)*(l-2))*((l*(a2+yv2)*(1-e2)+2*c*d*e*l*yv2-l*yv2*(c2+d2)).*D(:,ii,jj,ll+2)-(yv*(d-ce)+m*(1-e2)).*As2(:,ii,ll) ...
    -(h*(1-e2))*As2(:,ii+1,ll)-(yv*(d-ce)+s1*(1-e2)).*As1(:,ii,ll)-(yv*(c-de)+g*(1-e2)).*Br2(:,jj,ll)-(f*(1-e2))*Br2(:,jj+1,ll) ...
    -(yv*(c-de)+r1*(1-e2)).*Br1(:,jj,ll));

clear test 
test(1:8,1)= (l*(a2+yv2(1:8))*(1-e2)+2*c*d*e*l*yv2(1:8)-l*yv2(1:8)*(c2+d2)).*D(1:8,ii,jj,ll+2);
test(:,2)= -(yv*(d-ce)+m*(1-e2)).*As2(:,ii,ll);
test(:,3)= -(h*(1-e2)).*As2(:,ii+1,ll);
test(:,4)= -(yv*(d-ce)+sv*(1-e2)).*As1(:,ii,ll);
test(:,5)= -(yv*(c-de)+g*(1-e2)).*Br2(:,jj,ll);
test(:,6)= -(f*(1-e2)).*Br2(:,jj+1,ll);
test(:,7)= -(yv*(c-de)+r1*(1-e2)).*Br1(:,jj,ll);
vtest(1:8) = 1/((1-e2)*(l-2))*SmartSum(test);
D(:,ii,jj,ll) = vtest(1:8);

%EF[0,0,-3] %QUAD
i=0 ; j=0 ; k=0 ; l=-3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
E03s1(1:8) = 1/(1-c2)/(l-2-i-k)*((l*(a2+sv2)*(1-c2)+2*c*d*e*l*sv2-l*sv2*(e2+d2)).*Es1(:,ii,kk,ll+2)+(-rfs.^i.*(rfs*(1-c2)-d*sv*c+e*sv).*C03s1(:)-yv.^k.*(yv*(1-c2)+d*sv-e*c*sv).*A03s1(:)));
tp1=1+h2+2*eh-(dh+c)^2;
E03s2(1:8) = 1/(l-i-k-2)*(((hm+em)*(dh+c)-dm*(1+h2+2*eh))/tp1*A03s2(:)-yv.^(k+1).*A03s2(:) + (dm*(dh+c)-(hm+em))/tp1*C03s2(:)-rv.^(i+1).*C03s2(:) +l*(a2+m2)*Es2(:,ii,kk,ll+2) - l*((hm+em)*(hm+em-dm*(dh+c))+dm*(dm*(1+h2+2*eh)-(dh+c)*(hm+em)))/tp1*Es2(:,ii,kk,ll+2));

F03r1(1:8) = 1/(1-d2)/(l-2-j-k)*((l*(a2+rv2)*(1-d2)+2*c*d*e*l*rv2-l*rv2*(e2+c2)).*Fr1(:,jj,kk,ll+2)+(-sfr.^j.*(sfr*(1-d2)-c*rv*d+e*rv).*C03r1(:)-yv.^k.*(yv*(1-d2)+c*rv-e*d*rv).*B03r1(:)));
tp1=1+f2+2*ef-(cf+d)^2;
F03r2(1:8) = 1/(l-j-k-2)*(((fg+eg)*(cf+d)-cg*(1+f2+2*ef))/tp1*B03r2(:)-yv.^(k+1).*B03r2(:) + (cg*(cf+d)-(fg+eg))/tp1*C03r2(:)-sv.^(j+1).*C03r2(:) +l*(a2+g2)*Fr2(:,jj,kk,ll+2) - l*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(1+f2+2*ef)-(cf+d)*(fg+eg)))/tp1*Fr2(:,jj,kk,ll+2));

%D03(1:8)
%=1/((1-e2)*(l-2-i-j))*((l*(a2+yv2)*(1-e2)+2*c*d*e*l*yv2-l*yv2*(c2+d2)).*D(:,ii,jj,ll+2)-(yv*(d-e*c)+m.^(j+1)*(1-e2)).*A03s2(:)-(h*(1-e2)).*A13s2(:)-(yv*(d-e*c)+sv.^(j+1)*(1-e2)).*A03s1(:)-(yv*(c-e*d)+g.^(i+1)*(1-e2)).*B03r2(:)-(f*(1-e2)).*B13r2(:)-(yv*(c-e*d)+rv.^(i+1)*(1-e2)).*B03r1(:));
%D03(1:8)
%=1/((1-e2)*(l-2-i-j))*(l*(a2*(1-e2)+(1-c2-d2-e2+2*c*d*e)*yv2).*D(:,ii,jj,ll+2)-(m.^(j+1)*(1-e2)+yv*(d-ce)).*A03s2(:)-(h*(1-e2)).*A13s2(:)-(sv.^(j+1)*(1-e2)+yv*(d-ce)).*A03s1(:)-(g.^(i+1)*(1-e2)+yv*(c-de)).*B03r2(:)-(f*(1-e2)).*B13r2(:)-(rv.^(i+1)*(1-e2)+yv*(c-de)).*B03r1(:));
clear test 
test(1:8,1)= l*(a2*(1-e2)+(1-c2-d2-e2+2*c*d*e)*yv2).*D(:,ii,jj,ll+2);
test(:,2)=-(m.^(j+1)*(1-e2)+yv*(d-ce)).*A03s2(:);
test(:,3)=-(h*(1-e2)).*A13s2(:);
test(:,4)=-(sv.^(j+1)*(1-e2)+yv*(d-ce)).*A03s1(:);
test(:,5)=-(g.^(i+1)*(1-e2)+yv*(c-de)).*B03r2(:);
test(:,6)=-(f*(1-e2)).*B13r2(:);
test(:,7)=-(rv.^(i+1)*(1-e2)+yv*(c-de)).*B03r1(:);
D03(1:8) =1/((1-e2)*(l-2-i-j))*SmartSum(test);

% F03r1(:,jj,kk,ll) = 1/(1-d2)/(l-2-j-k)*((l*(a2+rv2)*(1-d2)+2*c*d*e*l*rv2-l*rv2*(e2+c2)).*Fr1(:,jj,kk,ll+2)+(-sfr.^j.*(sfr*(1-d2)-d*rv*c+e*rv).*Cr1(:,kk,ll)-yv.^k.*(yv*(1-d2)+c*rv-e*d*rv).*Br1(:,jj,ll)));
% tp1=1+f2+2*ef-(cf+d)^2;
% F03r2(:,jj,kk,ll) = 1/(l-j-k-2)*(((fg+eg)*(cf+d)-cg*(1+f2+2*ef))/tp1*Br2(:,jj,ll)-yv.^(k+1).*Br2(:,jj,ll) + (cg*(cf+d)-(fg+eg))/tp1*Cr2(:,kk,ll)-sv.^(j+1).*Cr2(:,kk,ll) +l*(a2+g2)*Fr2(:,jj,kk,ll+2) - l*((fg+eg)*(fg+eg-cg*(cf+d))+cg*(cg*(1+f2+2*ef)-(cf+d)*(fg+eg)))/tp1*Fr2(:,jj,kk,ll+2));

%EF[1,0,-1] %QUAD
i=0 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(-yv.^k.*A03s1(:)+c*r1.^i.*C03s1(:))+(ce-d)*sv.*Es1(:,ii,kk,ll));
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(+c*yv.^k.*A03s1(:)-rfs.^i.*C03s1(:))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*A03s2(:)-rv.^i.*C03s2(:)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );

%Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(-yv.^k.*B03r1(:)+d*s1.^j.*C03r1(:))+(de-c)*rv.*Fr1(:,jj,kk,ll));
Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(+d*yv.^k.*B03r1(:)-sfr.^j.*C03r1(:))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*B03r2(:)-sv.^j.*C03r2(:)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );

A03=A03s1+A03s2;
B03=B03r1+B03r2;
D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(+e*A03(:)-B03(:))-(c-d*e)*yv.*D(:,ii,jj,ll));
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(+e*B03(:)-A03(:))-(d-c*e)*yv.*D(:,ii,jj,ll));

%E[2,0,1] E[1,1,1] QUAD
i=1 ; j=1 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(i*Es1(:,ii-1,kk,ll-2)+c*yv.^k.*As1(:,ii,ll-2)-rfs.^i.*Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*As2(:,ii,ll-2)-rv.^i.*Cs2(:,kk,ll-2)+i*Es2(:,ii-1,kk,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );

Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(j*Fr1(:,jj-1,kk,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-sfr.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)+j*Fr2(:,jj-1,kk,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,ii,kk,ll) );

i=0 ; j=0 ; k=1 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(-c*k*Es1(:,ii,kk-1,ll-2)+c*yv.^k.*As1(:,ii,ll-2)-rfs.^i.*Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*As2(:,ii,ll-2)-rv.^i.*Cs2(:,kk,ll-2)-k*(dh+c)*Es2(:,ii,kk-1,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );

Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(-d*k*Fr1(:,jj,kk-1,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-sfr.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)-k*(cf+d)*Fr2(:,jj,kk-1,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );

%D[2,0,1] D[0,2,1] D[1,1,1] QUAD
i=1 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Bri(:,2,1,ll-2)= f*Br2(:,2,ll-2) + g*Br2(:,1,ll-2) + r1*Br1(:,1,ll-2);
Bri(:,2,2,ll-2)= f*Br2(:,3,ll-2) + g*Br2(:,2,ll-2) + r1*Br1(:,2,ll-2);
Bri(:,1,2,ll-2)= B(:,2,ll-2);
Asi(:,1,2,ll-2)= h*As2(:,2,ll-2) + m*As2(:,1,ll-2) + s1*As1(:,1,ll-2);
Asi(:,2,2,ll-2)= h*As2(:,3,ll-2) + m*As2(:,2,ll-2) + s1*As1(:,2,ll-2);
Asi(:,2,1,ll-2)= A(:,2,ll-2);
%D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)+e*A(:,ii,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));

loverl= 1/(l-2);

clear test;
test(1:8,1)=(loverl*(i*D(:,ii-1,jj,ll-2)));
test(:,2)=(1/(l-2)*e*A(:,ii,ll-2));
test(:,3)=(-1/(l-2)*Bri(:,ii,jj,ll-2));
test(:,4)=(-(c-de)*yv.*D(:,ii,jj,ll));
vtest(1:8) = 1/(1-e2)*SmartSum(test);
D(:,ii+1,jj,ll) = vtest(:);

%D111_1(1:8)= 1/(1-e2)*(1/(l-2)*(-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll))

clear test;
test(1:8,1)=(loverl*(-e*i*D(:,ii-1,jj,ll-2)));
test(:,2)=(1/(l-2)*e*Bri(:,ii,jj,ll-2));
test(:,3)=(-1/(l-2)*Asi(:,ii,jj,ll-2));
test(:,4)=(-(d-ce)*yv.*D(:,ii,jj,ll));
vtest(1:8) = 1/(1-e2)*SmartSum(test);
D111_1(1:8) = vtest(1:8);

i=0 ; j=1 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)+e*B(:,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));
%D111_2(1:8) =1/(1-e2)*(1/(l-2)*(-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll))

clear test;
test(1:8,1)=(loverl*(j*D(:,ii,jj-1,ll-2)));
test(:,2)=(1/(l-2)*e*B(:,jj,ll-2));
test(:,3)=(-1/(l-2)*Asi(:,ii,jj,ll-2));
test(:,4)=(-(d-ce)*yv.*D(:,ii,jj,ll));
vtest(1:8) = 1/(1-e2)*SmartSum(test);
D(:,ii,jj+1,ll) = vtest(:);

clear test;
test(1:8,1)=(loverl*(-e*j*D(:,ii,jj-1,ll-2)));
test(:,2)=(1/(l-2)*e*Asi(:,ii,jj,ll-2));
test(:,3)=(-1/(l-2)*Bri(:,ii,jj,ll-2));
test(:,4)=(-(c-de)*yv.*D(:,ii,jj,ll));
vtest(1:8) = 1/(1-e2)*SmartSum(test);
D111_2(1:8)= vtest(1:8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D(:,2,2,ll)=0.5*D111_1(1:8)+0.5*D111_2(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%EF[2,0,3]
i=1 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(i*Es1(:,ii-1,kk,ll-2)+c*As1(:,ii,ll-2)-rv.^i.*Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*As2(:,ii,ll-2)-rv.^i.*Cs2(:,kk,ll-2)+i*Es2(:,ii-1,kk,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );
Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(j*Fr1(:,jj-1,kk,ll-2)+d*Br1(:,jj,ll-2)-sfr.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)+j*Fr2(:,jj-1,kk,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );

%D[2,0,3] D[0,2,3]
i=1 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Bri(:,2,1,ll-2)= f*Br2(:,2,ll-2) + g*Br2(:,1,ll-2) + r1*Br1(:,1,ll-2);
Bri(:,2,2,ll-2)= f*Br2(:,3,ll-2) + g*Br2(:,2,ll-2) + r1*Br1(:,2,ll-2);
Bri(:,1,2,ll-2)= B(:,2,ll-2);
Asi(:,1,2,ll-2)= h*As2(:,2,ll-2) + m*As2(:,1,ll-2) + s1*As1(:,1,ll-2);
Asi(:,2,2,ll-2)= h*As2(:,3,ll-2) + m*As2(:,2,ll-2) + s1*As1(:,2,ll-2);
Asi(:,2,1,ll-2)= A(:,2,ll-2);
D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)+e*A(:,ii,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear test
test(:,1)= 1/(l-2)*(i*D(:,ii-1,jj,ll-2));
test(:,2)= 1/(l-2)*(+e*A(:,ii,ll-2));
test(:,3)= 1/(l-2)*(-Bri(:,ii,jj,ll-2));
test(:,4)= -(c-de)*yv.*D(:,ii,jj,ll);
vtest= 1/(1-e2)*SmartSum(test);
D(:,ii+1,jj,ll) = vtest;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=0 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)+e*B(:,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));

%%%%%%%%%%%%%%%%%%%%%%%
clear test
test(:,1)= 1/(l-2)*(j*D(:,ii,jj-1,ll-2));
test(:,2)= 1/(l-2)*(e*B(:,jj,ll-2));
test(:,3)= 1/(l-2)*(-Asi(:,ii,jj,ll-2));
test(:,4)= -(d-ce)*yv.*D(:,ii,jj,ll);
vtest= 1/(1-e2)*SmartSum(test);
D(:,ii,jj+1,ll) = vtest;
%%%%%%%%%%%%%%%%%%%%%%%

%E[0,2,3]
i=0 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii,kk+1,ll)=1/(1-c2)*(1/(l-2)*(k*Es1(:,ii,kk-1,ll-2)+c*Cs1(:,kk,ll-2)-yv.^k.*As1(:,ii,ll-2))+(c*e-d)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii,kk+1,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*Cs2(:,kk,ll-2)-(1+h2+2*eh)*yv.^k.*As2(:,ii,ll-2)+k*(1+h2+2*eh)*Es2(:,ii,kk-1,ll-2)+(l-2)*((dh+c)*(hm+em)-(1+h2+2*eh)*dm)*Es2(:,ii,kk,ll));
Fr1(:,jj,kk+1,ll)=1/(1-d2)*(1/(l-2)*(k*Fr1(:,jj,kk-1,ll-2)+d*Cr1(:,kk,ll-2)-yv.^k.*Br1(:,jj,ll-2))+(d*e-c)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj,kk+1,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*Cr2(:,kk,ll-2)-(1+f2+2*ef)*yv.^k.*Br2(:,jj,ll-2)+k*(1+f2+2*ef)*Fr2(:,jj,kk-1,ll-2)+(l-2)*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll));

%E[1,1,3] F[1,1,3] D[1,1,3]
i=0 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(-c*k*Es1(:,ii,kk-1,ll-2)+c*yv.^k.*As1(:,ii,ll-2)-Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*As2(:,ii,ll-2)-Cs2(:,kk,ll-2)-k*(dh+c)*Es2(:,ii,kk-1,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );
Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(-d*k*Fr1(:,jj,kk-1,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-Cr2(:,kk,ll-2)-k*(cf+d)*Fr2(:,jj,kk-1,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );
i=0 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-d*e)*yv.*D(:,ii,jj,ll));

i=0 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-d*e)*yv.*D(:,ii,jj,ll));
%%%%
D113_1(1:8) =1/(1-e2)*(1/(l-2)*(-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-d*e)*yv.*D(:,ii,jj,ll));


i=1 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-c*e)*yv.*D(:,ii,jj,ll));
D113_2(1:8) =1/(1-e2)*(1/(l-2)*(-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-c*e)*yv.*D(:,ii,jj,ll));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D(:,ii,jj+1,ll)= 0.5*D113_1(1:8) + 0.5*D113_2(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEF[2,1,3] DEF[1,2,3]
i=1 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(i*Es1(:,ii-1,kk,ll-2)-c*k*Es1(:,ii,kk-1,ll-2)+c*yv.^k.*As1(:,ii,ll-2)-rfs.^i.*Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es1(:,ii,kk+1,ll)=1/(1-c2)*(1/(l-2)*(k*Es1(:,ii,kk-1,ll-2)-c*i*Es1(:,ii-1,kk,ll-2)+c*rfs.^i.*Cs1(:,kk,ll-2)-yv.^k.*As1(:,ii,ll-2))+(c*e-d)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*As2(:,ii,ll-2)-rv.^i.*Cs2(:,kk,ll-2)+i*Es2(:,ii-1,kk,ll-2)-k*(dh+c)*Es2(:,ii,kk-1,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );
Es2(:,ii,kk+1,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*rv.^i.*Cs2(:,kk,ll-2)-(1+h2+2*eh)*yv.^k.*As2(:,ii,ll-2)+k*(1+h2+2*eh)*Es2(:,ii,kk-1,ll-2)-i*(dh+c)*Es2(:,ii-1,kk,ll-2)+(l-2)*((dh+c)*(hm+em)-(1+h2+2*eh)*dm)*Es2(:,ii,kk,ll));
Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(j*Fr1(:,jj-1,kk,ll-2)-d*k*Fr1(:,jj,kk-1,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-sv.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)+j*Fr2(:,jj-1,kk,ll-2)-k*(cf+d)*Fr2(:,jj,kk-1,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );
Fr1(:,jj,kk+1,ll)=1/(1-d2)*(1/(l-2)*(k*Fr1(:,jj,kk-1,ll-2)-d*j*Fr1(:,jj-1,kk,ll-2)+d*sv.^j.*Cr1(:,kk,ll-2)-yv.^k.*Br1(:,jj,ll-2))+(d*e-c)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj,kk+1,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*sv.^j.*Cr2(:,kk,ll-2)-(1+f2+2*ef)*yv.^k.*Br2(:,jj,ll-2)+k*(1+f2+2*ef)*Fr2(:,jj,kk-1,ll-2)-j*(cf+d)*Fr2(:,jj-1,kk,ll-2)+(l-2)*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll));
D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-d*e)*yv.*D(:,ii,jj,ll));
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-c*e)*yv.*D(:,ii,jj,ll));

%DEF[3,1,3] DEF[2,2,3] %QUAD
i=2 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii+1,kk,ll)=1/(1-c2)*(1/(l-2)*(i*Es1(:,ii-1,kk,ll-2)-c*k*Es1(:,ii,kk-1,ll-2)+c*yv.^k.*As1(:,ii,ll-2)-rfs.^i.*Cs1(:,kk,ll-2))+(c*d-e)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii+1,kk,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*yv.^k.*As2(:,ii,ll-2)-rv.^i.*Cs2(:,kk,ll-2)+i*Es2(:,ii-1,kk,ll-2)-k*(dh+c)*Es2(:,ii,kk-1,ll-2)-(l-2)*(hm+em-dm*(dh+c))*Es2(:,ii,kk,ll) );

Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(j*Fr1(:,jj-1,kk,ll-2)-d*k*Fr1(:,jj,kk-1,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-sfr.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)+j*Fr2(:,jj-1,kk,ll-2)-k*(cf+d)*Fr2(:,jj,kk-1,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );

i=2 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Es1(:,ii,kk+1,ll)=1/(1-c2)*(1/(l-2)*(k*Es1(:,ii,kk-1,ll-2)-c*i*Es1(:,ii-1,kk,ll-2)+c*rfs.^i.*Cs1(:,kk,ll-2)-yv.^k.*As1(:,ii,ll-2))+(c*e-d)*sv.*Es1(:,ii,kk,ll));
Es2(:,ii,kk+1,ll)=1/(l-2)/(1+h2+2*eh-(dh+c)^2)*((dh+c)*rv.^i.*Cs2(:,kk,ll-2)-(1+h2+2*eh)*yv.^k.*As2(:,ii,ll-2)+k*(1+h2+2*eh)*Es2(:,ii,kk-1,ll-2)-i*(dh+c)*Es2(:,ii-1,kk,ll-2)+(l-2)*((dh+c)*(hm+em)-(1+h2+2*eh)*dm)*Es2(:,ii,kk,ll));

Fr1(:,jj,kk+1,ll)=1/(1-d2)*(1/(l-2)*(k*Fr1(:,jj,kk-1,ll-2)-d*j*Fr1(:,jj-1,kk,ll-2)+d*sfr.^j.*Cr1(:,kk,ll-2)-yv.^k.*Br1(:,jj,ll-2))+(d*e-c)*rv.*Fr1(:,jj,kk,ll));
Fr2(:,jj,kk+1,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*sv.^j.*Cr2(:,kk,ll-2)-(1+f2+2*ef)*yv.^k.*Br2(:,jj,ll-2)+k*(1+f2+2*ef)*Fr2(:,jj,kk-1,ll-2)-j*(cf+d)*Fr2(:,jj-1,kk,ll-2)+(l-2)*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll));

i=2 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Asi(:,3,2,3)= h*As2(:,4,3) + m*As2(:,3,3) + s1*As1(:,3,3);
Bri(:,3,2,3)= f2*Br2(:,4,3) + 2*fg*Br2(:,3,3) + g2*Br2(:,2,3) + r1^2*Br1(:,2,3);
D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));

%%%%%
%D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));
D223_1(1:8)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));

clear test
test(:,1)= 1/(l-2)*(j*D(:,ii,jj-1,ll-2));
test(:,2)= -1/(l-2)*e*i*D(:,ii-1,jj,ll-2);
test(:,3)= +1/(l-2)*e*Bri(:,ii,jj,ll-2) ;
test(:,4)= -1/(l-2)*Asi(:,ii,jj,ll-2);
test(:,5)= -(d-ce)*yv.*D(:,ii,jj,ll);
vtest(1:8) = 1/(1-e2)*SmartSum(test);
D(:,ii,jj+1,ll) = vtest(1:8);

i=1 ; j=2 ; k=1 ; l=3 ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Asi(:,2,3,3)= h2*As2(:,4,3) + 2*hm*As2(:,3,3) + m2*As2(:,2,3) + s1^2*As1(:,2,3);
Bri(:,2,3,3)= f*Br2(:,4,3) + g*Br2(:,3,3) + r1*Br1(:,3,3);
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)-e*i*D(:,ii-1,jj,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));

%%%%%
D223_2(1:8)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)-e*j*D(:,ii,jj-1,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D(:,3,3,ll) = 0.5*D223_1(1:8) +0.5*D223_2(1:8);
%D(:,3,3,ll) = D223_1(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%D[3,0,3] %QUAD
i=2 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Asi(:,3,1,3)= As2(:,3,3) + As1(:,3,3);
Bri(:,3,1,3)= f2*Br2(:,3,3) + 2*fg*Br2(:,2,3) + g2*Br2(:,1,3) + r1^2*Br1(:,1,3);

%>>>>>>>>>>>>
% temp(1:8)=Asi(1:8,ii,jj,ll-2);
% testAsi=temp(1:8)*signv'
% F22=@(r) r.^(ii-1).*(h*r+m).^(jj-1).*(r.^2+(h*r+m).^2+y2^2+2*r.*(c*y2+e*(h*r+m))+2*d*(h*r+m)*y2+a2).^(-1/2);
% F21=@(r) r.^(ii-1).*(h*r+m).^(jj-1).*(r.^2+(h*r+m).^2+y1^2+2*r.*(c*y1+e*(h*r+m))+2*d*(h*r+m)*y1+a2).^(-1/2);
% F12=@(r) r.^(ii-1).*s1^(jj-1).*(r.^2+s1^2+y2^2+2*r*(c*y2+e*s1)+2*d*s1*y2+a2).^(-1/2);
% F11=@(r) r.^(ii-1).*s1^(jj-1).*(r.^2+s1^2+y1^2+2*r*(c*y1+e*s1)+2*d*s1*y1+a2).^(-1/2);
% NIAsi=quad(F22,r1,r2)-quad(F21,r1,r2)-quad(F12,r1,r2)+quad(F11,r1,r2)
% temp(1:8)=Bri(1:8,ii,jj,ll-2);
% testBri=temp(1:8)*signv'
% F22=@(s) (f*s+g).^(ii-1).*s.^(jj-1).*(s.^2+(f*s+g).^2+y2^2+2*s.*(d*y2+e*(f*s+g))+2*c*(f*s+g)*y2+a2).^(-(l-2)/2);
% F21=@(s) (f*s+g).^(ii-1).*s.^(jj-1).*(s.^2+(f*s+g).^2+y1^2+2*s.*(d*y1+e*(f*s+g))+2*c*(f*s+g)*y1+a2).^(-(l-2)/2);
% F12=@(s) r1.^(ii-1).*s.^(jj-1).*(s.^2+r1^2+y2^2+2*s*(d*y2+e*r1)+2*c*r1*y2+a2).^(-(l-2)/2);
% F11=@(s) r1.^(ii-1).*s.^(jj-1).*(s.^2+r1^2+y1^2+2*s*(d*y1+e*r1)+2*c*r1*y1+a2).^(-(l-2)/2);
% NIBri=quad(F22,s1,s2)-quad(F21,s1,s2)-quad(F12,s1,s2)+quad(F11,s1,s2)
%%pause

D(:,ii+1,jj,ll)=1/(1-e2)*(1/(l-2)*(i*D(:,ii-1,jj,ll-2)+e*Asi(:,ii,jj,ll-2)-Bri(:,ii,jj,ll-2))-(c-de)*yv.*D(:,ii,jj,ll));
i=0 ; j=2 ; k=1 ; l=3 ;  ii=i+1;jj=j+1;kk=k+1;ll=l+2;
Asi(:,1,3,3)= h2*As2(:,3,3) + 2*hm*As2(:,2,3) + m2*As2(:,1,3) + s1^2*As1(:,1,3);
Bri(:,1,3,3)= Br2(:,3,3) + Br1(:,3,3);
D(:,ii,jj+1,ll)=1/(1-e2)*(1/(l-2)*(j*D(:,ii,jj-1,ll-2)+e*Bri(:,ii,jj,ll-2)-Asi(:,ii,jj,ll-2))-(d-ce)*yv.*D(:,ii,jj,ll));

% % DEBUGGING D INTEGRALS
% vD00m3(1:8)=D(1:8,1,1,5);
% vD10m3(1:8)=D(1:8,2,1,5);
% vD01m3(1:8)=D(1:8,1,2,5);
% vD00m1(1:8)=D(1:8,1,1,3);
% vD10m1(1:8)=D(1:8,2,1,3);
% vD01m1(1:8)=D(1:8,1,2,3);
% vD20m3(1:8)=D(1:8,3,1,5);
% vD02m3(1:8)=D(1:8,1,3,5);
% vD11m3(1:8)=D(1:8,2,2,5);
% vD21m3(1:8)=D(1:8,3,2,5);
% vD12m3(1:8)=D(1:8,2,3,5);
% vD001(1:8)=D(1:8,1,1,1);
% vD003(1:8)=D03(1:8);
% vD101(1:8)=D(1:8,2,1,1);
% vD011(1:8)=D(1:8,1,2,1);
% vD20m1(1:8)=D(1:8,3,1,3);
% vD11m1(1:8)=D(1:8,2,2,3);
% vD02m1(1:8)=D(1:8,1,3,3);
% vD31m3(1:8)=D(1:8,4,2,5);
% vD22m3(1:8)=D(1:8,3,3,5);
% vD13m3(1:8)=D(1:8,2,4,5);
% vD30m3(1:8)=D(1:8,4,1,5);
% vD03m3(1:8)=D(1:8,1,4,5);

% scD00m3=SmartDot(vD00m3, signv)
% scD10m3=SmartDot(vD10m3, signv)
% scD01m3=SmartDot(vD01m3, signv)
% scD00m1=SmartDot(vD00m1, signv)
% scD10m1=SmartDot(vD10m1, signv)
% scD01m1=SmartDot(vD01m1, signv)
% scD20m3=SmartDot(vD20m3, signv)
% scD02m3=SmartDot(vD02m3, signv)
% scD11m3=SmartDot(vD11m3, signv)
% scD21m3=SmartDot(vD21m3, signv)
% scD12m3=SmartDot(vD12m3, signv)
% scD001=SmartDot(vD001, signv)
% scD003=SmartDot(D03, signv)
% scD101=SmartDot(vD101, signv)
% scD011=SmartDot(vD011, signv)
% scD20m1=SmartDot(vD20m1, signv)
% scD11m1=SmartDot(vD11m1, signv)
% scD02m1=SmartDot(vD02m1, signv)
% scD31m3=SmartDot(vD31m3, signv)
% scD22m3=SmartDot(vD22m3, signv)
% scD13m3=SmartDot(vD13m3, signv)
% scD30m3=SmartDot(vD30m3, signv)
% scD03m3=SmartDot(vD03m3, signv)

% nbin=3000;
% [xi,wi] = gqwp(nbin);
% AA=(r2-r1)/2;
% BB=(r2+r1)/2;
% CC=(s2-s1)/2;
% DD=(s2+s1)/2;
% drodx=AA;
% dsodx=CC;
% NIDy2(1:5,1:5,1:7)=0;
% NIDy1(1:5,1:5,1:7)=0;
% for ir=1:nbin
%     r=AA*xi(ir)+BB;
%     sf=-lqlp*r+s0;
%     df2(1:5,1:5,1:7)=0;
%     df1(1:5,1:5,1:7)=0;
%     for is=1:nbin
%         s=CC*xi(is)+DD;
%         if (s <= sf)
%             l=1; ipjmax=2;
%             for (i=2:2)
%                 for   j=1:1
%                     ii=i+1;jj=j+1;
%                     F1(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y1^2+2*e*r*s+2*y1*(c*r+d*s)+a2)^(-(l)/2));
%                     df1(ii,jj,l+2)=df1(ii,jj,l+2)+wi(is)*F1(ii,jj,l+2)*dsodx;
%                     F2(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y2^2+2*e*r*s+2*y2*(c*r+d*s)+a2)^(-(l)/2));
%                     df2(ii,jj,l+2)=df2(ii,jj,l+2)+wi(is)*F2(ii,jj,l+2)*dsodx;
%                 end
%             end
%         end
%     end
%     NIDy2(:,:,:)=NIDy2(:,:,:)+wi(ir)*df2(:,:,:)*drodx;
%     NIDy1(:,:,:)=NIDy1(:,:,:)+wi(ir)*df1(:,:,:)*drodx;
% end
% NID(:,:,:)=NIDy2(:,:,:)-NIDy1(:,:,:);
% NID(2,2,3)
% vD111(1:8) = D(:,2,2,3);
% vD111(1:8) *signv'
% 
% %pause

% %We use our own Gausse Quadrature, the native Matlab quadrature cannot be
% %used as integral bounds of the second integral depends on the first
% %variable
% clear F1 F2
% nbin=1000;
% [xi,wi] = gqwp(nbin);
% AA=(r2-r1)/2;
% BB=(r2+r1)/2;
% CC=(s2-s1)/2;
% DD=(s2+s1)/2;
% drodx=AA;
% dsodx=CC;
% NIDy2(1:5,1:5,1:7)=0;
% NIDy1(1:5,1:5,1:7)=0;
% for ir=1:nbin
%     r=AA*xi(ir)+BB;
%     sf=-lqlp*r+s0;
%     df2(1:5,1:5,1:7)=0;
%     df1(1:5,1:5,1:7)=0;
%     for is=1:nbin
%         s=CC*xi(is)+DD;
%         if (s <= sf)
%             l=3; ipjmax=4;
%             for (i=0:ipjmax)
%                 for   j=0:ipjmax-i
%                     ii=i+1;jj=j+1;
%                     F1(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y1^2+2*e*r*s+2*y1*(c*r+d*s)+a2)^(-(l)/2));
%                     df1(ii,jj,l+2)=df1(ii,jj,l+2)+wi(is)*F1(ii,jj,l+2)*dsodx;
%                     F2(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y2^2+2*e*r*s+2*y2*(c*r+d*s)+a2)^(-(l)/2));
%                     df2(ii,jj,l+2)=df2(ii,jj,l+2)+wi(is)*F2(ii,jj,l+2)*dsodx;
%                 end
%             end
%             l=1; ipjmax=2;
%             for (i=0:ipjmax)
%                 for   j=0:ipjmax-i
%                     ii=i+1;jj=j+1;
%                     F1(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y1^2+2*e*r*s+2*y1*(c*r+d*s)+a2)^(-(l)/2));
%                     df1(ii,jj,l+2)=df1(ii,jj,l+2)+wi(is)*F1(ii,jj,l+2)*dsodx;
%                     F2(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y2^2+2*e*r*s+2*y2*(c*r+d*s)+a2)^(-(l)/2));
%                     df2(ii,jj,l+2)=df2(ii,jj,l+2)+wi(is)*F2(ii,jj,l+2)*dsodx;
%                 end
%             end
%             l=-1; ipjmax=2;
%             for (i=0:ipjmax)
%                 for   j=0:ipjmax-i
%                     ii=i+1;jj=j+1;
%                     F1(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y1^2+2*e*r*s+2*y1*(c*r+d*s)+a2)^(-(l)/2));
%                     df1(ii,jj,l+2)=df1(ii,jj,l+2)+wi(is)*F1(ii,jj,l+2)*dsodx;
%                     F2(ii,jj,l+2)= r^(i)*s^(j)*((r^2+s^2+y2^2+2*e*r*s+2*y2*(c*r+d*s)+a2)^(-(l)/2));
%                     df2(ii,jj,l+2)=df2(ii,jj,l+2)+wi(is)*F2(ii,jj,l+2)*dsodx;
%                 end
%             end
%             l=-3; ipjmax=0;
%             for (i=0:ipjmax)
%                 for   j=0:ipjmax-i
%                     ii=i+1;jj=j+1;
%                     F1(ii,jj,7)= r^(i)*s^(j)*((r^2+s^2+y1^2+2*e*r*s+2*y1*(c*r+d*s)+a2)^(-(l)/2));
%                     df1(ii,jj,7)=df1(ii,jj,7)+wi(is)*F1(ii,jj,7)*dsodx;
%                     F2(ii,jj,7)= r^(i)*s^(j)*((r^2+s^2+y2^2+2*e*r*s+2*y2*(c*r+d*s)+a2)^(-(l)/2));
%                     df2(ii,jj,7)=df2(ii,jj,7)+wi(is)*F2(ii,jj,7)*dsodx;
%                 end
%             end
%         end
%     end
%     NIDy2(:,:,:)=NIDy2(:,:,:)+wi(ir)*df2(:,:,:)*drodx;
%     NIDy1(:,:,:)=NIDy1(:,:,:)+wi(ir)*df1(:,:,:)*drodx;
% end
% NID(:,:,:)=NIDy2(:,:,:)-NIDy1(:,:,:);
% TestD(:,1)= [SmartDot(vD00m3,signv) SmartDot(vD10m3,signv) SmartDot(vD01m3,signv) SmartDot(vD00m1,signv) SmartDot(vD10m1,signv)...
%     SmartDot(vD01m1,signv) SmartDot(vD20m3,signv) SmartDot(vD02m3,signv) SmartDot(vD11m3,signv) SmartDot(vD21m3,signv) SmartDot(vD12m3,signv)...
%     SmartDot(vD001,signv) SmartDot(vD003,signv) SmartDot(vD101,signv) SmartDot(vD011,signv)  SmartDot(vD11m1,signv) SmartDot(vD20m1,signv) SmartDot(vD02m1,signv)...
%     SmartDot(vD30m3,signv) SmartDot(vD03m3,signv) SmartDot(vD22m3,signv) SmartDot(vD31m3,signv)  SmartDot(vD13m3,signv)];
% TestD(:,2)= [NID(1,1,5) NID(2,1,5) NID(1,2,5) NID(1,1,3) NID(2,1,3) ...
%     NID(1,2,3) NID(3,1,5) NID(1,3,5) NID(2,2,5) NID(3,2,5) NID(2,3,5)...
%     NID(1,1,1) NID(1,1,7) NID(2,1,1) NID(1,2,1) NID(2,2,3) NID(3,1,3) NID(1,3,3)...
%     NID(4,1,5) NID(1,4,5) NID(3,3,5) NID(4,2,5) NID(2,4,5)];
% TestD(:,3)=(TestD(:,1)-TestD(:,2));
% 
% nt=size(TestD(:,1)); TEST_D(3)= 0;
% for i=1:nt
%     TEST_D(:)= TestD(i,:)
% end

%pause

Es1(1:2,:,:)=0;Es2(3:4,:,:)=0;Es1(5:6,:,:)=0;Es2(7:8,:,:)=0;
E03s1(1:2)=0;E03s2(3:4)=0;E03s1(5:6)=0;E03s2(7:8)=0;
E=Es1+Es2;
% vE00m3(1:8)=E(1:8,1,1,5);
% vE10m3(1:8)=E(1:8,2,1,5);
% vE01m3(1:8)=E(1:8,1,2,5);
% vE00m1(1:8)=E(1:8,1,1,3);
% vE10m1(1:8)=E(1:8,2,1,3);
% vE01m1(1:8)=E(1:8,1,2,3);
% vE20m3(1:8)=E(1:8,3,1,5);
% vE02m3(1:8)=E(1:8,1,3,5);
% vE11m3(1:8)=E(1:8,2,2,5);
% vE21m3(1:8)=E(1:8,3,2,5);
% vE12m3(1:8)=E(1:8,2,3,5);
% vE001(1:8)=E(1:8,1,1,1);
% vE003(1:8)=E03s2(1:8)+E03s1(1:8);
% vE101(1:8)=E(1:8,2,1,1);
% vE20m1(1:8)=E(1:8,3,1,3);
% vE11m1(1:8)=E(1:8,2,2,3);
% vE31m3(1:8)=E(1:8,4,2,5);
% vE22m3(1:8)=E(1:8,3,3,5);

% scE00m3=SmartDot(vE00m3, signv)
% scE10m3=SmartDot(vE10m3, signv)
% scE01m3=SmartDot(vE01m3, signv)
% scE00m1=SmartDot(vE00m1, signv)
% scE10m1=SmartDot(vE10m1, signv)
% scE01m1=SmartDot(vE01m1, signv)
% scE20m3=SmartDot(vE20m3, signv)
% scE02m3=SmartDot(vE02m3, signv)
% scE11m3=SmartDot(vE11m3, signv)
% scE21m3=SmartDot(vE21m3, signv)
% scE12m3=SmartDot(vE12m3, signv)
% scE001=SmartDot(vE001, signv)
% scE003=SmartDot(vE003, signv)
% scE101=SmartDot(vE101, signv)
% scE20m1=SmartDot(vE20m1, signv)
% scE11m1=SmartDot(vE11m1, signv)
% scE31m3=SmartDot(vE31m3, signv)
% scE22m3=SmartDot(vE22m3, signv)

% F1=@(r,y) 1.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) 1.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE00m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% F1=@(r,y) r.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) r.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE10m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% F1=@(r,y) y.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) y.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE01m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% 
% F12=@(r) log(2*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1);
% F11=@(r) log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1);
% F22=@(r) log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m));
% F21=@(r) log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m));
% NIE00m1=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F12=@(r) r.*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1);
% F11=@(r) r.*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1);
% F22=@(r) r.*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m));
% F21=@(r) r.*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m));
% NIE10m1=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F12=@(r) sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2) - (c*r+d*s1).*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1);
% F11=@(r) sqrt(r.^2+y1^2+s1^2+2*c*r*y1+2*s1*(d*y1+e*r)+a2) - (c*r+d*s1).*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1);
% F22=@(r) sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) - (c*r+d*(h*r+m)).*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m));
% F21=@(r) sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) - (c*r+d*(h*r+m)).*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m));
% NIE01m1=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F1=@(r,y) r.*r.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) r.*r.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE20m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% 
% 
% 
% F11=@(r) log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1)- 2*(4*(c*r +d*s1).^2*y1+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y1)./sqrt(r.^2+y1^2+s1^2+2*c*r*y1+2*s1*(d*y1+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F12=@(r) log(2*sqrt(r.^2+y2^2+s1^2+2*c*r.*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1)- 2*(4*(c*r +d*s1).^2*y2+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y2)./sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F22=@(r) log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r.*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m))- 2*(4*(c*r +d*(h*r+m)).^2*y2+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y2)./sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% F21=@(r) log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r.*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m))- 2*(4*(c*r +d*(h*r+m)).^2*y1+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y1)./sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% NIE02m3=(quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16));
% 
% F1=@(r,y) y.*r.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) y.*r.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE11m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% F1=@(r,y) y.*r.*r.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) y.*r.*r.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE21m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% 
% F11=@(r) r.*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1)- 2*r.*(4*(c*r +d*s1).^2*y1+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y1)./sqrt(r.^2+y1^2+s1^2+2*c*r*y1+2*s1*(d*y1+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F12=@(r) r.*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r.*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1)- 2*r.*(4*(c*r +d*s1).^2*y2+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y2)./sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F22=@(r) r.*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r.*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m))- 2*r.*(4*(c*r +d*(h*r+m)).^2*y2+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y2)./sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% F21=@(r) r.*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r.*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m))- 2*r.*(4*(c*r +d*(h*r+m)).^2*y1+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y1)./sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% NIE12m3=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% % bb= (2*c*r +2*d*s)
% % cc= (r.^2 +s.^2 +2*e*r.*s +a2)
% 
% F11=@(r) +1/4*((2*c*r +2*d*s1) +2*y1).*sqrt(r.^2+y1^2+s1.^2+2*c*r*y1+2*s1.*(d*y1+e*r)+a2) -1/8* ((2*c*r +2*d*s1).^2 -4*(r.^2 +s1.^2 +2*e*r.*s1 +a2)).*log(2*sqrt(r.^2+y1^2+s1.^2+2*c*r.*y1+2*s1.*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1) ;
% F12=@(r) +1/4*((2*c*r +2*d*s1) +2*y2).*sqrt(r.^2+y2^2+s1.^2+2*c*r*y2+2*s1.*(d*y2+e*r)+a2) -1/8* ((2*c*r +2*d*s1).^2 -4*(r.^2 +s1.^2 +2*e*r.*s1 +a2)).*log(2*sqrt(r.^2+y2^2+s1.^2+2*c*r.*y2+2*s1.*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1) ;
% F21=@(r) +1/4*((2*c*r +2*d*(h*r+m)) +2*y1).*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) -1/8* ((2*c*r +2*d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)).*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r.*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m)) ;
% F22=@(r) +1/4*((2*c*r +2*d*(h*r+m)) +2*y2).*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) -1/8* ((2*c*r +2*d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)).*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r.*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m)) ;
% NIE001=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% % F1=@(r,y) (r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(3/2);
% % F2=@(r,y) (r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(3/2);
% % NIE003=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% F11=@(r) r.^2.*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1)- 2*r.^2.*(4*(c*r +d*s1).^2*y1+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y1)./sqrt(r.^2+y1^2+s1^2+2*c*r*y1+2*s1*(d*y1+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F12=@(r) r.^2.*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r.*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1)- 2*r.^2.*(4*(c*r +d*s1).^2*y2+ 2*(c*r +d*s1).*(r.^2 +s1.^2 +2*e*r.*s1 +a2) -2*(r.^2 +s1.^2 +2*e*r*s1 +a2)*y2)./sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2) ./(4*(c*r+d*s1).^2 -4*(r.^2 +s1^2 +2*e*r.*s1 +a2) ) ;
% F22=@(r) r.^2.*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r.*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m))- 2*r.^2.*(4*(c*r +d*(h*r+m)).^2*y2+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y2)./sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% F21=@(r) r.^2.*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r.*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m))- 2*r.^2.*(4*(c*r +d*(h*r+m)).^2*y1+ 2*(c*r +d*(h*r+m)).*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) -2*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)*y1)./sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) ./(4*(c*r+d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2) ) ;
% NIE22m3=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F11=@(r) +1/4*r.*((2*c*r +2*d*s1) +2*y1).*sqrt(r.^2+y1^2+s1.^2+2*c*r*y1+2*s1.*(d*y1+e*r)+a2) -1/8* r.*((2*c*r +2*d*s1).^2 -4*(r.^2 +s1.^2 +2*e*r.*s1 +a2)).*log(2*sqrt(r.^2+y1^2+s1.^2+2*c*r.*y1+2*s1.*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1) ;
% F12=@(r) +1/4*r.*((2*c*r +2*d*s1) +2*y2).*sqrt(r.^2+y2^2+s1.^2+2*c*r*y2+2*s1.*(d*y2+e*r)+a2) -1/8* r.*((2*c*r +2*d*s1).^2 -4*(r.^2 +s1.^2 +2*e*r.*s1 +a2)).*log(2*sqrt(r.^2+y2^2+s1.^2+2*c*r.*y2+2*s1.*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1) ;
% F21=@(r) +1/4*r.*((2*c*r +2*d*(h*r+m)) +2*y1).*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) -1/8* r.*((2*c*r +2*d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)).*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r.*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m)) ;
% F22=@(r) +1/4*r.*((2*c*r +2*d*(h*r+m)) +2*y2).*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) -1/8* r.*((2*c*r +2*d*(h*r+m)).^2 -4*(r.^2 +(h*r+m).^2 +2*e*r.*(h*r+m) +a2)).*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r.*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m)) ;
% NIE101=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F12=@(r) r.^2.*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1);
% F11=@(r) r.^2.*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1);
% F22=@(r) r.^2.*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m));
% F21=@(r) r.^2.*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m));
% NIE20m1=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F12=@(r) r.*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2) - r.*(c*r+d*s1).*log(2*sqrt(r.^2+y2^2+s1^2+2*c*r*y2+2*s1*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*s1);
% F11=@(r) r.*sqrt(r.^2+y1^2+s1^2+2*c*r*y1+2*s1*(d*y1+e*r)+a2) - r.*(c*r+d*s1).*log(2*sqrt(r.^2+y1^2+s1^2+2*c*r.*y1+2*s1*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*s1);
% F22=@(r) r.*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2) - r.*(c*r+d*(h*r+m)).*log(2*sqrt(r.^2+y2^2+(h*r+m).^2+2*c*r*y2+2*(h*r+m).*(d*y2+e*r)+a2)+2*y2 +2*c*r+2*d*(h*r+m));
% F21=@(r) r.*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2) - r.*(c*r+d*(h*r+m)).*log(2*sqrt(r.^2+y1^2+(h*r+m).^2+2*c*r*y1+2*(h*r+m).*(d*y1+e*r)+a2)+2*y1 +2*c*r+2*d*(h*r+m));
% NIE11m1=quad(F22,r1,r2,1e-16)-quad(F12,r1,r2,1e-16)+quad(F11,r1,r2,1e-16)-quad(F21,r1,r2,1e-16);
% 
% F1=@(r,y) r.*r.*r.*y.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) r.*r.*r.*y.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIE31m3=dblquad(F2,r1,r2,y1,y2,1e-16)-dblquad(F1,r1,r2,y1,y2,1e-16);
% 
% TestE(:,1)= [vE00m3*signv' vE10m3*signv' vE01m3*signv' vE00m1*signv' vE10m1*signv' vE01m1*signv'  vE20m3*signv' vE02m3*signv' vE11m3*signv' vE21m3*signv' vE12m3*signv' vE101*signv' vE001*signv'   vE20m1*signv' vE11m1*signv' vE31m3*signv' vE22m3*signv']; % vE003*signv']%
% TestE(:,2)= [NIE00m3 NIE10m3 NIE01m3 NIE00m1 NIE10m1 NIE01m1  NIE20m3 NIE02m3 NIE11m3 NIE21m3 NIE12m3 NIE101 NIE001   NIE20m1 NIE11m1 NIE31m3 NIE22m3]; % NIE003]% TestE(:,2)= [ NIE00m3 NIE10m3 NIE01m3 NIE00m1  NIE10m1 NIE01m1  NIE20m3 NIE02m3 NIE11m3 NIE21m3 NIE12m3 NIE001 0 NIE101 NIE20m1 NIE11m1 NIE31m3 NIE22m3];
% TestE(:,3)=(TestE(:,1)-TestE(:,2));
% nt=size(TestE(:,1)); TEST_E(3)= 0;
% for i=1:nt
%     TEST_E(:)= TestE(i,:)
% end

%pause

% DEBUGGING F INTEGRALS
Fr1(1:4,:,:)=0;Fr2(5:8,:,:)=0;
F03r1(1:4)=0; F03r2(5:8)=0;
F=Fr1+Fr2;
% vF00m3(1:8)=F(1:8,1,1,5);
% vF10m3(1:8)=F(1:8,2,1,5);
% vF01m3(1:8)=F(1:8,1,2,5);
% vF00m1(1:8)=F(1:8,1,1,3);
% vF10m1(1:8)=F(1:8,2,1,3);
% vF01m1(1:8)=F(1:8,1,2,3);
% vF20m3(1:8)=F(1:8,3,1,5);
% vF02m3(1:8)=F(1:8,1,3,5);
% vF11m3(1:8)=F(1:8,2,2,5);
% vF21m3(1:8)=F(1:8,3,2,5);
% vF12m3(1:8)=F(1:8,2,3,5);
% vF001(1:8)=F(1:8,1,1,1);
% vF003(1:8)=F03r2(1:8)+F03r1(1:8);
% vF101(1:8)=F(1:8,2,1,1);
% vF20m1(1:8)=F(1:8,3,1,3);
% vF11m1(1:8)=F(1:8,2,2,3);
% vF31m3(1:8)=F(1:8,4,2,5);
% vF22m3(1:8)=F(1:8,3,3,5);

% scF00m3=SmartDot(vF00m3, signv)
% scF10m3=SmartDot(vF10m3, signv)
% scF01m3=SmartDot(vF01m3, signv)
% scF00m1=SmartDot(vF00m1, signv)
% scF10m1=SmartDot(vF10m1, signv)
% scF01m1=SmartDot(vF01m1, signv)
% scF20m3=SmartDot(vF20m3, signv)
% scF02m3=SmartDot(vF02m3, signv)
% scF11m3=SmartDot(vF11m3, signv)
% scF21m3=SmartDot(vF21m3, signv)
% scF12m3=SmartDot(vF12m3, signv)
% scF001=SmartDot(vF001, signv)
% scF003=SmartDot(vF003, signv)
% scF101=SmartDot(vF101, signv)
% scF20m1=SmartDot(vF20m1, signv)
% scF11m1=SmartDot(vF11m1, signv)
% scF31m3=SmartDot(vF31m3, signv)
% scF22m3=SmartDot(vF22m3, signv)

% F1=@(s,y) 1./(s.^2+r1^2+y.^2+2*c*r1*y+2*s*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) 1./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF00m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16)
% F1=@(s,y) s./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) s./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF10m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF01m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) 1./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(1/2);
% F2=@(s,y) 1./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(1/2);
% NIF00m1=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(1/2);
% F2=@(s,y) s./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(1/2);
% NIF10m1=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(1/2);
% F2=@(s,y) y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(1/2);
% NIF01m1=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*s./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) s.*s./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF20m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) y.*y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) y.*y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF02m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) s.*y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF11m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*s.*y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) s.*s.*y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF21m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*y.*y./(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) s.*y.*y./(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF12m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(1/2);
% F2=@(s,y) s.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(1/2);
% NIF101=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) (s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(1/2);
% F2=@(s,y) (s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(1/2);
% NIF001=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*s.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s*(d*y+e*r1)+a2).^(-1/2);
% F2=@(s,y) s.*s.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(-1/2);
% NIF20m1=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.*y.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(-1/2);
% F2=@(s,y) s.*y.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(-1/2);
% NIF11m1=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.^3.*y.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(-3/2);
% F2=@(s,y) s.^3.*y.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(-3/2);
% NIF31m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) s.^2.*y.^2.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(-3/2);
% F2=@(s,y) s.^2.*y.^2.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(-3/2);
% NIF22m3=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% F1=@(s,y) (s.^2+r1^2+y.^2+2*c*r1*y+2*s.*(d*y+e*r1)+a2).^(3/2);
% F2=@(s,y) (s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g).*y+2*s.*(d*y+e*(f*s+g))+a2).^(3/2);
% NIF003=dblquad(F2,s1,s2,y1,y2,1e-16)-dblquad(F1,s1,s2,y1,y2,1e-16);
% 
% TestF(:,1)= [SmartDot(vF00m3,signv)  SmartDot(vF10m3,signv) SmartDot(vF01m3,signv)  SmartDot(vF00m1,signv)  SmartDot(vF10m1,signv) SmartDot(vF01m1,signv) SmartDot(vF20m3,signv)  SmartDot(vF02m3,signv) SmartDot(vF11m3,signv) SmartDot(vF21m3,signv) SmartDot(vF12m3,signv) SmartDot(vF101,signv) SmartDot(vF001,signv) SmartDot(vF20m1,signv) SmartDot(vF11m1,signv) SmartDot(vF31m3,signv) SmartDot(vF22m3,signv)]; % vF003*signv']
% TestF(:,2)= [NIF00m3  NIF10m3  NIF01m3  NIF00m1  NIF10m1 NIF01m1 NIF20m3 NIF02m3 NIF11m3 NIF21m3 NIF12m3 NIF101 NIF001 NIF20m1 NIF11m1 NIF31m3 NIF22m3 ]; %NIF003];
% TestF(:,3)=(TestF(:,1)-TestF(:,2));
% nt=size(TestF(:,1)); TEST_F(3)= 0;
% for i=1:nt
%     TEST_F(:)= TestF(i,:)
% end

%pause


%TRIPLE INTEGRALS
% SEED INTEGRAL H[000-3]
% Dr1s1=@(y) 2./sqrt((1-e2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*y.^2).*atan(((1-e)*(sqrt(y.^2+r1^2+s1^2+2*y*(c*r1+d*s1)+2*e*r1*s1+a2)-s1-d*y-e*r1)+(c-d*e)*y+(1-e2)*r1)./sqrt((1-e2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*y.^2));
% Dr1s2=@(y) 2./sqrt((1-e2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*y.^2).*atan(((1-e)*(sqrt(y.^2+r1^2+s2^2+2*y*(c*r1+d*s2)+2*e*r1*s2+a2)-s2-d*y-e*r1)+(c-d*e)*y+(1-e2)*r1)./sqrt((1-e2)*a2+((1-d2)*(1-e2)-(c-d*e)^2)*y.^2));

% tp1=sqrt(1+f2+2*ef);
% tp3=sqrt(tp1-f-e);
% Dr2s1=@(y) 2/tp3./(sqrt(-((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3).^2-(tp1+f+e)*((y*(cf+d)+fg+eg)/tp1).^2+(tp1+f+e)*(a2+g2+y.^2+2*cg*y))).*atan((tp3*(sqrt(y.^2+(1+f2+2*e*f)*s1^2+2*s1*y*(cf+d)+2*s1*(fg+eg)+2*cg*y+g2+a2)-tp1*s1-((y*(cf+d)+fg+eg)/tp1))+((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3))./(sqrt(-((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3).^2-(tp1+f+e)*((y*(cf+d)+fg+eg)/tp1).^2+(tp1+f+e)*(a2+g2+y.^2+2*cg*y))));
% Dr2s2=@(y) 2/tp3./(sqrt(-((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3).^2-(tp1+f+e)*((y*(cf+d)+fg+eg)/tp1).^2+(tp1+f+e)*(a2+g2+y.^2+2*cg*y))).*atan((tp3*(sqrt(y.^2+(1+f2+2*e*f)*s2^2+2*s2*y*(cf+d)+2*s2*(fg+eg)+2*cg*y+g2+a2)-tp1*s2-((y*(cf+d)+fg+eg)/tp1))+((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3))./(sqrt(-((tp1*(c*y+g)-((y*(cf+d)+fg+eg)/tp1)*(f+e))/tp3).^2-(tp1+f+e)*((y*(cf+d)+fg+eg)/tp1).^2+(tp1+f+e)*(a2+g2+y.^2+2*cg*y))));

% nbin=10000;
% [xi,wi] = gqwp(nbin);
% AA=(r2-r1)/2;
% BB=(r2+r1)/2;
% CC=(s2-s1)/2;
% DD=(s2+s1)/2;
% drodx=AA;
% dsodx=CC;
% NID2=0;
% NID1=0;
% for i=1:nbin
%     r=AA*xi(i)+BB;
%     sf=-lqlp*r+s0;
%     df2=0;
%     df1=0;
%     for j=1:nbin
%         s=CC*xi(j)+DD;
%         if (s <= sf)
%             F1=-1./sqrt(y1^2+r.^2+s.^2+2*y1*(c*r+d*s)+2*e*r.*s+a2) ./(sqrt(y1^2+r.^2+s.^2+2*y1*(c*r+d*s)+2*e*r.*s+a2)+y1+c*r+d*s);
%             df1=df1+wi(j)*F1*dsodx;
%             F2=-1./sqrt(y2^2+r.^2+s.^2+2*y2*(c*r+d*s)+2*e*r.*s+a2) ./(sqrt(y2^2+r.^2+s.^2+2*y2*(c*r+d*s)+2*e*r.*s+a2)+y2+c*r+d*s);
%             df2=df2+wi(j)*F2*dsodx;
%         end
%     end
%     NID1=NID1+wi(i)*df1*drodx;
%     NID2=NID2+wi(i)*df2*drodx;
% end
% NIH00m3=NID2-NID1

%%%%%%%
H(1:8,1,1,1,5)=0;
% H(1,1,1,1,5)=quad(Dr2s2,y1,y2,1E-16);
% H(3,1,1,1,5)=quad(Dr2s1,y1,y2,1E-16);
% H(5,1,1,1,5)=quad(Dr1s2,y1,y2,1E-16);
% H(7,1,1,1,5)=quad(Dr1s1,y1,y2,1E-16);

for ll=1:5
    for kk=1:4
        for ii=1:4
            for jj=1:4
                if ii==1
                    Fri(1:8,ii,jj,kk,ll)=F(1:8,jj,kk,ll);
                elseif ii==2
                    Fri(1:8,ii,jj,kk,ll)=f*Fr2(:,jj+1,kk,ll)+g*Fr2(:,jj,kk,ll)+r1*Fr1(:,jj,kk,ll);
                elseif ii==3
                    Fri(1:8,ii,jj,kk,ll)=f2*Fr2(:,jj+2,kk,ll)+2*fg*Fr2(:,jj+1,kk,ll)+g2*Fr2(:,jj,kk,ll)+r1^2*Fr1(:,jj,kk,ll);
                elseif ii==4
                    Fri(1:8,ii,jj,kk,ll)=f3*Fr2(:,jj+3,kk,ll)+3*f2g*Fr2(:,jj+2,kk,ll)+3*fg2*Fr2(:,jj+1,kk,ll)+g3*Fr2(:,jj,kk,ll)+r1^3*Fr1(:,jj,kk,ll);
                end
            end
        end
    end
end
for ll=1:5
    for kk=1:4
        for jj=1:4
            for ii=1:4
                if jj==1
                    Esj(1:8,ii,jj,kk,ll)=E(1:8,ii,kk,ll);
                elseif jj==2
                    Esj(1:8,ii,jj,kk,ll)=h*Es2(:,ii+1,kk,ll)+m*Es2(:,ii,kk,ll)+s1*Es1(:,ii,kk,ll);
                elseif jj==3
                    Esj(1:8,ii,jj,kk,ll)=h2*Es2(:,ii+2,kk,ll)+2*hm*Es2(:,ii+1,kk,ll)+m2*Es2(:,ii,kk,ll)+s1^2*Es1(:,ii,kk,ll);
                elseif jj==4
                    Esj(1:8,ii,jj,kk,ll)=h3*Es2(:,ii+3,kk,ll)+3*h2m*Es2(:,ii+2,kk,ll)+3*hm2*Es2(:,ii+1,kk,ll)+m3*Es2(:,ii,kk,ll)+s1^3*Es1(:,ii,kk,ll);
                end
            end
        end
    end
end

% temp(1:8)=Fri(1:8,ii,jj,kk,ll);
% testFri=temp*signv'
% temp(1:8)=Esj(1:8,ii,jj,kk,ll);
% testEsj=temp*signv'
% F1=@(r,y) r.^i.*s1^j.*y.^k.*(r.^2+y.^2+s1^2+2*c*r.*y+2*s1*(d*y+e*r)+a2).^(-3/2);
% F2=@(r,y) r.^i.*(h*r+m).^j.*y.^k.*(r.^2+y.^2+(h*r+m).^2+2*c*r.*y+2*(h*r+m).*(d*y+e*r)+a2).^(-3/2);
% NIEsj=dblquad(F2,r1,r2,y1,y2)-dblquad(F1,r1,r2,y1,y2)
% F1=@(s,y) r1^i*s.^j.*y.^k.*(s.^2+r1^2+y.^2+2*c*r1*y+2*s*(d*y+e*r1)+a2).^(-3/2);
% F2=@(s,y) (f*s+g).^i.*s.^j.*y.^k.*(s.^2+(f*s+g).^2+y.^2+2*c*(f*s+g)*y+2*s.*(d*y+e*(f*s+g))+a2).^(-3/2);
% NIFri=dblquad(F2,s1,s2,y1,y2)-dblquad(F1,s1,s2,y1,y2)

% H[000-1]
i=0 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk,ll)=1/(l-i-j-k-3)*((l*a2)*H(:,ii,jj,kk,ll+2)-Fri(:,ii+1,jj,kk,ll)-Esj(:,ii,jj+1,kk,ll)-yv.^(k+1).*D(:,ii,jj,ll));

% H[0001] QUAD
i=0 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk,ll)=1/(l-i-j-k-3)*((l*a2)*H(:,ii,jj,kk,ll+2)-Fri(:,ii+1,jj,kk,ll)-Esj(:,ii,jj+1,kk,ll)-yv.^(k+1).*D(:,ii,jj,ll));

% H[001-1] H[100-1] H[010-1]
i=0 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*F(:,jj,kk,ll)-(c*e-d)*E(:,ii,kk,ll));
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[101-1] H[110-1] H[011-1] QUAD
i=1 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))-i*(c-de)*H(:,ii-1,jj,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll)-(ce-d)*Esj(:,ii,jj,kk,ll));
i=0 ; j=1 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))-j*(d-ce)*H(:,ii,jj-1,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll)-(ce-d)*Esj(:,ii,jj,kk,ll));

% !!!!!!!!
% ERREUR TROUVEE !!!! >>>>>>H101m1 et H011m1 re?cris de mani?re identique
% !!!!!!!!

%H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-e*H(:,ii,jj-1,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H110m1_1=1/l/(1-e2)*(-e*H(:,ii,jj-1,kk,ll) -Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=1 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-e*H(:,ii-1,jj,kk,ll) +e*Fri(:,ii,jj,kk,ll)-Esj(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
H110m1_2=1/l/(1-e2)*(-e*H(:,ii-1,jj,kk,ll) +e*Fri(:,ii,jj,kk,ll)-Esj(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H(:,2,2,1,3)= 0.5*H110m1_1 +0.5*H110m1_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!!!!!!!
% ERREUR TROUVEE !!!! >>>>>>H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
% !!!!!!!!

% H[200-1] H[200-1] QUAD
i=1 ; j=0 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=1 ; k=0 ; l=-1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[100-3] H[010-3] H[001-3]
i=0 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))+(c-de)*F(:,jj,kk,ll)-(ce-d)*E(:,ii,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*(smartsum([(1-e2)*(-yv.^k.*D(:,ii,jj,ll)) (c-de)*F(:,jj,kk,ll) -(ce-d)*E(:,ii,kk,ll)]));

% !!!!!!!!
% ERREUR TROUVEE !!!! round off H001m3 -> smartsum
% !!!!!!!!

H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[101-3] H[011-3] H[110-3]
i=1 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))+(c-de)*Fri(:,ii,jj,kk,ll)-i*(c-de)*H(:,ii-1,jj,kk,ll)-(ce-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

H110m3_1=H(:,ii,jj+1,kk,ll+2);

i=0 ; j=1 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))-(ce-d)*Esj(:,ii,jj,kk,ll) -j*(d-ce)*H(:,ii,jj-1,kk,ll) +(c-de)*Fri(:,ii,jj,kk,ll));
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));

Fri(:,ii,jj,kk,ll) = F(:,jj,kk,ll);
Esj(:,ii,jj,kk,ll) = h*Es2(:,ii+1,kk,ll) +m*Es2(:,ii,kk,ll) +s1*Es1(:,ii,kk,ll);

H110m3_2=H(:,ii+1,jj,kk,ll+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H(:,2,2,kk,ll+2) =  0.5*H110m3_1 + 0.5*H110m3_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H[002-3] H[200-3] H[020-3]
i=0 ; j=0 ; k=1 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
i=1 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=1 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[201-3] H[021-3] QUAD
i=2 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))-i*(c-de)*H(:,ii-1,jj,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll)-(ce-d)*Esj(:,ii,jj,kk,ll));
i=0 ; j=2 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))-j*(d-ce)*H(:,ii,jj-1,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll)-(ce-d)*Esj(:,ii,jj,kk,ll));


% H[111-3] H[210-3] H[120-3] QUAD
i=1 ; j=1 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*(SmartSum([(1-e2)*(-yv.^k.*D(:,ii,jj,ll))  -i*(c-de)*H(:,ii-1,jj,kk,ll)  +(c-de)*Fri(:,ii,jj,kk,ll)  -j*(d-ce)*H(:,ii,jj-1,kk,ll)  -(ce-d)*Esj(:,ii,jj,kk,ll)]));

H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));

clear test
test(:,1)= i*H(:,ii-1,jj,kk,ll);
test(:,2)= -Fri(:,ii,jj,kk,ll);
test(:,3)= -e*j*H(:,ii,jj-1,kk,ll);
test(:,4)= +e*Esj(:,ii,jj,kk,ll);
test(:,5)= -(c-de)*l*H(:,ii,jj,kk+1,ll+2);
vtest(1:8) = 1/l/(1-e2)*SmartSum(test);
H(:,ii+1,jj,kk,ll+2) = vtest(1:8);

H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

clear test
test(:,1)= j*H(:,ii,jj-1,kk,ll);
test(:,2)= -Esj(:,ii,jj,kk,ll);
test(:,3)= -e*i*H(:,ii-1,jj,kk,ll);
test(:,4)= +e*Fri(:,ii,jj,kk,ll);
test(:,5)= -(d-ce)*l*H(:,ii,jj,kk+1,ll+2);
vtest(1:8) = 1/l/(1-e2)*SmartSum(test);
H(:,ii,jj+1,kk,ll+2) = vtest(1:8);

%i=0 ; j=2 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll
%)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));

i=0 ; j=1 ; k=1 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll)) +(c-d*e)*Fri(:,ii,jj,kk,ll) -(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
i=1 ; j=0 ; k=1 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
% !!!!!!!!
% ERREUR TROUVEE !!!! >>>>>>H0123 missing
% !!!!!!!!
% !!!!!!!!
% ERREUR TROUVEE !!!! >>>>>>H1023 missing
% !!!!!!!!

% H[300-3] H[030-3] QUAD
i=2 ; j=0 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=2 ; k=0 ; l=1 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[001-5]
i=0 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*(SmartSum([(1-e2)*(-yv.^k.*D(:,ii,jj,ll))  +(c-de)*Fri(:,ii,jj,kk,ll)  -(ce-d)*Esj(:,ii,jj,kk,ll)]));
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[101-5] H[101-5] H[110-5]
i=1 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
i=0 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
i=0 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
scH110m5_1=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=1 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
scH110m5_2=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H(:,2,2,kk,ll+2) = 0.5*scH110m5_1 + 0.5*scH110m5_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H[002-5] H[200-5] H[020-5]
i=0 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
i=1 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[111-5]
i=1 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));

% H[210-5] H[120-5] H[102-5] H[012-5]
i=1 ; j=1 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

i=1 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-i*(c-de)*H(:,ii-1,jj,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll));
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-j*(d-ce)*H(:,ii,jj-1,kk,ll)+(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[202-5] H[022-5]
i=2 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*(smartsum([(1-e2)*(k*H(:,ii,jj,kk-1,ll) -yv.^k.*D(:,ii,jj,ll)) -i*(c-de)*H(:,ii-1,jj,kk,ll) +(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll)]));
i=0 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*(smartsum([(1-e2)*(k*H(:,ii,jj,kk-1,ll) -yv.^k.*D(:,ii,jj,ll)) -j*(d-ce)*H(:,ii,jj-1,kk,ll) +(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll)]));

% H[112-5] H[211-5] H[121-5]
i=1 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll) ...
    -(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));

clear test
test(:,1)= (1-e2)*(k*H(:,ii,jj,kk-1,ll));
test(:,2)= -(1-e2)*yv.^k.*D(:,ii,jj,ll);
test(:,3)= -(c*i-d*e*i)*H(:,ii-1,jj,kk,ll);
test(:,4)= +(c-d*e)*Fri(:,ii,jj,kk,ll);
test(:,5)= -(d*j-e*c*j)*H(:,ii,jj-1,kk,ll);
test(:,6)= -(c*e-d)*Esj(:,ii,jj,kk,ll);
vtest(1:8) = (1/l/(1-e2-c2-d2+2*e*c*d))*SmartSum(test);
H(:,ii,jj,kk+1,ll+2) = vtest(1:8);

H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[301-5] H[031-5]
i=2 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[300-5] H[030-5] QUAD
i=2 ; j=0 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=2 ; k=0 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% H[302-5] H[032-5] QUAD
i=3 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll) -yv.^k.*D(:,ii,jj,ll)) -i*(c-de)*H(:,ii-1,jj,kk,ll) +(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll));
i=0 ; j=3 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
%H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll) -yv.^k.*D(:,ii,jj,ll)) -j*(d-ce)*H(:,ii,jj-1,kk,ll) +(c-de)*Fri(:,ii,jj,kk,ll) -(ce-d)*Esj(:,ii,jj,kk,ll));

%H[311-5] H[131-5] QUAD
i=3 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=3 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=1 ; j=1 ; k=2 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));

% H[212-5] H[122-5] H[221-5] QUAD
i=2 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll) ...
-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));

clear test
test(:,1)= (1-e2)*(k*H(:,ii,jj,kk-1,ll));
test(:,2)= -(1-e2)*yv.^k.*D(:,ii,jj,ll);
test(:,3)= -(c*i-d*e*i)*H(:,ii-1,jj,kk,ll);
test(:,4)= +(c-d*e)*Fri(:,ii,jj,kk,ll);
test(:,5)= -(d*j-e*c*j)*H(:,ii,jj-1,kk,ll);
test(:,6)= -(c*e-d)*Esj(:,ii,jj,kk,ll);
vtest(1:8) = (1/l/(1-e2-c2-d2+2*e*c*d))*SmartSum(test);
H(:,ii,jj,kk+1,ll+2) = vtest(1:8);

i=1 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));

clear test
test(:,1)= (1-e2)*(k*H(:,ii,jj,kk-1,ll));
test(:,2)= -(1-e2)*yv.^k.*D(:,ii,jj,ll);
test(:,3)= -(c*i-d*e*i)*H(:,ii-1,jj,kk,ll);
test(:,4)= +(c-d*e)*Fri(:,ii,jj,kk,ll);
test(:,5)= -(d*j-e*c*j)*H(:,ii,jj-1,kk,ll);
test(:,6)= -(c*e-d)*Esj(:,ii,jj,kk,ll);
vtest(1:8) = (1/l/(1-e2-c2-d2+2*e*c*d))*SmartSum(test);
H(:,ii,jj,kk+1,ll+2) = vtest(1:8);

% i=1 ; j=1 ; k=2 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
% H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-
% e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
i=1 ; j=2 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));

clear test
test(:,1)= i*H(:,ii-1,jj,kk,ll);
test(:,2)= -Fri(:,ii,jj,kk,ll);
test(:,3)= -e*j*H(:,ii,jj-1,kk,ll);
test(:,4)= e*Esj(:,ii,jj,kk,ll);
test(:,5)= -(c-de)*l*H(:,ii,jj,kk+1,ll+2);
vtest(1:8) = 1/l/(1-e2)*SmartSum(test);
H221m5_1 = vtest(1:8);

i=2 ; j=1 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

clear test
test(:,1)= j*H(:,ii,jj-1,kk,ll);
test(:,2)= -Esj(:,ii,jj,kk,ll);
test(:,3)= -e*i*H(:,ii-1,jj,kk,ll);
test(:,4)= +e*Fri(:,ii,jj,kk,ll);
test(:,5)= -(d-ce)*l*H(:,ii,jj,kk+1,ll+2);
vtest(1:8) = 1/l/(1-e2)*SmartSum(test);
H221m5_2 = vtest(1:8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H(:,ii,jj+1,kk,ll+2)= 0.5*H221m5_1 + 0.5*H221m5_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H[401-5] H[041-5] QUAD
i=3 ; j=0 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
i=0 ; j=3 ; k=1 ; l=3 ; ii=i+1;jj=j+1;kk=k+1;ll=l+2;
H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));

% DEBUGGING TRIPPLE INTEGRALS
% nbin=100;
% nbiny=20
% [xi,wi] = gqwp(nbin);
% [xy,wy] = gqwp(nbiny);
% AA=(r2-r1)/2;
% BB=(r2+r1)/2;
% CC=(s2-s1)/2;
% DD=(s2+s1)/2;
% EE=(y2-y1)/2;
% FF=(y2+y1)/2;
% drodx=AA;
% dsodx=CC;
% dyodx=EE;
% NIH(1:6,1:6,1:6,1:7)=0;
% for iy=1:nbiny
%     iy
%     y=EE*xy(iy)+FF;
%     dFrs(1:6,1:6,1:6,1:7)=0;
%     for ir=1:nbin
%         r=AA*xi(ir)+BB;
%         sf=-lqlp*r+s0;
%         dFs(1:6,1:6,1:6,1:7)=0;
%         for is=1:nbin
%             s=CC*xi(is)+DD;
%             if (s <= sf)
%                 l=3; ijkmax=3; imax=3; jmax=3;
%                 for (i=0:ijkmax)
%                     for   j=0:ijkmax-i
%                         for   k=0:ijkmax-i-j
%                             ii=i+1;jj=j+1;kk=k+1;
%                             F1= r^(i)*s^(j)*y^(k)*(r^2+s^2+y^2+2*e*r*s+2*y*(c*r+d*s)+a2)^(-(l)/2);
%                             dFs(ii,jj,kk,l+2)=dFs(ii,jj,kk,l+2)+wi(is)*F1*dsodx;
%                         end
%                     end
%                 end
%                 l=5; ijkmax=5; imax=4; jmax=4;
%                 for (i=0:ijkmax)
%                     if(i < imax+1)
%                         for   (j=0:ijkmax-i)
%                             if (j < jmax+1)
%                                 for   k=0:ijkmax-i-j
%                                     ii=i+1;jj=j+1;kk=k+1;
%                                     F1= r^(i)*s^(j)*y^(k)*((r^2+s^2+y^2+2*e*r*s+2*y*(c*r+d*s)+a2)^(-(l)/2));
%                                     dFs(ii,jj,kk,l+2)=dFs(ii,jj,kk,l+2)+wi(is)*F1*dsodx;
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 
%             end
%         end
%         dFrs(:,:,:,:)=dFrs(:,:,:,:)+wi(ir)*dFs(:,:,:,:)*drodx;
%     end
%     NIH(:,:,:,:)=NIH(:,:,:,:)+wy(iy)*dFrs(:,:,:,:)*dyodx;
% end
% testH3(1:16,1)=[SmartDot(H(1:8,1,1,1,5),signv(:)) SmartDot(H(:,2,1,1,5),signv(:)) SmartDot(H(:,1,2,1,5),signv(:)) SmartDot(H(:,1,1,2,5),signv(:))...
%     SmartDot(H(:,2,2,1,5),signv(:)) SmartDot(H(:,1,2,2,5),signv(:)) SmartDot(H(:,2,1,2,5),signv(:))...
%     SmartDot(H(:,2,2,2,5),signv(:)) SmartDot(H(:,3,2,1,5),signv(:)) SmartDot(H(:,1,3,2,5),signv(:)) SmartDot(H(:,3,1,2,5),signv(:))...
%     SmartDot(H(:,2,3,1,5),signv(:)) SmartDot(H(:,1,2,3,5),signv(:)) SmartDot(H(:,2,1,3,5),signv(:))...
%     SmartDot(H(:,4,1,1,5),signv(:)) SmartDot(H(:,1,4,1,5),signv(:))];
% testH3(:,2)=[NIH(1,1,1,5) NIH(2,1,1,5) NIH(1,2,1,5) NIH(1,1,2,5)...
%     NIH(2,2,1,5) NIH(1,2,2,5) NIH(2,1,2,5)...
%     NIH(2,2,2,5) NIH(3,2,1,5) NIH(1,3,2,5) NIH(3,1,2,5)...
%     NIH(2,3,1,5) NIH(1,2,3,5) NIH(2,1,3,5)...
%     NIH(4,1,1,5) NIH(1,4,1,5)]
% testH5(1:4,1)=[SmartDot(H(:,1,1,1,7),signv(:)) SmartDot(H(:,2,1,1,7),signv(:)) SmartDot(H(:,1,2,1,7),signv(:)) SmartDot(H(:,1,1,2,7),signv(:))];
% testH5(1:4,2)=[NIH(1,1,1,7) NIH(2,1,1,7) NIH(1,2,1,7) NIH(1,1,2,7)]
% 
% testH5(1:3,1)=[ SmartDot(H(:,2,2,1,7),signv(:)) SmartDot(H(:,1,2,2,7),signv(:)) SmartDot(H(:,2,1,2,7),signv(:))];
% testH5(1:3,2)=[ NIH(2,2,1,7) NIH(1,2,2,7) NIH(2,1,2,7) ]
% 
% testH5(1:9,1)=[SmartDot(H(:,2,2,2,7),signv(:)) SmartDot(H(:,3,2,1,7),signv(:)) SmartDot(H(:,1,3,2,7),signv(:)) SmartDot(H(:,3,1,2,7),signv(:))...
%     SmartDot(H(:,2,3,1,7),signv(:)) SmartDot(H(:,1,2,3,7),signv(:)) SmartDot(H(:,2,1,3,7),signv(:))...
%     SmartDot(H(:,4,1,1,7),signv(:)) SmartDot(H(:,1,4,1,7),signv(:))];
% testH5(1:9,2)=[ NIH(2,2,2,7) NIH(3,2,1,7) NIH(1,3,2,7) NIH(3,1,2,7)...
%     NIH(2,3,1,7) NIH(1,2,3,7) NIH(2,1,3,7)...
%     NIH(4,1,1,7) NIH(1,4,1,7)]
% 
% testH5(1:9,1)=[SmartDot(H(:,3,3,1,7),signv(:)) SmartDot(H(:,3,1,3,7),signv(:)) SmartDot(H(:,1,3,3,7),signv(:))...
%     SmartDot(H(:,3,2,2,7),signv(:)) SmartDot(H(:,2,3,2,7),signv(:)) SmartDot(H(:,3,2,2,7),signv(:))...
%     SmartDot(H(:,2,3,2,7),signv(:)) SmartDot(H(:,2,2,3,7),signv(:)) SmartDot(H(:,2,2,3,7),signv(:))];
% testH5(1:9,2)=  [NIH(3,3,1,7) NIH(3,1,3,7) NIH(1,3,3,7)...
%     NIH(3,2,2,7) NIH(2,3,2,7) NIH(3,2,2,7)...
%     NIH(2,3,2,7) NIH(2,2,3,7) NIH(2,2,3,7)]
% 
% testH5(1:11,1)=[SmartDot(H(:,4,2,2,7),signv(:)) SmartDot(H(:,2,4,2,7),signv(:))...
%     SmartDot(H(:,4,1,3,7),signv(:)) SmartDot(H(:,1,4,3,7),signv(:)) SmartDot(H(:,4,3,1,7),signv(:)) SmartDot(H(:,3,4,1,7),signv(:))...
%     SmartDot(H(:,3,3,2,7),signv(:)) SmartDot(H(:,3,2,3,7),signv(:)) SmartDot(H(:,2,3,3,7),signv(:))...
%     SmartDot(H(:,5,1,2,7),signv(:)) SmartDot(H(:,1,5,2,7),signv(:))];
% testH5(1:11,2)=[ NIH(4,2,2,7) NIH(2,4,2,7)...
%     NIH(4,1,3,7) NIH(1,4,3,7) NIH(4,3,1,7) NIH(3,4,1,7)...
%     NIH(3,3,2,7) NIH(3,2,3,7) NIH(2,3,3,7)...
%     NIH(5,1,2,7) NIH(1,5,2,7)]

% nt=size(testH3(:,1)) 
% TEST_H3(1:3)= 0;
% for i=1:nt
%     TEST_H3(:)= testH3(i,:)
% end
% nt=size(testH5(:,1)) 
% TEST_H5(1:3)= 0;
% for i=1:nt
%     TEST_H5(:)= testH5(i,:)
% end

%pause

%##########################################################################
%

% %INTEGRALES DOUBLES
% %Integrales E sont fonctions seulement de s,y
% %Integrales F sont fonctions seulement de r,y
% %rien de change pour le moment pour les integrales D
% %Indices kk+2 can be shifter to kk as long as k is shifted as well
% D(:,ii+1,jj,ll+2)=1/(1-e2)*(1/l*(i*D(:,ii-1,jj,ll)-e*j*D(:,ii,jj-1,ll)+e*sv.^j.*A(:,ii,ll)-rv.^i.*B(:,jj,ll))-(c-d*e)*yv.*D(:,ii,jj,ll+2))
% D(:,ii,jj+1,ll+2)=1/(1-e2)*(1/l*(i*D(:,ii,jj-1,ll)-e*i*D(:,ii-1,jj,ll)+e*rv.^i.*B(:,jj,ll)-sv.^j.*A(:,ii,ll))-(d-c*e)*yv.*D(:,ii,jj,ll+2))
% D(:,ii,jj,ll+2)=1./(l*(a2+yv2).*(1-e2)+2*c*d*e*l*yv2-l*yv2*(c2+d2)).*(((1-e2)*(l-2-i-j)*D(:,ii,jj,ll)+i*c*(e-1)*yv.*D(:,ii-1,jj,ll)+j*d*yv.*D(:,ii,jj-1,ll)+(yv.*sv.^j*(d-e*c)+sv.^(j+1)*(1-e2)).*A(:,ii,ll)+(yv.*rv.^i*(c-e*d)+rv.^(i+1)*(1-e2)).*B(:,jj,ll)))
%
% D(:,ii,jj,ll) =1/((1-e2)*(l-2-i-j))*((l*(a2+yv2).*(1-e2)+2*c*d*e*l*yv2-l*yv2*(c2+d2)).*D(:,ii,jj,ll+2)-i*c*(e-1)*yv.*D(:,ii-1,jj,ll)-j*d*yv.*D(:,ii,jj-1,ll)-(yv.*sv.^j*(d-e*c)+sv.^(j+1)*(1-e2)).*A(:,ii,ll)-(yv.*rv.^i*(c-e*d)+rv.^(i+1)*(1-e2)).*B(:,jj,ll))
%
% E(:,ii+1,kk,ll+2)=1/(1-c2)*(1/l*(i*E(:,ii-1,kk,ll)-c*k*E(:,ii,kk-1,ll)+c*yv.^k*A(:,ii,ll)-rfs.^i.*sC(:,kk,ll))+(c*d-e)*sv.*E(:,ii,kk,ll+2))
% E(:,ii,kk+1,ll+2)=1/(1-c2)*(1/l*(k*E(:,ii,kk-1,ll)-c*i*E(:,ii-1,kk,ll)+c*rfs.^i*sC(:,kk,ll)-yv.^k.*A(:,ii,ll))+(c*e-d)*sv.*E(:,ii,kk,ll+2))
% E(:,ii,kk,ll+2)=1./(l*(a2+sv2)*(1-c2)+2*c*d*e*l*sv2-l*sv2*(e2+d2)).*(((1-c2)*(l-2-i-k)*E(:,ii,kk,ll)+i*(c*d-e)*sv.*E(:,ii-1,kk,ll)-k*(d-e*c)*sv.*E(:,ii,kk-1,ll)+rfs^i*(rfs*(1-c2)-d*sv*c+e*sv).*sC(:,kk,ll)+yv.^k.*(yv*(1-c2)+d*sv-e*c*sv).*A(:,ii,ll)))
%
% E(:,ii,kk,ll) = 1/(1-c2)/(l-2-i-k)*((l*(a2+sv2)*(1-c2)+2*c*d*e*l*sv2-l*sv2*(e2+d2)).*E(:,ii,kk,ll+2)+(-i*(c*d-e)*sv.*E(:,ii-1,kk,ll)+k*(d-e*c)*sv.*E(:,ii,kk-1,ll)-rfs^i*(rfs*(1-c2)-d*sv*c+e*sv).*sC(:,kk,ll)-yv.^k.*(yv*(1-c2)+d*sv-e*c*sv).*A(:,ii,ll)))
%
% Fr1(:,jj+1,kk,ll)=1/(1-d2)*(1/(l-2)*(j*Fr1(:,jj-1,kk,ll-2)-d*k*Fr1(:,jj,kk-1,ll-2)+d*yv.^k.*Br1(:,jj,ll-2)-sv.^j.*Cr1(:,kk,ll-2))+(c*d-e)*rv.*Fr1(:,jj,kk,ll));
% Fr2(:,jj+1,kk,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*yv.^k.*Br2(:,jj,ll-2)-sv.^j.*Cr2(:,kk,ll-2)+j*Fr2(:,jj-1,kk,ll-2)-k*(cf+d)*Fr2(:,jj,kk-1,ll-2)-(l-2)*(fg+eg-cg*(cf+d))*Fr2(:,jj,kk,ll) );
% Fr1(:,jj,kk+1,ll)=1/(1-d2)*(1/(l-2)*(k*Fr1(:,jj,kk-1,ll-2)-d*j*Fr1(:,jj-1,kk,ll-2)+d*sv.^j.*Cr1(:,kk,ll-2)-yv.^k.*Br1(:,jj,ll-2))+(d*e-c)*rv.*Fr1(:,jj,kk,ll));
% Fr2(:,jj,kk+1,ll)=1/(l-2)/(1+f2+2*ef-(cf+d)^2)*((cf+d)*sv.^j.*Cr2(:,kk,ll-2)-(1+f2+2*ef)*yv.^k.*Br2(:,jj,ll-2)+k*(1+f2+2*ef)*Fr2(:,jj,kk-1,ll-2)-j*(cf+d)*Fr2(:,jj-1,kk,ll-2)+(l-2)*((cf+d)*(fg+eg)-(1+f2+2*ef)*cg)*Fr2(:,jj,kk,ll));
%
% F(:,jj,kk,ll+2)=1./(l*(a2+rv2)*(1-d2)+2*c*d*e*l*rv2-l*rv2*(e2+c2)).*(((1-d2)*(l-2-j-k)*F(:,jj,kk,ll)+j*(c*d-e)*rv.*F(:,jj-1,kk,ll)-k*(c-e*d)*rv.*F(:,jj,kk-1,ll)+sfr.^j*(sfr*(1-d2)-d*rv*c+e*rv).*rC(:,kk,ll)+yv.^k*(yv*(1-d2)+c*rv-e*d*rv).*B(:,jj,ll)))
%
% %INTEGRALES TRIPPLES
%
% H(:,ii,jj,kk+1,ll+2)=(1/l/(1-e2-c2-d2+2*e*c*d))*((1-e2)*(k*H(:,ii,jj,kk-1,ll)-yv.^k.*D(:,ii,jj,ll))-(c*i-d*e*i)*H(:,ii-1,jj,kk,ll)+(c-d*e)*Fri(:,ii,jj,kk,ll)-(d*j-e*c*j)*H(:,ii,jj-1,kk,ll)-(c*e-d)*Esj(:,ii,jj,kk,ll));
% H(:,ii+1,jj,kk,ll+2)=1/l/(1-e2)*(i*H(:,ii-1,jj,kk,ll)-Fri(:,ii,jj,kk,ll)-e*j*H(:,ii,jj-1,kk,ll)+e*Esj(:,ii,jj,kk,ll)-(c-de)*l*H(:,ii,jj,kk+1,ll+2));
% H(:,ii,jj+1,kk,ll+2)=1/l/(1-e2)*(j*H(:,ii,jj-1,kk,ll)-Esj(:,ii,jj,kk,ll)-e*i*H(:,ii-1,jj,kk,ll)+e*Fri(:,ii,jj,kk,ll)-(d-ce)*l*H(:,ii,jj,kk+1,ll+2));
%
% H(:,ii,jj,kk,ll+2)=(1/(l*a2))*((l-i-j-k-3)*H(:,ii,jj,kk,ll)+rv.^(i+1).*F(:,jj,kk,ll)+sv.^(j+1).*E(:,ii,kk,ll)+yv.^(k+1).*D(:,ii,jj,ll));
%##########################################################################

% FINAL INTEGRAL EVALUATION
scH000m3= SmartDot(H(1:8,1,1,1,5), signv(1:8)');
scH000m1= SmartDot(H(1:8,1,1,1,3), signv(1:8)');
scH0001= SmartDot(H(1:8,1,1,1,1), signv(1:8)');

scH001m1= SmartDot(H(1:8,1,1,2,3), signv(1:8)');
scH100m1= SmartDot(H(1:8,2,1,1,3), signv(1:8)');
scH010m1= SmartDot(H(1:8,1,2,1,3), signv(1:8)');

scH200m1= SmartDot(H(1:8,3,1,1,3), signv(1:8)');
scH020m1= SmartDot(H(1:8,1,3,1,3), signv(1:8)');

scH101m1= SmartDot(H(1:8,2,1,2,3), signv(1:8)');
scH011m1= SmartDot(H(1:8,1,2,2,3), signv(1:8)');
scH110m1= SmartDot(H(1:8,2,2,1,3), signv(1:8)');

scH001m3= SmartDot(H(1:8,1,1,2,5), signv(1:8)');
scH100m3= SmartDot(H(1:8,2,1,1,5), signv(1:8)');
scH010m3= SmartDot(H(1:8,1,2,1,5), signv(1:8)');
%scH002m3= H(1:8,1,1,3,5), signv(1:8)');
scH200m3= SmartDot(H(1:8,3,1,1,5), signv(1:8)');
scH020m3= SmartDot(H(1:8,1,3,1,5), signv(1:8)');
scH101m3= SmartDot(H(1:8,2,1,2,5), signv(1:8)');
scH011m3= SmartDot(H(1:8,1,2,2,5), signv(1:8)');
scH110m3= SmartDot(H(1:8,2,2,1,5), signv(1:8)');

scH201m3= SmartDot(H(1:8,3,1,2,5), signv(1:8)') ; % QUAD
scH021m3= SmartDot(H(1:8,1,3,2,5), signv(1:8)') ; % QUAD
scH210m3= SmartDot(H(1:8,3,2,1,5), signv(1:8)') ; % QUAD
scH120m3= SmartDot(H(1:8,2,3,1,5), signv(1:8)') ; %; QUAD
scH111m3= SmartDot(H(1:8,2,2,2,5), signv(1:8)') ; % QUAD
scH300m3= SmartDot(H(1:8,4,1,1,5), signv(1:8)') ; % QUAD
scH030m3= SmartDot(H(1:8,1,4,1,5), signv(1:8)') ; % QUAD

scH001m5= SmartDot(H(1:8,1,1,2,7), signv(1:8)') ;
scH100m5= SmartDot(H(1:8,2,1,1,7), signv(1:8)') ;
scH010m5= SmartDot(H(1:8,1,2,1,7), signv(1:8)'); 
%scH002m5= H(1:8,1,1,3,7), signv(1:8)') 
scH200m5= SmartDot(H(1:8,3,1,1,7), signv(1:8)') ;
scH020m5= SmartDot(H(1:8,1,3,1,7), signv(1:8)') ;
scH101m5= SmartDot(H(1:8,2,1,2,7), signv(1:8)') ;
scH011m5= SmartDot(H(1:8,1,2,2,7), signv(1:8)') ;
scH110m5= SmartDot(H(1:8,2,2,1,7), signv(1:8)') ;
scH111m5= SmartDot(H(1:8,2,2,2,7), signv(1:8)') ;
scH112m5= SmartDot(H(1:8,2,2,3,7), signv(1:8)') ;
scH211m5= SmartDot(H(1:8,3,2,2,7), signv(1:8)') ;
scH121m5= SmartDot(H(1:8,2,3,2,7), signv(1:8)') ;
scH102m5= SmartDot(H(1:8,2,1,3,7), signv(1:8)') ;
scH012m5= SmartDot(H(1:8,1,2,3,7), signv(1:8)') ;
scH201m5= SmartDot(H(1:8,3,1,2,7), signv(1:8)') ;
scH021m5= SmartDot(H(1:8,1,3,2,7), signv(1:8)') ;
scH202m5= SmartDot(H(1:8,3,1,3,7), signv(1:8)') ;
scH022m5= SmartDot(H(1:8,1,3,3,7), signv(1:8)') ;
scH301m5= SmartDot(H(1:8,4,1,2,7), signv(1:8)') ;
scH031m5= SmartDot(H(1:8,1,4,2,7), signv(1:8)') ;

scH210m5= SmartDot(H(1:8,3,2,1,7), signv(1:8)') ;% QUAD
scH120m5= SmartDot(H(1:8,2,3,1,7), signv(1:8)') ;% QUAD
scH300m5= SmartDot(H(1:8,4,1,1,7), signv(1:8)') ;% QUAD
scH030m5= SmartDot(H(1:8,1,4,1,7), signv(1:8)') ;% QUAD
scH221m5= SmartDot(H(1:8,3,3,2,7), signv(1:8)') ;% QUAD
scH212m5= SmartDot(H(1:8,3,2,3,7), signv(1:8)') ;% QUAD
scH122m5= SmartDot(H(1:8,2,3,3,7), signv(1:8)') ;% QUAD
scH401m5= SmartDot(H(1:8,5,1,2,7), signv(1:8)') ;% QUAD
scH041m5= SmartDot(H(1:8,1,5,2,7), signv(1:8)') ;% QUAD
scH302m5= SmartDot(H(1:8,4,1,3,7), signv(1:8)') ;% QUAD
scH032m5= SmartDot(H(1:8,1,4,3,7), signv(1:8)') ;% QUAD
scH311m5= SmartDot(H(1:8,4,2,2,7), signv(1:8)') ;% QUAD
scH131m5= SmartDot(H(1:8,2,4,2,7), signv(1:8)') ;% QUAD

% FINAL SUM
% the force is now calculated for each of the 3 patch points
tcrossb= cross(t,b);
pcrossb= cross(p,b);
qcrossb= cross(q,b);
bcrosst= cross(b,t);
tdotn= dot(t,n);
ttbn= t*dot(tcrossb,n);
tbtn= t*dot(bcrosst,n);
tpbn= t*dot(pcrossb,n);
tqbn= t*dot(qcrossb,n);
factor=mu*alpha/4/pi/(1-nu);

I001m3=(tcrossb*tdotn+t*(dot(tcrossb,n)))*(1-nu) + (bcrosst*tdotn+tbtn);
I010m3=(qcrossb*tdotn+tqbn)*(1-nu)-(dot(qcrossb,t)*n)+(q*(dot(bcrosst,n)));
I100m3=(pcrossb*tdotn+tpbn)*(1-nu)-(dot(pcrossb,t)*n)+(p*(dot(bcrosst,n)));
I001m5=(tcrossb*tdotn+ttbn)*1.5*(1-nu)*a2;
I010m5=(qcrossb*tdotn+tqbn)*1.5*(1-nu)*a2-(dot(qcrossb,t)*n) * 3*a2;
I100m5=(pcrossb*tdotn+tpbn)*1.5*(1-nu)*a2-(dot(pcrossb,t)*n) * 3*a2;
I021m5=-(dot(qcrossb,t)*(tdotn))*q * 3;
I201m5=-(dot(pcrossb,t)*(tdotn))*p * 3;
I012m5=-(dot(qcrossb,t)*tdotn)*t * 3;
I102m5=-(dot(pcrossb,t)*tdotn)*t * 3;
I111m5=-((dot(pcrossb,t)*tdotn)*q +(dot(qcrossb,t)*tdotn)*p)*3;

aa=2;
bb=-(3*r1+r2);
cc=r1*r1+r1*r2;

F001m3_x4=aa*scH201m3+bb*scH101m3+cc*scH001m3 ;
F010m3_x4=aa*scH210m3+bb*scH110m3+cc*scH010m3 ;
F100m3_x4=aa*scH300m3+bb*scH200m3+cc*scH100m3 ;

F001m5_x4=aa*scH201m5+bb*scH101m5+cc*scH001m5 ;
F010m5_x4=aa*scH210m5+bb*scH110m5+cc*scH010m5 ;
F100m5_x4=aa*scH300m5+bb*scH200m5+cc*scH100m5 ;
F021m5_x4=aa*scH221m5+bb*scH121m5+cc*scH021m5 ;
F201m5_x4=aa*scH401m5+bb*scH301m5+cc*scH201m5 ;
F012m5_x4=aa*scH212m5+bb*scH112m5+cc*scH012m5 ;
F102m5_x4=aa*scH302m5+bb*scH202m5+cc*scH102m5 ;
F111m5_x4=aa*scH311m5+bb*scH211m5+cc*scH111m5 ;

fLLprime= I001m3 * F001m3_x4 + I100m3 * F100m3_x4 + I010m3 * F010m3_x4 ...
    + I001m5 * F001m5_x4 + I100m5 * F100m5_x4 + I010m5 * F010m5_x4 ...
    + I102m5 * F102m5_x4 + I012m5 * F012m5_x4 + I111m5 * F111m5_x4 ...
    + I201m5 * F201m5_x4 + I021m5 * F021m5_x4;

fx4= fLLprime*factor/Lp/Lp;

aa=2;
bb=-(3*s1+s2);
cc=s1*s1+s1*s2;

F001m3_x5=aa*scH021m3+bb*scH011m3+cc*scH001m3;
F010m3_x5=aa*scH030m3+bb*scH020m3+cc*scH010m3; 
F100m3_x5=aa*scH120m3+bb*scH110m3+cc*scH100m3; 

F001m5_x5=aa*scH021m5+bb*scH011m5+cc*scH001m5; 
F010m5_x5=aa*scH030m5+bb*scH020m5+cc*scH010m5; 
F100m5_x5=aa*scH120m5+bb*scH110m5+cc*scH100m5; 
F021m5_x5=aa*scH041m5+bb*scH031m5+cc*scH021m5; 
F201m5_x5=aa*scH221m5+bb*scH211m5+cc*scH201m5; 
F012m5_x5=aa*scH032m5+bb*scH022m5+cc*scH012m5; 
F102m5_x5=aa*scH122m5+bb*scH112m5+cc*scH102m5; 
F111m5_x5=aa*scH131m5+bb*scH121m5+cc*scH111m5; 

factor/Lp/Lp*I001m3*F001m3_x5;
factor/Lp/Lp*I100m3*F100m3_x5; 
factor/Lp/Lp*I010m3*F010m3_x5;

factor/Lp/Lp*I001m5*F001m5_x5;
factor/Lp/Lp*I100m5*F100m5_x5;
factor/Lp/Lp*I010m5*F010m5_x5;

factor/Lp/Lp*I102m5*F102m5_x5;
factor/Lp/Lp*I012m5 *F012m5_x5;
factor/Lp/Lp*I111m5 *F111m5_x5;

factor/Lp/Lp*I201m5 *F201m5_x5;
factor/Lp/Lp*I021m5 *F021m5_x5;
    

fLLprime= I001m3 * F001m3_x5 + I100m3 * F100m3_x5 + I010m3 * F010m3_x5...
    + I001m5 * F001m5_x5 + I100m5 * F100m5_x5 + I010m5 * F010m5_x5 ...
    + I102m5 * F102m5_x5 + I012m5 * F012m5_x5 + I111m5 * F111m5_x5 ...
    + I201m5 * F201m5_x5 + I021m5 * F021m5_x5;

fx5= fLLprime*factor/Lq/Lq;

aa=2/Lp/Lp;
bb=2/Lq/Lq;
cc=4/Lp/Lq;
dd=-4*r1/Lp/Lp-3/Lp-4*s1/Lp/Lq;
ee=-4*s1/Lq/Lq-3/Lq-4*r1/Lp/Lq;
ff=1+4*r1*s1/Lp/Lq+3*r1/Lp+3*s1/Lq+2*r1*r1/Lp/Lp+2*s1*s1/Lq/Lq;


F001m3_x3=aa*scH201m3 +dd*scH101m3 +bb*scH021m3 +ee*scH011m3 +cc*scH111m3 +ff*scH001m3;
F010m3_x3=aa*scH210m3 +dd*scH110m3 +bb*scH030m3 +ee*scH020m3 +cc*scH120m3 +ff*scH010m3;
F100m3_x3=aa*scH300m3 +dd*scH200m3 +bb*scH120m3 +ee*scH110m3 +cc*scH210m3 +ff*scH100m3;
 
F001m5_x3=aa*scH201m5 +dd*scH101m5 +bb*scH021m5 +ee*scH011m5 +cc*scH111m5 +ff*scH001m5;
F010m5_x3=aa*scH210m5 +dd*scH110m5 +bb*scH030m5 +ee*scH020m5 +cc*scH120m5 +ff*scH010m5; 
F100m5_x3=aa*scH300m5 +dd*scH200m5 +bb*scH120m5 +ee*scH110m5 +cc*scH210m5 +ff*scH100m5;

F021m5_x3=aa*scH221m5 +dd*scH121m5 +bb*scH041m5 +ee*scH031m5 +cc*scH131m5 +ff*scH021m5;
F201m5_x3=aa*scH401m5 +dd*scH301m5 +bb*scH221m5 +ee*scH211m5 +cc*scH311m5 +ff*scH201m5;
F012m5_x3=aa*scH212m5 +dd*scH112m5 +bb*scH032m5 +ee*scH022m5 +cc*scH122m5 +ff*scH012m5;
F102m5_x3=aa*scH302m5 +dd*scH202m5 +bb*scH122m5 +ee*scH112m5 +cc*scH212m5 +ff*scH102m5; 
F111m5_x3=aa*scH311m5 +dd*scH211m5 +bb*scH131m5 +ee*scH121m5 +cc*scH221m5 +ff*scH111m5;

clear test
test=[I001m3' * F001m3_x3   I100m3' * F100m3_x3   I010m3' * F010m3_x3 ...
      I001m5' * F001m5_x3   I100m5' * F100m5_x3   I010m5' * F010m5_x3 ...
      I102m5' * F102m5_x3   I012m5' * F012m5_x3   I111m5' * F111m5_x3 ...
      I201m5' * F201m5_x3   I021m5' * F021m5_x3];

fLLprime= smartsum(test)';

fx3= fLLprime*factor;

aa=-4/Lp/Lp;
bb=-4/Lp/Lq;
cc=4/Lp+8*r1/Lp/Lp+4*s1/Lp/Lq;
dd=4*r1/Lp/Lq;
ee=-4*r1/Lp-4*r1*r1/Lp/Lp-4*r1*s1/Lp/Lq;

F001m3_x6=aa*scH201m3+bb*scH111m3+cc*scH101m3+dd*scH011m3+ee*scH001m3; 
F010m3_x6=aa*scH210m3+bb*scH120m3+cc*scH110m3+dd*scH020m3+ee*scH010m3; 
F100m3_x6=aa*scH300m3+bb*scH210m3+cc*scH200m3+dd*scH110m3+ee*scH100m3; 

F001m5_x6=aa*scH201m5+bb*scH111m5+cc*scH101m5+dd*scH011m5+ee*scH001m5; 
F010m5_x6=aa*scH210m5+bb*scH120m5+cc*scH110m5+dd*scH020m5+ee*scH010m5; 
F100m5_x6=aa*scH300m5+bb*scH210m5+cc*scH200m5+dd*scH110m5+ee*scH100m5; 

F021m5_x6=aa*scH221m5+bb*scH131m5+cc*scH121m5+dd*scH031m5+ee*scH021m5; 
F201m5_x6=aa*scH401m5+bb*scH311m5+cc*scH301m5+dd*scH211m5+ee*scH201m5; 
F012m5_x6=aa*scH212m5+bb*scH122m5+cc*scH112m5+dd*scH022m5+ee*scH012m5; 
F102m5_x6=aa*scH302m5+bb*scH212m5+cc*scH202m5+dd*scH112m5+ee*scH102m5; 
F111m5_x6=aa*scH311m5+bb*scH221m5+cc*scH211m5+dd*scH121m5+ee*scH111m5; 

fLLprime= I001m3 * F001m3_x6 + I100m3 * F100m3_x6 + I010m3 * F010m3_x6 ...
    + I001m5 * F001m5_x6 + I100m5 * F100m5_x6 + I010m5 * F010m5_x6 ...
    + I102m5 * F102m5_x6 + I012m5 * F012m5_x6 + I111m5 * F111m5_x6 ...
    + I201m5 * F201m5_x6 + I021m5 * F021m5_x6;

fx6= fLLprime*factor;

aa=-4/Lq/Lq;
bb=-4/Lp/Lq;
cc=4*s1/Lp/Lq; 
dd=4/Lq+8*s1/Lq/Lq+4*r1/Lp/Lq;
ee=-4*s1/Lq-4*s1*s1/Lq/Lq-4*r1*s1/Lp/Lq;

F001m3_x7=aa*scH021m3+bb*scH111m3+cc*scH101m3+dd*scH011m3+ee*scH001m3; 
F010m3_x7=aa*scH030m3+bb*scH120m3+cc*scH110m3+dd*scH020m3+ee*scH010m3; 
F100m3_x7=aa*scH120m3+bb*scH210m3+cc*scH200m3+dd*scH110m3+ee*scH100m3; 

F001m5_x7=aa*scH021m5+bb*scH111m5+cc*scH101m5+dd*scH011m5+ee*scH001m5; 
F010m5_x7=aa*scH030m5+bb*scH120m5+cc*scH110m5+dd*scH020m5+ee*scH010m5; 
F100m5_x7=aa*scH120m5+bb*scH210m5+cc*scH200m5+dd*scH110m5+ee*scH100m5; 
F021m5_x7=aa*scH041m5+bb*scH131m5+cc*scH121m5+dd*scH031m5+ee*scH021m5; 
F201m5_x7=aa*scH221m5+bb*scH311m5+cc*scH301m5+dd*scH211m5+ee*scH201m5; 
F012m5_x7=aa*scH032m5+bb*scH122m5+cc*scH112m5+dd*scH022m5+ee*scH012m5; 
F102m5_x7=aa*scH122m5+bb*scH212m5+cc*scH202m5+dd*scH112m5+ee*scH102m5; 
F111m5_x7=aa*scH131m5+bb*scH221m5+cc*scH211m5+dd*scH121m5+ee*scH111m5; 

fLLprime= I001m3 * F001m3_x7 + I100m3 * F100m3_x7 + I010m3 * F010m3_x7 ...
    + I001m5 * F001m5_x7 + I100m5 * F100m5_x7 + I010m5 * F010m5_x7 ...
    + I102m5 * F102m5_x7 + I012m5 * F012m5_x7 + I111m5 * F111m5_x7 ...
    + I201m5 * F201m5_x7 + I021m5 * F021m5_x7;

fx7= fLLprime*factor;

bb=4;
cc=-4*s1;
dd=-4*r1;
ee=4*r1*s1;

F001m3_x8=bb*scH111m3+cc*scH101m3+dd*scH011m3+ee*scH001m3; 
F010m3_x8=bb*scH120m3+cc*scH110m3+dd*scH020m3+ee*scH010m3; 
F100m3_x8=bb*scH210m3+cc*scH200m3+dd*scH110m3+ee*scH100m3; 

F001m5_x8=bb*scH111m5+cc*scH101m5+dd*scH011m5+ee*scH001m5; 
F010m5_x8=bb*scH120m5+cc*scH110m5+dd*scH020m5+ee*scH010m5; 
F100m5_x8=bb*scH210m5+cc*scH200m5+dd*scH110m5+ee*scH100m5; 
F021m5_x8=bb*scH131m5+cc*scH121m5+dd*scH031m5+ee*scH021m5; 
F201m5_x8=bb*scH311m5+cc*scH301m5+dd*scH211m5+ee*scH201m5; 
F012m5_x8=bb*scH122m5+cc*scH112m5+dd*scH022m5+ee*scH012m5; 
F102m5_x8=bb*scH212m5+cc*scH202m5+dd*scH112m5+ee*scH102m5; 
F111m5_x8=bb*scH221m5+cc*scH211m5+dd*scH121m5+ee*scH111m5; 

fLLprime= I001m3 * F001m3_x8 + I100m3 * F100m3_x8 + I010m3 * F010m3_x8 ...
    + I001m5 * F001m5_x8 + I100m5 * F100m5_x8 + I010m5 * F010m5_x8 ...
    + I102m5 * F102m5_x8 + I012m5 * F012m5_x8 + I111m5 * F111m5_x8 ...
    + I201m5 * F201m5_x8 + I021m5 * F021m5_x8;

fx8= fLLprime*factor/Lp/Lq;

ftot=fx3+fx4+fx5+fx6+fx7+fx8;

end

function sum = SmartDot(u,v)

n = length(u);
w = u .* v;
w = bubblesortabsval(w);

sum=0;
for i = 1:n;
    sum=sum+w(i);
end

end

function sum = SmartSum(u)

ni = size(u(:,1),1);
nj = size(u(1,:),2);

sum(1:ni)= 0;
usort(1:ni,1:nj)= 0;
vecu(1:nj)= 0;

for i = 1:ni;
    vecu(:)= u(i,:);
    usort(i,:)= bubblesortabsval(vecu);
    
    for j =1:nj;
        sum(i)= sum(i) +usort(i,j);
    end
end

sum= sum';

end

function x = bubblesortabsval(x)
%--------------------------------------------------------------------------
% based on bubble sort by: Brian Moore
%               brimoor@umich.edu
%--------------------------------------------------------------------------

% Bubble sort
n = length(x);
while (n > 0)
    % Iterate through x
    nnew = 0;
    for i = 2:n
        % Swap elements in wrong order
        if (abs(x(i)) < abs(x(i - 1)))
            x = swap(x,i,i - 1);
            nnew = i;
        end
    end
    n = nnew;
end

end

function x = swap(x,i,j)
% Swap x(i) and x(j)
% Note: In practice, x xhould be passed by reference

val = x(i);
x(i) = x(j);
x(j) = val;

end
