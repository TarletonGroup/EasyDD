function [kee,kec,kcc]=intNs_gauss(segnode0,segnode1,elenode1,s0,s1,wx,wy,wz)

gauss=5;
%gauss points
if gauss == 2
    ss=[-1 1]/sqrt(3);
    wt=[1 1];
   
elseif gauss == 3
    ss = [-sqrt(3/5) 0 sqrt(3/5)];
    wt = [5 8 5]/9;
  
elseif gauss == 4
    ss = [   -sqrt((15+2*sqrt(30))/35),  -sqrt((15-2*sqrt(30))/35), ...
        sqrt((15-2*sqrt(30))/35),   sqrt((15+2*sqrt(30))/35)];
    wt = [  (90-5*sqrt(30))/180,    (90+5*sqrt(30))/180,...
        (90+5*sqrt(30))/180,    (90-5*sqrt(30))/180];
elseif gauss == 5
    ss = [ -sqrt(245-14*sqrt(70))/21,  -sqrt(245+14*sqrt(70))/21, 0, ...
        sqrt(245+14*sqrt(70))/21,   sqrt(245-14*sqrt(70))/21];
    wt = [  (322+13*sqrt(70))/900,    (322-13*sqrt(70))/900, 128/225,...
        (322-13*sqrt(70))/900,    (322+13*sqrt(70))/900];
end

s=(s1-s0)/2*ss +(s1+s0)/2;
J0=(s1-s0)/2;
w=[wx wy wz];
pm1 =[-1  1  1 -1 -1  1  1 -1];
pm2 =[-1 -1  1  1 -1 -1  1  1];
pm3 =[-1 -1 -1 -1  1  1  1  1];
%
Nes=zeros(1,8);
kee=zeros(8,8);
kec=zeros(8,2);
kcc=zeros(2,2);
 for mm=1:length(s)
     Ncs=[(1-s(mm))/2 (1+s(mm))/2];
     Xs=Ncs*[segnode0;segnode1];
     us=-1+2*(Xs-elenode1)./w;

     for i=1:8
         Nes(i)=1/8*(1+pm1(i)*us(1))*(1+pm2(i)*us(2))*(1+pm3(i)*us(3));
     end

     kes=Nes'*Nes;
     kecs=-Nes'*Ncs;
     kcs=Ncs'*Ncs;    
     
     kee= kee + wt(mm)*kes*J0;
     kec= kec + wt(mm)*kecs*J0;
     kcc= kcc + wt(mm)*kcs*J0;
 end


end
