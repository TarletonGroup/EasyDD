
function [finc0,finc1]=inc_nodalforce_gauss(MU,NU,segments)
% add stress due to inclusion 2019/11/12 fengxian Liu
% for reference: Xiang Y, Srolovitz D J. Philosophical magazine, 2006. (https://doi.org/10.1080/14786430600575427)

%integral along the segment, gauss inegration 2020/09/16 Fengxian Liu

global RR QQ strain0

%RR=[30 ]/bur; % radius of the spheres
%QQ=[50 0 0; ]/bur; %location of the centers of inclusion


LINKMAX=size(segments,1);
finc=zeros(1,3);
finc0=zeros(LINKMAX,3);
finc1=zeros(LINKMAX,3);

% syms s; %[-1,1]

% strain0=0.001; % dialatational misfit strain;
NU_inc = NU;
MU_inc = MU*1;
strain1=strain0*(1+NU_inc)/(1-NU_inc);


gauss=3;
%gauss points
if gauss == 2
    ss=[-1 1]/sqrt(3);
    wt=[1 1];
    %         fprintf('***************************  2     points gauss  *******')
elseif gauss == 3
    ss = [-sqrt(3/5) 0 sqrt(3/5)];
    wt = [5 8 5]/9;
    %         fprintf('***************************    3   points gauss  *******')
elseif gauss == 4
    %         fprintf('***************************    4   points gauss  *******')
    ss = [   -sqrt((15+2*sqrt(30))/35),  -sqrt((15-2*sqrt(30))/35), ...
        sqrt((15-2*sqrt(30))/35),   sqrt((15+2*sqrt(30))/35)];
    wt = [  (90-5*sqrt(30))/180,    (90+5*sqrt(30))/180,...
        (90+5*sqrt(30))/180,    (90-5*sqrt(30))/180];
elseif gauss == 5
    %         fprintf('***************************    5    points gauss *******')
    %         ss=zeros(1,gauss);
    %         ss(1)=.906179845938664 ; ss(2)=.538469310105683;
    %         ss(3)=.0;      ss(4)=-s(2) ; ss(5)=-s(1);
    ss = [ -sqrt(245-14*sqrt(70))/21,  -sqrt(245+14*sqrt(70))/21, 0, ...
        sqrt(245+14*sqrt(70))/21,   sqrt(245-14*sqrt(70))/21];
    wt = [  (322+13*sqrt(70))/900,    (322-13*sqrt(70))/900, 128/225,...
        (322-13*sqrt(70))/900,    (322+13*sqrt(70))/900];
    %         wt(1)=.236926885056189 ; wt(2)=.478628670499366;
    %         wt(3)=.568888888888889 ; wt(4)=wt(2) ; wt(5)=wt(1);
end

for m=1:length(ss)
    s=ss(m);
    for i=1:LINKMAX
        
        r0=segments(i,6:8);
        r1=segments(i,9:11);
        
        r0r1=r1-r0;
        magr0r1=norm(r0r1);
        r0r1_nomalize=r0r1/magr0r1;
        XX = r0*(1+s)/2 + r1*(1-s)/2;
        
        b=segments(i,3:5);
        Sigma=zeros(3,3);
        sigma_pt=zeros(3,3);
        
        %    rr=zeros(size(RR,2) ,3); % distance from XX to QQ
        
        for n=1:size(RR,2)
            
            rr(n,1:3)=XX-QQ(n,1:3);
            magrr=norm(rr);
            
            Sigma=[magrr^2-3*(rr(n,1))^2   -3*rr(n,1)*rr(n,2)           -3*rr(n,1)*rr(n,3);
                -3*rr(n,1)*rr(n,2)         (magrr)^2-3*(rr(n,2))^2    -3*rr(n,2)*rr(n,3);
                -3*rr(n,1)*rr(n,3)          -3*rr(n,2)*rr(n,3)          (magrr)^2-3*(rr(n,3))^2;];
            
            
            deltaV(n)=2*MU_inc*strain1*RR(n)^3;
            
            sigma_pt= sigma_pt + Sigma*deltaV(n)/magrr^5;
        end
        
        sigma = sigma_pt*b';
        finc= cross(sigma',r0r1_nomalize);
        
        J=1/2*magr0r1; % dx=L/2*ds
        
        fun0=finc*(1-s)/2*J;
        fun1=finc*(1+s)/2*J;
        
        finc0(i,1:3)= finc0(i,1:3)+ wt(m)*fun0;
        finc1(i,1:3)= finc1(i,1:3)+ wt(m)*fun1;
        
        
    end
end
end

