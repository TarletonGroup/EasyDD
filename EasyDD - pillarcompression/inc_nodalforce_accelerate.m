
function [finc0,finc1]=inc_nodalforce_accelerate(MU,NU,segments,linkid)
% add stress due to inclusion 2019/11/12 fengxian Liu
% for reference: Xiang Y, Srolovitz D J. Philosophical magazine, 2006. (https://doi.org/10.1080/14786430600575427)

global RR QQ strain0

LINKMAX=size(segments,1);
ninc=size(RR,2);

if linkid==0
    finc0=zeros(LINKMAX,3);
    finc1=finc0;
    QQ2r0=zeros(ninc,3);
    QQ2r0mag=zeros(ninc,1);
    QQ2r1=QQ2r0;
    QQ2r1mag=QQ2r0mag;
    
    Sigincr0=zeros(ninc,6);
    Sigincr1=Sigincr0;
    
    stressr0=zeros(3,3);
    stressr1=stressr0;
    
    deltaV=QQ2r0mag;
    
%     strain0=0.01; % dialatational misfit strain;
    NU_inc = NU;
    MU_inc = MU*1;
    strain1=strain0*(1+NU_inc)/(1-NU_inc);
    
    
    for i=1:LINKMAX
        r0=segments(i,6:8);
        r1=segments(i,9:11);
        
        rf0=((r0+r1)/2+r0)/2;
        rf1=((r0+r1)/2+r1)/2;
        
        rmid_r0=(r0+r1)/2-r0;
        r1_rmid=r1-(r0+r1)/2;
        b=segments(i,3:5);
        
        QQ2r0=rf0-QQ(:,1:3);
        QQ2r0mag=vecnorm(QQ2r0');
        QQ2r1=rf1-QQ(:,1:3);
        QQ2r1mag=vecnorm(QQ2r1');
        deltaV=2*MU_inc*strain1*RR.^3;
        deltaV=deltaV./(QQ2r0mag.^5);
        
        
        Sigincr0(:,1) = QQ2r0mag.^2' - 3*QQ2r0(:,1).^2; % sigma11
        Sigincr0(:,2) = QQ2r0mag.^2' - 3*QQ2r0(:,2).^2; %22
        Sigincr0(:,3) = QQ2r0mag.^2' - 3*QQ2r0(:,3).^2; %33
        Sigincr0(:,4) = -3*QQ2r0(:,1).*QQ2r0(:,2); %12
        Sigincr0(:,5) = -3*QQ2r0(:,1).*QQ2r0(:,3); %13
        Sigincr0(:,6) = -3*QQ2r0(:,2).*QQ2r0(:,3); %23
        
        Sigincr0=Sigincr0.*deltaV';
        
        Sigincr1(:,1) = QQ2r1mag.^2' - 3*QQ2r1(:,1).^2; % sigma11
        Sigincr1(:,2) = QQ2r1mag.^2' - 3*QQ2r1(:,2).^2; %22
        Sigincr1(:,3) = QQ2r1mag.^2' - 3*QQ2r1(:,3).^2; %33
        Sigincr1(:,4) = -3*QQ2r1(:,1).*QQ2r1(:,2); %12
        Sigincr1(:,5) = -3*QQ2r1(:,1).*QQ2r1(:,3); %13
        Sigincr1(:,6) = -3*QQ2r1(:,2).*QQ2r1(:,3); %23
        
        Sigincr1 =Sigincr1.*deltaV';
        
        sigincr0=sum(Sigincr0);
        sigincr1=sum(Sigincr1);
        
        stressr0(1,1)=sigincr0(1);
        stressr0(2,2)=sigincr0(2);
        stressr0(3,3)=sigincr0(3);
        stressr0(1,2)=sigincr0(4);
        stressr0(1,3)=sigincr0(5);
        stressr0(2,3)=sigincr0(6);
        
        stressr0(2,1)=stressr0(1,2);
        stressr0(3,1)=stressr0(1,3);
        stressr0(3,2)=stressr0(2,3);
        
        
        stressr1(1,1)=sigincr1(1);
        stressr1(2,2)=sigincr1(2);
        stressr1(3,3)=sigincr1(3);
        stressr1(1,2)=sigincr1(4);
        stressr1(1,3)=sigincr1(5);
        stressr1(2,3)=sigincr1(6);
        
        stressr1(2,1)=stressr1(1,2);
        stressr1(3,1)=stressr1(1,3);
        stressr1(3,2)=stressr1(2,3);
        
        sigb0=stressr0*b';
        sigb1=stressr1*b';
        
        finc0(i,:)= cross(sigb0',rmid_r0);
        finc1(i,:)= cross(sigb1',r1_rmid);
        
    end
else
    finc0=zeros(1,3);
    finc1=finc0;
    QQ2r0=zeros(ninc,3);
    QQ2r0mag=zeros(ninc,1);
    QQ2r1=QQ2r0;
    QQ2r1mag=QQ2r0mag;
    
    Sigincr0=zeros(ninc,6);
    Sigincr1=Sigincr0;
    
    stressr0=zeros(3,3);
    stressr1=stressr0;
    
    deltaV=QQ2r0mag;
    
    strain0=0.001; % dialatational misfit strain;
    NU_inc = NU;
    MU_inc = MU*1;
    strain1=strain0*(1+NU_inc)/(1-NU_inc);
    
    r0=segments(linkid,6:8);
    r1=segments(linkid,9:11);
    
    rf0=((r0+r1)/2+r0)/2;
    rf1=((r0+r1)/2+r1)/2;
    
    rmid_r0=(r0+r1)/2-r0;
    r1_rmid=r1-(r0+r1)/2;
    b=segments(linkid,3:5);
    
    QQ2r0=rf0-QQ(:,1:3);
    QQ2r0mag=vecnorm(QQ2r0');
    QQ2r1=rf1-QQ(:,1:3);
    QQ2r1mag=vecnorm(QQ2r1');
    deltaV=2*MU_inc*strain1*RR.^3;
    deltaV=deltaV./(QQ2r0mag.^5);
    
    
    Sigincr0(:,1) = QQ2r0mag.^2' - 3*QQ2r0(:,1).^2; % sigma11
    Sigincr0(:,2) = QQ2r0mag.^2' - 3*QQ2r0(:,2).^2; %22
    Sigincr0(:,3) = QQ2r0mag.^2' - 3*QQ2r0(:,3).^2; %33
    Sigincr0(:,4) = -3*QQ2r0(:,1).*QQ2r0(:,2); %12
    Sigincr0(:,5) = -3*QQ2r0(:,1).*QQ2r0(:,3); %13
    Sigincr0(:,6) = -3*QQ2r0(:,2).*QQ2r0(:,3); %23
    
    Sigincr0=Sigincr0.*deltaV';
    
    Sigincr1(:,1) = QQ2r1mag.^2' - 3*QQ2r1(:,1).^2; % sigma11
    Sigincr1(:,2) = QQ2r1mag.^2' - 3*QQ2r1(:,2).^2; %22
    Sigincr1(:,3) = QQ2r1mag.^2' - 3*QQ2r1(:,3).^2; %33
    Sigincr1(:,4) = -3*QQ2r1(:,1).*QQ2r1(:,2); %12
    Sigincr1(:,5) = -3*QQ2r1(:,1).*QQ2r1(:,3); %13
    Sigincr1(:,6) = -3*QQ2r1(:,2).*QQ2r1(:,3); %23
    
    Sigincr1 =Sigincr1.*deltaV';
    
    sigincr0=sum(Sigincr0);
    sigincr1=sum(Sigincr1);
    
    stressr0(1,1)=sigincr0(1);
    stressr0(2,2)=sigincr0(2);
    stressr0(3,3)=sigincr0(3);
    stressr0(1,2)=sigincr0(4);
    stressr0(1,3)=sigincr0(5);
    stressr0(2,3)=sigincr0(6);
    
    stressr0(2,1)=stressr0(1,2);
    stressr0(3,1)=stressr0(1,3);
    stressr0(3,2)=stressr0(2,3);
    
    
    stressr1(1,1)=sigincr1(1);
    stressr1(2,2)=sigincr1(2);
    stressr1(3,3)=sigincr1(3);
    stressr1(1,2)=sigincr1(4);
    stressr1(1,3)=sigincr1(5);
    stressr1(2,3)=sigincr1(6);
    
    stressr1(2,1)=stressr1(1,2);
    stressr1(3,1)=stressr1(1,3);
    stressr1(3,2)=stressr1(2,3);
    
    sigb0=stressr0*b';
    sigb1=stressr1*b';

    finc0(:)= cross(sigb0',rmid_r0);
    finc1(:)= cross(sigb1',r1_rmid);
end
    
end

