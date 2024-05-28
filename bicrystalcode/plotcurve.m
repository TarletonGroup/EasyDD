
MU=0.83e11;                      
NU=0.29; 
Modulus = 2*(1+NU)*MU;

mumag=8.3e4;
NU=0.29;
amag= 2.87e-4; 
amag=sqrt(3)/2*amag;

dx=2.0/amag; %1.2micron
dy=2.0/amag; %0.6micron
dz=4.0/amag; %0.6micron

load('.\output\output_result.dat')
x=output_result(:,1); %s


y3=output_result(:,4)./dx/dy*mumag; % Ftop  MPa
elasticstrain = y3./(2*(1+NU)*mumag);
y2=output_result(:,3)./dz - elasticstrain; % Uend/dz - elastic
y1=output_result(:,2) + elasticstrain ; %u/l 
y4=output_result(:,3)./dz*2 -elasticstrain  ; %u/l 

y5=output_result(:,5); %density
y6=output_result(:,6);%ratio climb steps/glide steps
y7=output_result(:,7);%ratio climb strain/glide strain


plotyy(x,-y2*100,x,y5);

plot(x,-y4*100,'LineWidth',2);
xlabel('Time (s)','FontSize',14)
ylabel('Plastic strain','FontSize',14)



plotyy(x,-y4,x,y5);

yyaxis left
plot(x,-y1*100,'k',x,-y4*100,'b','LineWidth',2);
xlabel('Time (s)','FontSize',14)
ylabel('Total strain (100%)','FontSize',14)
% legend({'\rho bv + elasticstrain','u/l',},'FontSize',14);


yyaxis right
plot(x,y5, '-.r','LineWidth',2)
ylabel('Density \times 10^{12}/m^2','FontSize',14)
legend({'\rho bv + elasticstrain','u/l','density'},'FontSize',14);

