% created by gaoyuan 09-10-14
% modified by gaoyuan 09-10-21

function [rn,links]=inicon_loop()

 radius = 5*[20 12 4]/0.247;
%  radius =160*ones(1,1)/0.247;
d=-15; % distance in the habit plane
distance=0/0.247; %distance in perpendicular 
% height=600;

% % beta=45/180*pi;
% % 
% % burg=[0 -sin(beta) cos(beta)];
% normplane = [0 0 1];
theta  =20;   % angle of each segment
nloop = size(radius,2); %  number of loops

% for i=1:nloop
%     beta(i)=(-1)^i*pi/4;
% 
% end

   beta=zeros(1,nloop);
%   beta=[0,pi/2];

for i=1:nloop
    burg(i,1:3)=[0 -sin(beta(i)) cos(beta(i))];
end



rn = zeros(nloop*360/theta,4);
links = zeros(nloop*360/theta,4);
%r_loop=radius;

x = 700;
y= 700;
z= 500;
% z=distance;
% z=0;
for iloop=1:nloop
    ii = 1 + (iloop-1)*360/theta;
    rn(ii,1:3) = [x + 1*radius(iloop), y+ 0*radius(iloop), z];
    rn(ii,1:3)=rn(ii,1:3)*[1 0 0;0 cos(beta(iloop)) sin(beta(iloop));0 -sin(beta(iloop)) cos(beta(iloop))];

    
    for i=ii:(ii+360/theta-1)

       
        if i==ii+360/theta-1
            links(i,1:2) = [i 1 + (iloop-1)*360/theta];
            vector=rn(ii,1:3)-rn(i,1:3);
            
%              burg*vector'
%              pause

            normplane=-cross(burg(iloop,:),vector)/norm(cross(vector, burg(iloop,:)));
            
        else
            rn(i+1,1:3) = [x + 1*radius(iloop)*cos(i*theta*pi/180),y+(radius(iloop)*1.0)*sin(i*theta*pi/180)*1,z];
            rn(i+1,1:3)=rn(i+1,1:3)*[1 0 0;0 cos(beta(iloop)) sin(beta(iloop));0 -sin(beta(iloop)) cos(beta(iloop))];
            vector=rn(i+1,1:3)-rn(i,1:3);
%             +height*sin(i*theta*pi/180)
 %           burg =burg*(-1)^iloop;
            normplane=-cross(burg(iloop,:),vector)/norm(cross(vector, burg(iloop,:)));
%          +height*sin(i*theta*pi/180)  
            links(i,1:2) = [i i+1];

        end

        links(i,3:5) = burg(iloop,:);
 %       links(i,5) = burg(3)*(-1)^(1+iloop);
        links(i,6:8) = normplane(1:3);        
    end
    
    if iloop ~= nloop
         x = x + 1.0*(radius(iloop) + radius(iloop+1))+d/0.247;
         y = y + 1.0*(radius(iloop) + radius(iloop+1))+d/0.247;
%          x = 2*x +radius(iloop+1)-40;
       z=z+distance;
      
    end
end
%  rn(:,1:3) = rn(:,1:3)/0.256;
% rn
% links
% fid = fopen('.\output\Initialall.plt','w');
% fprintf(fid,'TITLE = "Dislocation configuation of step%d"\n',1);
% fprintf(fid,'VARIABLES = "X", "Y","Z"\n');
% fprintf(fid,'ZONE  N=%d,E=%d,F=FEPOINT, ET=LINESEG\n',length(rn(:,1)),length(links(:,1)) );
% rnn=rn(:,1:3)';
% fprintf(fid,'%f     %f      %f\n',rnn);
% linksn=links(:,1:2)';
% fprintf(fid,'%d     %d\n',linksn);
% fclose(fid);

% 
hold on
plim =5000;
viewangle = [135 45];
plotnodes(rn,links,2*plim); 
view(viewangle);
% xlim([0 plim]);
% ylim([0 plim]); 
% zlim([0 plim]);
hold on


