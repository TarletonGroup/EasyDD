% clear all
% mx=1,dx=6,dy=1,dz=1,mu=1,nu=.3,loading=1; 

function [B,xnodes,mno,nc,n,ke,gammaj,gammaMiu,wx,wy,wz,mel,my,mz] = finiteElement3D_new(dx,dy,dz,mx,Dv)         

% fengxian Liu, 15/04/2020
%calculate the FEM elemental matrix ke, and the node set gammaJ and
%gammaMiu at boundary
% 3D FEM code using linear 8 node element with 8 integration pts (2x2x2) per element
% 4. ----- .3
%  �\       �\
%  � \      � \
% 1. -\---- .2 \ 
%   \  \     \  \
%    \ 8. ----\- .7
%     \ �      \ �
%      \�       \�
%      5. ----- .6
%
%
% rectangular domain.
% (y going in to the page) note this is rotated about x axis w.r.t. local (s1,s2,s3)
% system
% -------------------------------------------------
%           
% ^z                   (mx,my)
% |          
% |           
% |------>x-----------------------------------------
% not that this code uses a constant element size in base and beam: (wx,wy,wz)
% bottom surface x-y is fixed

loading =1; % different loading condition

wx = dx/mx; % elements width

my = round(mx*dy/dx); % # elements in y direction
my=max(my,1);
wy=dy/my; % element height

mz = round(mx*dz/dx); % elements in z direction
mz=max(mz,1);
wz = dz/mz; % element depth

mel = mx*my*mz; % number of element 
mno=(mx+1)*(my+1)*(mz+1);  % number of node

D = Dv; %diffusivity

zv = 1/sqrt(3);  %gauss quadrature
z = zeros(8,3);
z(1,1:3) = [-zv,-zv,-zv];
z(2,1:3) = [ zv,-zv,-zv];
z(3,1:3) = [ zv, zv,-zv];
z(4,1:3) = [-zv, zv,-zv];
z(5,1:3) = [-zv,-zv, zv];
z(6,1:3) = [ zv,-zv, zv];
z(7,1:3) = [ zv, zv, zv];
z(8,1:3) = [-zv, zv, zv];

xnodes= zeros(mno,3); %global nodal coordinates gn

%back face y = dy, C = (mx+1)*(mz+1)*my
%1+mz*(mx+1)+C ...          (mx+1)*(mz+1)+C 
% .          .                  .
% .            .                .
% 1+(mx+1)+C  2+(mx+1)+C  ...   2*(mx+1)+C
% 1+C         2+C         ...     mx+1 + C

%    front face y = 0
%    1+mz*(mx+1) ...          (mx+1)*(mz+1) 
%    .          .                  .
%    .            .                .
%    1+(mx+1)  2+(mx+1)  ...   2*(mx+1)
%    1         2         ...     mx+1 

for k = 1:mz+1
    for j =1:my+1 % nodes in x direction
        for i = 1:mx+1
            gn = i + (k-1)*(mx+1) + (j-1)*(mx+1)*(mz+1); %global node #
            xnodes(gn,1) = (i-1)*wx;
            xnodes(gn,2) = (j-1)*wy;
            xnodes(gn,3) = (k-1)*wz;
        end
    end

end
   
nc = zeros(mel,8); %element connectivity
% 4. ----- .3
%  �\       �\
%  � \      � \
% 1. -\---- .2 \ 
%   \  \     \  \
%    \ 8. ----\- .7
%     \ �      \ �
%      \�       \�
%      5. ----- .6



for k =1:mz
    for j =1:my % i = row index, j = column index
        for i = 1:mx %global node #(element #,local node #)
            ge = i + (k-1)*(mx) + (j-1)*mx*mz; % element number
            nc(ge,1) = i + (k-1)*(mx+1) + (mx+1)*(mz+1)*j;
            nc(ge,2) = nc(ge,1)+1;
            nc(ge,4) = i +  k*(mx+1) + (mx+1)*(mz+1)*j;
            nc(ge,3) = nc(ge,4) + 1;            
            nc(ge,5) = i + (k-1)*(mx+1) + (mx+1)*(mz+1)*(j-1);
            nc(ge,6) = nc(ge,5)+1;
            nc(ge,8) = i + k*(mx+1) + (mx+1)*(mz+1)*(j-1);
            nc(ge,7) = nc(ge,8) + 1;
        end
    end
end

% figure(1);clf;
% hold on;view(3)
% xlabel('x');ylabel('y');zlabel('z')
% for p =1:mel
% %     plot3(x(:,1),x(:,2),x(:,3),'.') % plot nodes
%     % plot elements
%     plot3(xnodes(nc(p,[1:4,1]),1),xnodes(nc(p,[1:4,1]),2),xnodes(nc(p,[1:4,1]),3),'b-') 
%     plot3(xnodes(nc(p,[5:8,5]),1),xnodes(nc(p,[5:8,5]),2),xnodes(nc(p,[5:8,5]),3),'b-') 
%     plot3(xnodes(nc(p,[1,5]),1),xnodes(nc(p,[1,5]),2),xnodes(nc(p,[1,5]),3),'b-') % 
%     plot3(xnodes(nc(p,[2,6]),1),xnodes(nc(p,[2,6]),2),xnodes(nc(p,[2,6]),3),'b-') % 
%     plot3(xnodes(nc(p,[3,7]),1),xnodes(nc(p,[3,7]),2),xnodes(nc(p,[3,7]),3),'b-') % 
%     plot3(xnodes(nc(p,[4,8]),1),xnodes(nc(p,[4,8]),2),xnodes(nc(p,[4,8]),3),'b-') % 
% end
% axis('equal')
% hold off

n = zeros(8,8); 
% n(q,a) is shape function Na(s1,s2,s3) evaluated at integration point z
for q =1:8
    n(q,1) = 1/8*(1-z(q,1))*(1-z(q,2))*(1-z(q,3));
    n(q,2) = 1/8*(1+z(q,1))*(1-z(q,2))*(1-z(q,3));
    n(q,3) = 1/8*(1+z(q,1))*(1+z(q,2))*(1-z(q,3));
    n(q,4) = 1/8*(1-z(q,1))*(1+z(q,2))*(1-z(q,3));
    n(q,5) = 1/8*(1-z(q,1))*(1-z(q,2))*(1+z(q,3));
    n(q,6) = 1/8*(1+z(q,1))*(1-z(q,2))*(1+z(q,3));
    n(q,7) = 1/8*(1+z(q,1))*(1+z(q,2))*(1+z(q,3));
    n(q,8) = 1/8*(1-z(q,1))*(1+z(q,2))*(1+z(q,3));
end


ns = zeros(8,8,3); 
% derivative of shape function Na(s1,s2,s3)  w.r.t. sj evaluated at 
%integration point q  ns(q,a,j) = (dNa/dsj)s=zq
for q =1:8
    ns(q,1,1) = -1/8*(1-z(q,2))*(1-z(q,3));
    ns(q,1,2) = -1/8*(1-z(q,1))*(1-z(q,3));
    ns(q,1,3) = -1/8*(1-z(q,1))*(1-z(q,2));
    
    ns(q,2,1) =  1/8*(1-z(q,2))*(1-z(q,3));
    ns(q,2,2) = -1/8*(1+z(q,1))*(1-z(q,3));
    ns(q,2,3) = -1/8*(1+z(q,1))*(1-z(q,2));
    
    ns(q,3,1) = 1/8*(1+z(q,2))*(1-z(q,3));
    ns(q,3,2) = 1/8*(1+z(q,1))*(1-z(q,3));
    ns(q,3,3) = -1/8*(1+z(q,1))*(1+z(q,2));
    
    ns(q,4,1) = -1/8*(1+z(q,2))*(1-z(q,3));
    ns(q,4,2) = 1/8*(1-z(q,1))*(1-z(q,3));
    ns(q,4,3) = -1/8*(1-z(q,1))*(1+z(q,2));
    
    ns(q,5,1) = -1/8*(1-z(q,2))*(1+z(q,3));
    ns(q,5,2) = -1/8*(1-z(q,1))*(1+z(q,3));
    ns(q,5,3) = 1/8*(1-z(q,1))*(1-z(q,2));
    
    ns(q,6,1) =  1/8*(1-z(q,2))*(1+z(q,3));
    ns(q,6,2) = -1/8*(1+z(q,1))*(1+z(q,3));
    ns(q,6,3) = 1/8*(1+z(q,1))*(1-z(q,2));
    
    ns(q,7,1) = 1/8*(1+z(q,2))*(1+z(q,3));
    ns(q,7,2) = 1/8*(1+z(q,1))*(1+z(q,3));
    ns(q,7,3) = 1/8*(1+z(q,1))*(1+z(q,2));
    
    ns(q,8,1) = -1/8*(1+z(q,2))*(1+z(q,3));
    ns(q,8,2) = 1/8*(1-z(q,1))*(1+z(q,3));
    ns(q,8,3) = 1/8*(1-z(q,1))*(1+z(q,2));
end



% validate ns：

pm1 =[-1  1  1 -1 -1  1  1 -1];
pm2 =[-1 -1  1  1 -1 -1  1  1];
pm3 =[-1 -1 -1 -1  1  1  1  1];

dNds=zeros(8,8,3); 
% Na(s1,s2,s3) = 1/8 *(1+pm1(a)s1)(1+pm2(a)s2)(1+pm3(a)s3)
for q=1:8
    for a=1:8
        dNds(q,a,1) = 1/8*pm1(a)*(1+pm2(a)*z(q,2))*(1+pm3(a)*z(q,3));
        dNds(q,a,2) = 1/8*(1+pm1(a)*z(q,1))*pm2(a)*(1+pm3(a)*z(q,3));
        dNds(q,a,3) = 1/8*(1+pm1(a)*z(q,1))*(1+pm2(a)*z(q,2))*pm3(a);
    end
end

if any(dNds~=ns) 
    'error'
    pause
end

J = zeros(3,3,mel,8); % 3Dx3D, mel elements, 8 quad pts/element
invJ= zeros(3,3); % inv(J)
detJ = zeros(mel,8); % det(J) at mel elements, 8 quad points/element
nx = zeros(mel,8,8,3); % derivative of shape functions w.r.t global x,y
% B = zeros(6,24,mel,8); % (# strain components, # dof/element, # elements, int pts/element)
B = zeros(3,8,mel,8);  % (# nodal potential, # dof/element, # elements, int pts/element)

for p = 1:mel % all elements
    for q = 1:8 % integration points per element
        
        for i =1:3 % DOF
            for j = 1:3 % DOF
                 J(i,j,p,q) = ns(q,1:8,j)*xnodes(nc(p,1:8),i); 
                 % implied sum over a: local shape function number                
                % sum over a of dNa(s1,s2,s3)/dsj at S=z(q) *xi(local node a of element p) 
                %Jij= dxi/dsj evaluated at element p integration point q                
            end
        end                            
        
        detJ(p,q) = det(J(:,:,p,q)); % det(J) evaluated at element p int point q

        invJ = inv(J(:,:,p,q)); % invJ_ij = dsi/dxj
        
        for a =1:8
            for j =1:3
                for i =1:3
                    nx(p,q,a,j) = nx(p,q,a,j)+ns(q,a,i)*invJ(i,j);
                end
            end
            
            B(1,a,p,q)=nx(p,q,a,1);

            B(2,a,p,q)=nx(p,q,a,2);

            B(3,a,p,q)=nx(p,q,a,3);
                        
        end
    end
end


disp('local K...');        
ke = zeros(8,8,mel); %local stiffness matrix
for p =1:mel %all elements 
    for q = 1:8 % int points per element
        ke(:,:,p) = ke(:,:,p)+B(:,:,p,q)'*D*B(:,:,p,q)*detJ(p,q);
    end
end

% ensure ke is symmetric eg remove any very small entries due to numerical
% error.

%------------------------------ boundary conditions ---------------------
if loading == 1 % bulk diffusion
    % define empty sets (not used but allows same code as for bending)       
    % priority for inclusion of shared nodes is in order below
    %S = [node #, nodal Area, outward normal]
    
    Sleft = zeros((my+1)*(mz+1),5);  % second priority
    Sright = zeros((my+1)*(mz+1),5); 
    
    Stop =  zeros((mx-1)*(my+1),5);  % third priority
    Sbot =  zeros((mx-1)*(my+1),5);  

    Sfront = zeros((mz-1)*(mx-1),5); % fourth priority
    Sback = zeros((mz-1)*(mx-1),5);
    
   
    m=1;
    for j =1:my+1
        for k =1:mz+1 
            % left surface of beam (u=0)
            gn = 1+(k-1)*(mx+1)+(mx+1)*(mz+1)*(j-1);
            Sleft(m,1)= gn; 
            
            if k == 1 || k == mz+1
                L3 = 0.5*wz;            
            else
                L3 = wz;
            end
            
            if j ==1 || j == my+1
                L2 = 0.5*wy;
            else
                L2 = wy;
            end 
            
            Sleft(m,2) = L2*L3;
            Sleft(m,3:5) = [-1,0,0];
            m=m+1;
        end
    end
    
    m=1;
    for j =1:my+1
        for k =1:mz+1
            % right surface of beam (t=0)
            gn = k*(mx+1)+(mx+1)*(mz+1)*(j-1);                 
            Sright(m,1)= gn;
                        
            if k == 1
                L3 = 0.5*wz;            
            else
                L3 = wz;
            end
            
            if j ==1 || j == my+1
                L2 = 0.5*wy;
            else
                L2 = wy;
            end 
            
            Sright(m,2) = L2*L3;
            Sright(m,3:5) = [1,0,0];
            m=m+1;
        end               
    end
       
    m=1; 
    for j = 1:my+1
        for i=2:mx
            %bottom surface of beam (t=0)
            gn = i+(mx+1)*(mz+1)*(j-1);
            Sbot(m,1) = gn;
            %top surface of beam (t=0)
            gn = i+(mx+1)*mz+(mx+1)*(mz+1)*(j-1);
            Stop(m,1) = gn;
            
            L1 = wx;
            if j == 1 || j == my+1
                L2 = 0.5*wy;
            else
                L2 = wy;
            end
            
            Sbot(m,2) = L1*L2;
            Sbot(m,3:5) = [0,0,-1];
            Stop(m,2) = L1*L2;
            Stop(m,3:5) = [0,0,1];            
            m=m+1;            
        end
    end
    
    m=1;
    for k=2:mz
        for i = 2:mx
            %front surface of beam (t=0)  
            gn = i + (k-1)*(mx+1);
            Sfront(m) = gn; %(y=0)
            
            gn =  i + (k-1)*(mx+1) + (mx+1)*(mz+1)*my;
            Sback(m) = gn; %(y=dy)
            
            L1 = wx;
            L3 = wz;
            Sfront(m,2) = L1*L3;
            Sfront(m,3:5) = [0,-1,0];
            Sback(m,2) = L1*L3;
            Sback(m,3:5) = [0,1,0];
            
            m=m+1;
        end
    end
    
else
    disp('error; only coded for end loadeded cantilever')   
    pause
end

if loading ==1 % bulk diffusion
           
    Ssurf=[Stop;Sbot;Sright;Sfront;Sback;Sleft]; % t=0         
    gammaj=Ssurf; %  flux across external surfaces 
    gammaMiu=[]; % chemical potential on internal surfaces
else
   disp('loading not specified correctly')
   pause
end


% {freeDofs,fixedDofs} should contain every degree of freedom on boundary +
% any internal Dofs with a force or displacement specified.
% figure(2);clf;hold on;view(3)
% xlabel('x');ylabel('y');zlabel('z')
% plot3(xnodes(Stop(:,1),1),xnodes(Stop(:,1),2),xnodes(Stop(:,1),3),'r*')
% plot3(xnodes(Sbot(:,1),1),xnodes(Sbot(:,1),2),xnodes(Sbot(:,1),3),'r*')
% plot3(xnodes(Sright(:,1),1),xnodes(Sright(:,1),2),xnodes(Sright(:,1),3),'b.')
% plot3(xnodes(Sleft(:,1),1),xnodes(Sleft(:,1),2),xnodes(Sleft(:,1),3),'b.')
% plot3(xnodes(Sfront(:,1),1),xnodes(Sfront(:,1),2),xnodes(Sfront(:,1),3),'k*')
% plot3(xnodes(Sback(:,1),1),xnodes(Sback(:,1),2),xnodes(Sback(:,1),3),'k*')
% % plot3(xnodes(Smixed(:,1),1),xnodes(Smixed(:,1),2),xnodes(Smixed(:,1),3),'g*')
% axis('equal')
% hold off


%HY20171206:********************************************************
disp('finished FEM')
