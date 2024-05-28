% clear all
% mx=1,dx=6,dy=1,dz=1,mu=1,nu=.3,loading=1; 

function [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D_haiyang2(dx,dy,dz,mx,mu,nu,loading)         

% E Tarleton edmund.tarleton@materials.ox.ac.uk
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
% rectangulare domain.
% (y going in to the page) note this is rotated about x axis w.r.t. local (s1,s2,s3)
% system
% -------------------------------------------------
%           
% ^z                   (mx,my)
% |          
% |           
% |------>x-----------------------------------------
% not that this code uses a constant element size in base and beam: (w,h,d)
% also gammau is the left surface (eg baseLength=0)
loading =1;

w = dx/mx; % elements width

my = round(mx*dy/dx); % # elements in y direction
my=max(my,1);
h=dy/my; % element height

mz = round(mx*dz/dx); % elements in z direction
mz=max(mz,1);
d = dz/mz; % element depth

mel = mx*my*mz;
mno=(mx+1)*(my+1)*(mz+1);

    
ey =mu*(2*(1+nu));
la =nu*ey/((1+nu)*(1-2*nu));

D = zeros(6,6);

for i =1:3
    D(i,i) = la+ 2*mu;
end
for i =4:6
    D(i,i) = mu ;
end
D(1,2) = la;
D(1,3) = la;
D(2,3) = la;

% plane stress
% D(1,1) = ey/(1-nu^2);
% D(1,2) = ey*nu/(1-nu^2);
% D(2,2) = D(1,1);
% D(3,3) = ey/2/(1+nu);
% 
D(2,1) = D(1,2);
D(3,1) = D(1,3);
D(3,2) = D(2,3);


%           
%  ^s2
%  |          
%  |           
%  |------>s1
% 
%(s3 out of page)
zv = 1/sqrt(3);
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
            xnodes(gn,1) = (i-1)*w;
            xnodes(gn,2) = (j-1)*h;
            xnodes(gn,3) = (k-1)*d;
        end
    end

end
   
nc = zeros(mel,4); %element connectivity
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

% figure(1);clf;hold on;view(3)
% xlabel('x');ylabel('y');zlabel('z')
% for p =1:mel
% %     plot3(x(:,1),x(:,2),x(:,3),'.') % plot nodes
%     % plot elements
%     plot3(xnodes(nc(p,[1:4,1]),1),xnodes(nc(p,[1:4,1]),2),xnodes(nc(p,[1:4,1]),3),'-') 
%     plot3(xnodes(nc(p,[5:8,5]),1),xnodes(nc(p,[5:8,5]),2),xnodes(nc(p,[5:8,5]),3),'-') 
%     plot3(xnodes(nc(p,[1,5]),1),xnodes(nc(p,[1,5]),2),xnodes(nc(p,[1,5]),3),'-') % 
%     plot3(xnodes(nc(p,[2,6]),1),xnodes(nc(p,[2,6]),2),xnodes(nc(p,[2,6]),3),'-') % 
%     plot3(xnodes(nc(p,[3,7]),1),xnodes(nc(p,[3,7]),2),xnodes(nc(p,[3,7]),3),'-') % 
%     plot3(xnodes(nc(p,[4,8]),1),xnodes(nc(p,[4,8]),2),xnodes(nc(p,[4,8]),3),'-') % 
% end
% axis('equal')
% hold off

n = zeros(8,8); 
% n(q,a) is shape function Na(s1,s2,s3) evaluated at integration point q
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
detJ = zeros(mel,8); % det(J) at mel elements, 9 quad points/element
nx = zeros(mel,8,8,3); % derivative of shape functions w.r.t global x,y
B = zeros(6,24,mel,8); % (# strain components, # dof/element, # elements, int pts/element)

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
            
            B(1,(a-1)*3+1,p,q)=nx(p,q,a,1);
            B(1,(a-1)*3+2,p,q)=0;
            B(1,(a-1)*3+3,p,q)=0;
            
            B(2,(a-1)*3+1,p,q)=0;
            B(2,(a-1)*3+2,p,q)=nx(p,q,a,2);
            B(2,(a-1)*3+3,p,q)=0;
            
            B(3,(a-1)*3+1,p,q)=0;
            B(3,(a-1)*3+2,p,q)=0;
            B(3,(a-1)*3+3,p,q)=nx(p,q,a,3);
            
            B(4,(a-1)*3+1,p,q)=nx(p,q,a,2);
            B(4,(a-1)*3+2,p,q)=nx(p,q,a,1);
            B(4,(a-1)*3+3,p,q)=0;
            
            B(5,(a-1)*3+1,p,q)=nx(p,q,a,3);
            B(5,(a-1)*3+2,p,q)=0;
            B(5,(a-1)*3+3,p,q)=nx(p,q,a,1);
            
            B(6,(a-1)*3+1,p,q)=0;
            B(6,(a-1)*3+2,p,q)=nx(p,q,a,3);
            B(6,(a-1)*3+3,p,q)=nx(p,q,a,2);
            
        end
    end
end


disp('local K...');        
ke = zeros(24,24,mel); %local stiffness matrix
for p =1:mel %all elements 
    for q = 1:8 % int points per element
        ke(:,:,p) = ke(:,:,p)+B(:,:,p,q)'*D*B(:,:,p,q)*detJ(p,q);
    end
end

% ensure ke is symmetric eg remove any very small entries due to numerical
% error.
% disp('global K...');
% 
% kg = sparse(mno*3,mno*3);
% tic
% a =1:8;
% b =1:8;
% for p =1:mel
% %     disp(['assembly ', num2str(100*p/mel), ' % complete'])
%      percentage = 100*p/mel;
%      if mod(percentage,1)==0
%          fprintf('Assembly %d percent complete \n',100*p/mel);
%      end
%     gna=nc(p,a);
%     gnb=nc(p,b);
%     for i =1:3
%         for j=1:3        
%             kg(3*(gna-1)+i,3*(gnb-1)+j)= kg(3*(gna-1)+i,3*(gnb-1)+j)+ke(3*(a-1)+i,3*(b-1)+j,p);
%         end
%     end
% end
% toc

% kg = sparse(mno*3,mno*3);
% tic;
% a =1:8;
% % b =1:8;
% for p =1:mel
%     percentage = 100*p/mel;
%     if mod(percentage,1)==0
%         fprintf('Assembly %d percent complete \n',100*p/mel);
%     end
%     gn=nc(p,a);
% %     gnb=nc(p,b);
%     dof(1:24)=[3*(gn-1)+1,3*(gn-1)+2,3*(gn-1)+3];
% %     dofb(1:24)=[3*(gnb-1)+1,3*(gnb-1)+2,3*(gnb-1)+3]
%     dofLocal=[3*(a-1)+1,3*(a-1)+2,3*(a-1)+3];
% %     dofbl=[3*(b-1)+1,3*(b-1)+2,3*(b-1)+3]
%     kg(dof,dof)= kg(dof,dof)+ke(dofLocal,dofLocal,p);
%     
% end
% toc

% %HY20171204: modified by Haiyang following Ed's instruction
% % see http://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/
% % instead of K(i,k) = K(i,j) + ... construct I J X tripletes and use sparse
% % as it's much faster O(nlog(n)) compared to O(n^2)
tic
disp('global stifness matrix assembly new (fast) method...')
a=1:8; %local node numbers 
dofLocal=[3*(a-1)+1,3*(a-1)+2,3*(a-1)+3];
ntriplets = 3*mno;
I = zeros(ntriplets,1);
J = zeros(ntriplets,1);
X = zeros(ntriplets,1);
ntriplets = 0;

for p =1:mel
    gn=nc(p,a); % global node numbers
    dof(1:24)=[3*(gn-1)+1,3*(gn-1)+2,3*(gn-1)+3]; % global degree of freedom
    for i =1:24
        for j =1:24
            ntriplets = ntriplets + 1;
            I(ntriplets) = dof(i);
            J(ntriplets) = dof(j);
            X(ntriplets) = ke(dofLocal(i),dofLocal(j),p);
        end
    end
  %  K{N}(dof,dof)= K{N}(dof,dof)+Ke(dofLocal,dofLocal,p); %the (full) global stiffness matrix

end

kg= sparse(I,J,X,3*mno,3*mno); %the (full) global stiffness matrix
toc

%------------------------------ boundary conditions ---------------------
disp('reformatting K');
% reformat kg due to applied displacements
%f=zeros(mno*2,1);
K=kg;
bcwt=mean(diag(kg));%=trace(K)/length(K)
bcwt = full(bcwt);

if loading == 1||loading==0 % cantilever bending
    % define empty sets (not used but allows same code as for bending)       
    % priority for inclusion of shared nodes is in order below
    %S = [node #, nodal Area, outward normal]
    Smixed = zeros(my+1,5); % top right edge (u3=U) first priority
    
    Sleft = zeros((my+1)*(mz+1),5);  % second priority
    Sright = zeros((my+1)*mz,3); 
    
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
                L3 = 0.5*d;            
            else
                L3 = d;
            end
            
            if j ==1 || j == my+1
                L2 = 0.5*h;
            else
                L2 = h;
            end 
            
            Sleft(m,2) = L2*L3;
            Sleft(m,3:5) = [-1,0,0];
            m=m+1;
        end
    end
    
    m=1;
    for j =1:my+1
        for k =1:mz
            % right surface of beam (t=0)
            gn = k*(mx+1)+(mx+1)*(mz+1)*(j-1);                 
            Sright(m,1)= gn;
                        
            if k == 1
                L3 = 0.5*d;            
            else
                L3 = d;
            end
            
            if j ==1 || j == my+1
                L2 = 0.5*h;
            else
                L2 = h;
            end 
            
            Sright(m,2) = L2*L3;
            Sright(m,3:5) = [1,0,0];
            m=m+1;
        end               
    end
    
    m=1;
    for j =1:my+1
        % top right loaded edge (t1=0, u2=-U) mixed surface
        gn = (mx+1)*(mz+1)*j;
        Smixed(m,1) = gn;
        
        if j ==1 || j == my+1
            L2 = 0.5*h;
        else
            L2 = h;
        end
        
        L3 = 0.5*d;
        Smixed(m,2) = L2*L3;
        Smixed(m,3:5) = [1,0,0];
        m=m+1;
    end
      
    
    m=1; 
    for j = 1:my+1
        for i=2:mx
            % bottom surface of beam (t=0)
            gn = i+(mx+1)*(mz+1)*(j-1);
            Sbot(m,1) = gn;
            % top surface of beam (t=0)
            gn = i+(mx+1)*mz+(mx+1)*(mz+1)*(j-1);
            Stop(m,1) = gn;
            
            L1 = w;
            if j == 1 || j == my+1
                L2 = 0.5*h;
            else
                L2 = h;
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
            
            L1 = w;
            L3 = d;
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

if loading ==1 % end loaded cantilever
           
    gammat=[Stop;Sbot;Sright;Sfront;Sback]; % t=0         
    gammau=Sleft; %  set containg fixed nodes u=0    
    gammaMixed=Smixed; % t1=t2=0, u3= U
    
    fixedDofs =[3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1); 3*gammaMixed(:,1)];
    freeDofs = [3*gammat(:,1)-2; 3*gammat(:,1)-1; 3*gammat(:,1); 3*gammaMixed(:,1)-2;...
        3*gammaMixed(:,1)-1];
   
elseif loading == 0 % used for debugging
    
     Sright = [Sright;Smixed]; 
%     
%    
%     gnr = (mx+1);                   % x(gnr) = [dx 0 0] u3 = 0 -> R(x2)=0; u2 = 0 -> R(x3)=0 
%     gnl = [(mx+1)*(mz+1)*my+1];     % x(gnl) = [0 dy 0] u3 = 0 -> R(x1)=0                                 
%     
%     Smixed = [Sright(ismember(Sright(:,1), gnr),:); Sleft(ismember(Sleft(:,1),gnl),:)];
%      
%     
%     gammat=[Stop;Sbot;Sfront;Sback;...
%         Sright(~ismember(Sright(:,1), gnr),:); Sleft(~ismember(Sleft(:,1),[1,gnl]),:)]; % t=0            
%    
%     gammau = Sleft(ismember(Sleft(:,1),1),:); 
%     
%     gammaMixed=Smixed;
% fixedDofs =[3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1);...
%         3*gammaMixed(1,1);  3*gammaMixed(1,1)-1; 3*gammaMixed(2,1)];    
%   freeDofs = [3*gammat(:,1)-2; 3*gammat(:,1)-1; 3*gammat(:,1);...
%         3*gammaMixed(1,1)-2; 3*gammaMixed(2,1)-2; 3*gammaMixed(2,1)-1];

    Smixed = [];
    gammaMixed=Smixed;
    gammat = [];
    gammau = [Stop;Sbot;Sleft;Sright;Sfront;Sback]; 
    
    fixedDofs =[3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1)];        
    
    freeDofs = [];
elseif loading == 4 % HY201130: modified for relaxation to get initial structure
    Smixed = [];
    gammaMixed=Smixed;
    gammat = [Stop;Sbot;Sleft;Sright;Sfront;Sback];
    gammau = []; 
    
    fixedDofs =[];        
    
    freeDofs = [3*gammau(:,1)-2; 3*gammau(:,1)-1; 3*gammau(:,1)];
else
   disp('loading not specified correctly')
   pause
end

% for m = 1:length(fixedDofs)    
%     i = fixedDofs(m);
%     K(:,i) = 0;
%     K(i,:) = 0;
%     K(i,i) = bcwt; 
% end
% 
%  if length([fixedDofs;freeDofs])>length(unique([fixedDofs;freeDofs]))
%         disp('error')
%         pause    
%  end

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
%pause
% disp('LU decomposition of K...')
% %L=0;U=0;
% tic;
% [L,U] = lu(K);    
% toc;

% disp('Cholesky Factorization of K...'); %should be symmetric!
% tic;
% U = chol(K);
% L = U';
% toc;
%pause

% if max(abs(diag( U\(L\K) )))-1 > 1000*eps
%     disp('Error in inverse K')
%     pause
% end

%HY20171206:********************************************************
%HY20171206: modified by HY to make the code cleaner by removing the
%equations related to the fixedDofs; since FreeDofs has been used to
%represent the free boundary nodes, a new term, unfixedDofs, is used to
%represent all the nodes other than the fixed ones. i.e.
%unfixedDofs=allDofs - fixedDofs

allDofs = [1:3*mno];
unfixedDofs = setdiff(allDofs,fixedDofs); %HY20171206: setdiff not only obtain the different elements but also sort them.
Kred = K(unfixedDofs,unfixedDofs);
% tic
% disp('LU decomposition of Kred')
% [Lred, Ured] = lu(Kred);
% toc
disp('Cholesky Factorization of Kred...'); %should be symmetric!
tic;
Ured = chol(Kred);
Lred = Ured';
toc;
L=[];
U=[];
%HY20171206:********************************************************
disp('finished FEM')
%-------------------------------------------------------------------
