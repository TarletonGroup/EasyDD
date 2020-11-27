% Finite element code for a rectangle subject to any combination of applied
% force and displacement boundary conditions  - agrees with Abaqus CPE4,
% CPS4 elements. Written by E Tarleton, edmund.tarleton@eng.ox.ac.uk
%
% /|---------------------------------------------------.u2 = -U
% /|                                                   |
% /|                                                   |
% /|                                                   |
% /|                                                   | Sright
% /|---------------------------------------------------|

clear all
close all
% input parameters
dx = 30; % x dimension
dy = 6;  % y dimension
mx = 30; % # elements in x direction
%elastic material model
mu = 100 % shear modulus GPa
nu = 0.3 % Poisson's ratio

w = dx/mx; % element width

my = 1+round(mx*dy/dx); % # elements in y direction
h=dy/(my); % element height

mel = mx*my; % # elements
mno=(mx+1)*(my+1); % # of nodes


ey =mu*(2*(1+nu)); % Young's modulus
la =nu*ey/((1+nu)*(1-2*nu)); % Lambda

D = zeros(3,3); 
% matrix of elastic stiffness constants
% plane strain 
D(1,1) = la+ 2*mu;
D(1,2) = la;
D(1,3) = 0;
D(2,2) = la + 2*mu;
D(2,3) = 0;
D(3,3) = mu;

% plane stress
% D(1,1) = ey/(1-nu^2);
% D(1,2) = ey*nu/(1-nu^2);
% D(2,2) = D(1,1);
% D(3,3) = ey/2/(1+nu);
% 
D(2,1) = D(1,2);
D(3,1) = D(1,3);
D(3,2) = D(2,3);


% Gauss point coordinates
zq = 1/sqrt(3);
z = zeros(4,2);
z(1,1:2) = [-zq,-zq];
z(2,1:2) = [zq,-zq];
z(3,1:2) = [zq,zq];
z(4,1:2) = [-zq,zq];
wq=1; % Gauss point weights

nc = zeros(mel,4); %element connectivity
% 4. ----- .3
%  ¦       ¦
%  ¦       ¦
% 1. ----- .2

for i =1:my % i = row index, j = column index
    for j = 1:mx %global node #(element #,local node #)
        ge = j + (i-1)*mx;
        nc(ge,1) = j + i*(mx+1);
        nc(ge,2) = nc(ge,1)+1;
        nc(ge,4) = j + (i-1)*(mx+1);
        nc(ge,3) = nc(ge,4) + 1;
    end
end

% generat the FE mesh
x= zeros(mno,2); %global nodal coordinates gn
% 1         2         ...   mx+1
% 1+(mx+1)  2+(mx+1)  ...   2*(mx+1)
% .          .                  .
% .            .                .
% 1+my*(mx+1) ...          (my+1)*(mx+1)

% region 1
for j =1:mx+1 % nodes in x direction
    for i = 1:my+1
        gn = j + (i-1)*(mx+1); %global node #
        x(gn,1) = (j-1)*w;
        x(gn,2) = dy-(i-1)*h;
    end
end


N = zeros(4,4); 
% N(q,a) is shape function Na(s1,s2) evaluated at integration point q
for q =1:4
    N(q,1) = .25*(1-z(q,1))*(1-z(q,2));
    N(q,2) = .25*(1+z(q,1))*(1-z(q,2));
    N(q,3) = .25*(1+z(q,1))*(1+z(q,2));
    N(q,4) = .25*(1-z(q,1))*(1+z(q,2));
end

Ns = zeros(4,4,2); 
% derivative of shape function Na(s1,s2)  w.r.t. sj evaluated at 
%integration point q  Ns(q,a,j) = (dNa/dsj)s=zq
for q =1:4
    Ns(q,1,1) = -.25*(1-z(q,2));
    Ns(q,1,2) = -.25*(1-z(q,1));
    Ns(q,2,1) =  .25*(1-z(q,2));
    Ns(q,2,2) = -.25*(1+z(q,1));
    Ns(q,3,1) = .25*(1+z(q,2));
    Ns(q,3,2) = .25*(1+z(q,1));
    Ns(q,4,1) = -.25*(1+z(q,2));
    Ns(q,4,2) = .25*(1-z(q,1));
end

J = zeros(2,2,mel,4);
invJ= zeros(2,2); % inv(J)
detJ = zeros(mel,4); % det(J)
Nx = zeros(mel,4,4,2); % derivative of shape functions w.r.t global x,y
B = zeros(3,8,mel,4); % (1:3,# dof/element, # elements, int pts/element)

for p = 1:mel % all elements
    for q = 1:4 % integration points per element
        
        for i =1:2 % DOF
            for j = 1:2 % DOF
                for a = 1:4 %local shape function number
                 J(i,j,p,q) = J(i,j,p,q) + Ns(q,a,j)*x(nc(p,a),i);               
                % sum over a of dNa(s1,s2)/dsj at S=z(q) *xi(local node a of element p) 
                %Jij= dxi/dsj evaluated at element p integration point q                
                end    
            end
        end                            
        
        detJ(p,q) = det(J(:,:,p,q)); % det(J) evaluated at element p int point q

        invJ = inv(J(:,:,p,q));
        
        for a =1:4
            for j =1:2
                for i =1:2
                    Nx(p,q,a,j) = Nx(p,q,a,j)+Ns(q,a,i)*invJ(i,j);
                end
            end
            
            B(1,(a-1)*2+1,p,q)=Nx(p,q,a,1);
            B(1,(a-1)*2+2,p,q)=0;
            
            B(2,(a-1)*2+1,p,q)=0;
            B(2,(a-1)*2+2,p,q)=Nx(p,q,a,2);
            
            B(3,(a-1)*2+1,p,q)=Nx(p,q,a,2);
            B(3,(a-1)*2+2,p,q)=Nx(p,q,a,1);
            
        end
    end
end


disp('local K...')        
Ke = zeros(8,8,mel); %local stiffness matrix for each element
for p =1:mel %all elements 
    for q = 1:4 % int points per element
        Ke(:,:,p) = Ke(:,:,p)+B(:,:,p,q)'*D*B(:,:,p,q)*detJ(p,q)*wq;
    end
end


disp('global K...') 
K = sparse(2*mno,2*mno); % global stiffness matrix

a=1:4; %local node numbers 
dofLocal=[2*(a-1)+1,2*(a-1)+2];
for p =1:mel
    gn=nc(p,a); % global node numbers
    dof(1:8)=[2*(gn-1)+1,2*(gn-1)+2];
    
    K(dof,dof)= K(dof,dof)+Ke(dofLocal,dofLocal,p); %the (full) global stiffness matrix
end


%% -----------------  apply boundary conditions and solution ---------------------
% % reformat K due to applied displacements
f=zeros(mno*2,1);
u = zeros(mno*2,1);
% bending
% Ftotal = -0.1;
% f((mx+1)*2) = Ftotal; % 

Utotal = -0.1780;
u((mx+1)*2) = Utotal; % 

allDofs = [1:2*mno]';
Sleft = find(x(:,1) == 0); % set of nodes on left edge
fixedDofs = [2*Sleft-1;2*Sleft;(mx+1)*2]; % degrees of freedom where u is specified
freeDofs = setdiff(allDofs,fixedDofs);

f = f - K*u;
% rf = K(:,fixedDofs)*u(fixedDofs);
% f = f - rf;
u(freeDofs) = K(freeDofs,freeDofs)\f(freeDofs); % Solve system of equations for nodal values

rf = K*u;
rf(2*(mx+1))
%% postprocessing: plotting of deformed mesh
scaleFactor = 10
xnew(:,1) = x(:,1) + scaleFactor*u([1:mno]*2-1);
xnew(:,2) = x(:,2) + scaleFactor*u([1:mno]*2);

%plot deformed mesh
% plot mesh to check it looks ok
figure(1);clf
axis equal 
hold on;
for p =1:mel
    plot(xnew(:,1),xnew(:,2),'ko') % plot nodes
    plot(xnew(nc(p,[1:4,1]),1),xnew(nc(p,[1:4,1]),2),'k-') % plot elements
end
hold off

%% postprocessing: plotting stress field
figure(2)
clf
axis equal 
colorbar
%caxis([-0.4,0.4])
hold on;
for p =1:mel
    % returns s11,s22,s12 at nodes by interpolating back from Gauss points
    
    strain = zeros(3,4); % strain at Gauss points
    for i =1:3
        for q = 1:4 % Gauss points
            strain(i,q) = B(i,[1,3,5,7],p,q)*u(2*nc(p,:)-1) + B(i,[2,4,6,8],p,q)*u(2*nc(p,:));
        end
    end
    
    stress = D*strain; % Stress at Gasuss points
    
    Z = mean(stress(1,:))*ones(2); % average stress in the element
    surf(x(nc(p,[1,3]),1),x(nc(p,[1,3]),2), Z )

end

hold off

%%
% 1a) nodal displacements at node mx+1: [ux; uy] = u((mx+1)*2-1: (mx+1)*2)
% 1b) comment out line 191-192
% 1c)Utotal = -0.1780; u((mx+1)*2) = Utotal; fixedDofs = [2*Sleft-1;2*Sleft;(mx+1)*2];