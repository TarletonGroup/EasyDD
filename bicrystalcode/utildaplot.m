%%
clear all %#ok<CLALL>
close

%code rn and links
rn=[35  10e5  -10e5  0;
    35  -10e5  10e5  0;
    35  -10e5  10e5+10  0;
    35  10-10e5  10e5+10  0;
    35  10e5+10  10-10e5  0;
    35  10e5+10  -10e5  0];

b  = [0 1 -1]/norm([0 1 -1]);

n  = [1 0 0]/norm([1 0 0]);
% rn=[A1';B1';C1'];
% A1=[2.5;2.5;-3];
% B1=[7.5;5;-3];
% C1=[2.5;7.5;-3];
% links=[ 1 2 b n;
%         2 3 b n;
%         3 1 b n];
links=[1  2  b  n;
       2  3  b  n;
       3  4  b  n;
       4  5  b  n;
       5  6  b  n;
       6  1  b  n];
llinks=size(links,1);
cm=[35;5;0];
nu=0.3;
ny = 11; %number of points in y 
nx = 71; % number of points in x
nz = 11;
nnodes=ny*nx*nz;

Lx = 70;
Ly = 10;
Lz = 10;

dx = Lx/(nx-1);
dy = Ly/(ny-1);
dz = Lz/(nz-1);
xnodes = zeros(nnodes,3);
n=0;
for i =1:ny
    for j =1:nx
        for k = 1:nz
            n = n+1;
            xnodes(n,1) = (j-1)*dx;
            xnodes(n,2) = (i-1)*dy;
            xnodes(n,3) = (k-1)*dz;               
        end
    end
end

utilda=zeros(nnodes,3);

for n=1:nnodes
    p = xnodes(n,:)';
    for m=1:llinks
        A=rn(links(m,1),1:3)';
        B=rn(links(m,2),1:3)';
        C=cm;
        utilda(n,:) = utilda(n,:) + displacement_et(p,A,B,C,b',nu);
        %utilda(n,:) = utilda(n,:) + displacement_bb(p,A,B,C,b',nu);
    end
end

        xp = xnodes+1e1*utilda;
 
    amag=1;
    xp = amag*xp;
    figure;clf;hold on;view(0,0)
    xlabel('x');ylabel('y');zlabel('z')
    %style='-k';
       plot3(xp(:,1),xp(:,2),xp(:,3),'.','markers',18) % plot nodes
      
    axis equal
    zlabel('z-direction (b)');
    xlabel('x-direction (b)');
    ylabel('y-direction (b)');
     %zlim([-10 15])
     %xlim([0.5 50])
    title('$\tilde{u}$','FontSize',14,'Interpreter','Latex');
