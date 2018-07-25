% This is badly coded as there is a lot of hard code but it is a proof of 
% concept for an idealised case.

mu = 1;
nu = 0.29;
a = 1E-6;

%  Same numbering scheme sas the analytical solution.
%  x5------x6
%  |       |
%  |       |
%  x3------x4

x3 = [-1 -1 0];
x4 = [ 1 -1 0];
x5 = [-1  1 0];
x6 = [ 1  1 0];

% There seems to be a factor of 3 missing from the analytical solutions.
x1 = [0 0 1];
x2 = [0 0 1E6];
b  = [1 0 0];

% Gauss-Legrende quadrature seems to be numerically unstable w.r.t. 
% the number of quadrature points for this problem. I've tried with other 
% dislocation lines and they have different behaviours.

% There seems to be a factor of 3 missing from the analytical solutions.
pts = 101;
[fx3n, fx4n, fx5n, fx6n, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, pts);

[fx3a, fx4a, ...
 fx5a, fx6a, ...
 ftota] = ...
   nodal_surface_force_linear_rectangle_mex(x1, x2, ...
      x3, x4, x5, x6, b, mu, nu, a, 1, 1);
 fx3a = fx3a';
 fx4a = fx4a';
 fx5a = fx5a';
 fx6a = fx6a';
 ftota = ftota'
 ftotn
 
 ex3 = (fx3n - fx3a)./fx3a;
 ex4 = (fx4n - fx4a)./fx4a;
 ex5 = (fx5n - fx5a)./fx5a;
 ex6 = (fx6n - fx6a)./fx6a;
 extot = (ftotn - ftota)./ftota;
 
 ex3(isnan(ex3)) = 0;
 ex4(isnan(ex4)) = 0;
 ex5(isnan(ex5)) = 0;
 ex6(isnan(ex6)) = 0;
 extot(isnan(extot)) = 0;
 
 ex3(isinf(ex3)) = 0;
 ex4(isinf(ex4)) = 0;
 ex5(isinf(ex5)) = 0;
 ex6(isinf(ex6)) = 0;
 extot(isinf(extot)) = 0;
%  
%  ex3
%  ex4
%  ex5
%  ex6
%  extot

% pts = 1;
% [fx3, fx4, fx5, fx6, ftot] = f_num(x1, x2, x3, x4, x5, x6, b, n, mu, nu, a, pts)

function [N1, N2, N3, N4] = shape_func(s1, s2)
    N1 = 0.25*(1-s1)*(1-s2);
    N2 = 0.25*(1+s1)*(1-s2);
    N3 = 0.25*(1+s1)*(1+s2);
    N4 = 0.25*(1-s1)*(1+s2);
end %function

function [fx1, fx2, fx3, fx4, ftot] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, n_q)
     n = cross(x4-x3, x5-x3);
     n = n/norm(n);

     Lx = (x4(1)-x3(1));
     Ly = (x5(2)-x3(2));
     
     LxLy = Lx*Ly;
     factor = 1/4;
     
     fx1 = zeros(1,3);
     fx2 = zeros(1,3);
     fx3 = zeros(1,3);
     fx4 = zeros(1,3);
     ftot = zeros(1,3);
     
     [q, w] = lgwt(n_q, -1, 1);
     
     for i = 1: n_q
         y = q(i);
         wy = w(i);
         for j = 1: n_q
             x = q(j);
             wx = w(j);
             [N1, N2, N3, N4] = shape_func(x, y);
             
             [sxx, syy, szz,...
              sxy, syz, sxz] = StressDueToSegs(1, 1, x, y, 0,...
                                              x1(1), x1(2), x1(3),...
                                              x2(1), x2(2), x2(3),...
                                              b(1),b(2),b(3),...
                                              a,mu,nu);
                                          
            T = [sxx.*n(1,1) + sxy.*n(1,2) + sxz.*n(1,3) ...
                 sxy.*n(1,1) + syy.*n(1,2) + syz.*n(1,3) ...
                 sxz.*n(1,1) + syz.*n(1,2) + szz.*n(1,3)];
            
            ft = wy*wx*T*LxLy*factor;
             
            fx1 = fx1 + ft*N1;
            fx2 = fx2 + ft*N2;
            fx3 = fx3 + ft*N3;
            fx4 = fx4 + ft*N4;
         end %for
     end %for
     ftot = fx1 + fx2 + fx3 + fx4;
     
    
%     gamma = gammat(gammat(:, idx) == val,:);
%     area = zeros(dim(1),1) + gamma(1,2);
%     gamma(1,2)
%     normal = gamma(1,3:5);
%     lseg = size(segments,1);
%     lgrid = dim(1);
%     p1x = segments(:,6);
%     p1y = segments(:,7);
%     p1z = segments(:,8);
%     p2x = segments(:,9);
%     p2y = segments(:,10);
%     p2z = segments(:,11);
%     bx = segments(:,3);
%     by = segments(:,4);
%     bz = segments(:,5);
%     x = xmid(:,1);
%     y = xmid(:,2);
%     z = xmid(:,3);
%     [sxx, syy, szz, sxy, syz, sxz] = StressDueToSegs(lgrid, lseg,...
%                                                   x, y, z,...
%                                                   p1x,p1y,p1z,...
%                                                   p2x,p2y,p2z,...
%                                                   bx,by,bz,...
%                                                   a,MU,NU); 
%     Tx = sxx*normal(1) + sxy*normal(2) + sxz*normal(3);
%     Ty = sxy*normal(1) + syy*normal(2) + syz*normal(3);
%     Tz = sxz*normal(1) + syz*normal(2) + szz*normal(3);
%     ATx = area.*Tx;
%     ATy = area.*Ty;
%     ATz = area.*Tz;
%     f_dln = [ATx,ATy,ATz];
end %function

function [x, w] = lgwt(N,a,b)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end