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
x1i = [0 0 0];
x2i = [0 0 1];
b  = [1 0 0];

% Gauss-Legrende quadrature seems to be numerically unstable w.r.t. 
% the number of quadrature points for this problem. I've tried with other 
% dislocation lines and they have different behaviours.



% There seems to be a factor of 3 missing from the analytical solutions.

n_q  = [1; 2; 10; 11; 100; 101; 500; 501; 1000; 1001; 1500; 1501];
dist = [0; 1E-6; 1E-5; 1E-4; 1E-3; 1E-2; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6];
save("parameters.mat")
% Cycle gauss points
for i = 1: size(n_q)
    [q, w] = lgwt(n_q(i), -1, 1);
    quad = [q, w];
    for j = 1: size(dist)
        x1 = x1i;
        x1(1,3) = x1(1,3) + dist(j);
        x2i = x1;
        for k = 1: size(dist)
            x2(1,3) = x2i(1,3) + dist(k);
            if dist(k) == 0
                x2(1,3) = x2(1,3) + a;
            end %if
            [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad);
            if i == 1
                [~, ~, ~, ~, ftota] = ...
                    nodal_surface_force_linear_rectangle_mex(x1, x2, ...
                                    x3, x4, x5, x6, b, mu, nu, a, 1, 1);
                ftota = ftota';
                save(sprintf("a_%f_%f.mat", x1(1,3), x2(1,3)), "ftota")
            end %if
            save(sprintf("n_%d_%f_%f.mat", n_q(i), x1(1,3), x2(1,3)), "ftotn")
        end %for
    end %for
end %for
% cntr = 0;
% for i = 1: size(dist)
%     x1(1, 3) = x1i(1, 3) + dist(i);
%     for j = 1: size(dist)
%         x2(1, 3) = x2i(1, 3) + dist(j);
%         if (x1(1,3) == x2(1,3))
%             continue
%         end %if
%         [fx3a, fx4a, fx5a, fx6a, ftota] = ...
%             nodal_surface_force_linear_rectangle_mex(x1, x2, ...
%                             x3, x4, x5, x6, b, mu, nu, a, 1, 1);
%             fx3a = fx3a';
%             fx4a = fx4a';
%             fx5a = fx5a';
%             fx6a = fx6a';
%             ftota = ftota';
%         for k = 1: size(n_q)
%             [q, w] = lgwt(n_q(k), -1, 1);
%             quad = [q, w];
%             [fx3n, fx4n, fx5n, fx6n, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad);
%             
%             ex3 = error(fx3a, fx3n);
%             ex4 = error(fx4a, fx4n);
%             ex5 = error(fx5a, fx5n);
%             ex6 = error(fx6a, fx6n);
%             extot = error(ftota, ftotn);
%             
%             cntr = cntr + 1;
%             save(sprintf("%d", cntr));
%         end %for
%     end %for
% end %for

function err = error(analytical, numerical)
    diff = numerical - analytical;
    diff(diff < 1E-10) = 0;
    err = diff./analytical;
    err(isnan(err)) = 0;
end %function

function [N1, N2, N3, N4] = shape_func(s1, s2)
    N1 = 0.25*(1-s1)*(1-s2);
    N2 = 0.25*(1+s1)*(1-s2);
    N3 = 0.25*(1+s1)*(1+s2);
    N4 = 0.25*(1-s1)*(1+s2);
end %function

function [fx1, fx2, fx3, fx4, ftot] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad)
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
     
     n_q = size(quad, 1);
     
     for i = 1: n_q
         y  = quad(i, 1);
         wy = quad(i, 2);
         for j = 1: n_q
             x  = quad(j, 1);
             wx = quad(j, 2);
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