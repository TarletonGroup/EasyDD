% This is badly coded as there is a lot of hard code but it is a proof of
% concept for an idealised case.
clear all;
close all;
mu = 1;
nu = 0.29;

%  Same numbering scheme sas the analytical solution.
%  x5------x6
%  |       |
%  |       |
%  x3------x4
%%
x3 = [-0.5E3 -0.5E3 0];
x4 = [ 0.5E3 -0.5E3 0];
x5 = [-0.5E3  0.5E3 0];
x6 = [ 0.5E3  0.5E3 0];
%
n_q  = [1; 2; 10; 11; 100; 101; 500; 501; 1000; 1001; 1500; 1501];
%%
% Dislocation passes over x3
x1i = [-0.5E3 -0.5E3 0];
x2i = [-0.5E3 -0.5E3 1];
b   = [1 0 0];
a   = 5*norm(b);
dist = a + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];

save("./worst_case/eperp_params_x3.mat")
for i = 1: size(n_q)
    [q, w] = lgwt(n_q(i), -0.5E3, 0.5E3);
    quad = [q, w];
    for j = 1: size(dist)
        x1 = x1i;
        x1(1,3) = x1(1,3) + dist(j);
        x2i = x1;
        for k = 1: size(dist)
            x2(1,3) = x2i(1,3) + dist(k);
            if x1(1,3) == x2(1,3)
                x2(1,3) = x2(1,3) + 1;
            end %if
            [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, 1);
            if i == 1
                [~, ~, ~, ~, ftota] = ...
                    nodal_surface_force_linear_rectangle_mex(x1, x2, ...
                                    x3, x4, x5, x6, b, mu, nu, a, 1, 1, 0.0);
                ftota = ftota';
                save(sprintf("./worst_case/aeperp_%f_%f_x3.mat", x1(1,3), x2(end,3)), "ftota")
            end %if
            save(sprintf("./worst_case/neperp_%d_%f_%f_x3.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
        end %for
%         sprintf("q = %d, dist = %f, err = %f, ", n_q(i), dist(j), abs(ftotn(2)-ftota(2))./ftota(2))
    end %for
end %for
%%
n_dln =1000;
x1i = zeros(n_dln,3);
x2i = zeros(n_dln,3);
b = zeros(n_dln,3);
% Dislocation passes over x3 and x6
x1i(:,1) = linspace(-0.5E6,0.5E6-(0.5E6+0.5E6)/(n_dln),n_dln);
x2i(:,1) = linspace(-0.5E6+(0.5E6+0.5E6)/(n_dln),0.5E6,n_dln);
x1i(:,2) = linspace(-0.5E6,0.5E6-(0.5E6+0.5E6)/(n_dln),n_dln);
x2i(:,2) = linspace(-0.5E6+(0.5E6+0.5E6)/(n_dln),0.5E6,n_dln);
b(:,3)   = 1;
a   = 5*norm(b(1,:));
dim = size(x1i,1)*size(x1i,2);
dist = a + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];

% save("./worst_case/epar_params_x3x6.mat")
for i = 1: size(n_q)
    [q, w] = lgwt(n_q(i), -0.5E3, 0.5E3);
    quad = [q, w];
    for j = 1: size(dist)
        x1 = x1i;
        x2 = x2i;
        x1(:,3) = x1(:,3) + dist(j);
        x2(:,3) = x2(:,3) + dist(j);
        [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, n_dln);
        if i == 1
            [~, ~, ~, ~, ftota] = ...
                nodal_surface_force_linear_rectangle_mex(reshape(x1',dim,1), reshape(x2',dim,1), ...
                                x3, x4, x5, x6, b, mu, nu, a, 1, n_dln, 1e-6);
            ftota = ftota';
            save(sprintf("./worst_case/aepar_%f_%f_x3x6.mat", x1(1,3), x2(end,3)), "ftota")
        end %if
        save(sprintf("./worst_case/nepar_%d_%f_%f_x3x6.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
    end %for
end %for

%%
% x1i = [-0.5E9 0 0];
% x2i = [ 0.5E9 0 0];
% b   = [1 0 0];
% a   = 5*norm(b);
% dist = a + [0; 1E-6; 1E-5; 1E-4; 1E-3; 1E-2; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
% n_dln = 10000;
% x1i = zeros(n_dln,3);
% x2i = zeros(n_dln,3);
% b = zeros(n_dln,3);
% x1i(:,1) = linspace(-0.5E6,0.5E6-(0.5E6+0.5E6)/(n_dln),n_dln);
% x2i(:,1) = linspace(-0.5E6+(0.5E6+0.5E6)/(n_dln),0.5E6,n_dln);
% b(:,1)   = 1;
% a   = 5*norm(b(1,:));
% dim = size(x1i,1)*size(x1i,2);
% dist = a + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
%
% save("spar_params.mat")
% for i = 1: size(n_q)
%     [q, w] = lgwt(n_q(i), -0.5E3, 0.5E3);
%     quad = [q, w];
%     for j = 1: size(dist)
%         x1 = x1i;
%         x2 = x2i;
%         x1(:,3) = x1(:,3) + dist(j);
%         x2(:,3) = x2(:,3) + dist(j);
%         [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, n_dln);
%         if i == 1
%             [~, ~, ~, ~, ftota] = ...
%                 nodal_surface_force_linear_rectangle_mex(reshape(x1',dim,1), reshape(x2',dim,1), ...
%                                 x3, x4, x5, x6, b, mu, nu, a, 1, n_dln, 1e-6);
%             ftota = ftota';
%             save(sprintf("aspar_%f_%f.mat", x1(1,3), x2(end,3)), "ftota")
%         end %if
%         save(sprintf("nspar_%d_%f_%f.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
%     end %for
% end %for
%%

% %%
% x3 = [-0.5E3 -0.5E3 0];
% x4 = [ 0.5E3 -0.5E3 0];
% x5 = [-0.5E3  0.5E3 0];
% x6 = [ 0.5E3  0.5E3 0];
%
% n_q  = [1; 2; 10; 11; 100; 101; 500; 501; 1000; 1001; 1500; 1501];
%
% % There seems to be a factor of 3 missing from the analytical solutions.
% % x1i = [0.5E3 0.5E3 0];
% % x2i = [0.5E3 0.5E3 1];
% % b   = [1 0 0];
% % a   = 5*norm(b);
% % dist = a + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
% %
% % save("aeperp_params.mat")
% % for i = 1: size(n_q)
% %     [q, w] = lgwt(n_q(i), -0.5E3, 0.5E3);
% %     quad = [q, w];
% %     for j = 1: size(dist)
% %         x1 = x1i;
% %         x1(1,3) = x1(1,3) + dist(j);
% %         x2i = x1;
% %         for k = 1: size(dist)
% %             x2(1,3) = x2i(1,3) + dist(k);
% %             if x1(1,3) == x2(1,3)
% %                 x2(1,3) = x2(1,3) + 1;
% %             end %if
% %             [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, 1);
% %             if i == 1
% %                 [~, ~, ~, ~, ftota] = ...
% %                     nodal_surface_force_linear_rectangle_mex(x1, x2, ...
% %                                     x3, x4, x5, x6, b, mu, nu, a, 1, 1);
% %                 ftota = ftota';
% %                 save(sprintf("aaeperp_%f_%f.mat", x1(1,3), x2(end,3)), "ftota")
% %             end %if
% %             save(sprintf("aneperp_%d_%f_%f.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
% %         end %for
% %     end %for
% % end %for
%
% n_dln = 100;
% x1i = zeros(n_dln,3);
% x2i = zeros(n_dln,3);
% b = zeros(n_dln,3);
% x1i(:,1) = linspace(-0.5E9,0.5E9-(0.5E9+0.5E9)/(n_dln),n_dln);
% x2i(:,1) = linspace(-0.5E9+(0.5E9+0.5E9)/(n_dln),0.5E9,n_dln);
% x1i(:,2) = 0.5E3;
% x2i(:,2) = 0.5E3;
% b(:,2)   = 1;
% a   = 5*norm(b(1,:));
% dim = size(x1i,1)*size(x1i,2);
% dist = a;% + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
%
% % save("aepar_params.mat")
% for i = 1: size(n_q)
%     [q, w] = lgwt(n_q(i), -0.5E3, 0.E3);
%     quad = [q, w];
%     for j = 1: size(dist)
%         x1 = x1i;
%         x2 = x2i;
%         x1(:,3) = x1(:,3) + dist(j);
%         x2(:,3) = x2(:,3) + dist(j);
%         for k = 1: size(dist)
%             [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, n_dln);
%             if i == 1
%                 [~, ~, ~, ~, ftota] = ...
%                     nodal_surface_force_linear_rectangle_mex(reshape(x1',dim,1), reshape(x2',dim,1), ...
%                                     x3, x4, x5, x6, b, mu, nu, a, 1, n_dln);
%                 ftota = ftota';
% %                 save(sprintf("aaepar_%f_%f.mat", x1(1,3), x2(end,3)), "ftota")
%             end %if
% %             save(sprintf("anepar_%d_%f_%f.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
%         end %for
%     end %for
% end %for
%
% % x1i = [-0.5E9 0 0];
% % x2i = [ 0.5E9 0 0];
% % b   = [1 0 0];
% % a   = 5*norm(b);
% % dist = a + [0; 1E-6; 1E-5; 1E-4; 1E-3; 1E-2; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
% n_dln = 100;
% x1i = zeros(n_dln,3);
% x2i = zeros(n_dln,3);
% b = zeros(n_dln,3);
% x1i(:,1) = linspace(-0.5E9,0.5E9-(0.5E9+0.5E9)/(n_dln),n_dln);
% x2i(:,1) = linspace(-0.5E9+(0.5E9+0.5E9)/(n_dln),0.5E9,n_dln);
% x1i(:,2) = 0.5E3;
% x2i(:,2) = 0.5E3;
% b(:,1)   = 1;
% a   = 5*norm(b(1,:));
% dim = size(x1i,1)*size(x1i,2);
% dist = a;% + [0; 1E-1; 1E0; 1E1; 1E2; 1E3; 1E4; 1E5; 1E6; 1E7; 1E8; 1E9];
%
% % save("aaspar_params.mat")
% for i = 1: size(n_q)
%     [q, w] = lgwt(n_q(i), -0.5E3, 0.E3);
%     quad = [q, w];
%     for j = 1: size(dist)
%         x1 = x1i;
%         x2 = x2i;
%         x1(:,3) = x1(:,3) + dist(j);
%         x2(:,3) = x2(:,3) + dist(j);
%         for k = 1: size(dist)
%             [~, ~, ~, ~, ftotn] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, 1, n_dln);
%             if i == 1
%                 [~, ~, ~, ~, ftota] = ...
%                     nodal_surface_force_linear_rectangle_mex(reshape(x1',dim,1), reshape(x2',dim,1), ...
%                                     x3, x4, x5, x6, b, mu, nu, a, 1, n_dln);
%                 ftota = ftota';
% %                 save(sprintf("aaspar_%f_%f.mat", x1(1,3), x2(end,3)), "ftota")
%             end %if
% %             save(sprintf("anspar_%d_%f_%f.mat", n_q(i), x1(1,3), x2(end,3)), "ftotn")
%         end %for
%     end %for
% end %for
% %%

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

function [fx1, fx2, fx3, fx4, ftot] = f_num(x1, x2, x3, x4, x5, x6, b, mu, nu, a, quad, n_se, n_dln)
     n = cross(x4-x3, x5-x3);
     n = n/norm(n);

     Lx = (x4(1)-x3(1));
     Ly = (x5(2)-x3(2));
     LxLy = Lx*Ly;

     factor = 1/LxLy;
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
              sxy, syz, sxz] = StressDueToSegs(n_se, n_dln, x, y, 0,...
                                              x1(:,1), x1(:,2), x1(:,3),...
                                              x2(:,1), x2(:,2), x2(:,3),...
                                              b(:,1),b(:,2),b(:,3),...
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
