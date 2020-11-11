% n-dimensional Gauss-Legendre in orthogonal coordinates.
% coords = [-3,5; -1,3; -2,2];
% dx = coords(:,2) - coords(:,1);
% ds = limits(:,2)-limits(:,1);
% InvJ = prod(ds./dx);

%%
limits = [-1, 1; -1, 1; -1, 1];
ndim = 3;
npoints = 2;
ivol = 1 / (prod(limits(:, 2) - limits(:, 1)));
[q, w] = gen_gauss_quad(ndim, npoints, limits);
v = threed_gauss_mesh(q, w, npoints);
% zv = 1/sqrt(3);
% z = zeros(8,3);
% z(1,1:3) = [-zv,-zv,-zv];
% z(2,1:3) = [ zv,-zv,-zv];
% z(3,1:3) = [ zv, zv,-zv];
% z(4,1:3) = [-zv, zv,-zv];
% z(5,1:3) = [-zv,-zv, zv];
% z(6,1:3) = [ zv,-zv, zv];
% z(7,1:3) = [ zv, zv, zv];
% z(8,1:3) = [-zv, zv, zv];
N = threed_shape_func(v(:, 1:3));

n_dln = 1;
tpoints = size(v, 1);
% x1 = [-0.75, 0, 0];
% x2 = [ 0.25, 0, 0];
x1 = [-0.75 0, 0];
x2 = [0.25, 0, 0];
b = [0, 1, 0];
a = 0.2;
mu = 1;
nu = 0.29;

fx_fem = zeros(8, 3);

[sxx, syy, szz, ...
        sxy, syz, sxz] = StressDueToSegs(tpoints, n_dln, v(:, 1), v(:, 2), v(:, 3), ...
    x1(:, 1), x1(:, 2), x1(:, 3), ...
    x2(:, 1), x2(:, 2), x2(:, 3), ...
    b(:, 1), b(:, 2), b(:, 3), ...
    a, mu, nu);

sigmap = [sxx, syy, szz, ...
        sxy, syz, sxz];

sigmav = zeros(n_dln, size(sigmap, 2));

% \sigma_v = v^{-1} \int_v \sigma(\vec{x}) dv
% \approx v^{-1} \sum_i=1^I \sum_j=1^J \sum_k=1^K wi wj wk \sigma(xi,xj,xk)
sigmav(:, :) = sum(v(:, 4) .* sigmap(:, :)) * ivol;

% Fn = \int_v \bm{\sigmav} \cdot \nabla N dv
% \approx \sum_i=1^I \sum_j=1^J \sum_k=1^K wi wj wk \bm{\sigmav} \cdot
% \nabla N(xi,yj,zk)
svx = [sigmav(:, 1), sigmav(:, 4), sigmav(:, 6)];
svy = [sigmav(:, 4), sigmav(:, 2), sigmav(:, 5)];
svz = [sigmav(:, 6), sigmav(:, 5), sigmav(:, 3)];

npoints = 1;
[q, w] = gen_gauss_quad(ndim, npoints, limits);
v = threed_gauss_mesh(q, w, npoints);

ns = dthreed_shape_func_dxi(v(:, 1:3));

fx_fem = zeros(8, 3);

for j = 1:8

    for i = 1:3
        fx_fem(j, 1) = fx_fem(j, 1) + sum(v(:, 4) .* ns(:, j, i)) .* svx(:, i);
        fx_fem(j, 2) = fx_fem(j, 2) + sum(v(:, 4) .* ns(:, j, i)) .* svy(:, i);
        fx_fem(j, 3) = fx_fem(j, 3) + sum(v(:, 4) .* ns(:, j, i)) .* svz(:, i);
    end

end

fx_fem
sum(fx_fem)

% fx1(1,1) = fx1(1,1) + sum(v(:,4).*ns(:,1,1).*svx(:,1));
% fx1(1,1) = fx1(1,1) + sum(v(:,4).*ns(:,1,2).*svx(:,2));
% fx1(1,1) = fx1(1,1) + sum(v(:,4).*ns(:,1,3).*svx(:,3));
%
% fx1(1,2) = fx1(1,2) + sum(v(:,4).*ns(:,1,1).*svy(:,1));
% fx1(1,2) = fx1(1,2) + sum(v(:,4).*ns(:,1,2).*svy(:,2));
% fx1(1,2) = fx1(1,2) + sum(v(:,4).*ns(:,1,3).*svy(:,3));
%
% fx1(1,3) = fx1(1,3) + sum(v(:,4).*ns(:,1,1).*svz(:,1));
% fx1(1,3) = fx1(1,3) + sum(v(:,4).*ns(:,1,2).*svz(:,2));
% fx1(1,3) = fx1(1,3) + sum(v(:,4).*ns(:,1,3).*svz(:,3));

% fx1_fem = zeros(8,3);
% for i = 1:3 % x, y, z
%     for j = 1:tpoints
%         for k = 1:8 %nodes
%             fx1_fem(k,1) = fx1_fem(k,1) + v(j,4)*ns(j,k,i)*svx(1,i);
%             fx1_fem(k,2) = fx1_fem(k,2) + v(j,4)*ns(j,k,i)*svy(1,i);
%             fx1_fem(k,3) = fx1_fem(k,3) + v(j,4)*ns(j,k,i)*svz(1,i);
%         end
%     end
% end

%%
z = zeros(8, 3);
z(1, 1:3) = [-1, -1, -1];
z(2, 1:3) = [1, -1, -1];
z(3, 1:3) = [1, 1, -1];
z(4, 1:3) = [-1, 1, -1];
z(5, 1:3) = [-1, -1, 1];
z(6, 1:3) = [1, -1, 1];
z(7, 1:3) = [1, 1, 1];
z(8, 1:3) = [-1, 1, 1];

z = threed_gauss_mesh([1, 1, 1; -1, -1, -1], [1, 1, 1; 1, 1, 1], 2)

% x3 = [z(1,:)'; z(6,:)'; z(2,:)'; z(5,:)'; z(4,:)'; z(5,:)'];
% x4 = [z(2,:)'; z(5,:)'; z(6,:)'; z(1,:)'; z(3,:)'; z(6,:)'];
% x5 = [z(4,:)'; z(7,:)'; z(3,:)'; z(8,:)'; z(8,:)'; z(1,:)'];
% x6 = [z(3,:)'; z(8,:)'; z(7,:)'; z(4,:)'; z(7,:)'; z(2,:)'];

x3 = [z(1, :); z(6, :); z(2, :); z(5, :); z(4, :); z(5, :)];
x4 = [z(2, :); z(5, :); z(6, :); z(1, :); z(3, :); z(6, :)];
x5 = [z(4, :); z(7, :); z(3, :); z(8, :); z(8, :); z(1, :)];
x6 = [z(3, :); z(8, :); z(7, :); z(4, :); z(7, :); z(2, :)];

n_se = 6;
n_dln = 1;
eps = 1e-6;
fx_n = zeros(n_se * 3, 4);
fx_se = zeros(n_se * 3);
[fx_n(:, 1), fx_n(:, 2), fx_n(:, 3), fx_n(:, 4), ...
        fx_se] = nodal_surface_force_linear_rectangle_mex(...
    x1', x2', ...
    reshape(x3', [], 1), reshape(x4', [], 1), ...
    reshape(x5', [], 1), reshape(x6', [], 1), ...
    b, mu, nu, a, n_se, n_dln, eps);
fx3 = reshape(fx_n(:, 1), 3, [])';
fx4 = reshape(fx_n(:, 2), 3, [])';
fx5 = reshape(fx_n(:, 3), 3, [])';
fx6 = reshape(fx_n(:, 4), 3, [])';

fx_trac = zeros(8, 3);

for i = 1:8
    idx3 = find(sum(x3 == z(i, :), 2) == 3);
    idx4 = find(sum(x4 == z(i, :), 2) == 3);
    idx5 = find(sum(x5 == z(i, :), 2) == 3);
    idx6 = find(sum(x6 == z(i, :), 2) == 3);
    fx_trac(i, :) = sum(fx3(idx3, :), 1) + sum(fx4(idx4, :), 1) + sum(fx5(idx5, :), 1) + sum(fx6(idx6, :), 1);
end

fx_trac
sum(fx_trac)

diff = fx_trac - fx_fem
sum(diff)

[points, weights] = gen_gauss_quad(2, 1, limits);
twod_shape_func(points);

function mesh = threed_gauss_mesh(points, weights, npoints)
    mesh = zeros(npoints^3, 3 + 1);

    idx = 0;
    cntr = 0;
    % Loop through z.
    for i = npoints:-1:1

        for j = npoints:-1:1

            if mod(cntr, 2) == 0

                for k = npoints:-1:1
                    idx = idx + 1;
                    n = npoints - mod(idx - 1, npoints);
                    mesh(idx, 1:3) = [points(n, 1), points(j, 2), points(i, 3)];
                    mesh(idx, 4) = weights(n, 1) * weights(j, 2) * weights(i, 3);
                end

            else

                for k = 1:npoints
                    idx = idx + 1;
                    n = k;
                    mesh(idx, 1:3) = [points(n, 1), points(j, 2), points(i, 3)];
                    mesh(idx, 4) = weights(n, 1) * weights(j, 2) * weights(i, 3);
                end

            end

            cntr = cntr + 1;
        end

    end

end

% 8 integrals, one per node (column).
% Each row represents the contribution of a given point to the node.
function n = threed_shape_func(s)
    n = zeros(size(s, 1), 8);
    n(:, 1) = 0.125 * (1 - s(:, 1)) .* (1 - s(:, 2)) .* (1 - s(:, 3));
    n(:, 2) = 0.125 * (1 + s(:, 1)) .* (1 - s(:, 2)) .* (1 - s(:, 3));
    n(:, 3) = 0.125 * (1 + s(:, 1)) .* (1 + s(:, 2)) .* (1 - s(:, 3));
    n(:, 4) = 0.125 * (1 - s(:, 1)) .* (1 + s(:, 2)) .* (1 - s(:, 3));
    n(:, 5) = 0.125 * (1 - s(:, 1)) .* (1 - s(:, 2)) .* (1 + s(:, 3));
    n(:, 6) = 0.125 * (1 + s(:, 1)) .* (1 - s(:, 2)) .* (1 + s(:, 3));
    n(:, 7) = 0.125 * (1 + s(:, 1)) .* (1 + s(:, 2)) .* (1 + s(:, 3));
    n(:, 8) = 0.125 * (1 - s(:, 1)) .* (1 + s(:, 2)) .* (1 + s(:, 3));
end %function

function ns = dthreed_shape_func_dxi(s)
    ns = zeros(size(s, 1), 8, 3);

    ns(:, 1, 1) = -0.125 * (1 - s(:, 2)) .* (1 - s(:, 3));
    ns(:, 1, 2) = -0.125 * (1 - s(:, 1)) .* (1 - s(:, 3));
    ns(:, 1, 3) = -0.125 * (1 - s(:, 1)) .* (1 - s(:, 2));

    ns(:, 2, 1) = 0.125 * (1 - s(:, 2)) .* (1 - s(:, 3));
    ns(:, 2, 2) = -0.125 * (1 + s(:, 1)) .* (1 - s(:, 3));
    ns(:, 2, 3) = -0.125 * (1 + s(:, 1)) .* (1 - s(:, 2));

    ns(:, 3, 1) = 0.125 * (1 + s(:, 2)) .* (1 - s(:, 3));
    ns(:, 3, 2) = 0.125 * (1 + s(:, 1)) .* (1 - s(:, 3));
    ns(:, 3, 3) = -0.125 * (1 + s(:, 1)) .* (1 + s(:, 2));

    ns(:, 4, 1) = -0.125 * (1 + s(:, 2)) .* (1 - s(:, 3));
    ns(:, 4, 2) = 0.125 * (1 - s(:, 1)) .* (1 - s(:, 3));
    ns(:, 4, 3) = -0.125 * (1 - s(:, 1)) .* (1 + s(:, 2));

    ns(:, 5, 1) = -0.125 * (1 - s(:, 2)) .* (1 + s(:, 3));
    ns(:, 5, 2) = -0.125 * (1 - s(:, 1)) .* (1 + s(:, 3));
    ns(:, 5, 3) = 0.125 * (1 - s(:, 1)) .* (1 - s(:, 2));

    ns(:, 6, 1) = 0.125 * (1 - s(:, 2)) .* (1 + s(:, 3));
    ns(:, 6, 2) = -0.125 * (1 + s(:, 1)) .* (1 + s(:, 3));
    ns(:, 6, 3) = 0.125 * (1 + s(:, 1)) .* (1 - s(:, 2));

    ns(:, 7, 1) = 0.125 * (1 + s(:, 2)) .* (1 + s(:, 3));
    ns(:, 7, 2) = 0.125 * (1 + s(:, 1)) .* (1 + s(:, 3));
    ns(:, 7, 3) = 0.125 * (1 + s(:, 1)) .* (1 + s(:, 2));

    ns(:, 8, 1) = -0.125 * (1 + s(:, 2)) .* (1 + s(:, 3));
    ns(:, 8, 2) = 0.125 * (1 - s(:, 1)) .* (1 + s(:, 3));
    ns(:, 8, 3) = 0.125 * (1 - s(:, 1)) .* (1 + s(:, 2));
end %function

function n = twod_shape_func(s)
    n = zeros(size(s, 1), 4);
    n(:, 1) = 0.25 * (1 - s(:, 1)) .* (1 - s(:, 2));
    n(:, 2) = 0.25 * (1 + s(:, 1)) .* (1 - s(:, 2));
    n(:, 3) = 0.25 * (1 + s(:, 1)) .* (1 + s(:, 2));
    n(:, 4) = 0.25 * (1 - s(:, 1)) .* (1 + s(:, 2));
end %function

function [points, weights] = gen_gauss_quad(n_coords, n_nodes, limits)
    % Weights and points for gauss quadrature for the number of coordinates
    % and nodes.
    weights = zeros(n_nodes, n_coords);
    points = zeros(n_nodes, n_coords);

    for i = 1:n_coords
        [points(:, i), weights(:, i)] = lgwt(n_nodes, limits(i, 1), limits(i, 2));
    end %for coord

end %function 3d_gausss_quad

function limits = gen_limits(limits)
    weights = zeros(n_nodes, n_coords);
    points = zeros(n_nodes, n_coords);

    for i = 1:n_coords
        [points(:, i), weights(:, i)] = lgwt(n_nodes, limits(i, 1), limits(i, 2));
    end %for coord

end %function gen_limits

function [x, w] = lgwt(N, a, b)

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
    N = N - 1;
    N1 = N + 1; N2 = N + 2;

    xu = linspace(-1, 1, N1)';

    % Initial guess
    y = cos((2 * (0:N)' + 1) * pi / (2 * N + 2)) + (0.27 / N1) * sin(pi * xu * N / N2);

    % Legendre-Gauss Vandermonde Matrix
    L = zeros(N1, N2);

    % Derivative of LGVM
    Lp = zeros(N1, N2);

    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method

    y0 = 2;

    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y - y0)) > eps
        L(:, 1) = 1;
        Lp(:, 1) = 0;

        L(:, 2) = y;
        Lp(:, 2) = 1;

        for k = 2:N1
            L(:, k + 1) = ((2 * k - 1) * y .* L(:, k) - (k - 1) * L(:, k - 1)) / k;
        end

        Lp = (N2) * (L(:, N1) - y .* L(:, N2)) ./ (1 - y.^2);

        y0 = y;
        y = y0 - L(:, N2) ./ Lp;

    end

    % Linear map from[-1,1] to [a,b]
    x = (a * (1 - y) + b * (1 + y)) / 2;

    % Compute the weights
    w = (b - a) ./ ((1 - y.^2) .* Lp.^2) * (N2 / N1)^2;

end
