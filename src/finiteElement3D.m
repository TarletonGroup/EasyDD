function [vertices, B, xnodes, mno, nc, n, D, kg, w, h, d, mx, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, my, mz, mu, nu)
    %=========================================================================%
    % E Tarleton edmund.tarleton@materials.ox.ac.uk
    % 3D FEM code using linear 8 node element with 8 integration pts (2x2x2)
    % per element
    %
    % 4. ----- .3
    %  |\       |\
    %  | \      | \
    % 1. -\---- .2 \
    %   \  \     \  \
    %    \ 8. ----\- .7
    %     \ |      \ |
    %      \|       \|
    %      5. ----- .6
    %
    %
    % rectangulare domain.
    % (y going in to the page) note this is rotated about x axis w.r.t. local
    % (s1,s2,s3) system.
    % -------------------------------------------------
    %
    % ^z                   (mx,my)
    % |
    % |
    % |------>x-----------------------------------------
    %=========================================================================%

    fprintf('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.\n');

    % Simulation domain vertices. Required for surface remeshing.
    vertices = [0, 0, 0; ...
                dx, 0, 0; ...
                0, dy, 0; ...
                dx, dy, 0; ...
                0, 0, dz; ...
                dx, 0, dz; ...
                0, dy, dz; ...
                dx, dy, dz];
    % not that this code uses a constant element size in base and beam: (w,h,d)
    % also gammau is the left surface (eg baseLength=0)
    % loading = 1;

    w = dx / mx; % elements width
    h = dy / my; % element height
    d = dz / mz; % element depth

    mel = mx * my * mz;
    mno = (mx + 1) * (my + 1) * (mz + 1);

    ey = mu * (2 * (1 + nu));
    la = nu * ey / ((1 + nu) * (1 - 2 * nu));

    D = zeros(6, 6);

    for i = 1:3
        D(i, i) = la + 2 * mu;
    end

    for i = 4:6
        D(i, i) = mu;
    end

    D(1, 2) = la;
    D(1, 3) = la;
    D(2, 3) = la;

    % plane stress
    % D(1,1) = ey/(1-nu^2);
    % D(1,2) = ey*nu/(1-nu^2);
    % D(2,2) = D(1,1);
    % D(3,3) = ey/2/(1+nu);
    %
    D(2, 1) = D(1, 2);
    D(3, 1) = D(1, 3);
    D(3, 2) = D(2, 3);

    %
    %  ^s2
    %  |
    %  |
    %  |------>s1
    %
    %(s3 out of page)
    zv = 1 / sqrt(3);
    z = zeros(8, 3);
    z(1, 1:3) = [-zv, -zv, -zv];
    z(2, 1:3) = [zv, -zv, -zv];
    z(3, 1:3) = [zv, zv, -zv];
    z(4, 1:3) = [-zv, zv, -zv];
    z(5, 1:3) = [-zv, -zv, zv];
    z(6, 1:3) = [zv, -zv, zv];
    z(7, 1:3) = [zv, zv, zv];
    z(8, 1:3) = [-zv, zv, zv];

    xnodes = zeros(mno, 3); %global nodal coordinates gn

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

    for k = 1:mz + 1

        for j = 1:my + 1% nodes in x direction

            for i = 1:mx + 1
                gn = i + (k - 1) * (mx + 1) + (j - 1) * (mx + 1) * (mz + 1); %global node #
                xnodes(gn, 1) = (i - 1) * w;
                xnodes(gn, 2) = (j - 1) * h;
                xnodes(gn, 3) = (k - 1) * d;
            end

        end

    end

    nc = zeros(mel, 4); %element connectivity
    % 4. ----- .3
    %  �\       �\
    %  � \      � \
    % 1. -\---- .2 \
    %   \  \     \  \
    %    \ 8. ----\- .7
    %     \ �      \ �
    %      \�       \�
    %      5. ----- .6

    for k = 1:mz

        for j = 1:my% i = row index, j = column index

            for i = 1:mx%global node #(element #,local node #)
                ge = i + (k - 1) * (mx) + (j - 1) * mx * mz; % element number
                nc(ge, 1) = i + (k - 1) * (mx + 1) + (mx + 1) * (mz + 1) * j;
                nc(ge, 2) = nc(ge, 1) + 1;
                nc(ge, 4) = i + k * (mx + 1) + (mx + 1) * (mz + 1) * j;
                nc(ge, 3) = nc(ge, 4) + 1;
                nc(ge, 5) = i + (k - 1) * (mx + 1) + (mx + 1) * (mz + 1) * (j - 1);
                nc(ge, 6) = nc(ge, 5) + 1;
                nc(ge, 8) = i + k * (mx + 1) + (mx + 1) * (mz + 1) * (j - 1);
                nc(ge, 7) = nc(ge, 8) + 1;
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

    n = zeros(8, 8);
    % n(q,a) is shape function Na(s1,s2,s3) evaluated at integration point q
    for q = 1:8
        n(q, 1) = 1/8 * (1 - z(q, 1)) * (1 - z(q, 2)) * (1 - z(q, 3));
        n(q, 2) = 1/8 * (1 + z(q, 1)) * (1 - z(q, 2)) * (1 - z(q, 3));
        n(q, 3) = 1/8 * (1 + z(q, 1)) * (1 + z(q, 2)) * (1 - z(q, 3));
        n(q, 4) = 1/8 * (1 - z(q, 1)) * (1 + z(q, 2)) * (1 - z(q, 3));
        n(q, 5) = 1/8 * (1 - z(q, 1)) * (1 - z(q, 2)) * (1 + z(q, 3));
        n(q, 6) = 1/8 * (1 + z(q, 1)) * (1 - z(q, 2)) * (1 + z(q, 3));
        n(q, 7) = 1/8 * (1 + z(q, 1)) * (1 + z(q, 2)) * (1 + z(q, 3));
        n(q, 8) = 1/8 * (1 - z(q, 1)) * (1 + z(q, 2)) * (1 + z(q, 3));
    end

    ns = zeros(8, 8, 3);
    % derivative of shape function Na(s1,s2,s3)  w.r.t. sj evaluated at
    %integration point q  ns(q,a,j) = (dNa/dsj)s=zq
    for q = 1:8
        ns(q, 1, 1) = -1/8 * (1 - z(q, 2)) * (1 - z(q, 3));
        ns(q, 1, 2) = -1/8 * (1 - z(q, 1)) * (1 - z(q, 3));
        ns(q, 1, 3) = -1/8 * (1 - z(q, 1)) * (1 - z(q, 2));
        ns(q, 2, 1) = 1/8 * (1 - z(q, 2)) * (1 - z(q, 3));
        ns(q, 2, 2) = -1/8 * (1 + z(q, 1)) * (1 - z(q, 3));
        ns(q, 2, 3) = -1/8 * (1 + z(q, 1)) * (1 - z(q, 2));
        ns(q, 3, 1) = 1/8 * (1 + z(q, 2)) * (1 - z(q, 3));
        ns(q, 3, 2) = 1/8 * (1 + z(q, 1)) * (1 - z(q, 3));
        ns(q, 3, 3) = -1/8 * (1 + z(q, 1)) * (1 + z(q, 2));
        ns(q, 4, 1) = -1/8 * (1 + z(q, 2)) * (1 - z(q, 3));
        ns(q, 4, 2) = 1/8 * (1 - z(q, 1)) * (1 - z(q, 3));
        ns(q, 4, 3) = -1/8 * (1 - z(q, 1)) * (1 + z(q, 2));

        ns(q, 5, 1) = -1/8 * (1 - z(q, 2)) * (1 + z(q, 3));
        ns(q, 5, 2) = -1/8 * (1 - z(q, 1)) * (1 + z(q, 3));
        ns(q, 5, 3) = 1/8 * (1 - z(q, 1)) * (1 - z(q, 2));
        ns(q, 6, 1) = 1/8 * (1 - z(q, 2)) * (1 + z(q, 3));
        ns(q, 6, 2) = -1/8 * (1 + z(q, 1)) * (1 + z(q, 3));
        ns(q, 6, 3) = 1/8 * (1 + z(q, 1)) * (1 - z(q, 2));
        ns(q, 7, 1) = 1/8 * (1 + z(q, 2)) * (1 + z(q, 3));
        ns(q, 7, 2) = 1/8 * (1 + z(q, 1)) * (1 + z(q, 3));
        ns(q, 7, 3) = 1/8 * (1 + z(q, 1)) * (1 + z(q, 2));
        ns(q, 8, 1) = -1/8 * (1 + z(q, 2)) * (1 + z(q, 3));
        ns(q, 8, 2) = 1/8 * (1 - z(q, 1)) * (1 + z(q, 3));
        ns(q, 8, 3) = 1/8 * (1 - z(q, 1)) * (1 + z(q, 2));
    end

    pm1 = [-1 1 1 -1 -1 1 1 -1];
    pm2 = [-1 -1 1 1 -1 -1 1 1];
    pm3 = [-1 -1 -1 -1 1 1 1 1];

    % pm1 =[-1  1  1 -1 -1  1  1 -1]
    % pm2 =[-1 -1  1  1 -1 -1  1  1]
    % pm3 =[-1 -1 -1 -1  1  1  1  1]
    dNds = zeros(8, 8, 3);
    % Na(s1,s2,s3) = 1/8 *(1+pm1(a)s1)(1+pm2(a)s2)(1+pm3(a)s3)
    for q = 1:8

        for a = 1:8
            dNds(q, a, 1) = 1/8 * pm1(a) * (1 + pm2(a) * z(q, 2)) * (1 + pm3(a) * z(q, 3));
            dNds(q, a, 2) = 1/8 * (1 + pm1(a) * z(q, 1)) * pm2(a) * (1 + pm3(a) * z(q, 3));
            dNds(q, a, 3) = 1/8 * (1 + pm1(a) * z(q, 1)) * (1 + pm2(a) * z(q, 2)) * pm3(a);
        end

    end

    if any(dNds ~= ns)
        fprintf('error\n')
        pause
    end

    J = zeros(3, 3, mel, 8); % 3Dx3D, mel elements, 8 quad pts/element
    detJ = zeros(mel, 8); % det(J) at mel elements, 9 quad points/element
    nx = zeros(mel, 8, 8, 3); % derivative of shape functions w.r.t global x,y
    B = zeros(6, 24, mel, 8); % (# strain components, # dof/element, # elements, int pts/element)

    for p = 1:mel% all elements

        for q = 1:8% integration points per element

            for i = 1:3% DOF

                for j = 1:3% DOF
                    J(i, j, p, q) = ns(q, 1:8, j) * xnodes(nc(p, 1:8), i);
                    % implied sum over a: local shape function number
                    % sum over a of dNa(s1,s2,s3)/dsj at S=z(q) *xi(local node a of element p)
                    %Jij= dxi/dsj evaluated at element p integration point q
                end

            end

            detJ(p, q) = det(J(:, :, p, q)); % det(J) evaluated at element p int point q

            invJ = inv(J(:, :, p, q)); % invJ_ij = dsi/dxj

            for a = 1:8

                for j = 1:3

                    for i = 1:3
                        nx(p, q, a, j) = nx(p, q, a, j) + ns(q, a, i) * invJ(i, j);
                    end

                end

                B(1, (a - 1) * 3 + 1, p, q) = nx(p, q, a, 1);
                B(1, (a - 1) * 3 + 2, p, q) = 0;
                B(1, (a - 1) * 3 + 3, p, q) = 0;

                B(2, (a - 1) * 3 + 1, p, q) = 0;
                B(2, (a - 1) * 3 + 2, p, q) = nx(p, q, a, 2);
                B(2, (a - 1) * 3 + 3, p, q) = 0;

                B(3, (a - 1) * 3 + 1, p, q) = 0;
                B(3, (a - 1) * 3 + 2, p, q) = 0;
                B(3, (a - 1) * 3 + 3, p, q) = nx(p, q, a, 3);

                B(4, (a - 1) * 3 + 1, p, q) = nx(p, q, a, 2);
                B(4, (a - 1) * 3 + 2, p, q) = nx(p, q, a, 1);
                B(4, (a - 1) * 3 + 3, p, q) = 0;

                B(5, (a - 1) * 3 + 1, p, q) = nx(p, q, a, 3);
                B(5, (a - 1) * 3 + 2, p, q) = 0;
                B(5, (a - 1) * 3 + 3, p, q) = nx(p, q, a, 1);

                B(6, (a - 1) * 3 + 1, p, q) = 0;
                B(6, (a - 1) * 3 + 2, p, q) = nx(p, q, a, 3);
                B(6, (a - 1) * 3 + 3, p, q) = nx(p, q, a, 2);

            end

        end

    end

    fprintf('local K...\n');
    ke = zeros(24, 24, mel); %local stiffness matrix

    for p = 1:mel%all elements

        for q = 1:8% int points per element
            ke(:, :, p) = ke(:, :, p) + B(:, :, p, q)' * D * B(:, :, p, q) * detJ(p, q);
        end

    end

    % ensure ke is symmetric eg remove any very small entries due to numerical
    % error.
    fprintf('global K...\n');

    %% kg Haiyang
    tic
    a = 1:8; %local node numbers
    dofLocal = [3 * (a - 1) + 1, 3 * (a - 1) + 2, 3 * (a - 1) + 3];
    ntriplets = 3 * mno;
    I = zeros(ntriplets, 1);
    J = zeros(ntriplets, 1);
    X = zeros(ntriplets, 1);
    ntriplets = 0;

    for p = 1:mel
        gn = nc(p, a); % global node numbers
        dof(1:24) = [3 * (gn - 1) + 1, 3 * (gn - 1) + 2, 3 * (gn - 1) + 3]; % global degree of freedom

        for i = 1:24

            for j = 1:24
                ntriplets = ntriplets + 1;
                I(ntriplets) = dof(i);
                J(ntriplets) = dof(j);
                X(ntriplets) = ke(dofLocal(i), dofLocal(j), p);
            end

        end

        %  K{N}(dof,dof)= K{N}(dof,dof)+Ke(dofLocal,dofLocal,p); %the (full) global stiffness matrix

    end

    kg = sparse(I, J, X, 3 * mno, 3 * mno); %the (full) global stiffness matrix

    toc
end
