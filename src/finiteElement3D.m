function [S, vertices, B, xnodes, mno, nc, n, D, kg, w, h, d, mx, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, my, mz, mu, nu)
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
    
    % @dhortela made this, i'm backporting.
    globnidx = zeros(mx + 1, my + 1, mz + 1);
    for k = 1:mz + 1

        for j = 1:my + 1% nodes in x direction

            for i = 1:mx + 1
                gn = i + (k - 1) * (mx + 1) + (j - 1) * (mx + 1) * (mz + 1); %global node #
                xnodes(gn, 1) = (i - 1) * w;
                xnodes(gn, 2) = (j - 1) * h;
                xnodes(gn, 3) = (k - 1) * d;
                globnidx(i,j,k) = gn;
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
    
    % Backporting @dhortela's node set.
    S = surfaceCuboid(w, h, d, mx, my, mz, globnidx, [1;0;0], [0;1;0], [0;0;1]);

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
    %%
    % These signs were used in hatStress, they're wrong, they error out with
    % the explicit equations. This is why we need unit tests!
    % I don't know why these are used there.
    %     pm1 = [-1 1 1 -1 -1 1 1 -1];
    %     pm2 = [1 1 1 1 -1 -1 -1 -1];
    %     pm3 = [-1 -1 1 1 -1 -1 1 1];
    %%
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

    J = zeros(3, 3, 8); % 3Dx3D, mel elements, 8 quad pts/element
    detJ = zeros(8, 1); % det(J) at mel elements, 9 quad points/element
    nx = zeros(8, 8, 3); % derivative of shape functions w.r.t global x,y
    B = zeros(6, 24, 8); % (# strain components, # dof/element, # elements, int pts/element)

    for q = 1:8% integration points per element

        for i = 1:3% DOF

            for j = 1:3% DOF
                J(i, j, q) = ns(q, 1:8, j) * xnodes(nc(1, 1:8), i);
                % implied sum over a: local shape function number
                % sum over a of dNa(s1,s2,s3)/dsj at S=z(q) *xi(local node a of element p)
                %Jij= dxi/dsj evaluated at element p integration point q
            end

        end

        detJ(q) = det(J(:, :, q)); % det(J) evaluated at element p int point q

        invJ = inv(J(:, :, q)); % invJ_ij = dsi/dxj

        for a = 1:8

            for j = 1:3

                for i = 1:3
                    nx(q, a, j) = nx(q, a, j) + ns(q, a, i) * invJ(i, j);
                end

            end

            B(1, (a - 1) * 3 + 1, q) = nx(q, a, 1);
            B(1, (a - 1) * 3 + 2, q) = 0;
            B(1, (a - 1) * 3 + 3, q) = 0;

            B(2, (a - 1) * 3 + 1, q) = 0;
            B(2, (a - 1) * 3 + 2, q) = nx(q, a, 2);
            B(2, (a - 1) * 3 + 3, q) = 0;

            B(3, (a - 1) * 3 + 1, q) = 0;
            B(3, (a - 1) * 3 + 2, q) = 0;
            B(3, (a - 1) * 3 + 3, q) = nx(q, a, 3);

            B(4, (a - 1) * 3 + 1, q) = nx(q, a, 2);
            B(4, (a - 1) * 3 + 2, q) = nx(q, a, 1);
            B(4, (a - 1) * 3 + 3, q) = 0;

            B(5, (a - 1) * 3 + 1, q) = nx(q, a, 3);
            B(5, (a - 1) * 3 + 2, q) = 0;
            B(5, (a - 1) * 3 + 3, q) = nx(q, a, 1);

            B(6, (a - 1) * 3 + 1, q) = 0;
            B(6, (a - 1) * 3 + 2, q) = nx(q, a, 3);
            B(6, (a - 1) * 3 + 3, q) = nx(q, a, 2);

        end

    end

    fprintf('local K...\n');
    ke = zeros(24, 24); %local stiffness matrix

    for q = 1:8% int points per element
        ke(:, :) = ke(:, :) + B(:, :, q)' * D * B(:, :, q) * detJ(q);
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
                X(ntriplets) = ke(dofLocal(i), dofLocal(j));
            end

        end

        %  K{N}(dof,dof)= K{N}(dof,dof)+Ke(dofLocal,dofLocal,p); %the (full) global stiffness matrix

    end

    kg = sparse(I, J, X, 3 * mno, 3 * mno); %the (full) global stiffness matrix

    toc
end


function [S] = surfaceCuboid(w, h, d, mx, my, mz, globnidx, nx, ny, nz)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Defines the geometry of the simulation domain surfaces.
    % This function specifically defines a cantilever (cuboid) with
    % orthorhombic axes.
    %===============================================================%
    
    %% Extract data from structures
    
%     % FEM:
%     w = FEM.w; h = FEM.h; d = FEM.d; % Finite element dimensions
%     mx = FEM.mx; my = FEM.my; mz = FEM.mz; % # FEs per dimension
%     globnidx = FEM.globnidx; % Matrix of global FE node indices
%     nx = FEM.normals(1,:); ny = FEM.normals(2,:); nz = FEM.normals(3,:);
    
    %% Define surface sets
    
    S = struct; % S contains all surface sets
    % S.s = [global node idx, nodal area, outward normal]
    % Each surface set corresponds to a subset of the FE nodes.
    % Nodes are never shared between surface sets.
    
    % Faces do not contain edges nor corners.
    % Edges do not contain corners nor faces.
    
    %%% Faces
    
    % Constant x:
    S.left = zeros((my-1)*(mz-1), 5);
    S.right = zeros((my-1)*(mz-1), 5);
    
    % Constant y:
    S.front = zeros((mx-1)*(mz-1), 5);
    S.back = zeros((mx-1)*(mz-1), 5);
    
    % Constant z:
    S.top = zeros((mx-1)*(my-1), 5);
    S.bot = zeros((mx-1)*(my-1), 5);
    
    %%% Edges
    
    % Top edges:
    S.topleft = zeros(my-1, 5); % Top left edge
    S.topright = zeros(my-1, 5); % Top right edge
    S.topfront = zeros(mx-1, 5); % Top front edge (& front top edge)
    S.topback = zeros(mx-1, 5); % Top back edge (& back top edge)
    
    % Bottom edges:
    S.botleft = zeros(my-1, 5);
    S.botright = zeros(my-1, 5);
    S.botfront = zeros(mx-1, 5); % Bot front edge (& front bot edge)
    S.botback = zeros(mx-1, 5); % Bot front edge (& front bot edge)
    
    % Front edges:
    S.frontleft = zeros(mz-1, 5);
    S.frontright = zeros(mz-1, 5);
    % Back edges:
    S.backleft = zeros(mz-1, 5);
    S.backright = zeros(mz-1, 5);
    
    %%% Corners
    
    S.corners = zeros(8,5);
    
    % Corner indexing:
%     vertices = [0, 0, 0; ...
%                 dx, 0, 0; ...
%                 0, dy, 0; ...
%                 dx, dy, 0; ...
%                 0, 0, dz; ...
%                 dx, 0, dz; ...
%                 0, dy, dz; ...
%                 dx, dy, dz];
    
    %% Fill surface sets %%
    
    
    %% Faces
    
    % Constant x
    fxA = h*d;
    S.left(:, 2) = fxA;%S.left(:, 3:5) = -nx;
    S.right(:, 2) = fxA;%S.right(:, 3:5) = nx;
    
    % Constant y
    fyA = w*d;
    S.front(:, 2) = fyA;%S.front(:, 3:5) = -ny;
    S.back(:, 2) = fyA;%S.back(:, 3:5) = ny;
    
    % Constant z
    fzA = w*h;
    S.bot(:, 2) = fzA;%S.bot(:, 3:5) = -nz;
    S.top(:, 2) = fzA;%S.top(:, 3:5) = nz;
    
    %%% Indexing
    
    m = 1;
    % Constant x
    for j = 2:my % Nodes in y direction
        for k = 2:mz % Nodes in z direction
            
            % Left face (x = 0)
            i = 1;
            S.left(m, 1) = globnidx(i,j,k);
            S.left(m, 3:5) = -nx;
            
            % Right face (x = dx)
            i = mx + 1;
            S.right(m, 1) = globnidx(i,j,k);
            S.right(m, 3:5) = nx;
            
            m = m + 1; % Increment node index
        end
    end
    
    m = 1;
    % Constant y
    for i = 2:mx % Nodes in x direction
        for k = 2:mz % Nodes in z direction
            
            % Front face (y = 0)
            j = 1;
            S.front(m, 1) = globnidx(i,j,k);
            S.front(m, 3:5) = -ny;
            
            % Back face (y = dy)
            j = my + 1;
            S.back(m, 1) = globnidx(i,j,k);
            S.back(m, 3:5) = ny;
            
            m = m + 1;
        end
    end
    
    m = 1;
    % Constant z
    for i = 2:mx
        for j = 2:my
            
            % Bottom face (z = 0)
            k = 1;
            S.bot(m, 1) =  globnidx(i,j,k);
            S.bot(m, 3:5) = -nz;
            
            % Top face (z = dz)
            k = mz + 1;
            S.top(m, 1) =  globnidx(i,j,k);
            S.top(m, 3:5) = nz;
            
            m = m + 1;
        end
    end
        
    %% Edges
    
    % Constant x
    exA = (0.5*w + 0.5*d)*h;
    S.topleft(:, 2) = exA;% S.topleft(:, 3:5) = (-nx + nz) / norm(nx + nz);
    S.topright(:, 2) = exA;% S.topright(:, 3:5) = (nx + nz) / norm(nx + nz);
    S.botleft(:, 2) = exA;% S.botleft(:, 3:5) = (-nx - nz) / norm(nx + nz);
    S.botright(:, 2) = exA;% S.botright(:, 3:5) = (nx - nz) / norm(nx + nz);
    
    % Constant y
    eyA = (0.5*h + 0.5*d)*w;
    S.topfront(:, 2) = eyA;% S.topfront(:, 3:5) = (-ny + nz) / norm(ny + nz);
    S.topback(:, 2) = eyA;% S.topback(:, 3:5) = (ny + nz) / norm(ny + nz);
    S.botfront(:, 2) = eyA;% S.botfront(:, 3:5) = (-ny - nz) / norm(ny + nz);
    S.botback(:, 2) = eyA;% S.botback(:, 3:5) = (ny - nz) / norm(ny + nz);
    
    % Constant z
    ezA = (0.5*w + 0.5*h)*d;
    S.frontleft(:, 2) = ezA;% S.frontleft(:, 3:5) = (-nx - ny) / norm(nx + ny);
    S.frontright(:, 2) = ezA;% S.frontright(:, 3:5) = (nx - ny) / norm(nx + ny);
    S.backleft(:, 2) = ezA;% S.backleft(:, 3:5) = (-nx + ny) / norm(nx + ny);
    S.backright(:, 2) = ezA;% S.backright(:, 3:5) = (nx + ny) / norm(nx + ny);
    
    %%% Indexing
    
    m = 1;
    % Constant x
    for j = 2:my
        
        % Top left edge
        i = 1;
        k = mz + 1;
        S.topleft(m, 1) = globnidx(i,j,k);
        S.topleft(m, 3:5) = (-nx + nz) / norm(nx + nz);
        
        % Top right edge
        i = mx + 1;
        k = mz + 1;
        S.topright(m, 1) = globnidx(i,j,k);
        S.topright(m, 3:5) = (nx + nz) / norm(nx + nz);
        
        % Bottom left edge
        i = 1;
        k = 1;
        S.botleft(m, 1) = globnidx(i,j,k);
        S.botleft(m, 3:5) = (-nx - nz) / norm(nx + nz);
        
        % Bottom right edge
        i = mx + 1;
        k = 1;
        S.botright(m, 1) = globnidx(i,j,k);
        S.botright(m, 3:5) = (nx - nz) / norm(nx + nz);
        
        m = m + 1;
    end
    
    m = 1;
    % Constant y
    for i = 2:mx
        
        % Top front edge
        j = 1;
        k = mz + 1;
        S.topfront(m, 1) = globnidx(i,j,k);
        S.topfront(m, 3:5) = (-ny + nz) / norm(ny + nz);
        
        % Top back edge
        j = my + 1;
        k = mz + 1;
        S.topback(m, 1) = globnidx(i,j,k);
        S.topback(m, 3:5) = (ny + nz) / norm(ny + nz);
        
        % Bottom front edge
        j = 1;
        k = 1;
        S.botfront(m, 1) = globnidx(i,j,k);
        S.botfront(m, 3:5) = (-ny - nz) / norm(ny + nz);
        
        % Bottom back edge
        j = my + 1;
        k = 1;
        S.botback(m, 1) = globnidx(i,j,k);
        S.botback(m, 3:5) = (ny - nz) / norm(ny + nz);
        
        m = m + 1;
    end
    
    m = 1;
    % Constant z
    for k = 2:mz
        
        % Front left edge
        i = 1;
        j = 1;
        S.frontleft(m, 1) = globnidx(i,j,k);
        S.frontleft(m, 3:5) = (-nx - ny) / norm(nx + ny);
        
        % Front right edge
        i = mx + 1;
        j = 1;
        S.frontright(m, 1) = globnidx(i,j,k);
        S.frontright(m, 3:5) = (nx - ny) / norm(nx + ny);
        
        % Back left edge
        i = 1;
        j = my + 1;
        S.backleft(m, 1) = globnidx(i,j,k);
        S.backleft(m, 3:5) = (-nx + ny) / norm(nx + ny);
        
        % Back right edge
        i = mx + 1;
        j = my + 1;
        S.backright(m, 1) = globnidx(i,j,k);
        S.backright(m, 3:5) = (nx + ny) / norm(nx + ny);
        
        m = m + 1;
    end
    
    %% Corners

    S.corners(1,1) = globnidx(1,1,1);
    S.corners(2,1) = globnidx(mx+1,1,1);
    S.corners(3,1) = globnidx(1,my+1,1);
    S.corners(4,1) = globnidx(mx+1,my+1,1);
    S.corners(5,1) = globnidx(1,1,mz+1);
    S.corners(6,1) = globnidx(mx+1,1,mz+1);
    S.corners(7,1) = globnidx(1,my+1,mz+1);
    S.corners(8,1) = globnidx(mx+1,my+1,mz+1);
    
    cA = (h*d + w*d + w*h) / 4; % Nodal area of corners
    S.corners(:,2) = cA;
    
    S.corners(1,3:5) = (-nx - ny - nz);
    S.corners(2,3:5) = (nx - ny - nz);
    S.corners(3,3:5) = (-nx + ny - nz);
    S.corners(4,3:5) = (nx + ny - nz);
    S.corners(5,3:5) = (-nx - ny + nz);
    S.corners(6,3:5) = (nx - ny + nz);
    S.corners(7,3:5) = (-nx + ny + nz);
    S.corners(8,3:5) = (nx + ny + nz);
    
    S.corners(:,3:5) = S.corners(:,3:5) / norm(nx + ny + nz);
    
    %% Concatenate surface sets
    
    % Convert surface sets to a cell structure:
    Sc = struct2cell(S);
    % Concatenate surface sets vertically:
    S.cat = cat(1,Sc{:,:}); % size (M,5), M = # surface nodes
    
end