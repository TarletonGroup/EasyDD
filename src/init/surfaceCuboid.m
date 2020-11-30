function [S] = surfaceCuboid(FEM)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Defines the geometry of the simulation domain surfaces.
    % This function specifically defines a cantilever (cuboid) with
    % orthorhombic axes.
    %===============================================================%
    
    %% Extract data from structures
    
    % FEM:
    w = FEM.w; h = FEM.h; d = FEM.d; % Finite element dimensions
    mx = FEM.mx; my = FEM.my; mz = FEM.mz; % # FEs per dimension
    globnidx = FEM.globnidx; % Matrix of global FE node indices
    nx = FEM.normals(1,:); ny = FEM.normals(2,:); nz = FEM.normals(3,:);
    
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
