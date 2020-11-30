function [surfmesh] = MeshSurfaceTriangulation(S, FEM)
    %===============================================================%
    % Oxford Materials (11/11/2020)
    
    % Beginning of surface remeshing for surface nodes. Creates mesh 
    % of adjacent triangles on the geometry surface. The surface FEM
    % nodes are the vertices of the mesh triangles.
    %===============================================================%
    
    %% Extraction
    
    % Surface sets:
    Scat = S.cat;
    
    % FEM:
    xnodes = FEM.xnodes;
    
    %% Meshing
    
    % Create triangulation (uses all surface nodes from FEM code):
    pts = xnodes(Scat(:, 1), :);
    dt = delaunayTriangulation(pts); % Fine triangulation

    [tri, Xb] = freeBoundary(dt);
    TR = triangulation(tri, Xb);

    % Calculate centroids of triangulated surfaces:
    P = incenter(TR);

    % Calculate normals of triangulated surfaces:
    fn = faceNormal(TR);
    
    %% Storage
    
    surfmesh = struct; % Stores surface triangulation mesh data
    
    surfmesh.TriangleCentroids = P;
    surfmesh.TriangleNormals = fn;
    surfmesh.tri = tri;
    surfmesh.Xb = Xb;
    
end
