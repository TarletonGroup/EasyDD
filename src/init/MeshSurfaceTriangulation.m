function [P, fn, tri, Xb] = MeshSurfaceTriangulation(xnodes, Scat)
        
    % Beginning of surface remeshing for surface node. %%%%%%%%%%%%%%%%%%%%

    %create triangulation (uses all surface nodes from FEM code)
    pts = xnodes(Scat(:, 1), :);
    dt = delaunayTriangulation(pts);
    %dt = delaunayTriangulation(vertices); %coarse triangulation

    [tri, Xb] = freeBoundary(dt);
    TR = triangulation(tri, Xb);

    %calculate centroids of triangulated surfaces
    P = incenter(TR);

    %calculate normals of triangulated surfaces
    fn = faceNormal(TR);

end
