run Thesis_Zinc_100_tensile.m
run inputCompletion.m

% Compile mex files.
CUDA_flag = compileCode(CUDA_flag);

% Cleanup input structures.
[rn, links] = cleanupnodes(rn, links);

% Generate connectivity of inputs.
[connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

% Check input consistency.
consistencycheck(rn, links, connectivity, linksinconnect);

% Construct stiffeness matrix K and pre-compute L,U decompositions.
[S, vertices, B, xnodes, mno, nc, n, D, kg, w, h, d, mx, my, mz, mel] = ...
    finiteElement3D(dx, dy, dz, mx, my, mz, MU, NU);

[K, L, U, P_l, P_u, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, ...
        gammau, gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = ...
    simType(kg, w, h, d, mno, mx, my, mz, S);

plotFEMDomain(Stop, Sbot, Sright, Sleft, Sfront, Sback, Smixed, xnodes)

% Construct data structures needed for analytic tractions.
[f, f_hat, para_tol, x3x6, n_se, gamma_dln, f_tilda_node, f_tilda_se, f_tilda, ...
        idxi, n_nodes_t, n_threads, para_scheme, gamma_disp, u_tilda_0, ...
        u, u_hat, u_tilda] = AuxFEMCoupler(mno, dx, dy, dz, mx, my, mz, ...
    xnodes, nc, gammat, gammau, gammaMixed, calculateTractions, CUDA_flag, ...
    n_threads, para_scheme);

% Use Delaunay triangulation to create surface mesh, used for visualisation
% dislocation remeshing algorithm.
[TriangleCentroids, TriangleNormals, tri, Xb] = ...
    MeshSurfaceTriangulation(xnodes, Stop, Sbot, Sfront, Sback, Sleft, ...
    Sright, gammaMixed);

%Remesh considering surfaces in case input file incorrect.
[rn, links, connectivity, linksinconnect] = remesh_surf(rn, links, ...
    connectivity, linksinconnect, vertices, TriangleCentroids, ...
    TriangleNormals);

u_tilda_0 = calculateUtilda(rn, links, gamma_disp, NU, xnodes, dx, ...
    dy, dz, mx, my, mz, u_tilda_0);
close all

fseg = segforcevec(MU, NU, a, Ec, rn, links, 0, vertices, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, CUDA_flag);

    %mobility function
[vn, fn] = mobfcc0(fseg, rn, links, connectivity, [], [], Bcoeff, rotMatrix);