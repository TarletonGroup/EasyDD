%=========================================================================%
% Developed by Oxford Materials (as of 11/11/2020)

% Main EasyDD script.
%=========================================================================%

%% Initialisation

% Add relative paths.
run addPaths.m

% Check for missing variables.
run inputCompletion.m

% Compile mex files.
CUDA_flag = compileCode(CUDA_flag);

% Cleanup input structures.
[rn, links] = cleanupnodes(rn, links);

% Generate connectivity of inputs.
[connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

% Check input consistency.
consistencycheck(rn, links, connectivity, linksinconnect);

% Construct stiffness matrix kg.
[B, xnodes, mno, nc, n, D, kg, w, h, d, my, mz, mel, globnidx] = ...
    finiteElement3D(dx, dy, dz, mx, MU, NU);

% Define domain surfaces.
[S, Scat] = cantilever(w, h, d, mx, my, mz, globnidx);

% Prescribe boundary conditions on degrees of freedom.
[gamma, fixedDofs, freeDofs] = cantileverBC(S, maxstepsBC);

% Reformat stiffness matrix into K and pre-compute L,U decompositions.
[K, L, U, bcwt] = cantileverK(kg, fixedDofs, freeDofs);

fprintf('Finished FEM.\n')

plotFEMDomain(S, xnodes)

% Construct data structures needed for analytic tractions.
[f, f_hat, para_tol, x3x6, n_se, gamma_dln, f_tilda_node, f_tilda_se, f_tilda, ...
        idxi, n_nodes_t, n_threads, para_scheme, gamma_disp, u_tilda_0, ...
        u, u_hat, u_tilda] = AuxFEMCoupler(mno, dx, dy, dz, mx, my, mz, ...
    xnodes, nc, gamma, a_trac, CUDA_flag, ...
    n_threads, para_scheme);

% Use Delaunay triangulation to create surface mesh, used for visualisation
% dislocation remeshing algorithm.
[TriangleCentroids, TriangleNormals, tri, Xb] = ...
    MeshSurfaceTriangulation(xnodes, Scat);

% Remesh considering surfaces in case input file is incorrect.
[rn, links, connectivity, linksinconnect] = remesh_surf(rn, links, ...
    connectivity, linksinconnect, vertices, TriangleCentroids, ...
    TriangleNormals);

u_tilda_0 = calculateUtilda(rn, links, gamma_disp, NU, xnodes, dx, ...
    dy, dz, mx, my, mz, u_tilda_0);

fprintf('Initialisation complete.\n');

%% Simulation
while simTime < totalSimTime

    % DDD+FEM coupling
    [f, f_hat, f_tilda, u, u_hat, u_tilda, r_hat, BCstore] = FEM_DDD_Superposition(...
        rn, links, a, MU, NU, xnodes, kg, L, U, bcwt, gamma_disp, gamma, S, ...
        fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, diffBC, ...
        u_tilda_0, u, u_hat, u_tilda, loading, a_trac, gamma_dln, x3x6, 4, ...
        n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, f_hat, CUDA_flag, ...
        n_threads, para_scheme, para_tol);

    %integrating equation of motion
    [rnnew, vn, dt, fn, fseg] = feval(integrator, rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, mobstruct, CUDA_flag);

    % plastic strain and plastic spin calculations
    [ep_inc, wp_inc] = calcPlasticStrainIncrement(rnnew, rn, links, (2 * plim)^3);

    plotSimulation(rn, links, plim, vertices, plotFreq, viewangle, plotForceDisp, BCstore, simscale, curstep);

    [planeindex] = outofplanecheck(rn, links);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = updateMatricesForward(rnnew, vn, links, ...
        connectivity, linksinconnect, fseg);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remeshPreCollision(rnnew, linksnew, ...
        connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, ...
        mobility, doremesh, dovirtmesh, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, ...
        TriangleCentroids, TriangleNormals, CUDA_flag, mobstruct);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideNodesAndSegments(docollision, ...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, lmin, CUDA_flag, mobstruct, curstep);

    rnnew = fixBlockadingNodes(rnnew, connectivitynew);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
        separation(doseparation, rnnew, linksnew, connectivitynew, linksinconnectnew, ...
        fsegnew, mobility, MU, NU, a, Ec, 2 * rann, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, CUDA_flag, mobstruct);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, ...
        connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, ...
        mobility, doremesh, 0, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, ...
        TriangleNormals, CUDA_flag, mobstruct);

    [rn, vn, links, connectivity, linksinconnect, fseg] = updateMatricesBackward(rnnew, ...
        linksnew, connectivitynew, linksinconnectnew, fsegnew);

    [curstep, simTime] = updateTime(curstep, simTime, dt);

    saveSimulation(simName, curstep, saveFreq)

    if ~any(rn(:, 4) == 0)
        fprintf('No more real segments, ending simulation.\n')
        saveSimulation(simName, curstep, 1)
        return;
    end

end

fprintf('Simulation complete.\n')
saveSimulation(simName, curstep, 1)
