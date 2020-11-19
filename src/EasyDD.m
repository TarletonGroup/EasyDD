%=========================================================================%
% EasyDD.m v2.0
%
% 3D Discrete Dislocation Plasticity. Is currently capable of simulating
% nanoindentation, micropillar compression and cantilever bending.
% Mixed-language model C, CUDA and Matlab so requires at least a C/C++
% compatible compiler, CUDA computation is optional. Explicitly calculates
% dislocation-dislocation interactions O(N^2).
%
%-------------------------------------------------------------------------%
%
% Data structure:
% rn: (:,4) array of nodal positions (x, y, z, label)
% vn: (:,3) array of nodal velocities (vx, vy, vz)
% fn: (:,3) array of nodal forces (fx, fy, fz)
% links: (:,8) array of links (idx1, idx2, bx, by, bz, nx, ny, nz)
%
%-------------------------------------------------------------------------%
%
% Compatibility:
% As long as Matlab and a compatible C/C++ and/or CUDA compiler are
% supported, the code will be compatible.
%
%-------------------------------------------------------------------------%
%
% Known issues and Improvement wish list:
% * Memory model is highly inefficient. There is a lot of naive dynamic
%   resizing of arrays.
% * Some matrices variables change in their number of columns are are hard
%   to track/debug (rn, connectivity).
% * Some variables hang around unused.
% * Collision, Separation and Remeshing can be antagonistic at higher
%   dislocation densities and may end up repeating processes until the
%   simulation progresses enough.
% * FEM boundary conditions need reworking so they can be passed in as
%   inputs rather than soft-coded into the FEM mesh builder.
% * There is no sanity check on input parameters.
% * Some performance-critical functions have no Matlab equivalent so a
%   C/C++ compiler is strictly required.
% * There are potential opportunities for further parallelisation in
%   collision, separation and integration processes.
% * Better modularisation would help increase the usability of the code.
% * Slip system in links(:, 6:8) may be inaccurate for some dislocations.
%
%=========================================================================%

%% Initialisation
% Fengxian's input reader.

% Check for missing variables.
% TODO: #10 make vertices and faces an argument, if they are not defined by the
% input provide a default.
% TODO: #11 in remesh_surf, make it so dislocations do not leave the domain via
% the fixed end depending on the simulation type.

% [fList, pList] = matlab.codetools.requiredFilesAndProducts('EasyDD.m');
% fList(:)
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
[vertices, B, xnodes, mno, nc, n, D, kg, w, h, d, mx, my, mz, mel] = ...
    finiteElement3D(dx, dy, dz, mx, my, mz, MU, NU);

[K, L, U, Sleft, Sright, Stop, Sbot, Sfront, Sback, Smixed, gammat, ...
        gammau, gammaMixed, fixedDofs, freeDofs, processForceDisp, plotForceDisp] = ...
    feval(simType, kg, w, h, d, mx, my, mz);

plotFEMDomain(Stop, Sbot, Sright, Sleft, Sfront, Sback, Smixed, xnodes)

% Construct data structures needed for analytic tractions.
[f, f_hat, para_tol, x3x6, n_se, gamma_dln, f_tilda_node, f_tilda_se, f_tilda, ...
        idxi, n_nodes_t, n_threads, para_scheme, gamma_disp, u_tilda_0, ...
        u, u_hat, u_tilda] = AuxFEMCoupler(mno, dx, dy, dz, mx, my, mz, ...
    xnodes, nc, gammat, gammau, gammaMixed, a_trac, CUDA_flag, ...
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

fprintf('Initialisation complete.\n');

while simTime < totalSimTime

    % DDD+FEM coupling
    [f_bar, f_hat, f_tilda, u_bar, u_hat, u_tilda, r_hat] = FEM_DDD_Superposition(...
        rn, links, a, MU, NU, xnodes, kg, L, U, gamma_disp, gammaMixed, fixedDofs, ...
        freeDofs, dx, dy, dz, simTime, mx, my, mz, sign_u_dot, u_dot, sign_f_dot, ...
        f_dot, u_tilda_0, u, u_hat, u_tilda, loading, a_trac, gamma_dln, x3x6, 4, ...
        n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, f_hat, CUDA_flag, ...
        n_threads, para_scheme, para_tol);

    [Fsim, Usim, t] = feval(processForceDisp, Fsim, f_bar, f_hat, f_tilda, Usim, u_bar, u_hat, u_tilda, ...
        r_hat, gammaMixed, fixedDofs, freeDofs, curstep, simTime);

    %integrating equation of motion
    [rnnew, vn, dt, fn, fseg] = feval(integrator, rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, Bcoeff, CUDA_flag);

    % plastic strain and plastic spin calculations
    [ep_inc, wp_inc] = calcPlasticStrainIncrement(rnnew, rn, links, (2 * plim)^3);

    plotSimulation(Usim, Fsim, rn, links, plim, vertices, plotFreq, viewangle, plotForceDisp, amag, mumag, curstep);

    [planeindex] = outofplanecheck(rn, links);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = updateMatricesForward(rnnew, vn, links, ...
        connectivity, linksinconnect, fseg);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remeshPreCollision(rnnew, linksnew, ...
        connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, ...
        mobility, doremesh, dovirtmesh, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, ...
        TriangleCentroids, TriangleNormals, CUDA_flag, Bcoeff);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideNodesAndSegments(docollision, ...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, lmin, CUDA_flag, Bcoeff, curstep);

    rnnew = fixBlockadingNodes(rnnew);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
        separation(doseparation, rnnew, linksnew, connectivitynew, linksinconnectnew, ...
        fsegnew, mobility, MU, NU, a, Ec, 2 * rann, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, CUDA_flag, Bcoeff);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, ...
        connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, ...
        mobility, doremesh, 0, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, ...
        TriangleNormals, CUDA_flag, Bcoeff);

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

fprintf('Simulation completed.\n')
saveSimulation(simName, curstep, 1)
