%=========================================================================%
% Developed by Oxford Materials (as of 11/11/2020)

% Main EasyDD script.
%=========================================================================%

%% Initialisation

% Add relative paths.
[mainpath] = addPaths('../', ... % Main path
    'src/modules/','src/init/','output/'); % Selected paths

% Check for missing variables.
run inputCompletion.m

% Compile mex files.
[flags] = compileCode(flags);

% Cleanup input structures.
[rn, links] = cleanupnodes(rn, links);

% Generate connectivity of inputs.
[connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

% Check input consistency.
consistencycheck(rn, links, connectivity, linksinconnect);

% Construct stiffness matrix kg.
[FEM] = finiteElement3DCuboid(FEM, MU, NU);

% Define domain surface sets S.
[S] = surfaceCuboid(FEM);

% Prescribe boundary conditions on degrees of freedom.
[gamma, dofs] = feval(mods.prescribeDofs, S);

% Reformat stiffness matrix into K and pre-compute L,U decompositions.
[decompK] = reformatStiffnessMatrix(FEM, dofs);

fprintf('Finished FEM.\n')

% Plot initial surfaces.
plotFEMDomain(S, FEM)

% Initialise solution data as zeros or empty containers.
[u, u_hat, u_tilda, u_tilda_0, f, f_hat, fseg] = ...
    initialiseSolutions(FEM);

% Construct data structures needed for analytic tractions.
[f_tilda, tract] = ...
    AuxFEMCoupler(FEM, gamma, tract, flags);

% Use Delaunay triangulation to create surface mesh, used for visualisation
% dislocation remeshing algorithm.
[surfmesh] = ...
    MeshSurfaceTriangulation(S, FEM);

% Remesh considering surfaces in case input file is incorrect.
[rn, links, connectivity, linksinconnect, ~] = remesh_surf(...
    rn, links, connectivity, linksinconnect, fseg, ...
    FEM, surfmesh);

% Compute initial dislocation displacement field.
u_tilda_0 = calculateUtilda(rn, links, gamma, FEM, NU, u_tilda_0);

fprintf('Initialisation complete.\n');

%% Simulation
while simTime < totalSimTime

    % Perform DDD+FEM coupling.
    [u, u_hat, u_tilda, f, f_hat, f_tilda, r_hat] = FEM_DDD_Superposition(...
        rn, links, ...
        simTime, ...
        u, u_hat, u_tilda, u_tilda_0, ...
        f, f_hat, f_tilda, ...
        matpara, mods, flags, decompK, FEM, gamma, dofs, tract, diffBC);
    
    % Integrate equation of motion.
    [rnnew, vn, dt, fn, fseg] = feval(mods.integrator, ...
        rn, links, connectivity, ...
        u_hat, dt, dt0, ...
        matpara, mods, flags, FEM, Bcoeff);

    % Calculate plastic strain and plastic spin.
    [ep_inc, wp_inc] = calcPlasticStrainIncrement(rnnew, rn, links, FEM);
    
    % Print and store relevant BCs for plotting.
    [saveBC] = feval(mods.storeBC, ...
        u, u_hat, u_tilda, f, f_hat, f_tilda, r_hat, ...
        curstep, simTime, ...
        saveBC, gamma, dofs);
    
    % Plot simulation.
    plotSimulation(rn, links, plotFreq, viewangle, curstep, ...
        FEM, saveBC);

    [planeindex] = outofplanecheck(rn, links);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = updateMatricesForward(...
        rnnew, vn, links, connectivity, linksinconnect, fseg);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remeshPreCollision(...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
        u_hat, ...
        matpara, mods, flags, surfmesh, FEM, Bcoeff);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideNodesAndSegments(...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
        u_hat, curstep, ...
        matpara, mods, flags, FEM, Bcoeff);

    rnnew = fixBlockadingNodes(rnnew, connectivitynew);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = separation(...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
        u_hat, ...
        matpara, mods, flags, FEM, Bcoeff);

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_final(...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
        u_hat, ...
        matpara, mods, flags, FEM, Bcoeff);

    [rn, vn, links, connectivity, linksinconnect, fseg] = updateMatricesBackward(...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew);

    [curstep, simTime] = updateTime(curstep, simTime, dt);

    saveSimulation(simName, curstep, saveFreq)
    
    % Check that there exists nodes in the domain.
    if ~any(rn(:, 4) == 0)
        fprintf('No more real segments. Ending simulation.\n')
        saveSimulation(simName, curstep, 1)
        return;
    end

end

fprintf('Simulation complete.\n')
saveSimulation(simName, curstep, 1)
