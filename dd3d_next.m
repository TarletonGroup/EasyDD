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
% TODO: make vertices and faces an argument, if they are not defined by the
% input provide a default.
% TODO: in remesh_surf, make it so dislocations do not leave the domain via
% the fixed end depending on the simulation type.

% Compile mex files.
CUDA_flag = compileCode(CUDA_flag);

% Cleanup input structures.
[rn, links] = cleanupnodes(rn, links);

% Generate connectivity of inputs.
[connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

% Check input consistency.
consistencycheck(rn, links, connectivity, linksinconnect);

% Construct stiffeness matrix K and pre-compute L,U decompositions.
% gammat -> nodes where tractions = 0
% gammau -> nodes where displacements = 0
% gammaLoad -> nodes where tractions != 0
% gammaDisp -> nodes where displacements != 0
[vertices, B, xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, ...
        Sbot, Sfront, Sback, gammat, gammau, gammaMixed, fixedDofs, ...
        freeDofs, w, h, d, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, ...
    MU, NU, loading);

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

disp('Initialisation complete.');
%%
if (~exist('dt'))
    dt = dt0;
end

dt = min(dt, dt0);
plotCounter = 1;
close all

Fend = zeros(1e6, 1); fend = [];
U_bar = zeros(1e6, 1); Ubar = [];
t = zeros(1e6, 1); simTime = 0;
sign_u_dot = -1;
sign_f_dot = -1;
u_dot = 100 * 1E3 * dx * (1E-4/160E9) * 100; % Marielle
f_dot = 0;
simType = 1;

while simTime < totalSimTime

    % DDD+FEM coupling

    [f_hat, u_hat, r_hat] = FEM_DDD_Superposition(rn, links, a, MU, NU, ...
        xnodes, kg, L, U, gamma_disp, gammaMixed, fixedDofs, ...
        freeDofs, dx, dy, dz, simTime, mx, my, mz, sign_u_dot, u_dot, sign_f_dot, f_dot, u_tilda_0, u, ...
        u_hat, u_tilda, simType, a_trac, gamma_dln, x3x6, 4, ...
        n_nodes_t, n_se, idxi, f, f_tilda_node, f_tilda_se, f_tilda, f_hat, CUDA_flag, ...
        n_threads, para_scheme, para_tol);

    if a_trac == 0
        [u_hat, fend, Ubar] = FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
            gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, u_tilda_0);
    else
        [u_hat, fend, Ubar] = analytic_FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
            gamma_disp, gammat, gammaMixed, fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, u_tilda_0, ...
            gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, f, ...
            f_tilda_node, f_tilda_se, f_tilda, f_hat, use_gpu, n_threads, para_scheme, para_tol);
    end

    Fend(curstep + 1) = fend;
    U_bar(curstep + 1) = Ubar;
    t(curstep + 1) = simTime;

    fprintf('fend = %d, Ubar = %d, simTime = %d \n', fend, Ubar, simTime);

    %integrating equation of motion
    [rnnew, vn, dt, fn, fseg] = feval(integrator, rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d);
    % plastic strain and plastic spin calculations
    [ep_inc, wp_inc] = calcplasticstrainincrement(rnnew, rn, links, (2 * plim)^3);

    if (mod(curstep, printfreq) == 0) && not(isempty(vn))
        fprintf('step%3d dt=%e v%d=(%e,%e,%e) \n', ...
            curstep, dt, printnode, vn(printnode, 1), vn(printnode, 2), vn(printnode, 3));
        close all;
        save restart_temp
    elseif (mod(curstep, printfreq) == 0)
        fprintf('step%3d dt=%e v%d=(%e,%e,%e) \n', ...
            curstep, dt, printnode);
    end

    if (mod(curstep, plotfreq) == 0)
        plotnodes(rn, links, plim, vertices);
        view(viewangle);
        %schematicplot2(rn,links,vertices,U_bar,Fend,amag,dx,totalSimTime);
        drawnow
        pause(0.01);
        %
        %         figure(2);
        %          clf
        %          plot(U_bar*bmag,-Fend*bmag^2*mumag)
        %          pause(0.01);
    end

    [planeindex] = outofplanecheck(rn, links);

    rnnew = [rnnew(:, 1:3) vn rnnew(:, 4)];
    linksnew = links;
    connectivitynew = connectivity;
    linksinconnectnew = linksinconnect;
    fsegnew = fseg;

    if (doremesh)%do virtual re-meshing first
        %remeshing virtual dislocation structures
        if (dovirtmesh)
            %[rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen_mex(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,DIST_SOURCE*0.49,dx,MU,NU,a,Ec);
            %[rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,DIST_SOURCE*0.49,dx,MU,NU,a,Ec);
            [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = virtualmeshcoarsen3(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, MU, NU, a, Ec, dx, dy, dz);
            %[rnnew,linksnew,connectivitynew,linksinconnectnew] = virtualmeshcoarsen2(rnnew,linksnew,maxconnections,10*lmin);
        end

        %remeshing internal dislocation structures
        [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, mobility, doremesh, dovirtmesh, vertices, ...
            u_hat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, TriangleNormals);
    end

    %save restart.mat
    if (docollision)
        s1_old = 0;
        s2_old = 0;
        %collision detection and handling

        %COLLISION GPU / GPUPLUS  Marielle  + CollisionCheckerMexMarielle
        %it requires a GPU device
        %it requires CollisionCheckerMexMarielle, collisionGPUplus (or collision GPU), mindistcalcGPU1, mindistcalcGPU2,CreateInputMex
        colliding_segments = 1;

        while colliding_segments == 1
            [colliding_segments, n1s1, n2s1, n1s2, n2s2, floop, s1, s2, segpair] = CollisionCheckerMex(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
                rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann);

            %                   if colliding_segments == 1 %scan and update dislocation structure.
            %                         reset(gpuDevice)
            %                         [rnnew,linksnew,~,~,fsegnew]=...
            %                         collisionGPUplus(rnnew,linksnew,connectivitynew,linksinconnectnew,...
            %                         fsegnew,rann,MU,NU,a,Ec,mobility,vertices,u_hat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair);
            %                   end

            %              INTIAL COLLISION
            %              [colliding_segments]=CollisionCheckerMex(rnnew(:,1),rnnew(:,2),rnnew(:,3),rnnew(:,end),...
            %              rnnew(:,4),rnnew(:,5),rnnew(:,6),linksnew(:,1),linksnew(:,2),connectivitynew,rann);
            if colliding_segments == 1%scan and update dislocation structure.

                if floop == 1
                    fprintf("Step %d. Unconnected links found. Links %d and %d are colliding.\n", curstep, s1, s2)
                elseif floop == 2
                    fprintf("Step %d. Links %d and %d colliding by hinge condition.\n", curstep, s1, s2)
                end

                %                 [rnnew,linksnew,~,~,fsegnew]=...
                %                 collision(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                %                 fsegnew,rann,MU,NU,a,Ec,mobility,vertices,u_hat,nc,xnodes,D,mx,mz,w,h,d);

                %                 fprintf(' Segment %d and segment %d are colliding\n',s1,s2);
                if colliding_segments == 1
                    [rnnew, linksnew, ~, ~, fsegnew, colliding_segments] = collision_basic(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, ...
                        u_hat, nc, xnodes, D, mx, mz, w, h, d, floop, n1s1, n2s1, n1s2, n2s2, s1, s2, segpair, lmin);

                    %removing links with effective zero Burgers vectors
                    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = cleanupsegments(rnnew, linksnew, fsegnew);

                    %                 [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,0,vertices,...
                    %     u_hat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);

                    %                 if (doseparation)
                    %                     if max(connectivitynew(:,1))>3
                    %                         %spliting of nodes with 4 or more connections
                    %                         [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=...
                    %                             separation(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                    %                             fsegnew,mobility,MU,NU,a,Ec,2*rann,vertices,u_hat,nc,xnodes,D,mx,mz,w,h,d);
                    %                     end
                    %                 end
                end

                %                 [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,0,vertices,...
                %         u_hat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
            end

        end

    end

    %find blockading fixed nodes and allow code to deal with them
    for i = 1:size(rnnew, 1)

        if rnnew(i, end) ~= 7
            continue
        end

        if connectivitynew(i, 1) > 7
            rnnew(i, end) = 0;
        end

    end

    if (doseparation)

        if max(connectivitynew(:, 1)) > 3
            %spliting of nodes with 4 or more connections
            [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
                separation(rnnew, linksnew, connectivitynew, linksinconnectnew, ...
                fsegnew, mobility, MU, NU, a, Ec, 2 * rann, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d);
        end

    end

    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, mobility, doremesh, 0, vertices, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, TriangleNormals);

    rn = [rnnew(:, 1:3) rnnew(:, 7)];
    vn = rnnew(:, 4:6);
    links = linksnew;
    connectivity = connectivitynew;
    linksinconnect = linksinconnectnew;
    fseg = fsegnew;

    %store run time information
    %time step
    curstep = curstep + 1;
    %data(curstep,1)=dt;
    simTime = simTime + dt;

    %    save restart;
    %     if all(rn(:,4)==67) %no more real segments, stop simulation
    %         disp('Dislocation-free real domain. Only virtual segments remain!');
    %         disp('Computing distorted mesh. Please wait...');
    %         [u_tilda]=visualise(rn,links,NU,D,mx,my,mz,mel,mno,xnodes,nc,...
    %             dx,dy,dz,w,h,d,vertices,u_hat);
    %         return;
    %     end

end

save restart;
disp('completed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
