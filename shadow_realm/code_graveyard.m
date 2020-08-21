%%=======================================================================%%
% This file unused functions and versions of functions. We keep them in
% case bits and pieces prove useful in the future. However, anything in
% this file should be treated as a repository of deprecated and outdated
% functionality whose only use is for reference.
%%=======================================================================%%

function EasyDD()
    %%=======================================================================%%
    % One of the many old versions of EasyDD.
    %%=======================================================================%%
    % %Dislocation Dynamics simulation in 3-dimension
    % % Features:
    % % mixed boundary conditions cantilever.
    % % linear mobility law (mobfcc0,mobfcc1)
    % % N^2 interaction (no neighbor list, no fast multipole)
    %
    % %Data structure:
    % %NMAX:    maximum number of nodes (including disabled ones)
    % %LINKMAX: maximum number of links (including disabled ones)
    % %rn: (NMAX,4) array of nodal positions (last column is flag: -1 means disabled)
    % %vn: (NMAX,3) array of nodal velocities
    % %fn: (NMAX,3) array of nodal forces
    % %links: (LINKMAX,8) array of links (id1,id2,bx,by,bz,nx,ny,nz)
    %
    % compile the c source code for seg-seg force evaluation and makes a dynamic linked library
    %  disp('Compiling MEX files');
    %  mex SegSegForcesMex.c
    %  mex StressDueToSegsMex.c
    %  mex UtildaMex.c
    %  mex MinDistCalcMex.c
    %  mex  CollisionCheckerMex.c
    %  mex mobbcc1mex.c
    %  mex displacementmex_et.c
    %  mex CreateInputMex.c %CollisionMarielle
    %  mex CollisionCheckerMexMarielle.c %CollisionMarielle
    %  mex ./nodalForce/nodal_surface_force_linear_rectangle_mex.c
    %  mexcuda -v COMPFLAGS="-Xptxas -v -arch=compute_60 -code=sm_60 -use_fast_math" ./nodalForce/nodal_surface_force_linear_rectangle_mex_cuda.cu
    % disp('Done!');
    %  default value if run by itself (e.g. not through "rundd3d")
    %  cleanup the empty node and link entries at the end of the initial data structures
    [rn, links] = cleanupnodes(rn, links);

    % genererate the connectivity list from the list of links
    disp('Initiliazing connectivity list. Please wait.');
    [connectivity, linksinconnect] = genconnectivity(rn, links, maxconnections);

    consistencycheck(rn, links, connectivity, linksinconnect);
    disp('Consistencycheck : Done!');

    %construct stiffeness matrix K and pre-compute L,U decompositions.
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
    [vertices, B, xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, Sbot, ...
            Sfront, Sback, gammat, gammau, gammaMixed, fixedDofs, freeDofs, ...
            w, h, d, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, MU, NU, loading); % Haiyang's addition

    % [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    %     Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    %     w,h,d,my,mz,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D_haiyang2(dx,dy,dz,mx,MU,NU,loading);

    %% Addition by Daniel Celis Garza 12/19/2017.
    % Precomputing variables needed to couple tractions induced by dislocations
    % to FEM. Generate the list of planes to be extracted. The FEM coupler is
    % hardcoded so that the min(x) yz-plane does not have traction boundary
    % conditions.
    f_hat = zeros(3 * mno, 1);

    planes = (1:1:6)';
    [x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx; my; mz], ...
        planes, 4);
    gamma_dln = [gammat(:, 1); gammaMixed(:, 1)];
    [f_dln_node, f_dln_se, ...
            f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
    para_tol = dx / 1e7;
    % Parallel CUDA C flags.
    if use_gpu == 1
        % Provide a default number of threads in case none is given.
        if ~exist('n_threads', 'var')
            n_threads = 256;
        end %if

        % Provide a default parallelisaion scheme in case none is given.
        if ~exist('para_scheme', 'var')
            % Parallelise over dislocations.
            para_scheme = 1;
        end %if

    else
        n_threads = 0;
        para_scheme = 0;
    end %if

    gamma_disp = [gammau; gammaMixed];

    disp('Done! Initializing simulation.');

    global USE_GPU;
    USE_GPU = 0; %0 if CPU only.

    if (USE_GPU == 1)
        disp('Going to use GPU as well...'); %setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\amd64']);
        system('nvcc -ptx -m64 -arch sm_35 SegForceNBodyCUDADoublePrecision.cu');
    end

    %Use Delaunay triangulation to create surface mesh, used for visualisation

    %and for dislocation remeshing algorithm.
    [TriangleCentroids, TriangleNormals, tri, Xb] = ...
        MeshSurfaceTriangulation(xnodes, Stop, Sbot, Sfront, Sback, Sleft, Sright, gammaMixed);
    %Remesh considering surfaces in case input file incorrect.
    %disp('Remeshing...');
    [rn, links, connectivity, linksinconnect] = remesh_surf(rn, links, connectivity, linksinconnect, vertices, TriangleCentroids, TriangleNormals);

    %plot dislocation structure
    % figure(1);
    % plotHandle = plotnodes(rn,links,plim,vertices); view(viewangle);
    % drawnow

    % data=zeros(totalsteps,1);
    if (~exist('dt', 'var'))
        dt = dt0;
    end

    dt = min(dt, dt0);
    plotCounter = 1;
    close all

    Fend = zeros(1e6, 1); fend = [];
    U_bar = zeros(1e6, 1); Ubar = [];
    t = zeros(1e6, 1); simTime = 0;
    %%
    gamma = [gammau; gammaMixed];
    utilda_0 = zeros(3 * mno, 1);
    gn = gamma(:, 1);
    [Ux, Uy, Uz] = Utilda_bb3(rn, links, gn, NU, xnodes, dx, dy, dz, mx, my, mz);

    utilda_0(3 * gn -2) = Ux;
    utilda_0(3 * gn -1) = Uy;
    utilda_0(3 * gn) = Uz;
    %%
    while simTime < totalSimTime

        % frame recording
        intSimTime = intSimTime + dt;

        if intSimTime > dtplot && doplot == 1
            %plotHandle=plotnodes(rn,links,plim,vertices);view(viewangle);
            plotCounter = plotCounter + 1;
            plotCounterString = num2str(plotCounter, '%03d');
            %saveas(plotHandle,plotCounterString,'png')
            %save(plotCounterString,'rn','links','fend','Ubar','simTime');
            %plotHandle=plotnodes(rn,links,plim,vertices);view(viewangle);
            %plotnodes(rn,links,plim,vertices);view(viewangle);
            %zlim([-100 100])
            %xlim([-100 100])
            %ylim([-100 100])
            intSimTime = intSimTime - dtplot;
        end

        %DDD+FEM coupling
        if a_trac == 0
            [uhat, fend, Ubar] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
                gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, utilda_0);
        else
            [uhat, fend, Ubar] = analytic_FEMcoupler(rn, links, a, MU, NU, xnodes, mno, kg, L, U, ...
                gamma_disp, gammat, gammaMixed, fixedDofs, freeDofs, dx, dy, dz, simTime, mx, my, mz, utilda_0, ...
                gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, ...
                f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, para_tol);
        end

        Fend(curstep + 1) = fend;
        U_bar(curstep + 1) = Ubar;
        t(curstep + 1) = simTime;
        fprintf('fend = %d, Ubar = %d, simTime = %d \n', fend, Ubar, simTime);

        %     if (dovirtmesh)
        %         %remeshing virtual dislocation structures
        %         %[rn,links,connectivity,linksinconnect]=remesh_surf(rn,links,connectivity,linksinconnect,vertices,TriangleCentroids,TriangleNormals);
        %         %[rn,links,connectivity,linksinconnect] = virtualmeshcoarsen2(rn,links,maxconnections,10*lmin);
        %     end

        %integrating equation of motion
        [rnnew, vn, dt, fn, fseg] = feval(integrator, rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
            rmax, rntol, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);

        % Haiyang's screw dislocation correction.
        [rn, links] = cross_slip_BCCrotate5(fseg, rn, links, connectivity, 0, ...
            0, curstep, 0, NU, a, Ec, 0, ...
            uhat, nc, xnodes, D, mx, mz, w, h, d);

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
            drawnow
            pause(0.01);
            %
            %         figure(2);
            %          clf
            %          plot(U_bar*bmag,-Fend*bmag^2*mumag)
            %          pause(0.01);
        end

        rnnew = [rnnew(:, 1:3) vn rnnew(:, 4)];
        linksnew = links;
        connectivitynew = connectivity;
        linksinconnectnew = linksinconnect;
        fsegnew = fseg;

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
                [colliding_segments, n1s1, n2s1, n1s2, n2s2, floop, s1, s2, segpair] = CollisionCheckerMexMariellebis(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
                    rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann);

                if floop == 2
                    colliding_segments = 0;
                end

                %                   if colliding_segments == 1 %scan and update dislocation structure.
                %                         reset(gpuDevice)
                %                         [rnnew,linksnew,~,~,fsegnew]=...
                %                         collisionGPUplus(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                %                         fsegnew,rann,MU,NU,a,Ec,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair);
                %                   end

                %              INTIAL COLLISION
                %              [colliding_segments]=CollisionCheckerMex(rnnew(:,1),rnnew(:,2),rnnew(:,3),rnnew(:,end),...
                %              rnnew(:,4),rnnew(:,5),rnnew(:,6),linksnew(:,1),linksnew(:,2),connectivitynew,rann);
                if colliding_segments == 1%scan and update dislocation structure.

                    %                 [rnnew,linksnew,~,~,fsegnew]=...
                    %                 collision(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                    %                 fsegnew,rann,MU,NU,a,Ec,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);

                    if floop == 2
                        %                     if s1==s1_old && s2==s2_old
                        colliding_segments = 0;
                        %                     else
                        %                         s1_old=s1;
                        %                         s2_old=s2;
                        %                     end
                    end

                    if mod(curstep, 100) ~= 0 && floop == 1
                        nodeids1 = connectivitynew(n1s1, 1);
                        nodeids1 = 1:nodeids1;
                        nodeids1 = 2 * nodeids1;
                        nodeids2 = connectivitynew(n2s1, 1);
                        nodeids2 = 1:nodeids2;
                        nodeids2 = 2 * nodeids2;
                        nodeids3 = connectivitynew(n1s2, 1);
                        nodeids3 = 1:nodeids3;
                        nodeids3 = 2 * nodeids3;
                        nodeids4 = connectivitynew(n2s2, 1);
                        nodeids4 = 1:nodeids4;
                        nodeids4 = 2 * nodeids4;
                        alt1 = connectivitynew(n1s1, nodeids1);
                        alt2 = connectivitynew(n2s1, nodeids2);
                        alt3 = connectivitynew(n1s2, nodeids3);
                        alt4 = connectivitynew(n2s2, nodeids4);
                        other1 = any(alt1 == alt3');
                        other1 = sum(other1);
                        other2 = any(alt1 == alt4');
                        other2 = sum(other2);
                        other3 = any(alt2 == alt3');
                        other3 = sum(other3);
                        other4 = any(alt2 == alt4');
                        other4 = sum(other4);

                        if other1 > eps || other3 > eps || other2 > eps || other4 > eps
                            fprintf('Remeshing conflict detected. Cancelling collision\n')
                            colliding_segments = 0;
                        end

                    end

                    fprintf(' Segment %d and segment %d are colliding\n', s1, s2);

                    if colliding_segments == 1
                        [rnnew, linksnew, ~, ~, fsegnew, colliding_segments] = collision_basic(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, ...
                            uhat, nc, xnodes, D, mx, mz, w, h, d, floop, n1s1, n2s1, n1s2, n2s2, s1, s2, segpair);
                    end

                end

                %removing links with effective zero Burgers vectors
                [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = cleanupsegments(rnnew, linksnew, fsegnew);
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

        if (doseparation) && max(connectivitynew(:, 1)) > 3
            %spliting of nodes with 4 or more connections
            [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = ...
                separation(rnnew, linksnew, connectivitynew, linksinconnectnew, ...
                fsegnew, mobility, MU, NU, a, Ec, 2 * rann, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
        end

        if (doremesh)%do virtual re-meshing first
            %remeshing virtual dislocation structures
            if (dovirtmesh)
                %[rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen_mex(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,DIST_SOURCE*0.49,dx,MU,NU,a,Ec);
                [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = virtualmeshcoarsen(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, DIST_SOURCE * 0.49, dx, MU, NU, a, Ec);
                %[rnnew,linksnew,connectivitynew,linksinconnectnew] = virtualmeshcoarsen2(rnnew,linksnew,maxconnections,10*lmin);
            end

            %remeshing internal dislocation structures
            [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = remesh_all(rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, lmin, lmax, areamin, areamax, MU, NU, a, Ec, mobility, doremesh, dovirtmesh, vertices, ...
                uhat, nc, xnodes, D, mx, mz, w, h, d, TriangleCentroids, TriangleNormals);
        end

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
        %         [utilda]=visualise(rn,links,NU,D,mx,my,mz,mel,mno,xnodes,nc,...
        %             dx,dy,dz,w,h,d,vertices,uhat);
        %         return;
        %     end
        %save restart;
        %       if (curstep >= 2004)
        %           fprintf("\n")
        %           save(sprintf('20180925_gpu_debug_%d_%d', n_threads, curstep));
        %       end
        %     if mod(curstep, 1000)==0
        %         close all
        %         plotnodes(rn,links,plim,vertices);
        %         pause
        % % %         save(sprintf('./mat_files/debris_%d', curstep));
        %     end
        %     if -Fend(curstep)*amag^2*mumag < -2580
        %         save(sprintf('./mat_files/20200128_%d_%d_%d_HYmob', curstep, a_trac, simple));
        %         pause
        %     end
        if (mod(curstep, 1000) == 0)
            close all
            save(sprintf('./mat_files/HaiMesha_%d', curstep));
            %         save(sprintf('./mat_files/20200128_%d_%d_%d_HYmob', curstep, a_trac, simple));
            %         save(sprintf('./mat_files/20191216_%d_%d_%d_HYmob', curstep, a_trac, simple));
            %         save(sprintf('./mat_files/20180928_%d', curstep));
            %         save(sprintf('20181001_gpu_%d_%d', n_threads, curstep));
        end

    end

    disp('completed')
end

%%
function compile(CUDA_flag)
    file = dir(sprintf("UtildaMex.%s", ext));

    if ~isfile("UtildaMex.c") ||~isfile(file.name) && days(file.date - datetime('now')) > 30
        mex -v COPTIMFLAGS = "-o3 -oy -DNDEBUG" UtildaMex.c
    end

end

%%
function visualise
    %%

    % % code below is just to visulise the defomed mesh

    clear all
    close all
    disp('Loading restart file');
    load restart.mat%restart.50sources.75744.mat

    plotmesh = 1;
    plotstress = 0;

    %refine mesh
    mx = mx * 2;

    %update results
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
    [B, xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, Sbot, ...
            Sfront, Sback, gammat, gammau, gammaMixed, fixedDofs, freeDofs, ...
            w, h, d, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, MU, NU, loading);
    disp('Done!');

    [TriangleCentroids, TriangleNormals, tri, Xb] = ...
        MeshSurfaceTriangulation(xnodes, Stop, Sbot, Sfront, Sback, Sleft, Sright);

    [uhat, fend, Ubar] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, simTime);

    segments = constructsegmentlist(rn, links);
    utilda = zeros(3 * mno, 1);
    gn = 1:mno; % global node number
    x0 = xnodes(gn, 1:3); % field point
    point_array_length = size(x0, 1);
    segments_array_length = size(segments, 1);
    %Full C version (including combinatorials, checks, corrections etc.)

    disp('Calculating displacements');
    [Ux, Uy, Uz] = UtildaMex(x0(:, 1), x0(:, 2), x0(:, 3), ...%coordinates
    segments(:, 3), segments(:, 4), segments(:, 5), ...%burgers vector
    segments(:, 6), segments(:, 7), segments(:, 8), ...%start node segs
    segments(:, 9), segments(:, 10), segments(:, 11), ...%end node segs
    segments(:, 12), segments(:, 13), segments(:, 14), ...%slip plane
    NU, point_array_length, segments_array_length);
    %[Uxf, Uyf, Uzf] =displacement_fivel(x0,segments,NU); %gives same answer

    utilda(3 * gn -2) = Ux;
    utilda(3 * gn -1) = Uy;
    utilda(3 * gn) = Uz;

    if plotmesh
        disp('Plotting mesh');
        xp = zeros(mno, 3);

        for j = 1:mno;
            xp(j, 1:3) = xnodes(j, 1:3) ...
                + 1e3 * utilda(3 * j - 2:3 * j)'; %+ 0e4*uhat(3*j-2:3*j)';
        end

        %     amag=1;
        %     xp = amag*xp;
        figure; clf; hold on; view(0, 0)
        xlabel('x'); ylabel('y'); zlabel('z')
        style = '-k';

        for p = 1:mel
            %      plot3(xp(:,1),xp(:,2),xp(:,3),'.') % plot nodes
            % plot elements
            plot3(xp(nc(p, [1:4, 1]), 1), xp(nc(p, [1:4, 1]), 2), xp(nc(p, [1:4, 1]), 3), style)
            plot3(xp(nc(p, [5:8, 5]), 1), xp(nc(p, [5:8, 5]), 2), xp(nc(p, [5:8, 5]), 3), style)
            plot3(xp(nc(p, [1, 5]), 1), xp(nc(p, [1, 5]), 2), xp(nc(p, [1, 5]), 3), style)%
            plot3(xp(nc(p, [2, 6]), 1), xp(nc(p, [2, 6]), 2), xp(nc(p, [2, 6]), 3), style)%
            plot3(xp(nc(p, [3, 7]), 1), xp(nc(p, [3, 7]), 2), xp(nc(p, [3, 7]), 3), style)%
            plot3(xp(nc(p, [4, 8]), 1), xp(nc(p, [4, 8]), 2), xp(nc(p, [4, 8]), 3), style)%
        end

        axis equal
        zlabel('z-direction (\mum)');
        xlabel('x-direction (\mum)');
        %     zlim([-6 6])
        %     xlim([0 9e4])
        title('$\tilde{u}$ scaled', 'FontSize', 14, 'Interpreter', 'Latex');
    end

    %-------------------generate stress--------------------------
    if plotstress
        X = linspace(0, dx, 2 * mx)';
        Z = linspace(0, dy, 2 * my)';
        Y = 0.5 * dy; %middle slide

        X_size = length(X);
        Y_size = length(Y);
        Z_size = length(Z);

        Sxxu = zeros(X_size, Y_size);
        Syyu = zeros(X_size, Y_size);
        Szzu = zeros(X_size, Y_size);
        Sxyu = zeros(X_size, Y_size);
        Sxzu = zeros(X_size, Y_size);
        Syzu = zeros(X_size, Y_size);
        Sxx = zeros(X_size, Y_size);
        Syy = zeros(X_size, Y_size);
        Szz = zeros(X_size, Y_size);
        Sxy = zeros(X_size, Y_size);
        Sxz = zeros(X_size, Y_size);
        Syz = zeros(X_size, Y_size);
        p1x = segments(:, 6);
        p1y = segments(:, 7);
        p1z = segments(:, 8);
        p2x = segments(:, 9);
        p2y = segments(:, 10);
        p2z = segments(:, 11);
        bx = segments(:, 3);
        by = segments(:, 4);
        bz = segments(:, 5);
        x1 = [p1x, p1y, p1z];
        x2 = [p2x, p2y, p2z];
        b = [bx, by, bz];

        disp('Calculating stresses');

        for i = 1:X_size;

            for j = 1:Z_size
                x0 = [X(i) Y Z(j)]; % field point
                sigmahat = hatStress(uhat, nc, xnodes, D, mx, mz, w, h, d, x0);
                Sxxu(i, j) = sigmahat(1, 1);
                Syyu(i, j) = sigmahat(2, 2);
                Szzu(i, j) = sigmahat(3, 3);

                Sxyu(i, j) = sigmahat(1, 2); %isotropic
                Sxzu(i, j) = sigmahat(1, 3); %isotropic
                Syzu(i, j) = sigmahat(2, 3); %isotropic

                sigmatilde = FieldPointStress(x0, x1, x2, b, a, MU, NU);
                Sxx(i, j) = sigmatilde(1);
                Syy(i, j) = sigmatilde(2);
                Szz(i, j) = sigmatilde(3);
                Sxy(i, j) = sigmatilde(4);
                Syz(i, j) = sigmatilde(5);
                Sxz(i, j) = sigmatilde(6);

            end

        end

        figure; clf
        subplot(3, 1, 1)
        surf(X * amag, Z * amag, mumag * (Sxxu + Sxx)', 'EdgeColor', 'none');
        view(2)
        axis equal;
        axis([0 dx * amag 0 dz * amag])
        xlabel('x-direction (\mum)')
        ylabel('z-direction (\mum)')
        title('$$\hat{\sigma}_{xx}$$+$$\tilde{\sigma}_{xx}$$', 'Interpreter', 'Latex');
        grid off
        h = colorbar;
        xlabel(h, 'MPa');

        subplot(3, 1, 2)
        surf(X * amag, Z * amag, mumag * Sxx', 'EdgeColor', 'none');
        view(2)
        axis equal;
        axis([0 dx * amag 0 dz * amag])
        xlabel('x-direction (\mum)')
        ylabel('z-direction (\mum)')
        title('$$\tilde{\sigma}_{xx}$$', 'Interpreter', 'Latex');
        grid off
        h = colorbar;
        xlabel(h, 'MPa');

        subplot(3, 1, 3)
        surf(X * amag, Z * amag, mumag * Sxxu', 'EdgeColor', 'none');
        view(2)
        axis equal;
        axis([0 dx * amag 0 dz * amag])
        xlabel('x-direction (\mum)')
        ylabel('z-direction (\mum)')
        title('$$\hat{\sigma}_{xx}$$', 'Interpreter', 'Latex');
        grid off
        %saveas(gcf,'sxx','epsc')
    end

end

%%
function tractionCheck
    %load a test condition
    run inputCheck.m

    segments = constructsegmentlist(rn, links);

    %construct finite element arrays
    %construct stiffeness matrix K and pre-compute L,U decompositions.
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
    [B, xnodes, mno, nc, n, D, kg, K, L, U, Sleft, Sright, Stop, Sbot, ...
            Sfront, Sback, gammat, gammau, gammaMixed, fixedDofs, freeDofs, ...
            w, h, d, my, mz, mel] = finiteElement3D(dx, dy, dz, mx, MU, NU, loading);

    X = linspace(0, dx, 2 * mx)';
    Z = linspace(0, dy, 2 * my)';
    Y = 0.5 * dy;

    X_size = length(X);
    Y_size = length(Y);
    Z_size = length(Z);

    Sxxuhat = zeros(X_size, Z_size);
    Syyuhat = zeros(X_size, Z_size);
    Szzuhat = zeros(X_size, Z_size);
    Sxyuhat = zeros(X_size, Z_size);
    Sxzuhat = zeros(X_size, Z_size);
    Syzuhat = zeros(X_size, Z_size);

    Sxxutilde = zeros(X_size, Z_size);
    Syyutilde = zeros(X_size, Z_size);
    Szzutilde = zeros(X_size, Z_size);
    Sxyutilde = zeros(X_size, Z_size);
    Sxzutilde = zeros(X_size, Z_size);
    Syzutilde = zeros(X_size, Z_size);

    Sxx = zeros(X_size, Z_size);
    Syy = zeros(X_size, Z_size);
    Szz = zeros(X_size, Z_size);
    Sxy = zeros(X_size, Z_size);
    Sxz = zeros(X_size, Z_size);
    Syz = zeros(X_size, Z_size);

    %%
    %---------------------------------------------------------------
    % evaluate tractions on gammat and solve for uhat  to recover
    % dislocation stress field

    gn = 1:mno; %gammau(:,1);% global node number
    x0 = xnodes(gn, 1:3); % field point
    point_array_length = size(x0, 1);
    segments_array_length = size(segments, 1);

    [Ux, Uy, Uz] = UtildaMex(x0(:, 1), x0(:, 2), x0(:, 3), ...%coordinates
    segments(:, 3), segments(:, 4), segments(:, 5), ...%burgers vector
    segments(:, 6), segments(:, 7), segments(:, 8), ...%start node segs
    segments(:, 9), segments(:, 10), segments(:, 11), ...%end node segs
    segments(:, 12), segments(:, 13), segments(:, 14), ...%slip plane
    NU, point_array_length, segments_array_length);

    utilda(3 * gn - 2) = Ux;
    utilda(3 * gn - 1) = Uy;
    utilda(3 * gn) = Uz;

    %% ---------------------------------------------------------------
    for i = 1:X_size;

        for j = 1:Z_size

            x0 = [X(i) Y Z(j)]; % field point

            %sigmahatFEM = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
            %Sxxuhat(i,j) = sigmahatFEM(1,1);
            %Syyuhat(i,j) = sigmahatFEM(2,2);
            %Szzuhat(i,j) = sigmahatFEM(3,3);
            %Sxyuhat(i,j) = sigmahatFEM(1,2); %isotropic
            %Sxzuhat(i,j) = sigmahatFEM(1,3); %isotropic
            %Syzuhat(i,j) = sigmahatFEM(2,3); %isotropic

            sigmatildeFEM = hatStress(utilda, nc, xnodes, D, mx, mz, w, h, d, x0);
            Sxxutilde(i, j) = sigmatildeFEM(1, 1);
            Syyutilde(i, j) = sigmatildeFEM(2, 2);
            Szzutilde(i, j) = sigmatildeFEM(3, 3);
            Sxyutilde(i, j) = sigmatildeFEM(1, 2); %isotropic
            Sxzutilde(i, j) = sigmatildeFEM(1, 3); %isotropic
            Syzutilde(i, j) = sigmatildeFEM(2, 3); %isotropic

            x1 = segments(:, 6:8);
            x2 = segments(:, 9:11);
            b = segments(:, 3:5);
            sigmatilde = FieldPointStress(x0, x1, x2, b, a, MU, NU);
            Sxx(i, j) = sigmatilde(1);
            Syy(i, j) = sigmatilde(2);
            Szz(i, j) = sigmatilde(3);
            Sxy(i, j) = sigmatilde(4);
            Syz(i, j) = sigmatilde(5);
            Sxz(i, j) = sigmatilde(6);
        end

    end

    %%
    subplot(5, 2, 1)
    contourf(X, Z, Sxx');
    axis([0 dx / 2 0 dz])
    axis equal
    xlabel('x-direction', 'FontSize', 14)
    %ylabel('z-direction','FontSize',14)
    title(['\sigma_{xx}^{Analytical}'], 'FontSize', 14)
    caxis([-5e-45e-4]);
    subplot(5, 2, 2)
    contourf(X, Z, Sxxutilde');
    axis([0 dx / 2 0 dz])
    axis equal;
    xlabel('x-direction', 'FontSize', 14)
    ylabel('z-direction', 'FontSize', 14)
    title(['\sigma_{xx}^{FEM}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 3)
    contourf(X, Z, Syy');
    axis([0 dx / 2 0 dy])
    axis equal
    xlabel('x-direction', 'FontSize', 14)
    %ylabel('z-direction','FontSize',14)
    title(['\sigma_{yy}^{Analytical}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 4)
    %c=1E-4*linspace(-1,1,100);
    % contourf(X,Z,Sxx',c,'EdgeColor','none');
    contourf(X, Z, Syyutilde');
    axis([0 dx / 2 0 dy])
    axis equal;
    xlabel('x-direction', 'FontSize', 14)
    ylabel('z-direction', 'FontSize', 14)
    title(['\sigma_{yy}^{FEM}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 5)
    contourf(X, Z, Szz');
    axis([0 dx / 2 0 dy])
    axis equal
    xlabel('x-direction', 'FontSize', 14)
    %ylabel('z-direction','FontSize',14)
    title(['\sigma_{zz}^{Analytical}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 6)
    %c=1E-4*linspace(-1,1,100);
    contourf(X, Z, Szzutilde');
    %surf(X,Z,Szzutilde','EdgeColor','none');
    axis([0 dx / 2 0 dy])
    axis equal;
    xlabel('x-direction', 'FontSize', 14)
    ylabel('z-direction', 'FontSize', 14)
    title(['\sigma_{zz}^{FEM}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 7)
    contourf(X, Z, Sxy');
    axis([0 dx / 2 0 dy])
    axis equal
    xlabel('x-direction', 'FontSize', 14)
    %ylabel('z-direction','FontSize',14)
    title(['\sigma_{xy}^{Analytical}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 8)
    %c=1E-4*linspace(-1,1,100);
    % contourf(X,Z,Sxx',c,'EdgeColor','none');
    contourf(X, Z, Sxyutilde');
    axis([0 dx / 2 0 dy])
    axis equal;
    xlabel('x-direction', 'FontSize', 14)
    ylabel('z-direction', 'FontSize', 14)
    title(['\sigma_{xy}^{FEM}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 9)
    contourf(X, Z, Sxz');
    axis([0 dx / 2 0 dy])
    axis equal
    xlabel('x-direction', 'FontSize', 14)
    %ylabel('z-direction','FontSize',14)
    title(['\sigma_{xz}^{Analytical}'], 'FontSize', 14)
    caxis([-5e-45e-4]);

    subplot(5, 2, 10)
    %c=1E-4*linspace(-1,1,100);
    % contourf(X,Z,Sxx',c,'EdgeColor','none');
    contourf(X, Z, Sxzutilde');
    axis([0 dx / 2 0 dy])
    axis equal;
    xlabel('x-direction', 'FontSize', 14)
    ylabel('z-direction', 'FontSize', 14)
    title(['\sigma_{xz}^{FEM}'], 'FontSize', 14)
    caxis([-5e-45e-4]);
end

%%
function displacement_test

    segments = constructsegmentlist(rn, links);

    x0 = [linspace(0, dx)', 2.5e4 * ones(100, 1), 4.5e4 * ones(100, 1)];

    tic;
    [Ux_f, Uy_f, Uz_f] = displacement_fivel(x0, segments, NU);
    toc;

    point_array_length = size(x0, 1);
    segments_array_length = size(segments, 1);
    tic;
    [Ux, Uy, Uz] = UtildaMex(x0(:, 1), x0(:, 2), x0(:, 3), ...%coordinates
    segments(:, 3), segments(:, 4), segments(:, 5), ...%burgers vector
    segments(:, 6), segments(:, 7), segments(:, 8), ...%start node segs
    segments(:, 9), segments(:, 10), segments(:, 11), ...%end node segs
    segments(:, 12), segments(:, 13), segments(:, 14), ...%slip plane
    NU, point_array_length, segments_array_length);
    utilda2 = horzcat(Ux, Uy, Uz);
    toc;

    norm([Ux - Ux_f, Uy - Uy_f, Uz - Uz_f]);

    clear utilda3

    for i = 1:size(x0, 1)

        uspf = displacement_spf(x0(i, :), rn(:, 1:3), segments(1, 3:5), NU);
        utilda3(i, 1:3) = uspf;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether virtual nodes can be eliminated
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% February 2014
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = virtualmeshcoarsen(rn, links, connectivity, linksinconnect, fseg, lmin, lmax, MU, NU, a, Ec)
    rnnew = rn;
    %[lrn, lrn2]=size(rn);
    areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
    linksnew = links;
    fsegnew = fseg;
    connectivitynew = connectivity;
    linksinconnectnew = linksinconnect;
    areamin2 = areamin * areamin;
    delta = 1e-16;
    i = 1;

    while i <= length(rnnew(:, 1))%the loop looks each nodes M
        %fprintf('%i \n',i);
        if (connectivitynew(i, 1) == 2) && rnnew(i, end) == 67%virtual nodes
            % the discretization node is normal so set up the conditions to check for whether is should be coarsened away
            % looking for the environnment of 'node i': 2 links and 2 other nodes M.
            link1 = connectivitynew(i, 2);
            link2 = connectivitynew(i, 4);
            posi_link1 = connectivitynew(i, 3);
            posi_link2 = connectivitynew(i, 5);
            posnoti_link1 = 3 - posi_link1;
            posnoti_link2 = 3 - posi_link2;
            link1_nodenoti = linksnew(link1, posnoti_link1);
            link2_nodenoti = linksnew(link2, posnoti_link2);

            if rnnew(link1_nodenoti, end) ~= 67 || rnnew(link2_nodenoti, end) ~= 67
                % if neighbour nodes are not virtual nodes then continue M
                % if at least one neighbour node is a virtual node : see after the loop
                i = i + 1;
                continue;
            end

            vec1 = rnnew(link1_nodenoti, 1:3) - rnnew(i, 1:3);
            vec2 = rnnew(link2_nodenoti, 1:3) - rnnew(i, 1:3);
            vec3 = vec2 - vec1;
            r1 = sqrt(vec1 * vec1');
            r2 = sqrt(vec2 * vec2');
            r3 = sqrt(vec3 * vec3');

            if r3 < lmax%if the coarsening would result in a link length larger than lmax don't coarsen
                s = 0.5 * (r1 + r2 + r3);
                area2 = (s * (s - r1) * (s - r2) * (s - r3));
                dvec1dt = rnnew(link1_nodenoti, 4:6) - rnnew(i, 4:6);
                dvec2dt = rnnew(link2_nodenoti, 4:6) - rnnew(i, 4:6);
                dvec3dt = dvec2dt - dvec1dt;
                dr1dt = (vec1 * dvec1dt') / (r1 + delta);
                dr2dt = (vec2 * dvec2dt') / (r2 + delta);
                dr3dt = (vec3 * dvec3dt') / (r3 + delta);
                dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);
                darea2dt = dsdt * (s - r1) * (s - r2) * (s - r3);
                darea2dt = darea2dt + s * (dsdt - dr1dt) * (s - r2) * (s - r3);
                darea2dt = darea2dt + s * (s - r1) * (dsdt - dr2dt) * (s - r3);
                darea2dt = darea2dt + s * (s - r1) * (s - r2) * (dsdt - dr3dt);

                if ((area2 < areamin2) && (darea2dt < 0.0d0)) || ((r1 < lmin) || (r2 < lmin))%remove single node critierion
                    %the area is less than minimum and shrinking or one of the arms is less than lmin and shrinking
                    [rnnew, connectivitynew, linksnew, linksinconnectnew, fsegnew, ~] = mergenodes(rnnew, connectivitynew, linksnew, linksinconnectnew, fsegnew, link2_nodenoti, i, MU, NU, a, Ec);
                else %((area2>=areamin2)|(dareadt>=0.0d0)) & ((r1>=lmin)|(dr1dt>=lmin)) & ((r2>=lmin)|(dr2dt>=lmin))
                    i = i + 1;
                end

            else % r3>=lmax
                i = i + 1;
            end

        else % connectivitynew(i,1)>2
            i = i + 1;
        end

    end %while loop

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether virtual nodes can be eliminated based on:
% 1) Segment length (i.e. between lmin lmax)
% Simpler/faster than virtualmeshcoarsen.m
%
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% February 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew, linksnew, connectivitynew, linksinconnectnew] = virtualmeshcoarsen2(rn, links, maxconnections, lmin)

    %| rn(loop_list(:,1),4)==6
    % N.B. loop_list gives a list of nodes that below in a loop, in sequence.
    % However, it does not give the "direction" of this sequence since segments
    % are vectors, hence you need to check this based on input links array.

    lrn = size(rn, 1); %total number of nodes M
    [connectivity, linksinconnect, ~] = genconnectivity(rn, links, maxconnections);
    counter = 1;
    fseg = []; %empty

    while counter > 0;
        counter = 0;
        i = 1;

        while i < lrn% i<total numbers of nodes M

            if rn(i, end) ~= 67 || connectivity(i, 1) ~= 2
                %ignore domain/fixed/surface nodes
                %ignore junctions
                i = i + 1;
                continue;
            end

            node1 = links(connectivity(i, 2), logicswap(connectivity(i, 3)));
            node2 = links(connectivity(i, 4), logicswap(connectivity(i, 5)));
            seg1 = rn(i, 1:3) - rn(node1, 1:3);
            seg2 = rn(i, 1:3) - rn(node2, 1:3);
            length = sqrt(sum(seg1 .* seg1)) + sqrt(sum(seg2 .* seg2)); %length=seg1+seg2 ? M

            if length < lmin
                %check angle
                %angle=acos(dot(seg1,seg2)/norm(seg1)/norm(seg2));
                %if angle < (pi/10) || angle > (9*pi/10) %if the curvature is very high, then don't delete.
                % remove node and remesh
                counter = counter + 1;
                [rn, connectivity, links, linksinconnect, ~, ~] = ...
                    mergenodes(rn, connectivity, links, linksinconnect, fseg, node1, i);
                %end
            end

            lrn = size(rn, 1);
            i = i + 1;
        end

        %clean-up orphan nodes flagged as superfluous and re-index nodal structure.
        %[rn,links] = clearorphannodes(rn,links);
        %update connectivity and linksinconnect.
        %[~,~,~] = genconnectivity(rn,links,maxconnections);
    end

    %clean-up orphan nodes flagged as superfluous and re-index nodal structure.
    [rnnew, linksnew] = clearorphannodes(rn, links);
    %update connectivity and linksinconnect.
    [connectivitynew, linksinconnectnew] = genconnectivity(rn, links, maxconnections);

end

%%%%%%%%

function [rnnew, linksnew] = clearorphannodes(rnnew, linksnew)

    orphan_nodes = rnnew(:, end) == -1;

    while any(orphan_nodes)

        for i = 1:size(rnnew, 1)

            if rnnew(i, 4) == -1
                [row, col] = find(linksnew(:, 1:2) == i);
                %check it's a pair!
                if length(row) ~= 2
                    disp('WARNING! Not analyzing pairs! See virtualmeshcoarsen.m (row 79)')
                end

                if col(1) == 2
                    %merge
                    linksnew(row(1), 2) = linksnew(row(2), 2);
                    linksnew(row(2), :) = [];
                elseif col(1) == 1
                    %merge
                    linksnew(row(1), 1) = linksnew(row(2), 1);
                    linksnew(row(2), :) = [];
                end

                rnnew(i, end) = 1; %change flag
                break;
            end

        end

        %update orphan_nodes list
        orphan_nodes = (rnnew(:, end) == -1);
    end

    orphan_nodes = (rnnew(:, end) == 1); %look for flagged rn.
    index_corr = cumsum(orphan_nodes);
    linksnew(:, 1) = linksnew(:, 1) - index_corr(linksnew(:, 1)); %update numbering
    linksnew(:, 2) = linksnew(:, 2) - index_corr(linksnew(:, 2)); %update numbering
    rnnew(orphan_nodes, :) = []; %finally get rid of nodes

end

%%%%%%%%

function out = logicswap(in)

    if in == 1
        out = 2;
    elseif in == 2;
        out = 1;
    else
        disp('Incorrect use of logicswap function');
    end

end

%%%%%%%%

function [rn, links, connectivity, linksinconnect, fseg] = collision(rn, links, connectivity, linksinconnect, fseg, mindist, MU, NU, a, Ec, mobility, vertices, ...
        uhat, nc, xnodes, D, mx, mz, w, h, d)
    % this subroutine goes through the existing links and checks for collisions
    % it first checks through unconnected links
    % it then checks for hinges that are coming within the mininum distance
    mindist2 = mindist * mindist;
    lrn2 = length(rn(1, :));
    lrn3 = lrn2 - 1;
    % eps is the error factor of the calculation
    eps = 1e-12;
    % check for two links colliding
    i = 1;

    while i < size(links, 1)%for every link M
        %ignore virtual segments - FF
        if rn(links(i, 1), end) == 67 || rn(links(i, 2), end) == 67
            i = i + 1;
            continue;
        end

        j = i + 1;

        while j <= size(links, 1)
            n1s1 = links(i, 1);

            n2s1 = links(i, 2);
            n1s2 = links(j, 1);
            n2s2 = links(j, 2);

            if (n1s1 ~= n1s2) && (n1s1 ~= n2s2) && (n2s1 ~= n1s2) && (n2s1 ~= n2s2)

                %[dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(n1s2,1:lrn3),rn(n2s2,1:lrn3));
                [dist2, ddist2dt, L1, L2] = MinDistCalcMex(rn(n1s1, 1:lrn3), rn(n2s1, 1:lrn3), rn(n1s2, 1:lrn3), rn(n2s2, 1:lrn3));
                collision_condition_is_met = ((dist2 < mindist2) & (ddist2dt <- eps)) | (dist2 < eps);
                % there are two conditions here the first condition handles non planar collisions
                % the second conditions looks for coplanar collisions
                if collision_condition_is_met
                    % links are unconnected and colliding
                    % identify the first node to be merged
                    vec = rn(n1s1, 1:3) - rn(n2s1, 1:3);
                    close_to_n1s1 = ((L1 * L1 * (vec * vec')) < mindist2);
                    close_to_n2s1 = (((1 - L1) * (1 - L1) * (vec * vec')) < mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    if close_to_n1s1
                        mergenode1 = n1s1;
                    elseif close_to_n2s1
                        mergenode1 = n2s1;
                    else
                        spnode = n1s1;
                        splitconnection = linksinconnect(i, 1); % linki=s1 M
                        posvel = rn(n1s1, 1:lrn3) .* (1 - L1) + rn(n2s1, 1:lrn3) .* L1;
                        [rn, links, connectivity, linksinconnect] = splitnode(rn, links, connectivity, linksinconnect, spnode, splitconnection, posvel);
                        mergenode1 = length(rn(:, 1)); %nodeid of mergenode1 M
                        linknew = length(links(:, 1)); %linkid of linknew M
                        links(linknew, 6:8) = links(i, 6:8); %glide plane M
                        fseg = [fseg; zeros(1, 6)];
                    end

                    % identify the second node to be merged
                    vec = rn(n1s2, 1:3) - rn(n2s2, 1:3);
                    close_to_n1s2 = ((L2 * L2 * (vec * vec')) < mindist2);
                    close_to_n2s2 = (((1 - L2) * (1 - L2) * (vec * vec')) < mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    if close_to_n1s2
                        mergenode2 = n1s2;
                    elseif close_to_n2s2
                        mergenode2 = n2s2;
                    else
                        spnode = n1s2;
                        splitconnection = linksinconnect(j, 1); %linkj=s2 M
                        posvel = rn(n1s2, 1:lrn3) .* (1 - L2) + rn(n2s2, 1:lrn3) .* L2;
                        [rn, links, connectivity, linksinconnect] = splitnode(rn, links, connectivity, linksinconnect, spnode, splitconnection, posvel);
                        mergenode2 = length(rn(:, 1));
                        linknew = length(links(:, 1));
                        links(linknew, 6:8) = links(j, 6:8);
                        fseg = [fseg; zeros(1, 6)];
                    end

                    % merge the two colliding nodes
                    %                 disp(sprintf('node %d and node %d are colliding by two line collision',mergenode1,mergenode2))
                    collisionpoint = findcollisionpoint(mergenode1, mergenode2, rn, connectivity, links);
                    rn(mergenode1, 1:lrn2) = [collisionpoint 0 0 0 max(rn(mergenode1, lrn2), rn(mergenode2, lrn2))];
                    [rn, connectivity, links, linksinconnect, fseg, mergednodeid] = mergenodes(rn, connectivity, links, linksinconnect, fseg, mergenode1, mergenode2, MU, NU, a, Ec);

                    if mergednodeid > 0

                        for k = 1:connectivity(mergednodeid, 1)
                            linkid = connectivity(mergednodeid, 2 * k);
                            fseg(linkid, :) = segforcevec(MU, NU, a, Ec, rn(:, [1 2 3 lrn2]), links, linkid, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
                            othernode = links(linkid, 3 - connectivity(mergednodeid, 2 * k + 1)); % 3-connectivity(mergednodeid,2*k+1) = 1 or 2, it corresponds to the position of the other node of the link ( the one which is not mergenode ) M
                            clist = [connectivity(othernode, 1) linspace(1, connectivity(othernode, 1), connectivity(othernode, 1))];
                            [rn(othernode, 4:6), ~] = feval(mobility, fseg, rn, links, connectivity, othernode, clist);
                        end

                        numbcon = connectivity(mergednodeid, 1);
                        conlist = [numbcon linspace(1, numbcon, numbcon)];
                        [rn(mergednodeid, 4:6), ~] = feval(mobility, fseg, rn, links, connectivity, mergednodeid, conlist);
                    end

                end

            end

            j = j + 1;
        end

        i = i + 1;
    end

    % check for a hinge condition
    i = 1;

    while i <= length(rn(:, 1))
        %ignore virtual nodes - FF
        if rn(i, end) == 67
            i = i + 1;
            continue;
        end

        j = 1;

        while j <= connectivity(i, 1)
            nodenoti = links(connectivity(i, 2 * j), 3 - connectivity(i, 2 * j + 1));
            k = 1;

            while k <= connectivity(i, 1)
                linkid = connectivity(i, 2 * k);
                % if node is on the link do not check for collision
                if j ~= k
                    % identify the nodes on the link
                    n1s1 = links(linkid, 1);
                    n2s1 = links(linkid, 2);
                    %[dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(nodenoti,1:lrn3),rn(nodenoti,1:lrn3));
                    [dist2, ddist2dt, L1, ~] = MinDistCalcMex(rn(n1s1, 1:lrn3), rn(n2s1, 1:lrn3), rn(nodenoti, 1:lrn3), rn(nodenoti, 1:lrn3));
                    %dist2
                    %ddist2dt
                    collision_condition_is_met = (dist2 < mindist2) & (ddist2dt <- eps);

                    if collision_condition_is_met
                        % identify the first node to be merged
                        mergenode1 = nodenoti;
                        % identify the second node to be merged
                        vec = rn(n1s1, 1:3) - rn(n2s1, 1:3);
                        close_to_n1s1 = ((L1 * L1 * (vec * vec')) < mindist2);
                        close_to_n2s1 = (((1 - L1) * (1 - L1) * (vec * vec')) < mindist2);
                        % if collision point is close to one of the existing nodes use that node
                        if close_to_n1s1
                            mergenode2 = n1s1;
                        elseif close_to_n2s1
                            mergenode2 = n2s1;
                        else
                            spnode = n1s1;
                            splitconnection = linksinconnect(linkid, 1);
                            posvel = rn(n1s1, 1:lrn3) .* (1 - L1) + rn(n2s1, 1:lrn3) .* L1;
                            [rn, links, connectivity, linksinconnect] = splitnode(rn, links, connectivity, linksinconnect, spnode, splitconnection, posvel);
                            mergenode2 = length(rn(:, 1));
                            newlink = length(links(:, 1));
                            links(newlink, 6:8) = links(linkid, 6:8);
                            fseg = [fseg; zeros(1, 6)];
                        end

                        %merge the two nodes
                        %                     disp(sprintf('node %d and node %d are colliding by hinge condition',mergenode2,mergenode1))
                        collisionpoint = findcollisionpoint(mergenode1, mergenode2, rn, connectivity, links);
                        rn(mergenode1, 1:lrn2) = [collisionpoint 0 0 0 max(rn(mergenode1, lrn2), rn(mergenode2, lrn2))];
                        [rn, connectivity, links, linksinconnect, fseg, mergednodeid] = mergenodes(rn, connectivity, links, linksinconnect, fseg, mergenode1, mergenode2, MU, NU, a, Ec);

                        if mergednodeid > 0

                            for k = 1:connectivity(mergednodeid, 1)
                                linkid = connectivity(mergednodeid, 2 * k);
                                fseg(linkid, :) = segforcevec(MU, NU, a, Ec, rn(:, [1 2 3 lrn2]), links, linkid, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
                                othernode = links(linkid, 3 - connectivity(mergednodeid, 2 * k + 1));
                                clist = [connectivity(othernode, 1) linspace(1, connectivity(othernode, 1), connectivity(othernode, 1))];
                                [rn(othernode, 4:6), ~] = feval(mobility, fseg, rn, links, connectivity, othernode, clist);
                            end

                            numbcon = connectivity(mergednodeid, 1);
                            conlist = [numbcon linspace(1, numbcon, numbcon)];
                            [rn(mergednodeid, 4:6), ~] = feval(mobility, fseg, rn, links, connectivity, mergednodeid, conlist);
                        end

                        %there has been a connectivity change in node i start the search through node i's connections from the beginning
                        if i > size(rn, 1)
                            % this is a rare but possible case.
                            % for this condition to be satisfied the last node was being checked for closed hinge condition and it merged with another node
                            % since this was the last node being checked exit the function
                            return;
                        else
                            j = 0;
                            k = connectivity(i, 1);
                        end

                    end

                end

                k = k + 1;
            end

            j = j + 1;
        end

        i = i + 1;
    end

end

% function [dist2,ddist2dt,L1,L2]=mindistcalc(x0vx0,x1vx1,y0vy0,y1vy1)
% % this function finds the minimum distance bewtween two line segments
% % seg1=x0->x1 seg2=y0->y1
% % dist2 = square of the minimum distance between the two points
% % L1 = normalize position on seg1 that is closest to seg2
% % L2 = normalized position on seg2 that is closest to seg1
% % ddist2dt = time rate of change of the distance between L1 and L2
% x0=x0vx0(1:3); %seg1 M
% x1=x1vx1(1:3); %seg1 M
% y0=y0vy0(1:3); %seg2 M
% y1=y1vy1(1:3); %seg2 M
% if length(x0vx0)==6 %if there are velocities M
%     vx0=x0vx0(4:6);
%     vx1=x1vx1(4:6);
%     vy0=y0vy0(4:6);
%     vy1=y1vy1(4:6);
% else %if there are no velocities M
%     vx1=zeros(1,3);
%     vx0=zeros(1,3);
%     vy1=zeros(1,3);
%     vy0=zeros(1,3);
% end
%
% seg1=x1-x0; %vector seg1 M
% seg2=y1-y0; %vector seg2 M
% vseg1=vx1-vx0; %vector velocity seg1 M
% vseg2=vy1-vy0; %vector velocity seg2 M
%
% A=seg1*seg1'; %��seg1��^2 M
% B=2*seg1*(x0'-y0');
% C=2*seg1*seg2';
% D=2*seg2*(y0'-x0');
% E=seg2*seg2'; %��seg2��^2 M
% F=x0*x0'+y0*y0';
% G=C*C-4*A*E;
% eps=1e-12;
% if A<eps % seg1 is just a point
%     L1=0;
%     if E<eps
%         L2=0;
%     else
%         L2=-0.5*D/E;
%     end
% elseif E<eps % seg2 is just a point
%     L2=0;
%     if A<eps
%         L1=0;
%     else
%         L1=-0.5*B/A;
%     end
% elseif abs(G)<eps % lines are parallel
%     dist2=[(y0-x0)*(y0-x0)' (y1-x0)*(y1-x0)' (y0-x1)*(y0-x1)' (y1-x1)*(y1-x1)'];
%     [mindist2,pos]=min(dist2);
%     L1=floor(pos/2);
%     L2=mod(pos-1,2);
% else
%     L2=(2*A*D+B*C)/G;
%     L1=0.5*(C*L2-B)/A;
% end
%
% % now check to make sure that L2 and L1 are betwen 0 and 1
% L1=min(max([L1,0]),1);
% L2=min(max([L2,0]),1);
%
% % now calculate the distance^2 and the time rate of change of the distance between the points at L1 and L2
% dist2=(x0+seg1.*L1-y0-seg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)';
% ddist2dt=2*((vx0+vseg1.*L1-vy0-vseg2.*L2)*(x0+seg1.*L1-y0-seg2.*L2)');

function collisionpoint = findcollisionpoint(mergenode1, mergenode2, rn, connectivity, links)
    % this subroutine finds the collision point of two nodes given that there are strict glide plane constraints
    eps = 1e-12;
    newplanecondition = 0.875;
    p1 = rn(mergenode1, 1:3);
    p2 = rn(mergenode2, 1:3);
    Nmat = zeros(3, 3);
    Nsize = 0;
    vector = zeros(3, 1);
    s = size(rn, 2);

    if rn(mergenode1, s) == 7
        collisionpoint = rn(mergenode1, 1:3);
        return;
    elseif rn(mergenode2, s) == 7
        collisionpoint = rn(mergenode2, 1:3);
        return;
    end

    for i = 1:connectivity(mergenode1, 1)

        if Nsize < 3
            linkid = connectivity(mergenode1, 2 * i);
            connode = links(connectivity(mergenode1, 2 * i), 3 - connectivity(mergenode1, 2 * i + 1));
            rt = rn(mergenode1, 1:3) - rn(connode, 1:3);
            L = norm(rt);
            linedir = rt ./ L;
            n1 = cross(linedir, links(linkid, 3:5));
            n2 = cross(linedir, rn(mergenode1, 4:6));

            if n1 * n1' > eps
                plane = n1 ./ norm(n1);
            elseif n2 * n2' > eps
                plane = n2 ./ norm(n2);
            end

            if ((n1 * n1' > eps) || (n2 * n2' > eps))

                if Nsize == 0
                    conditionismet = 1;
                elseif Nsize == 1
                    conditionismet = ((Nmat(1, :) * plane')^2 < newplanecondition * newplanecondition);
                else
                    detN = det([Nmat(1:2, :); plane]);
                    conditionismet = detN * detN > (1 - newplanecondition)^4;
                end

                if conditionismet
                    Nsize = Nsize + 1;
                    Nmat(Nsize, :) = plane;
                    vector(Nsize) = plane * p1';
                end

            end

        end

    end

    for i = 1:connectivity(mergenode2, 1)

        if Nsize < 3
            linkid = connectivity(mergenode2, 2 * i);
            connode = links(connectivity(mergenode2, 2 * i), 3 - connectivity(mergenode2, 2 * i + 1));
            rt = rn(mergenode2, 1:3) - rn(connode, 1:3);
            L = norm(rt);
            linedir = rt ./ L;
            n1 = cross(linedir, links(linkid, 3:5));
            n2 = cross(linedir, rn(mergenode2, 4:6));

            if n1 * n1' > eps
                plane = n1 ./ norm(n1);
            elseif n2 * n2' > eps
                plane = n2 ./ norm(n2);
            end

            if ((n1 * n1' > eps) || (n2 * n2' > eps))

                if Nsize == 1
                    conditionismet = ((Nmat(1, :) * plane')^2 < newplanecondition * newplanecondition);
                else
                    detN = det([Nmat(1:2, :); plane]);
                    conditionismet = detN * detN > (1 - newplanecondition)^4;
                end

                if conditionismet
                    Nsize = Nsize + 1;
                    Nmat(Nsize, :) = plane;
                    vector(Nsize) = plane * p2';
                end

            end

        end

    end

    Matrix = [eye(3) Nmat(1:Nsize, :)'; Nmat(1:Nsize, :) zeros(Nsize, Nsize)];
    V = [(rn(mergenode1, 1:3)' + rn(mergenode2, 1:3)') ./ 2; vector(1:Nsize)];
    res = Matrix \ V;
    collisionpoint = res(1:3)';
end

%%%%%%
function [uhat, fend, Ubar] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, dx, t, unfixedDofs, Kred, Lred, Ured, curstep, dt)

    %HY20190307: added global unloading flag and counter
    global unloadflag unloadcount Ubarglobal

    %Coupling of FEM and DDD
    % u = uhat + utilda
    % f = fhat + ftilda

    unloadbegin = 1000;

    segments = constructsegmentlist(rn, links);

    %HY20190307: modified by HY for cyclic loading

    % Udot = 100*1E3*dx*(1E-4/160E3); %for tungsten...
    % Udot =100*1E3*dx*(1E-4/mumag)*100 ; % Marielle
    Udot = dx * 1E-5; % HY20180119

    if curstep <= unloadbegin
        Ubarglobal = Udot * t;
    end

    %Ubar = 0.1*1E4; for debuggin

    % if curstep>unloadbegin
    %     unloadflag = 1;
    % end

    if unloadflag == 1
        Udot = -1.0 * dx * 1E-4;
    end

    if curstep > unloadbegin
        Ubarglobal = Ubarglobal + Udot * dt;
    end

    Ubar = Ubarglobal;

    %HY20180119
    % if Ubar>0.03*dx
    %     Ubar = Ubar*1E-2;
    % end

    u = zeros(3 * (mno), 1);
    gamma = [gammau; gammaMixed];

    u(3 * gammaMixed(:, 1)) = -Ubar; %applied displacements in z at right edge nodes

    uhat = zeros(3 * mno, 1);
    utilda = zeros(3 * mno, 1);

    gn = gamma(:, 1); % global node number
    x0 = xnodes(gn, 1:3); % field point
    point_array_length = size(x0, 1);
    segments_array_length = size(segments, 1);

    %Matlab wrapper
    % tic;
    % displacements = displacement_fivel(x0,segments,NU); %utilda array ordered by x0
    % toc;

    %Full C version (including combinatorials, checks, corrections etc.)
    % tic;
    % [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
    %                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
    %                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
    %                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
    %                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
    %                        NU,point_array_length,segments_array_length);
    % % displacements = horzcat(Ux,Uy,Uz);
    % % toc;
    % % disp(displacementsMEX-displacements);
    % %  [Ux, Uy, Uz]  = displacement_fivel(x0,segments,NU);
    %
    % utilda(3*gn -2) = Ux;
    % utilda(3*gn -1) = Uy;
    % utilda(3*gn   ) = Uz;

    if any(isnan(utilda))
        disp('some entries of utilda are NaN')
        pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
    end

    %HY20180122
    % if Ubar<dx*1E-4
    %     utilda= zeros(3*mno,1);
    % end

    uhat(fixedDofs) = u(fixedDofs) - utilda(fixedDofs);

    fhat = zeros(3 * (mno), 1);

    gamma = [gammat; gammaMixed];

    %ftilda = zeros(3*mno,1);
    ftilda = traction(gamma, segments, xnodes, mno, a, MU, NU);
    %ftilda=zeros(3*mno,1); %ET uncomment later!

    %HY20180122
    % if Ubar<dx*1E-4
    %     ftilda = zeros(3*mno,1);
    % end

    fhat(freeDofs) = -ftilda(freeDofs); % no applied forces

    %f=zeros(2*(mno),1);
    f = fhat - kg(:, fixedDofs) * uhat(fixedDofs);

    %HY20171206:********************************************************
    %HY20171206: modified by HY to make the code cleaner by removing the
    %equations related to the fixedDofs; since FreeDofs has been used to
    %represent the free boundary nodes, a new term, unfixedDofs, is used to
    %represent all the nodes other than the fixed ones. i.e.
    %unfixedDofs=allDofs - fixedDofs
    % tic;
    % disp('setdiff')
    uhat = uhat;
    %uhatHY = uhat;
    fred = f(unfixedDofs);
    u_new = Ured \ (Lred \ fred);
    uhat(unfixedDofs) = u_new;
    %uhatHY(unfixedDofs) = u_new;
    % toc;
    % pause
    %HY20171206:********************************************************
    % tic;
    % disp('none setdiff')
    % bcwt=mean(diag(kg));%=trace(K)/length(K)
    % bcwt = full(bcwt);
    %
    % f(fixedDofs) = bcwt*uhat(fixedDofs);
    % uhat = U\(L\f); %using LU decomposition
    % % uhat2=K\f;
    % toc;
    % pause

    % uhatdiff=uhatHY-uhat;
    % maxdiff=max(abs(uhatdiff))
    % pause

    rhat = kg * uhat; % reaction force

    fend = rhat(3 * gammaMixed(:, 1)) + ftilda(3 * gammaMixed(:, 1));
    fend = sum(fend);

    % fprintf('fend=%3d\n',fend);

    if fend > 0 && curstep > 1000
        unloadflag = 0
        unloadcount = unloadcount + 1;
    end

end

function [uhat, fend, Ubar] = FEMcoupler(rn, links, maxconnections, a, MU, NU, xnodes, mno, kg, L, U, ...
        gammau, gammat, gammaMixed, fixedDofs, freeDofs, unfixedDofs, dx, dy, dz, t, mx, my, mz, utilda_0)

    %Coupling of FEM and DDD
    % u = uhat + utilda
    % f = fhat + ftilda

    segments = constructsegmentlist(rn, links);

    % Udot = 1E3*dx*(1E-4/160E9); %test by BB
    Udot = 100 * 1E3 * dx * (1E-4/160E9); %for tungsten...
    % Udot =100*1E3*dx*(1E-4/160E9)*100 ; % Marielle
    % Udot = dx*1e-5; %HY

    Ubar = Udot * t;
    %Ubar = 0.1*1E4; for debuggin
    u = zeros(3 * (mno), 1);
    gamma = [gammau; gammaMixed];

    u(3 * gammaMixed(:, 1)) = -Ubar; %applied displacements in z at right edge nodes

    uhat = zeros(3 * mno, 1);
    utilda = zeros(3 * mno, 1);

    gn = gamma(:, 1); % global node number
    % x0 = xnodes(gn,1:3); % field point
    % point_array_length = size(x0,1);
    % segments_array_length = size(segments,1);

    %Matlab wrapper
    % tic;
    % displacements = displacement_fivel(x0,segments,NU); %utilda array ordered by x0
    % toc;

    %Full C version (including combinatorials, checks, corrections etc.)
    % tic;
    % [Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
    %                        segments(:,3), segments(:,4), segments(:,5),... %burgers vector
    %                        segments(:,6), segments(:,7), segments(:,8),... %start node segs
    %                        segments(:,9), segments(:,10), segments(:,11),... %end node segs
    %                        segments(:,12), segments(:,13), segments(:,14),... %slip plane
    %                        NU,point_array_length,segments_array_length);
    % displacements = horzcat(Ux,Uy,Uz);
    % toc;
    % disp(displacementsMEX-displacements);
    %  [Ux, Uy, Uz]  = displacement_fivel(x0,segments,NU);

    % % Displacement calculation Bruce Bromage
    % [Ux, Uy, Uz] = Utilda_bb3_vec(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);
    %
    % utilda(3*gn -2) = Ux;
    % utilda(3*gn -1) = Uy;
    % utilda(3*gn   ) = Uz;
    %
    % utilda = utilda - utilda_0;
    %
    % if any(isnan(utilda))
    %     disp('some entries of utilda are NaN')
    %     pause; %some entries of utilda are NaN -> leads to issues in hatStress routine
    % end
    %
    uhat(fixedDofs) = u(fixedDofs); % - utilda(fixedDofs);

    fhat = zeros(3 * (mno), 1);

    gamma = [gammat; gammaMixed];

    %ftilda = zeros(3*mno,1);
    ftilda = traction(gamma, segments, xnodes, mno, a, MU, NU);

    %%
    %ftilda=zeros(3*mno,1); %ET uncomment later!

    %% we think the plane normals are wrong, lets experiment.
    %fhat(freeDofs) = -ftilda(freeDofs);% no applied forces
    fhat(freeDofs) = ftilda(freeDofs); % no applied forces
    %%

    %f=zeros(2*(mno),1);
    f = fhat - kg(:, fixedDofs) * uhat(fixedDofs);

    %HY20171206:********************************************************
    %HY20171206: modified by HY to make the code cleaner by removing the
    %equations related to the fixedDofs; since FreeDofs has been used to
    %represent the free boundary nodes, a new term, unfixedDofs, is used to
    %represent all the nodes other than the fixed ones. i.e.
    %unfixedDofs=allDofs - fixedDofs
    %
    uhat(unfixedDofs) = U \ (L \ f(unfixedDofs));

    %HY20171206:********************************************************

    %HY20171206: commented by HY
    % bcwt=mean(diag(kg));%=trace(K)/length(K)
    % bcwt = full(bcwt);
    %
    % f(fixedDofs) = bcwt*uhat(fixedDofs);
    % uhat = U\(L\f); %using LU decomposition
    % uhat2=K\f;

    rhat = kg * uhat; % reaction force

    fend = rhat(3 * gammaMixed(:, 1)) + ftilda(3 * gammaMixed(:, 1));
    fend = sum(fend);

end

%%%%%%%%%%%%%%%%%%%%
function [rn, vn, dt, fn, fseg] = int_trapezoid_bb_old(rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d)

    %Implicit numerical integrator using the Euler-trapezoid method adapted
    %from [Cai & Bulatov, Algorithm 10.2, pg. 216]. The tiemstep size is
    %controlled so that it can increase suitably quickly and not decrease
    %excessively whilst remaining within acceptable tolerence limits
    %Written by B.Bromage and D.Celis-Garza 05/11/2018

    %Convert rn into a single column of coordinates and store the node flags
    rnvec0 = [rn(:, 1); rn(:, 2); rn(:, 3)]; flag = rn(:, 4);

    %Calculate the current nodal velocities
    [vnvec0, fn, fseg] = drndt(rnvec0, flag, MU, NU, a, Ec, links, connectivity, ...
        mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);

    maxiter = 10; %Maximum number of times the timestep can be increased
    counter = 1; %Counter variable for the while loop
    dt_old = 0; %Variable for the maximumum acceptable timestep
    maxchange = 1.2; %Maximum factor the timestep can be increased by
    exponent = 20; %Variable that controls the degree that the calculated error affects the cange in timestep
    dt_old_good = 0; %Logical which flags whether an acceptable timestep has been calculated
    convergent = 0; %Logical for the operation of the while loop

    while (~convergent)

        rnvec1 = rnvec0 + vnvec0 * dt; %Euler forward method [Cai & Bulatov, eq. 10.43]

        if isempty(rnvec1)%If there are no dislocations then use the maximum possible timestep
            convergent = 1;
        end

        %Calculate the nodal velocities for the next timestep accordin to
        %Euler forward method
        [vnvec1, fn, fseg] = drndt(rnvec1, flag, MU, NU, a, Ec, links, connectivity, ...
            mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);

        distvec = rnvec1 - rnvec0;
        distmag = max(abs(distvec)); %Largest distance any node moves

        err = distvec - (vnvec1 + vnvec0) / 2 * dt; %Euler-trapezoid method
        errmag = max(abs(err)); %Largest difference (error) calculated by the Euler-trapezoid method

        if isempty(errmag)%If the Euler-trapzoid method yields no results use maximum time step and end loop
            dt = dt0;
            break
        end

        if (errmag < rntol) && (distmag < rmax)%If error and max distance move are in acceptable limits
            dt_old = dt; %Store current timestep as maximum acceptable timestep
            factor = maxchange * (1 / (1 + (maxchange^exponent - 1) * (errmag / rntol)))^(1 / exponent);
            dt = min(dt * factor, dt0); %Increase timestep depending on the magnitude of the error
            dt_old_good = 1; %Flag acceptable timestep calculated
            counter = counter + 1; %Proceed to next iteration
        else

            if dt_old_good == 1%If current timestep is too large, use largest acceptable timestep
                dt = dt_old;
                counter = maxiter;
            else
                dt = dt / 2; %If no acceptable timestep has been calculated, halve timestep and try again
            end

        end

        if counter > maxiter || dt == dt0%End loop if maximum number of iterations is reached or curren timestep is maximum timestep
            convergent = 1;
        end

    end

    %Rearrange rnvec and vnvec into rn and vn matrices
    rn = [reshape(rnvec1, length(rnvec1) / 3, 3), flag];
    vn = reshape(vnvec1, length(vnvec1) / 3, 3);
end

%%%%%%%%%%%%%%%%%%%%
function [rn, vn, dt, fn, fseg] = int_trapezoid(rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d)
    %implicit numerical integrator using the Trapezoid method
    %dt: suggested timestep (usually from previous iteration)
    %dt0: maximum allowed timestep

    %dummy variable
    t = 0;
    rnold = rn;

    %scramble rn into a single column vector
    rnvec = [rn(:, 1); rn(:, 2); rn(:, 3)]; flag = rn(:, 4);

    %Backward Euler
    rnvec0 = rnvec;

    [vnvec0, fn, fseg] = drndt(rnvec0, flag, MU, NU, a, Ec, links, connectivity, ...
        mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
    %dt=1/max(1/dt0,max(vnvec0)/rmax);
    %dt=dt0;

    dt1 = dt;
    maxiter = 1;
    convergent = 0;

    % This algorithm [Cai & Bulatov, Algorithm 10.2, pg. 216] is a simple
    % routine to dynamically change the time-step. Whenever the difference
    % between the predicted and corrected positions for any node exceeds
    % threshold err, the timestep is decreased by half and new positions for
    % all nodes are recomputed with the reduced time step.
    %   1. Initialize time step dt = dt_max
    %   2. dt_0 = dt
    %   3. Compute predictor rn_P(t+dt) and corrector rn(t+dt) from forward
    %   Euler and trapezoid method respectively.
    %   4. If max(||rn_P(t+dt) -  rn(t+dt)||) > e, reduce timestep by half
    %   5. t = t + dt
    %   6. If dt = dt_0, increase the step to dt=min(1.2dt, dt_max)
    %   7. Return to 2, unless total number of cycles is reached

    while (~convergent)

        rnvec1 = rnvec0 + vnvec0 * dt; %Euler forward method [Cai & Bulatov, eq. 10.43]

        if isempty(rnvec1)
            convergent = 1;
        end

        % ET - rnvec1 can contain a node outside the domain!
        for iter = 1:maxiter,
            [vnvec, fn, fseg] = drndt(rnvec1, flag, MU, NU, a, Ec, links, connectivity, ...
                mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
            %err=rnvec1-rnvec0-vnvec.*dt;          %backward Euler
            err = rnvec1 - rnvec0 - (vnvec + vnvec0) / 2 * dt; %trapzoid
            errmag = max(abs(err));
            %disp(sprintf('iter=%d err=%e',iter,errmag));
            if (errmag < rntol)
                convergent = 1;
                break;
            else
                rnvec1 = rnvec1 - err;
            end

        end

        if (convergent)
            break;
        else
            dt = dt / 2;
        end

    end

    %unscramble rn and vn vectors
    rn = [reshape(rnvec1, length(rnvec1) / 3, 3), flag];
    vn = reshape(vnvec, length(vnvec) / 3, 3); % trapezoidal rule modification

    %When no time step reduction is necessary, the algorithm attempts to
    %increase dt by 20% for the next cycle, but not to exceed a preset value of
    %dtmax.

    %automatically adjust time step
    if isempty(errmag)
        dt = dt0;
    elseif ((dt == dt1) && (iter == 1))
        maxchange = 1.2;
        exponent = 20;
        factor = maxchange * (1 / (1 + (maxchange^exponent - 1) * (errmag / rntol)))^(1 / exponent);
        dt = min(dt1 * factor, dt0);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rn, vn, dt, fn, fseg] = int_trapezoid_stoc(rn, dt, dt0, MU, NU, a, Ec, links, connectivity, ...
        ~, rntol, mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d)
    % implicit numerical integrator using the Trapezoid method
    % dt: suggested timestep (usually from previous iteration)
    % dt0: maximum allowed timestep
    %
    % includes simplectic algorithm for stochastic noise in langevin simulation
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NB NOISE IS NOT CONSTRAINED BASED ON CRYSTALLOGRAPHY!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global LangevinDiffusionCoeff;

    %scramble rn into a single column vector
    rnvec = [rn(:, 1); rn(:, 2); rn(:, 3)]; flag = rn(:, 4);

    %Generate thermal noise (nodal based), from normal distribution with mean=0
    %scaled with sqrt(2*D), where D is the diffusion coefficient.

    % S.L.Dudarev et al./Journal of Nuclear Materials 455 (2014) 16-20
    % D(T) = (Dbar/2rho) exp(-Ea/kT), where Dbar = 5.7x10^13 nm^3/s, rho is the
    % loop radius (in nm units), Ea=1.3eV is the activation for loop migration,
    % and Eb approx. 0.4eV
    %what is D in DDLab units? Scales with loop radius apparently?
    vn_langevin = sqrt(2 * LangevinDiffusionCoeff) .* randn(size(rnvec, 1), size(rnvec, 2));

    %Backward Euler
    %rnvec0=rnvec;
    rnvec0 = rnvec + vn_langevin * 0.5 * dt; %Step 1 of simplectic algorithm

    %predict velocities from elastic interactions (predictor step)
    [vnvec0, fn, fseg] = drndt(rnvec0, flag, MU, NU, a, Ec, links, connectivity, ...
        mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);

    dt1 = dt;
    maxiter = 1;
    convergent = 0;

    % This algorithm [Cai & Bulatov, Algorithm 10.2, pg. 216] is a simple
    % routine to dynamically change the time-step. Whenever the difference
    % between the predicted and corrected positions for any node exceeds
    % threshold err, the timestep is decreased by half and new positions for
    % all nodes are recomputed with the reduced time step.
    %   1. Initialize time step dt = dt_max
    %   2. dt_0 = dt
    %   3. Compute predictor rn_P(t+dt) and corrector rn(t+dt) from forward
    %   Euler and trapezoid method respectively.
    %   4. If max(||rn_P(t+dt) -  rn(t+dt)||) > e, reduce timestep by half
    %   5. t = t + dt
    %   6. If dt = dt_0, increase the step to dt=min(1.2dt, dt_max)
    %   7. Return to 2, unless total number of cycles is reached

    while (~convergent)

        rnvec1 = rnvec0 + vnvec0 * dt; %Euler forward method [Cai & Bulatov, eq. 10.43]
        %This is also Step 2 of simplectic algorithm.

        for iter = 1:maxiter,
            [vnvec, fn, fseg] = drndt(rnvec1, flag, MU, NU, a, Ec, links, connectivity, ...
                mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
            %err=rnvec1-rnvec0-vnvec.*dt;          %backward Euler
            err = rnvec1 - rnvec0 - (vnvec + vnvec0) / 2 * dt; %trapzoid
            errmag = max(abs(err));
            %disp(sprintf('iter=%d err=%e',iter,errmag));
            if (errmag < rntol)
                convergent = 1;
                break;
            else
                rnvec1 = rnvec1 - err;
            end

        end

        if (convergent)
            %Step 3 of the simplectic algorithm, i.e. 2nd part of noise
            vnvec = vnvec + vn_langevin * 0.5 * dt;
            break;
        else
            dt = dt / 2;
            %Update Step 1 from simplectic algorithm...
            %This is important since dt needs to be the same for all three
            %steps of the simplectic algorithm. So, if we incorporate this in
            %the Euler forward method, we need to update rnvec0 to account for
            %the new dt being applied to the noise part.
            rnvec0 = rnvec + vn_langevin * 0.5 * dt;
            [vnvec0, fn, fseg] = drndt(rnvec0, flag, MU, NU, a, Ec, links, connectivity, ...
                mobility, vertices, uhat, nc, xnodes, D, mx, mz, w, h, d);
        end

    end

    %unscramble rn and vn vectors
    rn = [reshape(rnvec1, length(rnvec1) / 3, 3), flag];
    vn = reshape(vnvec, length(vnvec) / 3, 3); % trapezoidal rule modification

    %When no time step reduction is necessary, the algorithm attempts to
    %increase dt by 20% for the next cycle, but not to exceed a preset value of
    %dtmax.

    %automatically adjust time step
    if ((dt == dt1) && (iter == 1))
        maxchange = 1.2;
        exponent = 20;
        factor = maxchange * (1 / (1 + (maxchange^exponent - 1) * (errmag / rntol)))^(1 / exponent);
        dt = min(dt1 * factor, dt0);
    end

end
