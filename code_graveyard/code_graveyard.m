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
[rn,links]=cleanupnodes(rn,links);

% genererate the connectivity list from the list of links
disp('Initiliazing connectivity list. Please wait.');
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);

consistencycheck(rn,links,connectivity,linksinconnect);
disp('Consistencycheck : Done!');

%construct stiffeness matrix K and pre-compute L,U decompositions.
disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
[vertices,B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading); % Haiyang's addition

% [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
%     Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
%     w,h,d,my,mz,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D_haiyang2(dx,dy,dz,mx,MU,NU,loading);

%% Addition by Daniel Celis Garza 12/19/2017.
% Precomputing variables needed to couple tractions induced by dislocations
% to FEM. Generate the list of planes to be extracted. The FEM coupler is
% hardcoded so that the min(x) yz-plane does not have traction boundary
% conditions.
f_hat = zeros(3*mno, 1);

planes = (1:1:6)';
[x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
    planes, 4);
gamma_dln = [gammat(:,1); gammaMixed(:,1)];
[f_dln_node, f_dln_se,...
    f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
para_tol = dx/1e7;
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
USE_GPU=0; %0 if CPU only.

if (USE_GPU==1)
    disp('Going to use GPU as well...'); %setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\amd64']);
    system('nvcc -ptx -m64 -arch sm_35 SegForceNBodyCUDADoublePrecision.cu');
end


%Use Delaunay triangulation to create surface mesh, used for visualisation

%and for dislocation remeshing algorithm.
[TriangleCentroids,TriangleNormals,tri,Xb] = ...
    MeshSurfaceTriangulation(xnodes,Stop,Sbot,Sfront,Sback,Sleft,Sright,gammaMixed);
%Remesh considering surfaces in case input file incorrect.
%disp('Remeshing...');
[rn,links,connectivity,linksinconnect]=remesh_surf(rn,links,connectivity,linksinconnect,vertices,TriangleCentroids,TriangleNormals);

%plot dislocation structure
% figure(1);
% plotHandle = plotnodes(rn,links,plim,vertices); view(viewangle);
% drawnow

% data=zeros(totalsteps,1);
if(~exist('dt','var'))
    dt=dt0;
end
dt=min(dt,dt0);
plotCounter=1;
close all

Fend=zeros(1e6,1); fend=[];
U_bar=zeros(1e6,1); Ubar=[];
t=zeros(1e6,1); simTime=0;
%%
gamma=[gammau;gammaMixed];
utilda_0=zeros(3*mno,1);
gn = gamma(:,1);
[Ux, Uy, Uz] = Utilda_bb3(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);

utilda_0(3*gn -2) = Ux;
utilda_0(3*gn -1) = Uy;
utilda_0(3*gn   ) = Uz;
%%
while simTime < totalSimTime
    
    % frame recording
    intSimTime=intSimTime+dt;
    if intSimTime > dtplot && doplot == 1
        %plotHandle=plotnodes(rn,links,plim,vertices);view(viewangle);
        plotCounter=plotCounter+1;
        plotCounterString=num2str(plotCounter,'%03d');
        %saveas(plotHandle,plotCounterString,'png')
        %save(plotCounterString,'rn','links','fend','Ubar','simTime');
        %plotHandle=plotnodes(rn,links,plim,vertices);view(viewangle);
        %plotnodes(rn,links,plim,vertices);view(viewangle);
        %zlim([-100 100])
        %xlim([-100 100])
        %ylim([-100 100])
        intSimTime=intSimTime-dtplot;
    end
    
    %DDD+FEM coupling
    if a_trac == 0
        [uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
            gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,dy,dz,simTime,mx,my,mz,utilda_0);
    else
        [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
            gamma_disp, gammat, gammaMixed, fixedDofs,freeDofs,dx,dy,dz,simTime,mx,my,mz,utilda_0,...
            gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, ...
            f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, para_tol);
    end
    Fend(curstep+1) = fend;
    U_bar(curstep+1) = Ubar;
    t(curstep+1) = simTime;
    fprintf('fend = %d, Ubar = %d, simTime = %d \n',fend,Ubar,simTime);
    
    %     if (dovirtmesh)
    %         %remeshing virtual dislocation structures
    %         %[rn,links,connectivity,linksinconnect]=remesh_surf(rn,links,connectivity,linksinconnect,vertices,TriangleCentroids,TriangleNormals);
    %         %[rn,links,connectivity,linksinconnect] = virtualmeshcoarsen2(rn,links,maxconnections,10*lmin);
    %     end
    
    
    %integrating equation of motion
    [rnnew,vn,dt,fn,fseg]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
        rmax,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
    
    % Haiyang's screw dislocation correction.
    [rn,links] = cross_slip_BCCrotate5(fseg,rn,links,connectivity,0,...
        0,curstep,0,NU,a,Ec,0,...
        uhat,nc,xnodes,D,mx,mz,w,h,d);
    
    % plastic strain and plastic spin calculations
    [ep_inc,wp_inc]=calcplasticstrainincrement(rnnew,rn,links,(2*plim)^3);
    
    if(mod(curstep,printfreq)==0) && not(isempty(vn))
        fprintf('step%3d dt=%e v%d=(%e,%e,%e) \n',...
            curstep,dt,printnode,vn(printnode,1),vn(printnode,2),vn(printnode,3));
        close all;
        save restart_temp
    elseif (mod(curstep,printfreq)==0)
        fprintf('step%3d dt=%e v%d=(%e,%e,%e) \n',...
            curstep,dt,printnode);
    end
    
    if(mod(curstep,plotfreq)==0)
        plotnodes(rn,links,plim,vertices);
        view(viewangle);
        drawnow
        pause(0.01);
        %
        %         figure(2);
        %          clf
        %          plot(U_bar*bmag,-Fend*bmag^2*mumag)
        %          pause(0.01);
    end
    
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    
    %save restart.mat
    if (docollision)
        s1_old=0;
        s2_old=0;
        %collision detection and handling
        
        %COLLISION GPU / GPUPLUS  Marielle  + CollisionCheckerMexMarielle
        %it requires a GPU device
        %it requires CollisionCheckerMexMarielle, collisionGPUplus (or collision GPU), mindistcalcGPU1, mindistcalcGPU2,CreateInputMex
        colliding_segments=1;
        while colliding_segments==1
            [colliding_segments,n1s1,n2s1,n1s2,n2s2,floop,s1,s2,segpair]=CollisionCheckerMexMariellebis(rnnew(:,1),rnnew(:,2),rnnew(:,3),rnnew(:,end),...
                rnnew(:,4),rnnew(:,5),rnnew(:,6),linksnew(:,1),linksnew(:,2),connectivitynew,rann);
            if floop==2
                colliding_segments=0;
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
            if colliding_segments == 1 %scan and update dislocation structure.
                
                %                 [rnnew,linksnew,~,~,fsegnew]=...
                %                 collision(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                %                 fsegnew,rann,MU,NU,a,Ec,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                
                if floop==2
                    %                     if s1==s1_old && s2==s2_old
                    colliding_segments=0;
                    %                     else
                    %                         s1_old=s1;
                    %                         s2_old=s2;
                    %                     end
                end
                
                if mod(curstep,100)~=0 && floop==1
                    nodeids1=connectivitynew(n1s1,1);
                    nodeids1=1:nodeids1;
                    nodeids1=2*nodeids1;
                    nodeids2=connectivitynew(n2s1,1);
                    nodeids2=1:nodeids2;
                    nodeids2=2*nodeids2;
                    nodeids3=connectivitynew(n1s2,1);
                    nodeids3=1:nodeids3;
                    nodeids3=2*nodeids3;
                    nodeids4=connectivitynew(n2s2,1);
                    nodeids4=1:nodeids4;
                    nodeids4=2*nodeids4;
                    alt1=connectivitynew(n1s1,nodeids1);
                    alt2=connectivitynew(n2s1,nodeids2);
                    alt3=connectivitynew(n1s2,nodeids3);
                    alt4=connectivitynew(n2s2,nodeids4);
                    other1=any(alt1==alt3');
                    other1=sum(other1);
                    other2=any(alt1==alt4');
                    other2=sum(other2);
                    other3=any(alt2==alt3');
                    other3=sum(other3);
                    other4=any(alt2==alt4');
                    other4=sum(other4);
                    if other1>eps || other3>eps || other2>eps || other4>eps
                        fprintf('Remeshing conflict detected. Cancelling collision\n')
                        colliding_segments=0;
                    end
                end
                
                
                fprintf(' Segment %d and segment %d are colliding\n',s1,s2);
                if colliding_segments==1
                    [rnnew,linksnew,~,~,fsegnew,colliding_segments]=collision_basic(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,vertices,...
                        uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair);
                end
            end
            %removing links with effective zero Burgers vectors
            [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew] = cleanupsegments(rnnew,linksnew,fsegnew);
        end
    end
    
    %find blockading fixed nodes and allow code to deal with them
    for i=1:size(rnnew,1)
        if rnnew(i,end)~=7
            continue
        end
        if connectivitynew(i,1)>7
            rnnew(i,end)=0;
        end
    end
    
    if (doseparation) && max(connectivitynew(:,1))>3
        %spliting of nodes with 4 or more connections
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=...
            separation(rnnew,linksnew,connectivitynew,linksinconnectnew,...
            fsegnew,mobility,MU,NU,a,Ec,2*rann,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
    end
    
    if (doremesh) %do virtual re-meshing first
        %remeshing virtual dislocation structures
        if (dovirtmesh)
            %[rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen_mex(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,DIST_SOURCE*0.49,dx,MU,NU,a,Ec);
            [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,DIST_SOURCE*0.49,dx,MU,NU,a,Ec);
            %[rnnew,linksnew,connectivitynew,linksinconnectnew] = virtualmeshcoarsen2(rnnew,linksnew,maxconnections,10*lmin);
        end
        %remeshing internal dislocation structures
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,dovirtmesh,vertices,...
            uhat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
    end
    
    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    
    %store run time information
    %time step
    curstep = curstep + 1;
    %data(curstep,1)=dt;
    simTime = simTime+dt;
    
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
    if ~isfile("UtildaMex.c") || ~isfile(file.name) && days(file.date - datetime('now')) > 30
        mex -v COPTIMFLAGS="-o3 -oy -DNDEBUG" UtildaMex.c   
    end
end

%%
function visualise
%%



% % code below is just to visulise the defomed mesh

clear all
close all
disp('Loading restart file');
load restart.mat %restart.50sources.75744.mat

plotmesh=1;
plotstress=0;

%refine mesh
mx=mx*2;

%update results
disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
     Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
     w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);    
disp('Done!');

[TriangleCentroids,TriangleNormals,tri,Xb] = ...
     MeshSurfaceTriangulation(xnodes,Stop,Sbot,Sfront,Sback,Sleft,Sright);

[uhat,fend,Ubar] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
                    gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,simTime);

segments = constructsegmentlist(rn,links);
utilda = zeros(3*mno,1);
gn = 1:mno; % global node number
x0 = xnodes(gn,1:3); % field point
point_array_length = size(x0,1);
segments_array_length = size(segments,1);
%Full C version (including combinatorials, checks, corrections etc.)

disp('Calculating displacements');
[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);                       
%[Uxf, Uyf, Uzf] =displacement_fivel(x0,segments,NU); %gives same answer

utilda(3*gn -2) = Ux;
utilda(3*gn -1) = Uy;
utilda(3*gn   ) = Uz;

if plotmesh
    disp('Plotting mesh');
    xp = zeros(mno,3);
    for j =1:mno;
         xp(j,1:3) = xnodes(j,1:3)...
             + 1e3*utilda(3*j-2:3*j)'; %+ 0e4*uhat(3*j-2:3*j)';
    end
%     amag=1;
%     xp = amag*xp;
    figure;clf;hold on;view(0,0)
    xlabel('x');ylabel('y');zlabel('z')
    style='-k';
    for p =1:mel
    %      plot3(xp(:,1),xp(:,2),xp(:,3),'.') % plot nodes
        % plot elements
        plot3(xp(nc(p,[1:4,1]),1),xp(nc(p,[1:4,1]),2),xp(nc(p,[1:4,1]),3),style) 
        plot3(xp(nc(p,[5:8,5]),1),xp(nc(p,[5:8,5]),2),xp(nc(p,[5:8,5]),3),style) 
        plot3(xp(nc(p,[1,5]),1),xp(nc(p,[1,5]),2),xp(nc(p,[1,5]),3),style) % 
        plot3(xp(nc(p,[2,6]),1),xp(nc(p,[2,6]),2),xp(nc(p,[2,6]),3),style) % 
        plot3(xp(nc(p,[3,7]),1),xp(nc(p,[3,7]),2),xp(nc(p,[3,7]),3),style) % 
        plot3(xp(nc(p,[4,8]),1),xp(nc(p,[4,8]),2),xp(nc(p,[4,8]),3),style) % 
    end
    axis equal
    zlabel('z-direction (\mum)');
    xlabel('x-direction (\mum)');
%     zlim([-6 6])
%     xlim([0 9e4])
    title('$\tilde{u}$ scaled','FontSize',14,'Interpreter','Latex');
end

%-------------------generate stress--------------------------
if plotstress
    X = linspace(0,dx,2*mx)';
    Z = linspace(0,dy,2*my)';
    Y = 0.5*dy; %middle slide

    X_size = length(X);
    Y_size = length(Y);
    Z_size = length(Z);

    Sxxu = zeros(X_size,Y_size);
    Syyu = zeros(X_size,Y_size);
    Szzu = zeros(X_size,Y_size);
    Sxyu = zeros(X_size,Y_size);
    Sxzu = zeros(X_size,Y_size);
    Syzu = zeros(X_size,Y_size);
    Sxx = zeros(X_size,Y_size);
    Syy = zeros(X_size,Y_size);
    Szz = zeros(X_size,Y_size);
    Sxy = zeros(X_size,Y_size);
    Sxz = zeros(X_size,Y_size);
    Syz = zeros(X_size,Y_size);
    p1x = segments(:,6);
    p1y = segments(:,7);
    p1z = segments(:,8);
    p2x = segments(:,9);
    p2y = segments(:,10);
    p2z = segments(:,11);
    bx = segments(:,3);
    by = segments(:,4);
    bz = segments(:,5);
    x1 =  [ p1x,p1y,p1z];
    x2 = [p2x,p2y,p2z];
    b=[bx, by, bz];
    
    disp('Calculating stresses');
    for i= 1:X_size;
        for j = 1:Z_size
            x0 = [X(i) Y Z(j)]; % field point
            sigmahat = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
            Sxxu(i,j) = sigmahat(1,1);
            Syyu(i,j) = sigmahat(2,2);
            Szzu(i,j) = sigmahat(3,3);

            Sxyu(i,j) = sigmahat(1,2); %isotropic
            Sxzu(i,j) = sigmahat(1,3); %isotropic
            Syzu(i,j) = sigmahat(2,3); %isotropic

            sigmatilde=FieldPointStress(x0,x1,x2,b,a,MU,NU);
            Sxx(i,j) = sigmatilde(1);
            Syy(i,j) = sigmatilde(2);
            Szz(i,j) = sigmatilde(3);
            Sxy(i,j) = sigmatilde(4);
            Syz(i,j) = sigmatilde(5);
            Sxz(i,j) = sigmatilde(6);

        end
    end
    
    figure; clf
    subplot(3,1,1)
    surf(X*amag,Z*amag,mumag*(Sxxu+Sxx)','EdgeColor','none'); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)')
    ylabel('z-direction (\mum)')
    title('$$\hat{\sigma}_{xx}$$+$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    h=colorbar;
    xlabel(h,'MPa');
    
    subplot(3,1,2)
    surf(X*amag,Z*amag,mumag*Sxx','EdgeColor','none'); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)')
    ylabel('z-direction (\mum)')
    title('$$\tilde{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    h=colorbar;
    xlabel(h,'MPa');

    subplot(3,1,3)
    surf(X*amag,Z*amag,mumag*Sxxu','EdgeColor','none'); 
    view(2)
    axis equal;
    axis([0 dx*amag 0 dz*amag])
    xlabel('x-direction (\mum)')
    ylabel('z-direction (\mum)')
    title('$$\hat{\sigma}_{xx}$$','Interpreter','Latex');
    grid off
    %saveas(gcf,'sxx','epsc')
end

end

%%
function tractionCheck
%load a test condition
run inputCheck.m

segments = constructsegmentlist(rn,links);

%construct finite element arrays
%construct stiffeness matrix K and pre-compute L,U decompositions.
disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
[B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
    Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
    w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);    

X = linspace(0,dx,2*mx)';
Z = linspace(0,dy,2*my)';
Y = 0.5*dy; 

X_size = length(X);
Y_size = length(Y);
Z_size = length(Z);

Sxxuhat = zeros(X_size,Z_size);
Syyuhat = zeros(X_size,Z_size);
Szzuhat = zeros(X_size,Z_size);
Sxyuhat = zeros(X_size,Z_size);
Sxzuhat = zeros(X_size,Z_size);
Syzuhat = zeros(X_size,Z_size);

Sxxutilde = zeros(X_size,Z_size);
Syyutilde = zeros(X_size,Z_size);
Szzutilde = zeros(X_size,Z_size);
Sxyutilde = zeros(X_size,Z_size);
Sxzutilde = zeros(X_size,Z_size);
Syzutilde = zeros(X_size,Z_size);

Sxx = zeros(X_size,Z_size);
Syy = zeros(X_size,Z_size);
Szz = zeros(X_size,Z_size);
Sxy = zeros(X_size,Z_size);
Sxz = zeros(X_size,Z_size);
Syz = zeros(X_size,Z_size);

%%
%---------------------------------------------------------------
% evaluate tractions on gammat and solve for uhat  to recover
% dislocation stress field 

gn = 1:mno; %gammau(:,1);  % global node number
x0 = xnodes(gn,1:3); % field point
point_array_length = size(x0,1);
segments_array_length = size(segments,1);

[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);

utilda(3*gn - 2) = Ux;
utilda(3*gn - 1) = Uy;
utilda(3*gn   ) = Uz;

%% ---------------------------------------------------------------
for i= 1:X_size;
    
    for j = 1:Z_size
        
        x0 = [X(i) Y Z(j)]; % field point
        
        %sigmahatFEM = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,x0);
        %Sxxuhat(i,j) = sigmahatFEM(1,1);
        %Syyuhat(i,j) = sigmahatFEM(2,2);
        %Szzuhat(i,j) = sigmahatFEM(3,3);
        %Sxyuhat(i,j) = sigmahatFEM(1,2); %isotropic
        %Sxzuhat(i,j) = sigmahatFEM(1,3); %isotropic
        %Syzuhat(i,j) = sigmahatFEM(2,3); %isotropic
        
        sigmatildeFEM = hatStress(utilda,nc,xnodes,D,mx,mz,w,h,d,x0);
        Sxxutilde(i,j) = sigmatildeFEM(1,1);
        Syyutilde(i,j) = sigmatildeFEM(2,2);
        Szzutilde(i,j) = sigmatildeFEM(3,3);
        Sxyutilde(i,j) = sigmatildeFEM(1,2); %isotropic
        Sxzutilde(i,j) = sigmatildeFEM(1,3); %isotropic
        Syzutilde(i,j) = sigmatildeFEM(2,3); %isotropic
        
        x1=segments(:,6:8);
        x2=segments(:,9:11);
        b=segments(:,3:5);
        sigmatilde=FieldPointStress(x0,x1,x2,b,a,MU,NU);
        Sxx(i,j) = sigmatilde(1);
        Syy(i,j) = sigmatilde(2);
        Szz(i,j) = sigmatilde(3);
        Sxy(i,j) = sigmatilde(4);
        Syz(i,j) = sigmatilde(5);
        Sxz(i,j) = sigmatilde(6);           
    end
end

%%
subplot(5,2,1)
contourf(X,Z,Sxx'); 
axis([0 dx/2 0 dz])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xx}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);
subplot(5,2,2)
contourf(X,Z,Sxxutilde'); 
axis([0 dx/2 0 dz])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xx}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,3)
contourf(X,Z,Syy'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{yy}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,4)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Syyutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{yy}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,5)
contourf(X,Z,Szz'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{zz}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,6)
%c=1E-4*linspace(-1,1,100);
contourf(X,Z,Szzutilde'); 
%surf(X,Z,Szzutilde','EdgeColor','none'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{zz}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,7)
contourf(X,Z,Sxy'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xy}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,8)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Sxyutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xy}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,9)
contourf(X,Z,Sxz'); 
axis([0 dx/2 0 dy])
axis equal
xlabel('x-direction','FontSize',14)
%ylabel('z-direction','FontSize',14)
title(['\sigma_{xz}^{Analytical}'],'FontSize',14)
caxis([-5e-4 5e-4]);

subplot(5,2,10)
%c=1E-4*linspace(-1,1,100);
% contourf(X,Z,Sxx',c,'EdgeColor','none'); 
contourf(X,Z,Sxzutilde'); 
axis([0 dx/2 0 dy])
axis equal;
xlabel('x-direction','FontSize',14)
ylabel('z-direction','FontSize',14)
title(['\sigma_{xz}^{FEM}'],'FontSize',14)
caxis([-5e-4 5e-4]);
end

%%
function displacement_test


segments = constructsegmentlist(rn,links);

x0 = [linspace(0,dx)', 2.5e4*ones(100,1), 4.5e4*ones(100,1)];

tic;
[Ux_f,Uy_f,Uz_f]  = displacement_fivel(x0,segments,NU);
toc;

point_array_length = size(x0,1);
segments_array_length = size(segments,1);
tic;
[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);
utilda2 = horzcat(Ux,Uy,Uz);
toc;

norm([Ux-Ux_f,Uy-Uy_f,Uz-Uz_f]);

clear utilda3
for i =1:size(x0,1)
    
uspf = displacement_spf(x0(i,:),rn(:,1:3),segments(1,3:5),NU);
utilda3(i,1:3) = uspf;
end

end