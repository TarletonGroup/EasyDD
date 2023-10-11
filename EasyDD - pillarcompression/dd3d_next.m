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
%   disp('Compiling MEX files');
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" SegSegForcesMex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" StressDueToSegs.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" UtildaMex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" mindistcalcmex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG"  CollisionCheckerMex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" mobbcc1mex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" displacementmex_et.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" CreateInputMex.c %CollisionMarielle
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" CollisionCheckerMexMariellebis.c %CollisionMarielle
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" ./nodalForce/nodal_surface_force_linear_rectangle_mex.c
%  mexcuda -v COMPFLAGS="-Xptxas -v -arch=compute_60 -code=sm_60 -use_fast_math" ./nodalForce/nodal_surface_force_linear_rectangle_mex_cuda.cu
%   disp('Done!');
%  default value if run by itself (e.g. not through "rundd3d")
%  cleanup the empty node and link entries at the end of the initial data structures
clc;
clear all;
global tem amag RR QQ do_inclusion mumag SimBox
%tem: temperature
%amag: magnitude of lattice parameter a
%RR: list of the radius of particles
%QQ: list of centers of the particles
%do_inclusion: flag for particles
%mumag: magnitude of modulus MU
%SimBox: size of the simulation box

restart = 0; % flag for restart;

global USE_GPU;
USE_GPU=1; %0 if CPU only.

if (USE_GPU==1)
    disp('Going to use GPU as well...');  setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\amd64']);
    system('nvcc -ptx -m64 -arch sm_35 SegForceNBodyCUDADoublePrecision.cu');
end

% inital setup
if restart ==1
    load restart.mat % load the restart.mat file
    
    %construct stiffeness matrix K and pre-compute L,U decompositions.
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,mx,my,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D(dx,dy,dz,mz,MU,NU);
    disp('Done! Initializing simulation.');
else
    
    initial  % inital configuration and parameters
    
    [rn,links]=cleanupnodes(rn,links);
    
    % genererate the connectivity list from the list of links
    disp('Initiliazing connectivity list. Please wait.');
    [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
    
    consistencycheck(rn,links,connectivity,linksinconnect);
    disp('Consistencycheck : Done!');
        
    %construct stiffeness matrix K and pre-compute L,U decompositions.
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,mx,my,mel,unfixedDofs,Kred,Lred,Ured] = finiteElement3D(dx,dy,dz,mz,MU,NU);
%     disp('Done! Initializing simulation.');
    
    %%%% Addition by Daniel Celis Garza 12/19/2017.
    % Precomputing variables needed to couple tractions induced by dislocations
    % to FEM. Generate the list of planes to be extracted. The FEM coupler is
    % hardcoded so that the min(x) yz-plane does not have traction boundary
    % conditions.
    if(~exist('a_trac','var'))
        a_trac=0;
    end
    if (~exist('use_gpu', 'var'))
        use_gpu = 0;
    end
    if a_trac ~= 0
        f_hat = zeros(3*mno, 1);
        para_tol = dx/1e7;
        
        planes = (1:1:6)';
        [x3x6_lbl, x3x6, n_se] = extract_surface_nodes(xnodes, nc, [mx;my;mz],...
            planes, 4);
        gamma_dln = [gammat(:,1); gammaMixed(:,1)];
        [f_dln_node, f_dln_se,...
            f_dln, idxi, n_nodes_t] = nodal_force_map(x3x6_lbl, gamma_dln, 4, n_se, mno);
        tolerance = dx/10^7;
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
    end
    
    disp('Done! Initializing simulation.');
    
    
    %Use Delaunay triangulation to create surface mesh, used for visualisation
    
    %and for dislocation remeshing algorithm.
    [TriangleCentroids,TriangleNormals,tri,Xb] = ...
        MeshSurfaceTriangulation(xnodes,Stop,Sbot,Sfront,Sback,Sleft,Sright,gammaMixed);
    %Remesh considering surfaces in case input file incorrect.
    %disp('Remeshing...');
    [rn,links,connectivity,linksinconnect]=remesh_surf(rn,links,connectivity,linksinconnect,vertices,TriangleCentroids,TriangleNormals);
    
    % data=zeros(totalsteps,1);
    if(~exist('dt','var'))
        dt=dt0;
    end
    dt=min(dt,dt0);
    %plotCounter=1;
    close all
    
    
    gamma=[gammau;gammat];
    utilda_0=zeros(3*mno,1);
    gn = gamma(:,1);
    [Ux, Uy, Uz] = Utilda_bb3(rn,links,gn,NU,xnodes,dx,dy,dz,mx,my,mz);
    
    utilda_0(3*gn -2) = Ux;
    utilda_0(3*gn -1) = Uy;
    utilda_0(3*gn   ) = Uz;
    
    % parameters for the force loading condition
    F_c = 400/mumag*(dx*dy); % critial force applied on top
    simTime0=0;
    plasticstrain=0;
    
    %number of glide and climb steps
    step_glide=0;
    step_climb=0;
    %strain of glide ad climb
    strain_glide=0;
    strain_climb=0;
end


%% start simulation
while simTime < totalSimTime*1e6
    
    curstep=curstep+1;
    
    disp(curstep);
    
    if (length(links(:,1))<1)
        fprintf('there is no link\n');
        break;
    end
    
    %DDD+FEM coupling
    if a_trac == 0
        [uhat,Ftop,Uend] = FEMcoupler(rn,links,maxconnections,a,MU,NU,xnodes,mno,kg,L,U,...
            Stop,gammau,gammat,gammaMixed,fixedDofs,freeDofs,dx,dy,dz,simTime0,mx,my,mz,unfixedDofs,Kred,Lred,Ured,utilda_0);
    else
        [uhat,fend,Ubar] = analytic_FEMcoupler(rn,links,a,MU,NU,xnodes,mno,kg,L,U,...
            gamma_disp, gammat, gammaMixed, fixedDofs,freeDofs,dx,dy,dz,simTime,mx,my,mz,utilda_0,...
            gamma_dln, x3x6, 4, n_nodes_t, n_se, idxi, ...
            f_dln_node, f_dln_se, f_dln, f_hat, use_gpu, n_threads, para_scheme, para_tol);
    end
    
    fprintf('fend = %d, Ubar = %d, simTime = %d\n', Ftop/dx/dy*mumag, Uend, simTime);
    
    %%%%%%%%%%%%%%%%%%%%%checke the stress state
    if abs(Ftop) > F_c
        Ftop_c=Ftop;
        if do_climb==1
            strat_change=1;  %start to check the change of type of motion,
        end
    else
        simTime0=curstep;
    end
    
    %-----dislocation density
    Volume = SimBox(1)*SimBox(2)*SimBox(3);
    volume = Volume*(amag)^2;
    
   
    %remove the virtual nodes
    virtual = [find(rn(links(:,1),end)==67);find(rn(links(:,2),end)==67);];
    virtual = unique(virtual,'stable');
    
    dx1=rn(links(:,2),1:3)-rn(links(:,1),1:3);
    dx1(virtual,:) = 0;
    disLength=sum(vecnorm(dx1',2));
    
    disDensity = disLength/volume;  %1/um^2 or 10^12/m^2
    
    
    %integrating equation of motion
    %%%% determine the mobility law, glide/climb
    if flagc
        mobility = 'mobbcc_linear_full_climb';        
        %%integrating equation of motion
        %%calculate the the nodal force, nodal velocity, update rn, links
        [rnnew,vn,dt,fseg,links]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
            rmax,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);  
        
        real_dt=dt; %s        
    else
        mobility = 'mobbcc_bb1b_origin';        
        %%% calculate the the nodal force, nodal velocity, update rn, links
        [rnnew,vn,dt,fseg,links]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
            rmax,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);    
        
        real_dt=dt*1e-2/mumag/1e6; %s        
    end
       
    %--------plot the configuration and save figures, please refer to plot_gif to make videos or .gif
    if mod(curstep,plotfreq)==1
        
        if do_visualisestress==1
            visualisedislocationstress(SimBox,curstep,rnnew,links,a,MU,NU);
        end
        
        PlotBox_b(SimBox);
        plotnodes_climb(connectivity,rnnew,links,SimBox,simTime,dt)
        
        strtime = ['Time = ', num2str(simTime),'s' ];
        text(-1*plim,0.5*plim, 2.0*plim,strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'Color','b');
        if do_inclusion
            plot_sphere(RR,QQ,amag)
        end
       
        viewangle = [50 10];
        view(viewangle);
        
        print(1,'-dbmp',sprintf('%d',num_fig))
        num_fig=num_fig+1;
        
    end    
    
    rnnew=[rnnew(:,1:3) vn rnnew(:,4)];
    linksnew=links;
    connectivitynew=connectivity;
    linksinconnectnew=linksinconnect;
    fsegnew=fseg;
    
    % plastic strain and plastic spin calculations
    [ep_inc,~]=calcplasticstrainincrement_HY(rnnew,rn,links,Volume);
    plasticstrain = plasticstrain + ep_inc(3,3);
    
    % calculate the amount of glide and climb;
    if flagc
        step_climb = step_climb + 1;
        strain_climb = strain_climb+ ep_inc(3,3);
    else
        step_glide = step_glide + 1;
        strain_glide = strain_glide + ep_inc(3,3);
    end
    ratio_step = step_climb/step_glide;
    ratio_strain = strain_climb/strain_glide;
    
    %%------output stress-strain curve and dislcoation density
    string1=strcat('.\output\','output');
    if (curstep==1)
        fid = fopen(strcat(string1,'_result.dat'),'w');
    else
        fid = fopen(strcat(string1,'_result.dat'),'a');
    end
    fprintf(fid,'%e    %e     %e    %e    %e\n', simTime, plasticstrain, Uend,  Ftop, disDensity);
    fclose(fid);
    
    %% topologic evolution
    if flagc  % topological evolution during climb steps
        if (doremesh) %do virtual re-meshing first
            %remeshing virtual dislocation structures
            if (dovirtmesh)
                [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen3(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,MU,NU,a,Ec);
               
            end
            %remeshing internal dislocation structures
            [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all_c(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,dovirtmesh,vertices,...
                uhat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
        end
               
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
                
                if colliding_segments == 1 %scan and update dislocation structure.
                   
                    fprintf(' Segment %d and segment %d are colliding\n',s1,s2);
                    if colliding_segments==1
                        [rnnew,linksnew,~,~,fsegnew,colliding_segments]=collision_basic_c(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,vertices,...
                            uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair,lmin);
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
        
        if (doseparation)
            if max(connectivitynew(:,1))>3
                %spliting of nodes with 4 or more connections
                [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=...
                    separation_c(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                    fsegnew,mobility,MU,NU,a,Ec,2*rann,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
            end
        end
        
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all_c(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,0,vertices,...
            uhat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
        if any(any(isnan(rnnew)))
            disp('The new extended node is NaN! See Line 46 in remesh_surf.m');
            pause;
        end
        
        
    else
        
        if (doremesh) %do virtual re-meshing first
            %remeshing virtual dislocation structures
            if (dovirtmesh)
                [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen3(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,MU,NU,a,Ec);
            end
            %remeshing internal dislocation structures
            [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,dovirtmesh,vertices,...
                uhat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
        end
        
        if any(any(isnan(rnnew)))
            disp('The new extended node is NaN! See Line 46 in remesh_surf.m');
            pause;
        end
        
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

                if colliding_segments == 1 %scan and update dislocation structure.
                    
                    fprintf(' Segment %d and segment %d are colliding\n',s1,s2);
                    if colliding_segments==1
                        [rnnew,linksnew,~,~,fsegnew,colliding_segments]=collision_basic(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,vertices,...
                            uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair,lmin);
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
        
        if (doseparation)
            if max(connectivitynew(:,1))>3
                %spliting of nodes with 4 or more connections
                [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=...
                    separation(rnnew,linksnew,connectivitynew,linksinconnectnew,...
                    fsegnew,mobility,MU,NU,a,Ec,2*rann,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
            end
        end
        
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh_all(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,mobility,doremesh,0,vertices,...
            uhat,nc,xnodes,D,mx,mz,w,h,d,TriangleCentroids,TriangleNormals);
        if any(any(isnan(rnnew)))
            disp('The new extended node is NaN! See Line 46 in remesh_surf.m');
            pause;
        end
    end
    
    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    fseg=fsegnew;
    
    %%  reset the type of mobility
    if strat_change==1
        if flagc == 0
            dotarea = norm(ep_inc(3,3))/real_dt;  % strain rate            
            if(dotarea< crit_glide)                
                num_glide=num_glide+1;
                if num_glide > numofstep   % since v_g*t_c is very large, 20 steps are set to gurantee that the time step has recovered to t_g.
                    disp('climb is activated');
                    flagc=1;
                    num_climb=0;
                end                
            end
            
        elseif  (max(vecnorm(vn,2,2))*real_dt > crit_climb)
            
            num_climb=num_climb+1;
            if num_climb > numofstep  % 20 steps are set to gurantee that the time step has recovered to t_g.
                flagc=0;
                disp('glide is activated');
                num_glide=0;
            end
        end
        
    end
    
    %store run time information
    simTime = simTime+real_dt;
    if mod(curstep, 5*savefreq)==1
        fname=num2str(curstep);
        save(fname,'-regexp','^(?!(K|kg|Kred|Lred|Ured)$).')
    end
    

end

save restart;
disp('completed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%