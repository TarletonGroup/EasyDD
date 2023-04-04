%Features:
% free boundary condition (no pbc)
% linear mobility law for glide 
% N^2 interaction (no neighbor list, no fast multipole)

%Data structure:
% NMAX:    maximum number of nodes (including disabled ones)
% LINKMAX: maximum number of links (including disabled ones)
% rn: (NMAX,4) array of nodal positions (last column is flag: -1 means disabled)
% vn: (NMAX,3) array of nodal velocities
% fn: (NMAX,3) array of nodal forces
% links: (LINKMAX,8) array of links (id1,id2,bx,by,bz,nx,ny,nz)
function dd3d_couple( )


clc; clear all;
global maxconnections 
%flag_b; % flag_b=1 if lattice diffusion is considered, otherwise flag_b=0;
restart=0; 


if (restart)
    %read in the restart information
    load('restart.mat');
    
else
    % generate a new initial dislcoation structure  
    initial      
    kk=0;
    totstrain = 0;
    flagc=0;   %  0 for glide; 1 for climb
    tottime = 0;
        
end


[rn,links]=cleanupnodes(rn,links);
[connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);

 consistencycheck(rn,links,connectivity,linksinconnect,SimBox);
    disp('Consistencycheck : Done!');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag_b
        %construct stiffeness matrix K and pre-compute L,U decompositions.
        disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.');
        [B,xnodes,mno,nc,n,ke,gammaj,gammaU,wx,wy,wz,mel,my,mz] = finiteElement3D_new(dx,dy,dz,mx,Dv);
        disp('Done! Initializing simulation.');
    else
        xnodes=[];
        mno=[];
        nc=[];
        ke=[];
        wx=[];
        wy=[];
        wz=[];
        mel=[];
        mx=[];
        my=[];
        mz=[];
        Dv=[];
        ke=[];
    end

data=zeros(totalsteps,1);
if(~exist('dt'))
    dt=dt0;
end
dt=min(dt,dt0);


global USE_GPU;
USE_GPU=0; %0 if CPU only.

if (USE_GPU==1)
    disp('Going to use GPU as well...'); %
    setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\amd64']);
    system('nvcc -ptx -m64 -arch sm_35 SegForceNBodyCUDADoublePrecision.cu');
end  

    
%%
telapse=0;
mmmm=1;

simtime=0; 
rn0=rn;
links0=links;
 for curstep=kk+1:totalsteps
        
    if mod(curstep,1)==0
        disp(curstep);
       
    end
   
   
    if (length(links(:,1))<1)
        fprintf('there is no link\n');
        break;
    end

    %clean up nodes
    llinks = size(links,1);
    if (links(llinks,1)==0)
        [rn,links]=cleanupnodes(rn,links);
        % genererate the connectivity list from the list of links
        [connectivity,linksinconnect]=genconnectivity(rn,links,maxconnections);
    end
         
    %plot & save config
   if mod(curstep,200)==1
               
        plim =SimBox(1);
        plotnodes_initial(rn0,links0,plim);
        plotnodes(rn,links,plim);

        PlotBox_b(SimBox);
        xlim([0 plim]);
        ylim([0 plim]);
        zlim([0 plim]);
        view([0 0 1]);
        
        time=simtime;
        strtime = ['Time = ', num2str(time),'s' ];
        text(0.2*plim,1.0*plim,0*plim,strtime,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',16,'Color','b');
        
%         print(1,'-dbmp',sprintf('%d',num_fig))
%         num_fig=num_fig+1;
       
    end
    
%calculate nodal velocity
    [rnnew,vn,dt,fseg]=feval(integrator,rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
    appliedstress,rmax,rntol,mobility,SimBox,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv);
    
     simtime = simtime +dt;

    
    rnnew = [rnnew(:,1:3) vn rnnew(:,4)];
    linksnew = links;
    
    connectivitynew = connectivity;
    linksinconnectnew = linksinconnect;
    fsegnew = fseg;
 
    
    if(doseparation)
        %spliting of nodes with 4 or more connections
         [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mmmm]=separation(mmmm,rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mobility,MU,NU,a,Ec,2*rann,appliedstress,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv);
    end

    if(docollision)
        %collision detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,mmmm]=collision(mmmm,rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv);
    end
    
    if(Boundary(1))
        %surface interaction detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=touchsur(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,SimBox(1),1);
    end
    if(Boundary(2))
        %surface interaction detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=touchsur(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,SimBox(2),2);
    end
    if(Boundary(3))
        %surface interaction detection and handling
        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=touchsur(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,rann,MU,NU,a,Ec,mobility,appliedstress,SimBox(3),3);
    end
    
    if(doremesh)
        %remesh
       [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv);
    end
%     if(curstep==1)
%         %remesh
%        [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=remesh(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,lmin,lmax,areamin,areamax,MU,NU,a,Ec,appliedstress,mobility);
%     end
    
%     
 

    rn=[rnnew(:,1:3) rnnew(:,7)];
    vn=rnnew(:,4:6);
    links=linksnew;
    connectivity=connectivitynew;
    linksinconnect=linksinconnectnew;
    %store run time information
    %time step
    data(curstep,1)=dt;
    %save restart
    if mod(curstep, 1*savefreq)==1
        fname=num2str(curstep);
        save(fname,'-regexp','^(?!(K|kg|Kred|Lred|Ured)$).')
    end

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%