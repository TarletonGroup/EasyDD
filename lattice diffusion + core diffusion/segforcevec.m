
function [fseg]=segforcevec(MU,NU,a,Ec,rn,links,sigext,linkid)
%compute nodal driving force of dislocation network by stress formula
%(vectorized version)
%rn: nodal position
%links: dislocation segments (connectivity, Burgers vector)
%sigext: external stress (applied)

%segment contribution to node force

%NMAX: max number of nodes[NMAX,m]=size(rn);
[NMAX,m]=size(rn);
if(m~=4)
    disp('rn should have 4 columns!');
    return;
end
[LINKMAX,m]=size(links);

%construct segment list
segments=constructsegmentlist(rn,links);
if linkid==0
    %PK force due to applied stress
    fpk=pkforcevec(sigext,segments);

    %self force due to self stress
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments);

    %remote force due to remote stress
    [fr0,fr1]=remoteforcevec(MU,NU,a,segments,0);
    
    %add force contributions together
    fseg=[fpk, fpk]*0.5+[fs0, fs1]+[fr0, fr1];
  

else
    %PK force due to applied stress
    fpk=pkforcevec(sigext,segments(linkid,:));

    %self force due to self stress
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments(linkid,:));

    %remote force due to remote stress
    [fr0,fr1]=remoteforcevec(MU,NU,a,segments,linkid);
   
    %add force contributions together
      fseg=[fpk, fpk]*0.5+[fs0, fs1]+[fr0, fr1];

  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segments=constructsegmentlist(rn,links)

[LINKMAX,m]=size(links);

segments=zeros(LINKMAX,14);
nseg=0;
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&(n1~=0))
        nseg=nseg+1;
        segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
    end
end
segments=segments(1:nseg,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=pkforcevec(sigext,segments)
%nodal force on dislocation segments due to applied stress sigext
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,nx,ny,nz)
[nseg,m] = size(segments);
f = zeros(nseg,3);

b = segments(:,3:5);
sigb = sigext*b';
r0 = segments(:,6:8);
r1 = segments(:,9:11);
r01 = r1-r0;

for i=1:nseg
    f(i,:) = cross(sigb(:,i)',r01(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1,f2]=selfforcevec(mu,nu,a,Ec,segments)
%self force (due to self stress)
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1)

    [nseg,m]=size(segments);
    Diff=segments(:,9:11)-segments(:,6:8);
 
    L=sqrt(sum(Diff.*Diff,2));
  if(min(L)>0)
    Linv=1./L;
    
    La=sqrt(L.*L+a*a);
    Lainv=1./La;
    
    t=zeros(nseg,3); %unite vector of segiment deirection
    for i=1:nseg
        t(i,:)=Diff(i,:)./norm(Diff(i,:));
    end
    
    omninv=1/(1-nu);
    
    bs=sum(segments(:,3:5).*t,2);
    bs2=bs.*bs;
    bev=segments(:,3:5)-[bs bs bs].*t;
    be2=sum(bev.*bev,2);
    
    % Elastic Self Interaction Force - Torsional Contribution
    S=(0.25*mu/pi).*bs.*((nu*omninv).*( log((La+L)./a)- 2.*(La-a).*Linv)- 0.5.*(La-a).*(La-a).*Linv.*Lainv);
    % Core Self Interaction Force - Torsional Contribution
    Score=2.*nu*omninv*Ec.*bs;
    Stot=S+Score;
    f2=[Stot Stot Stot].*bev;
    
    % Core Self Interaction Force - Longitudinal Component
    LTcore=(bs2 + be2.*omninv).*Ec;
    f2=f2-[LTcore LTcore LTcore].*t;
    f1=-f2;
  else
    for i=1:nseg
      f1(i,1:3)=0;
    end
    f2=-f1;
  end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f0,f1]=remoteforcevec(MU,NU,a,segments,linkid)

%nodal force on dislocation segment 01 due to another segment 23, for all
%segments

global USE_GPU;
S = size(segments,1);

%prepare inputs for MEX code.
bx = segments(:,3);
by = segments(:,4);
bz = segments(:,5);

p1x = segments(:,6); 
p1y = segments(:,7);
p1z = segments(:,8);

p2x = segments(:,9);
p2y = segments(:,10);
p2z = segments(:,11);

if (S<300 && USE_GPU==1) || (USE_GPU==0)
%MEX implementation if GPU option is OFF or if number of segments below 300.
%In this case, it's usually faster to run on a CPU...
    tic;[f0x, f0y, f0z, f1x, f1y, f1z] = SegSegForcesMex(p1x,p1y,p1z,...
                                                     p2x,p2y,p2z,...
                                                     bx,by,bz,...
                                                     a,MU,NU,...
                                                     linkid,S); 
else
    % Structure of Arrays; faster access data structure;
        SoA = reshape((segments(:,3:11))',[],1);
        f0x = zeros(S,1);
        f0y = zeros(S,1);
        f0z = zeros(S,1);
        f1x = zeros(S,1);
        f1y = zeros(S,1);
        f1z = zeros(S,1);
        
        % Initialize GPU kernel
        kern = parallel.gpu.CUDAKernel('SegForceNBodyCUDADoublePrecision.ptx', 'SegForceNBodyCUDADoublePrecision.cu');
        kern.ThreadBlockSize = 256;
        kern.GridSize = ceil(S/256);
     
        % Evaluate kernel on GPU device
        [f0x, f0y, f0z, f1x, f1y, f1z] = feval(kern,SoA,...
                                               a,MU,NU,...
                                               S,...
                                               f0x, f0y, f0z,...
                                               f1x, f1y, f1z);
        
        % Transfer results from device memory to host
        f0x=gather(f0x);
        f0y=gather(f0y);
        f0z=gather(f0z);
        f1x=gather(f1x);
        f1y=gather(f1y);
        f1z=gather(f1z);
end
    
if linkid==0                                             
    f0 = horzcat(f0x,f0y,f0z);
    f1 = horzcat(f1x,f1y,f1z);
else
    %when linkid is non-zero, SegSegForcesMex will swap linkid row with
    %first row of array. Therefore, the fseg of interest will be the first.
    f0 = [f0x(1),f0y(1),f0z(1)];
    f1 = [f1x(1),f1y(1),f1z(1)];
    
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%