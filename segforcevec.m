function [fseg_tot]=segforcevec(MU,NU,a,Ec,rn,links,linkid,~,...
    uhat,nc,xnodes,D,mx,mz,w,h,d)
%compute nodal driving force of dislocation network by stress formula
%(vectorized version)
%rn: nodal position
%links: dislocation segments (connectivity, Burgers vector)
%sigext: external stress (applied)
%[fs0,fs1,fr0,fr1,fimg,fseg,segments]
%segment contribution to node force

%NMAX: max number of nodes[NMAX,m]=size(rn);
[~,m]=size(rn);
if(m~=4)
    disp('rn should have 4 columns!');
    return;
end

%construct segment list
segments=constructsegmentlist(rn,links);

%find segments that are virtual with a surface node.
[virtindex,~,realindex,realID,~] = virtualsegfinder(rn,links,segments);
virtsegs = segments(virtindex,:);

S=size(segments,1);
%get only "real" segments, not virtual ones.
index=true(S,1);
for i=1:S
    if rn(segments(i,1),4) == 67 || rn(segments(i,2),4) == 67
        index(i) = false;
    end
end
segments = segments(index,:);

if any(any(isnan(segments)))
    disp('YDFUS, see line 36 segforcevec.m')
end

if linkid==0 %calculate forces on all segments
    %PK force due to applied stress
    %t0=clock;

    fpk=pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments);
    %t=etime(clock,t0); disp(sprintf('pkforcevec time = %6.2f seconds\n',t));

    %self force due to self stress
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments);
%     if not(isempty(virtsegs))
%         %zero out the self-force on surface nodes
%         realtosurfsegs=find(realindex==1);
%         for i=1:length(realtosurfsegs)
%             realtosurfsegs(i) = sum(index(1:realtosurfsegs(i))); %update for local indexing
%             if realID(i) == 1 %surface node is in first column of links entry
%                fs0(realtosurfsegs(i),:) = [0, 0, 0];
%             else %surface node is in second column of links entry
%                fs1(realtosurfsegs(i),:) = [0, 0, 0];
%             end
%         end
%     end

%     if isempty(virtsegs)
        %remote force due to remote stress, real segments
        [fr0,fr1]=remoteforcevec(MU,NU,a,segments,0);
%     else
%         %if virtual nodes present one must
%         %correct surface nodes by removing virtual contribution
%         %and delete appended forces on virtual nodes
%         segments_tot = vertcat(segments,virtsegs);
%         %force on every segment due to every other segment (virtual included)
%         [fr0,fr1]=remoteforcevec(MU,NU,a,segments_tot,0);
%         %subtract contribution in surface nodes of connected virtual nodes...
%         for vs = 1:size(segpairs,1)*0.5
%             segment_pair = segpairs(2*vs-1:2*vs,:); %first segment is real, second is virtual
%             [fr0v,fr1v]=remoteforcevec(MU,NU,a,segment_pair,1);
%             %decide which is surface node, and subtract "virtual" force
%             %from segment
%             rowpos=find(segments_tot(:,1)==segpairs(1,1) & segments_tot(:,2)==segpairs(1,2));
%             if realID(vs)==1 %first node is at surface
%                 fr0(rowpos,:) = fr0(rowpos,:) - fr0v;
%             elseif realID(vs)==2 %second node is at surface
%                 fr1(rowpos,:) = fr1(rowpos,:) - fr1v;
%             end
%         end
%         %delete appended forces on virtual nodes
%         fr0 = fr0(1:size(segments,1),:);
%         fr1 = fr1(1:size(segments,1),:);
%     end

    %PK force due to image stress
    %fimg = imageforce(segments,a,MU,NU,vertices);

    fseg=[fpk, fpk]*0.5+[fs0, fs1]+[fr0, fr1];   %combine all force contributions

%     if ~isempty(virtsegs)
%         %project forces on surface nodes into the plane of the surface
%         realtosurfsegs=find(realindex==1);
%         for i=1:length(realtosurfsegs)
%             realtosurfsegs(i) = sum(index(1:realtosurfsegs(i))); %update for local indexing
%             relevantsegs1=virtsegs(virtsegs(:,1)==segments(realtosurfsegs(i),realID(i)),:);
%             relevantsegs2=virtsegs(virtsegs(:,2)==segments(realtosurfsegs(i),realID(i)),:);
%             planenormal1=relevantsegs1(:,9:11)-relevantsegs1(:,6:8);
%             planenormal2=relevantsegs2(:,6:8)-relevantsegs2(:,9:11);
%             planenormal=[planenormal1;planenormal2];
%             planenormal=sum(planenormal,1);
%             planenormal=planenormal/norm(planenormal);
%             if realID(i) == 1 %surface node is in first column of links entry
%                 fseg(realtosurfsegs(i),1:3) = 2*(fseg(realtosurfsegs(i),1:3)-planenormal*dot(fseg(realtosurfsegs(i),1:3),planenormal));
%             else %surface node is in second column of links entry
%                 fseg(realtosurfsegs(i),4:6) = 2*(fseg(realtosurfsegs(i),4:6)-planenormal*dot(fseg(realtosurfsegs(i),4:6),planenormal));
%             end
%         end
%     end

    if any(any(isnan(fseg)))
        fseg(isnan(fseg))=0; %for when the collision code creates weird surface nodes
    end

    %fseg is calculated for real segments. The total fseg list used
    %throughout DDLab also contains virtual segments (flagged as 67).
    %Therefore, one must assign a nil force to the virtual segments (as
    %they are slaved to real ones) and remember to index the fseg elements
    %in the global list correctly. This is easy and fast since we have a
    %logic index list of real and virtual segments!
    fseg_tot = zeros(S,6);
    fseg_tot(index,1:6) = fseg;

else  %calculate force on segment specified by linkid

    %update linkid based on "real" index, rather than global indexing.
    if sum(index(1:linkid))==0
        fseg_tot=[0, 0, 0, 0, 0, 0]; %remesh.m wants to access virtual seg
        return;
    else
       linkid=sum(index(1:linkid));
    end

    %PK force due to applied stress
    %t0=clock;
    fpk=pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments(linkid,:));
    %t=etime(clock,t0); disp(sprintf('pkforcevec time = %6.2f seconds\n',t));

    %self force due to self stress
    [fs0,fs1]=selfforcevec(MU,NU,a,Ec,segments(linkid,:));

    %remote force due to remote stress
    %I think you do not need to include virtseg here because this is only called
    %by remesh.m to calculate the force on newly created nodes in real
    %segments when a node is created. length of [fr0,fr1]=1
    [fr0,fr1]=remoteforcevec(MU,NU,a,segments,linkid);

    fseg=[fpk, fpk]*0.5 + [fs0, fs1] + [fr0, fr1]; %sum force contributions

    %project forces on surface nodes into plane of surface
%     if rn(segments(linkid,1),end)==6
%        relevantsegs1=virtsegs(virtsegs(:,1)==segments(linkid,1),:);
%        relevantsegs2=virtsegs(virtsegs(:,2)==segments(linkid,1),:);
%        planenormal1=relevantsegs1(:,9:11)-relevantsegs1(:,6:8);
%        planenormal2=relevantsegs2(:,6:8)-relevantsegs2(:,9:11);
%        planenormal=[planenormal1;planenormal2];
%        planenormal=sum(planenormal,1);
%        planenormal=planenormal/norm(planenormal,1);
%        fseg(1:3) = 2*(fseg(1:3)-planenormal*dot(fseg(1:3),planenormal));
%     elseif rn(segments(linkid,2),end)==6
%        relevantsegs1=virtsegs(virtsegs(:,1)==segments(linkid,2),:);
%        relevantsegs2=virtsegs(virtsegs(:,2)==segments(linkid,2),:);
%        planenormal1=relevantsegs1(:,9:11)-relevantsegs1(:,6:8);
%        planenormal2=relevantsegs2(:,6:8)-relevantsegs2(:,9:11);
%        planenormal=[planenormal1;planenormal2];
%        planenormal=sum(planenormal,1);
%        planenormal=planenormal/norm(planenormal);
%        fseg(4:6) = 2*(fseg(4:6)-planenormal*dot(fseg(4:6),planenormal));
%     end

    if any(any(isnan(fseg)))
        fseg(isnan(fseg))=0; %for when the collision code creates weird surface nodes
    end

    fseg_tot = fseg; %output variable
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=pkforcevec(uhat,nc,xnodes,D,mx,mz,w,h,d,segments)
%nodal force on dislocation segments due to applied stress sigext
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1,nx,ny,nz)
[nseg,~]=size(segments);
f=zeros(nseg,3);

b=segments(:,3:5);
r0=segments(:,6:8);
r1=segments(:,9:11);
r01=r1-r0;
rmid = 0.5*(r0+r1);
for i =1:nseg % evaluate sigmahat at midpoint of each segment (do something better later)
    xi = rmid(i,1:3);

    if any(isnan(xi))
        disp('Segment midpoint is undefined! See segforcevec.m');
        pause;
    end

    sigext = hatStress(uhat,nc,xnodes,D,mx,mz,w,h,d,xi); %must not pass virtsegs!

    sigb=sigext*b(i,1:3)';
    l = r01(i,:);
    f(i,:)=cross(sigb',l);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1,f2]=selfforcevec(mu,nu,a,Ec,segments)
%self force (due to self stress)
%(vectorized version)
%format of segments array:
% (n0,n1,bx,by,bz,x0,y0,z0,x1,y1,z1)

    Diff=segments(:,9:11)-segments(:,6:8);

    L=sqrt(sum(Diff.*Diff,2));
    Linv=1./L;

    La=sqrt(L.*L+a*a);
    Lainv=1./La;

    t=Diff.*[Linv Linv Linv];

    omninv=1/(1-nu);

    bs=sum(segments(:,3:5).*t,2); %bs = b.lhat --> component of screw b
    bs2=bs.*bs;
    bev=segments(:,3:5)-[bs bs bs].*t; %bev = b - (b.lhat)lhat --> component edge b
    be2=sum(bev.*bev,2);

    % A. Arsenlis et al, Modelling Simul. Mater. Sci. Eng. 15 (2007)
    % 553?595: gives this expression in appendix A p590
    % f^s_43 = -(mu/4pi) [ t cross (t cross b)](t dot b) { v/(1-v) ( ln[
    % (L_a + L)/a] - 2*(L_a - a)/L ) - (L_a - a)^2/(2La*L) }

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
