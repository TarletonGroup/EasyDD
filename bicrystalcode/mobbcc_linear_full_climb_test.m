

function [v_c,fn_c] = mobbcc_linear_full_climb_test(fseg,rn,links,connectivity)

% test the nodal climb velocity, to find the nodes penetrating the inclusions

% vn_c climb velocity

% fn1 total nodal force

%mobility law function
% if nodelist and conlist are empty then vn is the length of rn otherwise
% vn is the length of nodelist


global amag tem mumag ;

vmax  = 3.9e11;        % unit:b/s,the max velocity is 100 m/s
eps   = 1e-12;
% bur=amag*sqrt(3)/2;
bur=amag;

% parameters for climb

% Coeff_core0 = 1.18e-5;
% formeng  = 0.67*1.6022e-19;%%%%%   formation energy of vacancy
% migeng   =0.6* 0.68*1.6022e-19; %        migration energy of vacancy
% Boltz    = 1.3807e-23;
% Coeff_core =Coeff_core0*exp(-migeng/Boltz/tem);  % core diffusivity;
%
% Atomvol  = 16.3e-30;  %%%%%%%       volum of the atom
% rc       = 5*bur*1e-9;   % radius of the dislocation core area

% Dc = 1e-3*Coeff_core*(2*pi*rc^2)*Atomvol/Boltz/tem/(bur*1e-9)^5; % effective diffusivity;


%%
pre_exponential = 1e-23; % m^4/s
active_energy = 174*1e3; % kJ/mol  1.8034eV
R = 8.314; % J/mol/K
Boltz    = 1.3807e-23; % J/K

D_core = pre_exponential*exp(-active_energy/R/tem);  % core diffusivity
Atomvol  = (2.856)^3*1e-30;  %%%%%%%       volum of the atom
% Dc= 1.0156e-43/(bur*1e-9)^5; % parameter from scientific report?
% Dc=D_core*Atomvol/Boltz/tem/(bur*1e-9)^5;
Dc=D_core*Atomvol/Boltz/tem;
% Dc=7.0488e-39;
%%

L1 = size(rn,1);    % when mobility law is call form remesh fucnction, L1~=L0;
nodelist = linspace(1,L1,L1)';
% flag=(rn(:,end)==67);
% nodelist(flag)=[];
[L2,L3] = size(connectivity);
conlist = zeros(L2,(L3-1)/2+1);
conlist(:,1) = connectivity(:,1);
for i=1:L2
    connumb = conlist(i,1);
    conlist(i,2:connumb+1) = linspace(1,connumb,connumb);
end

%% dislocation climb dominated by core diffusion  ( edited by Ally  12/11/2018 )
% flag = 1; % vacancy (1) or interstitial (-1)

v_c=zeros(L1,3);  % nodal climb velocity
fn_c=zeros(L1,3); % nodal climb force
lamla_c=zeros(L1,1); %chemical potential at each node

%% node force
nseg=size(links,1);
nnode=size(rn,1);
C=zeros(nnode,nnode+nseg);

k=zeros(nseg+nnode,nseg+nnode);
fnc=zeros(nnode,1);
fn_c=zeros(nnode,3);
K=zeros(nnode*2+nseg,nnode*2+nseg);
l=zeros(nseg,3);
nc=zeros(nnode,3);
fn=zeros(nnode,3);
Fnc=zeros(nnode+nseg,1);
F=zeros(nnode*2+nseg);

linvec=zeros(nnode,3);
linel=zeros(nnode,1);
vt=zeros(nnode,3); % line velocity;


cc=[-3/8,-1/8, -1;-1/8,-3/8, 1];
%%%%%%calculate the climb direction and climb force for each node
L=0;
for n=1:L1
    n0=nodelist(n);
    numNbrs=conlist(n,1);
    rt=zeros(numNbrs,3);
    
    for i=1:numNbrs
        ii=conlist(n,i+1);
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=(rn(n1,1:3)-rn(n0,1:3))*(-1)^(3-posinlink);
        
        
        linvec(n,:)= linvec(n,:)+rt/norm(rt)*(-1)^(3-posinlink);
        l(linkid)=norm(rt);
        linel(n)=linel(n)+l(linkid)*(-1)^(3-posinlink);
        
        
        %              normplane= links(linkid,6:8);
        normplane=cross(rt,links(linkid,3:5));
        if norm(normplane) > eps
            normplane=normplane/norm(normplane);
        else
            normplane=[0 0 0];
        end
        
        nc(n,:)=nc(n,:)+ normplane*l(linkid);
        %             L =L+ l(linkid);
        
        if l(linkid)>0
            %   fsegn0 =fseg(linkid,3*(posinlink-1)+[1:3])/l(linkid)*2;
            fsegn0 =fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:) = fn(n,:) + fsegn0;
        end
        
    end
    
    if norm(nc(n,:))> eps
        nc(n,:)=nc(n,:)/norm(nc(n,:));
    else

        nc(n,:)=0;
    end
    
    fnc(n,:)=dot(fn(n,:),nc(n,:));
    fn_c=fnc.*nc;
    linvec(n,:)=linvec(n,:)/norm(linvec(n,:));

    vt(n,:)=linel(n)*linvec(n,:);
    
end

%%%assemble the matrix
m=1; % the index of degrees
z=zeros(nseg+nnode,1); % track the index of the degrees;
loc_n0=0;
loc_n1=0;
loc_j=0;
k_e=zeros(3,3);
rnmax=size(rn,2);

for i=1:nseg
    n0=links(i,1);
    n1=links(i,2);
    
    
    r1=rn(links(i,1),1:3);
    r2=rn(links(i,2),1:3);
    rt=r2-r1;
    l(i)=norm(rt);
    burg=links(i,3:5);
    
    %  normplane= links(i,6:8);
    normplane=cross(rt,links(i,3:5));
    
    if norm(normplane) > eps
        normplane=normplane/norm(normplane);
    else
        normplane=[0 0 0];
    end
    
    an0= nc(n0,:)*normplane';
    an1= nc(n1,:)*normplane';
    
    if l(i)>0
        
        theta=dot(rt,burg)/norm(burg)/l(i);
        theta=acos(theta);
        sintheta=sin(theta);
    else
        sintheta=0;
    end
    
    if sintheta < 0.1  %  cutt off angle 10;
        
        sintheta=0.1;
        
    end
    
    
    %%% find the location of n0 n1 and middle node in the matrix;
    if ismember(n0,z)
        loc_n0=find(z==n0);
    else
        z(m)=n0;
        loc_n0=m;
        m=m+1;
    end
    z(m)=-i;
    loc_j=m;
    m=m+1;
    
    if ismember(n1,z)
        loc_n1=find(z==n1);
    else
        z(m)=n1;
        loc_n1=m;
        m=m+1;
    end
    
    k_e = l(i)*(amag*1e-6)/Dc* [23/960*(an0*l(i)*sintheta)^2*(amag*1e-6)^4,      17/960*(l(i)*sintheta)^2*an0*an1*(amag*1e-6)^4,   1/24*(an0*l(i)*sintheta)^1*(amag*1e-6)^2;
                               17/960*(l(i)*sintheta)^2*an0*an1*(amag*1e-6)^4,  23/960*(l(i)*sintheta*an1)^2*(amag*1e-6)^4,       -1/24*(an1*l(i)*sintheta)^1*(amag*1e-6)^2;
                               1/24*(an0*l(i)*sintheta)^1*(amag*1e-6)^2,        -1/24*(an1*l(i)*sintheta)^1*(amag*1e-6)^2,         1                           ];
    B_e = -[3/8*l(i)*sintheta*an0*(amag*1e-6)^2, 1/8*l(i)*sintheta*an1*(amag*1e-6)^2, 1;
            1/8*l(i)*sintheta*an0*(amag*1e-6)^2, 3/8*l(i)*sintheta*an1*(amag*1e-6)^2, -1];
    
    
    %        if sintheta < 0.1  %  cutt off angle 10;
    %            k_e=k_e*1e5;
    %            B_e=B_e*1e5;
    %        end
    % BC on pinned point
    %        if rn(n0,rnmax)~=0
    %            k_e(1,:)=k_e(1,:)*1e8;
    %            k_e(:,1)=k_e(:,1)*1e8;
    %        end
    %        if rn(n1,rnmax)~=0
    %            k_e(3,:)=k_e(3,:)*1e8;
    %            k_e(:,3)=k_e(:,3)*1e8;
    %        end
    % BC around the inclusions
    %      if nearinc_flag(n0)==1
    %          k_e(1,:)=k_e(1,:)*1e10;
    %          k_e(:,1)=k_e(:,1)*1e10;
    %      end
    %
    %      if nearinc_flag(n1)==1
    %          k_e(3,:)=k_e(3,:)*1e10;
    %          k_e(:,3)=k_e(:,3)*1e10;
    %      end
    
    
    
    
    % assemble elementary k
    k(loc_n0, loc_n0) = k(loc_n0, loc_n0) + k_e(1,1);
    k(loc_n0, loc_n1) = k(loc_n0, loc_n1) + k_e(1,2);
    k(loc_n0, loc_j)  = k(loc_n0, loc_j)  + k_e(1,3);
    k(loc_n1, loc_n0) = k(loc_n1, loc_n0) + k_e(2,1);
    k(loc_n1, loc_n1) = k(loc_n1, loc_n1) + k_e(2,2);
    k(loc_n1, loc_j)  = k(loc_n1, loc_j)  + k_e(2,3);
    k(loc_j,  loc_n0) = k(loc_j, loc_n0)  + k_e(3,1);
    k(loc_j,  loc_n1) = k(loc_j, loc_n1)  + k_e(3,2);
    k(loc_j,  loc_j)  = k(loc_j, loc_j)   + k_e(3,3);
    
    C(n0,loc_n0) =  C(n0,loc_n0) + B_e(1,1);
    C(n0,loc_n1) =  C(n0,loc_n1) + B_e(1,2);
    C(n0,loc_j)  =  C(n0,loc_j)  + B_e(1,3);
    C(n1,loc_n0) =  C(n1,loc_n0) + B_e(2,1);
    C(n1,loc_n1) =  C(n1,loc_n1) + B_e(2,2);
    C(n1,loc_j)  =  C(n1,loc_j)  + B_e(2,3);
    
    Fnc(loc_n0) = fnc(n0);
    Fnc(loc_n1) = fnc(n1);
    Fnc(loc_j)  = 0;
    
    
end


% reshape the matrix
lamla=zeros(nnode,nnode);
%     K=[k C';C lamla];
K =vertcat(horzcat(k,C'), horzcat(C,lamla));
F=[Fnc*(mumag*1e6)*(amag*1e-6)^1;zeros(nnode,1)];

%     tic
KK=sparse(K);
FF=sparse(F);
vc=KK\FF;  % mumag*amag a/s
%     toc

% tic
%      [vc,R]=linsolve(K,F); %% solve K*vc=F to get vc, R should be <<0 to ensure the accuracy of the solution;
%      toc
%    spy(K)  % show the structure of K
if sum( any(isnan(vc)==1))
    disp('g sum( any(isnan(vc)==1)) in climb mobility law line 300');
    %        pause
end
% %
%     string=[1:nnode];
%     string1=strtrim(cellstr(num2str(string'))');
%  plim = 1500;
% viewangle = [135 45];
% plotnodes(rn,links,plim);
% view(viewangle);
% xlim([-plim plim]);
% ylim([-plim plim]);
% zlim([-plim plim]);
% hold on
for n=1:nnode
    nnn=nc(n,:)/norm(nc(n,:));
    v_c(n,:)= vc(find(z==n))*nc(n,:);
    %
    %        lamda(n,:)=vc(size(z,1)+n)*nc(n,:)*(bur*1e-9)^2;  % chemical force; N/m^2
    %        fc(n,:)=fn_c(n,:)*(bur*1e-9)^2;
    %        text(rn(n,1),rn(n,2),rn(n,3),string1(n));
    %        burgg = links(n,3:5);
    %        normplane = links(n,6:8);
    %        quiver3(rn(n,1),rn(n,2),rn(n,3),burgg(1),burgg(2),burgg(3),1e2,'-m','LineWidth',2); % burgers
    %        quiver3(rn(n,1),rn(n,2),rn(n,3),normplane(1),normplane(2),normplane(3),1e2,'-c','LineWidth',2);
    %        quiver3(rn(n,1),rn(n,2),rn(n,3),lamda(n,1),lamda(n,2),lamda(n,3),1e3,'-b','LineWidth',2); %chemical
    %        quiver3(rn(n,1),rn(n,2),rn(n,3),fc(n,1),fc(n,2),fc(n,3),1e11,'-r','LineWidth',2); % mechanical
    %        quiver3(rn(n,1),rn(n,2),rn(n,3),v_c(n,1),v_c(n,2),v_c(n,3),1e3,'-k','LineWidth',2); % mechanical
    
end


%                 quiver3(rn(n,1),rn(n,2),rn(n,3),burgg(1),burgg(2),burgg(3),1e2,'-m','LineWidth',2);
%                 quiver3(rn(n,1),rn(n,2),rn(n,3),normplane(1),normplane(2),normplane(3),1e2,'-c','LineWidth',2);
%                 quiver3(rn(n,1),rn(n,2),rn(n,3),rt(1),rt(2),rt(3),1,'-b','LineWidth',2);

%                 text(rn(n,1),rn(n,2),rn(n,3),string1(n));


%     lamla_c=vc(3*nseg+1:size(F,1),1); % chemical force

rnmax=size(rn,2);
for i=1:size(rn,1)
    if rn(i,rnmax)~=0
        v_c(i,1:3)=0;
    end
end

%       v_c=v_c + norm(max(v_c))/norm(max(vt))*1*vt;
% %   v_c=vt;






